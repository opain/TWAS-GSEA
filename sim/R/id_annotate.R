# sim/R/id_annotate.R
#
# Adds `Symbol` and `Entrez` columns to a TWAS data.frame keyed on `ID`
# (unversioned ENSG). Deduped 1:1 maps loaded once and cached in the calling
# process via `local()`.

suppressMessages({
  library(data.table)
})

.load_annotate_maps <- local({
  cache <- NULL
  function(ensg2entrez_path, biomart_path){
    key <- paste(ensg2entrez_path, biomart_path, sep = '|')
    if(!is.null(cache) && identical(cache$key, key)) return(cache$maps)

    # ENSG -> Entrez (from NCBI gene2ensembl, pre-filtered).
    e2e <- fread(ensg2entrez_path)
    stopifnot(all(c('ENSG','Entrez') %in% names(e2e)))
    e2e <- unique(e2e[, .(ENSG, Entrez)])
    setorder(e2e, ENSG, Entrez)
    e2e <- e2e[, .(Entrez = Entrez[1]), by = ENSG]

    # ENSG -> Symbol (from biomart TSV).
    bm <- fread(biomart_path)
    stopifnot(all(c('ensembl_gene_id','external_gene_name') %in% names(bm)))
    bm <- unique(bm[, .(ensembl_gene_id, external_gene_name)])
    bm <- bm[external_gene_name != '' & !is.na(external_gene_name)]
    setorder(bm, ensembl_gene_id, external_gene_name)
    sym <- bm[, .(Symbol = external_gene_name[1]), by = ensembl_gene_id]
    setnames(sym, 'ensembl_gene_id', 'ENSG')

    maps <- list(entrez = e2e, symbol = sym)
    cache <<- list(key = key, maps = maps)
    maps
  }
})

annotate_twas <- function(twas_df,
                          ensg2entrez_path = 'sim/output/id_map/ensg2entrez.tsv',
                          biomart_path     = '/data/biomart/biomart_genes_grch37.tsv') {
  # twas_df must have an `ID` column of unversioned ENSG (matches Whole_Blood.pos).
  stopifnot('ID' %in% names(twas_df))
  maps <- .load_annotate_maps(ensg2entrez_path, biomart_path)

  dt <- as.data.table(twas_df)
  # Left join preserving twas_df row order.
  dt[, `:=`(row_ord = seq_len(.N))]
  dt <- maps$entrez[dt, on = c(ENSG = 'ID')]
  setnames(dt, 'ENSG', 'ID')
  dt <- maps$symbol[dt, on = c(ENSG = 'ID')]
  setnames(dt, 'ENSG', 'ID')
  setorder(dt, row_ord)
  dt[, row_ord := NULL]
  as.data.frame(dt)
}
