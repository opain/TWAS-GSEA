# sim/R/feature_builders.R
#
# Build the aligned feature matrix (rows = TWAS FILE ids, cols = gene sets or
# properties) that feeds sim/R/in_session_core.R::score_dataset(). Reproduces
# the ID-matching + filtering logic in TWAS-GSEA-fast.R:200-245.

suppressMessages({
  library(data.table)
  library(qusage)
})

# From a gmt file + a TWAS data.frame (with an Alt_ID column carrying the same
# ID space as the gmt genes), build a 0/1 membership matrix with rows named by
# TWAS$FILE and columns named by gene set (punctuation-normalised).
build_membership_matrix <- function(twas_df, gmt_path, alt_id_col) {
  stopifnot(alt_id_col %in% names(twas_df))
  gs <- read.gmt(gmt_path)
  names(gs) <- gsub('[[:punct:]]', '.', names(gs))
  alt_ids <- as.character(twas_df[[alt_id_col]])
  mem <- vapply(gs, function(members) alt_ids %in% as.character(members),
                logical(length(alt_ids)))
  storage.mode(mem) <- 'double'
  rownames(mem) <- twas_df$FILE
  attr(mem, 'set_sizes') <- vapply(gs, length, integer(1))
  mem
}

# Build a numeric property matrix aligned to TWAS FILE ids. `prop_mat` is
# rownames-keyed by the same ID space as `twas_df[[alt_id_col]]` (e.g. Symbol
# or Entrez).
build_property_matrix <- function(twas_df, prop_mat, alt_id_col) {
  stopifnot(alt_id_col %in% names(twas_df))
  alt_ids <- as.character(twas_df[[alt_id_col]])
  join_idx <- match(alt_ids, rownames(prop_mat))
  keep <- !is.na(join_idx)
  if(sum(keep) == 0) stop('no overlap between twas_df$', alt_id_col, ' and prop_mat rownames')
  out <- prop_mat[join_idx[keep], , drop = FALSE]
  storage.mode(out) <- 'double'
  out[!is.finite(out)] <- 0
  rownames(out) <- twas_df$FILE[keep]
  out
}

# Load an .rds or text property file into a matrix, matching the CLI logic at
# TWAS-GSEA-fast.R:215-224.
load_prop_file <- function(path) {
  if(grepl('\\.rds$', path, ignore.case = TRUE)){
    prop_mat <- readRDS(path)
  } else {
    dt <- fread(path)
    id <- dt[[1]]
    prop_mat <- as.matrix(dt[, -1, with = FALSE])
    rownames(prop_mat) <- id
  }
  storage.mode(prop_mat) <- 'double'
  prop_mat
}
