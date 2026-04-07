#!/usr/bin/Rscript
# build_cor_matrix.R
#
# Precompute the block-diagonal sparse gene-gene correlation matrix for a TWAS
# expression panel. The result is reusable across many TWAS-GSEA-fast.R runs
# (different GWAS, different gene-set / property files) on the same panel.

start.time <- Sys.time()
suppressMessages(library(optparse))

option_list <- list(
	make_option('--expression_ref', action='store', default=NA, type='character',
		help='Predicted expression reference file (gz allowed) [required]'),
	make_option('--pos', action='store', default=NA, type='character',
		help='FUSION .pos file with WGT, ID, CHR, P0, P1 [required]'),
	make_option('--cor_window', action='store', default=5e6, type='numeric',
		help='Window for retaining gene-gene correlations [default 5e6]'),
	make_option('--min_r2', action='store', default=0.0001, type='numeric',
		help='r^2 threshold below which correlations are zeroed [default 1e-4]'),
	make_option('--max_r2', action='store', default=1, type='numeric',
		help='r^2 threshold above which collinear genes are pruned [default 1]'),
	make_option('--n_cores', action='store', default=1, type='numeric',
		help='Cores for foreach over blocks [default 1]'),
	make_option('--output', action='store', default=NA, type='character',
		help='Output prefix; writes <output>.CorMat.RDS and <output>.CorMat.log [required]')
)
opt <- parse_args(OptionParser(option_list = option_list))

for(req in c('expression_ref','pos','output')){
	if(is.na(opt[[req]])) stop('--', req, ' is required')
}

if(dirname(opt$output) != '.') system(paste0('mkdir -p ', dirname(opt$output)))

LOG_FILE <- paste0(opt$output, '.CorMat.log')
log_msg  <- function(...) cat(..., file = LOG_FILE, append = TRUE)
cat('', file = LOG_FILE)
log_progress <- function(i, n){
	pct <- floor(i / n * 100)
	if(pct %in% seq(10, 100, 10) && i == floor(n / 100 * pct)) log_msg(pct, '% ', sep = '')
}

suppressMessages(library(data.table))
sink('/dev/null'); suppressMessages(library(WGCNA, quietly = TRUE)); sink()
suppressMessages(library(Matrix))
suppressMessages(library(matrixcalc))
suppressMessages(library(foreach))
suppressMessages(library(doMC))
registerDoMC(opt$n_cores)

source(file.path(dirname(sub('--file=', '', grep('--file=', commandArgs(trailingOnly = FALSE), value = TRUE)[1])), 'R', 'build_cor_matrix_helper.R'))

log_msg('build_cor_matrix.R\nOptions:\n')
log_msg(paste(capture.output(print(opt)), collapse = '\n'), '\n')
log_msg('Started at ', as.character(start.time), '\n')

# ---------------------------------------------------------------------------
# Load expression matrix.
# ---------------------------------------------------------------------------
if(substr(opt$expression_ref, nchar(opt$expression_ref) - 2, nchar(opt$expression_ref)) == '.gz'){
	GeneX_all <- data.frame(fread(cmd = paste0('zcat ', opt$expression_ref)))
} else {
	GeneX_all <- data.frame(fread(opt$expression_ref))
}
GeneX_all <- GeneX_all[-1:-2]
GeneX_all <- GeneX_all[, apply(GeneX_all, 2, function(x) !(var(x) == 0 | all(is.na(x))))]
names(GeneX_all) <- gsub(':', '.', names(GeneX_all))
names(GeneX_all) <- gsub('-', '.', names(GeneX_all))
log_msg('Expression panel:', ncol(GeneX_all), 'non-zero-variance features,', nrow(GeneX_all), 'individuals.\n')

# ---------------------------------------------------------------------------
# Load .pos and derive FILE column the same way TWAS-GSEA.V1.2.R does (so the
# rownames of the cor matrix match what TWAS-GSEA-fast.R will see in TWAS
# results FILE columns).
# ---------------------------------------------------------------------------
pos <- data.frame(fread(opt$pos))
if(!all(c('WGT','CHR','P0','P1') %in% names(pos))) stop('.pos file must contain WGT, CHR, P0, P1 columns')
pos$FILE <- sub('.wgt.RDat', '', sub('.*/', '', pos$WGT))
pos$FILE <- gsub(':', '.', pos$FILE)
pos$FILE <- gsub('-', '.', pos$FILE)

# Apply the same +/-5e5 SNP-window padding V1.2 uses (V1.2:248-250) so block
# construction is identical between the two pipelines.
pos$P0 <- pmax(pos$P0 - 5e5, 0)
pos$P1 <- pos$P1 + 5e5

# Intersect with expression panel.
genes_overlap <- intersect(pos$FILE, names(GeneX_all))
if(length(genes_overlap) == 0) stop('No overlap between .pos WGT names and expression matrix columns.')
genes_df <- pos[pos$FILE %in% genes_overlap, c('FILE','CHR','P0','P1')]
genes_df <- genes_df[!duplicated(genes_df$FILE), ]
genes_df <- genes_df[order(genes_df$CHR, genes_df$P0, genes_df$P1), ]
GeneX_all <- GeneX_all[, match(genes_df$FILE, names(GeneX_all)), drop = FALSE]
log_msg(nrow(genes_df), 'genes overlap between .pos and expression panel.\n')

# ---------------------------------------------------------------------------
# Build cor matrix and save.
# ---------------------------------------------------------------------------
cor_result <- build_cor_matrix(
	genes_df     = genes_df,
	GeneX_all    = GeneX_all,
	cor_window   = opt$cor_window,
	min_r2       = opt$min_r2,
	max_r2       = opt$max_r2,
	log_msg      = log_msg,
	log_progress = log_progress)
cor_block_all <- cor_result$K

saveRDS(cor_result, paste0(opt$output, '.CorMat.RDS'))

end.time <- Sys.time()
log_msg('Wrote ', opt$output, '.CorMat.RDS (', dim(cor_block_all)[1], ' x ', dim(cor_block_all)[2], ').\n', sep = '')
log_msg('Finished at ', as.character(end.time), ' (elapsed ', round(as.numeric(difftime(end.time, start.time, units = 'secs')), 1), 's)\n', sep = '')
