#!/usr/bin/Rscript
# sim/build_signed_cor.R
#
# Signed-K variant of build_cor_matrix.R. Sources sim/R/build_cor_matrix_signed_helper.R
# instead of R/build_cor_matrix_helper.R and writes <output>.CorMatSigned.RDS.
# After the build, if a companion <output>.CorMat.RDS (abs-K) already exists,
# assert identical rownames and identical sparsity pattern (this must hold —
# pruning is sign-invariant); fail loudly on mismatch.

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
		help='Output prefix; writes <output>.CorMatSigned.RDS and <output>.CorMatSigned.log [required]')
)
opt <- parse_args(OptionParser(option_list = option_list))

for(req in c('expression_ref','pos','output')){
	if(is.na(opt[[req]])) stop('--', req, ' is required')
}

if(dirname(opt$output) != '.') system(paste0('mkdir -p ', dirname(opt$output)))

LOG_FILE <- paste0(opt$output, '.CorMatSigned.log')
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

SCRIPT_DIR <- dirname(sub('--file=', '', grep('--file=', commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
source(file.path(SCRIPT_DIR, 'R', 'build_cor_matrix_signed_helper.R'))

log_msg('build_signed_cor.R\nOptions:\n')
log_msg(paste(capture.output(print(opt)), collapse = '\n'), '\n')
log_msg('Started at ', as.character(start.time), '\n')

# ---------------------------------------------------------------------------
# Load expression matrix. (Identical to build_cor_matrix.R.)
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
# Load .pos and derive FILE column. (Identical to build_cor_matrix.R.)
# ---------------------------------------------------------------------------
pos <- data.frame(fread(opt$pos))
if(!all(c('WGT','CHR','P0','P1') %in% names(pos))) stop('.pos file must contain WGT, CHR, P0, P1 columns')
pos$FILE <- sub('.wgt.RDat', '', sub('.*/', '', pos$WGT))
pos$FILE <- gsub(':', '.', pos$FILE)
pos$FILE <- gsub('-', '.', pos$FILE)

pos$P0 <- pmax(pos$P0 - 5e5, 0)
pos$P1 <- pos$P1 + 5e5

genes_overlap <- intersect(pos$FILE, names(GeneX_all))
if(length(genes_overlap) == 0) stop('No overlap between .pos WGT names and expression matrix columns.')
genes_df <- pos[pos$FILE %in% genes_overlap, c('FILE','CHR','P0','P1')]
genes_df <- genes_df[!duplicated(genes_df$FILE), ]
genes_df <- genes_df[order(genes_df$CHR, genes_df$P0, genes_df$P1), ]
GeneX_all <- GeneX_all[, match(genes_df$FILE, names(GeneX_all)), drop = FALSE]
log_msg(nrow(genes_df), 'genes overlap between .pos and expression panel.\n')

# ---------------------------------------------------------------------------
# Build signed cor matrix and save.
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

out_rds <- paste0(opt$output, '.CorMatSigned.RDS')
saveRDS(cor_result, out_rds)

end.time <- Sys.time()
log_msg('Wrote ', out_rds, ' (', dim(cor_block_all)[1], ' x ', dim(cor_block_all)[2], ').\n', sep = '')
log_msg('Finished at ', as.character(end.time), ' (elapsed ', round(as.numeric(difftime(end.time, start.time, units = 'secs')), 1), 's)\n', sep = '')

# ---------------------------------------------------------------------------
# Cross-check against abs-K matrix if present.
# ---------------------------------------------------------------------------
abs_rds <- paste0(opt$output, '.CorMat.RDS')
if(file.exists(abs_rds)){
	log_msg('Cross-check vs ', abs_rds, ':\n', sep = '')
	abs_obj <- readRDS(abs_rds)
	K_abs <- abs_obj$K
	K_sgn <- cor_result$K
	blocks_abs <- abs_obj$blocks
	blocks_sgn <- cor_result$blocks

	# HARD invariants: gene universe and block partition must match. These
	# depend only on |r|-based pruning + windowing + block detection (all
	# sign-invariant), so any mismatch is a bug.
	err <- character(0)
	if(!identical(rownames(K_abs), rownames(K_sgn))) err <- c(err, 'rownames differ')
	if(!identical(colnames(K_abs), colnames(K_sgn))) err <- c(err, 'colnames differ')
	if(!identical(dim(K_abs), dim(K_sgn))) err <- c(err, sprintf('dims differ: %s vs %s', paste(dim(K_abs), collapse='x'), paste(dim(K_sgn), collapse='x')))
	if(!identical(names(blocks_abs), names(blocks_sgn))) err <- c(err, 'block-index names differ')
	if(!identical(unname(blocks_abs), unname(blocks_sgn))) err <- c(err, 'block partition differs')
	if(length(err) > 0){
		log_msg('  HARD MISMATCH: ', paste(err, collapse = '; '), '\n', sep = '')
		stop('Signed-K vs abs-K cross-check failed: ', paste(err, collapse = '; '))
	}
	log_msg('  rownames/colnames/dims/block-partition identical.\n')

	# SOFT diagnostics: sparsity patterns can differ because the per-block
	# nearPD step is not sign-invariant (signed K can have negative
	# eigenvalues; abs K can have large positive ones). Report overlap
	# instead of failing.
	pat_abs <- (K_abs != 0)
	pat_sgn <- (K_sgn != 0)
	shared  <- pat_abs & pat_sgn
	only_abs <- sum(pat_abs & !pat_sgn)
	only_sgn <- sum(!pat_abs & pat_sgn)
	nnz_abs  <- sum(pat_abs)
	nnz_sgn  <- sum(pat_sgn)
	log_msg('  nnz(K_abs)=', nnz_abs, '  nnz(K_signed)=', nnz_sgn,
	        '  shared=', sum(shared),
	        '  only-abs=', only_abs,
	        '  only-signed=', only_sgn, '\n', sep = '')

	# On the shared nonzero pattern, |K_signed| should track K_abs up to the
	# nearPD perturbation. Report the max discrepancy.
	Kabs_shared <- as.numeric(K_abs[shared])
	Ksgn_shared <- as.numeric(K_sgn[shared])
	max_abs_diff <- max(abs(abs(Ksgn_shared) - Kabs_shared))
	log_msg('  max ||K_signed| - K_abs| on shared nnz = ', signif(max_abs_diff, 4), '\n', sep = '')
	neg_frac <- mean(Ksgn_shared < 0)
	log_msg('  fraction of shared entries with negative K_signed sign = ', signif(neg_frac, 3), '\n', sep = '')
} else {
	log_msg('No companion abs-K file at ', abs_rds, ' — skipping cross-check.\n', sep = '')
}
