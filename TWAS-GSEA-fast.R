#!/usr/bin/Rscript
# TWAS-GSEA-fast.R
#
# Fast standalone competitive gene-set / gene-property analysis for TWAS results.
#
# Statistical model is identical to TWAS-GSEA.V1.2.R:
#     y = X beta + g + e,  g ~ N(0, sigma_u^2 K),  e ~ N(0, sigma_e^2 I)
# where K is the block-diagonal predicted-expression correlation matrix
# precomputed once per expression panel by build_cor_matrix.R.
#
# Speed comes from two changes vs V1.2:
#   1. K is loaded from a precomputed .CorMat.RDS rather than rebuilt per run.
#   2. (sigma_u^2, sigma_e^2) are fitted by a custom REML routine that exploits
#      K's block-diagonal structure via per-block eigendecomposition + a 1-D
#      profiled likelihood (Brent's method). Mathematically equivalent to
#      lme4qtl::relmatLmer on the same model — not an approximation — but much
#      faster, with no lme4qtl dependency.
#
# After the null fit, the per-set GLS pipeline is the same as V1.2's
# --fast_competitive T path: sparse Cholesky of V_hat, whiten y / X / Z_gs,
# QR-residualise against the null design, vectorised crossprod across all
# gene sets.

start.time <- Sys.time()
suppressMessages(library(optparse))

option_list <- list(
	make_option('--twas_results', action='store', default=NA, type='character',
		help='TWAS results file (gz allowed) [required]'),
	make_option('--input_CorMat', action='store', default=NA, type='character',
		help='.CorMat.RDS produced by build_cor_matrix.R [required]'),
	make_option('--pos', action='store', default=NA, type='character',
		help='Optional FUSION .pos file. When provided, P0/P1 in --twas_results are overridden from .pos (FUSION bug fix), the +/-5e5 SNP-window padding is applied, and a GeneLength = P1-P0 column is computed so it can be passed to --covar [default NA]'),
	make_option('--gmt_file',  action='store', default=NA, type='character',
		help='Gene set file in gmt format [optional]'),
	make_option('--prop_file', action='store', default=NA, type='character',
		help='Gene property file (first column ID) [optional]'),
	make_option('--covar', action='store', default='none', type='character',
		help='Comma-separated covariate columns from --twas_results [default none]'),
	make_option('--use_alt_id', action='store', default=NA, type='character',
		help='Alt ID column in --twas_results to match gene set / property IDs (e.g. ID) [optional]'),
	make_option('--gene_id_map', action='store', default=NA, type='character',
		help='Gene symbol -> Entrez ID map (skips BioMart) [optional]'),
	make_option('--directional', action='store', default=FALSE, type='logical',
		help='Use signed TWAS.Z as outcome [default FALSE]'),
	make_option('--probit_P_as_Z', action='store', default=TRUE, type='logical',
		help='Use probit(1-TWAS.P) as outcome (overridden by --directional) [default TRUE]'),
	make_option('--two_sided', action='store', default=FALSE, type='logical',
		help='Two-sided p-values [default FALSE]'),
	make_option('--outlier_threshold', action='store', default='-3,6', type='character',
		help='Lower,upper z-score truncation [default -3,6]'),
	make_option('--min_Ngenes', action='store', default=2, type='numeric',
		help='Minimum number of available genes per set/property [default 2]'),
	make_option('--allow_duplicate_ID', action='store', default=FALSE, type='logical',
		help='Keep duplicate IDs (otherwise retain best MODELCV.R2) [default FALSE]'),
	make_option('--h_max', action='store', default=100, type='numeric',
		help='Upper bound on the variance ratio h = sigma_u^2/sigma_e^2 in REML search [default 100]'),
	make_option('--reml_tol', action='store', default=1e-6, type='numeric',
		help='Tolerance for the 1-D Brent search in the REML fit [default 1e-6]'),
	make_option('--twas_p_thresh', action='store', default='1', type='character',
		help='Comma-separated TWAS p-value thresholds for subsetting genes; 1 = all genes [default 1]'),
	make_option('--p_cor_method', action='store', default='fdr', type='character',
		help='Multiple testing correction method (passed to p.adjust) [default fdr]'),
	make_option('--n_cores', action='store', default=1, type='numeric',
		help='Cores [default 1]'),
	make_option('--output', action='store', default=NA, type='character',
		help='Output prefix; writes <output>.competitive.txt and <output>.log [required]')
)
opt <- parse_args(OptionParser(option_list = option_list))

for(req in c('twas_results','input_CorMat','output')){
	if(is.na(opt[[req]])) stop('--', req, ' is required')
}
if(is.na(opt$gmt_file) && is.na(opt$prop_file)) stop('Either --gmt_file or --prop_file must be specified')
if(!is.na(opt$gmt_file) && !is.na(opt$prop_file)) stop('Specify only one of --gmt_file / --prop_file')

if(dirname(opt$output) != '.') system(paste0('mkdir -p ', dirname(opt$output)))

opt$covar             <- as.character(unlist(strsplit(opt$covar, ',')))
opt$outlier_threshold <- as.numeric(unlist(strsplit(opt$outlier_threshold, ',')))
opt$twas_p_thresh     <- sort(unique(as.numeric(unlist(strsplit(opt$twas_p_thresh, ',')))))

LOG_FILE <- paste0(opt$output, '.log')
log_msg  <- function(...) cat(..., file = LOG_FILE, append = TRUE)
cat('', file = LOG_FILE)

suppressMessages(library(data.table))
suppressMessages(library(Matrix))
suppressMessages(library(VGAM))
suppressMessages(library(qusage))
suppressMessages(library(foreach))
suppressMessages(library(doMC))
registerDoMC(opt$n_cores)

# Sourced helpers (kept alongside the script).
SCRIPT_DIR <- dirname(sub('--file=', '', grep('--file=', commandArgs(trailingOnly=FALSE), value=TRUE)[1]))
source(file.path(SCRIPT_DIR, 'R', 'reml_blockdiag.R'))

# Reuse compute_pval helper (mirror of TWAS-GSEA.V1.2.R lines 175-183).
compute_pval <- function(t_val, two_sided){
	if(two_sided) 2 * pnorm(-abs(t_val)) else 1 - pnorm(t_val)
}

log_msg(
'#################################################################
# TWAS-GSEA-fast: blockwise-REML competitive GSEA
#################################################################
Options:\n')
log_msg(paste(capture.output(print(opt)), collapse = '\n'), '\n')
log_msg('Started at ', as.character(start.time), '\n')

# ---------------------------------------------------------------------------
# 1. Load TWAS results (mirrors TWAS-GSEA.V1.2.R:230-314).
# ---------------------------------------------------------------------------
TWAS <- data.frame(fread(opt$twas_results))
log_msg('TWAS results: ', nrow(TWAS), ' rows.\n', sep = '')

# Stage 1: normalise FILE to "Panel/Gene.wgt.RDat" form (last two path
# components). This matches the WGT column in FUSION .pos files.
file_list <- strsplit(as.character(TWAS$FILE), '/')
file_tab  <- lapply(file_list, function(x) x[(length(x) - 1):length(x)])
tmp <- data.frame(do.call(rbind, file_tab))
TWAS$FILE <- do.call(paste, c(tmp[, (ncol(tmp) - 1):ncol(tmp)], sep = '/'))

# Optional .pos merge: overrides P0/P1 (FUSION bug fix), pads by +/-5e5 (the
# FUSION SNP window), and creates GeneLength so it can be used as a covariate.
# Mirrors TWAS-GSEA.V1.2.R:241-253 exactly. Must happen BEFORE further FILE
# normalisation because pos$WGT is in the "Panel/Gene.wgt.RDat" form.
if(!is.na(opt$pos)){
	pos <- data.frame(fread(opt$pos))
	if(!all(c('WGT','P0','P1') %in% names(pos))) stop('.pos file must contain WGT, P0, P1 columns')
	TWAS$P0 <- NULL
	TWAS$P1 <- NULL
	TWAS <- merge(TWAS, pos[, c('WGT','P0','P1')], by.x = 'FILE', by.y = 'WGT')
	log_msg('Positional information available for ', nrow(TWAS), ' TWAS features after .pos merge.\n', sep = '')
	TWAS$P0 <- TWAS$P0 - 5e5
	TWAS$P0[TWAS$P0 < 0] <- 0
	TWAS$P1 <- TWAS$P1 + 5e5
	TWAS$GeneLength <- TWAS$P1 - TWAS$P0
}

# Stage 2: strip the path prefix and .wgt.RDat suffix and normalise punctuation
# so FILE matches the gene IDs in the cor matrix and gene-set/property files.
TWAS$FILE <- sub('.*/', '', TWAS$FILE)
TWAS$FILE <- sub('.wgt.RDat', '', TWAS$FILE)
TWAS$FILE <- gsub(':', '.', TWAS$FILE)
TWAS$FILE <- gsub('-', '.', TWAS$FILE)

TWAS <- TWAS[!is.na(TWAS$TWAS.P), ]
TWAS <- TWAS[!duplicated(TWAS$FILE), ]

if(!all(c('FILE','ID','TWAS.Z','TWAS.P') %in% names(TWAS))){
	stop('--twas_results must contain columns FILE, ID, TWAS.Z, TWAS.P')
}

# Alt ID handling.
if(!is.na(opt$use_alt_id)){
	if(opt$use_alt_id == 'ID'){
		TWAS$Alt_ID <- TWAS$ID
	} else {
		if(!(opt$use_alt_id %in% names(TWAS))) stop('--use_alt_id column not found')
		names(TWAS)[names(TWAS) == opt$use_alt_id] <- 'Alt_ID'
	}
}

# Drop duplicate IDs by best MODELCV.R2.
if(!opt$allow_duplicate_ID){
	if(!('MODELCV.R2' %in% names(TWAS))) stop('TWAS results must contain MODELCV.R2 column')
	TWAS <- TWAS[order(TWAS$MODELCV.R2), ]
	if(is.na(opt$use_alt_id)) TWAS <- TWAS[!duplicated(TWAS$ID), ] else TWAS <- TWAS[!duplicated(TWAS$Alt_ID), ]
}

# ZSCORE outcome.
if(opt$directional){
	TWAS$ZSCORE <- TWAS$TWAS.Z
	TWAS$ZSCORE[TWAS$ZSCORE > opt$outlier_threshold[2]] <- opt$outlier_threshold[2]
	TWAS$ZSCORE[TWAS$ZSCORE < opt$outlier_threshold[1]] <- opt$outlier_threshold[1]
} else if(opt$probit_P_as_Z){
	TWAS$ZSCORE <- probitlink(1 - TWAS$TWAS.P)
	TWAS$ZSCORE[TWAS$ZSCORE > opt$outlier_threshold[2]] <- opt$outlier_threshold[2]
	TWAS$ZSCORE[TWAS$ZSCORE < opt$outlier_threshold[1]] <- opt$outlier_threshold[1]
} else {
	TWAS$ZSCORE <- abs(TWAS$TWAS.Z)
	TWAS$ZSCORE[TWAS$ZSCORE > opt$outlier_threshold[2]] <- opt$outlier_threshold[2]
}

# Symbol -> entrez via gene_id_map (only used when --use_alt_id is not set).
if(is.na(opt$use_alt_id)){
	if(is.na(opt$gene_id_map)) stop('When --use_alt_id is not set, --gene_id_map is required (live BioMart not supported in fast script)')
	Genes <- read.table(opt$gene_id_map, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
	Genes <- Genes[!is.na(Genes$entrezgene_id) & !is.na(Genes$external_gene_name), ]
	Genes <- Genes[!duplicated(Genes$entrezgene_id), ]
	Genes <- Genes[!duplicated(Genes$external_gene_name), ]
	TWAS  <- merge(TWAS, Genes, by.x = 'ID', by.y = 'external_gene_name')
}

# ---------------------------------------------------------------------------
# 2. Load gene sets / properties (mirrors TWAS-GSEA.V1.2.R:339-396).
# ---------------------------------------------------------------------------
if(!is.na(opt$gmt_file)){
	gene_sets <- read.gmt(opt$gmt_file)
	names(gene_sets) <- gsub('[[:punct:]]', '.', names(gene_sets))
	log_msg('Gene set file: ', length(gene_sets), ' sets.\n', sep = '')

	mem_mat_full <- foreach(i = seq_along(gene_sets), .combine = cbind) %dopar% {
		if(is.na(opt$use_alt_id)) TWAS$entrezgene_id %in% as.character(unlist(gene_sets[i]))
		else TWAS$Alt_ID %in% as.character(unlist(gene_sets[i]))
	}
	mem_mat_full <- as.matrix(mem_mat_full); colnames(mem_mat_full) <- names(gene_sets); storage.mode(mem_mat_full) <- 'double'
	rownames(mem_mat_full) <- TWAS$FILE
	TWAS_GS_full <- TWAS
	using_prop <- FALSE
} else {
	log_msg('Reading prop file... ')
	if(grepl('\\.rds$', opt$prop_file, ignore.case = TRUE)){
		prop_mat <- readRDS(opt$prop_file)
	} else {
		prop_dt <- fread(opt$prop_file)
		id_col  <- prop_dt[[1]]
		prop_mat <- as.matrix(prop_dt[, -1, with = FALSE])
		rownames(prop_mat) <- id_col
		rm(prop_dt); gc(verbose = FALSE)
	}
	log_msg('done.\n')
	n_props <- ncol(prop_mat)
	log_msg('Gene property file: ', n_props, ' properties.\n', sep = '')

	twas_ids  <- if(is.na(opt$use_alt_id)) TWAS$entrezgene_id else TWAS$Alt_ID
	keep_rows <- which(rownames(prop_mat) %in% twas_ids)
	if(length(keep_rows) == 0) stop('No overlap between TWAS gene IDs and --prop_file ID column.')
	log_msg('  ', length(keep_rows), ' / ', nrow(prop_mat), ' prop rows overlap TWAS.\n', sep = '')
	prop_mat <- prop_mat[keep_rows, , drop = FALSE]

	prop_mat[!is.finite(prop_mat)] <- 0

	# Join prop_mat rows to TWAS order; defer min_Ngenes filter and z-scoring
	# to inside the threshold loop.
	join_idx  <- match(twas_ids, rownames(prop_mat))
	keep_twas <- !is.na(join_idx)
	TWAS_GS_full <- TWAS[keep_twas, ]
	prop_mat_raw <- prop_mat[join_idx[keep_twas], , drop = FALSE]
	rownames(prop_mat_raw) <- TWAS_GS_full$FILE

	using_prop <- TRUE
}
if(!using_prop){
	log_msg(ncol(mem_mat_full), ' gene sets loaded.\n', sep = '')
}

if(opt$covar[1] != 'none'){
	missing_covars <- opt$covar[!(opt$covar %in% names(TWAS_GS_full))]
	if(length(missing_covars) > 0) stop('--covar not in TWAS data: ', paste(missing_covars, collapse = ', '))
}

# ---------------------------------------------------------------------------
# 3. Load precomputed cor matrix and align rows.
# ---------------------------------------------------------------------------
log_msg('Loading precomputed cor matrix... ')
cor_obj <- readRDS(opt$input_CorMat)
# build_cor_matrix.R now saves list(K, blocks); older runs saved a bare matrix.
if(is.list(cor_obj) && !is.null(cor_obj$K)){
	cor_block_all    <- cor_obj$K
	block_index_full <- cor_obj$blocks
} else {
	cor_block_all    <- cor_obj
	block_index_full <- NULL
}
log_msg('Done (', dim(cor_block_all)[1], ' x ', dim(cor_block_all)[2], ').\n', sep = '')

# ---------------------------------------------------------------------------
# 3b. Loop over TWAS p-value thresholds.
# ---------------------------------------------------------------------------
for(pT in opt$twas_p_thresh){

log_msg('\n--- TWAS p-value threshold = ', pT, ' ---\n', sep = '')

# Subset genes by TWAS p-value.
if(pT < 1){
	TWAS_GS_t <- TWAS_GS_full[TWAS_GS_full$TWAS.P <= pT, ]
} else {
	TWAS_GS_t <- TWAS_GS_full
}
log_msg(nrow(TWAS_GS_t), ' genes with TWAS.P <= ', pT, '.\n', sep = '')
if(nrow(TWAS_GS_t) < opt$min_Ngenes){
	log_msg('  Fewer than --min_Ngenes (', opt$min_Ngenes, '); skipping threshold.\n', sep = '')
	next
}

# Re-apply min_Ngenes filter on the threshold-subset of genes.
if(!using_prop){
	mem_mat_t <- mem_mat_full[TWAS_GS_t$FILE, , drop = FALSE]
	keep_gs   <- colnames(mem_mat_t)[colSums(mem_mat_t) >= opt$min_Ngenes]
	if(length(keep_gs) == 0){
		log_msg('  No gene sets with >= ', opt$min_Ngenes, ' genes; skipping threshold.\n', sep = '')
		next
	}
	mem_mat_t <- mem_mat_t[, keep_gs, drop = FALSE]
	TWAS_GS_t <- cbind(TWAS_GS_t, as.data.frame(mem_mat_t))
	gene_sets_clean <- keep_gs
} else {
	prop_mat_t <- prop_mat_raw[TWAS_GS_t$FILE, , drop = FALSE]
	nz_per_prop <- colSums(prop_mat_t != 0)
	keep_gs     <- colnames(prop_mat_t)[nz_per_prop >= opt$min_Ngenes]
	if(length(keep_gs) == 0){
		log_msg('  No properties with >= ', opt$min_Ngenes, ' non-zero genes; skipping threshold.\n', sep = '')
		next
	}
	prop_mat_t <- prop_mat_t[, keep_gs, drop = FALSE]
	gene_sets_clean <- keep_gs
}
log_msg(length(gene_sets_clean), ' gene sets/properties retained after --min_Ngenes filter.\n', sep = '')

# Intersect with correlation matrix.
genes_overlap <- intersect(TWAS_GS_t$FILE, colnames(cor_block_all))
if(length(genes_overlap) < opt$min_Ngenes){
	log_msg('  Fewer than --min_Ngenes genes overlap cor matrix; skipping threshold.\n')
	next
}
TWAS_GS_t <- TWAS_GS_t[TWAS_GS_t$FILE %in% genes_overlap, ]
TWAS_GS_t <- TWAS_GS_t[order(match(TWAS_GS_t$FILE, colnames(cor_block_all))), ]
idx <- match(TWAS_GS_t$FILE, colnames(cor_block_all))
cor_block_t <- cor_block_all[idx, idx]
N_genes <- nrow(TWAS_GS_t)
log_msg(N_genes, ' genes used after intersecting TWAS / cor matrix / gene sets.\n', sep = '')

# Align and z-score prop_mat for this threshold.
if(using_prop){
	prop_mat_t <- prop_mat_t[TWAS_GS_t$FILE, , drop = FALSE]
	N_Mem_Avail_prop <- colSums(prop_mat_t != 0)
	col_means <- colMeans(prop_mat_t, na.rm = TRUE)
	col_sds   <- sqrt(colSums((prop_mat_t - rep(col_means, each = N_genes))^2, na.rm = TRUE) / max(N_genes - 1L, 1L))
	for(j in seq_len(ncol(prop_mat_t))){
		if(!is.na(col_sds[j]) && col_sds[j] > 0){
			prop_mat_t[, j] <- (prop_mat_t[, j] - col_means[j]) / col_sds[j]
		} else {
			prop_mat_t[, j] <- 0
		}
	}
	prop_mat_t[is.na(prop_mat_t)] <- 0
}

# Recover the per-row block index for the subsetted K.
if(!is.null(block_index_full)){
	block_index <- block_index_full[match(TWAS_GS_t$FILE, names(block_index_full))]
	block_index <- as.integer(factor(block_index))   # contiguous 1..B after subset
} else {
	K_pat <- as(cor_block_t != 0, 'lgCMatrix')
	parent <- seq_len(N_genes)
	find <- function(i){ while(parent[i] != i){ parent[i] <<- parent[parent[i]]; i <- parent[i] }; i }
	tri <- which(K_pat & upper.tri(K_pat), arr.ind = TRUE)
	for(k in seq_len(nrow(tri))){
		ra <- find(tri[k, 1]); rb <- find(tri[k, 2])
		if(ra != rb) parent[ra] <- rb
	}
	roots <- vapply(seq_len(N_genes), find, integer(1))
	block_index <- as.integer(factor(roots))
}
B <- max(block_index)
log_msg('Cor matrix decomposes into ', B, ' independent blocks.\n', sep = '')

K_blocks <- lapply(seq_len(B), function(b){
	idx_b <- which(block_index == b)
	as.matrix(cor_block_t[idx_b, idx_b, drop = FALSE])
})

# ---------------------------------------------------------------------------
# 4. REML fit of (sigma_u^2, sigma_e^2) on the null model.
# ---------------------------------------------------------------------------
y <- TWAS_GS_t$ZSCORE
if(opt$covar[1] != 'none'){
	X_null <- cbind(1, as.matrix(TWAS_GS_t[, opt$covar, drop = FALSE]))
} else {
	X_null <- matrix(1, nrow = N_genes, ncol = 1)
}

log_msg('Fitting variance components by block-diagonal REML... ')
reml_t0 <- Sys.time()
reml <- fit_reml_blockdiag(y = y, X = X_null,
                            K_blocks = K_blocks, block_index = block_index,
                            h_max = opt$h_max, tol = opt$reml_tol)
reml_secs <- as.numeric(difftime(Sys.time(), reml_t0, units = 'secs'))
log_msg('Done in ', round(reml_secs, 2), 's.\n', sep = '')
log_msg(sprintf('  sigma_u^2 = %.6f   sigma_e^2 = %.6f   h_hat = %.6f\n',
                reml$sigma_u2, reml$sigma_e2, reml$h_hat))

# ---------------------------------------------------------------------------
# 5. Build V_hat = sigma_u^2 K + sigma_e^2 I and sparse Cholesky once.
# ---------------------------------------------------------------------------
log_msg('Cholesky-decomposing V_hat... ')
V_hat  <- reml$sigma_u2 * cor_block_t + reml$sigma_e2 * Diagonal(N_genes)
chol_V <- Cholesky(forceSymmetric(V_hat), perm = FALSE, LDL = FALSE)
log_msg('Done.\n')

# ---------------------------------------------------------------------------
# 6. Whiten y and null design; QR-residualise.
# ---------------------------------------------------------------------------
y_wh      <- as.numeric(solve(chol_V, y, system = 'L'))
X_null_wh <- as.matrix(solve(chol_V, X_null, system = 'L'))
Q_null    <- qr.Q(qr(X_null_wh))
y_wh_r    <- y_wh - Q_null %*% crossprod(Q_null, y_wh)

# ---------------------------------------------------------------------------
# 7. Whiten + residualise all gene-set vectors in one batched solve.
# ---------------------------------------------------------------------------
log_msg('Running vectorised GLS over ', length(gene_sets_clean), ' gene sets/properties... ', sep = '')
Z_gs   <- if(using_prop) prop_mat_t else as.matrix(TWAS_GS_t[, gene_sets_clean, drop = FALSE])
storage.mode(Z_gs) <- 'double'
Z_wh   <- as.matrix(solve(chol_V, Z_gs, system = 'L'))
rm(Z_gs); gc(verbose = FALSE)
Z_wh   <- Z_wh - Q_null %*% crossprod(Q_null, Z_wh)

# Whitened residuals are unit-variance by construction (V_hat absorbs the full
# Var(y)), so SE = 1/sqrt(denom) — same as V1.2's --fast_competitive T path.
denom    <- colSums(Z_wh^2)
numer    <- drop(crossprod(Z_wh, y_wh_r))
beta_hat <- numer / denom
SE_hat   <- 1 / sqrt(denom)
t_stat   <- numer / sqrt(denom)
p_val    <- compute_pval(t_stat, opt$two_sided)
log_msg('Done.\n')

# ---------------------------------------------------------------------------
# 8. Assemble + write results (column layout matches TWAS-GSEA.V1.2.R).
# ---------------------------------------------------------------------------
if(using_prop){
	N_Mem_Avail <- N_Mem_Avail_prop[gene_sets_clean]
	Results <- data.frame(
		GeneSet     = gene_sets_clean,
		Estimate    = beta_hat,
		SE          = SE_hat,
		T           = t_stat,
		N_Mem_Avail = N_Mem_Avail,
		P           = p_val,
		stringsAsFactors = FALSE)
} else {
	N_Mem_Avail <- colSums(TWAS_GS_t[, gene_sets_clean, drop = FALSE] != 0)
	N_Mem       <- vapply(gene_sets_clean, function(gs) length(gene_sets[[gs]]), integer(1))
	Results <- data.frame(
		GeneSet     = gene_sets_clean,
		Estimate    = beta_hat,
		SE          = SE_hat,
		T           = t_stat,
		N_Mem_Avail = N_Mem_Avail,
		N_Mem       = N_Mem,
		P           = p_val,
		stringsAsFactors = FALSE)
}
Results <- Results[is.finite(Results$T), ]
Results$P.CORR <- p.adjust(Results$P, method = opt$p_cor_method)
Results <- Results[order(Results$P), ]

out_suffix <- if(pT == 1) '' else paste0('.pT', pT)
out_file   <- paste0(opt$output, out_suffix, '.competitive.txt')
write.table(Results, out_file, col.names = TRUE, row.names = FALSE, quote = FALSE)
log_msg('Wrote ', out_file, ' (', nrow(Results), ' rows).\n', sep = '')

}  # end for(pT ...)

end.time <- Sys.time()
log_msg('\nFinished at ', as.character(end.time), ' (elapsed ', round(as.numeric(difftime(end.time, start.time, units = 'secs')), 1), 's)\n', sep = '')
