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

	mem_mat <- foreach(i = seq_along(gene_sets), .combine = cbind) %dopar% {
		if(is.na(opt$use_alt_id)) TWAS$entrezgene_id %in% as.character(unlist(gene_sets[i]))
		else TWAS$Alt_ID %in% as.character(unlist(gene_sets[i]))
	}
	mem_mat <- as.matrix(mem_mat); colnames(mem_mat) <- names(gene_sets); storage.mode(mem_mat) <- 'double'
	keep_gs <- colnames(mem_mat)[colSums(mem_mat) >= opt$min_Ngenes]
	mem_mat <- mem_mat[, keep_gs, drop = FALSE]
	TWAS_GS <- cbind(TWAS, as.data.frame(mem_mat))
	gene_sets_clean <- keep_gs
	using_prop <- FALSE
} else {
	gene_prop <- data.frame(fread(opt$prop_file))
	log_msg('Gene property file: ', ncol(gene_prop) - 1, ' properties.\n', sep = '')
	if(is.na(opt$use_alt_id)){
		TWAS_GS <- merge(TWAS, gene_prop, by.x = 'entrezgene_id', by.y = 'ID')
	} else {
		TWAS_GS <- merge(TWAS, gene_prop, by.x = 'Alt_ID', by.y = 'ID')
	}
	prop_cols <- names(gene_prop)[-1]
	prop_cols <- prop_cols[prop_cols %in% names(TWAS_GS)]
	keep_gs <- prop_cols[colSums(abs(TWAS_GS[, prop_cols, drop = FALSE])) >= opt$min_Ngenes]
	for(i in keep_gs) TWAS_GS[[i]] <- as.numeric(scale(TWAS_GS[[i]]))
	gene_sets_clean <- keep_gs
	using_prop <- TRUE
}
log_msg(length(gene_sets_clean), ' gene sets/properties retained after --min_Ngenes filter.\n', sep = '')

if(opt$covar[1] != 'none'){
	missing_covars <- opt$covar[!(opt$covar %in% names(TWAS_GS))]
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

genes_overlap <- intersect(TWAS_GS$FILE, colnames(cor_block_all))
if(length(genes_overlap) == 0) stop('No overlap between TWAS FILE column and cor matrix rows.')
TWAS_GS <- TWAS_GS[TWAS_GS$FILE %in% genes_overlap, ]
TWAS_GS <- TWAS_GS[order(match(TWAS_GS$FILE, colnames(cor_block_all))), ]
idx <- match(TWAS_GS$FILE, colnames(cor_block_all))
cor_block_all <- cor_block_all[idx, idx]
N_genes <- nrow(TWAS_GS)
log_msg(N_genes, ' genes used after intersecting TWAS / cor matrix / gene sets.\n', sep = '')

# Recover the per-row block index for the subsetted K. If the precomputed file
# carries one we use it (renumbered after subset); otherwise derive it from the
# sparsity pattern via connected components.
if(!is.null(block_index_full)){
	block_index <- block_index_full[match(TWAS_GS$FILE, names(block_index_full))]
	block_index <- as.integer(factor(block_index))   # contiguous 1..B after subset
} else {
	# Derive blocks as connected components of K's non-zero pattern.
	K_pat <- as(cor_block_all != 0, 'lgCMatrix')
	parent <- seq_len(N_genes)
	find <- function(i){ while(parent[i] != i){ parent[i] <<- parent[parent[i]]; i <- parent[i] }; i }
	tri <- which(K_pat & upper.tri(K_pat), arr.ind = TRUE)
	for(k in seq_len(nrow(tri))){
		ra <- find(tri[k, 1]); rb <- find(tri[k, 2])
		if(ra != rb) parent[ra] <<- rb
	}
	roots <- vapply(seq_len(N_genes), find, integer(1))
	block_index <- as.integer(factor(roots))
}
B <- max(block_index)
log_msg('Cor matrix decomposes into ', B, ' independent blocks.\n', sep = '')

# Per-block dense submatrices for the REML fitter.
K_blocks <- lapply(seq_len(B), function(b){
	idx_b <- which(block_index == b)
	as.matrix(cor_block_all[idx_b, idx_b, drop = FALSE])
})

# ---------------------------------------------------------------------------
# 4. REML fit of (sigma_u^2, sigma_e^2) on the null model.
# ---------------------------------------------------------------------------
y <- TWAS_GS$ZSCORE
if(opt$covar[1] != 'none'){
	X_null <- cbind(1, as.matrix(TWAS_GS[, opt$covar, drop = FALSE]))
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
V_hat  <- reml$sigma_u2 * cor_block_all + reml$sigma_e2 * Diagonal(N_genes)
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
Z_gs   <- as.matrix(TWAS_GS[, gene_sets_clean, drop = FALSE])
storage.mode(Z_gs) <- 'double'
Z_wh   <- as.matrix(solve(chol_V, Z_gs, system = 'L'))
Z_wh_r <- Z_wh - Q_null %*% crossprod(Q_null, Z_wh)

# Whitened residuals are unit-variance by construction (V_hat absorbs the full
# Var(y)), so SE = 1/sqrt(denom) — same as V1.2's --fast_competitive T path.
denom    <- colSums(Z_wh_r^2)
numer    <- drop(crossprod(Z_wh_r, y_wh_r))
beta_hat <- numer / denom
SE_hat   <- 1 / sqrt(denom)
t_stat   <- numer / sqrt(denom)
p_val    <- compute_pval(t_stat, opt$two_sided)
log_msg('Done.\n')

# ---------------------------------------------------------------------------
# 7. Assemble + write results (column layout matches TWAS-GSEA.V1.2.R).
# ---------------------------------------------------------------------------
if(using_prop){
	Results <- data.frame(
		GeneSet  = gene_sets_clean,
		Estimate = beta_hat,
		SE       = SE_hat,
		T        = t_stat,
		P        = p_val,
		stringsAsFactors = FALSE)
} else {
	N_Mem_Avail <- colSums(TWAS_GS[, gene_sets_clean, drop = FALSE] != 0)
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
write.table(Results, paste0(opt$output, '.competitive.txt'),
            col.names = TRUE, row.names = FALSE, quote = FALSE)

end.time <- Sys.time()
log_msg('Wrote ', opt$output, '.competitive.txt (', nrow(Results), ' rows).\n', sep = '')
log_msg('Finished at ', as.character(end.time), ' (elapsed ', round(as.numeric(difftime(end.time, start.time, units = 'secs')), 1), 's)\n', sep = '')
