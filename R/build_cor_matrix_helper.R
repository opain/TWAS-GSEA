# build_cor_matrix_helper.R
#
# Construct a block-diagonal sparse gene-gene correlation matrix from a
# predicted-expression reference panel. Logic factored verbatim from
# TWAS-GSEA.V1.2.R lines 484-602 so the two callers (TWAS-GSEA.V1.2.R and
# build_cor_matrix.R) produce identical matrices.
#
# Required packages on the caller side: data.table, Matrix, WGCNA, foreach,
# matrixcalc (nearPD).

build_cor_matrix <- function(genes_df, GeneX_all, cor_window, min_r2, max_r2,
                             log_msg = function(...) cat(...),
                             log_progress = function(i, n) invisible(NULL)) {
	# genes_df: data.frame with columns FILE, CHR, P0, P1 (one row per gene),
	#           ordered as desired (caller is responsible for ordering by CHR/P0/P1).
	# GeneX_all: data.frame of expression values, one column per gene, column
	#            names matching genes_df$FILE in the same order.

	stopifnot(all(c('FILE','CHR','P0','P1') %in% names(genes_df)))
	stopifnot(nrow(genes_df) == ncol(GeneX_all))
	stopifnot(all(genes_df$FILE == names(GeneX_all)))

	# Determine gene blocks (vectorised cumsum, same logic as V1.2).
	n_genes <- nrow(genes_df)
	same_chr <- genes_df$CHR[-1] == genes_df$CHR[-n_genes]
	overlaps  <- genes_df$P1[-1] > (genes_df$P0[-n_genes] - cor_window) &
	             genes_df$P0[-1] < (genes_df$P1[-n_genes] + cor_window)
	genes_df$Block <- cumsum(c(TRUE, !(same_chr & overlaps)))

	log_msg('The genes could be separated into', length(unique(genes_df$Block)), 'blocks.\n')

	log_msg('Creating correlation matrix... ')
	cor_blocks_list <- foreach(i = unique(genes_df$Block)) %dopar% {
		frob_diff <- 0
		if(sum(genes_df$Block == i) == 1){
			cor_block_2 <- Matrix(1, nrow = 1, ncol = 1, sparse = TRUE)
			colnames(cor_block_2) <- genes_df$FILE[genes_df$Block == i]
			rownames(cor_block_2) <- genes_df$FILE[genes_df$Block == i]
		} else {
			cor_block <- abs(WGCNA::cor(as.matrix(GeneX_all[(names(GeneX_all) %in% genes_df$FILE[genes_df$Block == i])]), method = 'pearson'))
			tmp <- cor_block
			tmp[!lower.tri(tmp)] <- 0
			keep <- colnames(cor_block)[!apply(tmp, 2, function(x) any(abs(x) > sqrt(max_r2)))]
			cor_block_2 <- cor_block[(colnames(cor_block) %in% keep), (colnames(cor_block) %in% keep)]
			genes_df_Block <- genes_df[which(genes_df$Block == i), ]
			genes_df_Block <- genes_df_Block[(genes_df_Block$FILE %in% keep), ]
			if(length(cor_block_2) == 1){
				cor_block_2 <- Matrix(1, nrow = 1, ncol = 1, sparse = TRUE)
				colnames(cor_block_2) <- genes_df_Block$FILE
				rownames(cor_block_2) <- genes_df_Block$FILE
			} else {
				p0_b <- genes_df_Block$P0
				p1_b <- genes_df_Block$P1
				sparse_struc <- Matrix(
				    outer(p1_b, p0_b - cor_window, `>`) & outer(p0_b, p1_b + cor_window, `<`),
				    sparse = TRUE
				)
				cor_block_2[(sparse_struc[, ] != 1)@x] <- 0
				cor_block_2[abs(cor_block_2) < sqrt(min_r2)] <- 0
				pd_ok <- tryCatch({ chol(as.matrix(cor_block_2)); TRUE }, error = function(e) FALSE)
				if(!pd_ok){
					cor_block_before <- as.matrix(cor_block_2)
					cor_block_2 <- nearPD(cor_block_2, corr = TRUE)$mat
					frob_diff <- norm(as.matrix(cor_block_2) - cor_block_before, "F")
				}
				cor_block_2 <- Matrix(cor_block_2, sparse = TRUE)
			}
		}
		log_progress(i, length(unique(genes_df$Block)))
		list(mat = cor_block_2, frob_diff = frob_diff)
	}

	frob_diffs <- sapply(cor_blocks_list, `[[`, 'frob_diff')
	cor_blocks_list <- lapply(cor_blocks_list, `[[`, 'mat')
	cor_block_all <- bdiag(cor_blocks_list)
	all_block_names <- unlist(lapply(cor_blocks_list, rownames))
	rownames(cor_block_all) <- colnames(cor_block_all) <- all_block_names

	log_msg('Done!\n')
	n_repaired <- sum(frob_diffs > 0)
	if(n_repaired > 0){
		log_msg('WARNING:', n_repaired, 'of', length(frob_diffs), 'genomic blocks required positive-definiteness repair (nearPD). Max Frobenius norm distortion:', round(max(frob_diffs), 4), '.\n')
	}

	prop_sparse <- sum(cor_block_all == 0) / (dim(cor_block_all)[1] * dim(cor_block_all)[2])
	log_msg('The correlation matrix of gene expression is ', prop_sparse * 100, '% sparse.\n', sep = '')
	log_msg('After pruning ', dim(cor_block_all)[1], ' features remain.\n', sep = '')

	# Per-row block index aligned with rownames(cor_block_all). Used downstream
	# by the per-block REML fitter (see R/reml_blockdiag.R) so block structure
	# does not have to be re-derived from the sparsity pattern.
	block_sizes <- vapply(cor_blocks_list, nrow, integer(1))
	block_index <- rep(seq_along(cor_blocks_list), times = block_sizes)
	names(block_index) <- all_block_names

	list(K = cor_block_all, blocks = block_index)
}
