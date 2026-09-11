# sim/R/in_session_core.R
#
# In-session port of TWAS-GSEA-fast.R's core pipeline (lines 174-419 of that
# file). Lets the runner call the same statistical core thousands of times
# from one R process, skipping R startup, .CorMat.RDS load, TWAS file parsing,
# and per-block eigen recomputation.
#
# Only the numerical pipeline is reproduced. Set/property construction,
# ID mapping, and gmt reading live in the caller (sim/R/feature_builders.R).
# This keeps score_dataset() a pure function of aligned matrices, which is
# what the parity check needs.
#
# `feature_mat` (rows = TWAS FILE ids, cols = gene sets or properties):
#   feature_type = 'prop': numeric matrix; will be z-scored over the retained
#                          universe with (N-1)-denominator, matching
#                          TWAS-GSEA-fast.R:326-340. N_Mem_Avail = colSums(mat != 0).
#   feature_type = 'set' : 0/1 membership matrix; NOT z-scored.
#                          N_Mem_Avail = colSums(mat) counting non-zero entries.
#
# `mode` (character): 'directional' | 'probit' | 'magnitude', matching
#   TWAS-GSEA-fast.R:174-185. `outlier_threshold` = c(lower, upper).

suppressMessages({
  library(Matrix)
  library(VGAM)
})

compute_pval <- function(t_val, two_sided){
  if(two_sided) 2 * pnorm(-abs(t_val)) else 1 - pnorm(t_val)
}

# Apply the same clip as TWAS-GSEA-fast.R:174-185. Returns a numeric vector.
zscore_outcome <- function(twas_z, twas_p, mode, outlier_threshold){
  z <- switch(mode,
    directional = twas_z,
    probit      = probitlink(1 - twas_p),
    magnitude   = abs(twas_z),
    stop('unknown mode: ', mode))
  lo <- outlier_threshold[1]; hi <- outlier_threshold[2]
  z[z > hi] <- hi
  if(mode != 'magnitude') z[z < lo] <- lo
  z
}

# score_dataset(prep, twas_df, feature_mat, feature_type, K,
#               mode = 'probit', outlier_threshold = c(-3, 6),
#               two_sided = FALSE, min_Ngenes = 2,
#               h_max = 100, reml_tol = 1e-6,
#               p_cor_method = 'fdr')
#
# `prep`         : list from panel_prep(), containing gene_universe, block_index
#                  (used ONLY to look up universe order + block partition; NOT to
#                  reconstruct K, since panel_prep symmetrises via (K+t(K))/2 and
#                  the CLI uses two different triangles at two different stages).
# `twas_df`      : data.frame with FILE, TWAS.Z, TWAS.P (any subset of universe).
# `feature_mat`  : numeric matrix; rownames must intersect twas_df$FILE.
# `K`            : dgCMatrix as read from .CorMat.RDS / .CorMatSigned.RDS,
#                  RAW (no symmetrisation). fit_reml_blockdiag's internal
#                  eigen(symmetric=TRUE) reads lower tri; Cholesky(forceSymmetric(...))
#                  reads upper tri — same as the CLI, so results match bit for bit.
#
# Returns data.frame with the same columns as TWAS-GSEA-fast.R's
# <output>.competitive.txt: GeneSet, Estimate, SE, T, N_Mem_Avail, [N_Mem], P, P.CORR.
score_dataset <- function(prep, twas_df, feature_mat, feature_type, K,
                          mode = 'probit',
                          outlier_threshold = c(-3, 6),
                          two_sided = FALSE,
                          min_Ngenes = 2,
                          h_max = 100, reml_tol = 1e-6,
                          p_cor_method = 'fdr') {
  stopifnot(feature_type %in% c('prop','set'))
  stopifnot(all(c('FILE','TWAS.Z','TWAS.P') %in% names(twas_df)))
  stopifnot(inherits(K, 'Matrix'))
  # 1. Compute outcome.
  twas_df$ZSCORE <- zscore_outcome(twas_df$TWAS.Z, twas_df$TWAS.P, mode, outlier_threshold)

  # 2. Restrict to universe.
  universe <- intersect(twas_df$FILE, prep$gene_universe)
  universe <- intersect(universe, rownames(feature_mat))
  if(length(universe) < min_Ngenes) stop('universe smaller than min_Ngenes')
  # Preserve prep$gene_universe order (matches K row order).
  universe <- prep$gene_universe[prep$gene_universe %in% universe]

  # 3. Align twas_df to universe order.
  twas_df <- twas_df[match(universe, twas_df$FILE), , drop = FALSE]

  # 4. Align feature_mat and apply min_Ngenes filter.
  fm <- feature_mat[universe, , drop = FALSE]
  N_genes <- length(universe)
  N_Mem_Avail_all <- colSums(fm != 0)
  keep_gs <- names(N_Mem_Avail_all)[N_Mem_Avail_all >= min_Ngenes]
  if(length(keep_gs) == 0) stop('no sets/properties meet min_Ngenes')
  fm <- fm[, keep_gs, drop = FALSE]
  N_Mem_Avail <- N_Mem_Avail_all[keep_gs]

  # 5. If property: z-score over the retained universe with (N-1) denom,
  # matching TWAS-GSEA-fast.R:326-340.
  if(feature_type == 'prop'){
    col_means <- colMeans(fm, na.rm = TRUE)
    col_sds   <- sqrt(colSums((fm - rep(col_means, each = N_genes))^2, na.rm = TRUE) / max(N_genes - 1L, 1L))
    for(j in seq_len(ncol(fm))){
      if(!is.na(col_sds[j]) && col_sds[j] > 0){
        fm[, j] <- (fm[, j] - col_means[j]) / col_sds[j]
      } else {
        fm[, j] <- 0
      }
    }
    fm[is.na(fm)] <- 0
  }

  # 6. Assemble K subset — RAW, no symmetrisation. See the CLI mirror in
  # TWAS-GSEA-fast.R:322-364: K is passed asymmetric to both fit_reml_blockdiag
  # (which reads lower tri via eigen(symmetric=TRUE)) and to Cholesky+
  # forceSymmetric (which reads upper tri). Matching this exactly is what
  # keeps the parity gate satisfied.
  stopifnot(all(universe %in% rownames(K)))
  cor_block_t <- K[universe, universe]

  # Recover per-row block index for the retained universe (mirrors
  # TWAS-GSEA-fast.R:343-357 when block_index_full is supplied).
  block_index_universe <- prep$block_index[universe]
  block_index <- as.integer(factor(block_index_universe))
  B <- max(block_index)
  K_blocks <- lapply(seq_len(B), function(b){
    idx_b <- which(block_index == b)
    as.matrix(cor_block_t[idx_b, idx_b, drop = FALSE])
  })

  # 7. Fit REML.
  y      <- twas_df$ZSCORE
  X_null <- matrix(1, nrow = N_genes, ncol = 1L)
  reml <- fit_reml_blockdiag(y = y, X = X_null,
                              K_blocks = K_blocks, block_index = block_index,
                              h_max = h_max, tol = reml_tol)

  # 8. V_hat, Cholesky, whiten.
  V_hat  <- reml$sigma_u2 * cor_block_t + reml$sigma_e2 * Diagonal(N_genes)
  chol_V <- Cholesky(forceSymmetric(V_hat), perm = FALSE, LDL = FALSE)

  y_wh      <- as.numeric(solve(chol_V, y, system = 'L'))
  X_null_wh <- as.matrix(solve(chol_V, X_null, system = 'L'))
  Q_null    <- qr.Q(qr(X_null_wh))
  y_wh_r    <- y_wh - Q_null %*% crossprod(Q_null, y_wh)

  # 9. Whiten feature_mat in one batched solve, then residualise.
  Z_gs <- fm; storage.mode(Z_gs) <- 'double'
  Z_wh <- as.matrix(solve(chol_V, Z_gs, system = 'L'))
  Z_wh <- Z_wh - Q_null %*% crossprod(Q_null, Z_wh)

  denom    <- colSums(Z_wh^2)
  numer    <- drop(crossprod(Z_wh, y_wh_r))
  beta_hat <- numer / denom
  SE_hat   <- 1 / sqrt(denom)
  t_stat   <- numer / sqrt(denom)
  p_val    <- compute_pval(t_stat, two_sided)

  # 10. Assemble output. For gene sets we also emit N_Mem (set size in the
  # original gmt), which the CLI has via `length(gene_sets[[gs]])`. In the
  # in-session path the caller supplies mem_mat already, so N_Mem is only
  # available if we passed it through — expose it via an attribute for the
  # scan test.
  Results <- data.frame(
    GeneSet     = colnames(fm),
    Estimate    = as.numeric(beta_hat),
    SE          = as.numeric(SE_hat),
    T           = as.numeric(t_stat),
    N_Mem_Avail = as.integer(N_Mem_Avail),
    P           = as.numeric(p_val),
    stringsAsFactors = FALSE)
  Results <- Results[is.finite(Results$T), ]
  Results$P.CORR <- p.adjust(Results$P, method = p_cor_method)
  Results <- Results[order(Results$P), ]
  attr(Results, 'reml') <- reml
  attr(Results, 'N_genes') <- N_genes
  Results
}
