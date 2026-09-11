# sim/R/sim_generator.R
#
# Generative model for TWAS-GSEA-fast Tier-A simulations. Draws latent signed-z
# from a block-diagonal MVN(mu, Sigma_signed) using the SIGNED correlation
# matrix produced by sim/build_signed_cor.R. Provides injection modes (null,
# set-based, property-based) that shift mu. Emits a TWAS results data.frame
# / TSV compatible with TWAS-GSEA-fast.R.
#
# All routines are pure R (I/O only in the CLI wrapper). Panel preparation is
# cached once per (panel, seed = NA) and shared by all downstream sim runs.

suppressMessages({
  library(data.table)
  library(Matrix)
})

# -------------------------------------------------------------------------
# panel_prep(K_signed_path, cache_path = NULL)
#
# Loads the signed .CorMatSigned.RDS produced by sim/build_signed_cor.R and
# precomputes per-block Cholesky factors (for MVN sampling) and full
# eigendecompositions (for the in-session GLS core in Phase 3). Both derive
# from a nearPD-safe copy of the block, so downstream users never need to
# repeat the safety projection.
#
# Returns list with:
#   gene_universe : character vector, rownames of K in the same order as blocks
#   block_index   : named integer vector, gene -> block id (1..B)
#   block_names   : character vector of length B, each entry is FILE ids in block
#   chol_list     : list of B upper-triangular chol factors (matrix, dense)
#   eigen_list    : list of B eigen decompositions ($vectors, $values with d>=0)
#   frob_repair   : numeric vector length B, Frobenius distortion from nearPD (0 if none)
# -------------------------------------------------------------------------
panel_prep <- function(K_signed_path, cache_path = NULL,
                       log_msg = function(...) cat(...)) {
  if(!is.null(cache_path) && file.exists(cache_path)){
    log_msg('panel_prep: cache hit at ', cache_path, '\n', sep = '')
    return(readRDS(cache_path))
  }

  log_msg('panel_prep: loading ', K_signed_path, '\n', sep = '')
  obj <- readRDS(K_signed_path)
  stopifnot(is.list(obj) && !is.null(obj$K) && !is.null(obj$blocks))
  K       <- obj$K
  blocks  <- obj$blocks
  stopifnot(identical(names(blocks), rownames(K)))

  B <- max(blocks)
  gene_universe <- rownames(K)

  chol_list  <- vector('list', B)
  eigen_list <- vector('list', B)
  block_names <- vector('list', B)
  frob_repair <- numeric(B)

  for(b in seq_len(B)){
    idx <- which(blocks == b)
    block_names[[b]] <- gene_universe[idx]
    Kb <- as.matrix(K[idx, idx, drop = FALSE])

    # Ensure symmetric to machine precision.
    Kb <- (Kb + t(Kb)) / 2

    # Try chol; if not PD, run nearPD once more (belt & braces; the build step
    # already applied per-block nearPD but we saw non-trivial repair for some
    # blocks and want to guarantee samplers succeed).
    ch <- tryCatch(chol(Kb), error = function(e) NULL)
    if(is.null(ch)){
      before <- Kb
      Kb <- as.matrix(Matrix::nearPD(Kb, corr = TRUE)$mat)
      Kb <- (Kb + t(Kb)) / 2
      frob_repair[b] <- norm(Kb - before, 'F')
      ch <- chol(Kb)
    }
    chol_list[[b]] <- ch

    if(nrow(Kb) == 1){
      eigen_list[[b]] <- list(vectors = matrix(1, 1, 1), values = as.numeric(Kb[1,1]))
    } else {
      eg <- eigen(Kb, symmetric = TRUE)
      eg$values[eg$values < 0] <- 0     # clamp fp noise
      eigen_list[[b]] <- list(vectors = eg$vectors, values = eg$values)
    }
  }

  prep <- list(
    gene_universe = gene_universe,
    block_index   = blocks,
    block_names   = block_names,
    chol_list     = chol_list,
    eigen_list    = eigen_list,
    frob_repair   = frob_repair,
    K_source      = normalizePath(K_signed_path, mustWork = TRUE))

  n_rep <- sum(frob_repair > 0)
  log_msg('panel_prep: ', length(gene_universe), ' genes, ', B, ' blocks; ',
          n_rep, ' blocks needed a second nearPD pass',
          if(n_rep > 0) sprintf(' (max Frobenius %.3g)', max(frob_repair)) else '',
          '.\n', sep = '')

  if(!is.null(cache_path)){
    dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(prep, cache_path)
    log_msg('panel_prep: wrote cache ', cache_path, '\n', sep = '')
  }
  prep
}

# -------------------------------------------------------------------------
# simulate_z(prep, mu = NULL, seed = NULL)
#
# Draws z ~ MVN(mu, block-diag(K)) where each block is sampled independently
# using the cached Cholesky factor:
#     z_b = mu_b + t(L_b) %*% rnorm(n_b),   L_b upper-tri chol of K_b.
#
# If mu is NULL, draws N(0, K) directly.
# Returns a numeric vector aligned with prep$gene_universe.
# -------------------------------------------------------------------------
simulate_z <- function(prep, mu = NULL, seed = NULL) {
  if(!is.null(seed)) set.seed(seed)
  N <- length(prep$gene_universe)
  if(is.null(mu)) mu <- numeric(N)
  stopifnot(length(mu) == N)

  z <- numeric(N)
  for(b in seq_along(prep$chol_list)){
    idx <- which(prep$block_index == b)
    n_b <- length(idx)
    L_b <- prep$chol_list[[b]]         # upper tri, t(L)%*%L = K_b
    e   <- rnorm(n_b)
    z[idx] <- mu[idx] + as.numeric(crossprod(L_b, e))
  }
  z
}

# -------------------------------------------------------------------------
# Injection helpers. Each returns a numeric vector length = |universe|,
# aligned with prep$gene_universe.
# -------------------------------------------------------------------------
mu_null <- function(prep) numeric(length(prep$gene_universe))

# mu_set: pick rho * length(target_members) members as causal, assign
# s_g in {+1,-1} with prob pi_pos, mu_g = Delta * s_g. Non-causal members and
# non-members are zero.
#   target_members: character vector of FILE ids in prep$gene_universe.
#   Delta:  effect magnitude on the latent Z scale.
#   rho:    fraction of target members that carry signal (0..1).
#   pi_pos: fraction of causal genes with positive sign.
mu_set <- function(prep, target_members, Delta, rho = 1, pi_pos = 0.5) {
  mu <- numeric(length(prep$gene_universe))
  members <- intersect(target_members, prep$gene_universe)
  if(length(members) == 0) return(mu)
  n_causal <- max(1L, round(rho * length(members)))
  causal <- sample(members, n_causal)
  signs  <- ifelse(runif(n_causal) < pi_pos, 1, -1)
  idx    <- match(causal, prep$gene_universe)
  mu[idx] <- Delta * signs
  mu
}

# mu_property: mu_g = beta * t_g. Caller passes an already-standardised t on
# prep$gene_universe (or a named vector we align).
mu_property <- function(prep, t_vec, beta) {
  N <- length(prep$gene_universe)
  if(is.null(names(t_vec))){
    stopifnot(length(t_vec) == N)
    t_aligned <- as.numeric(t_vec)
  } else {
    t_aligned <- numeric(N)
    hit <- match(prep$gene_universe, names(t_vec))
    t_aligned[!is.na(hit)] <- t_vec[hit[!is.na(hit)]]
  }
  beta * t_aligned
}

# -------------------------------------------------------------------------
# emit_twas(z, template_twas, universe = NULL)
#
# Constructs a TWAS results data.frame usable as --twas_results input for
# TWAS-GSEA-fast.R. Copies FILE, ID, CHR, MODELCV.R2, P0, P1 (and any other
# columns present) from template rows matching `universe`, overwrites TWAS.Z
# with z and TWAS.P with 2*pnorm(-abs(z)). Preserves the row order of the
# universe.
# -------------------------------------------------------------------------
emit_twas <- function(z, template_twas, universe) {
  stopifnot('FILE' %in% names(template_twas))
  stopifnot(length(z) == length(universe))
  tpl <- as.data.table(template_twas)
  keep <- tpl[FILE %in% universe]
  keep <- keep[match(universe, FILE)]
  keep[, TWAS.Z := z]
  keep[, TWAS.P := 2 * pnorm(-abs(z))]
  as.data.frame(keep)
}

# -------------------------------------------------------------------------
# load_template_twas(pattern_or_dir, universe = NULL)
#
# Reads /data/twas_results which is a per-chromosome collection. Concatenates
# all 22+ files. Returns a data.frame with all columns and one row per gene.
# If universe is non-NULL, subsets to those FILE names.
# -------------------------------------------------------------------------
load_template_twas <- function(path, universe = NULL) {
  if(dir.exists(path)){
    files <- list.files(path, full.names = TRUE)
    dts <- lapply(files, function(f) tryCatch(fread(f), error = function(e) NULL))
    dts <- Filter(Negate(is.null), dts)
    tpl <- rbindlist(dts, use.names = TRUE, fill = TRUE)
  } else {
    tpl <- fread(path)
  }
  # Normalise FILE column the same way TWAS-GSEA-fast.R does.
  if('FILE' %in% names(tpl)){
    tpl[, FILE := sub('.wgt.RDat', '', sub('.*/', '', FILE))]
    tpl[, FILE := gsub(':', '.', FILE)]
    tpl[, FILE := gsub('-', '.', FILE)]
  }
  if(!is.null(universe)) tpl <- tpl[FILE %in% universe]
  as.data.frame(tpl)
}
