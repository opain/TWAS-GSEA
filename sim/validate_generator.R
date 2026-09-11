#!/usr/bin/Rscript
# sim/validate_generator.R
#
# Phase 2 gate. Produces:
#   (1) Null marginal QQ + KS: TWAS.Z from mu=0 simulations should be N(0,1).
#   (2) Block correlation recovery: pooled empirical cor(z) across many null
#       reps should match K_signed on a handful of representative blocks.
#   (3) Injection sanity: with a strong signed signature t, the ordinary
#       correlation cor(z, t) should recover beta positively.
# Writes a multi-page PDF + a text log.

suppressMessages({
  library(optparse)
  library(data.table)
  library(Matrix)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option('--panel_prep', default = 'sim/output/panel_cache/Whole_Blood.prep.RDS', type = 'character'),
  make_option('--K_signed',   default = 'sim/output/panel_cache/Whole_Blood.CorMatSigned.RDS', type = 'character'),
  make_option('--n_null',     default = 500L, type = 'integer',
              help = 'Null reps for QQ + block-corr recovery [default 500]'),
  make_option('--n_signal',   default = 50L,  type = 'integer',
              help = 'Reps for injection sanity check [default 50]'),
  make_option('--outdir',     default = 'sim/output/validate', type = 'character')
)))

SCRIPT_DIR <- dirname(sub('--file=', '', grep('--file=', commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
source(file.path(SCRIPT_DIR, 'R', 'sim_generator.R'))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
LOG <- file.path(opt$outdir, 'generator.log')
PDF <- file.path(opt$outdir, 'generator.pdf')
log_msg <- function(...) {
  cat(..., '\n', sep = '')
  cat(..., '\n', sep = '', file = LOG, append = TRUE)
}
cat('', file = LOG)

prep    <- panel_prep(opt$panel_prep, cache_path = opt$panel_prep)
K_obj   <- readRDS(opt$K_signed)
K       <- K_obj$K
stopifnot(identical(rownames(K), prep$gene_universe))
N       <- length(prep$gene_universe)
blocks  <- prep$block_index

log_msg('validate_generator.R')
log_msg('universe = ', N, ' genes over ', max(blocks), ' blocks')
log_msg('n_null = ', opt$n_null, '   n_signal = ', opt$n_signal)

# -------------------------------------------------------------------------
# Draw the null replicates (in memory: N x n_null; ~7813 x 500 = 30 MB double)
# -------------------------------------------------------------------------
set.seed(1L)
t0 <- Sys.time()
Z_null <- matrix(0, nrow = N, ncol = opt$n_null)
for(r in seq_len(opt$n_null)) Z_null[, r] <- simulate_z(prep, mu = NULL)
log_msg('null sim: ', opt$n_null, ' reps in ',
        round(as.numeric(difftime(Sys.time(), t0, units = 'secs')), 2), ' s')

pdf(PDF, width = 8, height = 6)
op <- par(no.readonly = TRUE)

# --- (1) Marginal null: pool all Z into a single sample
z_all <- as.numeric(Z_null)
ks <- suppressWarnings(ks.test(sample(z_all, min(50000, length(z_all))), 'pnorm'))
log_msg(sprintf('marginal null: pooled mean=%.4f  sd=%.4f  KS p=%.3g',
                mean(z_all), sd(z_all), ks$p.value))
qq_lim <- range(qnorm(ppoints(min(20000, length(z_all)))))
qq_sample <- sample(z_all, min(20000, length(z_all)))
qqnorm(qq_sample, main = sprintf('Null pooled Z (subsample n=%d, KS p=%.3g)',
                                 length(qq_sample), ks$p.value),
       pch = '.', xlim = qq_lim, ylim = qq_lim)
abline(0, 1, col = 'red')

# --- (2) Block correlation recovery
# Pick a small representative block, a medium one, and the largest.
block_sizes <- as.numeric(table(blocks))
targets <- c(
  small_idx  = which(block_sizes >= 5 & block_sizes <= 20)[1],
  medium_idx = which.min(abs(block_sizes - median(block_sizes[block_sizes > 1]))),
  large_idx  = which.max(block_sizes)
)
log_msg('block-corr targets: sizes = ',
        paste(block_sizes[targets], collapse = ', '))

# Under the null MVN, each empirical off-diagonal cor is approximately
# N(rho_true, (1 - rho_true^2)^2 / (n_null - 1)). We z-score the residuals
# and check the fraction within +/- 2 (expected ~95%) and RMSE (expected ~SE).
# max|diff| is misleading because it scales with the number of off-diagonal
# entries (i.e. block size squared).
pass_flags <- logical(length(targets))
for(k in seq_along(targets)){
  b   <- targets[k]
  idx <- which(blocks == b)
  n_b <- length(idx)
  if(n_b < 2) next
  Kb_true <- as.matrix(K[idx, idx])
  Kb_true <- (Kb_true + t(Kb_true)) / 2   # see note in earlier commit
  Zb <- Z_null[idx, , drop = FALSE]
  Kb_emp  <- cor(t(Zb))
  # Off-diagonal, upper triangle only (each pair counted once)
  off <- upper.tri(Kb_true, diag = FALSE)
  d <- Kb_emp[off] - Kb_true[off]
  se <- (1 - Kb_true[off]^2) / sqrt(opt$n_null - 1)
  z  <- d / pmax(se, 1e-8)
  frac_within_2 <- mean(abs(z) <= 2)
  frac_within_3 <- mean(abs(z) <= 3)
  rmse <- sqrt(mean(d^2))
  mean_se <- mean(se)
  pass_flags[k] <- (frac_within_2 >= 0.90) && (frac_within_3 >= 0.985)
  log_msg(sprintf('  block %d (size %d, %d pairs): RMSE=%.4f  mean(SE)=%.4f  |z|<=2: %.1f%%  |z|<=3: %.2f%%  PASS=%s',
                  b, n_b, sum(off), rmse, mean_se, 100*frac_within_2, 100*frac_within_3, pass_flags[k]))
  plot(Kb_true[off], Kb_emp[off], pch = '.', cex = 1.4,
       xlab = 'K_signed (true)', ylab = 'cor(Z_null) (empirical)',
       main = sprintf('Block %d  n=%d  RMSE=%.3f  |z|<=2: %.1f%%',
                      b, n_b, rmse, 100*frac_within_2))
  abline(0, 1, col = 'red')
}

# --- (3) Injection sanity: strong beta, random signed signature
set.seed(7L)
sparsity <- 0.05
n_nz  <- max(1L, round(sparsity * N))
sig_idx <- sample.int(N, n_nz)
t_vec <- numeric(N); t_vec[sig_idx] <- rnorm(n_nz)
t_vec <- (t_vec - mean(t_vec)) / sd(t_vec)
names(t_vec) <- prep$gene_universe

betas <- c(0, 0.2, 0.5, 1.0)
recov <- data.frame(beta = betas, mean_cor = NA_real_, se_cor = NA_real_)
for(bi in seq_along(betas)){
  bt <- betas[bi]
  cors <- numeric(opt$n_signal)
  for(r in seq_len(opt$n_signal)){
    mu <- mu_property(prep, t_vec, bt)
    z  <- simulate_z(prep, mu = mu)
    cors[r] <- cor(z, t_vec)
  }
  recov$mean_cor[bi] <- mean(cors)
  recov$se_cor[bi]   <- sd(cors) / sqrt(length(cors))
  log_msg(sprintf('injection beta=%.2f: mean cor(z, t) = %.3f  (SE %.3f)',
                  bt, recov$mean_cor[bi], recov$se_cor[bi]))
}

par(op)
plot(recov$beta, recov$mean_cor,
     xlab = 'beta (injection scale)', ylab = 'mean cor(z, t) over reps',
     main = 'Injection sanity',
     pch = 19, ylim = c(-0.05, max(recov$mean_cor + 3*recov$se_cor)))
arrows(recov$beta, recov$mean_cor - 2*recov$se_cor,
       recov$beta, recov$mean_cor + 2*recov$se_cor, angle = 90, length = 0.05, code = 3)
abline(h = 0, lty = 2)

dev.off()
log_msg('wrote ', PDF)
log_msg('wrote ', LOG)

# ---------- Gate: pass / fail ----------
pass_marginal <- ks$p.value > 0.001
pass_blockcor <- all(pass_flags, na.rm = TRUE)
pass_injection <- with(recov, mean_cor[betas == 0] < 0.05 & mean_cor[betas == max(betas)] > 0.2)

log_msg('GATE marginal-null KS p > 1e-3:  ', pass_marginal)
log_msg('GATE block-corr z-scored (|z|<=2 >=90% AND |z|<=3 >=98.5%) across chosen blocks:  ', pass_blockcor)
log_msg('GATE injection: cor(z,t) at beta=0 <0.05 AND at beta=', max(betas), ' >0.2:  ', pass_injection)
log_msg('OVERALL: ', if(pass_marginal && pass_blockcor && pass_injection) 'PASS' else 'FAIL')
