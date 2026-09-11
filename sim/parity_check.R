#!/usr/bin/Rscript
# sim/parity_check.R
#
# HARD GATE for Phase 3. Generates a handful of simulated TWAS datasets,
# runs both:
#   (a) TWAS-GSEA-fast.R via Rscript (subprocess);
#   (b) sim/R/in_session_core.R::score_dataset() in the current R process;
# and asserts entrywise agreement of P, Estimate, and N_Mem_Avail per set /
# property. Sensitive to any drift in the in-session port.
#
# Coverage:
#   * mu=0 (null) x 3 reps: probit mode, --outlier_threshold '-3,6', abs-K
#   * property injection x 2 reps: directional two-sided, '-6,6', signed-K
# Feature: a synthetic property matrix built from a random signature and a
# small membership matrix built from a sampled slice of the c2 gmt (entrez).

suppressMessages({
  library(optparse)
  library(data.table)
  library(Matrix)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option('--panel_prep', default = 'sim/output/panel_cache/Whole_Blood.prep.RDS',   type = 'character'),
  make_option('--K_abs',      default = 'sim/output/panel_cache/Whole_Blood.CorMat.RDS', type = 'character'),
  make_option('--K_signed',   default = 'sim/output/panel_cache/Whole_Blood.CorMatSigned.RDS', type = 'character'),
  make_option('--template_twas', default = '/data/twas_results',           type = 'character'),
  make_option('--gmt',        default = '/data/c2.all.v2026.1.Hs.entrez.gmt', type = 'character'),
  make_option('--tool',       default = 'TWAS-GSEA-fast.R',                 type = 'character'),
  make_option('--rscript',    default = 'sim/run_r.sh',                    type = 'character'),
  make_option('--outdir',     default = 'sim/output/parity',               type = 'character'),
  make_option('--n_null',     default = 3L, type = 'integer'),
  make_option('--n_signal',   default = 2L, type = 'integer'),
  # P-values are defined only to REML tolerance (default 1e-6 in
  # TWAS-GSEA-fast.R). Empirically, worst |dP| across our scenarios is O(1e-7)
  # on the one seed where Brent's search lands close to a boundary. 1e-5 is a
  # safe gate; results with |dP|>1e-5 would indicate a real bug, not tolerance.
  make_option('--tol_p',      default = 1e-5, type = 'numeric'),
  make_option('--tol_est',    default = 1e-6, type = 'numeric')
)))

SCRIPT_DIR <- dirname(sub('--file=', '', grep('--file=', commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
source(file.path(SCRIPT_DIR, 'R', 'sim_generator.R'))
source(file.path(SCRIPT_DIR, 'R', 'id_annotate.R'))
source(file.path(SCRIPT_DIR, 'R', 'feature_builders.R'))
source(file.path(SCRIPT_DIR, 'R', 'in_session_core.R'))
source(file.path(dirname(SCRIPT_DIR), 'R', 'reml_blockdiag.R'))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
LOG <- file.path(opt$outdir, 'parity.log')
log_msg <- function(...) { cat(..., '\n', sep = ''); cat(..., '\n', sep = '', file = LOG, append = TRUE) }
cat('', file = LOG)

prep     <- panel_prep(opt$panel_prep, cache_path = opt$panel_prep)
universe <- prep$gene_universe
tpl      <- load_template_twas(opt$template_twas, universe = universe)
tpl      <- tpl[match(universe, tpl$FILE), , drop = FALSE]

# --- Build a small sampled subset of the c2 gmt: keep 500 sets to keep the
# CLI-run wall clock manageable during the parity check.
gs_all <- read.gmt(opt$gmt)
set.seed(42L)
gs_small <- gs_all[sort(sample(seq_along(gs_all), min(500, length(gs_all))))]

gmt_small_path <- file.path(opt$outdir, 'gmt_small.gmt')
sink(gmt_small_path)
for(nm in names(gs_small)){
  cat(nm, 'sim', paste(gs_small[[nm]], collapse = '\t'), sep = '\t'); cat('\n')
}
sink()

# --- Build a synthetic property file (1000 signed columns keyed by Symbol).
set.seed(1L)
N <- length(universe)
n_props <- 200
prop_syn <- matrix(rnorm(N * n_props, sd = 1), nrow = N, ncol = n_props)
# Add sparsity to a subset — some columns are ternary +/- 1.
for(j in seq_len(n_props)){
  if(j %% 3 == 0) prop_syn[, j] <- sign(prop_syn[, j]) * (abs(prop_syn[, j]) > 0.5)
}
colnames(prop_syn) <- paste0('drug_', seq_len(n_props))
# Row-key: universe ENSG unversioned (matches --use_alt_id ID)
rownames(prop_syn) <- sub('\\..*$', '', universe)
prop_syn_path <- file.path(opt$outdir, 'prop_synthetic.rds')
saveRDS(prop_syn, prop_syn_path)

# ---------------------------------------------------------------------------
# Helper: run one simulation + both scoring paths + compare.
# ---------------------------------------------------------------------------
run_one <- function(scenario_name, seed, mu_fn, mode, two_sided, outlier_thr,
                    input_mode, cor_rds_path) {
  # input_mode: 'prop_syn' or 'gmt_small'

  # 1. Simulate one TWAS dataset (in memory + on disk).
  z <- simulate_z(prep, mu = mu_fn(prep), seed = seed)
  tw <- emit_twas(z, tpl, universe)
  tw <- annotate_twas(tw)
  # For --use_alt_id ID (unversioned ENSG) the prop file must key on that.
  tw$AltID_ENSG <- sub('\\..*$', '', tw$FILE)

  # CLI expects FILE in WGT-path form; add the prefix so its normalisation
  # (lines 118-147) reduces it back to the same string the in-session core sees.
  # R.utils not installed in this env, so avoid .gz -> plain .tsv is fine here.
  tw_out <- tw
  tw_out$FILE <- paste0('Whole_Blood/', tw$FILE, '.wgt.RDat')
  twas_tsv <- file.path(opt$outdir, sprintf('twas_%s_seed%d.tsv', scenario_name, seed))
  fwrite(tw_out, twas_tsv, sep = '\t')

  # 2. CLI run.
  cli_out <- file.path(opt$outdir, sprintf('cli_%s_seed%d', scenario_name, seed))
  cli_flags <- c(
    '--twas_results', twas_tsv,
    '--input_CorMat', cor_rds_path,
    '--output', cli_out,
    # Use = form to prevent optparse from treating negative outlier bounds as flags.
    sprintf('--outlier_threshold=%g,%g', outlier_thr[1], outlier_thr[2]),
    '--use_alt_id',
    if(input_mode == 'gmt_small') 'Entrez' else 'AltID_ENSG',
    '--two_sided', if(two_sided) 'T' else 'F')
  if(mode == 'directional'){
    cli_flags <- c(cli_flags, '--directional', 'T')
  } else if(mode == 'magnitude'){
    cli_flags <- c(cli_flags, '--probit_P_as_Z', 'F', '--directional', 'F')
  } # else probit is the default
  if(input_mode == 'gmt_small') cli_flags <- c(cli_flags, '--gmt_file', gmt_small_path)
  else                          cli_flags <- c(cli_flags, '--prop_file', prop_syn_path)

  cli_cmd <- c(opt$rscript, opt$tool, cli_flags)
  t0 <- Sys.time()
  status <- system2(cli_cmd[1], args = cli_cmd[-1], stdout = FALSE, stderr = FALSE)
  cli_secs <- as.numeric(difftime(Sys.time(), t0, units = 'secs'))
  if(status != 0) stop('CLI failed: ', paste(cli_cmd, collapse = ' '))
  cli_res <- fread(paste0(cli_out, '.competitive.txt'))

  # 3. In-session run using the SAME inputs and the same K. Mirror the CLI's
  # dedup-on-Alt_ID step (TWAS-GSEA-fast.R:167-171): sort by MODELCV.R2 asc,
  # then drop duplicate Alt_ID (this quietly drops rows where the alt-id column
  # has NA — ~986 of 7813 for Entrez on Whole_Blood).
  alt_col <- if(input_mode == 'gmt_small') 'Entrez' else 'AltID_ENSG'
  tw_dedup <- tw[order(tw$MODELCV.R2), ]
  tw_dedup <- tw_dedup[!duplicated(tw_dedup[[alt_col]]), ]

  if(input_mode == 'gmt_small'){
    feature_mat <- build_membership_matrix(tw_dedup, gmt_small_path, alt_id_col = 'Entrez')
    feature_type <- 'set'
  } else {
    prop_mat_loaded <- load_prop_file(prop_syn_path)
    feature_mat <- build_property_matrix(tw_dedup, prop_mat_loaded, alt_id_col = 'AltID_ENSG')
    feature_type <- 'prop'
  }
  tw <- tw_dedup
  # For parity we must use the SAME K (asymmetric, raw) as the CLI.
  K_here <- readRDS(cor_rds_path)$K

  t0 <- Sys.time()
  ins_res <- score_dataset(prep, tw, feature_mat, feature_type, K = K_here,
                            mode = mode, outlier_threshold = outlier_thr,
                            two_sided = two_sided,
                            p_cor_method = 'fdr')
  ins_secs <- as.numeric(difftime(Sys.time(), t0, units = 'secs'))

  # 4. Compare on shared GeneSet ids.
  both <- merge(as.data.table(cli_res)[, .(GeneSet, Estimate_cli = Estimate, P_cli = P,
                                            N_Mem_Avail_cli = N_Mem_Avail)],
                as.data.table(ins_res)[, .(GeneSet, Estimate_ins = Estimate, P_ins = P,
                                            N_Mem_Avail_ins = N_Mem_Avail)],
                by = 'GeneSet')
  d_p    <- max(abs(both$P_cli - both$P_ins))
  d_est  <- max(abs(both$Estimate_cli - both$Estimate_ins))
  d_nm   <- max(abs(both$N_Mem_Avail_cli - both$N_Mem_Avail_ins))
  n_only_cli <- nrow(cli_res) - nrow(both)
  n_only_ins <- nrow(ins_res) - nrow(both)

  log_msg(sprintf(
    '[%s seed=%d input=%s mode=%s two=%s outlier=%s cor=%s]  shared=%d  only-cli=%d  only-ins=%d  max|dP|=%.3g  max|dEst|=%.3g  max|dN|=%d  CLI=%.1fs  in=%.1fs',
    scenario_name, seed, input_mode, mode, if(two_sided) 'T' else 'F',
    paste(outlier_thr, collapse=','), basename(cor_rds_path),
    nrow(both), n_only_cli, n_only_ins, d_p, d_est, d_nm, cli_secs, ins_secs))

  list(name = scenario_name, seed = seed, d_p = d_p, d_est = d_est, d_nm = d_nm,
       n_only_cli = n_only_cli, n_only_ins = n_only_ins,
       cli_secs = cli_secs, ins_secs = ins_secs)
}

# ---------------------------------------------------------------------------
# Scenarios.
# ---------------------------------------------------------------------------
mu_null_fn     <- function(prep) mu_null(prep)
signature_seed <- 111L

# random signed property signature (fixed across signal runs)
set.seed(signature_seed)
n_nz <- max(1L, round(0.05 * N))
sig_idx <- sample.int(N, n_nz)
t_vec <- numeric(N); t_vec[sig_idx] <- rnorm(n_nz)
t_vec <- (t_vec - mean(t_vec)) / sd(t_vec)
names(t_vec) <- universe
mu_signal_fn <- function(prep) mu_property(prep, t_vec, beta = 0.7)

results <- list()

# Null x N_null: probit + abs-K + property
for(s in seq_len(opt$n_null)){
  results[[length(results) + 1]] <- run_one(
    scenario_name = 'null_prop_probit_absK',
    seed         = 1000 + s,
    mu_fn        = mu_null_fn,
    mode         = 'probit',
    two_sided    = FALSE,
    outlier_thr  = c(-3, 6),
    input_mode   = 'prop_syn',
    cor_rds_path = opt$K_abs)
}

# Null x 1: probit + abs-K + gmt
results[[length(results) + 1]] <- run_one(
  scenario_name = 'null_gmt_probit_absK',
  seed         = 2001,
  mu_fn        = mu_null_fn,
  mode         = 'probit',
  two_sided    = FALSE,
  outlier_thr  = c(-3, 6),
  input_mode   = 'gmt_small',
  cor_rds_path = opt$K_abs)

# Signal x N_signal: directional two-sided + signed-K + property
for(s in seq_len(opt$n_signal)){
  results[[length(results) + 1]] <- run_one(
    scenario_name = 'signal_prop_dir2s_signedK',
    seed         = 3000 + s,
    mu_fn        = mu_signal_fn,
    mode         = 'directional',
    two_sided    = TRUE,
    outlier_thr  = c(-6, 6),
    input_mode   = 'prop_syn',
    cor_rds_path = opt$K_signed)
}

# Signal x 1: directional one-sided + abs-K + gmt
results[[length(results) + 1]] <- run_one(
  scenario_name = 'signal_gmt_dir1s_absK',
  seed         = 4001,
  mu_fn        = mu_signal_fn,
  mode         = 'directional',
  two_sided    = FALSE,
  outlier_thr  = c(-6, 6),
  input_mode   = 'gmt_small',
  cor_rds_path = opt$K_abs)

# ---------------------------------------------------------------------------
# Gate.
# ---------------------------------------------------------------------------
res_df <- do.call(rbind, lapply(results, function(r) data.frame(
  name = r$name, seed = r$seed, d_p = r$d_p, d_est = r$d_est, d_nm = r$d_nm,
  n_only_cli = r$n_only_cli, n_only_ins = r$n_only_ins,
  cli_secs = r$cli_secs, ins_secs = r$ins_secs)))
fwrite(res_df, file.path(opt$outdir, 'parity_summary.tsv'), sep = '\t')

pass_p   <- all(res_df$d_p   <= opt$tol_p)
pass_est <- all(res_df$d_est <= opt$tol_est)
pass_nm  <- all(res_df$d_nm  == 0)
pass_row <- all(res_df$n_only_cli == 0 & res_df$n_only_ins == 0)

log_msg('---')
log_msg(sprintf('GATE max|dP|  <= %g  : %s   (worst=%.3g)', opt$tol_p,   pass_p,   max(res_df$d_p)))
log_msg(sprintf('GATE max|dEst|<= %g  : %s   (worst=%.3g)', opt$tol_est, pass_est, max(res_df$d_est)))
log_msg(sprintf('GATE N_Mem_Avail identical : %s   (worst=%d)', pass_nm, max(res_df$d_nm)))
log_msg(sprintf('GATE row sets match       : %s', pass_row))
log_msg('OVERALL PARITY: ', if(pass_p && pass_est && pass_nm && pass_row) 'PASS' else 'FAIL')
