#!/usr/bin/Rscript
# sim/run_grid.R
#
# Grid runner: expands sim/scenarios.tsv x n_reps into a task list and executes
# each task in a forked worker (mclapply). Each task simulates one dataset,
# scores it via sim/R/in_session_core.R::score_dataset(), and returns tidy long
# rows.
#
# Parallelism strategy:
#   * one level only — across replicates, at concurrency = --n_cores (default 10)
#   * each task is SINGLE-THREADED (BLAS threads capped via env vars from run_r.sh)
#   * panel_prep, K matrices, gmt membership, prop matrices loaded ONCE in the
#     parent process and inherited via copy-on-write after fork
#
# Scenarios.tsv columns:
#   scenario_id, arm, mode, two_sided, cor_matrix, outlier_lo, outlier_hi,
#   min_Ngenes, injection_type, feature_source, feature_alt_id,
#   target_size, Delta, rho, pi_pos, sparsity, sigma_obs, beta,
#   signature_seed, n_reps, notes
#
# Feature sources (fixed set, resolved once at grid start):
#   'synthetic_prop'  — random signed matrix (see build_synthetic_prop_mat())
#   'gmt_c2_full'     — /data/c2.all.v2026.1.Hs.entrez.gmt (entrez, all sets)
#   'gmt_c2_sampled'  — a fixed random subset of c2 sets (for cheap tests)
#   'cmap_a375_full'  — /data/lvl5Allcompounds_...A375_24h_10uM.rds (symbol)
#   'ternary_full'    — /data/wholedatabase_for_targetor_directional.prop (symbol)

suppressMessages({
  library(optparse)
  library(data.table)
  library(Matrix)
  library(parallel)
  library(qusage)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option('--scenarios', default = NA, type = 'character', help = 'TSV with scenario grid [required]'),
  make_option('--panel_prep', default = 'sim/output/panel_cache/Whole_Blood.prep.RDS', type = 'character'),
  make_option('--K_abs',      default = 'sim/output/panel_cache/Whole_Blood.CorMat.RDS', type = 'character'),
  make_option('--K_signed',   default = 'sim/output/panel_cache/Whole_Blood.CorMatSigned.RDS', type = 'character'),
  make_option('--template_twas', default = '/data/twas_results', type = 'character'),
  make_option('--ensg2entrez',  default = 'sim/output/id_map/ensg2entrez.tsv', type = 'character'),
  make_option('--biomart',      default = '/data/biomart/biomart_genes_grch37.tsv', type = 'character'),
  make_option('--outdir',       default = 'sim/output/grid', type = 'character'),
  make_option('--n_cores',      default = 10L, type = 'integer'),
  make_option('--master_seed',  default = 20260911L, type = 'integer'),
  make_option('--only_scenario', default = '', type = 'character', help = 'Comma-separated scenario_ids to run [default all]')
)))
if(is.na(opt$scenarios)) stop('--scenarios is required')
SCRIPT_DIR <- dirname(sub('--file=', '', grep('--file=', commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
source(file.path(SCRIPT_DIR, 'R', 'sim_generator.R'))
source(file.path(SCRIPT_DIR, 'R', 'id_annotate.R'))
source(file.path(SCRIPT_DIR, 'R', 'feature_builders.R'))
source(file.path(SCRIPT_DIR, 'R', 'in_session_core.R'))
source(file.path(dirname(SCRIPT_DIR), 'R', 'reml_blockdiag.R'))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
GRID_LOG <- file.path(opt$outdir, 'run_grid.log')
log_msg <- function(...) { msg <- paste0(..., '\n'); cat(msg); cat(msg, file = GRID_LOG, append = TRUE) }
cat('', file = GRID_LOG)
log_msg('run_grid.R starting at ', as.character(Sys.time()), '  n_cores=', opt$n_cores)

# ---------------------------------------------------------------------------
# Read scenarios, filter if requested.
# ---------------------------------------------------------------------------
scen <- fread(opt$scenarios)
if(nzchar(opt$only_scenario)){
  wanted <- strsplit(opt$only_scenario, ',')[[1]]
  scen <- scen[scenario_id %in% wanted]
}
log_msg(nrow(scen), ' scenarios in this run.')

# ---------------------------------------------------------------------------
# Load panel once.
# ---------------------------------------------------------------------------
prep     <- panel_prep(opt$panel_prep, cache_path = opt$panel_prep, log_msg = function(...) log_msg(...))
universe <- prep$gene_universe
N        <- length(universe)
K_abs    <- readRDS(opt$K_abs)$K
K_signed <- readRDS(opt$K_signed)$K
stopifnot(identical(rownames(K_abs), universe), identical(rownames(K_signed), universe))

# Template TWAS (positional + MODELCV.R2 metadata reused per replicate)
tpl <- load_template_twas(opt$template_twas, universe = universe)
tpl <- tpl[match(universe, tpl$FILE), , drop = FALSE]

# ID annotation once — cached via id_annotate.R
tpl_ann <- annotate_twas(tpl, ensg2entrez_path = opt$ensg2entrez, biomart_path = opt$biomart)
tpl_ann$AltID_ENSG <- sub('\\..*$', '', tpl_ann$FILE)
alt_id_by_source <- list(
  synthetic_prop = 'AltID_ENSG',
  gmt_c2_full    = 'Entrez',
  gmt_c2_sampled = 'Entrez',
  cmap_a375_full = 'Symbol',
  ternary_full   = 'Symbol',
  # Signature-derived feature sources (single column per matrix). Keyed on
  # AltID_ENSG so we don't have to reconcile with real drug ID space.
  sig_cont       = 'AltID_ENSG',
  sig_tern       = 'AltID_ENSG',
  sig_bin        = 'AltID_ENSG'
)

# Deterministic signature vector used both for injection (mu_property) and
# for the sig_cont/sig_tern/sig_bin observed feature columns.
# Preserves zeros (unstandardised in raw form) so tau=0 gives a clean ternary.
get_signature_vector <- function(signature_seed, sparsity) {
  stopifnot(!is.na(signature_seed), !is.na(sparsity))
  set.seed(as.integer(signature_seed))
  n_nz <- max(1L, round(sparsity * N))
  idx  <- sample.int(N, n_nz)
  t_raw <- numeric(N)
  t_raw[idx] <- rnorm(n_nz)
  # Standardise WITHOUT shifting: divide by empirical sd only. Preserves zero
  # entries (sign(0) == 0, so ternary stays clean).
  t_std <- t_raw / sd(t_raw)
  names(t_std) <- universe
  attr(t_std, 'signature_seed') <- as.integer(signature_seed)
  attr(t_std, 'sparsity')       <- sparsity
  t_std
}

# ---------------------------------------------------------------------------
# Feature-source registry. Each entry returns a numeric matrix keyed by
# TWAS FILE rows (already aligned to the dedup'd TWAS in the caller).
# ---------------------------------------------------------------------------
feature_cache <- new.env(parent = emptyenv())

get_feature_mat <- function(source_id, tw_dedup, scen_row = NULL) {
  # sig_* sources depend on scenario-level params (signature_seed, sparsity,
  # sigma_obs) — key on those too so the cache hits precisely.
  extra_key <- ''
  if(!is.null(scen_row) && grepl('^sig_', source_id)){
    extra_key <- sprintf('|sig=%s|sp=%s|sig_obs=%s',
                          scen_row$signature_seed, scen_row$sparsity, scen_row$sigma_obs)
  }
  key <- paste(source_id, nrow(tw_dedup), extra_key, sep = '|')
  if(!is.null(feature_cache[[key]])) return(feature_cache[[key]])

  alt_col <- alt_id_by_source[[source_id]]
  fmat <- switch(source_id,
    synthetic_prop = {
      # 200 columns; some ternary. Fixed rownames = universe unversioned ENSG.
      set.seed(1L)
      n_props <- 200L
      p <- matrix(rnorm(N * n_props), nrow = N, ncol = n_props)
      for(j in seq_len(n_props)) if(j %% 3 == 0) p[, j] <- sign(p[, j]) * (abs(p[, j]) > 0.5)
      colnames(p) <- paste0('drug_', seq_len(n_props))
      rownames(p) <- sub('\\..*$', '', universe)
      build_property_matrix(tw_dedup, p, alt_id_col = alt_col)
    },
    gmt_c2_full = build_membership_matrix(tw_dedup, '/data/c2.all.v2026.1.Hs.entrez.gmt', alt_id_col = alt_col),
    gmt_c2_sampled = {
      gs_all <- read.gmt('/data/c2.all.v2026.1.Hs.entrez.gmt')
      set.seed(42L)
      gs_sm  <- gs_all[sort(sample(seq_along(gs_all), min(500, length(gs_all))))]
      sm_path <- file.path(opt$outdir, 'gmt_c2_sampled.gmt')
      sink(sm_path); for(nm in names(gs_sm)) { cat(nm, 'sim', paste(gs_sm[[nm]], collapse='\t'), sep='\t'); cat('\n') } ; sink()
      build_membership_matrix(tw_dedup, sm_path, alt_id_col = alt_col)
    },
    cmap_a375_full = {
      prop <- load_prop_file('/data/lvl5Allcompounds_population_TWAS-GSEA_prop.A375_24h_10uM.rds')
      build_property_matrix(tw_dedup, prop, alt_id_col = alt_col)
    },
    ternary_full = {
      prop <- load_prop_file('/data/wholedatabase_for_targetor_directional.prop')
      build_property_matrix(tw_dedup, prop, alt_id_col = alt_col)
    },
    sig_cont = {
      t_vec <- get_signature_vector(scen_row$signature_seed, scen_row$sparsity)
      # Noise seeded off the signature seed so it's reproducible but independent
      # of any particular replicate.
      set.seed(as.integer(scen_row$signature_seed) + 88711L)
      noisy <- t_vec + rnorm(N, sd = scen_row$sigma_obs)
      m <- matrix(noisy, ncol = 1)
      colnames(m) <- 'target_cont'
      rownames(m) <- sub('\\..*$', '', universe)
      build_property_matrix(tw_dedup, m, alt_id_col = alt_col)
    },
    sig_tern = {
      t_vec <- get_signature_vector(scen_row$signature_seed, scen_row$sparsity)
      # Threshold at the (1-sparsity) |t| quantile so the ternary retains the
      # true nonzero signature entries. tau=0 works exactly when t is raw-sparse,
      # but with sd-normalisation there's tiny zero-drift, so quantile is safer.
      tau <- quantile(abs(t_vec), 1 - scen_row$sparsity)
      tern <- ifelse(abs(t_vec) > tau, sign(t_vec), 0)
      m <- matrix(tern, ncol = 1)
      colnames(m) <- 'target_tern'
      rownames(m) <- sub('\\..*$', '', universe)
      build_property_matrix(tw_dedup, m, alt_id_col = alt_col)
    },
    sig_bin = {
      t_vec <- get_signature_vector(scen_row$signature_seed, scen_row$sparsity)
      tau <- quantile(abs(t_vec), 1 - scen_row$sparsity)
      bin <- as.numeric(abs(t_vec) > tau)
      m <- matrix(bin, ncol = 1)
      colnames(m) <- 'target_bin'
      rownames(m) <- sub('\\..*$', '', universe)
      build_property_matrix(tw_dedup, m, alt_id_col = alt_col)
    },
    stop('unknown feature_source: ', source_id))

  feature_cache[[key]] <- fmat
  fmat
}

# ---------------------------------------------------------------------------
# Task helpers.
# ---------------------------------------------------------------------------
build_mu <- function(scen_row, seed) {
  set.seed(seed)
  switch(as.character(scen_row$injection_type),
    null = mu_null(prep),
    set = {
      target <- sample(universe, scen_row$target_size)
      mu_set(prep, target, scen_row$Delta,
             rho    = if(is.na(scen_row$rho))    1 else scen_row$rho,
             pi_pos = if(is.na(scen_row$pi_pos)) 0.5 else scen_row$pi_pos)
    },
    property_random = {
      sig_seed <- if(is.na(scen_row$signature_seed)) 111L else as.integer(scen_row$signature_seed)
      s        <- if(is.na(scen_row$sparsity))       0.05 else scen_row$sparsity
      t_vec <- get_signature_vector(sig_seed, s)
      set.seed(seed)
      mu_property(prep, t_vec, scen_row$beta)
    },
    stop('unknown injection_type: ', scen_row$injection_type))
}

# Deterministic per-scenario integer seed base (hash to keep unique across runs).
scenario_seed_base <- function(scenario_id, master_seed){
  # Sum of ASCII codes + master; small collision risk between similar scen ids,
  # but we combine with replicate below so per-task seeds are unique.
  as.integer((sum(utf8ToInt(as.character(scenario_id))) * 1009L + master_seed) %% .Machine$integer.max)
}

# Simulate one dataset, dedup, get feature_mat, run in-session core.
run_task <- function(scen_row, replicate) {
  seed <- as.integer((scenario_seed_base(scen_row$scenario_id, opt$master_seed) +
                       replicate * 7919L) %% .Machine$integer.max)
  mu <- build_mu(scen_row, seed)
  set.seed(seed + 1L)                        # separate stream for the MVN draw
  z <- simulate_z(prep, mu = mu)
  tw <- emit_twas(z, tpl_ann, universe)      # emit onto annotated template
  tw$AltID_ENSG <- tpl_ann$AltID_ENSG        # emit_twas may drop non-standard cols
  # replay Symbol/Entrez if emit dropped them; annotate_twas already ran on tpl.
  tw$Symbol <- tpl_ann$Symbol
  tw$Entrez <- tpl_ann$Entrez

  alt_col <- alt_id_by_source[[as.character(scen_row$feature_source)]]
  tw_dedup <- tw[order(tw$MODELCV.R2), ]
  tw_dedup <- tw_dedup[!duplicated(tw_dedup[[alt_col]]), ]

  fmat <- get_feature_mat(as.character(scen_row$feature_source), tw_dedup, scen_row = scen_row)
  ftype <- if(as.character(scen_row$feature_source) %in% c('gmt_c2_full','gmt_c2_sampled','sig_bin')) 'set' else 'prop'
  K_here <- if(scen_row$cor_matrix == 'signed') K_signed else K_abs

  t0 <- Sys.time()
  res <- tryCatch(
    score_dataset(prep, tw_dedup, fmat, ftype, K = K_here,
                  mode = as.character(scen_row$mode),
                  outlier_threshold = c(scen_row$outlier_lo, scen_row$outlier_hi),
                  two_sided  = as.logical(scen_row$two_sided),
                  min_Ngenes = as.integer(scen_row$min_Ngenes),
                  p_cor_method = 'fdr'),
    error = function(e) { attr(e, 'traceback') <- sys.calls(); e })
  wall <- as.numeric(difftime(Sys.time(), t0, units = 'secs'))
  if(inherits(res, 'error')){
    return(data.table(scenario_id = scen_row$scenario_id, replicate = replicate,
                      seed = seed, error = conditionMessage(res)))
  }
  # Long-format rows
  out <- as.data.table(res)
  out[, `:=`(scenario_id = scen_row$scenario_id,
             replicate   = replicate,
             seed        = seed,
             wall_secs   = wall)]
  # is_target: for set-based, mark the target set if the scenario names one.
  # (Left blank here; the aggregator resolves target sets by naming convention.)
  out
}

# ---------------------------------------------------------------------------
# Expand tasks and run.
# ---------------------------------------------------------------------------
tasks <- do.call(rbind, lapply(seq_len(nrow(scen)), function(i){
  sr <- scen[i]
  data.table(row = i, replicate = seq_len(sr$n_reps))
}))
log_msg('total tasks: ', nrow(tasks), '  (approx ', round(nrow(tasks)/opt$n_cores), ' per worker at ', opt$n_cores, '-way)')

t_run0 <- Sys.time()
out_all <- mclapply(seq_len(nrow(tasks)), function(k){
  scen_row <- scen[tasks$row[k]]
  res <- tryCatch(run_task(scen_row, replicate = tasks$replicate[k]),
                  error = function(e){
                    data.table(scenario_id = scen_row$scenario_id,
                               replicate = tasks$replicate[k],
                               error = conditionMessage(e))
                  })
  # Also wrap when mclapply itself returns a try-error object
  if(inherits(res, 'try-error')){
    res <- data.table(scenario_id = scen_row$scenario_id,
                       replicate = tasks$replicate[k],
                       error = attr(res, 'condition')$message)
  }
  res
}, mc.cores = opt$n_cores, mc.preschedule = FALSE, mc.silent = FALSE)
t_run1 <- Sys.time()
log_msg('mclapply done in ', round(as.numeric(difftime(t_run1, t_run0, units='secs')), 1), 's')

# Split into per-scenario tables and save.
combined <- rbindlist(out_all, use.names = TRUE, fill = TRUE)
# Report errors separately.
err_mask <- !is.na(combined$error) & is.na(combined$P)
if(sum(err_mask) > 0){
  errs <- combined[err_mask]
  fwrite(errs, file.path(opt$outdir, 'errors.tsv'), sep = '\t')
  log_msg('WARNING: ', nrow(errs), ' tasks failed. See errors.tsv')
  combined <- combined[!err_mask]
}
for(sid in unique(combined$scenario_id)){
  sub <- combined[scenario_id == sid]
  saveRDS(sub, file.path(opt$outdir, sprintf('%s.long.rds', sid)))
}
log_msg('wrote per-scenario long.rds files under ', opt$outdir)
log_msg('run_grid.R finished at ', as.character(Sys.time()))
