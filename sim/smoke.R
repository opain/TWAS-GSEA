#!/usr/bin/Rscript
# sim/smoke.R
#
# Phase 3 gate. Runs 20 null + 20 signal reps against the FULL c2 gmt and the
# FULL CMAP + ternary property panels, measuring per-task wall time and peak
# RSS. Projects full-grid wall clock at 10-way concurrency and reports.

suppressMessages({
  library(optparse)
  library(data.table)
  library(Matrix)
  library(parallel)
  library(ps)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option('--outdir',   default = 'sim/output/smoke', type = 'character'),
  make_option('--n_cores',  default = 10L, type = 'integer'),
  make_option('--n_null',   default = 20L, type = 'integer'),
  make_option('--n_signal', default = 20L, type = 'integer')
)))

SCRIPT_DIR <- dirname(sub('--file=', '', grep('--file=', commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Emit a smoke scenarios.tsv, hand it to run_grid.R, capture the log.
# ---------------------------------------------------------------------------
smoke_scen <- data.table(
  scenario_id  = c('smoke_null_cmap', 'smoke_null_gmt', 'smoke_null_ternary',
                    'smoke_sig_cmap',  'smoke_sig_gmt',  'smoke_sig_ternary'),
  arm          = c('magnitude_probit_absK', 'magnitude_probit_absK', 'magnitude_probit_absK',
                    'directional_2s_signedK','directional_1s_absK',   'directional_2s_signedK'),
  mode           = c('probit','probit','probit','directional','directional','directional'),
  two_sided      = c(FALSE, FALSE, FALSE, TRUE,  FALSE, TRUE),
  cor_matrix     = c('abs',   'abs',  'abs',   'signed','abs',   'signed'),
  outlier_lo     = c(-3, -3, -3, -6, -6, -6),
  outlier_hi     = c(6, 6, 6, 6, 6, 6),
  min_Ngenes     = 2L,
  injection_type = c('null','null','null','property_random','set','property_random'),
  feature_source = c('cmap_a375_full','gmt_c2_full','ternary_full',
                     'cmap_a375_full','gmt_c2_full','ternary_full'),
  feature_alt_id = c('Symbol','Entrez','Symbol','Symbol','Entrez','Symbol'),
  target_size    = c(NA, NA, NA, NA, 100L, NA),
  Delta          = c(NA, NA, NA, NA, 1.5,  NA),
  rho            = c(NA, NA, NA, NA, 0.5,  NA),
  pi_pos         = c(NA, NA, NA, NA, 0.6,  NA),
  sparsity       = c(NA, NA, NA, 0.05, NA, 0.05),
  sigma_obs      = NA_real_,
  beta           = c(NA, NA, NA, 0.4, NA, 0.4),
  signature_seed = c(NA, NA, NA, 111L, NA, 111L),
  n_reps         = c(opt$n_null, opt$n_null, opt$n_null,
                     opt$n_signal, opt$n_signal, opt$n_signal),
  notes          = c('null-cmap','null-gmt','null-ternary','signal-cmap','signal-gmt','signal-ternary'))

scen_path <- file.path(opt$outdir, 'smoke_scenarios.tsv')
fwrite(smoke_scen, scen_path, sep = '\t')

# ---------------------------------------------------------------------------
# Package check — mandatory per plan.
# ---------------------------------------------------------------------------
missing <- c()
for(p in c('data.table','Matrix','WGCNA','VGAM','qusage','foreach','doMC','matrixcalc','ps')){
  if(!requireNamespace(p, quietly = TRUE)) missing <- c(missing, p)
}
if(length(missing) > 0) stop('Missing R packages: ', paste(missing, collapse = ', '))
cat('Packages OK.\n')

# ---------------------------------------------------------------------------
# Run the grid, measuring peak RSS across workers.
# ---------------------------------------------------------------------------
run_grid_path <- file.path(SCRIPT_DIR, 'run_grid.R')
child_cmd <- c(
  file.path(SCRIPT_DIR, 'run_r.sh'),
  run_grid_path,
  '--scenarios', scen_path,
  '--outdir', file.path(opt$outdir, 'grid'),
  '--n_cores', opt$n_cores)

# Launch as a subprocess so we can sample RSS from the parent process tree.
smoke_log <- file.path(opt$outdir, 'smoke.log')
proc <- ps::ps()  # baseline
child <- processx::process$new(child_cmd[1], child_cmd[-1],
                                stdout = smoke_log, stderr = '2>&1',
                                supervise = TRUE)
peak_rss <- 0
t0 <- Sys.time()
while(child$is_alive()){
  Sys.sleep(2)
  # Poll the process tree from ps
  tryCatch({
    pids <- ps::ps_children(ps::ps_handle(child$get_pid()), recursive = TRUE)
    rss_now <- 0
    for(h in pids) rss_now <- rss_now + ps::ps_memory_info(h)$rss
    # Add parent
    rss_now <- rss_now + ps::ps_memory_info(ps::ps_handle(child$get_pid()))$rss
    if(rss_now > peak_rss) peak_rss <- rss_now
  }, error = function(e) NULL)
}
wall <- as.numeric(difftime(Sys.time(), t0, units = 'secs'))
child$wait()
status <- child$get_exit_status()

# NOTE: the process-tree RSS walk via ps_children often returns 0 when the
# child ran through the env -i wrapper (grandchildren exit before we sample).
# For a reliable per-worker peak, use /proc/*/status externally or set
# ulimit + RUsage. Kept for a rough upper bound only.
cat(sprintf('smoke run wall = %.1f s   exit = %d   peak RSS (process tree) = %.2f GB (may under-report)\n',
            wall, status, peak_rss / 2^30))

# ---------------------------------------------------------------------------
# Load per-scenario long.rds files, produce timing table + projection.
# ---------------------------------------------------------------------------
long_dir <- file.path(opt$outdir, 'grid')
long_files <- list.files(long_dir, pattern = '\\.long\\.rds$', full.names = TRUE)
if(length(long_files) == 0){
  cat('No long.rds files produced — check smoke.log\n')
  quit(status = 1)
}
long_all <- rbindlist(lapply(long_files, readRDS), use.names = TRUE, fill = TRUE)
# Timing summary (one row per replicate; wall_secs already recorded)
by_scen <- long_all[, .(reps = uniqueN(replicate),
                         mean_wall = mean(wall_secs, na.rm = TRUE),
                         med_wall  = median(wall_secs, na.rm = TRUE),
                         max_wall  = max(wall_secs, na.rm = TRUE)),
                     by = scenario_id]
print(by_scen)
fwrite(by_scen, file.path(opt$outdir, 'timing_by_scenario.tsv'), sep = '\t')

# Projection: for each scenario in a hypothetical FULL grid, extrapolate.
# Placeholder full-grid sizes: T1E null needs 1000 reps per arm/feature; power
# needs 300 reps per cell over a modest grid; directional headline 300 reps per
# cell over a signed-beta sweep. Report a coarse projection.
proj <- data.table(
  arm         = by_scen$scenario_id,
  reps_per_scenario = c(rep(1000, 3), rep(300, 3))[seq_len(nrow(by_scen))],
  # `arms x cells` multiplier for the fully-populated grid (rough).
  n_cells     = c(4, 4, 4, 8, 6, 8)[seq_len(nrow(by_scen))],
  mean_secs   = by_scen$mean_wall)
proj[, secs_total := reps_per_scenario * n_cells * mean_secs]
proj[, hours_wallclock_10way := secs_total / (opt$n_cores * 3600)]
print(proj)
fwrite(proj, file.path(opt$outdir, 'projection.tsv'), sep = '\t')

cat(sprintf('\nSUMMARY\n  smoke wall = %.1f s (%d cores)\n  peak RSS   = %.2f GB\n  projected full-grid at %d-way = %.1f hours\n',
            wall, opt$n_cores, peak_rss / 2^30, opt$n_cores, sum(proj$hours_wallclock_10way)))
