#!/usr/bin/Rscript
# sim/aggregate.R
#
# CLI wrapper over sim/R/aggregate.R. Writes tidy .tsv tables.

suppressMessages({
  library(optparse)
  library(data.table)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option('--grid_dir',      default = 'sim/output/grid', type = 'character'),
  make_option('--scenarios',     default = 'sim/scenarios.tsv', type = 'character'),
  make_option('--outdir',        default = 'sim/output/aggregate', type = 'character')
)))

SCRIPT_DIR <- dirname(sub('--file=', '', grep('--file=', commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
source(file.path(SCRIPT_DIR, 'R', 'aggregate.R'))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
res <- aggregate_grid(opt$grid_dir, opt$scenarios)
fwrite(res$scenario_summary, file.path(opt$outdir, 'scenario_summary.tsv'), sep = '\t')
fwrite(res$power_summary,    file.path(opt$outdir, 'power_summary.tsv'),    sep = '\t')

cat('scenario_summary rows:', nrow(res$scenario_summary), '\n')
cat('power_summary rows:',    nrow(res$power_summary),    '\n')
cat('wrote', opt$outdir, '\n')

# ---- Console preview: split by grouping ---
if(nrow(res$scenario_summary) > 0){
  cat('\n=== T1E null calibration (frac_p_lt_05 should be ~ 0.05) ===\n')
  print(res$scenario_summary[injection_type == 'null',
    .(scenario_id, arm, feature_source, n_reps, n_features,
      frac_p_lt_05, fwer_at_05, mean_n_fd, pooled_p_ks_p)])
}
if(nrow(res$power_summary) > 0){
  cat('\n=== Directional headline power (dh_*) ===\n')
  print(res$power_summary[order(feature_source, arm, cor_matrix, beta),
    .(scenario_id, feature_source, arm, cor_matrix, beta, n_reps,
      mean_est, sd_est, power_p05, power_p_corr05,
      frac_sign_pos, sign_recovery)])
}
