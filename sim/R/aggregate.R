# sim/R/aggregate.R
#
# Reads per-scenario long.rds outputs from sim/output/grid/, computes the
# study's estimands, and returns tidy tables.

suppressMessages({
  library(data.table)
})

# For each scenario, compute T1E and power summaries. Signals are recognised by
# scenario_id conventions (t1e_* = null; dh_* have beta in scenarios.tsv).
aggregate_grid <- function(grid_dir, scenarios_tsv) {
  scen <- fread(scenarios_tsv)
  files <- list.files(grid_dir, pattern = '\\.long\\.rds$', full.names = TRUE)
  if(length(files) == 0) stop('no long.rds under ', grid_dir)

  per_scen <- list()
  per_rep  <- list()
  for(f in files){
    sid_from_file <- sub('\\.long\\.rds$', '', basename(f))
    x <- as.data.table(readRDS(f))
    x[, scenario_id := sid_from_file]
    # attach scenario metadata
    s <- scen[scenario_id == sid_from_file]
    if(nrow(s) == 0) next
    x[, `:=`(arm            = s$arm,
             mode           = s$mode,
             two_sided      = s$two_sided,
             cor_matrix     = s$cor_matrix,
             injection_type = s$injection_type,
             feature_source = s$feature_source,
             beta           = s$beta,
             signature_seed = s$signature_seed,
             n_reps         = s$n_reps)]

    # Per-scenario summary (across replicates x features): pooled uniform-P,
    # per-rep FWER, mean number of false discoveries at 0.05 FDR.
    pooled_ks <- tryCatch(
      suppressWarnings(ks.test(x$P, 'punif')$p.value),
      error = function(e) NA_real_)
    frac_p_lt_alpha <- mean(x$P < 0.05, na.rm = TRUE)
    # per-rep FWER at raw P (dependent tests but useful for calibration diag)
    fwer_p <- x[, .(hit = any(P.CORR <= 0.05, na.rm = TRUE),
                    n_fd = sum(P.CORR <= 0.05, na.rm = TRUE)),
                by = replicate]
    per_scen[[length(per_scen) + 1L]] <- data.table(
      scenario_id = sid_from_file,
      arm = s$arm, mode = s$mode, two_sided = s$two_sided,
      cor_matrix = s$cor_matrix, injection_type = s$injection_type,
      feature_source = s$feature_source, beta = s$beta,
      n_reps = uniqueN(x$replicate),
      n_features = uniqueN(x$GeneSet),
      pooled_p_ks_p = pooled_ks,
      frac_p_lt_05  = frac_p_lt_alpha,
      fwer_at_05    = mean(fwer_p$hit),
      mean_n_fd     = mean(fwer_p$n_fd))

    # For directional headline (single-feature scenarios), also compute
    # signed-Estimate + power summaries.
    if(nrow(x) > 0 && s$injection_type != 'null' && !is.na(s$beta)){
      # 'target' is the single feature the injection was drawn from.
      pow <- x[, .(power_p05     = mean(P     < 0.05,  na.rm = TRUE),
                    power_p_corr05= mean(P.CORR < 0.05,  na.rm = TRUE),
                    mean_est      = mean(Estimate, na.rm = TRUE),
                    sd_est        = sd(Estimate, na.rm = TRUE),
                    n_reps        = uniqueN(replicate),
                    frac_sign_pos = mean(Estimate > 0, na.rm = TRUE),
                    sign_recovery = if(s$beta != 0) mean(sign(Estimate) == sign(s$beta) & P.CORR < 0.05, na.rm = TRUE) else NA_real_),
                by = GeneSet]
      pow[, `:=`(scenario_id = sid_from_file, arm = s$arm, mode = s$mode,
                 two_sided = s$two_sided, cor_matrix = s$cor_matrix,
                 feature_source = s$feature_source, beta = s$beta)]
      per_rep[[length(per_rep) + 1L]] <- pow
    }
  }

  scenario_summary <- rbindlist(per_scen, use.names = TRUE, fill = TRUE)
  power_summary    <- rbindlist(per_rep,  use.names = TRUE, fill = TRUE)
  list(scenario_summary = scenario_summary,
       power_summary    = power_summary)
}
