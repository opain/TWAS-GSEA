#!/usr/bin/Rscript
# sim/write_scenarios.R
#
# Generate sim/scenarios.tsv covering the three main experimental groups:
#   A) T1E on real drug property panels (cmap + ternary)
#   B) T1E on gmt (set-based)
#   C) Directional headline: sig_cont / sig_tern / sig_bin  ×  mode × K × beta

suppressMessages(library(data.table))

# -----------------------------------------------------------
# Base template with every column the runner recognises. NAs where unused.
# -----------------------------------------------------------
tmpl <- function(...) {
  base <- data.table(
    scenario_id = NA_character_, arm = NA_character_,
    mode = NA_character_, two_sided = FALSE, cor_matrix = NA_character_,
    outlier_lo = -3, outlier_hi = 6, min_Ngenes = 2L,
    injection_type = 'null', feature_source = NA_character_, feature_alt_id = NA_character_,
    target_size = NA_integer_, Delta = NA_real_, rho = NA_real_, pi_pos = NA_real_,
    sparsity = NA_real_, sigma_obs = NA_real_, beta = NA_real_,
    signature_seed = NA_integer_, n_reps = NA_integer_, notes = NA_character_)
  arg <- list(...)
  for(k in names(arg)) base[[k]] <- arg[[k]]
  base
}

rows <- list()

# ===========================================================
# GROUP A: T1E on real drug panels (property mode)
# ===========================================================
N_A <- 1000L
for(fs in c('cmap_a375_full', 'ternary_full')){
  alt <- 'Symbol'
  # mag (probit) + directional 2-sided abs + directional 2-sided signed
  rows[[length(rows)+1L]] <- tmpl(scenario_id = sprintf('t1e_%s_mag',           fs),
    arm = 'mag_probit_absK',    mode = 'probit',      two_sided = FALSE,
    cor_matrix = 'abs',    outlier_lo = -3, outlier_hi = 6,
    injection_type = 'null', feature_source = fs, feature_alt_id = alt,
    n_reps = N_A, notes = 'T1E null: magnitude/probit')
  rows[[length(rows)+1L]] <- tmpl(scenario_id = sprintf('t1e_%s_dir2s_abs',     fs),
    arm = 'dir2s_absK',         mode = 'directional', two_sided = TRUE,
    cor_matrix = 'abs',    outlier_lo = -6, outlier_hi = 6,
    injection_type = 'null', feature_source = fs, feature_alt_id = alt,
    n_reps = N_A, notes = 'T1E null: directional two-sided, abs-K')
  rows[[length(rows)+1L]] <- tmpl(scenario_id = sprintf('t1e_%s_dir2s_signed',  fs),
    arm = 'dir2s_signedK',      mode = 'directional', two_sided = TRUE,
    cor_matrix = 'signed', outlier_lo = -6, outlier_hi = 6,
    injection_type = 'null', feature_source = fs, feature_alt_id = alt,
    n_reps = N_A, notes = 'T1E null: directional two-sided, signed-K')
}

# ===========================================================
# GROUP B: T1E on the full c2 gmt (set-based)
# ===========================================================
N_B <- 500L
rows[[length(rows)+1L]] <- tmpl(scenario_id = 't1e_gmt_mag',
  arm = 'mag_probit_absK', mode = 'probit', two_sided = FALSE,
  cor_matrix = 'abs',    outlier_lo = -3, outlier_hi = 6,
  injection_type = 'null', feature_source = 'gmt_c2_full', feature_alt_id = 'Entrez',
  n_reps = N_B, notes = 'T1E null: gmt set-based, magnitude/probit')
rows[[length(rows)+1L]] <- tmpl(scenario_id = 't1e_gmt_dir1s_abs',
  arm = 'dir1s_absK',     mode = 'directional', two_sided = FALSE,
  cor_matrix = 'abs',    outlier_lo = -6, outlier_hi = 6,
  injection_type = 'null', feature_source = 'gmt_c2_full', feature_alt_id = 'Entrez',
  n_reps = N_B, notes = 'T1E null: gmt set-based, directional one-sided, abs-K')
rows[[length(rows)+1L]] <- tmpl(scenario_id = 't1e_gmt_dir2s_signed',
  arm = 'dir2s_signedK',  mode = 'directional', two_sided = TRUE,
  cor_matrix = 'signed', outlier_lo = -6, outlier_hi = 6,
  injection_type = 'null', feature_source = 'gmt_c2_full', feature_alt_id = 'Entrez',
  n_reps = N_B, notes = 'T1E null: gmt set-based, directional two-sided, signed-K')

# ===========================================================
# GROUP C: Directional headline
#   For each (form × mode × K × beta): 300 reps.
#   Fixed signature seed 111, sparsity 0.05, sigma_obs 1.0.
# ===========================================================
N_C <- 300L
sig_seed <- 111L
sparsity <- 0.05
sigma_obs <- 1.0

forms <- c('sig_cont', 'sig_tern', 'sig_bin')
modes <- list(
  dir2s     = list(mode = 'directional', two_sided = TRUE,  ohl = c(-6, 6)),
  dir1s     = list(mode = 'directional', two_sided = FALSE, ohl = c(-6, 6)),
  magnitude = list(mode = 'probit',      two_sided = FALSE, ohl = c(-3, 6))
)
Ks <- c('abs', 'signed')
betas <- c(-0.4, -0.2, 0, 0.2, 0.4)

for(form in forms){
  for(mn in names(modes)){
    m <- modes[[mn]]
    for(kk in Ks){
      # sig_bin is unsigned membership -> only meaningful under magnitude or
      # directional-two-sided (unsigned features can't drive directional 1-sided
      # correctly). Keep dir1s row but expect near-null power.
      for(bt in betas){
        scen_id <- sprintf('dh_%s_%s_%s_b%s', form, mn, kk, gsub('-', 'm', gsub('\\.', 'p', sprintf('%.2f', bt))))
        rows[[length(rows)+1L]] <- tmpl(
          scenario_id    = scen_id,
          arm            = sprintf('%s_%sK', mn, kk),
          mode           = m$mode,
          two_sided      = m$two_sided,
          cor_matrix     = kk,
          outlier_lo     = m$ohl[1], outlier_hi = m$ohl[2],
          injection_type = if(bt == 0) 'null' else 'property_random',
          feature_source = form,
          feature_alt_id = 'AltID_ENSG',
          sparsity       = sparsity,
          sigma_obs      = sigma_obs,
          beta           = bt,
          signature_seed = sig_seed,
          n_reps         = N_C,
          notes          = sprintf('DH form=%s mode=%s K=%s beta=%.2f', form, mn, kk, bt))
      }
    }
  }
}

scenarios <- rbindlist(rows)
fwrite(scenarios, 'sim/scenarios.tsv', sep = '\t')
cat('wrote sim/scenarios.tsv with', nrow(scenarios), 'scenarios.\n')

# Rough wall-clock projection using smoke timings.
per_rep_secs <- list(
  cmap_a375_full = 3.5,
  ternary_full   = 2.4,
  gmt_c2_full    = 9.8,
  sig_cont       = 1.0,
  sig_tern       = 1.0,
  sig_bin        = 1.0)
scenarios[, per_rep := unlist(per_rep_secs)[as.character(feature_source)]]
tot <- sum(scenarios$per_rep * scenarios$n_reps)
cat(sprintf('Total task-seconds (sequential): %.0f  |  @10-way parallelism: %.1f min\n',
            tot, tot / (10 * 60)))
cat(sprintf('T1E-real:      %5d tasks   ~%.1f min @10-way\n',
            sum(scenarios[grepl('^t1e_', scenario_id) & feature_source != 'gmt_c2_full']$n_reps),
            sum(scenarios[grepl('^t1e_', scenario_id) & feature_source != 'gmt_c2_full']$per_rep *
                scenarios[grepl('^t1e_', scenario_id) & feature_source != 'gmt_c2_full']$n_reps) / (10 * 60)))
cat(sprintf('T1E-gmt:       %5d tasks   ~%.1f min @10-way\n',
            sum(scenarios[grepl('^t1e_gmt_', scenario_id)]$n_reps),
            sum(scenarios[grepl('^t1e_gmt_', scenario_id)]$per_rep *
                scenarios[grepl('^t1e_gmt_', scenario_id)]$n_reps) / (10 * 60)))
cat(sprintf('Dir headline:  %5d tasks   ~%.1f min @10-way\n',
            sum(scenarios[grepl('^dh_', scenario_id)]$n_reps),
            sum(scenarios[grepl('^dh_', scenario_id)]$per_rep *
                scenarios[grepl('^dh_', scenario_id)]$n_reps) / (10 * 60)))
