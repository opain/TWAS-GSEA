#!/usr/bin/Rscript
# sim/plots.R
#
# Consumes aggregate.R's outputs + raw long.rds files to build the report
# figures: T1E QQ, power vs beta by arm, discretisation penalty, abs-K vs
# signed-K comparisons, sign recovery.

suppressMessages({
  library(optparse)
  library(data.table)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option('--grid_dir',   default = 'sim/output/grid',      type = 'character'),
  make_option('--agg_dir',    default = 'sim/output/aggregate', type = 'character'),
  make_option('--scenarios',  default = 'sim/scenarios.tsv',    type = 'character'),
  make_option('--outdir',     default = 'sim/output/plots',     type = 'character')
)))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
scen <- fread(opt$scenarios)
scen_sum <- fread(file.path(opt$agg_dir, 'scenario_summary.tsv'))
pow      <- fread(file.path(opt$agg_dir, 'power_summary.tsv'))

# ---------- Helpers ----------
qq_uniform <- function(pvals, main, sub = '') {
  p <- pvals[!is.na(pvals) & pvals > 0 & pvals < 1]
  if(length(p) < 10){
    plot.new(); title(main = main, sub = 'n < 10, skipped'); return(invisible())
  }
  n <- length(p)
  x <- -log10(seq_len(n) / (n + 1))
  y <- -log10(sort(p, decreasing = TRUE))
  # Use a subsample if huge
  if(length(x) > 20000){
    idx <- sample(seq_along(x), 20000)
    x <- x[idx]; y <- y[idx]
  }
  plot(x, y, pch = '.', cex = 1.5,
       xlab = '-log10 expected', ylab = '-log10 observed',
       main = main, sub = sub)
  abline(0, 1, col = 'red', lty = 2)
}

# ---------- Figure 1: T1E QQ per arm (real prop panels + gmt) ----------
t1e_scenarios <- scen_sum[injection_type == 'null']
if(nrow(t1e_scenarios) > 0){
  n_scen <- nrow(t1e_scenarios)
  ncol <- 3
  nrow <- ceiling(n_scen / ncol)
  pdf(file.path(opt$outdir, 'fig1_t1e_qq.pdf'), width = ncol * 4, height = nrow * 3.5)
  par(mfrow = c(nrow, ncol), mar = c(3, 3, 2.5, 0.6), mgp = c(1.8, 0.6, 0), cex = 0.7)
  for(i in seq_len(n_scen)){
    sid <- t1e_scenarios$scenario_id[i]
    f <- file.path(opt$grid_dir, paste0(sid, '.long.rds'))
    if(!file.exists(f)) next
    x <- as.data.table(readRDS(f))
    qq_uniform(x$P,
      main = sprintf('%s\nfrac<.05=%.3f  FWER=%.3f',
                     sid, t1e_scenarios$frac_p_lt_05[i], t1e_scenarios$fwer_at_05[i]))
  }
  dev.off()
  cat('wrote fig1_t1e_qq.pdf\n')
}

# ---------- Figure 2: Directional headline power vs beta ----------
if(nrow(pow) > 0){
  pdf(file.path(opt$outdir, 'fig2_dh_power_vs_beta.pdf'), width = 12, height = 8)
  par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
  # x = beta,  y = power_p05.  One panel per (feature_source, cor_matrix).
  fs_levels <- c('sig_cont', 'sig_tern', 'sig_bin')
  K_levels  <- c('abs', 'signed')
  arm_colors <- c(dir2s_absK = 'black', dir2s_signedK = 'blue',
                  dir1s_absK = 'purple', dir1s_signedK = 'orange',
                  magnitude_absK = 'red', magnitude_signedK = 'darkred')
  for(K in K_levels){
    for(fs in fs_levels){
      sub <- pow[feature_source == fs & cor_matrix == K]
      if(nrow(sub) == 0){ plot.new(); title(paste(fs, K)); next }
      plot(NA, xlim = range(sub$beta),
           ylim = c(0, 1),
           xlab = 'beta', ylab = 'power (frac P<0.05)',
           main = sprintf('%s  K=%s', fs, K))
      abline(h = 0.05, col = 'grey', lty = 2)
      arms <- unique(sub$arm)
      for(a in arms){
        s2 <- sub[arm == a][order(beta)]
        col <- arm_colors[a]; if(is.na(col)) col <- 'darkgreen'
        lines(s2$beta, s2$power_p05, type = 'b', pch = 19, col = col, lwd = 2)
      }
      legend('topleft', legend = arms, col = arm_colors[arms], lwd = 2, pch = 19,
             cex = 0.7, bty = 'n')
    }
  }
  dev.off()
  cat('wrote fig2_dh_power_vs_beta.pdf\n')

  # ---------- Figure 3: Sign recovery (fraction significant AND correct sign) ----------
  pdf(file.path(opt$outdir, 'fig3_dh_sign_recovery.pdf'), width = 12, height = 6)
  par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
  for(fs in fs_levels){
    sub <- pow[feature_source == fs & !is.na(sign_recovery)]
    if(nrow(sub) == 0){ plot.new(); title(fs); next }
    plot(NA, xlim = range(sub$beta), ylim = c(0, 1),
         xlab = 'beta', ylab = 'sign recovery (frac correct-sign & P.CORR<0.05)',
         main = sprintf('%s  sign recovery', fs))
    keys <- unique(sub[, .(arm, cor_matrix)])
    for(kk in seq_len(nrow(keys))){
      s2 <- sub[arm == keys$arm[kk] & cor_matrix == keys$cor_matrix[kk]][order(beta)]
      col <- if(keys$cor_matrix[kk] == 'signed') 'blue' else 'black'
      lty <- if(grepl('^dir', keys$arm[kk])) 1 else 2
      lines(s2$beta, s2$sign_recovery, col = col, lty = lty, lwd = 2, pch = 19, type = 'b')
    }
    legend('topleft', c('directional (solid)', 'magnitude (dashed)',
                        'abs K (black)', 'signed K (blue)'),
           cex = 0.7, bty = 'n')
  }
  dev.off()
  cat('wrote fig3_dh_sign_recovery.pdf\n')

  # ---------- Figure 4: Discretisation penalty (cont - tern power at matched beta) ----------
  pdf(file.path(opt$outdir, 'fig4_discretisation_penalty.pdf'), width = 9, height = 6)
  pen <- merge(
    pow[feature_source == 'sig_cont', .(arm, cor_matrix, beta, power_cont = power_p05)],
    pow[feature_source == 'sig_tern', .(arm, cor_matrix, beta, power_tern = power_p05)],
    by = c('arm', 'cor_matrix', 'beta'))
  pen[, penalty := power_cont - power_tern]
  par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
  plot(NA, xlim = range(pen$beta), ylim = range(pen$penalty, na.rm = TRUE),
       xlab = 'beta', ylab = 'power(cont) - power(tern)',
       main = 'Discretisation penalty (higher = ternary hurts more)')
  abline(h = 0, col = 'grey', lty = 2)
  arms <- unique(pen[, .(arm, cor_matrix)])
  for(i in seq_len(nrow(arms))){
    s2 <- pen[arm == arms$arm[i] & cor_matrix == arms$cor_matrix[i]][order(beta)]
    col <- if(arms$cor_matrix[i] == 'signed') 'blue' else 'black'
    lty <- if(grepl('^dir', arms$arm[i])) 1 else 2
    lines(s2$beta, s2$penalty, col = col, lty = lty, lwd = 2, pch = 19, type = 'b')
  }
  legend('topleft', c('directional (solid)', 'magnitude (dashed)',
                      'abs K (black)', 'signed K (blue)'),
         cex = 0.8, bty = 'n')
  dev.off()
  cat('wrote fig4_discretisation_penalty.pdf\n')

  # ---------- Figure 5: abs-K vs signed-K power delta ----------
  pdf(file.path(opt$outdir, 'fig5_absK_vs_signedK.pdf'), width = 12, height = 6)
  par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
  for(fs in fs_levels){
    ab <- pow[feature_source == fs & cor_matrix == 'abs',
              .(arm, beta, power_abs = power_p05)]
    sg <- pow[feature_source == fs & cor_matrix == 'signed',
              .(arm, beta, power_signed = power_p05)]
    m <- merge(ab, sg, by = c('arm', 'beta'))
    if(nrow(m) == 0){ plot.new(); title(fs); next }
    plot(NA, xlim = range(m$beta), ylim = c(-0.3, 0.3),
         xlab = 'beta', ylab = 'power(signed K) - power(abs K)',
         main = sprintf('%s: signed vs abs K', fs))
    abline(h = 0, col = 'grey', lty = 2)
    arms <- unique(m$arm)
    for(a in arms){
      s2 <- m[arm == a][order(beta)]
      lines(s2$beta, s2$power_signed - s2$power_abs, type = 'b', pch = 19, lwd = 2)
      text(max(s2$beta), tail(s2$power_signed - s2$power_abs, 1), a, pos = 4, cex = 0.7)
    }
  }
  dev.off()
  cat('wrote fig5_absK_vs_signedK.pdf\n')
}

cat('\nplots.R done.\n')
