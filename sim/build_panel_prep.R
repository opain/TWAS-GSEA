#!/usr/bin/Rscript
# sim/build_panel_prep.R
#
# One-shot: run panel_prep() on the signed-K matrix and cache the per-block
# chol/eigen decompositions. Downstream sim and in-session GLS reuse the
# result without repeating the eigen.

suppressMessages(library(optparse))

opt <- parse_args(OptionParser(option_list = list(
  make_option('--cor_signed', default = NA, type = 'character',
              help = 'Path to .CorMatSigned.RDS from sim/build_signed_cor.R'),
  make_option('--output', default = NA, type = 'character',
              help = 'Output path for the .prep.RDS cache')
)))
stopifnot(!is.na(opt$cor_signed), !is.na(opt$output))

SCRIPT_DIR <- dirname(sub('--file=', '', grep('--file=', commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
source(file.path(SCRIPT_DIR, 'R', 'sim_generator.R'))

prep <- panel_prep(opt$cor_signed, cache_path = opt$output)
cat('Wrote ', opt$output, '\n', sep = '')
