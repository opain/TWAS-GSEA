#!/usr/bin/Rscript
# sim/simulate_twas.R
#
# Thin CLI wrapper over sim/R/sim_generator.R. Emits ONE simulated TWAS TSV
# suitable for --twas_results in TWAS-GSEA-fast.R.
#
# Supported injection modes:
#   --injection null                          — mu = 0
#   --injection set --target_size N --Delta D [--rho R] [--pi_pos P]
#                                             — random causal set of size N
#   --injection property_random --sparsity S --beta B [--sd_signal 1]
#                                             — random signed signature
#
# For the fuller property mode (real drug signature), use the runner directly.

suppressMessages(library(optparse))

opt <- parse_args(OptionParser(option_list = list(
  make_option('--panel_prep',   default = NA, type = 'character', help = 'Path to <panel>.prep.RDS'),
  make_option('--template_twas', default = NA, type = 'character', help = 'Path to /data/twas_results (dir or single file)'),
  make_option('--ensg2entrez',  default = 'sim/output/id_map/ensg2entrez.tsv', type = 'character'),
  make_option('--biomart',      default = '/data/biomart/biomart_genes_grch37.tsv', type = 'character'),
  make_option('--annotate',     default = TRUE, type = 'logical', help = 'Add Symbol + Entrez columns'),
  make_option('--injection',    default = 'null', type = 'character'),
  make_option('--target_size',  default = 100L,  type = 'integer'),
  make_option('--Delta',        default = 2,     type = 'numeric'),
  make_option('--rho',          default = 1,     type = 'numeric'),
  make_option('--pi_pos',       default = 0.5,   type = 'numeric'),
  make_option('--sparsity',     default = 0.05,  type = 'numeric'),
  make_option('--beta',         default = 1,     type = 'numeric'),
  make_option('--sd_signal',    default = 1,     type = 'numeric'),
  make_option('--seed',         default = 1L,    type = 'integer'),
  make_option('--output',       default = NA,    type = 'character', help = 'Output TSV path')
)))
for(req in c('panel_prep','template_twas','output')){
  if(is.na(opt[[req]])) stop('--', req, ' is required')
}

SCRIPT_DIR <- dirname(sub('--file=', '', grep('--file=', commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
source(file.path(SCRIPT_DIR, 'R', 'sim_generator.R'))
source(file.path(SCRIPT_DIR, 'R', 'id_annotate.R'))

set.seed(opt$seed)

prep <- panel_prep(opt$panel_prep, cache_path = opt$panel_prep)
universe <- prep$gene_universe

tpl <- load_template_twas(opt$template_twas, universe = universe)
if(nrow(tpl) < length(universe)){
  missing <- setdiff(universe, tpl$FILE)
  stop(sprintf('template_twas is missing %d of %d universe genes (first: %s)',
               length(missing), length(universe), missing[1]))
}
tpl <- tpl[match(universe, tpl$FILE), , drop = FALSE]

mu <- switch(opt$injection,
  null = mu_null(prep),
  set = {
    target <- sample(universe, opt$target_size)
    mu_set(prep, target, opt$Delta, opt$rho, opt$pi_pos)
  },
  property_random = {
    N <- length(universe)
    n_nz <- max(1L, round(opt$sparsity * N))
    idx <- sample.int(N, n_nz)
    t_vec <- numeric(N)
    t_vec[idx] <- rnorm(n_nz, sd = opt$sd_signal)
    # Standardise over G
    t_vec <- (t_vec - mean(t_vec)) / sd(t_vec)
    names(t_vec) <- universe
    mu_property(prep, t_vec, opt$beta)
  },
  stop('Unknown --injection: ', opt$injection))

z <- simulate_z(prep, mu = mu)
out <- emit_twas(z, tpl, universe)

if(opt$annotate){
  out <- annotate_twas(out, ensg2entrez_path = opt$ensg2entrez, biomart_path = opt$biomart)
}

dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)
fwrite(out, opt$output, sep = '\t')
cat('Wrote', opt$output, '(', nrow(out), 'rows,', ncol(out), 'cols )\n')
