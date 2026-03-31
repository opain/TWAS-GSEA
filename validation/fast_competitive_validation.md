# Validation: `--fast_competitive T` vs `--fast_competitive F`

## What was compared

The competitive mixed-model analysis has two paths:

| Flag | Method |
|------|--------|
| `--fast_competitive T` (default) | Fixed-V GLS whitening: builds V = σ²_u·K + σ²_e·I from null-model estimates, computes sparse Cholesky L, whitens y and all gene-set indicators simultaneously, evaluates Wald statistics in one vectorised pass |
| `--fast_competitive F` | Original per-gene-set path: for each gene set, calls `merPredD` to add a fixed-effect column then `refit` to re-optimise variance components (σ²_u, σ²_e) under the alternative model before computing the Wald test |

The fast path is an **approximation** because variance components are fixed at null
estimates rather than re-optimised per gene set. The approximation is expected to
be tight when N >> gene-set size; it has only been validated on the single test
dataset described below.

## Dataset

| Quantity | Value |
|----------|-------|
| TWAS features in input | 561 |
| Features after filtering (unique, non-missing) | 539 |
| Features matched to expression reference | 509 |
| Expression reference individuals | 9 |
| Genomic blocks | 3 |
| Correlation matrix sparsity | 60.8 % |

## Concordance (21 real gene sets)

Run on `c2.all.v7.5.1.mini.entrez.gmt` (21 gene sets passing `--min_Ngenes 5`).
See `compare_competitive.R` for the reproducible calculation.

| Metric | Value |
|--------|-------|
| Pearson r (T-statistic) | 0.999996 |
| Spearman ρ (rank order) | 0.99870 |
| Max \|ΔT\| | 0.014 |
| Max \|ΔP\| | 0.0020 |
| Sign agreement | 21 / 21 (100 %) |
| Top gene set same? | Yes (DAVICIONI.MOLECULAR.ARMS.VS.ERMS.UP) |
| Rank changes | 2 / 21, both with \|ΔT\| < 0.005 (effectively tied estimates) |

## Benchmark (210 gene sets, pre-computed correlation matrix)

To isolate the competitive step from preprocessing noise, the correlation
matrix was pre-computed once (`--save_CorMat T`) and both modes were run
with `--input_CorMat` on a 210-gene-set fixture (`validation/bench.gmt`,
the 21 real sets replicated 10× with distinct names).

| Step | fast_competitive T | fast_competitive F |
|------|-------------------|--------------------|
| Full run (incl. preprocessing) | 63 s | 201 s |
| Preprocessing (shared, estimated as fast-mode total) | ~63 s | ~63 s |
| Competitive step alone (difference) | ~0 s | ~138 s |
| Per gene set (competitive loop) | O(1) total | ~0.66 s / gene set |
| Speedup (full run) | 1× | 0.31× (3.2× slower) |

The 138-second saving comes entirely from eliminating the per-gene-set REML
re-optimisation. The fast path does one sparse Cholesky and a batched matrix
solve regardless of the number of gene sets.

For reference, the total-run comparison on the 21-gene-set fixture was noisy
(observations: fast=41 s / orig=48 s; fast=73 s / orig=62 s) and should not
be used to draw conclusions about speed.

## Conclusion

**Concordance**: On this dataset, the fast and original paths are indistinguishable
in practice — Pearson r = 0.999996, max |ΔT| = 0.014, 100% sign agreement,
identical top hit. The only rank changes (2/21) are between near-tied estimates
(|ΔT| < 0.005).

**Speed**: The saving is large and scales linearly with number of gene sets:
~138 s for 210 gene sets, ~0.66 s per gene set eliminated. For a typical analysis
with 300–500 gene sets the competitive step would take ~3–5 minutes less.

**Caveats**:
- Validated on one real TWAS dataset (UKBiobank sleep duration, N=561 genes)
  and one expression reference (CMC brain, N=9 individuals). Approximation
  quality may differ for datasets with very large gene sets, higher LD structure,
  or more variable per-set θ estimates.
- The benchmark GMT (`validation/bench.gmt`) is synthetic (replicated sets),
  so concordance metrics on 210 sets are not independent of those on 21.
- No validation with covariates (`--covar`) or weights (`--weights`).

**Recommendation**: `--fast_competitive T` is reasonable as the default for
typical TWAS-GSEA usage. `--fast_competitive F` should be retained as a
sensitivity-check option.
