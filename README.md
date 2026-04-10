# TWAS-based Gene Set Enrichment Analysis (TWAS-GSEA)

TWAS-GSEA is a tool for performing gene set or gene property analysis on TWAS results. It uses a similar method to MAGMA, in that it fits a mixed model to test for enrichment, specifying a gene-gene correlation matrix as a random effect to avoid bias due to non-independent gene-level statistics.

The recommended workflow is `build_cor_matrix.R` + `TWAS-GSEA-fast.R`. The original `TWAS-GSEA.V1.2.R` is **deprecated** for competitive analyses but is retained for linear, self-contained and weighted analyses (see [Legacy](#legacy-twas-geseav12r) below).

TWAS-GSEA was written to analyse the output of FUSION's [**FUSION.assoc_test.R**](https://github.com/gusevlab/fusion_twas/blob/master/FUSION.assoc_test.R), though it can be used on any gene-level association results.

## Install

You can use conda/mamba to install the dependencies, or install them manually.

### Using conda

```bash
conda env create -f env.yaml
conda activate twas_gsea
```

`lme4qtl` is only required if you want to run the legacy `TWAS-GSEA.V1.2.R`. The recommended `TWAS-GSEA-fast.R` workflow does not depend on it.

```R
# Optional, only for legacy V1.2:
library(devtools)
install_github("variani/lme4qtl")
```

### Manual install

```R
install.packages(c('data.table','optparse','WGCNA','Matrix','VGAM','foreach','doMC','matrixcalc','qusage'))
# Optional, legacy V1.2 only:
install.packages(c('gdata','lme4','pbkrtest','speedglm'))
library(devtools); install_github("variani/lme4qtl")
```

### Required upstream tools

- Perform a TWAS using FUSION: instructions [here](http://gusevlab.org/projects/fusion/).
- Predict gene expression in a reference panel: instructions [here](https://github.com/opain/Predicting-TWAS-features).

---

## Recommended workflow

Two steps. Step 1 is per expression panel and is reused across many TWAS / gene-set runs. Step 2 is the per-analysis call.

### Step 1: precompute the gene-gene correlation matrix

```sh
Rscript build_cor_matrix.R \
  --expression_ref CMC.BRAIN.RNASEQ_GeneX_all_MINI.csv \
  --pos            CMC.BRAIN.RNASEQ.pos \
  --output         demo
```

Produces `demo.CorMat.RDS` (a sparse block-diagonal correlation matrix plus per-row block index) and `demo.CorMat.log`. Run this once per expression panel; reuse the `.RDS` for every subsequent analysis on the same panel.

#### `build_cor_matrix.R` options

| Flag | Default | Description |
|---|---|---|
| `--expression_ref` | required | Predicted-expression file (FeaturePred output). One column per gene. `.gz` allowed. |
| `--pos` | required | FUSION `.pos` file with `WGT`, `CHR`, `P0`, `P1`. |
| `--cor_window` | `5e6` | Window in bp for retaining gene-gene correlations. |
| `--min_r2` | `1e-4` | r² threshold below which correlations are zeroed (sparsification). |
| `--max_r2` | `1` | r² threshold above which one of a collinear gene pair is dropped. |
| `--n_cores` | `1` | Cores for the per-block correlation step. |
| `--output` | required | Output prefix. Writes `<output>.CorMat.RDS` and `<output>.CorMat.log`. |

### Step 2: run the analysis

```sh
Rscript TWAS-GSEA-fast.R \
  --twas_results ukbiobank-2017-1160-prePRS-fusion-mini.tsv.GW \
  --input_CorMat demo.CorMat.RDS \
  --gmt_file     c2.all.v7.5.1.mini.entrez.gmt \
  --gene_id_map  data/gene_id_map.tsv \
  --output       demo
```

Produces `demo.competitive.txt` and `demo.log`.

#### `TWAS-GSEA-fast.R` options

| Flag | Default | Description |
|---|---|---|
| `--twas_results` | required | TWAS results file (FUSION format or any file with `FILE`, `ID`, `TWAS.Z`, `TWAS.P`, `MODELCV.R2`). `.gz` allowed. |
| `--input_CorMat` | required | `.CorMat.RDS` from `build_cor_matrix.R`. |
| `--gmt_file` | NA | Gene-set file in `.gmt` format. Mutually exclusive with `--prop_file`. |
| `--prop_file` | NA | Gene-property file (first column `ID`, then one numeric column per property; `.rds` files containing a numeric matrix with row names are also accepted). Mutually exclusive with `--gmt_file`. |
| `--covar` | `none` | Comma-separated covariate columns from `--twas_results`, e.g. `NSNP,GeneLength`. |
| `--use_alt_id` | NA | Alternative ID column in `--twas_results` to match against the gene-set / property file (e.g. `ID`). When unset, gene symbols → Entrez IDs are mapped via `--gene_id_map`. |
| `--gene_id_map` | NA | Tab-delimited file with columns `external_gene_name`, `entrezgene_id`. Required when `--use_alt_id` is not set. A bundled GRCh37 copy is at [`data/gene_id_map.tsv`](data/gene_id_map.tsv). |
| `--directional` | `FALSE` | If `T`, use signed `TWAS.Z` as outcome (tests for direction-aware enrichment). |
| `--probit_P_as_Z` | `TRUE` | Use `probit(1 − TWAS.P)` as outcome (overridden when `--directional T`). |
| `--two_sided` | `FALSE` | Two-sided p-values (test enrichment **and** depletion). See note below. |
| `--outlier_threshold` | `-3,6` | Lower,upper truncation of the outcome z-score. |
| `--min_Ngenes` | `2` | Minimum gene-set / property size to test. |
| `--allow_duplicate_ID` | `FALSE` | Keep multiple TWAS rows per gene; otherwise retain best `MODELCV.R2`. |
| `--h_max` | `100` | Upper bound on the variance ratio h = σ²_u / σ²_e in the REML search. |
| `--reml_tol` | `1e-6` | Convergence tolerance for the 1-D Brent search. |
| `--twas_p_thresh` | `1` | Comma-separated TWAS p-value thresholds. The analysis is run separately on each subset of genes with `TWAS.P ≤ threshold`. `1` retains all genes (default behaviour). |
| `--p_cor_method` | `fdr` | Multiple-testing correction (passed to `p.adjust`). |
| `--n_cores` | `1` | Cores. |
| `--output` | required | Output prefix. Writes `<output>.competitive.txt` and `<output>.log`. When multiple `--twas_p_thresh` values are given, non-1 thresholds write `<output>.pT<value>.competitive.txt`. |

---

## How it works

`TWAS-GSEA-fast.R` fits the same statistical model as V1.2:

```
y = X β + g + ε     g ~ N(0, σ²_u K)     ε ~ N(0, σ²_e I)
V = σ²_u K + σ²_e I
```

where `K` is the block-diagonal sparse predicted-expression correlation matrix from `build_cor_matrix.R`.

V1.2 estimates `(σ²_u, σ²_e)` via `lme4qtl::relmatLmer`, which is a general sparse-mixed-model fitter and does not exploit `K`'s block-diagonal structure. `TWAS-GSEA-fast.R` instead uses a custom REML routine (`R/reml_blockdiag.R`) that:

1. Computes a per-block eigendecomposition of `K` once upfront. Cost is `Σ_b n_b³`, sub-second for typical TWAS panels (~50 blocks of size ~150).
2. Reparameterises `V = σ²(hK + I)` and profiles out `σ²` analytically.
3. Optimises the resulting 1-D profiled likelihood in `h` via Brent's method (`stats::optimise`).

This is the FastLMM / GCTA-LMM trick specialised to the block-diagonal case. It is **mathematically equivalent** to `lme4qtl::relmatLmer` on the same model — same MLE, not an approximation — but much faster and with no `lme4qtl` dependency on the hot path.

After the null fit, `TWAS-GSEA-fast.R` runs the same vectorised GLS pipeline as V1.2's `--fast_competitive T` path: sparse Cholesky of `V_hat`, whiten `y` / `X` / gene-set indicator matrix `Z`, residualise against the null design via QR, then a single batched `crossprod` across all gene sets.

**Concordance and speed** (from `validation/run_validation.sh` Stage 3, on the bench fixture of 210 gene sets / ~539 genes):

| | Wall time | Concordance vs V1.2 slow |
|---|---|---|
| V1.2 slow (`--fast_competitive F`) | 3m 14s | reference |
| V1.2 fast (`--fast_competitive T`) | 15.6s | r(T) = 0.999996, max\|ΔT\| = 0.014 |
| **TWAS-GSEA-fast.R** | **4.0s** | **r(T) = 0.999985, max\|ΔT\| = 0.018** |

The custom REML fit alone takes ~0.24s; the rest is I/O and the GLS pipeline.

---

## Input file formats

### `--twas_results`

The output of [**FUSION.assoc_test.R**](https://github.com/gusevlab/fusion_twas/blob/master/FUSION.assoc_test.R) or any file containing the columns `FILE`, `ID`, `P0`, `P1`, `TWAS.Z`, `TWAS.P`. Per-chromosome files should be combined into a single file. An example is at [ukbiobank-2017-1160-prePRS-fusion-mini.tsv.GW](ukbiobank-2017-1160-prePRS-fusion-mini.tsv.GW). Gene IDs are expected to be gene symbols (override with `--use_alt_id`). When `--allow_duplicate_ID F` (default), the file must also contain `MODELCV.R2` (used to retain the best feature per gene).

### `--pos`

A FUSION `.pos` file with `WGT`, `ID`, `CHR`, `P0`, `P1`. Used by `build_cor_matrix.R` to determine block boundaries.

### `--expression_ref`

The output of FeaturePred: first two columns `FID`/`IID`, then one column per gene. Column names must match the `FILE` column of `--twas_results` after stripping the path prefix and `.wgt.RDat` suffix. Whitespace- or comma-delimited; `.gz` allowed.

### `--gmt_file` (gene-set analysis)

Standard `.gmt` file: tab-delimited, one row per set, first column = set name, second = description (ignored), then a series of Entrez IDs. Example: [c2.all.v7.5.1.mini.entrez.gmt](c2.all.v7.5.1.mini.entrez.gmt).

### `--prop_file` (gene-property analysis)

Either a text file (first column header `ID` with Entrez IDs, then one numeric column per property) or an `.rds` file containing a numeric matrix with Entrez IDs as row names. Properties are z-scored within `TWAS-GSEA-fast.R` before testing.

### `--gene_id_map`

Two tab-delimited columns: `external_gene_name`, `entrezgene_id`. Used to map TWAS gene symbols to Entrez IDs when `--use_alt_id` is not set, avoiding a live BioMart query. A bundled GRCh37 copy is at [data/gene_id_map.tsv](data/gene_id_map.tsv).

---

## Output files

### `<output>.competitive.txt` / `<output>.pT<value>.competitive.txt`

Space-delimited results of the competitive test. When `--twas_p_thresh` includes multiple values, threshold 1 writes `<output>.competitive.txt` and other thresholds write `<output>.pT<value>.competitive.txt` (e.g. `<output>.pT0.05.competitive.txt`).

| Column | Meaning |
|---|---|
| `GeneSet` | Gene-set or property name |
| `Estimate` | GLS coefficient |
| `SE` | Standard error |
| `T` | Test statistic |
| `N_Mem_Avail` | Genes from the set with TWAS data (gmt mode only) |
| `N_Mem` | Total genes in the set (gmt mode only) |
| `P` | One-sided (or two-sided if `--two_sided T`) p-value |
| `P.CORR` | Multiple-testing-corrected p-value (`p.adjust` with `--p_cor_method`) |

### `<output>.log`

Run log with options, gene counts, the REML fit (`σ²_u`, `σ²_e`, `h_hat`), and timing.

---

## Statistical note: one-sided p-values

By default all p-values reported by `TWAS-GSEA-fast.R` (and V1.2) are **one-sided**, testing for positive enrichment (gene-set members have higher `ZSCORE` than background). This is designed to detect gene sets whose members show stronger TWAS associations than expected, but will not detect systematic depletion. If you need two-sided tests use `--two_sided T`, or convert manually via `p_two_sided = 2 * min(p, 1 - p)`.

---

## Legacy: TWAS-GSEA.V1.2.R

`TWAS-GSEA.V1.2.R` is **deprecated** for competitive analyses — `TWAS-GSEA-fast.R` is mathematically equivalent and ~50× faster. V1.2 prints a deprecation banner on every run. It remains available because it still provides several features `TWAS-GSEA-fast.R` does not:

- **Linear stage** — `*.linear.txt` from a fast standard linear model, plus `--linear_p_thresh` to gate which sets go to the mixed model.
- **Self-contained analysis** (`--self_contained T`) — `*.self_contained.txt`.
- **`--weights`** — variance-weighting via `lme4qtl`.
- **`*.sig.txt`** — gene-level breakdown for significant gene sets (linear and competitive).
- **`*.png`** QQ-plots (`--qqplot T`).
- **`--save_CorMat`** / `--input_CorMat` — V1.2 can save and reload its own correlation matrix (uses the same `R/build_cor_matrix_helper.R` as `build_cor_matrix.R`, so `.CorMat.RDS` files are interchangeable across the two scripts).
- Live BioMart fallback when `--gene_id_map` is not provided.

V1.2-only flags (in addition to those listed above): `--competitive`, `--self_contained`, `--linear_p_thresh`, `--qqplot`, `--save_CorMat`, `--input_CorMat`, `--fast_competitive`, `--weights`.

### V1.2 example

```sh
Rscript TWAS-GSEA.V1.2.R \
  --twas_results ukbiobank-2017-1160-prePRS-fusion-mini.tsv.GW \
  --pos          CMC.BRAIN.RNASEQ.pos \
  --gmt_file     c2.all.v7.5.1.mini.entrez.gmt \
  --expression_ref CMC.BRAIN.RNASEQ_GeneX_all_MINI.csv \
  --gene_id_map  data/gene_id_map.tsv \
  --linear_p_thresh 1 \
  --output       demo
```

For full V1.2 option documentation see the deprecation banner printed by the script and the inline `optparse` definitions at the top of `TWAS-GSEA.V1.2.R`.

---

## Help

If you have questions or comments, use the [Google group](https://groups.google.com/forum/#!forum/twas-related-r-scripts) or email oliver.pain@kcl.ac.uk.
