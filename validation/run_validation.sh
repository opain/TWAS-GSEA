#!/usr/bin/env bash
# run_validation.sh
# Validate --fast_competitive T against F and benchmark runtime.
#
# Usage (from repo root):
#   bash validation/run_validation.sh
#
# Requires: miniforge3 module, twas_gsea conda environment.
# Runtime: ~5 minutes on a single core.
#
# Two stages:
#   Stage 1 — Concordance: run both modes on the real 21-gene-set fixture.
#   Stage 2 — Benchmark: pre-compute CorMat once, run both modes on
#             validation/bench.gmt (210 gene sets) so the competitive step
#             dominates total runtime and timing is not confounded by noise.

set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
OUTDIR="$SCRIPT_DIR/output"
mkdir -p "$OUTDIR"

# ---- Conda setup ------------------------------------------------------------
set +u
module load miniforge3/24.1.2-0-gcc-13.2.0
eval "$(conda shell.bash hook)"
conda activate twas_gsea
set -u
export R_LIBS_USER=''
RSCRIPT="Rscript"
TWAS="$REPO_ROOT/TWAS-GSEA.V1.2.R"

BASE_ARGS=(
  --twas_results "$REPO_ROOT/ukbiobank-2017-1160-prePRS-fusion-mini.tsv.GW"
  --pos          "$REPO_ROOT/CMC.BRAIN.RNASEQ.pos"
  --competitive  T
  --self_contained F
  --linear_p_thresh 1
  --n_cores 1
  --qqplot F
)

# =============================================================================
# Stage 1: Concordance on real 21-gene-set fixture
# =============================================================================
echo "=== Stage 1: Concordance (21 gene sets) ==="

echo "  Run 1: --fast_competitive T"
$RSCRIPT "$TWAS" "${BASE_ARGS[@]}" \
  --gmt_file        "$REPO_ROOT/c2.all.v7.5.1.mini.entrez.gmt" \
  --expression_ref  "$REPO_ROOT/CMC.BRAIN.RNASEQ_GeneX_all_MINI.csv" \
  --fast_competitive T \
  --output "$OUTDIR/val_fast" \
  > "$OUTDIR/val_fast.stdout" 2>&1

echo "  Run 2: --fast_competitive F"
$RSCRIPT "$TWAS" "${BASE_ARGS[@]}" \
  --gmt_file        "$REPO_ROOT/c2.all.v7.5.1.mini.entrez.gmt" \
  --expression_ref  "$REPO_ROOT/CMC.BRAIN.RNASEQ_GeneX_all_MINI.csv" \
  --fast_competitive F \
  --output "$OUTDIR/val_orig" \
  > "$OUTDIR/val_orig.stdout" 2>&1

echo "  Comparing outputs..."
$RSCRIPT "$SCRIPT_DIR/compare_competitive.R" \
  "$OUTDIR/val_fast.competitive.txt" \
  "$OUTDIR/val_orig.competitive.txt" 2>/dev/null

# =============================================================================
# Stage 2: Benchmark on 210-gene-set fixture with pre-computed CorMat
# =============================================================================
echo ""
echo "=== Stage 2: Benchmark (210 gene sets, pre-computed CorMat) ==="
echo "  (bench.gmt contains 21 real gene sets replicated 10x)"

# Pre-compute the correlation matrix if not already present
if [ ! -f "$OUTDIR/bench_precomp.CorMat.RDS" ]; then
  echo "  Pre-computing correlation matrix (one-time, ~30s)..."
  $RSCRIPT "$TWAS" "${BASE_ARGS[@]}" \
    --gmt_file       "$SCRIPT_DIR/bench.gmt" \
    --expression_ref "$REPO_ROOT/CMC.BRAIN.RNASEQ_GeneX_all_MINI.csv" \
    --fast_competitive T \
    --save_CorMat T \
    --output "$OUTDIR/bench_precomp" \
    > "$OUTDIR/bench_precomp.stdout" 2>&1
fi

BENCH_ARGS=(
  --gmt_file     "$SCRIPT_DIR/bench.gmt"
  --input_CorMat "$OUTDIR/bench_precomp.CorMat.RDS"
)

echo "  Run 1: --fast_competitive T"
{ time $RSCRIPT "$TWAS" "${BASE_ARGS[@]}" "${BENCH_ARGS[@]}" \
    --fast_competitive T \
    --output "$OUTDIR/bench_fast" \
    > "$OUTDIR/bench_fast.stdout" 2>&1
} 2>&1 | grep -E "^real"

echo "  Run 2: --fast_competitive F"
{ time $RSCRIPT "$TWAS" "${BASE_ARGS[@]}" "${BENCH_ARGS[@]}" \
    --fast_competitive F \
    --output "$OUTDIR/bench_orig" \
    > "$OUTDIR/bench_orig.stdout" 2>&1
} 2>&1 | grep -E "^real"

echo ""
FAST_DUR=$(grep "Analysis duration" "$OUTDIR/bench_fast.log" 2>/dev/null | tail -1)
ORIG_DUR=$(grep "Analysis duration" "$OUTDIR/bench_orig.log" 2>/dev/null | tail -1)
echo "  $FAST_DUR  [fast]"
echo "  $ORIG_DUR  [orig]"
echo "  Difference = competitive-step cost at 210 gene sets"

echo ""
echo "  Concordance on 210-set outputs:"
$RSCRIPT "$SCRIPT_DIR/compare_competitive.R" \
  "$OUTDIR/bench_fast.competitive.txt" \
  "$OUTDIR/bench_orig.competitive.txt" 2>/dev/null
