#!/usr/bin/env bash
# run_validation.sh
# Validate --fast_competitive T against F on the bundled mini fixture data.
#
# Usage (from repo root):
#   bash validation/run_validation.sh
#
# Requires: miniforge3 module available, twas_gsea conda environment installed.
# Runtime: ~2 minutes on a single core.

set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
OUTDIR="$SCRIPT_DIR/output"
mkdir -p "$OUTDIR"

# ---- Conda setup ------------------------------------------------------------
# conda activate scripts may reference unbound variables; disable -u for this block
set +u
module load miniforge3/24.1.2-0-gcc-13.2.0
eval "$(conda shell.bash hook)"
conda activate twas_gsea
set -u
export R_LIBS_USER=''

RSCRIPT="Rscript"
TWAS_SCRIPT="$REPO_ROOT/TWAS-GSEA.V1.2.R"

COMMON_ARGS=(
  --twas_results "$REPO_ROOT/ukbiobank-2017-1160-prePRS-fusion-mini.tsv.GW"
  --gmt_file     "$REPO_ROOT/c2.all.v7.5.1.mini.entrez.gmt"
  --expression_ref "$REPO_ROOT/CMC.BRAIN.RNASEQ_GeneX_all_MINI.csv"
  --pos          "$REPO_ROOT/CMC.BRAIN.RNASEQ.pos"
  --competitive  T
  --self_contained F
  --linear_p_thresh 1
  --n_cores 1
  --qqplot F
)

# ---- Run fast mode ----------------------------------------------------------
echo "=== Run 1: --fast_competitive T ==="
START_FAST=$SECONDS
$RSCRIPT "$TWAS_SCRIPT" "${COMMON_ARGS[@]}" \
  --fast_competitive T \
  --output "$OUTDIR/val_fast" \
  2>&1 | tee "$OUTDIR/val_fast.log"
ELAPSED_FAST=$(( SECONDS - START_FAST ))
echo "Fast mode wall time: ${ELAPSED_FAST}s"

# ---- Run original mode ------------------------------------------------------
echo ""
echo "=== Run 2: --fast_competitive F ==="
START_ORIG=$SECONDS
$RSCRIPT "$TWAS_SCRIPT" "${COMMON_ARGS[@]}" \
  --fast_competitive F \
  --output "$OUTDIR/val_orig" \
  2>&1 | tee "$OUTDIR/val_orig.log"
ELAPSED_ORIG=$(( SECONDS - START_ORIG ))
echo "Original mode wall time: ${ELAPSED_ORIG}s"

# ---- Benchmark summary ------------------------------------------------------
echo ""
echo "=== Benchmark ==="
echo "  fast_competitive T: ${ELAPSED_FAST}s"
echo "  fast_competitive F: ${ELAPSED_ORIG}s"
if [ "$ELAPSED_FAST" -gt 0 ]; then
  python3 -c "print(f'  Speedup (total):    {$ELAPSED_ORIG / $ELAPSED_FAST:.2f}x')" 2>/dev/null || \
    echo "  Speedup: $ELAPSED_ORIG / $ELAPSED_FAST (total run)"
fi

# ---- Compare outputs --------------------------------------------------------
echo ""
echo "=== Concordance ==="
$RSCRIPT "$SCRIPT_DIR/compare_competitive.R" \
  "$OUTDIR/val_fast.competitive.txt" \
  "$OUTDIR/val_orig.competitive.txt"
