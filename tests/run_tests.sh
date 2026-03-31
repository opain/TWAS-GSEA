#!/usr/bin/env bash
# run_tests.sh — Regression test suite for TWAS-GSEA
#
# Runs the main script under several configurations against mini test data
# and compares outputs against saved reference files.
#
# Usage:
#   bash tests/run_tests.sh              # run tests (requires reference files)
#   bash tests/run_tests.sh --update     # generate/update reference files
#
# Requires: twas_gsea conda environment (see env.yaml).

set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
OUTDIR="$SCRIPT_DIR/output"
REFDIR="$SCRIPT_DIR/reference"
COMPARE="$SCRIPT_DIR/compare_outputs.R"
TWAS="$REPO_ROOT/TWAS-GSEA.V1.2.R"

UPDATE_MODE=false
if [[ "${1:-}" == "--update" ]]; then
  UPDATE_MODE=true
fi

mkdir -p "$OUTDIR"

# ---- Conda setup (same as validation) ----------------------------------------
set +u
if command -v module &>/dev/null; then
  module load miniforge3/24.1.2-0-gcc-13.2.0 2>/dev/null || true
fi
eval "$(conda shell.bash hook)"
conda activate twas_gsea
set -u
export R_LIBS_USER=''
RSCRIPT="Rscript"

# ---- Shared test data --------------------------------------------------------
TWAS_RESULTS="$REPO_ROOT/ukbiobank-2017-1160-prePRS-fusion-mini.tsv.GW"
POS_FILE="$REPO_ROOT/CMC.BRAIN.RNASEQ.pos"
GMT_FILE="$REPO_ROOT/c2.all.v7.5.1.mini.entrez.gmt"
EXPR_REF="$REPO_ROOT/CMC.BRAIN.RNASEQ_GeneX_all_MINI.csv"
GENE_ID_MAP="$REPO_ROOT/data/gene_id_map.tsv"

# Common args shared by all test configurations
BASE_ARGS=(
  --twas_results "$TWAS_RESULTS"
  --pos          "$POS_FILE"
  --n_cores      1
  --qqplot       F
  --gene_id_map  "$GENE_ID_MAP"
)

# ---- Test definitions --------------------------------------------------------
# Each test: name, extra args (as a string to eval), output suffixes to check.
# Output suffixes: "linear" always produced; "competitive" and "self_contained"
# depend on config.

declare -A TEST_ARGS
declare -A TEST_OUTPUTS

# Test 1: Linear only (no expression_ref → linear runs, mixed models skipped)
TEST_ARGS[linear_only]='--gmt_file "$GMT_FILE" --competitive T --self_contained F'
TEST_OUTPUTS[linear_only]="linear"

# Test 2: Competitive fast (default path)
TEST_ARGS[competitive_fast]='--gmt_file "$GMT_FILE" --expression_ref "$EXPR_REF" --competitive T --self_contained F --fast_competitive T --linear_p_thresh 1'
TEST_OUTPUTS[competitive_fast]="linear competitive"

# Test 3: Competitive original (merPredD+refit)
TEST_ARGS[competitive_orig]='--gmt_file "$GMT_FILE" --expression_ref "$EXPR_REF" --competitive T --self_contained F --fast_competitive F --linear_p_thresh 1'
TEST_OUTPUTS[competitive_orig]="linear competitive"

# Test 4: Self-contained
TEST_ARGS[self_contained]='--gmt_file "$GMT_FILE" --expression_ref "$EXPR_REF" --competitive F --self_contained T'
TEST_OUTPUTS[self_contained]="linear self_contained"

# Test 5: Directional Z-scores
TEST_ARGS[directional]='--gmt_file "$GMT_FILE" --competitive T --self_contained F --directional T'
TEST_OUTPUTS[directional]="linear"

# Test 6: abs(TWAS.Z) (probit_P_as_Z=F, directional=F)
TEST_ARGS[abs_z]='--gmt_file "$GMT_FILE" --competitive T --self_contained F --probit_P_as_Z F --directional F'
TEST_OUTPUTS[abs_z]="linear"

# Test 7: Two-sided p-values (competitive fast)
TEST_ARGS[two_sided]='--gmt_file "$GMT_FILE" --expression_ref "$EXPR_REF" --competitive T --self_contained F --fast_competitive T --linear_p_thresh 1 --two_sided T'
TEST_OUTPUTS[two_sided]="linear competitive"

# Test 8: Covariate (competitive fast)
TEST_ARGS[covar_fast]='--gmt_file "$GMT_FILE" --expression_ref "$EXPR_REF" --competitive T --self_contained F --fast_competitive T --linear_p_thresh 1 --covar NSNP'
TEST_OUTPUTS[covar_fast]="linear competitive"

# Test 9: Weights (competitive fast)
TEST_ARGS[weights_fast]='--gmt_file "$GMT_FILE" --expression_ref "$EXPR_REF" --competitive T --self_contained F --fast_competitive T --linear_p_thresh 1 --weights MODELCV.R2'
TEST_OUTPUTS[weights_fast]="linear competitive"

# Test 10: Covariate + Weights (competitive fast)
TEST_ARGS[covar_weights_fast]='--gmt_file "$GMT_FILE" --expression_ref "$EXPR_REF" --competitive T --self_contained F --fast_competitive T --linear_p_thresh 1 --covar NSNP --weights MODELCV.R2'
TEST_OUTPUTS[covar_weights_fast]="linear competitive"

# Test 11: Covariate + Weights (competitive orig, for cross-check)
TEST_ARGS[covar_weights_orig]='--gmt_file "$GMT_FILE" --expression_ref "$EXPR_REF" --competitive T --self_contained F --fast_competitive F --linear_p_thresh 1 --covar NSNP --weights MODELCV.R2'
TEST_OUTPUTS[covar_weights_orig]="linear competitive"

# Test 12: Pre-computed correlation matrix
TEST_ARGS[precomp_cormat]='--gmt_file "$GMT_FILE" --input_CorMat "$OUTDIR/competitive_fast.CorMat.RDS" --competitive T --self_contained F --fast_competitive T --linear_p_thresh 1'
TEST_OUTPUTS[precomp_cormat]="linear competitive"

# ---- Ordered test list (test 7 depends on test 2's CorMat) -------------------
TESTS_PHASE1=(linear_only competitive_fast competitive_orig self_contained directional abs_z two_sided covar_fast weights_fast covar_weights_fast covar_weights_orig)
TESTS_PHASE2=(precomp_cormat)

# ---- Run a single test -------------------------------------------------------
run_test() {
  local name="$1"
  local extra_args
  eval "extra_args=(${TEST_ARGS[$name]})"

  echo -n "  Running $name... "
  if $RSCRIPT "$TWAS" "${BASE_ARGS[@]}" "${extra_args[@]}" \
      --output "$OUTDIR/$name" \
      > "$OUTDIR/$name.stdout" 2>&1; then
    echo "OK"
    return 0
  else
    echo "CRASHED (exit code $?)"
    echo "    See: $OUTDIR/$name.stdout"
    return 1
  fi
}

# ---- Compare or update reference for a single test ---------------------------
check_test() {
  local name="$1"
  local suffixes="${TEST_OUTPUTS[$name]}"
  local test_pass=true

  for suffix in $suffixes; do
    local test_file="$OUTDIR/${name}.${suffix}.txt"
    local ref_file="$REFDIR/${name}.${suffix}.txt"

    if [[ ! -f "$test_file" ]]; then
      echo "    FAIL: $test_file not produced"
      test_pass=false
      continue
    fi

    if $UPDATE_MODE; then
      mkdir -p "$REFDIR"
      cp "$test_file" "$ref_file"
      echo "    Updated reference: ${name}.${suffix}.txt"
    else
      if [[ ! -f "$ref_file" ]]; then
        echo "    SKIP: no reference file for ${name}.${suffix}.txt (run with --update first)"
        test_pass=false
        continue
      fi
      echo "  ${name}.${suffix}.txt:"
      if ! $RSCRIPT "$COMPARE" "$test_file" "$ref_file" 2>/dev/null; then
        test_pass=false
      fi
    fi
  done

  $test_pass
}

# ---- Main --------------------------------------------------------------------
echo "=== TWAS-GSEA Regression Tests ==="
echo ""

if $UPDATE_MODE; then
  echo "MODE: Updating reference files"
else
  echo "MODE: Testing against reference files"
fi
echo ""

# Phase 1: independent tests
N_PASS=0
N_FAIL=0
N_TOTAL=0

# Need CorMat from competitive_fast for precomp_cormat test — save it
# Override competitive_fast to also save CorMat
TEST_ARGS[competitive_fast]='--gmt_file "$GMT_FILE" --expression_ref "$EXPR_REF" --competitive T --self_contained F --fast_competitive T --linear_p_thresh 1 --save_CorMat T'

echo "--- Phase 1: Core tests ---"
for name in "${TESTS_PHASE1[@]}"; do
  N_TOTAL=$((N_TOTAL + 1))
  if run_test "$name" && check_test "$name"; then
    N_PASS=$((N_PASS + 1))
  else
    N_FAIL=$((N_FAIL + 1))
  fi
  echo ""
done

echo "--- Phase 2: Dependent tests ---"
for name in "${TESTS_PHASE2[@]}"; do
  N_TOTAL=$((N_TOTAL + 1))
  if [[ ! -f "$OUTDIR/competitive_fast.CorMat.RDS" ]]; then
    echo "  SKIP $name: competitive_fast.CorMat.RDS not available"
    N_FAIL=$((N_FAIL + 1))
  elif run_test "$name" && check_test "$name"; then
    N_PASS=$((N_PASS + 1))
  else
    N_FAIL=$((N_FAIL + 1))
  fi
  echo ""
done

# ---- Cross-check: fast vs original competitive should agree ------------------
echo "--- Cross-checks: fast vs original competitive concordance ---"
for pair in "competitive_fast:competitive_orig" "covar_weights_fast:covar_weights_orig"; do
  fast_name="${pair%%:*}"
  orig_name="${pair##*:}"
  fast_file="$OUTDIR/${fast_name}.competitive.txt"
  orig_file="$OUTDIR/${orig_name}.competitive.txt"
  if [[ -f "$fast_file" && -f "$orig_file" ]]; then
    echo "  ${fast_name} vs ${orig_name} (tolerance 0.1):"
    if $RSCRIPT "$COMPARE" "$fast_file" "$orig_file" 0.1 2>/dev/null; then
      echo "  Concordance: PASS"
    else
      echo "  Concordance: FAIL"
    fi
  else
    echo "  SKIP: ${fast_name} or ${orig_name} competitive output missing"
  fi
done
echo ""

# ---- Summary -----------------------------------------------------------------
echo "=== Summary ==="
echo "  $N_PASS / $N_TOTAL tests passed"
if [[ $N_FAIL -gt 0 ]]; then
  echo "  $N_FAIL FAILED"
  exit 1
else
  if $UPDATE_MODE; then
    echo "  Reference files updated in: $REFDIR/"
  else
    echo "  All tests PASSED"
  fi
fi
