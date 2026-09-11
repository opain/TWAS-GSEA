#!/usr/bin/env bash
# sim/run_r.sh
#
# Wrapper for the twas_gsea conda-env Rscript. Uses `env -i` to strip the host
# R_LIBS_USER / R_LIBS_SITE, otherwise the container's R 4.2 user library
# leaks into R 4.0 and causes undefined-symbol errors (e.g. rlang.so).
#
# Usage:  sim/run_r.sh <script.R> [args...]

set -euo pipefail
ENV_BIN=/home/claude/micromamba/envs/twas_gsea/bin
# HOME must NOT be /home/claude — its ~/.Renviron hardcodes R_LIBS_USER to the
# R 4.2 user library, which crashes R 4.0 with ABI mismatch on shared objects.
FAKE_HOME=/tmp/twas_gsea_home
mkdir -p "$FAKE_HOME"
exec env -i \
    PATH=$ENV_BIN:/usr/bin \
    HOME=$FAKE_HOME \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 \
    R_LIBS_USER=/nonexistent-force-env-only \
    R_LIBS_SITE=/nonexistent-force-env-only \
    "$ENV_BIN/Rscript" "$@"
