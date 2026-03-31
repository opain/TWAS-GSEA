#!/usr/bin/env Rscript
# compare_outputs.R
# Compare two TWAS-GSEA output files and assert numerical agreement.
#
# Usage: Rscript compare_outputs.R <test_file> <reference_file> [tolerance]
#
# Compares all numeric columns (Estimate/Est, SE, T, P, P.CORR, N_Mem_Avail, N_Mem)
# between a test run and a reference file.  Exits with status 1 on failure.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript compare_outputs.R <test_file> <reference_file> [tolerance]")
}

test_file <- args[1]
ref_file  <- args[2]
tol       <- if (length(args) >= 3) as.numeric(args[3]) else 1e-6

test <- read.table(test_file, header = TRUE, stringsAsFactors = FALSE)
ref  <- read.table(ref_file,  header = TRUE, stringsAsFactors = FALSE)

# Align on GeneSet
both <- merge(test, ref, by = "GeneSet", suffixes = c(".test", ".ref"))

if (nrow(both) == 0) {
  cat("FAIL: no overlapping gene sets between test and reference\n")
  quit(status = 1)
}
if (nrow(both) != nrow(ref)) {
  cat(sprintf("FAIL: reference has %d gene sets but only %d matched in test output\n",
              nrow(ref), nrow(both)))
  quit(status = 1)
}

# Find numeric columns to compare (present in both with .test/.ref suffixes)
test_cols <- grep("\\.test$", names(both), value = TRUE)
num_cols  <- sub("\\.test$", "", test_cols)

all_pass <- TRUE
for (col in num_cols) {
  tc <- paste0(col, ".test")
  rc <- paste0(col, ".ref")
  if (!(tc %in% names(both) && rc %in% names(both))) next

  v_test <- both[[tc]]
  v_ref  <- both[[rc]]

  # Skip non-numeric
  if (!is.numeric(v_test) || !is.numeric(v_ref)) next

  max_diff <- max(abs(v_test - v_ref), na.rm = TRUE)
  ok <- max_diff <= tol

  cat(sprintf("  %-15s max|diff| = %.2e  %s\n", col, max_diff, if (ok) "PASS" else "FAIL"))
  if (!ok) all_pass <- FALSE
}

cat(sprintf("  Gene sets compared: %d\n", nrow(both)))

if (!all_pass) {
  cat("RESULT: FAIL\n")
  quit(status = 1)
} else {
  cat("RESULT: PASS\n")
}
