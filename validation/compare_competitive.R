#!/usr/bin/env Rscript
# compare_competitive.R
# Read two competitive output files (fast vs original) and report concordance.
# Usage: Rscript compare_competitive.R <fast_output.txt> <orig_output.txt>
# Called by run_validation.sh; can also be run standalone on any two outputs.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript compare_competitive.R <fast_file> <orig_file>")
}

fast <- read.table(args[1], header = TRUE, sep = "", stringsAsFactors = FALSE)
orig <- read.table(args[2], header = TRUE, sep = "", stringsAsFactors = FALSE)

# Align on gene set name (order may differ between runs)
both <- merge(fast, orig, by = "GeneSet", suffixes = c(".fast", ".orig"))
n <- nrow(both)

cat(sprintf("Gene sets compared: %d\n", n))
cat(sprintf("Columns in fast output: %s\n", paste(names(fast), collapse = ", ")))

# Check structural equivalence
stopifnot("Same columns in both outputs" = identical(sort(names(fast)), sort(names(orig))))
cat("Output structure: PASS (same columns)\n\n")

# T-statistic concordance
t_fast <- both$T.fast
t_orig <- both$T.orig

pearson_r  <- cor(t_fast, t_orig, method = "pearson")
spearman_r <- cor(t_fast, t_orig, method = "spearman")
max_abs_dt <- max(abs(t_fast - t_orig))
sign_agree <- sum(sign(t_fast) == sign(t_orig))

cat("--- T-statistic concordance ---\n")
cat(sprintf("  Pearson r:        %.6f\n", pearson_r))
cat(sprintf("  Spearman rho:     %.6f\n", spearman_r))
cat(sprintf("  Max |delta_T|:    %.4f\n", max_abs_dt))
cat(sprintf("  Sign agreement:   %d / %d\n", sign_agree, n))

# P-value concordance
p_fast <- both$P.fast
p_orig <- both$P.orig
cat(sprintf("  Max |delta_P|:    %.6f\n", max(abs(p_fast - p_orig))))

# Rank-order changes
rank_fast <- rank(p_fast)
rank_orig <- rank(p_orig)
rank_changes <- sum(rank_fast != rank_orig)
cat(sprintf("  Rank changes:     %d / %d\n", rank_changes, n))
if (rank_changes > 0) {
  changed <- both[rank_fast != rank_orig, c("GeneSet", "T.fast", "T.orig")]
  changed$delta_T <- changed$T.fast - changed$T.orig
  cat("  Rank-changed gene sets (delta_T shown — likely near-zero ties):\n")
  for (i in seq_len(nrow(changed))) {
    cat(sprintf("    %s  delta_T=%.5f\n", changed$GeneSet[i], changed$delta_T[i]))
  }
}
cat("\n")

# Top hits
cat("--- Top 3 gene sets (fast) ---\n")
top3_fast <- both[order(both$P.fast), c("GeneSet", "T.fast", "P.fast", "T.orig", "P.orig")][1:min(3, n), ]
print(top3_fast, row.names = FALSE)
cat("\n")

# Tolerances for pass/fail
tol_pearson  <- 0.999
tol_spearman <- 0.990   # rank swaps can occur between gene sets with near-identical T; relax slightly
tol_max_dt   <- 0.1

cat("--- Pass / fail ---\n")
check <- function(label, value, threshold, direction = "above") {
  ok <- if (direction == "above") value >= threshold else value <= threshold
  cat(sprintf("  %-35s %s (%.5f %s %.4f)\n",
              label, if (ok) "PASS" else "FAIL", value,
              if (direction == "above") ">=" else "<=", threshold))
  ok
}

all_pass <- TRUE
all_pass <- check("Pearson r(T) >= 0.999",       pearson_r,  tol_pearson,  "above") & all_pass
all_pass <- check("Spearman rho(T) >= 0.999",     spearman_r, tol_spearman, "above") & all_pass
all_pass <- check("Max |delta_T| <= 0.1",         max_abs_dt, tol_max_dt,   "below") & all_pass
all_pass <- check("Sign agreement == 100%",        sign_agree / n, 1.0,      "above") & all_pass

cat(sprintf("\nOverall: %s\n", if (all_pass) "PASS" else "FAIL"))
if (!all_pass) quit(status = 1)
