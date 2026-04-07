# reml_blockdiag.R
#
# Restricted-maximum-likelihood (REML) estimation of (sigma_u^2, sigma_e^2) in
# the 2-variance-component LMM
#
#     y = X beta + g + e,    g ~ N(0, sigma_u^2 K),    e ~ N(0, sigma_e^2 I)
#
# specialised to the case where K is BLOCK-DIAGONAL with blocks K_1, ..., K_B.
#
# Reparameterise V = sigma^2 (h K + I)  with  h = sigma_u^2 / sigma_e^2,
# profile out sigma^2, and optimise the 1-D profiled REML in h via Brent's
# method (stats::optimise). All per-iteration work is element-wise on cached
# per-block eigendecompositions, so the cost is dominated by the upfront
# eigendecomposition Σ_b n_b^3 (cheap when blocks are small).
#
# This is mathematically the same MLE as lme4qtl::relmatLmer on the same
# model — not an approximation.

# Build the cached per-block rotated quantities. K_blocks is a list of dense
# (small) symmetric block matrices. y is the response, X is the null design.
# Returns a list with everything fit_reml_profile() needs.
prepare_reml_blocks <- function(y, X, K_blocks, block_index) {
	stopifnot(length(y) == length(block_index))
	stopifnot(nrow(X) == length(block_index))
	N <- length(y); p <- ncol(X); B <- length(K_blocks)

	# Precompute per-block: eigenvalues d_b, rotated y_b_tilde, rotated X_b_tilde
	cached <- vector('list', B)
	for(b in seq_len(B)){
		idx <- which(block_index == b)
		stopifnot(length(idx) == nrow(K_blocks[[b]]))
		Kb <- as.matrix(K_blocks[[b]])
		yb <- y[idx]
		Xb <- X[idx, , drop = FALSE]
		if(nrow(Kb) == 1){
			# Singleton block: K_b = [1], U = [1], d = 1.
			cached[[b]] <- list(d = 1, yt = yb, Xt = Xb)
		} else {
			eg <- eigen(Kb, symmetric = TRUE)
			U  <- eg$vectors
			d  <- eg$values
			# Numerical: clamp tiny negative eigenvalues from finite-precision noise
			d[d < 0] <- 0
			cached[[b]] <- list(d = d, yt = crossprod(U, yb), Xt = crossprod(U, Xb))
		}
	}
	list(N = N, p = p, B = B, blocks = cached)
}

# Profiled REML log-likelihood at a given variance ratio h. Returns 2*ll up to
# an additive constant (sufficient for argmax).
reml_profile_ll <- function(h, prep) {
	N <- prep$N; p <- prep$p
	A     <- matrix(0, p, p)   # X' W^{-1} X
	bvec  <- numeric(p)         # X' W^{-1} y
	yWy   <- 0                  # y' W^{-1} y
	logdetW <- 0                # log |hK + I|

	for(blk in prep$blocks){
		lam <- h * blk$d + 1                # eigenvalues of W_b
		inv <- 1 / lam
		logdetW <- logdetW + sum(log(lam))
		yt  <- as.numeric(blk$yt)
		Xt  <- blk$Xt
		yWy   <- yWy   + sum(yt^2 * inv)
		bvec  <- bvec  + crossprod(Xt, yt * inv)
		A     <- A     + crossprod(Xt * inv, Xt)
	}

	# Solve A x = bvec robustly (A is p x p, p small)
	A_chol <- tryCatch(chol(A), error = function(e) NULL)
	if(is.null(A_chol)){
		# Degenerate design — return very negative ll so optimiser avoids it
		return(-Inf)
	}
	x       <- backsolve(A_chol, backsolve(A_chol, bvec, transpose = TRUE))
	logdetA <- 2 * sum(log(diag(A_chol)))
	Q       <- yWy - sum(bvec * x)
	if(Q <= 0) return(-Inf)

	# 2 * profiled REML ll (constants dropped)
	-(N - p) * log(Q) - logdetW - logdetA
}

# Public entry point.
#
# Returns:
#   sigma_u2, sigma_e2  - REML variance components
#   h_hat               - argmax of the profiled likelihood
#   ll                  - profiled REML log-lik (up to constant) at h_hat
#   converged           - logical
fit_reml_blockdiag <- function(y, X, K_blocks, block_index,
                                h_max = 100, tol = 1e-6) {
	prep <- prepare_reml_blocks(y, X, K_blocks, block_index)

	# 1-D Brent search. Minimise the negative log-likelihood.
	opt <- optimise(function(h) -reml_profile_ll(h, prep),
	                interval = c(0, h_max), tol = tol)
	h_hat <- opt$minimum

	# Recover sigma_e^2 from the closed-form profile: sigma_e^2 = Q(h) / (N - p)
	N <- prep$N; p <- prep$p
	A     <- matrix(0, p, p)
	bvec  <- numeric(p)
	yWy   <- 0
	for(blk in prep$blocks){
		lam <- h_hat * blk$d + 1
		inv <- 1 / lam
		yt  <- as.numeric(blk$yt)
		Xt  <- blk$Xt
		yWy  <- yWy  + sum(yt^2 * inv)
		bvec <- bvec + crossprod(Xt, yt * inv)
		A    <- A    + crossprod(Xt * inv, Xt)
	}
	x <- solve(A, bvec)
	Q <- yWy - sum(bvec * x)
	sigma_e2 <- Q / (N - p)
	sigma_u2 <- h_hat * sigma_e2

	list(sigma_u2 = sigma_u2,
	     sigma_e2 = sigma_e2,
	     h_hat    = h_hat,
	     ll       = -opt$objective,
	     converged = TRUE)
}
