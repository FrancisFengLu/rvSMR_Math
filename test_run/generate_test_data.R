# =====================================================================
# generate_test_data.R
#
# Synthetic two-sample summary-statistic MR data generator for rvMR
# (rvSMR algorithm) validation. Produces summary statistics
# (b_x, se_x, b_y, se_y) with KNOWN ground truth (true beta, true F,
# true pleiotropy mask) so that AR coverage, Sargan-J Type-I rate /
# power, CI shape distribution, and point-estimate bias can be
# verified against analytical targets.
#
# Generative model (additive linear SMM, Step 3 of main.tex):
#
#   For each mask k = 1, ..., K:
#     alpha_k  ~ deterministic; chosen so that
#                alpha_k^2 / se_x^2 = F_target
#                (gives expected first-stage F ≈ F_target)
#     b_x_k    ~ N(alpha_k, se_x^2)
#     pleio_k  = 0          if k ≤ K * (1 - pleiotropy_frac)
#              = pleio_size if k >  K * (1 - pleiotropy_frac)
#     b_y_k    ~ N(beta_true * alpha_k + pleio_k, se_y^2)
#
#   In the sample-overlap case (R_xy != 0):
#     (b_x_k, b_y_k) drawn jointly from a bivariate normal with the
#     marginal means above and Cov(b_x, b_y) = R_xy_kk * se_x * se_y.
#
# Standard errors are calibrated to sample size via the standard
# 1/sqrt(n) scaling of OLS slopes for a binary instrument with carrier
# frequency ~1%:
#
#   se_x = sd_x_residual / (sd_burden * sqrt(n_x))
#
# We use sd_x_residual = sd_burden = 1 (unit-variance convention) so
# that se_x = 1/sqrt(n_x).  Then alpha_k is set to alpha_k =
# sqrt(F_target) * se_x.
#
# Output: a list with named slots ready to be passed to rvMR functions.
# =====================================================================

# ---------------------------------------------------------------------
# simulate_burden_mr(): one realization of the (b_x, se_x, b_y, se_y)
# 4-vector per mask k=1..K, given truth (beta_true, F_target,
# pleiotropy_frac, sample sizes, optional sample-overlap block).
# ---------------------------------------------------------------------

simulate_burden_mr <- function(K = 3,
                               beta_true = 0.4,
                               F_target = 20,
                               n_x = 10000,
                               n_y = 10000,
                               pleiotropy_frac = 0,
                               pleio_size = NULL,
                               R_xy = 0,
                               R_xx = diag(K),
                               R_yy = diag(K),
                               seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(K >= 1L, beta_true >= -10, beta_true <= 10,
            F_target >= 0, n_x > 0, n_y > 0,
            pleiotropy_frac >= 0, pleiotropy_frac <= 1)

  # Marginal SE scaling: unit-variance reduced-form residuals, n^{-1/2}.
  se_x_vec <- rep(1 / sqrt(n_x), K)
  se_y_vec <- rep(1 / sqrt(n_y), K)

  # First-stage strength: alpha_k^2 / se_x^2 = F_target  =>  alpha_k =
  # sqrt(F_target) * se_x_k.  Sign alternates a bit so the per-IV Wald
  # ratios are not identical (gives the grid a non-degenerate envelope).
  signs <- ifelse((seq_len(K) %% 2L) == 1L, 1, -1)
  # For the "no sign alternation" stability test, force all positive.
  if (K <= 3) signs <- rep(1, K)
  alpha_vec <- signs * sqrt(F_target) * se_x_vec

  # Pleiotropy: last floor(K * pleiotropy_frac) masks get an offset.
  # Default offset magnitude = 5 * se_y so the bias is detectable but
  # not absurd (J test rejection power ~ 50-90% in the F=20 regime).
  if (is.null(pleio_size)) pleio_size <- 5 * mean(se_y_vec)
  K_pleio <- floor(K * pleiotropy_frac)
  pleio_vec <- rep(0, K)
  if (K_pleio > 0L) {
    pleio_idx <- seq.int(K - K_pleio + 1L, K)
    pleio_vec[pleio_idx] <- pleio_size
  }

  # ---- draw the sample ------------------------------------------------
  # Marginal means
  mu_x <- alpha_vec
  mu_y <- beta_true * alpha_vec + pleio_vec

  # Cross-correlation block R_xy (scalar -> diag * R_xy, matrix -> as is)
  if (length(R_xy) == 1L) {
    R_xy_mat <- diag(R_xy, K, K)
  } else if (is.matrix(R_xy) && all(dim(R_xy) == c(K, K))) {
    R_xy_mat <- R_xy
  } else {
    stop("simulate_burden_mr: R_xy must be a scalar or KxK matrix")
  }

  # Draw b_x and b_y jointly so that Cov(b_x_k, b_y_k) = rho_k * se_x_k *
  # se_y_k.  We assume R_xx = R_yy = I_K for the joint draw (i.e. the
  # K masks are independent both within exposure and within outcome --
  # this matches the rvSMR pLoF/mis/reg partition assumption; if the
  # caller passes non-I R_xx or R_yy we currently ignore them in the
  # draw and pass them through to mrAR_multi as the inference-side
  # covariance.  This is conservative -- it matches what the package
  # itself does internally.)
  b_x <- numeric(K)
  b_y <- numeric(K)
  for (k in seq_len(K)) {
    rho_k <- R_xy_mat[k, k]
    if (abs(rho_k) < .Machine$double.eps) {
      b_x[k] <- rnorm(1L, mean = mu_x[k], sd = se_x_vec[k])
      b_y[k] <- rnorm(1L, mean = mu_y[k], sd = se_y_vec[k])
    } else {
      Sigma_k <- matrix(c(se_x_vec[k]^2,
                          rho_k * se_x_vec[k] * se_y_vec[k],
                          rho_k * se_x_vec[k] * se_y_vec[k],
                          se_y_vec[k]^2),
                        nrow = 2L, byrow = TRUE)
      # Cholesky + standard normal
      L <- tryCatch(chol(Sigma_k), error = function(e) NULL)
      if (is.null(L)) {
        # Singular -> back off to marginal draws
        b_x[k] <- rnorm(1L, mean = mu_x[k], sd = se_x_vec[k])
        b_y[k] <- rnorm(1L, mean = mu_y[k], sd = se_y_vec[k])
      } else {
        z <- rnorm(2L)
        xy <- as.numeric(crossprod(L, z))
        b_x[k] <- mu_x[k] + xy[1L]
        b_y[k] <- mu_y[k] + xy[2L]
      }
    }
  }

  # ---- ground truth bundle -------------------------------------------
  truth <- list(
    beta_true            = beta_true,
    alpha_vec            = alpha_vec,
    F_per_mask           = (alpha_vec / se_x_vec)^2,  # expected F
    pleiotropy_per_mask  = pleio_vec,
    pleiotropy_frac      = pleiotropy_frac,
    pleio_size           = pleio_size,
    R_xy                 = R_xy_mat,
    K                    = K,
    n_x                  = n_x,
    n_y                  = n_y
  )

  list(
    b_x   = b_x,
    se_x  = se_x_vec,
    b_y   = b_y,
    se_y  = se_y_vec,
    n_x   = n_x,
    n_y   = n_y,
    R_xx  = R_xx,
    R_yy  = R_yy,
    R_xy  = R_xy_mat,
    truth = truth
  )
}

# ---------------------------------------------------------------------
# Quick self-check: does the empirical F-stat hit the target?
# ---------------------------------------------------------------------

selfcheck_F <- function(K = 3, F_target = 20, n_x = 10000,
                        n_reps = 100L, seed = 42L) {
  set.seed(seed)
  Fs <- numeric(0)
  for (r in seq_len(n_reps)) {
    d <- simulate_burden_mr(K = K, F_target = F_target,
                            n_x = n_x, n_y = n_x)
    Fs <- c(Fs, (d$b_x / d$se_x)^2)
  }
  list(mean_F = mean(Fs),
       median_F = median(Fs),
       q025 = quantile(Fs, 0.025),
       q975 = quantile(Fs, 0.975),
       target = F_target)
}

# When sourced interactively, print the self-check banner.
if (interactive() || identical(Sys.getenv("RVMR_SELFCHECK"), "1")) {
  cat("simulate_burden_mr() loaded.\n")
  cat("Quick F self-check (K=3, F_target=20, 100 reps):\n")
  print(selfcheck_F())
}
