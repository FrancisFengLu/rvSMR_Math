# =====================================================================
# generate_test_data_v3a.R  (Round 3 / Worker A)
#
# Round 3 extends the Round-2 generator with:
#
#   (5) `dgp = "signalt"`     : sign-alternated alpha_k. The instrument-
#                                on-exposure effects alternate sign
#                                (alpha_1=+, alpha_2=-, alpha_3=+, ...)
#                                so the pooled b_x summary statistic is
#                                near zero in expectation even when
#                                individual |alpha_k| is large. This is
#                                the canonical Wang-Kang 2022 / Fig 6
#                                style stress scenario for non-robust
#                                scalar IVW: the IVW numerator
#                                sum(b_x_k * w_k * b_y_k) is small while
#                                sum(b_x_k * w_k * b_x_k) is also small,
#                                producing heavy-tailed ratios and
#                                under-coverage. AR is sign-invariant in
#                                the moment m_k = b_y_k - beta * b_x_k
#                                and should be unaffected.
#
# All existing modes (ar / confounder / ld_xx / overlap) are preserved
# verbatim from v2; only signs (and the `signalt` branch) are new.
#
# CRITIQUE_v2 §S2.v2.1 (seed aliasing across adjacent scenarios) is
# fixed in run_tests_v3a.R, not here -- per-rep seeds are now computed
# via digest::digest2int(paste(scenario_name, master_seed, r)). This
# file's `selfcheck_F_v3a()` uses set.seed(seed) only for the F-axis
# sanity check (where independent seeds are not required).
#
# References:
#  - Wang & Kang 2022 Biometrics §3 (DOI 10.1111/biom.13524) for the
#    weak-IV / sign-alternated-alpha collapse logic in summary-IVW.
#  - Patel, Lane & Burgess 2024 (arXiv 2408.09868) Fig 2 for the
#    scalar-IVW-collapses comparator demonstration.
#  - Burgess, Davies, Thompson 2016 Genet Epidemiol on sample overlap.
# =====================================================================

simulate_burden_mr_v3a <- function(K               = 3,
                                   beta_true       = 0.4,
                                   F_target        = 20,
                                   n_x             = 10000,
                                   n_y             = 10000,
                                   pleiotropy_frac = 0,
                                   pleio_size      = NULL,
                                   dgp             = c("ar", "confounder",
                                                       "ld_xx", "overlap",
                                                       "signalt"),
                                   conf_strength   = 0,
                                   rho_xx          = 0,
                                   rho_yy          = 0,
                                   rho_xy_diag     = 0,
                                   sign_pattern    = NULL,
                                   seed            = NULL) {
  if (!is.null(seed)) set.seed(seed)
  dgp <- match.arg(dgp)
  stopifnot(K >= 1L, beta_true >= -10, beta_true <= 10,
            F_target >= 0, n_x > 0, n_y > 0,
            pleiotropy_frac >= 0, pleiotropy_frac <= 1)

  se_x_vec <- rep(1 / sqrt(n_x), K)
  se_y_vec <- rep(1 / sqrt(n_y), K)

  # Sign pattern: default = all +1. signalt overrides with alternating.
  # User may also supply an explicit `sign_pattern` (length-K +-1 vector).
  if (is.null(sign_pattern)) {
    if (dgp == "signalt") {
      signs <- rep(c(1, -1), length.out = K)
    } else {
      signs <- rep(1, K)
    }
  } else {
    stopifnot(length(sign_pattern) == K,
              all(abs(sign_pattern) == 1))
    signs <- as.numeric(sign_pattern)
  }
  alpha_vec <- signs * sqrt(F_target) * se_x_vec

  if (is.null(pleio_size)) pleio_size <- 5 * mean(se_y_vec)
  K_pleio <- floor(K * pleiotropy_frac)
  pleio_vec <- rep(0, K)
  if (K_pleio > 0L) {
    pleio_idx <- seq.int(K - K_pleio + 1L, K)
    pleio_vec[pleio_idx] <- pleio_size
  }

  mu_x_base <- alpha_vec
  mu_y_base <- beta_true * alpha_vec + pleio_vec

  R_xx <- diag(K)
  R_yy <- diag(K)
  R_xy <- matrix(0, K, K)

  b_x <- numeric(K)
  b_y <- numeric(K)

  if (dgp == "ar" || dgp == "signalt") {
    # AR-structured generator. signalt differs only in the sign pattern
    # of alpha_vec; the residual covariance is the same as `ar`.
    if (rho_xy_diag != 0) {
      R_xy <- diag(rho_xy_diag, K, K)
    }
    for (k in seq_len(K)) {
      rho_k <- R_xy[k, k]
      if (abs(rho_k) < .Machine$double.eps) {
        b_x[k] <- rnorm(1L, mean = mu_x_base[k], sd = se_x_vec[k])
        b_y[k] <- rnorm(1L, mean = mu_y_base[k], sd = se_y_vec[k])
      } else {
        Sigma_k <- matrix(c(se_x_vec[k]^2,
                            rho_k * se_x_vec[k] * se_y_vec[k],
                            rho_k * se_x_vec[k] * se_y_vec[k],
                            se_y_vec[k]^2),
                          2L, 2L, byrow = TRUE)
        L <- chol(Sigma_k)
        z <- rnorm(2L)
        xy <- as.numeric(crossprod(L, z))
        b_x[k] <- mu_x_base[k] + xy[1L]
        b_y[k] <- mu_y_base[k] + xy[2L]
      }
    }
  } else if (dgp == "confounder") {
    stopifnot(abs(conf_strength) <= 1)
    u <- rnorm(1L)
    eps_x <- rnorm(K)
    eps_y <- rnorm(K)
    cs <- conf_strength
    s_idio <- sqrt(max(0, 1 - cs^2))
    b_x <- mu_x_base + cs * se_x_vec * u + s_idio * se_x_vec * eps_x
    b_y <- mu_y_base + cs * se_y_vec * u + s_idio * se_y_vec * eps_y
    # Canonical stress test: user does NOT know about u. Inference call
    # receives R_xx = R_yy = I, R_xy = 0.
  } else if (dgp == "ld_xx") {
    stopifnot(abs(rho_xx) < 1)
    R_xx <- (1 - rho_xx) * diag(K) + rho_xx * matrix(1, K, K)
    Sigma_x <- diag(se_x_vec) %*% R_xx %*% diag(se_x_vec)
    Lx <- chol(Sigma_x)
    zx <- rnorm(K)
    b_x <- mu_x_base + as.numeric(crossprod(Lx, zx))
    b_y <- mu_y_base + rnorm(K, sd = se_y_vec)
  } else if (dgp == "overlap") {
    stopifnot(abs(rho_xy_diag) < 1)
    rho <- rho_xy_diag
    R_xx <- (1 - rho) * diag(K) + rho * matrix(1, K, K)
    R_yy <- (1 - rho) * diag(K) + rho * matrix(1, K, K)
    R_xy <- rho * diag(K) + 0.5 * rho * (matrix(1, K, K) - diag(K))
    Dxy <- c(se_x_vec, se_y_vec)
    Joint_corr <- rbind(cbind(R_xx, R_xy),
                        cbind(t(R_xy), R_yy))
    Sigma_joint <- diag(Dxy) %*% Joint_corr %*% diag(Dxy)
    eig <- eigen(Sigma_joint, symmetric = TRUE, only.values = TRUE)$values
    if (min(eig) < 1e-12) {
      Sigma_joint <- Sigma_joint + diag(max(1e-10, -min(eig) + 1e-12),
                                        nrow(Sigma_joint))
    }
    Lj <- chol(Sigma_joint)
    z <- rnorm(2 * K)
    draw <- as.numeric(crossprod(Lj, z))
    b_x <- mu_x_base + draw[1:K]
    b_y <- mu_y_base + draw[(K + 1):(2 * K)]
  }

  truth <- list(
    beta_true            = beta_true,
    alpha_vec            = alpha_vec,
    signs                = signs,
    F_per_mask           = (alpha_vec / se_x_vec)^2,
    pleiotropy_per_mask  = pleio_vec,
    pleiotropy_frac      = pleiotropy_frac,
    pleio_size           = pleio_size,
    R_xx_true            = R_xx,
    R_yy_true            = R_yy,
    R_xy_true            = R_xy,
    dgp                  = dgp,
    conf_strength        = conf_strength,
    rho_xx               = rho_xx,
    rho_xy_diag          = rho_xy_diag,
    K                    = K,
    n_x                  = n_x,
    n_y                  = n_y,
    lambda_joint         = K * F_target
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
    R_xy  = R_xy,
    truth = truth
  )
}

# Quick F-axis sanity check (R_xx=I baseline) -- not used in driver.
selfcheck_F_v3a <- function(K = 3, F_target = 20, n_x = 10000,
                            n_reps = 100L, seed = 42L,
                            dgp = "ar") {
  set.seed(seed)
  Fs <- numeric(0)
  for (r in seq_len(n_reps)) {
    d <- simulate_burden_mr_v3a(K = K, F_target = F_target,
                                n_x = n_x, n_y = n_x, dgp = dgp)
    Fs <- c(Fs, (d$b_x / d$se_x)^2)
  }
  list(mean_F = mean(Fs),
       median_F = median(Fs),
       q025 = quantile(Fs, 0.025),
       q975 = quantile(Fs, 0.975),
       target = F_target,
       lambda_joint_expected = K * F_target)
}

if (interactive() || identical(Sys.getenv("RVMR_SELFCHECK"), "1")) {
  cat("simulate_burden_mr_v3a() loaded.\n")
  print(selfcheck_F_v3a())
}
