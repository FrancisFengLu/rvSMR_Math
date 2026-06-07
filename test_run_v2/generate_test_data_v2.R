# =====================================================================
# generate_test_data_v2.R  (Round 2)
#
# Round 1 generator drew (b_x_k, b_y_k) ~ N(alpha_k, beta*alpha_k + plei,
# Sigma) under R_xx=R_yy=I, which exactly matches the AR statistic's
# moment structure (CRITIQUE S1.3, S2.1, S2.2). Round 2 extends with:
#
#   (1) `dgp = "ar"`           : Round-1 default (kept for anchor cells).
#   (2) `dgp = "confounder"`   : non-AR-structured DGP via an unobserved
#                                shared latent u (CRITIQUE S2.1, Tier 4).
#                                b_x_k = alpha_k + conf_strength * u + eps_x
#                                b_y_k = beta*alpha_k + conf_strength*u + eps_y
#                                The joint (b_x, b_y) is no longer the
#                                exact AR moment structure: residual u
#                                induces a *common* (correlated across
#                                masks) shift in both b_x and b_y.
#   (3) `dgp = "ld_xx"`        : R_xx != I, cross-mask correlation
#                                rho_xx between exposure burden estimates
#                                (CRITIQUE S2.2). The inference call is
#                                supplied with the correct R_xx.
#   (4) `dgp = "overlap"`      : honest sample-overlap analog -- shared
#                                per-mask cross correlation rho_xy on the
#                                bivariate draw, AND non-zero off-diagonal
#                                R_xx and R_yy from a shared individual-
#                                level latent (CRITIQUE S3, "honest"
#                                handling).
#
# All four modes share the (alpha, SE) calibration so the F-stat axis
# is comparable: SE_x = 1/sqrt(n_x), alpha_k = sign_k * sqrt(F_target) *
# SE_x, so per-mask E[F_k] = F_target + 1 and joint concentration
# parameter (R_xx=I) lambda = K * F_target (CRITIQUE S1.4).
#
# References:
#  - Burgess, Davies, Thompson 2016 Genet Epidemiol on sample overlap
#    (replaces the fabricated steps_5_to_9_logic.md citation, CRITIQUE
#    S1.5).
# =====================================================================

simulate_burden_mr_v2 <- function(K               = 3,
                                  beta_true       = 0.4,
                                  F_target        = 20,
                                  n_x             = 10000,
                                  n_y             = 10000,
                                  pleiotropy_frac = 0,
                                  pleio_size      = NULL,
                                  dgp             = c("ar", "confounder",
                                                      "ld_xx", "overlap"),
                                  conf_strength   = 0,
                                  rho_xx          = 0,
                                  rho_yy          = 0,
                                  rho_xy_diag     = 0,
                                  seed            = NULL) {
  if (!is.null(seed)) set.seed(seed)
  dgp <- match.arg(dgp)
  stopifnot(K >= 1L, beta_true >= -10, beta_true <= 10,
            F_target >= 0, n_x > 0, n_y > 0,
            pleiotropy_frac >= 0, pleiotropy_frac <= 1)

  se_x_vec <- rep(1 / sqrt(n_x), K)
  se_y_vec <- rep(1 / sqrt(n_y), K)

  signs <- rep(1, K)
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

  if (dgp == "ar") {
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
    # Shared latent u affects b_x AND b_y across masks; variance budget
    # preserved so empirical F is comparable. Within-X cross-mask cor =
    # conf_strength^2; same for Y; within-mask X-Y cor = conf_strength^2.
    stopifnot(abs(conf_strength) <= 1)
    u <- rnorm(1L)
    eps_x <- rnorm(K)
    eps_y <- rnorm(K)
    cs <- conf_strength
    s_idio <- sqrt(max(0, 1 - cs^2))
    b_x <- mu_x_base + cs * se_x_vec * u + s_idio * se_x_vec * eps_x
    b_y <- mu_y_base + cs * se_y_vec * u + s_idio * se_y_vec * eps_y
    # We do NOT propagate the induced cross-correlation into R_xx/R_yy/R_xy
    # because the canonical confounder scenario is: user does not know
    # about u. AR is given R_xx=R_yy=I, R_xy=0 (the default). The point
    # is to stress-test AR under misspecified covariance.
  } else if (dgp == "ld_xx") {
    # LD between masks: compound-symmetric R_xx with cross-mask cor rho_xx.
    stopifnot(abs(rho_xx) < 1)
    R_xx <- (1 - rho_xx) * diag(K) + rho_xx * matrix(1, K, K)
    Sigma_x <- diag(se_x_vec) %*% R_xx %*% diag(se_x_vec)
    Lx <- chol(Sigma_x)
    zx <- rnorm(K)
    b_x <- mu_x_base + as.numeric(crossprod(Lx, zx))
    b_y <- mu_y_base + rnorm(K, sd = se_y_vec)
  } else if (dgp == "overlap") {
    # Honest sample overlap: shared per-individual latent induces
    # non-zero off-diagonals in R_xx, R_yy, AND a full R_xy block.
    # Inference call receives the correct (non-diag) R_xx / R_yy / R_xy.
    stopifnot(abs(rho_xy_diag) < 1)
    rho <- rho_xy_diag
    # R_xx, R_yy compound-symmetric (shared-cohort: cross-mask cor=rho).
    # R_xy (cross-block exposure-outcome) = rho on the diagonal
    # (within-mask cross-sample cor), with off-diagonals = 0.5*rho
    # (cross-mask cross-sample cor, smaller because shared only via
    # the cohort latent, not via mask-specific residuals).
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

selfcheck_F_v2 <- function(K = 3, F_target = 20, n_x = 10000,
                           n_reps = 100L, seed = 42L,
                           dgp = "ar") {
  set.seed(seed)
  Fs <- numeric(0)
  for (r in seq_len(n_reps)) {
    d <- simulate_burden_mr_v2(K = K, F_target = F_target,
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
  cat("simulate_burden_mr_v2() loaded.\n")
  print(selfcheck_F_v2())
}
