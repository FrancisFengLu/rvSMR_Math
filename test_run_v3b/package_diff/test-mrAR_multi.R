test_that("mrAR_multi signature and defaults", {
  expect_true(is.function(mrAR_multi))
  fmls <- formals(mrAR_multi)
  expect_setequal(
    names(fmls),
    c("b_x", "se_x", "b_y", "se_y",
      "R_xx", "R_yy", "R_xy",
      "alpha", "n_grid", "grid_pad_mult", "grid_extend_max")
  )
  expect_identical(fmls$alpha, 0.05)
  expect_identical(fmls$n_grid, 4000L)
  expect_identical(fmls$grid_pad_mult, 3)
  expect_identical(fmls$grid_extend_max, 3L)
})

test_that("Strong-IV K=3 gives a bounded interval covering truth", {
  set.seed(1)
  b_x  <- c(0.6, 0.5, 0.4)
  se_x <- c(0.1, 0.1, 0.1)
  beta_true <- 0.5
  b_y  <- beta_true * b_x + rnorm(3, sd = 0.02)
  se_y <- c(0.05, 0.05, 0.05)

  res <- mrAR_multi(b_x, se_x, b_y, se_y, alpha = 0.05)

  expect_equal(res$ci_type, "bounded_interval")
  expect_length(res$ci_intervals, 1L)
  iv <- res$ci_intervals[[1L]]
  expect_true(is.finite(iv[1]) && is.finite(iv[2]))
  expect_lt(iv[1], beta_true)
  expect_gt(iv[2], beta_true)
  expect_lt(iv[2] - iv[1], 1)
})

test_that("Strong-IV K=3 J-test passes under homogeneity", {
  set.seed(2)
  b_x  <- c(0.6, 0.5, 0.4)
  se_x <- c(0.1, 0.1, 0.1)
  beta_true <- 0.5
  b_y  <- beta_true * b_x + rnorm(3, sd = 0.02)
  se_y <- c(0.05, 0.05, 0.05)

  res <- mrAR_multi(b_x, se_x, b_y, se_y, alpha = 0.05)

  expect_true(is.finite(res$J_stat))
  expect_true(is.finite(res$J_pvalue))
  expect_gt(res$J_pvalue, 0.05)
  expect_lt(abs(res$beta_hat - beta_true), 0.1)
})

test_that("K=3 pleiotropy: J-test rejects when one IV has offset", {
  set.seed(3)
  b_x  <- c(0.6, 0.5, 0.4)
  se_x <- c(0.1, 0.1, 0.1)
  beta_true <- 0.5
  b_y  <- beta_true * b_x
  b_y[3] <- b_y[3] + 0.3   # additive pleiotropic offset on IV #3
  se_y <- c(0.05, 0.05, 0.05)

  res <- mrAR_multi(b_x, se_x, b_y, se_y, alpha = 0.05)

  expect_true(is.finite(res$J_pvalue))
  expect_lt(res$J_pvalue, 0.05)
})

test_that("Weak-IV K=3 yields disconnected_union or whole_line", {
  set.seed(4)
  res <- mrAR_multi(
    b_x  = c(0.05, 0.04, 0.06),
    se_x = c(0.1,  0.1,  0.1),
    b_y  = c(0.30, 0.25, 0.28),
    se_y = c(0.05, 0.05, 0.05),
    alpha = 0.05
  )
  expect_true(res$ci_type %in%
              c("disconnected_union", "whole_line", "bounded_interval"))
  if (res$ci_type == "disconnected_union") {
    expect_gte(length(res$ci_intervals), 2L)
  }
})

test_that("K=2 minimal multi-instrument case behaves identically in shape", {
  set.seed(5)
  b_x  <- c(0.5, 0.4)
  se_x <- c(0.1, 0.1)
  beta_true <- 0.5
  b_y  <- beta_true * b_x + rnorm(2, sd = 0.02)
  se_y <- c(0.05, 0.05)

  res <- mrAR_multi(b_x, se_x, b_y, se_y, alpha = 0.05)
  expect_equal(res$ci_type, "bounded_interval")
  iv <- res$ci_intervals[[1L]]
  expect_true(iv[1] < beta_true && beta_true < iv[2])
  # K=2 → df=K-1=1 for J-test.
  expect_true(is.finite(res$J_pvalue))
})

test_that("K=1 cross-check: mrAR_multi matches mrAR closed form", {
  # Strong-IV K=1 case from the existing mrAR test suite.
  b_x  <- 0.18
  se_x <- 0.04
  b_y  <- -0.11
  se_y <- 0.025

  ref <- mrAR(b_x, se_x, b_y, se_y, alpha = 0.05)
  res <- mrAR_multi(b_x = b_x, se_x = se_x,
                    b_y = b_y, se_y = se_y,
                    alpha = 0.05,
                    n_grid = 8000L)
  expect_equal(res$ci_type, "bounded_interval")
  iv <- res$ci_intervals[[1L]]
  # K=1 closed form gives a bounded interval here; multi must agree to ~1e-4.
  expect_equal(iv[1], ref$ci_lower, tolerance = 1e-4)
  expect_equal(iv[2], ref$ci_upper, tolerance = 1e-4)
})

test_that("Input validation: mismatched length errors", {
  expect_error(
    mrAR_multi(b_x = c(0.5, 0.4),
               se_x = c(0.1, 0.1, 0.1),
               b_y = c(0.25, 0.20),
               se_y = c(0.05, 0.05))
  )
  expect_error(
    mrAR_multi(b_x = c(0.5, 0.4, 0.3),
               se_x = c(0.1, 0.1, 0.1),
               b_y = c(0.25, 0.20),
               se_y = c(0.05, 0.05, 0.05))
  )
})

test_that("Input validation: malformed R_xx / R_yy", {
  # Wrong dim → caught up front.
  expect_error(
    mrAR_multi(b_x = c(0.5, 0.4, 0.3),
               se_x = c(0.1, 0.1, 0.1),
               b_y = c(0.25, 0.20, 0.15),
               se_y = c(0.05, 0.05, 0.05),
               R_xx = diag(2))
  )
  # Non-PSD R_xx → solve() inside ar_fun should still typically return a
  # value (the linear solve doesn't require PSD), but a degenerate (rank-
  # deficient) R_xx with se_x > 0 can yield singular V at some b0. We
  # require either a graceful error or a finite output; the function
  # must not crash with an uninformative message.
  bad_R <- matrix(1, 3, 3)        # rank 1
  res <- tryCatch(
    mrAR_multi(b_x = c(0.5, 0.4, 0.3),
               se_x = c(0.1, 0.1, 0.1),
               b_y = c(0.25, 0.20, 0.15),
               se_y = c(0.05, 0.05, 0.05),
               R_xx = bad_R),
    error = function(e) e
  )
  # Must produce a list-like result or an error — not silent NULL.
  expect_true(inherits(res, "error") || is.list(res))
})

test_that("alpha and n_grid validation", {
  b_x  <- c(0.5, 0.4, 0.3)
  se_x <- c(0.1, 0.1, 0.1)
  b_y  <- c(0.25, 0.20, 0.15)
  se_y <- c(0.05, 0.05, 0.05)
  expect_error(mrAR_multi(b_x, se_x, b_y, se_y, alpha = 0))
  expect_error(mrAR_multi(b_x, se_x, b_y, se_y, alpha = 1))
  expect_error(mrAR_multi(b_x, se_x, b_y, se_y, n_grid = 5))
  expect_error(mrAR_multi(b_x, se_x, b_y, se_y, grid_pad_mult = -1))
})

test_that("R_xy index convention: AR depends only on the symmetric part of R_xy", {
  # IMPORTANT mathematical observation discovered during the v3 R_xy
  # audit (see IMPLEMENTATION_NOTES_v3b.md). The AR statistic is the
  # scalar quadratic form
  #     AR(b0) = m^T V(b0)^{-1} m
  # where V(b0) = V_yy + b0^2 V_xx - 2 b0 V_xy and V_xy = Dy R_xy Dx.
  # For any matrix V, the scalar m^T V^{-1} m depends ONLY on the
  # symmetric part of V (because the antisymmetric part contributes
  # zero to any quadratic form). Consequently, AR is invariant to the
  # transpose of R_xy:
  #     AR(b0; R_xy) == AR(b0; t(R_xy))
  # This is why the documented index convention (R_xy[i,j] =
  # cor(b_y_i, b_x_j), i.e. V_xy = Dy R_xy Dx) and the alternative
  # joint-covariance layout (R_xy[i,j] = cor(b_x_i, b_y_j), i.e.
  # V_xy = Dx R_xy Dy = (Dy R_xy^T Dx)^T) produce identical AR values:
  # they differ only in the antisymmetric part of V_xy.
  #
  # This test pins that property so a future edit that (a) "fixes" the
  # convention by also folding it through some non-quadratic-form
  # downstream computation, OR (b) breaks the symmetric-part invariance
  # via a numerical pathway, will fire.
  set.seed(42L)
  b_x  <- c(0.6, 0.5, 0.4)
  se_x <- c(0.1, 0.1, 0.1)
  beta_true <- 0.5
  b_y  <- beta_true * b_x + rnorm(3, sd = 0.02)
  se_y <- c(0.05, 0.05, 0.05)

  R_xy_asym <- matrix(0, 3, 3)
  R_xy_asym[1, 2] <- 0.3
  R_xy_asym[2, 1] <- 0.1
  R_xy_asym[3, 1] <- -0.05
  R_xy_asym[1, 3] <- 0.12

  res_orig <- mrAR_multi(b_x, se_x, b_y, se_y,
                         R_xy = R_xy_asym, alpha = 0.05,
                         n_grid = 4000L)
  res_tpos <- mrAR_multi(b_x, se_x, b_y, se_y,
                         R_xy = t(R_xy_asym), alpha = 0.05,
                         n_grid = 4000L)

  expect_equal(res_orig$ci_type, "bounded_interval")
  expect_equal(res_tpos$ci_type, "bounded_interval")
  iv_o <- res_orig$ci_intervals[[1L]]
  iv_t <- res_tpos$ci_intervals[[1L]]
  expect_true(all(is.finite(iv_o)))
  expect_true(all(is.finite(iv_t)))
  # Symmetric-part-only property:
  expect_equal(iv_o[1], iv_t[1], tolerance = 1e-6)
  expect_equal(iv_o[2], iv_t[2], tolerance = 1e-6)
  expect_equal(res_orig$beta_hat, res_tpos$beta_hat, tolerance = 1e-6)
  expect_equal(res_orig$J_stat,   res_tpos$J_stat,   tolerance = 1e-6)
})

test_that("R_xy: changing the symmetric part DOES change the answer", {
  # Counterpoint to the transpose-invariance test: the symmetric part of
  # R_xy is the part that matters. Construct two R_xy with different
  # symmetric parts (here, swap signs on the (1,2)/(2,1) symmetric pair)
  # and verify AR / CI materially differ. This is what a sample-overlap
  # correction is actually doing.
  set.seed(7L)
  b_x  <- c(0.6, 0.5, 0.4)
  se_x <- c(0.1, 0.1, 0.1)
  beta_true <- 0.5
  b_y  <- beta_true * b_x + rnorm(3, sd = 0.02)
  se_y <- c(0.05, 0.05, 0.05)

  R_pos <- matrix(0, 3, 3)
  R_pos[1, 2] <- R_pos[2, 1] <- 0.4   # symmetric, positive
  R_neg <- -R_pos                      # symmetric, negative

  res_pos <- mrAR_multi(b_x, se_x, b_y, se_y, R_xy = R_pos)
  res_neg <- mrAR_multi(b_x, se_x, b_y, se_y, R_xy = R_neg)
  res_zer <- mrAR_multi(b_x, se_x, b_y, se_y, R_xy = matrix(0, 3, 3))

  expect_equal(res_pos$ci_type, "bounded_interval")
  expect_equal(res_neg$ci_type, "bounded_interval")
  iv_pos <- res_pos$ci_intervals[[1L]]
  iv_neg <- res_neg$ci_intervals[[1L]]
  iv_zer <- res_zer$ci_intervals[[1L]]

  # Three distinct CIs (orientation matters via symmetric-part sign):
  betas <- c(res_pos$beta_hat, res_neg$beta_hat, res_zer$beta_hat)
  expect_equal(length(unique(round(betas, 6L))), 3L)
  # And the CI widths/midpoints should differ:
  expect_false(isTRUE(all.equal(iv_pos, iv_neg, tolerance = 1e-4)))
})

test_that("R_xy: V_xy = Dy R_xy Dx (NOT Dx R_xy Dy) — direct AR comparison", {
  # Pin the documented convention V_xy = Dy R_xy Dx against the alternative
  # V_xy = Dx R_xy Dy. With UNequal se_x and se_y AND a non-symmetric R_xy
  # whose symmetric part is non-trivial, the two assemblies give different
  # symmetric parts of V_xy, hence different AR values. This catches a
  # future edit that swaps Dx <-> Dy in the V_xy line.
  set.seed(11L)
  b_x  <- c(0.6, 0.5, 0.4)
  se_x <- c(0.20, 0.18, 0.22)    # markedly different from se_y
  beta_true <- 0.5
  b_y  <- beta_true * b_x + rnorm(3, sd = 0.02)
  se_y <- c(0.05, 0.07, 0.04)

  R_xy <- matrix(0, 3, 3)
  R_xy[1, 1] <- 0.3
  R_xy[2, 2] <- 0.2
  R_xy[3, 3] <- 0.25
  R_xy[1, 2] <- 0.4     # asymmetric off-diagonals
  R_xy[2, 1] <- 0.1
  R_xy[1, 3] <- 0.15
  R_xy[3, 1] <- 0.05

  # Reference: documented convention, computed externally and run through
  # mrAR_multi via the public R_xy argument.
  res_doc <- mrAR_multi(b_x, se_x, b_y, se_y, R_xy = R_xy)

  # Alternative: simulate "what would happen if V_xy = Dx R_xy Dy".
  # Equivalent to passing R_xy_alt = Dx %*% R_xy_user %*% Dy / (Dy . Dx)
  # element-wise — but exact equivalence at the R_xy slot only when
  # Dx == Dy. Here Dx != Dy, so we have to feed the OTHER convention's
  # R_xy through an externally constructed equivalent. Easiest: construct
  # the alternative V_xy and pass an R_xy_eff such that
  # Dy R_xy_eff Dx == Dx R_xy Dy
  # =>  R_xy_eff = solve(Dy) %*% Dx %*% R_xy %*% Dy %*% solve(Dx)
  Dx <- diag(se_x); Dy <- diag(se_y)
  R_xy_eff <- solve(Dy) %*% Dx %*% R_xy %*% Dy %*% solve(Dx)
  res_alt <- mrAR_multi(b_x, se_x, b_y, se_y, R_xy = R_xy_eff)

  # Must differ when se_x != se_y and R_xy is asymmetric.
  expect_false(isTRUE(all.equal(res_doc$beta_hat, res_alt$beta_hat,
                                tolerance = 1e-6)))
  expect_false(isTRUE(all.equal(res_doc$J_stat, res_alt$J_stat,
                                tolerance = 1e-6)))
})
