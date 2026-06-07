#' Multi-instrument Anderson-Rubin (AR) confidence set for two-sample
#' summary MR
#'
#' Constructs a weak-instrument-robust Anderson-Rubin confidence set
#' for the causal effect `beta = b_xy` using K >= 2 burden
#' instruments. Unlike the closed-form K=1 sibling [mrAR()] (which
#' reduces algebraically to inverted Fieller 1954), the multi-IV AR
#' statistic carries genuine over-identifying information and supports
#' a Sargan-style J-test for instrument-class homogeneity (no
#' horizontal pleiotropy).
#'
#' **Statistic.** With summary stats `b_x, se_x, b_y, se_y \in R^K`,
#' between-instrument correlations `R_xx, R_yy \in R^{K x K}`, and an
#' across-sample correlation block `R_xy \in R^{K x K}` (default 0 for
#' non-overlapping two-sample MR; index convention
#' `R_xy[i,j] = cor(b_y_i, b_x_j)` — outcome-indexed rows, exposure-
#' indexed cols), let
#' \deqn{m(\beta_0) = b_y - \beta_0 b_x}
#' \deqn{V(\beta_0) = D_y R_{yy} D_y + \beta_0^2 D_x R_{xx} D_x - 2 \beta_0 D_y R_{xy} D_x}
#' with `D_x = diag(se_x)`, `D_y = diag(se_y)`. The AR statistic is
#' \deqn{AR(\beta_0) = m(\beta_0)^T V(\beta_0)^{-1} m(\beta_0) \sim \chi^2_K}
#' under `H_0: beta = beta_0`. The 100(1-alpha)% confidence set is
#' \eqn{\{ \beta_0 : AR(\beta_0) \le q_{\chi^2_K, 1-\alpha} \}}.
#'
#' **Over-identification J-test.** At `\hat\beta = argmin AR`,
#' `J = AR(\hat\beta) ~ chi^2_{K-1}`. Reject homogeneity (i.e.
#' presence of horizontal pleiotropy across the K instruments) at
#' level alpha if `J > qchisq(1-alpha, K-1)`.
#'
#' **Geometry.** `AR(beta_0)` is a rational function of `beta_0`; its
#' level sets are unions of intervals. The implementation grid-
#' evaluates `AR`, locates sign changes of `AR - c_crit`, refines each
#' boundary with [stats::uniroot()], and classifies the set as
#' `bounded_interval`, `disconnected_union`, `whole_line`, or `empty`.
#'
#' @param b_x Numeric length-K. Burden-on-exposure effects.
#' @param se_x Numeric length-K. SE of `b_x`. Must all be > 0.
#' @param b_y Numeric length-K. Burden-on-outcome effects.
#' @param se_y Numeric length-K. SE of `b_y`. Must all be > 0.
#' @param R_xx K x K numeric. Between-instrument correlation of the
#'   exposure estimates. Default `diag(K)` (independence).
#' @param R_yy K x K numeric. Between-instrument correlation of the
#'   outcome estimates. Default `diag(K)`.
#' @param R_xy K x K numeric. Across-sample correlation block. Default
#'   `matrix(0, K, K)` for non-overlapping two-sample MR. Set
#'   non-zero for one-sample / overlapping-sample correction.
#'   **Index convention.** `R_xy[i, j] = cor(b_y_i, b_x_j)`: rows are
#'   indexed by the *outcome* effect components, columns by the
#'   *exposure* effect components. This matches the matrix product
#'   `V_xy = D_y R_xy D_x` used internally and in `main.tex` §Step 9 /
#'   `steps_5_to_9_logic.md`. It is the *transpose* of the more common
#'   joint-covariance layout `Cov([b_x; b_y])[1:K, (K+1):(2K)] =
#'   cor(b_x_i, b_y_j)`; if you estimated `R_xy` from a stacked joint
#'   block in that layout, pass `t(R_xy)` to this function.
#'   **Symmetric-part invariance.** The AR statistic
#'   `m^T V(\beta_0)^{-1} m` is a scalar quadratic form, so it depends
#'   only on the **symmetric part** of `V_xy` (the antisymmetric part
#'   contributes zero to any quadratic form). Consequently
#'   `AR(beta_0; R_xy) == AR(beta_0; t(R_xy))`: transposing R_xy makes
#'   no difference at the AR / J-test layer. The two index conventions
#'   (this package's vs. the joint-covariance layout) are therefore
#'   equivalent for AR inference. Orientation matters only if you
#'   reuse `R_xy` for non-quadratic purposes (e.g. simulating from a
#'   joint distribution). The regression tests in
#'   `tests/testthat/test-mrAR_multi.R` pin (a) the transpose
#'   invariance, (b) that the symmetric part DOES matter, and (c) that
#'   a future edit swapping `D_y R_xy D_x` to `D_x R_xy D_y` with
#'   unequal SEs would be detected.
#' @param alpha Numeric in `(0, 1)`. Significance level (default
#'   0.05).
#' @param n_grid Positive integer. Grid density for the initial
#'   sign-change scan (default 4000L).
#' @param grid_pad_mult Positive numeric. Padding multiplier in
#'   setting the initial grid: range covers
#'   `[min(Wald_k) - pad, max(Wald_k) + pad]` with
#'   `pad = grid_pad_mult * max(|Wald_k|)`.
#' @param grid_extend_max Non-negative integer. Maximum number of
#'   doubling extensions before declaring an unbounded
#'   `whole_line` set (default 3L).
#'
#' @return A `list` with elements:
#'   * `ci_type` — `"bounded_interval"`, `"disconnected_union"`,
#'     `"whole_line"`, or `"empty"`.
#'   * `ci_intervals` — list of length-2 numeric pairs
#'     `c(lower, upper)`. Single pair for `bounded_interval`,
#'     two-or-more for `disconnected_union`, `c(-Inf, Inf)` for
#'     `whole_line`, `list()` for `empty`.
#'   * `beta_hat` — argmin of `AR(beta_0)` (point estimate).
#'   * `J_stat` — `AR(beta_hat)`, the over-id statistic.
#'   * `J_pvalue` — `pchisq(J_stat, df = K - 1, lower.tail = FALSE)`;
#'     `NA` when `K == 1` (no over-id).
#'   * `ar_crit` — chi-square critical value `qchisq(1-alpha, K)`.
#'   * `grid_used` — `c(lo, hi)` final grid bounds.
#'   * `inputs` — list echoing inputs.
#'
#' @references
#' Anderson TW, Rubin H (1949). Estimation of the parameters of a
#' single equation in a complete system of stochastic equations.
#' *Annals of Mathematical Statistics* 20(1): 46-63.
#' doi:10.1214/aoms/1177730090
#'
#' Wang S, Kang H (2022). Weak-instrument-robust tests in two-sample
#' summary-data Mendelian randomization. *Biometrics* 78(4):
#' 1699-1713.
#'
#' Patel A, Lane J, Burgess S (2024). Anderson-Rubin tests for
#' Mendelian randomization with weak and possibly invalid instruments.
#' arXiv:2408.09868.
#'
#' Fieller EC (1954). Some problems in interval estimation.
#' *Journal of the Royal Statistical Society: Series B* 16(2):
#' 175-185. (For the K=1 algebraic equivalence checked by
#' [mrAR()].)
#'
#' @seealso [mrAR()] for the K=1 closed-form sibling.
#'
#' @examples
#' \dontrun{
#'   # K=3 strong-IV example: bounded interval expected, J-test passes.
#'   set.seed(1)
#'   b_x  <- c(0.6, 0.5, 0.4)
#'   se_x <- c(0.1, 0.1, 0.1)
#'   beta <- 0.5
#'   b_y  <- beta * b_x + rnorm(3, sd = 0.02)
#'   se_y <- c(0.05, 0.05, 0.05)
#'   res  <- mrAR_multi(b_x, se_x, b_y, se_y, alpha = 0.05)
#'   res$ci_type           # "bounded_interval"
#'   res$ci_intervals[[1]] # endpoints
#'   res$J_pvalue          # >> 0.05 under homogeneity
#'
#'   # K=3 weak-IV example: disconnected union possible.
#'   res_weak <- mrAR_multi(
#'     b_x  = c(0.05, 0.04, 0.06), se_x = c(0.1, 0.1, 0.1),
#'     b_y  = c(0.30, 0.25, 0.28), se_y = c(0.05, 0.05, 0.05))
#'   res_weak$ci_type
#' }
#' @export
mrAR_multi <- function(b_x, se_x, b_y, se_y,
                       R_xx = diag(length(b_x)),
                       R_yy = diag(length(b_y)),
                       R_xy = matrix(0, length(b_x), length(b_y)),
                       alpha = 0.05,
                       n_grid = 4000L,
                       grid_pad_mult = 3,
                       grid_extend_max = 3L) {
  # ---- input validation -------------------------------------------------
  validate_summary_input(b_x, se_x, b_y, se_y, n_x = 1L, n_y = 1L)
  K <- length(b_x)
  if (length(b_y) != K || length(se_x) != K || length(se_y) != K)
    stop("mrAR_multi: b_x, se_x, b_y, se_y must share length K")

  if (!is.matrix(R_xx) || any(dim(R_xx) != c(K, K)))
    stop("mrAR_multi: R_xx must be a K x K matrix")
  if (!is.matrix(R_yy) || any(dim(R_yy) != c(K, K)))
    stop("mrAR_multi: R_yy must be a K x K matrix")
  if (!is.matrix(R_xy) || any(dim(R_xy) != c(K, K)))
    stop("mrAR_multi: R_xy must be a K x K matrix")
  if (any(!is.finite(R_xx)) || any(!is.finite(R_yy)) || any(!is.finite(R_xy)))
    stop("mrAR_multi: R_xx / R_yy / R_xy contain non-finite values")

  if (!is.numeric(alpha) || length(alpha) != 1L ||
      !is.finite(alpha) || alpha <= 0 || alpha >= 1)
    stop("mrAR_multi: alpha must be a scalar in (0, 1)")
  if (!is.numeric(n_grid) || length(n_grid) != 1L ||
      !is.finite(n_grid) || n_grid < 10)
    stop("mrAR_multi: n_grid must be a positive integer >= 10")
  n_grid <- as.integer(n_grid)
  if (!is.numeric(grid_pad_mult) || length(grid_pad_mult) != 1L ||
      !is.finite(grid_pad_mult) || grid_pad_mult <= 0)
    stop("mrAR_multi: grid_pad_mult must be a positive scalar")
  if (!is.numeric(grid_extend_max) || length(grid_extend_max) != 1L ||
      !is.finite(grid_extend_max) || grid_extend_max < 0)
    stop("mrAR_multi: grid_extend_max must be a non-negative integer")
  grid_extend_max <- as.integer(grid_extend_max)

  # ---- precompute scaled covariance terms -------------------------------
  Dx <- diag(se_x, K, K)
  Dy <- diag(se_y, K, K)
  V_yy <- Dy %*% R_yy %*% Dy           # outcome covariance
  V_xx <- Dx %*% R_xx %*% Dx           # exposure covariance
  # Cross-covariance: V_xy[i,j] = cov(b_y_i, b_x_j) = se_y_i * R_xy[i,j] * se_x_j.
  # R_xy index convention: (outcome row, exposure col). See roxygen.
  # NOTE: AR depends only on the symmetric part of V_xy (quadratic-form
  # invariance), so this assembly is internally robust to a transposed
  # R_xy. The convention is still chosen to match main.tex §Step 9.
  V_xy <- Dy %*% R_xy %*% Dx           # cross-covariance (K x K)

  c_crit <- stats::qchisq(1 - alpha, df = K)

  # ---- AR(beta_0) closure ----------------------------------------------
  ar_fun <- function(b0) {
    m  <- b_y - b0 * b_x
    Vb <- V_yy + (b0 * b0) * V_xx - (2 * b0) * V_xy
    # Solve V b = m via a robust path; fall back to ginv-style on failure.
    sol <- tryCatch(solve(Vb, m), error = function(e) NULL)
    if (is.null(sol)) return(NA_real_)
    val <- as.numeric(crossprod(m, sol))
    if (!is.finite(val)) return(NA_real_)
    val
  }
  ar_vec <- function(b0_vec) vapply(b0_vec, ar_fun, numeric(1))

  # ---- initial grid bounds via per-IV Wald ratios -----------------------
  # Guard against b_x = 0 producing Inf Wald ratios.
  wald_per <- ifelse(abs(b_x) > .Machine$double.eps,
                     b_y / b_x, NA_real_)
  wald_per <- wald_per[is.finite(wald_per)]
  if (length(wald_per) == 0L) {
    centre <- 0
    pad    <- 1
  } else {
    centre_lo <- min(wald_per)
    centre_hi <- max(wald_per)
    pad <- grid_pad_mult * max(abs(wald_per), 1)
    centre <- c(centre_lo, centre_hi)
  }
  grid_lo <- min(centre) - pad
  grid_hi <- max(centre) + pad

  # ---- grid evaluation with up-to grid_extend_max doublings -------------
  build_grid <- function(lo, hi, n) {
    seq.int(lo, hi, length.out = n)
  }

  # We extend grid if both extreme grid points still satisfy AR <= c_crit;
  # that signals the level set leaks past the current envelope (unbounded
  # candidate). We extend a maximum of grid_extend_max times by doubling
  # the half-width on both sides.
  extensions <- 0L
  repeat {
    grid <- build_grid(grid_lo, grid_hi, n_grid)
    ar_grid <- ar_vec(grid)

    # If V(b0) hit a singularity somewhere in the grid, mark those as
    # +Inf so the sign-change logic treats them as "above c_crit".
    ar_grid_finite <- ar_grid
    ar_grid_finite[!is.finite(ar_grid_finite)] <- Inf

    leak_lo <- isTRUE(ar_grid_finite[1L]            <= c_crit)
    leak_hi <- isTRUE(ar_grid_finite[length(grid)]  <= c_crit)
    if ((leak_lo || leak_hi) && extensions < grid_extend_max) {
      span <- grid_hi - grid_lo
      grid_lo <- grid_lo - span
      grid_hi <- grid_hi + span
      extensions <- extensions + 1L
    } else {
      break
    }
  }

  # ---- detect sign changes of AR - c_crit -------------------------------
  diffs <- ar_grid_finite - c_crit
  # sign-change between consecutive grid points
  sc_idx <- which(diffs[-length(diffs)] * diffs[-1L] < 0)

  # ---- refine each boundary with uniroot --------------------------------
  refine <- function(x_lo, x_hi) {
    tryCatch(
      stats::uniroot(function(b) ar_fun(b) - c_crit,
                     lower = x_lo, upper = x_hi,
                     tol = 1e-8, maxiter = 200L)$root,
      error = function(e) NA_real_
    )
  }

  roots <- if (length(sc_idx) > 0L) {
    vapply(sc_idx, function(i) refine(grid[i], grid[i + 1L]),
           numeric(1))
  } else {
    numeric(0)
  }
  roots <- roots[is.finite(roots)]
  roots <- sort(unique(roots))

  # ---- classify CI ------------------------------------------------------
  # Decide acceptance at the two grid extremes and at midpoints between
  # consecutive roots (and beyond the outermost roots, between grid edge
  # and nearest root).
  accepts_at <- function(b0) {
    val <- ar_fun(b0)
    is.finite(val) && val <= c_crit
  }

  n_roots <- length(roots)

  # Helper to glue contiguous accepted segments separated only by roots.
  if (n_roots == 0L) {
    # No boundary crossings detected on grid → either entire grid is
    # accepted (whole_line / unbounded) or entirely rejected (empty).
    if (accepts_at((grid_lo + grid_hi) / 2)) {
      ci_type <- "whole_line"
      ci_intervals <- list(c(-Inf, Inf))
    } else {
      ci_type <- "empty"
      ci_intervals <- list()
    }
  } else {
    # Build candidate endpoints: -Inf, roots..., +Inf, then test interior
    # midpoints for acceptance.
    endpoints <- c(-Inf, roots, Inf)
    seg_accept <- logical(length(endpoints) - 1L)
    for (i in seq_along(seg_accept)) {
      a <- endpoints[i]
      b <- endpoints[i + 1L]
      mid <- if (is.infinite(a) && is.infinite(b)) {
        0
      } else if (is.infinite(a)) {
        b - 1
      } else if (is.infinite(b)) {
        a + 1
      } else {
        (a + b) / 2
      }
      seg_accept[i] <- accepts_at(mid)
    }
    # Collect accepted segments.
    accepted_intervals <- list()
    for (i in seq_along(seg_accept)) {
      if (seg_accept[i]) {
        accepted_intervals[[length(accepted_intervals) + 1L]] <-
          c(endpoints[i], endpoints[i + 1L])
      }
    }
    n_acc <- length(accepted_intervals)
    if (n_acc == 0L) {
      ci_type <- "empty"
      ci_intervals <- list()
    } else if (n_acc == 1L) {
      iv <- accepted_intervals[[1L]]
      if (is.infinite(iv[1]) && is.infinite(iv[2])) {
        ci_type <- "whole_line"
      } else if (is.infinite(iv[1]) || is.infinite(iv[2])) {
        # Single half-line. If grid extensions are exhausted and one side
        # never crosses, treat as whole_line (in practice for K=1 weak
        # IV this matches the disconnected case but here we have only
        # one root, so it really is a half-line "bounded_interval" with
        # one infinite endpoint — keep as "bounded_interval" with the
        # corresponding endpoint = +-Inf).
        ci_type <- "bounded_interval"
      } else {
        ci_type <- "bounded_interval"
      }
      ci_intervals <- accepted_intervals
    } else {
      ci_type <- "disconnected_union"
      ci_intervals <- accepted_intervals
    }
  }

  # ---- point estimate via optimize over the (possibly extended) grid ----
  opt <- tryCatch(
    stats::optimize(ar_fun, lower = grid_lo, upper = grid_hi,
                    tol = 1e-8),
    error = function(e) NULL
  )
  if (is.null(opt) || !is.finite(opt$objective)) {
    # Fallback: argmin from the dense grid.
    j <- which.min(ar_grid_finite)
    beta_hat <- grid[j]
    J_stat   <- ar_grid_finite[j]
  } else {
    beta_hat <- opt$minimum
    J_stat   <- opt$objective
  }

  J_pvalue <- if (K >= 2L) {
    stats::pchisq(J_stat, df = K - 1L, lower.tail = FALSE)
  } else {
    NA_real_
  }

  list(
    ci_type      = ci_type,
    ci_intervals = ci_intervals,
    beta_hat     = beta_hat,
    J_stat       = J_stat,
    J_pvalue     = J_pvalue,
    ar_crit      = c_crit,
    grid_used    = c(grid_lo, grid_hi),
    inputs       = list(b_x = b_x, se_x = se_x,
                        b_y = b_y, se_y = se_y,
                        R_xx = R_xx, R_yy = R_yy, R_xy = R_xy,
                        alpha = alpha, K = K,
                        n_grid = n_grid,
                        grid_pad_mult = grid_pad_mult,
                        grid_extend_max = grid_extend_max,
                        grid_extensions_used = extensions)
  )
}
