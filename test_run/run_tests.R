# =====================================================================
# run_tests.R
#
# Phase-3 driver: loads rvMR (devtools::load_all), the generator
# (generate_test_data.R) and the scenario grid (scenarios.R); loops
# over the 6 scenarios x n_reps replicates; collects per-replicate
# AR/J results plus a naive-Wald comparator; writes results.rds and
# results.md.
# =====================================================================

suppressMessages(suppressWarnings(devtools::load_all(
  "/home/francisfenglu4/rvSMR/May_30md/rvMR",
  quiet = TRUE)))

source("generate_test_data.R")
source("scenarios.R")

# ---- helpers --------------------------------------------------------

# Determine whether a value is inside the AR confidence set (handles all
# 4 mrAR_multi shapes).
in_AR_CI <- function(beta, ci_type, ci_intervals) {
  if (ci_type == "empty") return(FALSE)
  if (ci_type == "whole_line") return(TRUE)
  for (iv in ci_intervals) {
    if (beta >= iv[1] && beta <= iv[2]) return(TRUE)
  }
  FALSE
}

# Length of the (finite portion of the) accepted set.  Inf for
# whole_line, sum of finite-interval lengths otherwise.
ci_length <- function(ci_type, ci_intervals) {
  if (ci_type == "empty") return(0)
  if (ci_type == "whole_line") return(Inf)
  total <- 0
  for (iv in ci_intervals) {
    if (is.finite(iv[1]) && is.finite(iv[2])) {
      total <- total + (iv[2] - iv[1])
    } else {
      total <- Inf
    }
  }
  total
}

# Naive comparator: inverse-variance-weighted pooled burden + Wald CI.
#   pooled_b_x = sum(b_x_k / se_x_k^2) / sum(1 / se_x_k^2)
#   pooled_b_y = sum(b_y_k / se_y_k^2) / sum(1 / se_y_k^2)
#   beta_naive = pooled_b_y / pooled_b_x
#   se_naive   = delta method
#   95% Wald CI = beta_naive +/- 1.96 * se_naive
ivw_pooled_wald <- function(d) {
  w_x <- 1 / d$se_x^2
  w_y <- 1 / d$se_y^2
  bx_p <- sum(d$b_x * w_x) / sum(w_x)
  by_p <- sum(d$b_y * w_y) / sum(w_y)
  sx_p <- sqrt(1 / sum(w_x))
  sy_p <- sqrt(1 / sum(w_y))
  beta_n <- by_p / bx_p
  se_n   <- sqrt((sy_p / bx_p)^2 + (by_p^2 * sx_p^2) / (bx_p^4))
  c(beta = beta_n, se = se_n,
    lo = beta_n - 1.96 * se_n, hi = beta_n + 1.96 * se_n,
    F = (bx_p / sx_p)^2)
}

# ---- main loop ------------------------------------------------------

run_one_scenario <- function(scn, n_reps, master_seed) {
  # Pre-allocate results
  res <- data.frame(
    rep            = integer(n_reps),
    ci_type        = character(n_reps),
    ci_covers_beta = logical(n_reps),
    ci_covers_zero = logical(n_reps),
    ci_len_finite  = numeric(n_reps),
    beta_hat       = numeric(n_reps),
    bias           = numeric(n_reps),
    J_stat         = numeric(n_reps),
    J_pvalue       = numeric(n_reps),
    F_mean         = numeric(n_reps),
    F_max          = numeric(n_reps),
    naive_beta     = numeric(n_reps),
    naive_se       = numeric(n_reps),
    naive_covers_beta = logical(n_reps),
    naive_covers_zero = logical(n_reps),
    naive_len      = numeric(n_reps),
    stringsAsFactors = FALSE
  )

  beta_true <- scn$beta_true
  K         <- scn$K

  # Build the R_xy matrix as a DIAGONAL block (within-mask correlation
  # only; cross-mask cross-sample correlations are zero because masks are
  # disjoint variant sets -- see steps_5_to_9_logic.md §"defaults
  # legality" and HANDOVER §5).  This matches what simulate_burden_mr()
  # actually draws (which uses only R_xy[k,k]) and is the rvSMR
  # convention for the K disjoint annotation classes.
  R_xy <- if (scn$R_xy_block == 0) {
    matrix(0, K, K)
  } else {
    diag(scn$R_xy_block, K, K)
  }

  for (r in seq_len(n_reps)) {
    # Use modular seed to avoid integer overflow (R int is 32-bit).
    seed_r <- ((as.numeric(master_seed) %% 1000000L) * 1000 + r) %% .Machine$integer.max
    d <- simulate_burden_mr(
      K               = K,
      beta_true       = beta_true,
      F_target        = scn$F_target,
      n_x             = scn$n_x,
      n_y             = scn$n_y,
      pleiotropy_frac = scn$pleiotropy_frac,
      R_xy            = R_xy,
      seed            = seed_r
    )

    ar <- tryCatch(
      mrAR_multi(b_x = d$b_x, se_x = d$se_x,
                 b_y = d$b_y, se_y = d$se_y,
                 R_xx = d$R_xx, R_yy = d$R_yy,
                 R_xy = d$R_xy,
                 alpha = 0.05),
      error = function(e) NULL
    )

    if (is.null(ar)) {
      res$rep[r]    <- r
      res$ci_type[r] <- "error"
      next
    }

    res$rep[r]            <- r
    res$ci_type[r]        <- ar$ci_type
    res$ci_covers_beta[r] <- in_AR_CI(beta_true, ar$ci_type, ar$ci_intervals)
    res$ci_covers_zero[r] <- in_AR_CI(0,         ar$ci_type, ar$ci_intervals)
    res$ci_len_finite[r]  <- ci_length(ar$ci_type, ar$ci_intervals)
    res$beta_hat[r]       <- ar$beta_hat
    res$bias[r]           <- ar$beta_hat - beta_true
    res$J_stat[r]         <- ar$J_stat
    res$J_pvalue[r]       <- ar$J_pvalue
    res$F_mean[r]         <- mean((d$b_x / d$se_x)^2)
    res$F_max[r]          <- max((d$b_x / d$se_x)^2)

    # naive comparator
    nv <- ivw_pooled_wald(d)
    res$naive_beta[r]         <- nv["beta"]
    res$naive_se[r]           <- nv["se"]
    res$naive_covers_beta[r]  <- nv["lo"] <= beta_true && beta_true <= nv["hi"]
    res$naive_covers_zero[r]  <- nv["lo"] <= 0         && 0         <= nv["hi"]
    res$naive_len[r]          <- nv["hi"] - nv["lo"]
  }

  res
}

summarize_scenario <- function(scn, df) {
  ok <- df$ci_type != "error"
  n  <- sum(ok)
  cov_ar      <- mean(df$ci_covers_beta[ok])
  cov_naive   <- mean(df$naive_covers_beta[ok])
  rej_zero    <- mean(!df$ci_covers_zero[ok])
  rej_zero_nv <- mean(!df$naive_covers_zero[ok])
  F_mean      <- mean(df$F_mean[ok])
  bias        <- mean(df$bias[ok], na.rm = TRUE)
  bias_med    <- median(df$bias[ok], na.rm = TRUE)
  # J p-value Type-I/power: under H0 (no pleiotropy) we report Type-I,
  # under H1 (pleiotropy) we report power.  Same number, interpretation
  # flips.
  J_rej_rate  <- mean(df$J_pvalue[ok] < 0.05, na.rm = TRUE)
  shape_tab   <- prop.table(table(df$ci_type[ok]))

  list(
    label             = scn$label,
    n_ok              = n,
    coverage_AR       = cov_ar,
    coverage_naive    = cov_naive,
    rej_zero_AR       = rej_zero,
    rej_zero_naive    = rej_zero_nv,
    F_mean            = F_mean,
    bias_mean         = bias,
    bias_median       = bias_med,
    J_rej_rate        = J_rej_rate,
    shape_distribution = shape_tab
  )
}

# ---- run all scenarios ----------------------------------------------

run_all <- function(n_reps = n_reps_default, master_seed = 20260603L,
                    verbose = TRUE) {
  results <- list()
  summaries <- list()
  t0 <- Sys.time()
  for (nm in names(scenarios)) {
    scn <- scenarios[[nm]]
    if (verbose) {
      cat(sprintf("\n=== %s (n=%d reps) ===\n", scn$label, n_reps))
      tic <- Sys.time()
    }
    df  <- run_one_scenario(scn, n_reps, master_seed + match(nm, names(scenarios)))
    sm  <- summarize_scenario(scn, df)
    results[[nm]]   <- df
    summaries[[nm]] <- sm
    if (verbose) {
      toc <- Sys.time()
      cat(sprintf("  duration: %.1fs\n", as.numeric(toc - tic, units = "secs")))
      cat(sprintf("  Coverage AR = %.3f | Coverage naive = %.3f | F_mean = %.2f\n",
                  sm$coverage_AR, sm$coverage_naive, sm$F_mean))
      cat(sprintf("  Rej beta=0 (AR / naive) = %.3f / %.3f\n",
                  sm$rej_zero_AR, sm$rej_zero_naive))
      cat(sprintf("  J p<0.05 rate = %.3f\n", sm$J_rej_rate))
      cat(sprintf("  Bias mean / median = %+.4f / %+.4f\n",
                  sm$bias_mean, sm$bias_median))
      cat("  Shape: ")
      print(sm$shape_distribution)
    }
  }
  t1 <- Sys.time()
  if (verbose) {
    cat(sprintf("\nTotal duration: %.1fs\n",
                as.numeric(t1 - t0, units = "secs")))
  }
  list(results = results, summaries = summaries,
       n_reps = n_reps, master_seed = master_seed,
       duration_sec = as.numeric(t1 - t0, units = "secs"))
}

# ---- markdown table writer ------------------------------------------

write_results_md <- function(R, path) {
  con <- file(path, open = "w")
  on.exit(close(con))
  cat("# rvMR Phase-3 simulation results\n\n", file = con)
  cat(sprintf("- Replicates per scenario: **%d**\n", R$n_reps), file = con)
  cat(sprintf("- Master seed: %d\n", R$master_seed), file = con)
  cat(sprintf("- Total wall time: %.1fs\n", R$duration_sec), file = con)
  cat(sprintf("- R version: %s\n", R.version.string), file = con)
  cat(sprintf("- rvMR loaded from /home/francisfenglu4/rvSMR/May_30md/rvMR\n\n"), file = con)

  cat("## Headline coverage table\n\n", file = con)
  cat("| Scenario | n_ok | AR cov(β) | naive cov(β) | mean F | bias β̂ | rej β=0 (AR) | J<0.05 rate |\n",
      file = con)
  cat("|---|---:|---:|---:|---:|---:|---:|---:|\n", file = con)
  for (nm in names(R$summaries)) {
    s <- R$summaries[[nm]]
    cat(sprintf("| %s | %d | %.3f | %.3f | %.2f | %+.4f | %.3f | %.3f |\n",
                s$label, s$n_ok, s$coverage_AR, s$coverage_naive,
                s$F_mean, s$bias_mean, s$rej_zero_AR, s$J_rej_rate),
        file = con)
  }
  cat("\n", file = con)

  cat("## CI shape distribution per scenario\n\n", file = con)
  for (nm in names(R$summaries)) {
    s <- R$summaries[[nm]]
    cat(sprintf("### %s\n\n", s$label), file = con)
    cat("| CI shape | proportion |\n|---|---:|\n", file = con)
    sh <- s$shape_distribution
    for (k in names(sh)) {
      cat(sprintf("| %s | %.3f |\n", k, as.numeric(sh[k])), file = con)
    }
    cat("\n", file = con)
  }

  cat("## Pass / fail interpretation\n\n", file = con)
  for (nm in names(R$summaries)) {
    s <- R$summaries[[nm]]
    cat(sprintf("- **%s**: coverage_AR = %.3f, bias = %+.4f, J-rej = %.3f. ",
                s$label, s$coverage_AR, s$bias_mean, s$J_rej_rate),
        file = con)
    cat(interpret_one(nm, s), file = con)
    cat("\n", file = con)
  }
}

interpret_one <- function(nm, s) {
  switch(
    nm,
    A_null      = sprintf("Coverage of β=0 should be ≈ 0.95 (PASS if in [0.93, 0.97]). Type-I of rejecting β=0 is %.3f. %s",
                          s$rej_zero_AR,
                          if (abs(s$coverage_AR - 0.95) < 0.02) "**PASS**." else "**REVIEW** (coverage out of tolerance)."),
    B_strong    = sprintf("Coverage should be ≈ 0.95; bias should be ~0. %s",
                          if (abs(s$coverage_AR - 0.95) < 0.02 && abs(s$bias_mean) < 0.05) "**PASS**." else "**REVIEW**."),
    C_weak      = sprintf("AR's headline regime: coverage **must hold** at F≈1. %s",
                          if (abs(s$coverage_AR - 0.95) < 0.03) "**PASS** -- AR is doing its job." else "**FAIL** -- AR coverage degraded at low F."),
    D_very_weak = sprintf("F<1: expect many whole_line CIs, but coverage should still ≥ 0.95 (conservative). %s",
                          if (s$coverage_AR >= 0.93) "**PASS** (no under-coverage)." else "**FAIL**."),
    E_pleiotropy = sprintf("AR no longer valid (model misspecified); but Sargan-J should reject (high J p<0.05 rate). %s",
                          if (s$J_rej_rate > 0.5) "**PASS** (J detects pleiotropy)." else "**REVIEW** (J underpowered here)."),
    F_overlap   = sprintf("Coverage should remain ≈ 0.95 when R_xy is supplied. %s",
                          if (abs(s$coverage_AR - 0.95) < 0.03) "**PASS**." else "**REVIEW**.")
  )
}

# ---- entry point ----------------------------------------------------

if (!interactive()) {
  N <- as.integer(Sys.getenv("RVMR_NREPS", n_reps_default))
  cat(sprintf("Running run_all() with n_reps=%d...\n", N))
  R <- run_all(n_reps = N)
  saveRDS(R, file = "results.rds")
  write_results_md(R, "results.md")
  cat("\nWrote results.rds and results.md.\n")
}
