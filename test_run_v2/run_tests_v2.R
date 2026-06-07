# =====================================================================
# run_tests_v2.R  (Round 2)
#
# Round 2 driver. Fixes from CRITIQUE:
#   S1.1/S1.2: replace pool-then-divide naive Wald with the *classical*
#              IVW-of-ratios (per-mask Wald ratio, delta-method SE,
#              inverse-variance weighted mean, Wald CI).
#   S1.6:      report bias *conditional* on ci_type == "bounded_interval"
#              only; add n_bounded column. Whole-line / disconnected /
#              empty reps have undefined beta-hat.
#   S2.4:      multi-seed support (3 master seeds), report coverage +/- SE.
#   S1.5:      do NOT fabricate "steps_5_to_9_logic.md" citation; for
#              overlap scenario, pass the TRUE block R_xx, R_yy, R_xy
#              from the generator into mrAR_multi (Burgess 2016).
#
# Inputs:
#   - generate_test_data_v2.R    (simulate_burden_mr_v2)
#   - scenarios_v2.R             (scenarios_v2, master_seeds_v2,
#                                 n_reps_default_v2)
# =====================================================================

suppressMessages(suppressWarnings(devtools::load_all(
  "/home/francisfenglu4/rvSMR/May_30md/rvMR",
  quiet = TRUE)))

source("generate_test_data_v2.R")
source("scenarios_v2.R")

# ---- helpers --------------------------------------------------------

in_AR_CI <- function(beta, ci_type, ci_intervals) {
  if (ci_type == "empty") return(FALSE)
  if (ci_type == "whole_line") return(TRUE)
  for (iv in ci_intervals) {
    if (beta >= iv[1] && beta <= iv[2]) return(TRUE)
  }
  FALSE
}

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

# Classical IVW-of-ratios comparator (CRITIQUE S1.1).
# Per-mask Wald ratio + delta-method SE + inverse-variance weighted mean.
# Note: at weak IV the delta-method SE inflates the denominator (b_x^2 in
# v_k), giving WIDER CIs -- so this estimator can over-cover, NOT
# under-cover (this is the partial weak-IV-buffering that CRITIQUE S3.1
# warned about). We report its coverage honestly.
ivw_of_ratios <- function(d) {
  bx <- d$b_x; sx <- d$se_x
  by <- d$b_y; sy <- d$se_y
  ok <- abs(bx) > .Machine$double.eps^0.5
  if (!any(ok)) {
    return(c(beta = NA_real_, se = NA_real_,
             lo = NA_real_, hi = NA_real_, F = 0,
             n_used = 0))
  }
  bx <- bx[ok]; sx <- sx[ok]; by <- by[ok]; sy <- sy[ok]
  r  <- by / bx
  v  <- (sy / bx)^2 + (by * sx / bx^2)^2   # delta-method Var(r_k)
  w  <- 1 / v
  beta_h <- sum(w * r) / sum(w)
  se_h   <- sqrt(1 / sum(w))
  c(beta = beta_h, se = se_h,
    lo = beta_h - 1.96 * se_h, hi = beta_h + 1.96 * se_h,
    F = mean((bx / sx)^2),
    n_used = length(bx))
}

# TSLS-style summary comparator (the *canonical* weak-IV-collapse form).
# beta_hat = (b_x' W b_x)^{-1} b_x' W b_y, W = diag(1/se_y^2).
# SE drops the delta-method `b_x^4` denominator: SE_TSLS = (b_x' W b_x)^{-1/2}.
# This is what Wang-Kang Fig 6 / PLB Fig 2 typically show collapsing at
# low F: SE stays "small" while beta_hat is heavy-tailed via the
# `(b_x'Wb_x)^{-1}` factor near zero. The summary-stat TSLS is the
# scalar IVW from MR-Base / Burgess MendelianRandomization::mr_ivw().
tsls_summary <- function(d) {
  bx <- d$b_x; sx <- d$se_x
  by <- d$b_y; sy <- d$se_y
  W  <- 1 / sy^2
  bxWbx <- sum(bx * W * bx)
  if (!is.finite(bxWbx) || bxWbx <= 0) {
    return(c(beta = NA_real_, se = NA_real_,
             lo = NA_real_, hi = NA_real_, F = 0))
  }
  beta_h <- sum(bx * W * by) / bxWbx
  se_h   <- 1 / sqrt(bxWbx)
  c(beta = beta_h, se = se_h,
    lo = beta_h - 1.96 * se_h, hi = beta_h + 1.96 * se_h,
    F = mean((bx / sx)^2))
}

# ---- main per-scenario per-seed runner ------------------------------

run_one_scenario_seed <- function(scn, n_reps, master_seed,
                                  scenario_idx_offset) {
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
    naive_lo       = numeric(n_reps),
    naive_hi       = numeric(n_reps),
    naive_covers_beta = logical(n_reps),
    naive_covers_zero = logical(n_reps),
    naive_len      = numeric(n_reps),
    naive_n_used   = integer(n_reps),
    tsls_beta      = numeric(n_reps),
    tsls_se        = numeric(n_reps),
    tsls_lo        = numeric(n_reps),
    tsls_hi        = numeric(n_reps),
    tsls_covers_beta = logical(n_reps),
    tsls_covers_zero = logical(n_reps),
    tsls_len       = numeric(n_reps),
    stringsAsFactors = FALSE
  )

  beta_true <- scn$beta_true
  K         <- scn$K

  # Decide pleio_size from pleio_size_mult * SE_y if relevant.
  se_y_proxy <- 1 / sqrt(scn$n_y)
  pleio_size <- if (!is.null(scn$pleio_size_mult)) {
    scn$pleio_size_mult * se_y_proxy
  } else {
    NULL
  }

  for (r in seq_len(n_reps)) {
    seed_r <- ((as.numeric(master_seed) + scenario_idx_offset) %% 1e6) * 10000 + r
    seed_r <- seed_r %% .Machine$integer.max

    d <- simulate_burden_mr_v2(
      K               = K,
      beta_true       = beta_true,
      F_target        = scn$F_target,
      n_x             = scn$n_x,
      n_y             = scn$n_y,
      pleiotropy_frac = scn$pleiotropy_frac,
      pleio_size      = pleio_size,
      dgp             = scn$dgp,
      conf_strength   = scn$conf_strength,
      rho_xx          = scn$rho_xx,
      rho_xy_diag     = scn$rho_xy_diag,
      seed            = seed_r
    )

    # Inference covariance choice:
    #   - dgp == "overlap": pass true block R_xx, R_yy, R_xy from generator.
    #   - dgp == "ld_xx":   pass true R_xx; R_yy = I; R_xy = 0.
    #   - dgp == "confounder": pass DEFAULT R_xx = R_yy = I, R_xy = 0
    #     (this is the canonical confounder stress test: user doesn't
    #     know about u).
    #   - dgp == "ar": pass diag R_xy (if rho_xy_diag != 0) else 0;
    #     R_xx = R_yy = I.
    if (scn$dgp == "overlap") {
      R_xx_inf <- d$R_xx; R_yy_inf <- d$R_yy; R_xy_inf <- d$R_xy
    } else if (scn$dgp == "ld_xx") {
      R_xx_inf <- d$R_xx; R_yy_inf <- diag(K); R_xy_inf <- matrix(0, K, K)
    } else if (scn$dgp == "confounder") {
      R_xx_inf <- diag(K); R_yy_inf <- diag(K); R_xy_inf <- matrix(0, K, K)
    } else {
      # "ar"
      R_xx_inf <- diag(K); R_yy_inf <- diag(K)
      R_xy_inf <- if (scn$rho_xy_diag != 0) diag(scn$rho_xy_diag, K, K) else matrix(0, K, K)
    }

    ar <- tryCatch(
      mrAR_multi(b_x = d$b_x, se_x = d$se_x,
                 b_y = d$b_y, se_y = d$se_y,
                 R_xx = R_xx_inf, R_yy = R_yy_inf,
                 R_xy = R_xy_inf,
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
    # CRITIQUE S1.6: bias only meaningful when CI is bounded_interval.
    res$bias[r]           <- if (ar$ci_type == "bounded_interval") {
      ar$beta_hat - beta_true
    } else {
      NA_real_
    }
    res$J_stat[r]         <- ar$J_stat
    res$J_pvalue[r]       <- ar$J_pvalue
    res$F_mean[r]         <- mean((d$b_x / d$se_x)^2)
    res$F_max[r]          <- max((d$b_x / d$se_x)^2)

    nv <- ivw_of_ratios(d)
    res$naive_beta[r]         <- nv["beta"]
    res$naive_se[r]           <- nv["se"]
    res$naive_lo[r]           <- nv["lo"]
    res$naive_hi[r]           <- nv["hi"]
    res$naive_covers_beta[r]  <- isTRUE(nv["lo"] <= beta_true && beta_true <= nv["hi"])
    res$naive_covers_zero[r]  <- isTRUE(nv["lo"] <= 0         && 0         <= nv["hi"])
    res$naive_len[r]          <- nv["hi"] - nv["lo"]
    res$naive_n_used[r]       <- as.integer(nv["n_used"])

    ts <- tsls_summary(d)
    res$tsls_beta[r]          <- ts["beta"]
    res$tsls_se[r]            <- ts["se"]
    res$tsls_lo[r]            <- ts["lo"]
    res$tsls_hi[r]            <- ts["hi"]
    res$tsls_covers_beta[r]   <- isTRUE(ts["lo"] <= beta_true && beta_true <= ts["hi"])
    res$tsls_covers_zero[r]   <- isTRUE(ts["lo"] <= 0         && 0         <= ts["hi"])
    res$tsls_len[r]           <- ts["hi"] - ts["lo"]
  }

  res
}

# ---- per-scenario aggregator (across reps in a single seed) ---------

summarize_one_seed <- function(scn, df) {
  ok <- df$ci_type != "error"
  n  <- sum(ok)
  cov_ar      <- mean(df$ci_covers_beta[ok])
  cov_naive   <- mean(df$naive_covers_beta[ok], na.rm = TRUE)
  cov_tsls    <- mean(df$tsls_covers_beta[ok], na.rm = TRUE)
  rej_zero    <- mean(!df$ci_covers_zero[ok])
  rej_zero_nv <- mean(!df$naive_covers_zero[ok], na.rm = TRUE)
  rej_zero_tsls <- mean(!df$tsls_covers_zero[ok], na.rm = TRUE)
  F_mean      <- mean(df$F_mean[ok])
  bounded_mask <- ok & (df$ci_type == "bounded_interval")
  n_bounded   <- sum(bounded_mask)
  bias_mean   <- if (n_bounded > 0) mean(df$bias[bounded_mask]) else NA_real_
  bias_median <- if (n_bounded > 0) median(df$bias[bounded_mask]) else NA_real_
  J_rej_rate  <- mean(df$J_pvalue[ok] < 0.05, na.rm = TRUE)
  shape_tab   <- prop.table(table(df$ci_type[ok]))

  list(
    label             = scn$label,
    n_ok              = n,
    n_bounded         = n_bounded,
    coverage_AR       = cov_ar,
    coverage_naive    = cov_naive,
    coverage_tsls     = cov_tsls,
    rej_zero_AR       = rej_zero,
    rej_zero_naive    = rej_zero_nv,
    rej_zero_tsls     = rej_zero_tsls,
    F_mean            = F_mean,
    bias_mean         = bias_mean,
    bias_median       = bias_median,
    J_rej_rate        = J_rej_rate,
    shape_distribution = shape_tab
  )
}

# ---- multi-seed aggregator -----------------------------------------

aggregate_seeds <- function(per_seed) {
  cov_ar   <- sapply(per_seed, function(s) s$coverage_AR)
  cov_nv   <- sapply(per_seed, function(s) s$coverage_naive)
  cov_ts   <- sapply(per_seed, function(s) s$coverage_tsls)
  rej_zero <- sapply(per_seed, function(s) s$rej_zero_AR)
  rej_zero_nv <- sapply(per_seed, function(s) s$rej_zero_naive)
  rej_zero_ts <- sapply(per_seed, function(s) s$rej_zero_tsls)
  Fm       <- sapply(per_seed, function(s) s$F_mean)
  J_rej    <- sapply(per_seed, function(s) s$J_rej_rate)
  n_bd     <- sapply(per_seed, function(s) s$n_bounded)
  bias_bd  <- sapply(per_seed, function(s) s$bias_mean)
  list(
    coverage_AR_mean     = mean(cov_ar),
    coverage_AR_se       = sd(cov_ar) / sqrt(length(cov_ar)),
    coverage_AR_per_seed = cov_ar,
    coverage_naive_mean  = mean(cov_nv, na.rm = TRUE),
    coverage_naive_se    = sd(cov_nv, na.rm = TRUE) / sqrt(length(cov_nv)),
    coverage_naive_per_seed = cov_nv,
    coverage_tsls_mean   = mean(cov_ts, na.rm = TRUE),
    coverage_tsls_se     = sd(cov_ts, na.rm = TRUE) / sqrt(length(cov_ts)),
    coverage_tsls_per_seed = cov_ts,
    rej_zero_AR_mean     = mean(rej_zero),
    rej_zero_naive_mean  = mean(rej_zero_nv, na.rm = TRUE),
    rej_zero_tsls_mean   = mean(rej_zero_ts, na.rm = TRUE),
    F_mean               = mean(Fm),
    J_rej_rate_mean      = mean(J_rej),
    n_bounded_mean       = mean(n_bd),
    bias_mean_bounded    = mean(bias_bd, na.rm = TRUE)
  )
}

# ---- main loop ------------------------------------------------------

run_all_v2 <- function(n_reps = n_reps_default_v2,
                       master_seeds = master_seeds_v2,
                       verbose = TRUE) {
  results   <- list()  # results[[scenario]][[seed_idx]] = data.frame
  summaries <- list()  # summaries[[scenario]] = aggregated across seeds
  per_seed  <- list()  # per_seed[[scenario]][[seed_idx]] = one-seed summary
  t0 <- Sys.time()
  for (i in seq_along(scenarios_v2)) {
    nm  <- names(scenarios_v2)[i]
    scn <- scenarios_v2[[nm]]
    if (verbose) {
      cat(sprintf("\n=== %s ===\n", scn$label))
    }
    results[[nm]] <- list()
    per_seed[[nm]] <- list()
    for (s_idx in seq_along(master_seeds)) {
      ms <- master_seeds[s_idx]
      if (verbose) cat(sprintf("  seed %d (idx %d) ... ", ms, s_idx))
      tic <- Sys.time()
      df <- run_one_scenario_seed(scn, n_reps = n_reps,
                                  master_seed = ms,
                                  scenario_idx_offset = i)
      sm <- summarize_one_seed(scn, df)
      results[[nm]][[s_idx]]  <- df
      per_seed[[nm]][[s_idx]] <- sm
      toc <- Sys.time()
      if (verbose) {
        cat(sprintf("done (%.0fs): cov_AR=%.3f cov_IVW=%.3f cov_TSLS=%.3f F_mean=%.2f n_bd=%d\n",
                    as.numeric(toc - tic, units = "secs"),
                    sm$coverage_AR, sm$coverage_naive, sm$coverage_tsls,
                    sm$F_mean, sm$n_bounded))
      }
    }
    summaries[[nm]] <- aggregate_seeds(per_seed[[nm]])
  }
  t1 <- Sys.time()
  if (verbose) {
    cat(sprintf("\nTotal duration: %.1fs\n",
                as.numeric(t1 - t0, units = "secs")))
  }
  list(results = results, summaries = summaries, per_seed = per_seed,
       n_reps = n_reps, master_seeds = master_seeds,
       duration_sec = as.numeric(t1 - t0, units = "secs"))
}

# ---- markdown writer (results_v2.md) -------------------------------

write_results_v2_md <- function(R, path) {
  con <- file(path, open = "w")
  on.exit(close(con))
  cat("# rvMR Round-2 simulation results\n\n", file = con)
  cat(sprintf("- Replicates per scenario per seed: **%d**\n", R$n_reps), file = con)
  cat(sprintf("- Master seeds: %s\n",
              paste(R$master_seeds, collapse = ", ")), file = con)
  cat(sprintf("- Total wall time: %.1fs\n", R$duration_sec), file = con)
  cat(sprintf("- R version: %s\n", R.version.string), file = con)
  cat(sprintf("- rvMR loaded from /home/francisfenglu4/rvSMR/May_30md/rvMR\n\n"), file = con)

  # --- F sweep table (headline) ---
  cat("## Headline: coverage vs lambda (F-sweep)\n\n", file = con)
  cat("`lambda_joint = K * F_target` at K=3, R_xx=I (CRITIQUE S1.4).\n", file = con)
  cat("Two comparators are reported:\n", file = con)
  cat("  - **IVW-of-ratios** (CRITIQUE S1.1 canonical form): per-mask Wald ratio\n", file = con)
  cat("    + delta-method SE + IVW. Partially weak-IV-buffered via SE inflation.\n", file = con)
  cat("  - **TSLS** (summary form): (b_x' W b_x)^{-1} b_x' W b_y, W=diag(1/SE_y^2).\n", file = con)
  cat("    Equivalent to scalar IVW from `MendelianRandomization::mr_ivw`.\n", file = con)
  cat("    This is the comparator that Wang-Kang Fig 6 reports collapsing.\n\n", file = con)
  cat("| F | lambda | AR cov (mean +/- SE) | IVW cov (mean +/- SE) | TSLS cov (mean +/- SE) | F_mean | n_bounded | bias|bounded |\n", file = con)
  cat("|---:|---:|---|---|---|---:|---:|---|\n", file = con)
  for (nm in names(R$summaries)) {
    if (!grepl("^Fsweep_", nm)) next
    s <- R$summaries[[nm]]
    scn <- scenarios_v2[[nm]]
    cat(sprintf("| %g | %.1f | %.3f +/- %.3f | %.3f +/- %.3f | %.3f +/- %.3f | %.2f | %.0f | %s |\n",
                scn$F_target, 3 * scn$F_target,
                s$coverage_AR_mean, s$coverage_AR_se,
                s$coverage_naive_mean, s$coverage_naive_se,
                s$coverage_tsls_mean, s$coverage_tsls_se,
                s$F_mean, s$n_bounded_mean,
                if (is.finite(s$bias_mean_bounded))
                  sprintf("%+.4f", s$bias_mean_bounded) else "NA"),
        file = con)
  }
  cat("\n", file = con)

  # --- Pleio sweep ---
  cat("## Pleiotropy magnitude sweep (at F=20, 1/3 invalid)\n\n", file = con)
  cat("| pleio mult | AR cov mean | AR cov SE | J<0.05 rate | F_mean |\n", file = con)
  cat("|---:|---:|---:|---:|---:|\n", file = con)
  for (nm in names(R$summaries)) {
    if (!grepl("^Plei_mult", nm)) next
    s <- R$summaries[[nm]]
    scn <- scenarios_v2[[nm]]
    cat(sprintf("| %g | %.3f | %.4f | %.3f | %.2f |\n",
                scn$pleio_size_mult,
                s$coverage_AR_mean, s$coverage_AR_se,
                s$J_rej_rate_mean, s$F_mean),
        file = con)
  }
  cat("\n", file = con)

  # --- Other scenarios (anchors, confounder, LD, overlap) ---
  cat("## Anchor / confounder / LD / overlap cells\n\n", file = con)
  cat("| Scenario | AR cov | IVW cov | TSLS cov | rej_zero (AR) | J<0.05 | F_mean | n_bd |\n", file = con)
  cat("|---|---:|---:|---:|---:|---:|---:|---:|\n", file = con)
  for (nm in names(R$summaries)) {
    if (grepl("^Fsweep_", nm) || grepl("^Plei_mult", nm)) next
    s <- R$summaries[[nm]]
    scn <- scenarios_v2[[nm]]
    cat(sprintf("| %s | %.3f +/- %.3f | %.3f | %.3f | %.3f | %.3f | %.2f | %.0f |\n",
                scn$label,
                s$coverage_AR_mean, s$coverage_AR_se,
                s$coverage_naive_mean,
                s$coverage_tsls_mean,
                s$rej_zero_AR_mean,
                s$J_rej_rate_mean, s$F_mean,
                s$n_bounded_mean),
        file = con)
  }
  cat("\n", file = con)

  # --- Per-seed coverage breakdown ---
  cat("## Per-seed AR coverage (each row = scenario; columns = 3 seeds)\n\n", file = con)
  cat("| Scenario | seed1 | seed2 | seed3 |\n|---|---:|---:|---:|\n", file = con)
  for (nm in names(R$summaries)) {
    s <- R$summaries[[nm]]
    cps <- s$coverage_AR_per_seed
    cat(sprintf("| %s | %.3f | %.3f | %.3f |\n",
                scenarios_v2[[nm]]$label,
                cps[1], cps[2], cps[3]),
        file = con)
  }
  cat("\n", file = con)
}

# ---- entry point ----------------------------------------------------

if (!interactive() && identical(Sys.getenv("RVMR_RUN_MAIN"), "1")) {
  N <- as.integer(Sys.getenv("RVMR_NREPS", n_reps_default_v2))
  SEEDS_env <- Sys.getenv("RVMR_SEEDS", "")
  SEEDS <- if (nzchar(SEEDS_env)) as.integer(strsplit(SEEDS_env, ",")[[1]]) else master_seeds_v2
  cat(sprintf("Running run_all_v2() with n_reps=%d, seeds=%s ...\n",
              N, paste(SEEDS, collapse = ",")))
  R <- run_all_v2(n_reps = N, master_seeds = SEEDS)
  saveRDS(R, file = "results_v2.rds")
  write_results_v2_md(R, "results_v2.md")
  cat("\nWrote results_v2.rds and results_v2.md.\n")
}
