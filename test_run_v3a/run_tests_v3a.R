# =====================================================================
# run_tests_v3a.R  (Round 3 / Worker A)
#
# Changes from Round 2 driver (CRITIQUE_v2 follow-ups):
#
#   S1.v2.1 (run completion):     this driver writes results_v3a.rds
#                                  and results_v3a.md at end with an
#                                  explicit "Wrote ..." line.
#
#   S1.v2.2 / S2.v2.3 (TSLS name): tsls_summary() renamed to
#                                  ivw_summary(). All `tsls_*` result
#                                  columns are `ivw_summary_*`. Prose
#                                  does NOT claim Wang-Kang Fig 6
#                                  reproduction; the sign-alternated
#                                  cells (separate scenarios) are
#                                  where we make the Wang-Kang
#                                  comparison.
#
#   S2.v2.1 (seed aliasing):      per-rep seed is now a hash of
#                                  (scenario name, master seed, rep)
#                                  via digest::digest2int. This
#                                  eliminates the modular-linear seed
#                                  collision that aliased Plei_mult0
#                                  with Fsweep_F20 in Round 2.
#
#   S1.6 (bias|bounded):          preserved from Round 2.
#
# Comparators reported (per rep):
#   AR             : mrAR_multi (the method under test).
#   IVW-of-ratios  : per-mask Wald ratio + delta-method SE, IVW mean,
#                     Wald CI. CRITIQUE S1.1 canonical-form.
#   ivw_summary    : (b_x' W b_x)^{-1} b_x' W b_y with W = diag(1/SE_y^2),
#                     SE = (b_x' W b_x)^{-1/2}. Equivalent to
#                     MendelianRandomization::mr_ivw(method="default")
#                     (Burgess et al.). This is the scalar IVW form,
#                     NOT 2SLS. It is the canonical comparator that
#                     under-covers at weak IV under sign-alternated
#                     alpha (Wang-Kang 2022 §3) and also collapses
#                     in Patel-Lane-Burgess 2024 Fig 2.
# =====================================================================

suppressMessages(suppressWarnings(devtools::load_all(
  "/home/francisfenglu4/rvSMR/May_30md/rvMR",
  quiet = TRUE)))

suppressMessages(suppressWarnings(library(digest)))

source("generate_test_data_v3a.R")
source("scenarios_v3a.R")

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

# Classical IVW-of-ratios (CRITIQUE S1.1). Delta-method SE makes this
# partially weak-IV-buffered at low F (over-coverage), as Round 2
# empirically confirmed.
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

# Scalar IVW (a.k.a. summary IVW), Burgess MendelianRandomization::mr_ivw
# form. NOT 2SLS in the strict individual-level sense; the v2 driver
# called this "tsls_summary" which was non-standard (CRITIQUE_v2
# §S1.v2.2 / §S2.v2.3). We rename to ivw_summary throughout.
ivw_summary <- function(d) {
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

# Hash-based per-rep seed (CRITIQUE_v2 §S2.v2.1 fix).
# digest2int returns a signed 32-bit int in (-2^31, 2^31-1). We map it
# to a positive R-compatible integer for set.seed().
hash_seed <- function(scenario_name, master_seed, r) {
  key <- paste(scenario_name, as.integer(master_seed), as.integer(r),
               sep = "|")
  h <- digest::digest2int(key)
  # Map to non-negative range; set.seed accepts any integer but using
  # positive ints keeps log output readable.
  s <- as.integer(h %% .Machine$integer.max)
  if (s < 0L) s <- s + .Machine$integer.max
  if (s == 0L) s <- 1L
  s
}

# ---- main per-scenario per-seed runner ------------------------------

run_one_scenario_seed <- function(scn_name, scn, n_reps, master_seed) {
  res <- data.frame(
    rep            = integer(n_reps),
    seed_used      = integer(n_reps),
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
    ivw_summary_beta      = numeric(n_reps),
    ivw_summary_se        = numeric(n_reps),
    ivw_summary_lo        = numeric(n_reps),
    ivw_summary_hi        = numeric(n_reps),
    ivw_summary_covers_beta = logical(n_reps),
    ivw_summary_covers_zero = logical(n_reps),
    ivw_summary_len       = numeric(n_reps),
    stringsAsFactors = FALSE
  )

  beta_true <- scn$beta_true
  K         <- scn$K

  se_y_proxy <- 1 / sqrt(scn$n_y)
  pleio_size <- if (!is.null(scn$pleio_size_mult)) {
    scn$pleio_size_mult * se_y_proxy
  } else {
    NULL
  }

  for (r in seq_len(n_reps)) {
    seed_r <- hash_seed(scn_name, master_seed, r)

    d <- tryCatch(
      simulate_burden_mr_v3a(
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
      ),
      error = function(e) NULL
    )

    if (is.null(d)) {
      res$rep[r]     <- r
      res$seed_used[r] <- seed_r
      res$ci_type[r] <- "dgp_error"
      next
    }

    # Inference covariance choice
    if (scn$dgp == "overlap") {
      R_xx_inf <- d$R_xx; R_yy_inf <- d$R_yy; R_xy_inf <- d$R_xy
    } else if (scn$dgp == "ld_xx") {
      R_xx_inf <- d$R_xx; R_yy_inf <- diag(K); R_xy_inf <- matrix(0, K, K)
    } else if (scn$dgp == "confounder") {
      R_xx_inf <- diag(K); R_yy_inf <- diag(K); R_xy_inf <- matrix(0, K, K)
    } else {
      # "ar" or "signalt": default diagonals; signalt does NOT change
      # the residual covariance structure passed to AR, only alpha_k
      # signs (which AR is moment-invariant to).
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
      res$seed_used[r] <- seed_r
      res$ci_type[r] <- "error"
      next
    }

    res$rep[r]            <- r
    res$seed_used[r]      <- seed_r
    res$ci_type[r]        <- ar$ci_type
    res$ci_covers_beta[r] <- in_AR_CI(beta_true, ar$ci_type, ar$ci_intervals)
    res$ci_covers_zero[r] <- in_AR_CI(0,         ar$ci_type, ar$ci_intervals)
    res$ci_len_finite[r]  <- ci_length(ar$ci_type, ar$ci_intervals)
    res$beta_hat[r]       <- ar$beta_hat
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

    iv <- ivw_summary(d)
    res$ivw_summary_beta[r]          <- iv["beta"]
    res$ivw_summary_se[r]            <- iv["se"]
    res$ivw_summary_lo[r]            <- iv["lo"]
    res$ivw_summary_hi[r]            <- iv["hi"]
    res$ivw_summary_covers_beta[r]   <- isTRUE(iv["lo"] <= beta_true && beta_true <= iv["hi"])
    res$ivw_summary_covers_zero[r]   <- isTRUE(iv["lo"] <= 0         && 0         <= iv["hi"])
    res$ivw_summary_len[r]           <- iv["hi"] - iv["lo"]
  }

  res
}

# ---- per-scenario aggregator ---------------------------------------

summarize_one_seed <- function(scn, df) {
  ok <- !(df$ci_type %in% c("error", "dgp_error", ""))
  n  <- sum(ok)
  cov_ar      <- if (n > 0) mean(df$ci_covers_beta[ok]) else NA_real_
  cov_naive   <- if (n > 0) mean(df$naive_covers_beta[ok], na.rm = TRUE) else NA_real_
  cov_ivw_summary <- if (n > 0) mean(df$ivw_summary_covers_beta[ok], na.rm = TRUE) else NA_real_
  rej_zero    <- if (n > 0) mean(!df$ci_covers_zero[ok]) else NA_real_
  rej_zero_nv <- if (n > 0) mean(!df$naive_covers_zero[ok], na.rm = TRUE) else NA_real_
  rej_zero_iv <- if (n > 0) mean(!df$ivw_summary_covers_zero[ok], na.rm = TRUE) else NA_real_
  F_mean      <- if (n > 0) mean(df$F_mean[ok]) else NA_real_
  bounded_mask <- ok & (df$ci_type == "bounded_interval")
  n_bounded   <- sum(bounded_mask)
  bias_mean   <- if (n_bounded > 0) mean(df$bias[bounded_mask]) else NA_real_
  bias_median <- if (n_bounded > 0) median(df$bias[bounded_mask]) else NA_real_
  J_rej_rate  <- if (n > 0) mean(df$J_pvalue[ok] < 0.05, na.rm = TRUE) else NA_real_
  shape_tab   <- if (n > 0) prop.table(table(df$ci_type[ok])) else table(integer(0))
  n_dgp_err   <- sum(df$ci_type == "dgp_error")
  n_ar_err    <- sum(df$ci_type == "error")

  list(
    label             = scn$label,
    n_ok              = n,
    n_bounded         = n_bounded,
    n_dgp_err         = n_dgp_err,
    n_ar_err          = n_ar_err,
    coverage_AR       = cov_ar,
    coverage_naive    = cov_naive,
    coverage_ivw_summary = cov_ivw_summary,
    rej_zero_AR       = rej_zero,
    rej_zero_naive    = rej_zero_nv,
    rej_zero_ivw_summary = rej_zero_iv,
    F_mean            = F_mean,
    bias_mean         = bias_mean,
    bias_median       = bias_median,
    J_rej_rate        = J_rej_rate,
    shape_distribution = shape_tab
  )
}

aggregate_seeds <- function(per_seed) {
  cov_ar   <- sapply(per_seed, function(s) s$coverage_AR)
  cov_nv   <- sapply(per_seed, function(s) s$coverage_naive)
  cov_iv   <- sapply(per_seed, function(s) s$coverage_ivw_summary)
  rej_zero <- sapply(per_seed, function(s) s$rej_zero_AR)
  rej_zero_nv <- sapply(per_seed, function(s) s$rej_zero_naive)
  rej_zero_iv <- sapply(per_seed, function(s) s$rej_zero_ivw_summary)
  Fm       <- sapply(per_seed, function(s) s$F_mean)
  J_rej    <- sapply(per_seed, function(s) s$J_rej_rate)
  n_bd     <- sapply(per_seed, function(s) s$n_bounded)
  bias_bd  <- sapply(per_seed, function(s) s$bias_mean)
  n_dgperr <- sapply(per_seed, function(s) s$n_dgp_err)
  n_arerr  <- sapply(per_seed, function(s) s$n_ar_err)
  list(
    coverage_AR_mean     = mean(cov_ar, na.rm = TRUE),
    coverage_AR_se       = sd(cov_ar, na.rm = TRUE) / sqrt(sum(!is.na(cov_ar))),
    coverage_AR_per_seed = cov_ar,
    coverage_naive_mean  = mean(cov_nv, na.rm = TRUE),
    coverage_naive_se    = sd(cov_nv, na.rm = TRUE) / sqrt(sum(!is.na(cov_nv))),
    coverage_naive_per_seed = cov_nv,
    coverage_ivw_summary_mean   = mean(cov_iv, na.rm = TRUE),
    coverage_ivw_summary_se     = sd(cov_iv, na.rm = TRUE) / sqrt(sum(!is.na(cov_iv))),
    coverage_ivw_summary_per_seed = cov_iv,
    rej_zero_AR_mean     = mean(rej_zero, na.rm = TRUE),
    rej_zero_naive_mean  = mean(rej_zero_nv, na.rm = TRUE),
    rej_zero_ivw_summary_mean   = mean(rej_zero_iv, na.rm = TRUE),
    F_mean               = mean(Fm, na.rm = TRUE),
    J_rej_rate_mean      = mean(J_rej, na.rm = TRUE),
    n_bounded_mean       = mean(n_bd, na.rm = TRUE),
    bias_mean_bounded    = mean(bias_bd, na.rm = TRUE),
    n_dgp_err_total      = sum(n_dgperr),
    n_ar_err_total       = sum(n_arerr)
  )
}

# ---- main loop ------------------------------------------------------

run_all_v3a <- function(n_reps = n_reps_default_v3a,
                        master_seeds = master_seeds_v3a,
                        scenarios = scenarios_v3a,
                        verbose = TRUE) {
  results   <- list()
  summaries <- list()
  per_seed  <- list()
  t0 <- Sys.time()
  for (i in seq_along(scenarios)) {
    nm  <- names(scenarios)[i]
    scn <- scenarios[[nm]]
    if (verbose) {
      cat(sprintf("\n=== %s ===\n", scn$label))
    }
    results[[nm]] <- list()
    per_seed[[nm]] <- list()
    for (s_idx in seq_along(master_seeds)) {
      ms <- master_seeds[s_idx]
      if (verbose) cat(sprintf("  seed %d (idx %d) ... ", ms, s_idx))
      tic <- Sys.time()
      df <- run_one_scenario_seed(nm, scn, n_reps = n_reps,
                                  master_seed = ms)
      sm <- summarize_one_seed(scn, df)
      results[[nm]][[s_idx]]  <- df
      per_seed[[nm]][[s_idx]] <- sm
      toc <- Sys.time()
      if (verbose) {
        cat(sprintf(
          "done (%.0fs): cov_AR=%.3f cov_IVWratios=%.3f cov_IVWsummary=%.3f F_mean=%.2f n_bd=%d n_dgperr=%d\n",
          as.numeric(toc - tic, units = "secs"),
          ifelse(is.na(sm$coverage_AR), NaN, sm$coverage_AR),
          ifelse(is.na(sm$coverage_naive), NaN, sm$coverage_naive),
          ifelse(is.na(sm$coverage_ivw_summary), NaN, sm$coverage_ivw_summary),
          ifelse(is.na(sm$F_mean), NaN, sm$F_mean),
          sm$n_bounded, sm$n_dgp_err))
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
       scenarios = scenarios,
       duration_sec = as.numeric(t1 - t0, units = "secs"))
}

# ---- markdown writer (results_v3a.md) ------------------------------

write_results_v3a_md <- function(R, path) {
  con <- file(path, open = "w")
  on.exit(close(con))
  cat("# rvMR Round-3 / Worker A simulation results\n\n", file = con)
  cat(sprintf("- Replicates per scenario per seed: **%d**\n", R$n_reps), file = con)
  cat(sprintf("- Master seeds: %s\n",
              paste(R$master_seeds, collapse = ", ")), file = con)
  cat(sprintf("- Total wall time: %.1fs\n", R$duration_sec), file = con)
  cat(sprintf("- R version: %s\n", R.version.string), file = con)
  cat(sprintf("- rvMR loaded from /home/francisfenglu4/rvSMR/May_30md/rvMR (Worker B territory; not modified)\n\n"), file = con)

  cat("Comparators reported:\n", file = con)
  cat("  - **IVW-of-ratios** (delta-method SE): per-mask Wald ratio + IVW mean. Partially weak-IV-buffered (Round 2 confirmed over-cover at low lambda).\n", file = con)
  cat("  - **scalar IVW (`ivw_summary`)**: (b_x' W b_x)^{-1} b_x' W b_y with W = diag(1/SE_y^2). Equivalent to `MendelianRandomization::mr_ivw(method=\"default\")`. NOT 2SLS (CRITIQUE_v2 §S2.v2.3 rename from v2 `tsls_*`).\n\n", file = con)

  # --- F-sweep (homogeneous-sign) table ---
  cat("## (a) Homogeneous-sign weak-IV F-sweep (the cells missing from Round 2)\n\n", file = con)
  cat("alpha_k = +sqrt(F)*SE_x for all k. K=3, lambda_joint = K*F. This is the row that\n", file = con)
  cat("directly extends the Round-2 v2 headline table at F in {0.25, 0.5, 1}.\n\n", file = con)
  cat("| F | lambda | AR cov (mean +/- SE) | IVW-ratios cov | scalar-IVW cov | F_mean | n_bounded | bias|bounded |\n", file = con)
  cat("|---:|---:|---|---|---|---:|---:|---|\n", file = con)
  for (nm in names(R$summaries)) {
    if (!grepl("^Fsweep_", nm)) next
    s <- R$summaries[[nm]]
    scn <- R$scenarios[[nm]]
    cat(sprintf("| %g | %.2f | %.3f +/- %.3f | %.3f +/- %.3f | %.3f +/- %.3f | %.2f | %.0f | %s |\n",
                scn$F_target, 3 * scn$F_target,
                s$coverage_AR_mean, s$coverage_AR_se,
                s$coverage_naive_mean, s$coverage_naive_se,
                s$coverage_ivw_summary_mean, s$coverage_ivw_summary_se,
                s$F_mean, s$n_bounded_mean,
                if (is.finite(s$bias_mean_bounded))
                  sprintf("%+.4f", s$bias_mean_bounded) else "NA"),
        file = con)
  }
  cat("\n", file = con)

  # --- Sign-alternated F-sweep table ---
  cat("## (b) Sign-alternated alpha_k weak-IV sweep (Wang-Kang 2022 §3 style)\n\n", file = con)
  cat("alpha_k alternates sign across K=3 masks (+ - +) so the pooled b_x summary statistic\n", file = con)
  cat("is near zero in expectation. Scalar IVW (b_x' W b_x)^{-1} factor blows up at low F.\n", file = con)
  cat("AR is sign-invariant in the moment m_k = b_y_k - beta * b_x_k and should hold at nominal.\n\n", file = con)
  cat("Reference: Wang & Kang 2022 Biometrics 78(4):1699-1713, §3 (DOI 10.1111/biom.13524).\n\n", file = con)
  cat("| F | lambda | AR cov (mean +/- SE) | IVW-ratios cov | scalar-IVW cov | F_mean | n_bounded | bias|bounded |\n", file = con)
  cat("|---:|---:|---|---|---|---:|---:|---|\n", file = con)
  for (nm in names(R$summaries)) {
    if (!grepl("^SignAlt_", nm)) next
    s <- R$summaries[[nm]]
    scn <- R$scenarios[[nm]]
    cat(sprintf("| %g | %.2f | %.3f +/- %.3f | %.3f +/- %.3f | %.3f +/- %.3f | %.2f | %.0f | %s |\n",
                scn$F_target, 3 * scn$F_target,
                s$coverage_AR_mean, s$coverage_AR_se,
                s$coverage_naive_mean, s$coverage_naive_se,
                s$coverage_ivw_summary_mean, s$coverage_ivw_summary_se,
                s$F_mean, s$n_bounded_mean,
                if (is.finite(s$bias_mean_bounded))
                  sprintf("%+.4f", s$bias_mean_bounded) else "NA"),
        file = con)
  }
  cat("\n", file = con)

  # --- Confounder strength sweep ---
  cat("## (c) Confounder-strength sweep at F=20\n\n", file = con)
  cat("Non-AR DGP via shared latent u: b_x_k = alpha_k + cs*SE_x*u + s_idio*SE_x*eps_x,\n", file = con)
  cat("b_y_k = beta*alpha_k + cs*SE_y*u + s_idio*SE_y*eps_y; s_idio = sqrt(1-cs^2).\n", file = con)
  cat("Inference call uses default R_xx = R_yy = I, R_xy = 0 (canonical \"user does not know about u\" stress test).\n", file = con)
  cat("cs > 1 violates s_idio's domain and is reported as DGP-error (skipped).\n\n", file = con)
  cat("| cs | AR cov (mean +/- SE) | IVW-ratios cov | scalar-IVW cov | rej_zero (AR) | F_mean | n_bounded | n_dgp_err |\n", file = con)
  cat("|---:|---|---|---|---:|---:|---:|---:|\n", file = con)
  for (nm in names(R$summaries)) {
    if (!grepl("^ConfSweep_", nm)) next
    s <- R$summaries[[nm]]
    scn <- R$scenarios[[nm]]
    cat(sprintf("| %g | %.3f +/- %.3f | %.3f +/- %.3f | %.3f +/- %.3f | %.3f | %.2f | %.0f | %d |\n",
                scn$conf_strength,
                s$coverage_AR_mean, s$coverage_AR_se,
                s$coverage_naive_mean, s$coverage_naive_se,
                s$coverage_ivw_summary_mean, s$coverage_ivw_summary_se,
                s$rej_zero_AR_mean,
                s$F_mean, s$n_bounded_mean, s$n_dgp_err_total),
        file = con)
  }
  cat("\n", file = con)

  # --- Per-seed coverage breakdown ---
  cat("## Per-seed AR coverage (each row = scenario; columns = 3 seeds)\n\n", file = con)
  cat("Seeds: 20260603, 99887766, 1234567 (spread-out, CRITIQUE_v2 §S3.v2.3).\n\n", file = con)
  cat("| Scenario | seed1 (20260603) | seed2 (99887766) | seed3 (1234567) |\n|---|---:|---:|---:|\n", file = con)
  for (nm in names(R$summaries)) {
    s <- R$summaries[[nm]]
    cps <- s$coverage_AR_per_seed
    cat(sprintf("| %s | %.3f | %.3f | %.3f |\n",
                R$scenarios[[nm]]$label,
                cps[1], cps[2], cps[3]),
        file = con)
  }
  cat("\n", file = con)
}

# ---- seed-aliasing verification -------------------------------------
# Required by Round 3 brief: confirm two adjacent scenarios at seed 1
# are NOT bit-identical with the hash-based per-rep seed.
verify_seed_decoupling <- function(scn_a = "Fsweep_F1",
                                   scn_b = "Fsweep_F0.5",
                                   master_seed = master_seeds_v3a[1],
                                   r_grid = 1:5) {
  cat("\n--- seed decoupling check ---\n")
  cat(sprintf("scenario A: %s, scenario B: %s, master_seed: %d\n",
              scn_a, scn_b, master_seed))
  for (r in r_grid) {
    sa <- hash_seed(scn_a, master_seed, r)
    sb <- hash_seed(scn_b, master_seed, r)
    cat(sprintf("  rep %d: seed_A=%d  seed_B=%d  identical=%s\n",
                r, sa, sb, identical(sa, sb)))
  }
}

# ---- entry point ----------------------------------------------------

if (!interactive() && identical(Sys.getenv("RVMR_RUN_MAIN"), "1")) {
  N <- as.integer(Sys.getenv("RVMR_NREPS", n_reps_default_v3a))
  SEEDS_env <- Sys.getenv("RVMR_SEEDS", "")
  SEEDS <- if (nzchar(SEEDS_env)) as.integer(strsplit(SEEDS_env, ",")[[1]]) else master_seeds_v3a
  SCN_FILTER <- Sys.getenv("RVMR_SCN_FILTER", "")
  scns <- scenarios_v3a
  if (nzchar(SCN_FILTER)) {
    keep <- grep(SCN_FILTER, names(scns), value = TRUE)
    scns <- scns[keep]
    cat(sprintf("Scenario filter '%s' matches %d cells: %s\n",
                SCN_FILTER, length(scns), paste(keep, collapse = ", ")))
  }
  cat(sprintf("Running run_all_v3a() with n_reps=%d, seeds=%s ...\n",
              N, paste(SEEDS, collapse = ",")))
  verify_seed_decoupling()
  R <- run_all_v3a(n_reps = N, master_seeds = SEEDS, scenarios = scns)
  saveRDS(R, file = "results_v3a.rds")
  write_results_v3a_md(R, "results_v3a.md")
  cat("\nWrote results_v3a.rds and results_v3a.md.\n")
}
