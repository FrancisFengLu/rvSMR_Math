# =====================================================================
# bmi_sbp_sanity.R  (Round 3, Worker B)
#
# Real-data sanity check: run rvMR::mrAR (K=1, single-IV Anderson-Rubin)
# per-SNP on the BMI -> SBP dataset shipped with the mr.raps CRAN
# package, then meta-analyze. Compare to Wang & Kang 2022 Table 1
# (Biometrics 78(4): 1699-1713), which is itself a replication of
# Zhao et al. 2020 Annals of Stats.
#
# Wang-Kang Table 1 reports BMI -> SBP causal effect (mmHg per s.d. BMI,
# both standardized) under several estimators, using two IV sets:
#   - p-threshold = 5e-8  (genome-wide significant): nominally 26-28 IVs,
#     mr.raps bmi.sbp has 25 SNPs at < 5e-8 after the package's QC.
#   - p-threshold = 1e-4  (suggestive, the weak-IV regime that motivates
#     AR): 160 SNPs in mr.raps bmi.sbp.
#
# Wang-Kang Table 1 BMI->SBP point/CI (from their published numbers
# for the BMI-SBP application; mr.raps and various AR variants in their
# Tables; Zhao 2020 Table 5 originally reported ~0.31-0.40 mmHg/s.d.
# range across IVW/MR-Egger/RAPS; AR-type intervals on the 25-SNP set
# nest those Wald-style point estimates).
#
# We exercise rvMR::mrAR per SNP (the K=1 closed-form inverted Fieller
# 1954, our single-IV AR), and aggregate using inverse-variance
# meta-analysis of per-SNP Wald ratios for the point estimate, with
# an AR-style intersection of acceptance sets across SNPs as the
# WK-comparable CI (a Q-style version).
# =====================================================================

suppressPackageStartupMessages({
  .libPaths(c("/home/francisfenglu4/R/library", .libPaths()))
  library(mr.raps)
  library(devtools)
  load_all("/home/francisfenglu4/rvSMR/May_30md/rvMR")
})

# --- Load data ---------------------------------------------------------
data("bmi.sbp", package = "mr.raps")
df_all <- bmi.sbp

cat("Loaded mr.raps::bmi.sbp:\n")
cat("  rows (SNPs):", nrow(df_all), "\n")
cat("  cols:", ncol(df_all), "\n")
cat("  N at p_selection < 5e-8:", sum(df_all$pval.selection < 5e-8), "\n")
cat("  N at p_selection < 1e-4:", sum(df_all$pval.selection < 1e-4), "\n")

# Subset to the two Wang-Kang IV sets.
df_25  <- df_all[df_all$pval.selection < 5e-8, ]
df_160 <- df_all[df_all$pval.selection < 1e-4, ]
stopifnot(nrow(df_25) == 25L, nrow(df_160) == 160L)

# --- Per-SNP rvMR::mrAR (K=1) -----------------------------------------
run_per_snp_mrAR <- function(df, label) {
  cat("\n=== Per-SNP mrAR on", label, "(", nrow(df), "SNPs) ===\n")
  out <- vector("list", nrow(df))
  for (i in seq_len(nrow(df))) {
    bx  <- df$beta.exposure[i]
    sx  <- df$se.exposure[i]
    by_ <- df$beta.outcome[i]
    sy  <- df$se.outcome[i]
    ar  <- mrAR(b_x = bx, se_x = sx, b_y = by_, se_y = sy, alpha = 0.05)
    wald <- by_ / bx                       # per-SNP Wald ratio (point est)
    se_wald <- sqrt((sy / bx)^2 + (by_^2 * sx^2) / (bx^4))  # delta method
    F_stat <- (bx / sx)^2
    out[[i]] <- data.frame(
      SNP            = df$SNP[i],
      b_x            = bx, se_x = sx,
      b_y            = by_, se_y = sy,
      F_stat         = F_stat,
      wald_ratio     = wald,
      se_wald        = se_wald,
      ar_ci_type     = ar$ci_type,
      ar_ci_lower    = ar$ci_lower,
      ar_ci_upper    = ar$ci_upper,
      pval_selection = df$pval.selection[i],
      stringsAsFactors = FALSE
    )
  }
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  cat("Per-SNP F summary:\n");        print(summary(res$F_stat))
  cat("Per-SNP Wald summary:\n");     print(summary(res$wald_ratio))
  cat("CI-shape table:\n");           print(table(res$ar_ci_type))
  res
}

per25  <- run_per_snp_mrAR(df_25,  "25-SNP set (p<5e-8)")
per160 <- run_per_snp_mrAR(df_160, "160-SNP set (p<1e-4)")

# --- Inverse-variance meta-analysis of Wald ratios --------------------
ivw_meta <- function(per, label) {
  # Standard IVW: beta_hat = sum(w_k beta_k) / sum(w_k), w_k = 1/se^2.
  # Q-statistic = sum(w_k (beta_k - beta_hat)^2) ~ chi^2_{K-1} under
  # homogeneity (Cochran). Random-effects CI uses sqrt(1/sum(w_k))
  # with optional tau^2 inflation; here we report fixed-effect SE.
  w <- 1 / per$se_wald^2
  bh <- sum(w * per$wald_ratio) / sum(w)
  se <- sqrt(1 / sum(w))
  ci_lo <- bh - 1.96 * se
  ci_hi <- bh + 1.96 * se
  Q  <- sum(w * (per$wald_ratio - bh)^2)
  df <- nrow(per) - 1L
  pQ <- pchisq(Q, df = df, lower.tail = FALSE)
  list(label = label, beta_hat = bh, se = se,
       ci_lo = ci_lo, ci_hi = ci_hi,
       Q = Q, Q_df = df, Q_pvalue = pQ,
       K = nrow(per),
       mean_F = mean(per$F_stat))
}

meta25  <- ivw_meta(per25,  "25-SNP IVW meta")
meta160 <- ivw_meta(per160, "160-SNP IVW meta")

cat("\n=== Meta-analytic summaries ===\n")
for (m in list(meta25, meta160)) {
  cat(sprintf("[%s] K=%d  beta=%.4f  SE=%.4f  95%%CI=[%.4f, %.4f]  Q=%.2f df=%d p=%.3g  mean(F)=%.1f\n",
              m$label, m$K, m$beta_hat, m$se,
              m$ci_lo, m$ci_hi, m$Q, m$Q_df, m$Q_pvalue, m$mean_F))
}

# --- AR-style meta CI by intersecting per-SNP acceptance sets ---------
# For each SNP, AR accepts a set S_k. The K-AR-style joint acceptance
# set is the INTERSECTION (under the homogeneity assumption). For
# bounded per-SNP intervals, the intersection is [max(lo), min(hi)]
# (empty if max(lo) > min(hi)).
ar_intersection <- function(per) {
  # Drop SNPs with disconnected / whole_line / empty (handle separately).
  bnd <- per[per$ar_ci_type == "bounded", ]
  cnt <- table(per$ar_ci_type)
  if (nrow(bnd) == 0L) {
    return(list(lo = NA, hi = NA, n_bounded = 0L, ci_shape_table = cnt))
  }
  lo <- max(bnd$ar_ci_lower, na.rm = TRUE)
  hi <- min(bnd$ar_ci_upper, na.rm = TRUE)
  list(lo = lo, hi = hi, n_bounded = nrow(bnd),
       n_disc = sum(per$ar_ci_type == "disconnected"),
       n_whole = sum(per$ar_ci_type == "whole_line"),
       n_empty = sum(per$ar_ci_type == "empty"),
       ci_shape_table = cnt,
       feasible = is.finite(lo) && is.finite(hi) && lo <= hi)
}

ari25  <- ar_intersection(per25)
ari160 <- ar_intersection(per160)

cat("\n=== AR-intersection CIs (each SNP's bounded CI, intersected) ===\n")
cat(sprintf("[25-SNP]  bounded=%d disconnected=%d whole_line=%d empty=%d\n",
            ari25$n_bounded, ari25$n_disc, ari25$n_whole, ari25$n_empty))
cat(sprintf("          intersect = [%.4f, %.4f]  feasible=%s\n",
            ari25$lo, ari25$hi, ari25$feasible))
cat(sprintf("[160-SNP] bounded=%d disconnected=%d whole_line=%d empty=%d\n",
            ari160$n_bounded, ari160$n_disc, ari160$n_whole, ari160$n_empty))
cat(sprintf("          intersect = [%.4f, %.4f]  feasible=%s\n",
            ari160$lo, ari160$hi, ari160$feasible))

# --- mr.raps reference numbers (the package's own MR estimator) -------
# Run mr.raps's own estimator for direct cross-check (same dataset,
# different point-estimate method). This is the closest thing to a
# "ground truth" within the same data source.
cat("\n=== mr.raps native estimator (reference) ===\n")
raps_25  <- mr.raps(b_exp = df_25$beta.exposure,
                    b_out = df_25$beta.outcome,
                    se_exp = df_25$se.exposure,
                    se_out = df_25$se.outcome,
                    over.dispersion = TRUE, loss.function = "huber")
cat(sprintf("[25-SNP]  RAPS beta=%.4f  SE=%.4f  95%%CI=[%.4f, %.4f]\n",
            raps_25$beta.hat, raps_25$beta.se,
            raps_25$beta.hat - 1.96*raps_25$beta.se,
            raps_25$beta.hat + 1.96*raps_25$beta.se))
raps_160 <- mr.raps(b_exp = df_160$beta.exposure,
                    b_out = df_160$beta.outcome,
                    se_exp = df_160$se.exposure,
                    se_out = df_160$se.outcome,
                    over.dispersion = TRUE, loss.function = "huber")
cat(sprintf("[160-SNP] RAPS beta=%.4f  SE=%.4f  95%%CI=[%.4f, %.4f]\n",
            raps_160$beta.hat, raps_160$beta.se,
            raps_160$beta.hat - 1.96*raps_160$beta.se,
            raps_160$beta.hat + 1.96*raps_160$beta.se))

# --- Save artifacts ----------------------------------------------------
saveRDS(list(per25 = per25, per160 = per160,
             meta25 = meta25, meta160 = meta160,
             ari25 = ari25, ari160 = ari160,
             raps_25 = raps_25, raps_160 = raps_160),
        "/home/francisfenglu4/projects/rvSMR_Math/test_run_v3b/bmi_sbp_results.rds")

write.csv(per25,  "/home/francisfenglu4/projects/rvSMR_Math/test_run_v3b/per_snp_25.csv",
          row.names = FALSE)
write.csv(per160, "/home/francisfenglu4/projects/rvSMR_Math/test_run_v3b/per_snp_160.csv",
          row.names = FALSE)

cat("\nDONE. Results saved to test_run_v3b/.\n")
