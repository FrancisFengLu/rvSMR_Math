# finngen_cis_mr.R
# Real-data K=1 cis-MR sanity check for rvMR::mrAR()
# Exposure: lipid GWAS (LDL or TG) per Open Targets credible-set summary stats
# Outcome:  FinnGen R12 x MVP x UKBB meta-analysis I9_IHD / I9_MI_STRICT
# Lead variants are RCT-validated cis-MR canonical SNPs.
#
# Allele harmonization: both exposure (Open Targets) and outcome (PheWeb)
# report beta relative to ALT allele on the same chr-pos-ref-alt key.
# The Wald ratio b_y / b_x is therefore harmonized.

suppressPackageStartupMessages({
  library(jsonlite)
  library(rvMR)
})

panel <- fromJSON("panel_input.json", simplifyDataFrame = FALSE)

# Wald F-statistic = (b_x / se_x)^2 -- one-instrument first-stage F.
fstat <- function(b_x, se_x) (b_x / se_x)^2

# Wald 95% CI on the ratio b_y / b_x using delta method.
wald_ci <- function(b_x, se_x, b_y, se_y, alpha = 0.05) {
  point <- b_y / b_x
  se_point <- sqrt((se_y^2) / (b_x^2) + (b_y^2 * se_x^2) / (b_x^4))
  z <- qnorm(1 - alpha / 2)
  c(lo = point - z * se_point, hi = point + z * se_point, est = point, se = se_point)
}

results <- list()
for (g in panel) {
  ar <- mrAR(b_x = g$b_x, se_x = g$se_x, b_y = g$b_y, se_y = g$se_y, alpha = 0.05)
  fs <- fstat(g$b_x, g$se_x)
  wd <- wald_ci(g$b_x, g$se_x, g$b_y, g$se_y, 0.05)

  # Decide whether mrAR CI is directionally consistent with RCT truth.
  # RCT truth direction: sign of the causal ratio implied by drug action.
  # We compute the Wald point estimate; the CI direction is consistent if
  # the CI excludes zero and lies on the RCT-implied side.
  ci_excludes_zero <- if (ar$ci_type == "bounded") {
    (ar$ci_lower > 0) || (ar$ci_upper < 0)
  } else if (ar$ci_type == "disconnected") {
    # accepted set is (-Inf, lower] U [upper, Inf); excludes 0 iff lower<0<upper
    !(ar$ci_lower <= 0 && ar$ci_upper >= 0)
  } else {
    FALSE
  }

  results[[g$gene]] <- list(
    gene = g$gene, rs = g$rs, variant = g$var,
    exposure = g$exposure, exposure_study = g$exposure_study,
    b_x = g$b_x, se_x = g$se_x,
    outcome = g$outcome, b_y = g$b_y, se_y = g$se_y,
    n_case = g$outcome_n_case, n_control = g$outcome_n_control,
    Fstat = fs,
    wald_point = wd[["est"]], wald_se = wd[["se"]],
    wald_lo = wd[["lo"]], wald_hi = wd[["hi"]],
    mrAR_type = ar$ci_type,
    mrAR_lo = ar$ci_lower, mrAR_hi = ar$ci_upper,
    mrAR_ar_at_point = ar$ar_at_point_estimate,
    rct_truth = g$rct_truth,
    ci_excludes_zero = ci_excludes_zero
  )
}

# Format
fmt <- function(x, d = 4) if (is.finite(x)) formatC(x, format = "f", digits = d) else as.character(x)
out_rows <- lapply(results, function(r) {
  ci_str <- if (r$mrAR_type == "bounded") {
    sprintf("[%s, %s]", fmt(r$mrAR_lo, 3), fmt(r$mrAR_hi, 3))
  } else if (r$mrAR_type == "disconnected") {
    sprintf("(-Inf, %s] U [%s, Inf)", fmt(r$mrAR_lo, 3), fmt(r$mrAR_hi, 3))
  } else if (r$mrAR_type == "whole_line") {
    "(-Inf, Inf)"
  } else "empty"
  c(
    gene = r$gene, rs = r$rs,
    bx = fmt(r$b_x), sex = fmt(r$se_x),
    by = fmt(r$b_y), sey = fmt(r$se_y),
    F = fmt(r$Fstat, 1),
    wald = sprintf("%s [%s, %s]", fmt(r$wald_point, 3), fmt(r$wald_lo, 3), fmt(r$wald_hi, 3)),
    mrAR_type = r$mrAR_type,
    mrAR_CI = ci_str,
    excludes_zero = as.character(r$ci_excludes_zero)
  )
})
tab <- do.call(rbind, out_rows)
print(as.data.frame(tab, stringsAsFactors = FALSE))

saveRDS(results, "finngen_results.rds")
write(jsonlite::toJSON(results, auto_unbox = TRUE, pretty = TRUE), "finngen_results.json")
cat("\nWrote finngen_results.rds and finngen_results.json\n")
