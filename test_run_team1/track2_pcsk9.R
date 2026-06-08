# Track 2: PCSK9 -> CHD via common-variant cis-MR using mrAR_multi().
#
# rvSMR VALIDATION_PLAN.md §3 Track 2 spec: "real-data-sized end-to-end run"
# proving the pipeline handles real LD/annotation/I/O. Direction sanity:
# PCSK9 expression UP -> LDL UP -> CHD UP, so b_xy > 0 expected (Wald ratio
# in units of log-OR-CHD per SD-PCSK9-expression-residual).
#
# Substrate (substitutions documented in the report MD):
#   exposure b_x: GTEx v8 LIVER cis-eQTLs (eQTL Catalogue QTD000266, n=208)
#                 -- substitution for TenK10K Phase 1 because (a) PCSK9 is
#                 hepatocyte-expressed, not PBMC-expressed; (b) TenK10K
#                 rare-variant Zenodo zips are 214-260 byte placeholders.
#   outcome b_y:  FinnGen R12 x MVP x UKBB joint meta, endpoint I9_IHD
#                 (Worker B's validated route at mvp-ukbb.finngen.fi).
#
# Instruments: K = 7 lead variants (200 kb position-window LD prune) in the
# PCSK9 +- 1 Mb cis window with non-NA cis-eQTL and non-NA I9_IHD meta beta.
#
# Output:
#   * track2_results.json: structured output (CI, F per-IV, J-stat, J-p,
#     per-cell-type breakdown placeholder)
#   * console: pretty-printed table.

suppressPackageStartupMessages({
  library(jsonlite)
  library(rvMR)
})

set.seed(2026)

panel <- fromJSON("pcsk9_track2_panel.json", simplifyDataFrame = FALSE)
inst  <- panel$instruments
K     <- length(inst)
cat(sprintf("Track 2: PCSK9 -> CHD via mrAR_multi(K = %d)\n", K))
cat(sprintf("  exposure substrate: %s\n", panel$exposure_substrate))
cat(sprintf("  outcome substrate:  %s\n", panel$outcome_substrate))

# ---- assemble vectors --------------------------------------------------------
b_x  <- vapply(inst, function(r) as.numeric(r$eqtl_beta), numeric(1))
se_x <- vapply(inst, function(r) as.numeric(r$eqtl_se),   numeric(1))
b_y  <- vapply(inst, function(r) as.numeric(r$y_beta),    numeric(1))
se_y <- vapply(inst, function(r) as.numeric(r$y_se),      numeric(1))
rsid <- vapply(inst, function(r) as.character(r$rsid),    character(1))
pos  <- vapply(inst, function(r) as.integer(r$pos),       integer(1))

cat("\n--- per-IV instrument table ---\n")
tab <- data.frame(
  k        = seq_len(K),
  rsid     = rsid,
  pos      = pos,
  b_x      = round(b_x, 4),  se_x = round(se_x, 4),
  b_y      = round(b_y, 4),  se_y = round(se_y, 4),
  F        = round((b_x / se_x)^2, 2),
  wald_per = round(b_y / b_x, 3),
  stringsAsFactors = FALSE
)
print(tab, row.names = FALSE)

# ---- run mrAR_multi ----------------------------------------------------------
# R_xy is the across-sample correlation block. GTEx Liver and FinnGen meta
# are non-overlapping samples (US/Brazil GTEx vs Finland/MVP-US/UK Biobank
# meta); R_xy = 0 is appropriate.
cat("\n--- mrAR_multi (alpha = 0.05) ---\n")
res <- mrAR_multi(
  b_x = b_x, se_x = se_x,
  b_y = b_y, se_y = se_y,
  alpha = 0.05
)
cat("  ci_type:    ", res$ci_type, "\n")
cat("  ci_intervals:\n")
for (ix in seq_along(res$ci_intervals)) {
  iv <- res$ci_intervals[[ix]]
  cat(sprintf("    [%g, %g]\n", iv[1], iv[2]))
}
cat(sprintf("  beta_hat:   %.4f\n", res$beta_hat))
cat(sprintf("  J_stat:     %.4f  (df = K - 1 = %d)\n", res$J_stat, K - 1L))
cat(sprintf("  J_pvalue:   %.4g\n", res$J_pvalue))
cat(sprintf("  ar_crit:    %.4f\n", res$ar_crit))

# ---- repeat on strong-IV subset only -----------------------------------------
strong_idx <- which((b_x / se_x)^2 >= 10)
cat(sprintf("\n--- mrAR_multi on strong-IV subset (F >= 10, K_strong = %d) ---\n",
            length(strong_idx)))
if (length(strong_idx) >= 1L) {
  if (length(strong_idx) == 1L) {
    # Single instrument: use mrAR() (closed-form Fieller)
    j <- strong_idx[1]
    res_strong <- mrAR(b_x = b_x[j], se_x = se_x[j],
                       b_y = b_y[j], se_y = se_y[j], alpha = 0.05)
    cat("  K = 1 — using mrAR() instead of mrAR_multi():\n")
    cat("  ci_type:    ", res_strong$ci_type, "\n")
    cat(sprintf("  ci:         [%g, %g]\n", res_strong$ci_lower, res_strong$ci_upper))
    cat(sprintf("  J_pvalue:   NA (K = 1)\n"))
    # Wrap in compatible structure
    res_strong_unified <- list(
      ci_type = res_strong$ci_type,
      ci_intervals = list(c(res_strong$ci_lower, res_strong$ci_upper)),
      beta_hat = b_y[j] / b_x[j],
      J_stat = NA_real_, J_pvalue = NA_real_, K = 1L
    )
  } else {
    res_strong <- mrAR_multi(
      b_x = b_x[strong_idx], se_x = se_x[strong_idx],
      b_y = b_y[strong_idx], se_y = se_y[strong_idx], alpha = 0.05)
    cat("  ci_type:    ", res_strong$ci_type, "\n")
    cat("  ci_intervals:\n")
    for (ix in seq_along(res_strong$ci_intervals)) {
      iv <- res_strong$ci_intervals[[ix]]
      cat(sprintf("    [%g, %g]\n", iv[1], iv[2]))
    }
    cat(sprintf("  beta_hat:   %.4f\n", res_strong$beta_hat))
    cat(sprintf("  J_stat:     %.4f  (df = %d)\n", res_strong$J_stat, length(strong_idx) - 1L))
    cat(sprintf("  J_pvalue:   %.4g\n", res_strong$J_pvalue))
    res_strong_unified <- res_strong
    res_strong_unified$K <- length(strong_idx)
  }
} else {
  res_strong_unified <- NULL
  cat("  No strong-IV (F >= 10) instruments in panel.\n")
}

# ---- iv_partial_r2 + e_value per IV (rvSMR sensitivity scalars) -------------
cat("\n--- per-IV sensitivity scalars ---\n")
sens_rows <- list()
for (k in seq_len(K)) {
  r2 <- iv_partial_r2(b_x = b_x[k], se_x = se_x[k], n = panel$n_x)
  # E-value uses the Wald per-IV (rough proxy; the panel-level E-value uses
  # res$beta_hat below).
  wald_per <- b_y[k] / b_x[k]
  se_wald_per <- sqrt((se_y[k]^2) / (b_x[k]^2) +
                      (b_y[k]^2 * se_x[k]^2) / (b_x[k]^4))
  ev <- e_value(b_xy = wald_per, se_xy = se_wald_per)
  sens_rows[[k]] <- list(
    k = k, rsid = rsid[k],
    partial_r2 = r2$partial_r2, RV = r2$RV, F = (b_x[k]/se_x[k])^2,
    e_value_point = ev$e_value_point, e_value_ci = ev$e_value_ci,
    wald_per = wald_per, wald_se_per = se_wald_per
  )
}
sens_df <- do.call(rbind, lapply(sens_rows, function(r) {
  data.frame(k = r$k, rsid = r$rsid, F = round(r$F, 2),
             partial_r2 = signif(r$partial_r2, 3), RV = round(r$RV, 3),
             E_point = round(r$e_value_point, 2),
             E_ci = round(r$e_value_ci, 2), stringsAsFactors = FALSE)
}))
print(sens_df, row.names = FALSE)

# ---- panel-level E-value at mrAR_multi point estimate -----------------------
ev_panel <- tryCatch(
  e_value(b_xy = res$beta_hat,
          se_xy = (res$ci_intervals[[1]][2] - res$ci_intervals[[1]][1]) /
                  (2 * qnorm(0.975))),
  error = function(e) NULL
)
if (!is.null(ev_panel)) {
  cat("\n--- panel-level (mrAR_multi point) E-value ---\n")
  cat(sprintf("  beta_hat (logOR-CHD per SD-PCSK9-expr): %.4f\n", res$beta_hat))
  cat(sprintf("  implied RR (Swanson-VanderWeele):       %.3f\n", ev_panel$rr_scale))
  cat(sprintf("  E_value at point:                       %.2f\n", ev_panel$e_value_point))
  cat(sprintf("  E_value at CI bound nearest null:       %.2f\n", ev_panel$e_value_ci))
}

# ---- per-tissue cell_type_q analog ------------------------------------------
# Substitution rationale: TenK10K Phase 1 28-PBMC-cell-type substrate is
# unavailable (rare-variant Zenodo placeholders, common-variant files 14-23
# GB each but PCSK9 is hepatocyte-not-PBMC); we substitute *bulk GTEx
# tissues* as analog "cell types". This is bulk-tissue resolution, not
# single-cell, but exercises the cell_type_q() plumbing identically. The
# substitution will be flagged in the report.
#
# For each tissue, take the per-tissue cis-eQTL beta on the lead variant
# rs471705 (the strongest PCSK9 liver eQTL) and form the per-tissue Wald
# ratio with the FinnGen I9_IHD beta. Run cell_type_q() across the tissues.
cat("\n--- cell_type_q analog: per-tissue Wald ratios on rs471705 ---\n")
tissues_raw <- fromJSON("pcsk9_per_tissue.json", simplifyDataFrame = FALSE)
lead_rsid <- "rs471705"
lead_idx_in_panel <- which(tab$rsid == lead_rsid)
if (length(lead_idx_in_panel) == 1L) {
  by_lead   <- b_y[lead_idx_in_panel]
  sey_lead  <- se_y[lead_idx_in_panel]
  ct_input  <- list()
  ct_rows   <- list()
  for (tlab in names(tissues_raw)) {
    tis <- tissues_raw[[tlab]]
    lk_list <- tis$lead_lookups
    # find the row whose rsid matches lead
    hit <- NULL
    for (lk in lk_list) {
      if (isTRUE(lk$rsid == lead_rsid) && isTRUE(lk$found) &&
          !is.null(lk$beta) && !is.null(lk$se)) {
        hit <- lk; break
      }
    }
    if (is.null(hit)) {
      ct_rows[[length(ct_rows) + 1L]] <- data.frame(
        cell_type = tlab, n_donors = tis$n_donors,
        eqtl_beta = NA_real_, eqtl_se = NA_real_,
        wald = NA_real_, wald_se = NA_real_,
        included = FALSE, stringsAsFactors = FALSE)
      next
    }
    # nominal eQTL p < 0.05 filter so we only use tissues where PCSK9 is
    # actually expression-regulated by rs471705
    if (!isTRUE(hit$pvalue < 0.05)) {
      ct_rows[[length(ct_rows) + 1L]] <- data.frame(
        cell_type = tlab, n_donors = tis$n_donors,
        eqtl_beta = hit$beta, eqtl_se = hit$se,
        wald = NA_real_, wald_se = NA_real_,
        included = FALSE, stringsAsFactors = FALSE)
      next
    }
    wald   <- by_lead / hit$beta
    se_wald <- sqrt((sey_lead^2) / (hit$beta^2) +
                    (by_lead^2 * hit$se^2) / (hit$beta^4))
    ct_input[[tlab]] <- list(cell_type = tlab, b_xy = wald,
                             se_xy = se_wald, n_donors = tis$n_donors)
    ct_rows[[length(ct_rows) + 1L]] <- data.frame(
      cell_type = tlab, n_donors = tis$n_donors,
      eqtl_beta = hit$beta, eqtl_se = hit$se,
      wald = wald, wald_se = se_wald,
      included = TRUE, stringsAsFactors = FALSE)
  }
  ct_table <- do.call(rbind, ct_rows)
  ct_table$wald    <- round(ct_table$wald, 3)
  ct_table$wald_se <- round(ct_table$wald_se, 3)
  ct_table$eqtl_beta <- round(ct_table$eqtl_beta, 4)
  ct_table$eqtl_se   <- round(ct_table$eqtl_se,   4)
  print(ct_table, row.names = FALSE)
  if (length(ct_input) >= 2L) {
    cell_q <- cell_type_q(ct_input, min_donors = 0)
    cat(sprintf("\n  cell_type_q (C_used = %d): Q = %.3f  df = %d  p = %.4g  -> %s\n",
                cell_q$C_used, cell_q$Q_cell, cell_q$df, cell_q$p_value,
                cell_q$interpretation))
  } else {
    cell_q <- NULL
    cat("\n  cell_type_q: insufficient tissues with significant lead eQTL\n")
  }
} else {
  ct_table <- NULL
  cell_q   <- NULL
  cat("  lead variant rs471705 not in panel — skipping cell_type_q analog\n")
}

# ---- save full output --------------------------------------------------------
out <- list(
  gene = panel$gene, ensg = panel$ensg,
  outcome = panel$outcome,
  exposure_substrate = panel$exposure_substrate,
  outcome_substrate = panel$outcome_substrate,
  K = K, n_x = panel$n_x,
  instruments = tab,
  mrAR_multi = list(
    ci_type = res$ci_type,
    ci_intervals = res$ci_intervals,
    beta_hat = res$beta_hat,
    J_stat = res$J_stat, J_pvalue = res$J_pvalue, ar_crit = res$ar_crit
  ),
  mrAR_strong_subset = res_strong_unified,
  sensitivity_per_IV = sens_df,
  e_value_panel = ev_panel,
  cell_type_q_table = ct_table,
  cell_type_q = if (!is.null(cell_q)) list(
    Q_cell = cell_q$Q_cell, df = cell_q$df, p_value = cell_q$p_value,
    bar = cell_q$bar, C_used = cell_q$C_used,
    interpretation = cell_q$interpretation
  ) else NULL
)
write(toJSON(out, auto_unbox = TRUE, pretty = TRUE),
      "track2_results.json")
cat("\nWrote track2_results.json\n")
