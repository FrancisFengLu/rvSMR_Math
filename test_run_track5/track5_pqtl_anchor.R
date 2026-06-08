# Track 5: pQTL-anchor renormalization of PCSK9 per-tissue Wald ratios
#
# Hypothesis: Team 1's Track 2 cell_type_q reported p = 0.041 (mild
# discordance) across 4 GTEx tissues for PCSK9 -> I9_IHD. We hypothesize
# the discordance is a "secreted-to-plasma scale" artifact: each
# tissue has a different (b_eQTL -> plasma PCSK9 protein) ratio
# because PCSK9 is secreted by liver and only secondarily by other
# tissues. Dividing each tissue's per-tissue Wald ratio by its
# per-tissue (b_eQTL -> b_pQTL) anchor should collapse the
# heterogeneity.
#
# Algebra (VALIDATION_PLAN Track 5):
#   b_xy_unanchored_tissue = b_y / b_x_tissue
#   anchor_tissue          = b_pqtl_variant / b_x_tissue
#   b_xy_anchored_tissue   = b_xy_unanchored_tissue / anchor_tissue
#                          = b_y / b_pqtl_variant
#
# Note: when the SAME variant is used across all tissues (Team 1's
# rs471705), b_pqtl_variant is identical across tissues, so the
# anchored Wald collapses to a single value and Q = 0 by construction.
# This is the mathematically correct outcome of the universal-anchor
# fallback the spec describes.
#
# Variance: delta method.
#   Var(b_y / b_pqtl) ~ se_y^2 / b_pqtl^2 + b_y^2 * se_pqtl^2 / b_pqtl^4
#
# Data sources (data_pull_log.md):
#   - per-tissue b_x_tissue: test_run_team1/pcsk9_per_tissue.json
#   - b_y per IV: test_run_team1/track2_results.json (rs471705 -> I9_IHD)
#   - b_pqtl: Pott et al. 2024 PCSK9 meta-GWAS (PMID 38491180) for
#     rs11591147 (the canonical PCSK9 cis-pQTL lead). Sun BB 2023
#     UKB-PPP was the preferred source per spec, but is NOT in the
#     GWAS Catalog REST API as of the data pull; Synapse / AWS
#     downloads require authentication. See data_pull_log.md.
#
# Output:
#   - per_tissue_renormalized.json
#   - console table comparing unanchored vs anchored cell_type_q
#
# Hard constraints (per spec):
#   - DO NOT modify rvMR package
#   - DO NOT touch test_run_team1/

suppressPackageStartupMessages({
  library(jsonlite)
  library(rvMR)
})

set.seed(2026)

# ----------------------------------------------------------------------
# Data load
# ----------------------------------------------------------------------
TEAM1 <- "/home/francisfenglu4/projects/rvSMR_Math/test_run_team1"
OUT   <- "/home/francisfenglu4/projects/rvSMR_Math/test_run_track5"

panel <- fromJSON(file.path(TEAM1, "pcsk9_track2_panel.json"),
                  simplifyDataFrame = FALSE)
tissues_raw <- fromJSON(file.path(TEAM1, "pcsk9_per_tissue.json"),
                       simplifyDataFrame = FALSE)
track2_res <- fromJSON(file.path(TEAM1, "track2_results.json"),
                       simplifyDataFrame = FALSE)

# Pull b_y, se_y at rs471705 from Team 1's panel (the lead variant
# Team 1 used for cell_type_q)
lead_rsid <- "rs471705"
inst_lead <- NULL
for (it in panel$instruments) {
  if (it$rsid == lead_rsid) { inst_lead <- it; break }
}
stopifnot(!is.null(inst_lead))

b_y  <- as.numeric(inst_lead$y_beta)
se_y <- as.numeric(inst_lead$y_se)
b_x_liver  <- as.numeric(inst_lead$eqtl_beta)
se_x_liver <- as.numeric(inst_lead$eqtl_se)

cat(sprintf("Anchor variant: %s (PCSK9 cis-eQTL lead, GTEx Liver)\n",
            lead_rsid))
cat(sprintf("  b_y  (I9_IHD)         = %.4f  SE = %.4f\n", b_y, se_y))
cat(sprintf("  b_x  (GTEx liver eQTL)= %.4f  SE = %.4f\n",
            b_x_liver, se_x_liver))

# ----------------------------------------------------------------------
# pQTL anchor (universal anchor, per spec fallback)
# ----------------------------------------------------------------------
# Pott J et al. 2024 Hum Mol Genet (PMID 38491180), PCSK9 sex-stratified
# meta-GWAS (6 European studies, n approx 20k); rs11591147-T:
#   beta = 0.37189 (decrease in PCSK9), SE = 0.0145297, p = 2e-144
# Scale: log-transformed PCSK9 protein measurement (mixed Olink + ELISA
# platforms across studies); interpretable as "per-SD log-PCSK9" not
# raw mass units.
#
# SUBSTITUTION RATIONALE (data_pull_log.md):
#   Spec preferred Sun BB 2023 UKB-PPP (Nature 622:329). That source
#   is NOT indexed in the GWAS Catalog REST API for PCSK9 protein
#   measurement; Synapse syn51365301 and AWS S3 ukbiobank.opendata.
#   sagebase.org return 403 (auth required) for bulk-download paths.
#   Pott 2024 is the closest available substitute (also Olink scale,
#   European-ancestry, does not overlap UKB but has UKB-PPP-compatible
#   units).
b_pqtl_RAW <- 0.37189   # rs11591147-T, decrease
se_pqtl    <- 0.0145297
# Convention: Pott 2024 reported as "decrease in PCSK9 per T allele".
# T is the LOF allele; the magnitude is what we need for renormalization.
# Sign-coherent representation: T allele LOWERS PCSK9 (negative effect)
# AND lowers LDL/CHD risk. So b_pqtl on PCSK9 expression scale (same
# sign as eQTL) is +0.37189 if we represent rs11591147 as "alt allele
# T raises PCSK9" (which is FALSE biologically -- T lowers it), or
# -0.37189 if we represent "T raises PCSK9".
#
# For the spec's algebra (tilde_b_xy = b_y / b_pqtl), the sign of
# b_pqtl matters only to the sign of the renormalized Wald. We use
# |b_pqtl| = 0.37189 here because:
#   (a) the discordance test is in magnitude (Q is sum-of-squares),
#   (b) the canonical eQTL direction (rs471705-alt raises PCSK9
#       expression, beta positive) is the same biological direction
#       as rs11591147-T LOWERS PCSK9 (LOF) -- so when we anchor
#       cross-variant we flip sign accordingly.
# This is documented in WORKER_TRACK5_REPORT.md as a coarse
# normalization choice.
b_pqtl <- abs(b_pqtl_RAW)

cat(sprintf("\npQTL anchor (Pott 2024 PMID 38491180; rs11591147-T):\n"))
cat(sprintf("  |b_pqtl| = %.4f  SE = %.4f\n", b_pqtl, se_pqtl))
cat(sprintf("  scale: log-PCSK9 protein measurement (mixed Olink + ELISA)\n"))

# ----------------------------------------------------------------------
# Per-tissue Wald (unanchored = Team 1 reproduce) and anchored
# ----------------------------------------------------------------------
TISSUES_INCLUDED <- c("liver", "blood", "adipose_visceral", "artery_aorta")
# Team 1 used these 4 tissues (those where rs471705 has nominal eQTL
# p < 0.05). We reproduce that filter exactly so the unanchored Q
# matches Team 1's p = 0.041.

rows <- list()
for (tlab in names(tissues_raw)) {
  tis <- tissues_raw[[tlab]]
  hit <- NULL
  for (lk in tis$lead_lookups) {
    if (isTRUE(lk$rsid == lead_rsid) && isTRUE(lk$found) &&
        !is.null(lk$beta) && !is.null(lk$se)) {
      hit <- lk; break
    }
  }
  if (is.null(hit)) {
    next
  }
  included <- (tlab %in% TISSUES_INCLUDED) && isTRUE(hit$pvalue < 0.05)
  b_x_t  <- as.numeric(hit$beta)
  se_x_t <- as.numeric(hit$se)

  # Unanchored Wald (Team 1)
  wald_un    <- b_y / b_x_t
  se_wald_un <- sqrt(se_y^2 / b_x_t^2 + b_y^2 * se_x_t^2 / b_x_t^4)

  # Anchored Wald: tilde_b_xy = b_y / b_pqtl
  # (universal anchor; same value for all tissues)
  wald_anch  <- b_y / b_pqtl
  # delta method: Var(b_y / b_pqtl) = (1/b_pqtl)^2 * Var(b_y) +
  #              (b_y / b_pqtl^2)^2 * Var(b_pqtl)
  se_wald_anch <- sqrt(se_y^2 / b_pqtl^2 +
                       b_y^2 * se_pqtl^2 / b_pqtl^4)

  rows[[length(rows) + 1L]] <- list(
    tissue        = tlab,
    n_donors      = tis$n_donors,
    b_x_tissue    = b_x_t,
    se_x_tissue   = se_x_t,
    eqtl_p        = hit$pvalue,
    b_y           = b_y,
    se_y          = se_y,
    b_pqtl        = b_pqtl,
    se_pqtl       = se_pqtl,
    wald_unanchored    = wald_un,
    se_wald_unanchored = se_wald_un,
    wald_anchored      = wald_anch,
    se_wald_anchored   = se_wald_anch,
    included      = included
  )
}

# ----------------------------------------------------------------------
# Per-tissue table
# ----------------------------------------------------------------------
df <- do.call(rbind, lapply(rows, function(r) data.frame(
  tissue     = r$tissue,
  n_donors   = r$n_donors,
  b_x_tissue = round(r$b_x_tissue, 4),
  eqtl_p     = signif(r$eqtl_p, 3),
  wald_un    = round(r$wald_unanchored, 4),
  se_un      = round(r$se_wald_unanchored, 4),
  wald_anch  = round(r$wald_anchored, 4),
  se_anch    = round(r$se_wald_anchored, 4),
  included   = r$included,
  stringsAsFactors = FALSE
)))

cat("\n--- per-tissue table (anchor variant rs471705, universal pQTL anchor rs11591147 / Pott 2024) ---\n")
print(df, row.names = FALSE)

# ----------------------------------------------------------------------
# cell_type_q on unanchored and anchored
# ----------------------------------------------------------------------
ct_input_un   <- list()
ct_input_anch <- list()
for (r in rows) {
  if (!r$included) next
  ct_input_un[[r$tissue]] <- list(
    cell_type = r$tissue,
    b_xy      = r$wald_unanchored,
    se_xy     = r$se_wald_unanchored,
    n_donors  = r$n_donors
  )
  ct_input_anch[[r$tissue]] <- list(
    cell_type = r$tissue,
    b_xy      = r$wald_anchored,
    se_xy     = r$se_wald_anchored,
    n_donors  = r$n_donors
  )
}

cat(sprintf("\nIncluded tissues for cell_type_q: %d (target: 4 to match Team 1)\n",
            length(ct_input_un)))
stopifnot(length(ct_input_un) >= 2L)

cq_un <- cell_type_q(ct_input_un, min_donors = 0L)
cat(sprintf("\nUNANCHORED cell_type_q (reproduce Team 1):\n"))
cat(sprintf("  C_used = %d  Q = %.4f  df = %d  p = %.4g  -> %s\n",
            cq_un$C_used, cq_un$Q_cell, cq_un$df, cq_un$p_value,
            cq_un$interpretation))

cq_anch <- cell_type_q(ct_input_anch, min_donors = 0L)
cat(sprintf("\nANCHORED cell_type_q (pQTL-anchor universal renorm):\n"))
cat(sprintf("  C_used = %d  Q = %.4f  df = %d  p = %.4g  -> %s\n",
            cq_anch$C_used, cq_anch$Q_cell, cq_anch$df, cq_anch$p_value,
            cq_anch$interpretation))

# ----------------------------------------------------------------------
# SECONDARY ANALYSIS: per-tissue best variant + universal pQTL anchor
# ----------------------------------------------------------------------
# This is the spec's "use the same plasma PCSK9 pQTL lead variant
# (rs11591147) as the universal anchor across all tissues" fallback,
# combined with PER-TISSUE selection of the strongest eQTL variant.
#
# For each tissue, pick the panel variant with largest F-stat
# (b_x_tissue / se_x_tissue)^2. Use that variant's b_y and b_x_tissue
# from Team 1's panel. Renormalize each tissue's Wald by the SAME
# pQTL anchor b_pqtl. Because each tissue's b_y differs across variants
# (panel variants have different y_beta), this produces a NON-DEGENERATE
# anchored Wald distribution.
cat("\n--- SECONDARY ANALYSIS: per-tissue best variant + universal pQTL anchor ---\n")

# Map of panel variants -> b_y, se_y for quick lookup
panel_y <- list()
for (it in panel$instruments) {
  panel_y[[it$rsid]] <- list(
    b_y  = as.numeric(it$y_beta),
    se_y = as.numeric(it$y_se),
    pos  = as.integer(it$pos)
  )
}

rows_v2 <- list()
for (tlab in names(tissues_raw)) {
  tis <- tissues_raw[[tlab]]
  best <- NULL
  for (lk in tis$lead_lookups) {
    if (!isTRUE(lk$found)) next
    if (is.null(lk$beta) || is.null(lk$se)) next
    bx <- as.numeric(lk$beta); sx <- as.numeric(lk$se)
    if (!is.finite(bx) || !is.finite(sx) || sx <= 0) next
    F_stat <- (bx / sx)^2
    if (is.null(best) || F_stat > best$F_stat) {
      best <- list(rsid = lk$rsid, b_x = bx, se_x = sx,
                   eqtl_p = as.numeric(lk$pvalue), F_stat = F_stat)
    }
  }
  if (is.null(best)) next
  py <- panel_y[[best$rsid]]
  if (is.null(py)) next
  b_y_v <- py$b_y; se_y_v <- py$se_y

  # Include tissue if best variant has eQTL p < 0.05 (matches Team 1 filter)
  included <- isTRUE(best$eqtl_p < 0.05)

  wald_un      <- b_y_v / best$b_x
  se_wald_un   <- sqrt(se_y_v^2 / best$b_x^2 +
                       b_y_v^2 * best$se_x^2 / best$b_x^4)
  # Anchored: variant-specific b_y / universal b_pqtl
  wald_anch    <- b_y_v / b_pqtl
  se_wald_anch <- sqrt(se_y_v^2 / b_pqtl^2 +
                       b_y_v^2 * se_pqtl^2 / b_pqtl^4)

  rows_v2[[length(rows_v2) + 1L]] <- list(
    tissue    = tlab,
    best_rsid = best$rsid,
    n_donors  = tis$n_donors,
    b_x       = best$b_x,
    se_x      = best$se_x,
    eqtl_p    = best$eqtl_p,
    F_stat    = best$F_stat,
    b_y       = b_y_v,
    se_y      = se_y_v,
    wald_un       = wald_un,
    se_wald_un    = se_wald_un,
    wald_anch     = wald_anch,
    se_wald_anch  = se_wald_anch,
    included      = included
  )
}

df2 <- do.call(rbind, lapply(rows_v2, function(r) data.frame(
  tissue = r$tissue, best_rsid = r$best_rsid,
  n_donors = r$n_donors, F_stat = round(r$F_stat, 2),
  eqtl_p = signif(r$eqtl_p, 3),
  b_y = round(r$b_y, 4),
  wald_un = round(r$wald_un, 4), se_un = round(r$se_wald_un, 4),
  wald_anch = round(r$wald_anch, 4), se_anch = round(r$se_wald_anch, 4),
  included = r$included,
  stringsAsFactors = FALSE
)))
print(df2, row.names = FALSE)

ct_input_un2   <- list()
ct_input_anch2 <- list()
for (r in rows_v2) {
  if (!r$included) next
  ct_input_un2[[r$tissue]] <- list(
    cell_type = paste0(r$tissue, ":", r$best_rsid),
    b_xy = r$wald_un, se_xy = r$se_wald_un, n_donors = r$n_donors)
  ct_input_anch2[[r$tissue]] <- list(
    cell_type = paste0(r$tissue, ":", r$best_rsid),
    b_xy = r$wald_anch, se_xy = r$se_wald_anch, n_donors = r$n_donors)
}

cq_un2 <- if (length(ct_input_un2) >= 2L)
            cell_type_q(ct_input_un2, min_donors = 0L) else NULL
cq_anch2 <- if (length(ct_input_anch2) >= 2L)
              cell_type_q(ct_input_anch2, min_donors = 0L) else NULL

if (!is.null(cq_un2)) {
  cat(sprintf("\nSEC-UNANCHORED (per-tissue best var): C = %d  Q = %.4f  p = %.4g  -> %s\n",
              cq_un2$C_used, cq_un2$Q_cell, cq_un2$p_value, cq_un2$interpretation))
}
if (!is.null(cq_anch2)) {
  cat(sprintf("SEC-ANCHORED   (per-tissue best var): C = %d  Q = %.4f  p = %.4g  -> %s\n",
              cq_anch2$C_used, cq_anch2$Q_cell, cq_anch2$p_value, cq_anch2$interpretation))
}

# ----------------------------------------------------------------------
# Save outputs
# ----------------------------------------------------------------------
out <- list(
  anchor_variant_eqtl    = lead_rsid,
  anchor_variant_pqtl    = "rs11591147",
  pqtl_source            = list(
    paper = "Pott J et al. 2024",
    pmid  = "38491180",
    journal = "Hum Mol Genet",
    note  = paste("Substitution for Sun BB 2023 UKB-PPP",
                  "(Nature 622:329) which was not accessible",
                  "via public REST API. Pott 2024 is Olink-",
                  "scale-compatible European meta-GWAS,",
                  "n approx 20k; does NOT overlap UKB-PPP",
                  "but uses comparable log-PCSK9 units."),
    b_pqtl_used = b_pqtl,
    b_pqtl_raw  = b_pqtl_RAW,
    se_pqtl     = se_pqtl
  ),
  outcome_b_y       = b_y,
  outcome_se_y      = se_y,
  outcome           = "FinnGen R12 x MVP x UKBB I9_IHD",
  team1_cell_q_unanchored = list(
    Q = cq_un$Q_cell, df = cq_un$df, p_value = cq_un$p_value,
    C_used = cq_un$C_used, interpretation = cq_un$interpretation
  ),
  track5_cell_q_anchored = list(
    Q = cq_anch$Q_cell, df = cq_anch$df, p_value = cq_anch$p_value,
    C_used = cq_anch$C_used, interpretation = cq_anch$interpretation
  ),
  secondary_per_tissue_best_variant = list(
    description = paste("Per-tissue best (max F) panel variant +",
                        "universal pQTL anchor (rs11591147). Tests",
                        "whether variant heterogeneity drives",
                        "discordance separately from tissue-eQTL",
                        "scale heterogeneity."),
    unanchored = if (!is.null(cq_un2)) list(
      Q = cq_un2$Q_cell, df = cq_un2$df, p_value = cq_un2$p_value,
      C_used = cq_un2$C_used, interpretation = cq_un2$interpretation) else NULL,
    anchored = if (!is.null(cq_anch2)) list(
      Q = cq_anch2$Q_cell, df = cq_anch2$df, p_value = cq_anch2$p_value,
      C_used = cq_anch2$C_used, interpretation = cq_anch2$interpretation) else NULL,
    per_tissue = rows_v2
  ),
  per_tissue = rows
)

write(toJSON(out, auto_unbox = TRUE, pretty = TRUE),
      file.path(OUT, "per_tissue_renormalized.json"))
cat(sprintf("\nWrote %s\n", file.path(OUT, "per_tissue_renormalized.json")))

# ----------------------------------------------------------------------
# Interpretation block
# ----------------------------------------------------------------------
cat("\n=== INTERPRETATION ===\n")
cat(sprintf("Team 1 unanchored: p = %.4f (discordant_investigate, Q = %.2f, df = %d)\n",
            cq_un$p_value, cq_un$Q_cell, cq_un$df))
cat(sprintf("Track 5 anchored:  p = %.4f (Q = %.2g, df = %d)\n",
            cq_anch$p_value, cq_anch$Q_cell, cq_anch$df))
cat("\n")
cat("ANCHOR ALGEBRA NOTE: When the universal anchor variant (rs11591147)\n")
cat("is used across all 4 tissues, b_y / b_pqtl is identical for every\n")
cat("tissue -> Q collapses to 0 by construction. This is the MATHEMATICALLY\n")
cat("EXPECTED outcome of the spec's fallback algebra, NOT a finding of\n")
cat("biological homogeneity. The substantive test of the 'secretion-scale\n")
cat("artifact' hypothesis would require per-tissue lead variants WITH\n")
cat("matched per-variant pQTL anchors, which the 7-panel-variant pQTL\n")
cat("lookup did not return (none of the 7 panel variants are reported as\n")
cat("PCSK9 pQTLs in GWAS Catalog -- expected: they are intergenic eQTL\n")
cat("leads, not coding-variant pQTL leads). See WORKER_TRACK5_REPORT.md.\n")
