# WORKER_FINNGEN_REPORT — rvMR cis-MR sanity check on public FinnGen × MVP × UKBB meta

Worker: external real-data check. Date: 2026-06-08. Output dir: `test_run_finngen/`.

## What we did

Pulled K=1 cis-MR summary statistics for 5 RCT-validated lipid drug targets (PCSK9, HMGCR, ANGPTL3, APOC3, LPA) and ran `rvMR::mrAR()` to test whether the 95% Anderson–Rubin CI for the causal ratio recovers the RCT-known direction and (where possible) magnitude.

- **Outcome (b_y)**: FinnGen R12 × MVP × UKBB joint meta-analysis at `https://mvp-ukbb.finngen.fi` (no-login PheWeb). Endpoint `I9_IHD` ("Ischaemic heart disease, wide definition") for PCSK9, HMGCR, APOC3; substituted `I9_MI_STRICT` for ANGPTL3 and LPA because `I9_IHD` is NA at their lead variants in the meta.
- **Exposure (b_x)**: Open Targets Platform GraphQL API (`https://api.platform.opentargets.org/api/v4/graphql`) credible-set rows for LDL-C / triglyceride GWAS at each lead variant. Per-ALT-allele coding throughout; allele harmonization is automatic because both sources key on the same `chr:pos:ref:alt`.
- **rvMR engine**: `rvMR::mrAR(b_x, se_x, b_y, se_y, alpha=0.05)` installed from the local source at `/home/francisfenglu4/rvSMR/May_30md/rvMR` into a user R 4.5 library.

## Substitutions (honest reporting)

| Gene | Originally planned | Used | Reason |
|---|---|---|---|
| APOC3 | rs138326449 (R19X stop-gain) | rs964184 (APOA5–APOC3 cluster GWAS tag) | rs138326449 fg_af ≈ 3e-4 in the FinnGen exome arm; all circulatory endpoints NA in the cross-cohort meta |
| ANGPTL3 | I9_IHD | I9_MI_STRICT | I9_IHD NA at rs10889353 in the meta (no chr1:62.5-62.8 MB variants appear in the I9_IHD manhattan) |
| LPA | I9_IHD; Lp(a) GWAS for b_x | I9_MI_STRICT; LDL GWAS as Lp(a) proxy | I9_IHD NA at rs10455872; no Lp(a)-mg/dL GWAS in Open Targets credible-set rows for this variant |

## Results (5/5 directionally correct, 3/5 magnitudinally containing the published point)

| Gene | F | mrAR 95% CI (per SD exposure → log-OR) | RCT direction match |
|---|---:|---|---|
| PCSK9 | 173 | [0.355, 0.518] | yes |
| HMGCR | 2458 | [0.111, 0.292] | yes |
| ANGPTL3 | 344 | [0.023, 0.282] | yes |
| APOC3 | 287 | [0.219, 0.321] | yes |
| LPA | 70 | [1.117, 1.832] | yes |

All 5 mrAR CIs are bounded (every F is well above the weak-IV regime) and exclude 0. The mrAR upper bound exceeds the Wald upper bound by ~1 % in every row, consistent with the Fieller / AR construction being conservative at finite F. The `ar_at_point_estimate` field is exactly 0 in every row, confirming the AR statistic vanishes at the Wald ratio (built-in sanity check passes).

### What this proves and what it doesn't

- **Proves**: `mrAR()` runs end-to-end on real publicly accessible summary statistics; the closed-form quadratic root finder selects the `bounded` branch correctly when F is strong; the AR CI numerically agrees with Wald in the strong-F limit (as the theory predicts).
- **Does not prove**: weak-IV behaviour (every F here is ≥ 70). To exercise the `disconnected` and `whole_line` branches we would need genuinely weak instruments — outside the scope of the 5-gene RCT panel.

## Failures / honest gaps

- The FinnGen meta-PheWeb endpoint `I9_IHD` is NA for ANGPTL3 and LPA lead variants — not an rvMR limitation, a meta-coverage filter at the PheWeb side.
- The canonical APOC3 R19X stop-gain was not callable in the cross-cohort meta; we used the common-tag substitute.
- Lp(a) exposure GWAS is not in the Open Targets credible-set rows for rs10455872. The LDL proxy yields the right sign but an inflated magnitude (because Lp(a)-cholesterol is a fraction of LDL-C), so the LPA point estimate of 1.39 is not directly comparable to a Lp(a)-anchored published cis-MR estimate.
- I9_IHD is a wider endpoint than MI-strict; the PCSK9 mrAR CI [0.355, 0.518] is shifted slightly low relative to the canonical Ference 2016 NEJM per-SD-LDL log-OR of ~0.5–0.7, which is the typical wide-IHD-vs-MI dilution.

## Files produced

- `finngen_panel_results.md` — main table
- `finngen_pull_log.md` — URL/data-source log
- `finngen_cis_mr.R` — analysis script
- `finngen_results.rds`, `finngen_results.json` — mrAR output
- `panel_input.json` — harmonized 5-gene input
- `per_gene/*.json` — raw PheWeb and Open Targets dumps for reproducibility

## Verdict sentence

**rvMR's K=1 mrAR on the 5-gene panel matches RCT direction for 5/5 genes, with CI containing the published cis-MR estimate for 3/5 genes (HMGCR, ANGPTL3, APOC3); PCSK9 is in-direction but slightly low due to I9_IHD being broader than MI, and LPA is in-direction but uses an LDL proxy rather than Lp(a) so its numerical CI is not directly comparable to the Lp(a)-anchored published estimate.**
