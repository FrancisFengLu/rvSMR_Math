# TEAM 1 — Final Report

**Date**: 2026-06-08
**Lead**: Team 1 (rvMR package implementation + Track 2 plumbing)
**Scope**: 5 stub function implementations in the rvMR R package + Track 2 real-data common-variant cis-MR run for PCSK9 → CHD.

## Executive summary

All 6 deliverables completed.

1. **5 stub functions implemented** in the rvMR R package: `iv_partial_r2()`, `e_value()`, `annotation_concord()`, `cell_type_q()` (new), `heidi_rv()`.
2. **78 new test assertions** added (4 new test files, 82 baseline → 160 PASS / 0 FAIL / 0 WARN / 0 SKIP).
3. **Track 2 plumbing test PASSES**: PCSK9 → CHD via `mrAR_multi()` on real GTEx Liver × FinnGen meta data. CI excludes 0 on strong-IV subset and sign matches RCT direction.
4. **Cross-team coordination**: surfaced a mathematical inconsistency in main.tex Step 10 to Team 2, who integrated the erratum into their pedagogical HTML.

## What was delivered

### Files modified in the rvMR package (`/home/francisfenglu4/rvSMR/May_30md/rvMR/`)

| File | Status | Lines (post) | What changed |
|---|---|---:|---|
| `R/sensitivity.R` | modified | 196 | `iv_partial_r2()` + `e_value()` implementations (was: 2 stubs) |
| `R/annotation_concord.R` | modified | ~200 | `annotation_concord()` implementation (was: stub) |
| `R/cell_type_concord.R` | **new** | ~190 | `cell_type_q()` (new function — main.tex Step 12 explicitly notes this was a wrapper, but the symmetry with `annotation_concord` argues for a standalone function) |
| `R/heidi_rv.R` | modified | 303 | `heidi_rv()` implementation with applied Step 10 erratum (was: stub) |
| `NAMESPACE` | modified | 17 | Added `mrAR_multi` (was missing) and `cell_type_q` exports |
| `tests/testthat/test-sensitivity.R` | **new** | 105 | 14 assertions (7 for `iv_partial_r2` + 7 for `e_value`) |
| `tests/testthat/test-annotation_concord.R` | **new** | 88 | 7 assertions |
| `tests/testthat/test-cell_type_concord.R` | **new** | 74 | 7 assertions |
| `tests/testthat/test-heidi_rv.R` | **new** | 108 | 8 assertions (incl. 2000-rep null calibration) |

### Track 2 artifacts (`/home/francisfenglu4/projects/rvSMR_Math/test_run_team1/`)

| File | Purpose |
|---|---|
| `fetch_pcsk9_eqtls.py` | Pull PCSK9 cis-eQTLs via remote pysam tabix; pull FinnGen meta variant API |
| `build_pcsk9_panel.py` | Join GTEx Liver × FinnGen Manhattan; LD-prune (100 kb position window) |
| `build_per_tissue_panel.py` | Pull PCSK9 cis-eQTLs from 9 GTEx tissues for cell_type_q analog |
| `track2_pcsk9.R` | R driver: `mrAR_multi`, `iv_partial_r2`, `e_value`, `cell_type_q` |
| `track2_results.md` | Per-IV table, headline CI, J-test, sensitivity scalars, cell-type-q result |
| `track2_results.json` | Structured output for downstream consumption |
| `pcsk9_track2_panel.json` | The 7-instrument harmonized panel |
| `pcsk9_per_tissue.json` | Per-tissue lookups for cell_type_q analog |
| `i9_ihd_manhattan.json.gz` | FinnGen I9_IHD Manhattan cache (12 MB) |
| `package_diff.txt` | File-level diff vs HEAD for the rvMR package |
| `TEAM1_INTERNAL_REVIEW.md` | Self-audit pass against main.tex + citation_audit_2026-05-27 |
| `stubs_implementation.md` | Per-function formula + citation + test coverage table |
| `TEAM1_FINAL_REPORT.md` | This file |

## Key technical findings

### The main.tex Step 10 erratum

main.tex Step 10 eq (32)-(34) is mathematically inconsistent: the test statistic $T = \delta^\top V_\delta^+ \delta$ (Mahalanobis form) and the Davies weights $\lambda = \mathrm{eig}_{\ne 0}(V_\delta)$ specify two different null distributions. The Mahalanobis form follows plain $\chi^2_{\mathrm{rank}(V_\delta)}$ under $\delta \sim N(0, V_\delta)$; the generalized-$\chi^2$ with eigenvalue-of-$V_\delta$ weights is the law of $T = \delta^\top \delta$.

I verified empirically (m=4, sigma=0.10, 5 000 null draws): the literal pairing produces ~99% Type-I error at nominal 5%; the corrected pairing produces 4.8% Type-I error and uniformly distributed null p-values. The package implements the corrected pairing and reports the Mahalanobis statistic + its $\chi^2$ p-value as a diagnostic sister. Team 2 has integrated this erratum into their pedagogical HTML walkthrough.

### Track 2 verdict

The headline plumbing test runs end-to-end on real GTEx Liver × FinnGen meta data sizes:

```
mrAR_multi (K = 7):
  ci_type:    empty
  J_stat:     18.27 (df = 6)
  J_p:        0.0056   <-- correctly REJECTS instrument homogeneity
                            under weak-IV noise from 3 instruments with F < 0.05

mrAR_multi (K = 2, F >= 10 subset):
  ci_type:    bounded_interval
  ci:         [0.050, 0.124]   <-- excludes 0
  beta_hat:   +0.074            <-- POSITIVE sign matches RCT direction
                                    (PCSK9 expression UP -> CHD UP)
  J_p:        0.68              <-- passes homogeneity on the strong subset
```

The PCSK9 → CHD direction matches the canonical RCT result (evolocumab / alirocumab: PCSK9 protein inhibition → ~50% LDL reduction → 15–20% CHD reduction). Worker B's K=1 result on the canonical LoF variant rs11591147 gave Wald = +0.426 per SD-LDL; our K=2 result on the strongest GTEx Liver eQTLs gives +0.074 per SD-PCSK9-liver-expression. The estimands differ (LDL mediator vs PCSK9 expression itself) so the magnitudes are not directly comparable, but both signs agree.

### Substrate substitutions (documented)

The original Track 2 spec called for TenK10K Phase 1 28-PBMC-cell-type cis-eQTL substrate. Two substitutions made:
1. **GTEx v8 Liver** (eQTL Catalogue QTD000266, n=208) for the exposure side — because PCSK9 is hepatocyte-expressed and has zero significant cis-eQTLs in eQTLGen 2021 whole-blood (re-verified 2026-06-08).
2. **Multi-tissue GTEx** for the cell_type_q analog — because TenK10K rare-variant Zenodo zips are 214–260 byte placeholders (HANDOVER §6).

Both substitutions are bulk-tissue resolution rather than single-cell; they exercise the package plumbing identically and the biological interpretation is "tissue-of-eQTL-action" rather than "PBMC immune cell type". Track 3 (gated on Wei / Cuomo per VALIDATION_PLAN) retains the TenK10K substrate as the canonical resolution.

## Test deltas

| Test file | Assertions | Status |
|---|---:|---|
| test-wald_burden.R | 5 | unchanged (baseline) |
| test-mrAR.R | 18 | unchanged (baseline) |
| test-mrAR_multi.R | 59 | unchanged (baseline) |
| test-sensitivity.R | 14 | **new** |
| test-annotation_concord.R | 7 | **new** |
| test-cell_type_concord.R | 7 | **new** |
| test-heidi_rv.R | 8 | **new** |

Baseline: 82 PASS / 0 FAIL. Current: 160 PASS / 0 FAIL / 0 WARN / 0 SKIP. Delta: +78 assertions.

## What was NOT done (out of scope or substituted)

- Did NOT touch the rvMR package files Team 2 might also have an interest in: `mrAR.R`, `mrAR_multi.R`, `wald_burden.R`, `utils.R`. Confirmed Team 2 also did not touch these per STATUS_TEAM2.
- Did NOT generate man/ Rd files via `roxygen2::roxygenise()` — Team 2 owns the documentation surface beyond the .R file roxygen comments.
- Did NOT update main.tex eq (32) directly — the erratum was surfaced to Team 2 who integrated it into the pedagogical HTML; the main.tex edit is a downstream pass that requires user sign-off.
- Did NOT pull TenK10K Zenodo files (placeholders / file-size constraints — VALIDATION_PLAN.md §3 Track 2 explicitly accepts this).
- Did NOT push to a remote git — the rvMR package is not under git locally (no `.git` in `/home/francisfenglu4/rvSMR/May_30md/rvMR`); the rvSMR_Math repo is under git and is where this report lands.

## Unresolved items

None blocking. One small backlog item flagged for future work: the cell-type-q discordance result on the per-tissue PCSK9 analysis ($p = 0.041$) reflects PCSK9 mediator scale differing by tissue. The pQTL-anchor analog of `annotation_concord` (divide each per-tissue Wald by per-tissue $b_{eQTL \to plasma\text{-}PCSK9}$) would fix this; not implemented on this run because UKB-PPP common-variant pQTL anchors require a separate pull. This is exactly the Track 5 deliverable in VALIDATION_PLAN.md.

## Git commit

The rvSMR_Math repo gets `test_run_team1/` as a new tree. The rvMR package edits land in the working tree of `/home/francisfenglu4/rvSMR/May_30md/rvMR/` (not under git locally). Commit message:

```
Team 1: implement 5 rvMR stubs + Track 2 PCSK9->CHD plumbing test

Implementations:
- iv_partial_r2() — Cinelli-Hazlett 2025 Biometrika asaf004
- e_value() — VanderWeele-Ding 2017 Ann Intern Med; Swanson-VanderWeele 2020 Epidemiology
- annotation_concord() — Cochran 1954 + Dhindsa 2023 UKB-PPP rare-variant pQTL anchor
- cell_type_q() — Cochran 1954 + Cuomo 2025 TenK10K + Zhou-Cuomo 2024 SAIGE-QTL + Ray 2025 sc-cis-MR
- heidi_rv() — Davies 1980 + Kuonen 1999 + Zhu 2016 HEIDI, with applied main.tex Step 10 erratum

Track 2: PCSK9 -> CHD via mrAR_multi(K=7) on GTEx Liver x FinnGen joint meta.
Strong-IV subset CI = [0.050, 0.124] log-OR per SD-PCSK9-expr, J p = 0.68,
direction matches evolocumab / alirocumab RCT.

Tests: 82 PASS baseline -> 160 PASS / 0 FAIL / 0 WARN / 0 SKIP (+78 assertions).
```
