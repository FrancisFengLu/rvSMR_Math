# TEAM 1 STATUS

## Milestone 0 / 0% — Bootstrap (2026-06-08)

Read scope:
- `main.tex` Steps 10-13 — math walkthrough for HEIDI-rv, annotation-class Cochran-Q, cell-type Cochran-Q, IV partial R^2, E-value.
- `VALIDATION_PLAN.md` §3 Track 2 — common-variant cis-MR plumbing test.
- `test_run_finngen/finngen_panel_results.md` — Worker B's K=1 mrAR ran clean; PCSK9 mrAR CI = [0.355, 0.518] per SD LDL, sign matches RCT.
- `citation_audit_2026-05-27.md` — 13/19 cites confirmed; primary IV sensitivity = Cinelli-Hazlett 2025 Biometrika asaf004.

## Milestone 1 / 50% — All 5 stubs implemented + tests passing (2026-06-08)

### Implementations done

| Function | File | Lines added | Tests added | Notes |
|---|---|---|---|---|
| `iv_partial_r2()` | sensitivity.R | ~50 | 7 | Cinelli-Hazlett 2025 IV partial R^2 + RV |
| `e_value()` | sensitivity.R | ~55 | 7 | Swanson-VanderWeele 2020 + VanderWeele-Ding 2017 |
| `annotation_concord()` | annotation_concord.R | ~135 | 7 | Cochran 1954 with pQTL-anchor delta-method |
| `cell_type_q()` | **new** cell_type_concord.R | ~190 | 7 | Cochran 1954 across cell types |
| `heidi_rv()` | heidi_rv.R | ~135 | 8 | Davies 1980 generalized chi^2 |

Total: 82 baseline tests -> 160 tests (+78); FAIL = 0, WARN = 0, SKIP = 0.

### Math choice flagged for cross-team review

**main.tex Step 10 contains a mathematical inconsistency** in the HEIDI-rv specification. The text states:

> T = delta^T V_delta^+ delta ... follows a generalized chi^2 with weights = non-zero eigenvalues of V_delta itself.

These two statements cannot both be true. The Mahalanobis form T_M = delta^T V_delta^+ delta exactly follows chi^2_{rank(V_delta)} when delta ~ N(0, V_delta) — NOT a non-trivial generalized chi^2. The un-normalized form T = delta^T delta has the law sum_i lambda_i Z_i^2 where lambda_i are eigenvalues of V_delta — THIS is the generalized chi^2 form whose Davies weights are the eigenvalues of V_delta itself.

**Empirical verification (2026-06-08)**: simulating delta ~ N(0, V_delta) for m = 4 and applying the literal-spec pairing (T = delta^T V^+ delta + Davies w/ eigenvalues of V_delta) gives ~99% Type-I error at nominal 5%. The corrected pairing (T = delta^T delta + Davies w/ eigenvalues of V_delta) gives 4.8% Type-I error and uniformly distributed null p-values, as expected.

**Resolution**: I implemented T = delta^T delta with Davies weights = eigenvalues of V_delta — the pairing that aligns with Davies 1980 and matches the spec's cited Davies call. The Mahalanobis statistic and its chi^2_{m-1} p-value are also reported as a diagnostic sister (`T_mahalanobis`, `p_mahalanobis`).

**Team 2: please flag this in your HTML if you reference Step 10.** Documented in `heidi_rv.R` roxygen.

### Citations used (all verified against citation_audit_2026-05-27)

- iv_partial_r2: Cinelli-Hazlett 2025 *Biometrika* asaf004; Cinelli-Hazlett 2020 *JRSS-B* 82(1):39-67.
- e_value: VanderWeele-Ding 2017 *Ann Intern Med* 167(4):268-274; Swanson-VanderWeele 2020 *Epidemiology* 31(3):e23-e24.
- annotation_concord: Cochran 1954 *Biometrics* 10(1):101-129; Dhindsa 2023 *Nature* 622:339-347 (UKB-PPP rare-variant pQTL anchor).
- cell_type_q: Cochran 1954; Cuomo 2025 *medRxiv* 2025.03.20.25324352 (TenK10K Phase 1); Zhou-Cuomo 2024 *medRxiv* 2024.05.15.24307317 (SAIGE-QTL); Ray 2025 *AJHG* 112(7):1597 (cell-type cis-MR precedent).
- heidi_rv: Davies 1980 *Appl Stat* 29(3):323-333; Kuonen 1999 *Biometrika* 86(4):929-935; Zhu 2016 *Nat Genet* 48:481-487 (HEIDI origin).

## Milestone 2 / 100% — Track 2 + internal review complete (2026-06-08)

### Track 2 (PCSK9 → CHD via common-variant cis-MR)

Pipeline:
1. Pull 7 344 PCSK9 cis-eQTL records from GTEx v8 Liver (eQTL Catalogue QTD000266, 2.9 GB bgzip file) via remote `pysam.TabixFile()` slicing on the PCSK9 ±1 Mb cis window — no full download required.
2. Pull 231 PCSK9-region I9_IHD GW-sig variants from FinnGen R12 × MVP × UKBB joint meta Manhattan API.
3. Join by (chr, pos, ref, alt); LD-prune at 100 kb position window → 7 lead instruments.
4. Run `mrAR_multi(K = 7)`, `mrAR_multi(K = 2, F ≥ 10 subset)`, per-IV `iv_partial_r2()`, per-IV `e_value()`, panel-level `e_value()`, and a 9-tissue `cell_type_q()` analog.

Results:
- `mrAR_multi(K = 7)`: ci_type = empty, J = 18.27 (df = 6), J p = 0.0056 → correctly rejects homogeneity under weak-IV noise.
- `mrAR_multi(K = 2, F ≥ 10)`: bounded CI [0.050, 0.124] log-OR-CHD per SD-PCSK9-liver-expression, J p = 0.68, beta_hat = +0.074 → **excludes 0, positive sign matches RCT direction (PCSK9 expression UP → CHD UP, consistent with evolocumab/alirocumab)**.
- cell_type_q across 9 GTEx tissues (4 with significant lead-eQTL): Q = 8.28, df = 3, p = 0.041 → biologically interpretable cross-tissue Wald-scale discordance (PCSK9 mediator is secreted, scale varies by tissue).

Substrate substitutions documented in `test_run_team1/track2_results.md`:
1. GTEx Liver for the multi-IV exposure side (PCSK9 hepatocyte not PBMC; eQTLGen has zero significant PCSK9 cis-eQTLs).
2. GTEx multi-tissue for cell_type_q analog (TenK10K rare-variant Zenodo placeholders; common-variant 14-23 GB; PCSK9 PBMC expression too low).

### Internal review

Self-audit pass (sub-subagent unavailable in this environment, performed as explicit second-stage audit) wrote `TEAM1_INTERNAL_REVIEW.md`. Re-read each modified .R file end-to-end, cross-checked math against main.tex Steps 10-13, cross-checked every citation against citation_audit_2026-05-27 (all 13/19 confirmed citations used; no disputed citations introduced), re-ran `devtools::test()` (160 PASS / 0 FAIL), re-ran Track 2 driver, read STATUS_TEAM2.md.

Issues caught and fixed during review: zero — all issues were caught + fixed at implementation time.

### Team 2 cross-check

Read STATUS_TEAM2.md at MS2: Team 2 delivered `algorithm_paper_walkthrough.html` (1900 lines, 6 inline SVGs, 19 citation cards) and **integrated the Step 10 erratum I flagged at MS1**. Team 2 reports no outstanding questions for Team 1. Hard constraints respected on both sides: Team 2 did not touch the rvMR package, Team 1 did not touch the HTML.

## Final deliverable structure (in `test_run_team1/`)

```
test_run_team1/
├── TEAM1_FINAL_REPORT.md        # 2-page executive summary
├── TEAM1_INTERNAL_REVIEW.md     # internal audit pass
├── stubs_implementation.md       # per-function formula + citation + test coverage
├── package_diff.txt              # file-level diff vs HEAD for rvMR
├── track2_pcsk9.R                # the Track 2 R driver
├── track2_results.md             # Track 2 narrative + per-IV table
├── track2_results.json           # Track 2 structured output
├── fetch_pcsk9_eqtls.py          # GTEx Liver tabix + FinnGen variant API
├── build_pcsk9_panel.py          # GTEx Liver × FinnGen Manhattan join + LD prune
├── build_per_tissue_panel.py     # 9-tissue panel for cell_type_q analog
├── pcsk9_track2_panel.json       # 7-IV harmonized panel
├── pcsk9_per_tissue.json         # per-tissue lookup data
├── pcsk9_gtex_liver_cis.json     # 7344 raw cis-eQTL records
├── pcsk9_lead_eqtls.json         # 12 LD-pruned pre-FinnGen-join leads
└── i9_ihd_manhattan.json.gz      # FinnGen I9_IHD Manhattan cache (12 MB)
```

## STATUS = COMPLETE
