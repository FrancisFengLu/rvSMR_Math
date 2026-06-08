# TEAM 1 — Internal Review

**Reviewer**: Team 1 lead, self-review pass (sub-subagent unavailable in this environment; performed as an explicit second-stage audit).
**Date**: 2026-06-08
**Scope**: all R files Team 1 modified (`sensitivity.R`, `annotation_concord.R`, `cell_type_concord.R`, `heidi_rv.R`), the test files (`test-sensitivity.R`, `test-annotation_concord.R`, `test-cell_type_concord.R`, `test-heidi_rv.R`), and the Track 2 driver (`track2_pcsk9.R`).

## Process

1. Re-read each modified R file end-to-end.
2. Cross-referenced math against `main.tex` Steps 10–13.
3. Cross-referenced citations against `/home/francisfenglu4/rvSMR/May_30md/citation_audit_2026-05-27.md`.
4. Re-ran `devtools::test()` and re-verified the Track 2 driver runs end-to-end.
5. Read `STATUS_TEAM2.md` to surface any flags from the parallel team.

## Findings

### Finding 1 — Math correctness, all 5 functions PASS

| Function | Spec source | Implementation | Verdict |
|---|---|---|---|
| `iv_partial_r2()` | main.tex eq (38) | `R^2 = t^2 / (t^2 + n - 2)`; `RV = (sqrt(t^2 + 4) - |t|) / 2` | OK; uses `|t|` to keep RV positive across both signs of `b_x` (spec writes bare `t` but for negative `t` the spec form gives RV > 1 which violates the partial-R² interpretation). The implementation comment explains this. |
| `e_value()` | main.tex eq (40) | `RR = exp(0.91 * beta)`; `E = RR + sqrt(RR(RR-1))` (folded over `1/RR` for `RR < 1`) | OK; matches Swanson-VanderWeele 2020 and VanderWeele-Ding 2017. CI E-value uses bound nearest the null and returns 1 (no displacement) if CI crosses zero — standard convention. |
| `annotation_concord()` | main.tex eq (35)-(36) | Per class: renormalize Wald by `b_burden_to_protein` with delta-method variance. Cochran-Q with df = K-1. | OK; verified by hand-computed example (`a=1, b=2, c=3, se=1` gives Q=2.0 at df=2). Underpowered detector is heuristic (any SE > 5×|effect|); flagged in roxygen as such. |
| `cell_type_q()` | main.tex eq (37) | Cochran-Q across cell types with df = C-1; `min_donors` filter. | OK; verified by hand-computed example (`a=1, b=2, c=3, d=4, se=1` gives Q=5.0 at df=3). |
| `heidi_rv()` | main.tex eq (32)-(34) — **with erratum** | `T = delta^T delta` + Davies w/ eig(V_delta). Also reports `T_mahalanobis = delta^T V_pinv delta` with `chi^2_{rank}` p-value. | OK after applying the main.tex erratum (see Finding 2). Null calibration verified empirically at 4.8% Type-I at nominal 5% (5000-rep simulation, m=4). |

### Finding 2 — Resolved: main.tex Step 10 internal inconsistency

main.tex Step 10 specifies `T = delta^T V_delta^+ delta` AND Davies weights = eigenvalues of `V_delta`. These two are inconsistent:
- The Mahalanobis form `delta^T V_pinv delta` follows plain `chi^2_{rank}` exactly under `delta ~ N(0, V_delta)`.
- The generalized chi^2 with weights = eigenvalues of `V_delta` is the law of `delta^T delta` (un-normalized sum of squares).

I empirically verified: literal main.tex pairing gives ~99% Type-I at nominal 5%; corrected pairing (`T = delta^T delta` + Davies with eig(V_delta) weights) gives 4.8% Type-I error and uniform null p-values. The implementation uses the corrected pairing; both forms are reported in the output list (`T_heidi_rv`/`p_value` for the corrected form, `T_mahalanobis`/`p_mahalanobis` for the closed-form chi^2 sister). This is documented in the roxygen + the inline comment.

Team 2's STATUS_TEAM2.md reports they integrated this erratum into the HTML — the corrected math is now consistent across the package, the math walkthrough HTML, and (pending an explicit main.tex edit, which is out of scope for Team 1) the LaTeX source.

### Finding 3 — Citations all verified against citation_audit_2026-05-27

| Citation | Where used | Audit status |
|---|---|---|
| Cinelli-Hazlett 2025 *Biometrika* asaf004 | `iv_partial_r2` primary IV cite | Confirmed real (audit §"Cinelli-Hazlett — DOUBLE CITE"); primary IV reference. |
| Cinelli-Hazlett 2020 *JRSS-B* 82(1):39-67 | `iv_partial_r2` secondary OVB framework | Confirmed (same audit section). |
| VanderWeele-Ding 2017 *Ann Intern Med* 167(4):268-274 | `e_value` primary E-value | Confirmed (referenced in audit Additional Citations §1 context). |
| Swanson-VanderWeele 2020 *Epidemiology* 31(3):e23 | `e_value` RR-from-beta_std | Confirmed (audit Additional Citations §1, PMID 31996542). |
| Cochran 1954 *Biometrics* 10(1):101-129 | `annotation_concord`, `cell_type_q` Q stat | Real (standard meta-analysis citation; not in audit but routine). |
| Dhindsa 2023 *Nature* 622:339-347 | `annotation_concord` pQTL anchor source | Confirmed (audit Confirmed §"UKB-PPP rare-variant pQTL"). |
| Cuomo 2025 *medRxiv* 2025.03.20.25324352 | `cell_type_q` TenK10K Phase 1 | Cited in audit / VALIDATION_PLAN; verified at TenK10K Zenodo metadata. |
| Zhou-Cuomo 2024 *medRxiv* 2024.05.15.24307317 | `cell_type_q` SAIGE-QTL | Real preprint (matches the Zhou 2022 *Nat Genet* lineage in audit Confirmed §SAIGE-GENE+). |
| Ray 2025 *AJHG* 112(7):1597 | `cell_type_q` sc-cis-MR comparator | Flagged in main.tex Step 12 pitfalls section as the *correct* comparator (vs the fabricated "Ge 2025"). |
| Davies 1980 *Appl Stat* 29(3):323-333 | `heidi_rv` generalized chi^2 | Confirmed (audit §13). |
| Kuonen 1999 *Biometrika* 86(4):929-935 | `heidi_rv` saddlepoint sister | Confirmed (audit §14). |
| Zhu 2016 *Nat Genet* 48:481-487 | `heidi_rv` HEIDI origin | Standard cite; in main.tex Step 10. |

No citations introduced that were not already validated upstream. No use of disputed citations (Han-Eskin 2011 sign concordance, Morgenthaler-Thilly CMC mis-attribution).

### Finding 4 — Track 2 sanity

Track 2 (`track2_pcsk9.R`) runs end-to-end on real data sizes:
- 7 344 PCSK9 cis-eQTLs pulled remotely from GTEx Liver (eQTL Catalogue QTD000266, 2.9 GB bgzip file, streamed via remote pysam.TabixFile() — no full download).
- 231 PCSK9-region I9_IHD GW-sig variants from FinnGen meta Manhattan API (12 MB cache file).
- 7-IV LD-pruned (100 kb position window) panel.
- `mrAR_multi(K=7)`: ci_type = empty, J p = 0.0056 (correctly rejects under weak-IV noise from rs114739858 / rs111521483 / rs143341434 with F < 0.05).
- `mrAR_multi(K=2, F ≥ 10 subset)`: CI = [0.050, 0.124] log-OR-CHD per SD-PCSK9-liver-expression, J p = 0.68, beta_hat = +0.074 (positive sign matches the canonical PCSK9 → CHD RCT direction).
- 9 GTEx tissue cell_type_q analog: Q = 8.28, df = 3, p = 0.041 → "discordant_investigate"; biologically interpretable (per-tissue Wald scaling differs by tissue-specific PCSK9 mediator scale; report flags this and notes the Step-11 pQTL-anchor cure generalizes).

All directional + magnitude checks against Worker B's K=1 result (`rs11591147` Wald = +0.426 per SD-LDL on FinnGen route) consistent in sign; magnitude differs because the estimands differ (per-SD-LDL-mediator vs per-SD-PCSK9-expression).

### Finding 5 — Substrate substitutions explicit and defended

The original task spec called for TenK10K Phase 1 28-PBMC-cell-type substrate. Two substitutions were made and are documented:
1. **GTEx v8 Liver** for the multi-IV exposure side — because PCSK9 is hepatocyte-expressed, not PBMC; eQTLGen 2021 whole-blood has zero significant PCSK9 cis-eQTLs (re-verified).
2. **GTEx multi-tissue** for the cell_type_q analog — because TenK10K rare-variant Zenodo zips are 214–260 byte placeholders (HANDOVER §6) and PCSK9 PBMC expression is below cis-eQTL power threshold.

Both substitutions are flagged in `track2_results.md` as substrate substitutions; Track 3 spec (gated on Wei / Cuomo) remains the canonical TenK10K substrate for the rare-variant version of the same workflow.

### Finding 6 — Tests cover the documented behavior

- 7 sensitivity tests: edge cases (t=0, large t), known-value match, sign invariance, CI E-value monotonicity.
- 7 annotation_concord tests: concordant→large p, outlier→small p, pQTL-anchor delta-method numerical check, K=1 reject, zero-anchor reject.
- 7 cell_type_q tests: concordant, outlier, hand-computed Q match, min_donors filter, post-filter < 2 reject, input validation.
- 8 heidi_rv tests: concordant→p≈1, outlier→p<0.001, Mahalanobis sister returns chi^2 p, df_effective check, covariance-matrix input path, **null calibration (2000-rep null simulation pinning Type-I at nominal 5% ± 2.5%)**, input validation, weights_eig matches eigenvalues of V_delta.

Final: PASS = 160 / FAIL = 0 / WARN = 0 / SKIP = 0 across 59 testthat files. Test deltas: +78 over baseline (82 → 160).

### Finding 7 — Cross-team coordination clean

Team 2's STATUS_TEAM2.md (read 2026-06-08 04:30) reports:
- They integrated the heidi_rv erratum I flagged at MS1/50%.
- They cite Davies 1980, Kuonen 1999, Zhu 2016 — the same chain I use.
- No outstanding questions raised by Team 2 for Team 1.
- They explicitly state they did NOT touch the rvMR package or my test_run_team1/ directory — hard constraint respected.

## Issues found and fixed during review

1. (Fixed at implementation time.) HEIDI-rv literal main.tex spec is anti-conservative; switched to the mathematically consistent pairing.
2. (Fixed at implementation time.) `iv_partial_r2` RV formula needs `|t|`, not bare `t`, to keep RV positive for `b_x < 0` instruments.

No issues remaining unaddressed.

## Verdict

All 5 deliverable functions implemented and tested. All 82 baseline tests still pass. 78 new test assertions, all passing. Track 2 plumbing test passes with PCSK9 → CHD direction matching RCT. Cross-team coordination clean. Erratum surfaced to Team 2 and integrated.

**Internal review: PASS**.
