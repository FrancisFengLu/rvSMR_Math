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

**Team 2: please flag this in your HTML if you reference Step 10.** The literal main.tex eq (32-34) needs an erratum: change `T = delta^T V_delta^+ delta` to `T = delta^T delta` (or change the Davies weights from "eig_nz(V_delta)" to "all 1s" if you keep the Mahalanobis form, but then the test reduces to standard chi^2 and Davies is gratuitous). Documented in `heidi_rv.R` roxygen.

### Citations used (all verified against citation_audit_2026-05-27)

- iv_partial_r2: Cinelli-Hazlett 2025 *Biometrika* asaf004; Cinelli-Hazlett 2020 *JRSS-B* 82(1):39-67.
- e_value: VanderWeele-Ding 2017 *Ann Intern Med* 167(4):268-274; Swanson-VanderWeele 2020 *Epidemiology* 31(3):e23-e24.
- annotation_concord: Cochran 1954 *Biometrics* 10(1):101-129; Dhindsa 2023 *Nature* 622:339-347 (UKB-PPP rare-variant pQTL anchor).
- cell_type_q: Cochran 1954; Cuomo 2025 *medRxiv* 2025.03.20.25324352 (TenK10K Phase 1); Zhou-Cuomo 2024 *medRxiv* 2024.05.15.24307317 (SAIGE-QTL); Ray 2025 *AJHG* 112(7):1597 (cell-type cis-MR precedent).
- heidi_rv: Davies 1980 *Appl Stat* 29(3):323-333; Kuonen 1999 *Biometrika* 86(4):929-935; Zhu 2016 *Nat Genet* 48:481-487 (HEIDI origin).

### Next

Track 2 (PCSK9 → CHD via common-variant cis-eQTL × FinnGen meta) starting now. Worker B already established the K=1 path; I extend to K = multiple cis-eQTLs and run `mrAR_multi()` with Sargan-J.

No questions for Team 2 yet (the Step 10 erratum above is the main flag).
