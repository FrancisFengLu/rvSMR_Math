# rvMR Package Algorithm Audit — Phase 1

*Substrate: `/home/francisfenglu4/rvSMR/May_30md/rvMR/` (R 4.1.2, devtools loadable, 68/68 testthat assertions pass at baseline).*
*Output destination: input for Phase 2 (test-data generator) and Phase 3 (validation runs).*

---

## 1. Implemented functions

| Function | File:line | Signature | What it computes | testthat asserts |
|---|---|---|---|---|
| `validate_summary_input` | `utils.R:32` | `(b_x, se_x, b_y, se_y, n_x, n_y)` | Input sanitizer; matched lengths, positive SE, finite values, positive n | indirect (every caller) |
| `delta_method_ratio_se` | `utils.R:89` | `(b_x, se_x, b_y, se_y, cor_xy = 0)` | First-order delta-method SE of `b_y/b_x` with optional sample-overlap correction | 1 (hand-check in wald test) |
| `f_statistic` | `utils.R:133` | `(b_x, se_x)` | `(b_x/se_x)^2` first-stage F | 1 (in wald test) |
| `wald_burden` | `wald_burden.R:77` | `(b_x, se_x, b_y, se_y, n_x, n_y, cor_xy = 0)` | Wald ratio point estimate + delta SE + F | 7 |
| `mrAR` (K=1 closed-form) | `mrAR.R:88` | `(b_x, se_x, b_y, se_y, alpha = 0.05, cor_xy = 0)` | Inverts the scalar quadratic `A b0^2 + B b0 + C ≤ 0`; classifies CI as `bounded` / `disconnected` / `whole_line` / `empty` | 27 |
| `mrAR_multi` (K≥2 grid+uniroot) | `mrAR_multi.R:114` | `(b_x, se_x, b_y, se_y, R_xx=I_K, R_yy=I_K, R_xy=0, alpha=0.05, n_grid=4000, grid_pad_mult=3, grid_extend_max=3)` | `AR(β0)=m^T V^{-1} m`; numerically inverts level set; returns CI intervals, β̂, J, J p-value (df=K-1) | 34 |

**Total: 6 functions implemented, 68 testthat assertions pass (verified 2026-06-03 via `devtools::test()`).**

## 2. Stub functions

| Function | File:line | Planned signature | Blocker |
|---|---|---|---|
| `heidi_rv` | `heidi_rv.R:102` | `(b_xy_per_variant, se_per_variant, weights)` | Needs `CompQuadForm::davies()`; per-variant burden SAIGE-QTL re-run; math doc §7 not yet final |
| `annotation_concord` | `annotation_concord.R:85` | `(estimates_list, pqtl_anchor = NULL)` | Needs pQTL-anchor data (UKB-PPP / deCODE); Cochran-Q with mediator-scale normalization not yet coded |
| `iv_partial_r2` | `sensitivity.R:48` | `(b_x, se_x, n)` | Trivial to fill (Cinelli-Hazlett: `R²=t²/(t²+n-2)`, `RV=(√(t²+4)-t)/2`) — pending math-doc §10 |
| `e_value` | `sensitivity.R:99` | `(b_xy, se_xy)` | Trivial to fill (Swanson-VanderWeele: `RR=exp(0.91 β_std)`, `E=RR+√(RR(RR-1))`) — pending §10 |

**Status: 4 stubs remaining. 2 of them (`iv_partial_r2`, `e_value`) are one-liners; the other 2 (`heidi_rv`, `annotation_concord`) require either external R packages or external data and are real engineering tasks.**

## 3. Coverage matrix — main.tex steps ↔ R implementation

| main.tex step | Description | R file:line | Status |
|---|---|---|---|
| Step 0 | Two-sample summary-stat setup | n/a (notation) | – |
| Step 1 | Burden construction `Z=Σw_j G_j` | external (SAIGE-GENE+) | upstream |
| Step 2 | Per-mask `(b̂_x, SE_x)` from SAIGE-QTL | external | upstream |
| Step 3 | Per-mask `(b̂_y, SE_y)` from Genebass | external | upstream |
| Step 4 | K=1 closed-form AR CI | `mrAR.R:88` | ✅ implemented |
| Step 5 | Stack K masks into vectors | `mrAR_multi.R:153-157` | ✅ implemented |
| Step 6 | Multi-IV AR statistic `m^T V^{-1} m ~ χ²_K` | `mrAR_multi.R:162-171` | ✅ implemented |
| Step 7 | Grid + uniroot level-set inversion | `mrAR_multi.R:192-319` | ✅ implemented |
| Step 8 | Sargan-J = AR(β̂) ~ χ²_{K-1} | `mrAR_multi.R:322-341` | ✅ implemented |
| Step 9 | Sample-overlap cross-term `-2β D_y R_xy D_x` | `mrAR_multi.R:157,164` + `mrAR.R:116` | ✅ implemented |
| Step 10 | HEIDI-rv within-burden LOO | `heidi_rv.R:102` | 🔴 stub |
| Step 11 | Annotation-class Cochran-Q + pQTL anchor | `annotation_concord.R:85` | 🔴 stub |
| Step 12 | Cinelli-Hazlett partial R² / RV | `sensitivity.R:48` | 🔴 stub |
| Step 13 | Swanson-VanderWeele E-value | `sensitivity.R:99` | 🔴 stub |
| Step 14 | Cell-type concordance | (not in package yet) | 🔴 missing |

**Verdict on completeness for Phase 2/3 validation work**: Steps 4–9 (the AR causal-inference core) are fully implemented and testable. Stubs 10–14 are over-id diagnostics / sensitivity scalars — useful for production but not needed to validate AR coverage and J behavior, which is what Track 1 of `VALIDATION_PLAN.md` targets. **Phase 2/3 can proceed on `mrAR_multi` + `wald_burden` alone.**

## 4. Known correctness anchors

- **K=1 cross-check**: `test-mrAR_multi.R` line 98-115 verifies `mrAR_multi(K=1)` matches `mrAR()` closed form to `tolerance=1e-4` on a strong-IV scalar (a stricter ~1e-10 was noted in HANDOVER §4 informal smoke). The closed-form K=1 result is the rosetta stone of the package — every numerical regression in `mrAR_multi` would surface here.
- **Strong-IV bounded interval covers truth**: `test-mrAR_multi.R` line 16-33 verifies CI brackets `β_true=0.5` under K=3 and width<1.
- **J-test passes under homogeneity**: line 35-49 (K=3 pleiotropy-free, J p-value > 0.05).
- **J-test rejects under pleiotropic offset**: line 51-64 (one mask offset by 0.3, J p-value < 0.05).
- **Weak-IV K=3 yields disconnected or whole-line CI**: line 66-80 (the headline rvSMR regime).
- **Quadratic root algebra**: `test-mrAR.R` line 36-39 verifies AR at endpoint exactly equals χ²_{1,0.95} crit value.
- **Endpoints satisfy boundary def**: line 62-67 (disconnected case).

## 5. Test-coverage gaps for Phase 2/3 to fill

The 68 baseline assertions are *unit* tests on hand-picked points. They do NOT establish:

1. **95% coverage rate at α=0.05** under repeated sampling — only that a single CI contains a single true value.
2. **Type-I rate of Sargan-J** under simulated null pleiotropy (homogeneity) — only that one specific homogeneous example has p>0.05.
3. **Power of Sargan-J** under varying pleiotropy fraction.
4. **CI shape distribution** as a function of F-statistic (the rvSMR headline regime — F<10).
5. **Sample-overlap correction**: there is no test that supplying R_xy≠0 changes coverage in the right direction (only the K=1 `cor_xy` smoke `expect_false(all.equal)` test, line 48-56 of `test-wald_burden.R`).
6. **Boundary behavior**: empty CI at very small α; whole-line at very large pad; disconnected→bounded transitions as F crosses χ²_{1,0.95}≈3.84.

Phase 2/3 simulation grid (6 scenarios × 1000 reps) is designed to fill gaps 1–5.

## 6. Recommended testthat additions (not implemented in this run — flagged for Francis)

a. **Coverage regression test**: 100 reps at F_target=20, β=0.4; assert empirical 95% CI coverage ∈ [0.93, 0.97].
b. **J-test Type-I under homogeneity**: 200 reps with no pleiotropy; assert rejection rate ≤ 0.075 (one-tailed binomial).
c. **Sample-overlap correction**: with R_xy=0.3 block, assert CI shifts in the direction predicted by Step 9 sign-flip with β₀.
d. **Empty CI prevalence**: at F_target=20, no pleiotropy, β=0.4, empty CI rate should be ≤ α + 2% (Patel-Lane-Burgess 2024 §4 calibration).
e. **K=1 cross-check tightened to 1e-8**: relax tolerance back from 1e-4 → 1e-8 in `test-mrAR_multi.R` line 113-114 once grid density is tuned (current default `n_grid=8000` in that test gives ~1e-5 in practice).

## 7. Issues / risk notes carried into Phase 2/3

- **mrAR_multi rounding to `bounded_interval` when one root is at +∞**: lines 303-313 classify a single half-line as `bounded_interval` (not `disconnected_union`). For the simulation harness this is fine — it just means we should not interpret `bounded_interval` as "F > 10" without inspecting endpoints.
- **`grid_extend_max=3` ceiling**: very weak IV (F<0.5) can leak past 8× envelope; expect occasional `whole_line` classifications even when math wants `disconnected_union`. Acceptable for coverage purposes (whole_line is a valid 100% coverage instance) but worth recording in shape distribution.
- **No vectorized API**: `mrAR_multi` takes a single (gene, cell type, mask-vector) tuple per call. The Phase 3 driver must loop in R, not vectorize.

## 8. Summary verdict (for orchestrator)

Algorithm-completeness: **6/6 core (AR + Wald + utils) implemented**; **4/4 over-id/sensitivity stubs remain**. The implemented core is sufficient to validate the Track 1 simulation harness (coverage, J Type-I/power, CI shape distribution). Stubs do not block Phase 2/3. They block production use for HEIDI-rv (within-burden LOO) and annotation-class concordance, but those are over-id axes 2 and 3 — independent of the AR coverage proof.
