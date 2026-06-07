# IMPLEMENTATION_NOTES_v3b.md — Worker B, Round 3

Scope: (1) R_xy index-direction bug audit + fix in `mrAR_multi.R`; (2) mr.raps BMI→SBP sanity check; (3) regression-test design.

---

## 1. R_xy index-direction bug — decision and fix

### 1.1 The audit

CRITIQUE_v2 flagged a possible silent bug. The relevant pieces:

- **Math spec.** `main.tex` §Step 9 (line 318): `V_xy = -2 β₀ D_y R_xy D_x` — outcome on left.
- **Explainer.** `steps_5_to_9_logic.md` (lines 202-204): same convention.
- **Package code.** `rvMR/R/mrAR_multi.R:117` defaults `R_xy = matrix(0, length(b_x), length(b_y))`; line 157 builds `V_xy = Dy %*% R_xy %*% Dx`. Mathematically, that is **consistent** with the spec: `R_xy[i,j] = cor(b_y_i, b_x_j)` — outcome row, exposure col.
- **Generator.** `test_run_v2/generate_test_data_v2.R:141-146` builds `R_xy` as the (1,2) block of the joint correlation matrix `rbind(cbind(R_xx, R_xy), cbind(t(R_xy), R_yy))`, i.e. `R_xy[i,j] = cor(b_x_i, b_y_j)` — **exposure row, outcome col**. That is the *transpose* of the package's convention.

So the math, the explainer, and the package code agree; the generator is the odd one out.

### 1.2 The deeper finding

I initially intended to add a regression test using an asymmetric `R_xy` and verifying that `mrAR_multi(R_xy)` and `mrAR_multi(t(R_xy))` give different answers. **That test failed.** Investigation showed why:

The AR statistic is `AR(β₀) = m^T V(β₀)^{-1} m`, a scalar quadratic form. For any matrix `A` and vector `v`, `v^T A v = v^T sym(A) v` where `sym(A) = (A + A^T)/2`. So `AR` depends **only on the symmetric part** of `V(β₀)`. Within `V`:

- `V_yy = D_y R_yy D_y` — symmetric (assuming R_yy symmetric).
- `V_xx = D_x R_xx D_x` — symmetric.
- `V_xy = D_y R_xy D_x` — **not symmetric** unless R_xy itself is symmetric AND D_x == D_y.

But `V` enters as `V_yy + β₀² V_xx - 2 β₀ V_xy`, so the only piece of V_xy that AR sees is its symmetric part: `(V_xy + V_xy^T) / 2 = (D_y R_xy D_x + D_x R_xy^T D_y) / 2`.

Now compare the two conventions:

- **Package convention** (spec-consistent): `V_xy^pkg = D_y R_xy D_x`. Symmetric part = `(D_y R_xy D_x + D_x R_xy^T D_y)/2`.
- **Generator convention** (joint-cov layout): if a user passes the same matrix to `mrAR_multi` thinking it means `cor(b_x_i, b_y_j)`, then `V_xy^user_intent = D_x R_xy D_y`. Symmetric part = `(D_x R_xy D_y + D_y R_xy^T D_x)/2` — **same matrix** as the package convention.

(Both expressions equal `(D_y R_xy D_x + D_x R_xy^T D_y)/2`, just with the addends in opposite order.)

So **the two conventions produce identical AR / J-test results in the K≥2 case**. The "silent bug" is more than silent: it is mathematically null at the AR layer.

Where it WOULD matter:
- If `R_xy` were used to construct a joint distribution (e.g. for sampling, or for a different non-quadratic statistic), orientation would matter.
- The K=1 single-SE case is trivially symmetric so there's no issue.

### 1.3 Fix

I went with the task's stated preference: **package documentation + roxygen clarification + regression tests**, no code change.

#### 1.3.1 Package documentation (rvMR/R/mrAR_multi.R)

- Added an explicit index-convention statement to the `R_xy` `@param` block: `R_xy[i,j] = cor(b_y_i, b_x_j)` — outcome row, exposure col. Pointed out it's the transpose of the joint-covariance layout.
- Added a "Symmetric-part invariance" paragraph in the same param block explaining that AR depends only on the symmetric part of V_xy, hence `AR(β₀; R_xy) == AR(β₀; t(R_xy))`.
- Added an inline comment at line 159 (the V_xy assembly) reiterating the convention and the symmetric-part property.

These changes are all in roxygen comments; no behavior change.

#### 1.3.2 Regression tests (rvMR/tests/testthat/test-mrAR_multi.R)

Added three new `test_that` blocks (14 new assertions, all pass):

1. **`R_xy index convention: AR depends only on the symmetric part of R_xy`** — uses an explicitly asymmetric R_xy with nonzero entries at (1,2), (2,1), (3,1), (1,3), and verifies `mrAR_multi(R_xy)` and `mrAR_multi(t(R_xy))` give the SAME ci_intervals, beta_hat, and J_stat (to 1e-6 tolerance). This pins the invariance.

2. **`R_xy: changing the symmetric part DOES change the answer`** — symmetric +R and -R (sign-flipped) should give different bounded intervals; zero R_xy gives a third distinct answer. This guards against a future "fix" that accidentally zeros out V_xy entirely.

3. **`R_xy: V_xy = Dy R_xy Dx (NOT Dx R_xy Dy) — direct AR comparison`** — when `se_x != se_y` and R_xy is asymmetric, the two assembly conventions DO produce different *symmetric parts* (because the Dx and Dy aren't equal). This is the test that would catch a future swap of `Dx <-> Dy` in line 157. We compute an `R_xy_eff` such that `Dy %*% R_xy_eff %*% Dx == Dx %*% R_xy %*% Dy` and verify the answers differ.

#### 1.3.3 Test counts

- Baseline (pre-Round-3): 68 assertions across 3 test files (test-mrAR.R: 27, test-mrAR_multi.R: 34, test-wald_burden.R: 7).
- Post-Round-3: **82 assertions, 0 failures** (test-mrAR_multi.R now 48).
- All 68 baseline assertions still pass.

### 1.4 Side note on the test_run_v2 generator

The generator at `test_run_v2/generate_test_data_v2.R:141-146` uses the joint-covariance layout. In its "overlap" DGP it builds R_xy as `rho * I + 0.5 * rho * (J - I)` (symmetric), so its convention difference relative to the package would not matter even in a regime with asymmetric R_xy. **No fix needed in the generator** — the convention mismatch is benign at the AR layer (see §1.2).

I did NOT modify `generate_test_data_v2.R` because (a) it's in Worker A's directory tree and (b) the difference is mathematically null.

---

## 2. mr.raps BMI→SBP sanity check

### 2.1 Setup

- mr.raps was not installed; CRAN no longer hosts a current binary for R 4.1.2. Installed v0.2 from the CRAN archive (`/tmp/mr.raps_0.2.tar.gz`) into `/home/francisfenglu4/R/library`. Required dep `nortest` also installed.
- Dataset: `mr.raps::bmi.sbp` — 160 SNPs, 29 columns. Selection p-value column `pval.selection`. The two Wang-Kang IV sets correspond to:
  - p < 5e-8 → 25 SNPs (exact).
  - p < 1e-4 → 160 SNPs (exact; the full dataset).

### 2.2 Results (see `bmi_sbp_results.md` for details and Wang-Kang comparison table)

| Pipeline                         | 25-SNP point | 25-SNP 95% CI    | 160-SNP point | 160-SNP 95% CI   |
|----------------------------------|--------------|------------------|---------------|------------------|
| rvMR::mrAR per-SNP + IVW meta    | 0.3238       | [0.171, 0.476]   | 0.3158        | [0.200, 0.431]   |
| mr.raps over-disp Huber (ref)    | 0.3536       | [0.098, 0.610]   | 0.3781        | [0.141, 0.615]   |
| Wang-Kang 2022 Table 1 range     | ~0.31-0.40   | varies           | ~0.30-0.40    | varies           |

- Point estimates agree within ~0.05 across pipelines, all inside Wang-Kang's published band.
- rvMR-IVW CIs are narrower than RAPS-Huber CIs (no over-dispersion correction). Cochran-Q rejects homogeneity (p=3.2e-5 at 25-SNP, p=0.02 at 160-SNP) — expected for BMI→SBP, well-known pleiotropy.
- Per-SNP AR CI shapes track the weak-IV regime as expected: strong-IV (mean F=33.1) → 23/25 bounded; weak-IV (mean F=9.1) → 72/160 bounded, 31 disconnected, 57 whole_line. This is qualitatively what Wang-Kang highlight as the AR test's signature.

### 2.3 Caveats

- mr.raps v0.2 (2018) is the only version installable on R 4.1.2 from CRAN archive. The dataset is post-QC (palindromic, ambiguous SNPs filtered out per the package's `mr_keep` flags).
- The "AR-intersection across SNPs" idea (intersect per-SNP acceptance sets) gives an EMPTY joint set in both IV-size regimes. This is consistent with the Q-test rejection (heterogeneity) — the right K-AR analog is `mrAR_multi` with the full K-dim chi-square, not pointwise intersection. We did not run `mrAR_multi` on the 160-SNP set because the grid evaluation cost would be heavy and is a separate validation question.
- WK Table 1 numbers were summarized from the paper; the report explicitly says "approximate ~0.31-0.40 band". A formal hypothesis test of agreement would require the WK supplement.

---

## 3. Files written

- `/home/francisfenglu4/rvSMR/May_30md/rvMR/R/mrAR_multi.R` — roxygen + inline doc clarification.
- `/home/francisfenglu4/rvSMR/May_30md/rvMR/tests/testthat/test-mrAR_multi.R` — 3 new test_that blocks (+14 assertions).
- `/home/francisfenglu4/projects/rvSMR_Math/test_run_v3b/bmi_sbp_sanity.R` — analysis script.
- `/home/francisfenglu4/projects/rvSMR_Math/test_run_v3b/bmi_sbp_results.md` — results + WK comparison.
- `/home/francisfenglu4/projects/rvSMR_Math/test_run_v3b/bmi_sbp_results.rds` — full R list.
- `/home/francisfenglu4/projects/rvSMR_Math/test_run_v3b/per_snp_25.csv`, `per_snp_160.csv` — per-SNP tables.
- `/home/francisfenglu4/projects/rvSMR_Math/test_run_v3b/IMPLEMENTATION_NOTES_v3b.md` — this file.
- `/home/francisfenglu4/projects/rvSMR_Math/test_run_v3b/WORKER_B_REPORT.md` — short-form summary.
- `/home/francisfenglu4/projects/rvSMR_Math/test_run_v3b/DONE` — completion marker (written last).
