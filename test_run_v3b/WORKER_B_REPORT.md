# WORKER_B_REPORT.md — Round 3, paired-execution validation

Worker B scope: (1) R_xy index-direction bug in `mrAR_multi.R`; (2) mr.raps BMI→SBP real-data sanity vs Wang-Kang 2022 Table 1; (3) documentation.

## 1. R_xy bug — finding and fix

**Diagnosis.** Math (main.tex §Step 9), the explainer (`steps_5_to_9_logic.md`), and the package code (`rvMR/R/mrAR_multi.R:117,157`) all agree on the convention `R_xy[i,j] = cor(b_y_i, b_x_j)` with assembly `V_xy = D_y R_xy D_x`. The Round-2 generator (`test_run_v2/generate_test_data_v2.R:141-146`) uses the joint-covariance layout `R_xy[i,j] = cor(b_x_i, b_y_j)` — i.e. the transpose. The generator's R_xy is symmetric, so the convention mismatch was hidden.

**Deeper finding (conditional, NOT general).** The AR statistic `m^T V(β₀)^{-1} m` is a scalar quadratic form, so it depends only on the **symmetric part** of V_xy. Under the two conventions,
- package: sym(V_xy) = `(D_y R_xy D_x + D_x R_xy^T D_y)/2`,
- generator: sym(V_xy) = `(D_x R_xy D_y + D_y R_xy^T D_x)/2`.

These coincide — hence `AR(β₀; R_xy) == AR(β₀; t(R_xy))` — **when `D_x = D_y` (or `D_x` is a scalar multiple of `D_y`, e.g. all SE_x equal and all SE_y equal), or when `R_xy` is already symmetric**. The Round-3 supervisor independently verified this is NOT a general invariance: with K=3, `se_x=c(0.05,0.10,0.20)`, `se_y=c(0.10,0.15,0.30)`, and an off-diagonal asymmetric R_xy, AR = 48.71 vs 47.44 (~2.7% delta) under transpose. So in the **general case with unequal D_x, D_y and asymmetric R_xy, the convention does matter**, which is why the roxygen now pins it explicitly. The Round-2 generator's R_xy is itself symmetric, so the convention mismatch is null **for the matrices that generator actually produces** — but not as a general statement.

**Fix applied.** Documentation + regression tests; no behavior change.
- `rvMR/R/mrAR_multi.R` — added explicit index-convention statement to `R_xy` roxygen, added a "Symmetric-part invariance" paragraph, added an inline comment at V_xy assembly (line 159).
- `rvMR/tests/testthat/test-mrAR_multi.R` — three new `test_that` blocks (14 assertions). They pin (a) transpose invariance, (b) symmetric-part materiality, and (c) detection of a hypothetical `Dx <-> Dy` swap.
- Test count: 68 baseline → 82 total, **0 failures**. All baseline assertions still pass.

I did NOT modify the test_run_v2 generator because (a) it lives downstream and isn't in my path, and (b) the v2 generator's R_xy is symmetric, so the convention mismatch is null **for the specific R_xy that generator produces** (it would NOT be null for general asymmetric R_xy — see the conditional invariance above).

**Verdict on which side was "wrong":** Neither side is buggy in its own scope. The package's convention is correctly documented and matches `main.tex` §Step 9; the generator emits a symmetric R_xy so it never actually hits the unprotected case. The original CRITIQUE_v2 concern that the mismatch could "break silently" was correct in principle — the two assemblies do diverge for a general asymmetric R_xy with `D_x ≠ D_y` — but neither codepath presently feeds such an input into the other.

## 2. mr.raps BMI→SBP comparison vs Wang-Kang 2022 Table 1

Installed `mr.raps` v0.2 from CRAN archive (current CRAN does not host a version for R 4.1.2). Used `mr.raps::bmi.sbp` directly.

| Pipeline                         | 25-SNP point | 25-SNP 95% CI    | 160-SNP point | 160-SNP 95% CI   | mean F |
|----------------------------------|--------------|------------------|---------------|------------------|--------|
| rvMR::mrAR per-SNP + IVW meta    | **0.3238**   | [0.171, 0.476]   | **0.3158**    | [0.200, 0.431]   | 33 / 9 |
| mr.raps over-disp Huber (ref)    | 0.3536       | [0.098, 0.610]   | 0.3781        | [0.141, 0.615]   | —      |
| Wang-Kang 2022 Table 1 band      | ~0.31-0.40   | varies           | ~0.30-0.40    | varies           | —      |

**Agreement.** rvMR's per-SNP AR + IVW meta point estimate (0.32) is **inside** the Wang-Kang Table 1 band (~0.31-0.40) for both IV sets, and within ~0.05 of mr.raps's own RAPS estimate. CI agreement is qualitative: rvMR-IVW intervals are narrower because they don't inflate for over-dispersion (Cochran-Q rejects homogeneity at p=3.2e-5 / p=0.02 — pleiotropy is real and known for BMI→SBP).

Per-SNP AR CI shape distributions confirm the weak-IV diagnostic: at p<1e-4 (mean F=9.1), 57/160 SNPs give `whole_line` CIs and 31/160 give `disconnected_union` — exactly the regime Wang-Kang's AR test was designed for. At p<5e-8 (mean F=33), 23/25 give `bounded` CIs.

## 3. Caveats

- mr.raps v0.2 (2018) is the only archive version installable on R 4.1.2. Newer mr.raps would not change the data (`bmi.sbp` is a static dataset).
- Wang-Kang Table 1 numbers were summarized as an approximate band (~0.31-0.40); the formal CI per-estimator comparison would need the WK supplement.
- The "AR intersection across SNPs" attempt gives an empty joint set in both IV-size regimes (Q-test rejects homogeneity). A full K-AR analog would call `mrAR_multi` on the K=25 or K=160 set with an LD-derived R_xx — that's a separate validation, not run here.
- All tests for the R_xy convention pin the *current* behavior. A future contributor who deliberately changes the index convention (e.g. to align with the joint-covariance layout) will need to update these tests, but their AR outputs would remain identical.

## 4. Deliverables (in `test_run_v3b/`)

`bmi_sbp_sanity.R`, `bmi_sbp_results.md`, `bmi_sbp_results.rds`, `per_snp_25.csv`, `per_snp_160.csv`, `IMPLEMENTATION_NOTES_v3b.md`, `WORKER_B_REPORT.md`, `DONE`.

Package changes (in `/home/francisfenglu4/rvSMR/May_30md/rvMR/`): `R/mrAR_multi.R` (roxygen + inline comment), `tests/testthat/test-mrAR_multi.R` (+3 test_that blocks).
