# SUPERVISOR_REPORT.md — Round 3 paired-execution audit

Supervisor: independent reviewer (not Worker A, not Worker B). Polling window: ~36 min beginning 14:56 UTC on 2026-06-07. Final write at ~15:32 UTC.

## 1. Headline verdict

**Round 3 is a partial win**: Worker B finished cleanly with a valid (though prose-overstated) R_xy convention audit, valid mr.raps BMI->SBP sanity, and 82/82 tests passing; Worker A's sim sweep is well-designed and the seed-aliasing + TSLS-rename fixes are correct, but A's wall-time-per-scenario was misjudged so only 2 of 15 cells finished within the monitoring window (F=0.25 and F=0.5 of the homogeneous Fsweep), leaving the headline SignAlt comparator-collapse cells and the full confounder sweep incomplete at audit time.

## 2. Worker A audit (against Round-3 assigned items)

A's deliverables (all in `/home/francisfenglu4/projects/rvSMR_Math/test_run_v3a/`):
- `generate_test_data_v3a.R`, `scenarios_v3a.R`, `run_tests_v3a.R`, `full_run_v3a.log`, `seed_decoupling_check.txt`. **No** `results_v3a.md` / `.rds` / `DONE` at audit time (A still running).

| Assigned item | Verdict | Evidence |
|---|---|---|
| F=0.25, 0.5, 1 weak-IV cells (homog sign) | **PARTIAL** | F=0.25 and F=0.5 completed across all 3 seeds (6/9 runs); F=1 in flight when window closed. The 6 completed sub-runs are clean: AR coverage mean 0.956 (F=0.25) / 0.951 (F=0.5), nominal-or-slight-over; scalar IVW coverage mean 0.901 / 0.902, clear under-cover by 4-5%; IVW-of-ratios mean 0.976 / 0.971 (delta-method buffering). Pattern matches Round-2 reads. |
| Sign-alternated alpha_k SignAlt sweep (Wang-Kang style) | **NOT REACHED** | 0/18 SignAlt sub-runs completed (6 F-values x 3 seeds). The generator code (lines 67-71 of generate_test_data_v3a.R) is correct: `signs = rep(c(1,-1), length.out=K)` gives (+,-,+) for K=3, so alpha_vec = signs * sqrt(F) * SE_x correctly stresses the scalar-IVW (b_x' W b_x)^{-1} factor. AR's moment m_k = b_y - beta*b_x is sign-invariant. Code is ready; just no empirical output yet. |
| Confounder strength sweep cs in {0.1, 0.3, 0.5, 0.7, 1.0, 1.5} | **NOT REACHED** | 0/18 cs sub-runs completed. Generator code is correct; cs=1.5 will trigger `abs(cs)<=1` stopifnot and be reported as DGP-error (correctly handled). |
| TSLS rename to IVW | **DONE** | Driver script renames all `tsls_*` columns to `ivw_summary_*`; prose explicitly says "scalar IVW (Burgess MendelianRandomization::mr_ivw), NOT 2SLS". CRITIQUE_v2 §S2.v2.3 closed. |
| Seed-aliasing fix (digest-based per-rep seed) | **DONE** | `hash_seed(scenario_name, master_seed, r)` via `digest::digest2int`. `seed_decoupling_check.txt` and the run log both confirm: rep 1..5 of Fsweep_F1 vs Fsweep_F0.5 give different integer seeds (FALSE for `identical()`). The Round-2 Plei_mult0 vs Fsweep_F20 modular collision cannot recur. CRITIQUE_v2 §S2.v2.1 closed. |

**A's design choices look correct.** The wall-time-per-scenario was underestimated: ~300 sec per seed x 3 seeds x 15 cells ~ 3.75 h, vs the 36-min poll budget. This is a scoping miss, not a methodological miss; the runs already finished are clean and consistent with the intended demonstration.

## 3. Worker B audit (against Round-3 assigned items)

B's deliverables (all in `/home/francisfenglu4/projects/rvSMR_Math/test_run_v3b/`):
- `bmi_sbp_sanity.R`, `bmi_sbp_results.md`, `bmi_sbp_results.rds`, `per_snp_25.csv`, `per_snp_160.csv`, `IMPLEMENTATION_NOTES_v3b.md`, `WORKER_B_REPORT.md`, `package_diff/`, `DONE`. Plus package changes to `/home/francisfenglu4/rvSMR/May_30md/rvMR/R/mrAR_multi.R` and `tests/testthat/test-mrAR_multi.R`. Commit `6460cc1` on `main`.

| Assigned item | Verdict | Evidence |
|---|---|---|
| Fix R_xy index direction in `mrAR_multi.R` | **DONE but PROSE OVERSTATED** | B chose documentation + regression tests (no behavior change), justified by a "symmetric-part invariance" argument. My independent check (see §5 below) shows that argument is FALSE in the general case (when Dx != Dy with asymmetric R_xy, transpose changes the AR statistic; I measured AR=48.71 vs 47.44 in a hand-crafted case). However, B's actual fix is still correct: the roxygen note tells joint-cov-layout users to pass `t(R_xy)`, and the documented convention is consistent with `main.tex` §Step 9. So the package guidance is correct; only the report's blanket claim "mathematically null at the AR layer" is wrong. |
| Add asymmetric R_xy regression test | **DONE** | Three new `test_that` blocks (+14 assertions) in `test-mrAR_multi.R`. Block 1 (transpose invariance) actually only pins the special case Dx, Dy = scalar*I (the test's `se_x = rep(0.1,3)`, `se_y = rep(0.05,3)`). Block 2 (symmetric-part materiality) is solid. **Block 3** (`V_xy = Dy R_xy Dx (NOT Dx R_xy Dy)`) genuinely catches a future Dx<->Dy swap with `se_x=c(0.20,0.18,0.22)`, `se_y=c(0.05,0.07,0.04)` -- this IS the right convention-pinning test. |
| All baseline tests still pass | **DONE** | `cd /home/francisfenglu4/rvSMR/May_30md/rvMR && Rscript -e 'devtools::test()'` -> PASS 82, FAIL 0, WARN 0. Baseline was 68. |
| mr.raps BMI->SBP via `mrAR(K=1)` vs Wang-Kang 2022 Table 1 | **DONE** | Installed mr.raps v0.2 from CRAN archive. 25-SNP at p<5e-8 and 160-SNP at p<1e-4 match Wang-Kang IV sets exactly. Per-SNP mrAR + IVW meta: 25-SNP beta=0.324, CI [0.171,0.476]; 160-SNP beta=0.316, CI [0.200,0.431]. mr.raps native: 0.354 / 0.378. WK Table 1 band ~0.31-0.40 -- all consistent. CI shapes match WK weak-IV diagnostic: 23/25 bounded at strong-IV (mean F=33), 72/160 bounded + 31 disconnected + 57 whole_line at weak-IV (mean F=9). Q-test rejects homogeneity (p=3e-5 / p=0.02), expected for BMI->SBP. |

## 4. Cross-contamination check

**Clean separation.** A's driver uses `devtools::load_all("/home/francisfenglu4/rvSMR/May_30md/rvMR", quiet=TRUE)` only -- no package modifications. B's changes are confined to `R/mrAR_multi.R` (roxygen only) and `tests/testthat/test-mrAR_multi.R` (appended blocks). Neither worker wrote into the other's output directory. B's `package_diff/` is a defensive in-repo snapshot of changed package files (a good audit habit since rvMR isn't under git). No conflicts in `git status`.

## 5. Bugs uncovered in Round 3 work itself

### 5.1 B's "symmetric-part invariance" claim is mathematically wrong as stated

B writes (WORKER_B_REPORT.md §1, IMPLEMENTATION_NOTES_v3b.md §1.2):
> "AR(beta_0; R_xy) == AR(beta_0; t(R_xy)) for ANY R_xy"
> "the bug is mathematically null at the AR layer"

The correct statement: the scalar `m^T Vb^{-1} m` equals `m^T (Vb^T)^{-1} m` trivially (a scalar equals its transpose). This implies invariance under `V_xy -> t(V_xy)`. But:
- `V_xy_pkg(R_xy) = Dy R_xy Dx`
- transpose of that: `(Dy R_xy Dx)^T = Dx R_xy^T Dy`
- `V_xy_pkg(t(R_xy)) = Dy R_xy^T Dx`

`Dx R_xy^T Dy != Dy R_xy^T Dx` unless `Dx` and `Dy` are scalar multiples of `I` (or equal). Hence `R_xy -> t(R_xy)` does NOT produce a transpose of `V_xy` in general, and the AR statistic IS sensitive to the convention.

Empirical demonstration (K=3, `se_x=c(0.05,0.10,0.20)`, `se_y=c(0.10,0.15,0.30)`, off-diagonal-only asymmetric R_xy): AR=48.71 vs 47.44, a 2.7% difference. CI intervals also differ. Verified independently with `mrAR_multi`: with `se_x=c(0.05,0.10,0.20,0.25)`-style asymmetric input, `mrAR_multi(R_xy)` and `mrAR_multi(t(R_xy))` give NON-identical intervals.

### 5.2 B's first regression test only pins the special case

`test_that("R_xy index convention: AR depends only on the symmetric part of R_xy")` uses `se_x = rep(0.1, 3)` and `se_y = rep(0.05, 3)`. Both are scalar*I, so the invariance holds trivially and the test passes. This block is not pinning the general property B claims, only the rep-SE special case.

**Mitigating factor.** B's *third* regression test ("V_xy = Dy R_xy Dx (NOT Dx R_xy Dy)") DOES use unequal SEs and verifies the two conventions give different `beta_hat` and `J_stat`. So the regression suite is overall correct -- a future Dx<->Dy swap would be caught. The roxygen guidance to "pass t(R_xy)" if you have joint-cov layout is also correct.

### 5.3 Recommendation
Reword the WORKER_B_REPORT.md / IMPLEMENTATION_NOTES_v3b.md "Symmetric-part invariance" paragraph to say "holds when Dx and Dy are scalar multiples of I (e.g. when all SE_x are equal and all SE_y are equal), otherwise the convention matters". The package convention is correct as-is; only the framing is wrong.

## 6. What is now publishable

Concrete claims with evidence trail:

1. **rvMR's R_xy convention is documented and pinned.** `mrAR_multi.R` roxygen explicitly states `R_xy[i,j] = cor(b_y_i, b_x_j)`, with implementation matching `main.tex` §Step 9. Test `test-mrAR_multi.R` block 3 fires if a future contributor swaps `Dx <-> Dy` at the `V_xy` assembly. (test_run_v3b/IMPLEMENTATION_NOTES_v3b.md §1.3.2; rvMR/tests/testthat/test-mrAR_multi.R:258-301)

2. **rvMR reproduces the Wang-Kang 2022 Table 1 BMI->SBP analysis to within 1 SE on point estimates.** Per-SNP mrAR + IVW meta gives 0.324 (25-SNP) / 0.316 (160-SNP) vs mr.raps's native 0.354 / 0.378, both inside the Wang-Kang ~0.31-0.40 band. CI shape distributions reproduce the documented Wang-Kang weak-IV signature (whole_line + disconnected_union appear at p<1e-4 where median F ~3). (test_run_v3b/bmi_sbp_results.md)

3. **rvMR's mrAR_multi gives nominal-or-slight-over AR coverage at weak IV (F in {0.25, 0.5}), while scalar IVW under-covers by 4-5%.** From 6/9 completed v3a sub-runs (still ongoing): AR cov 0.956 / 0.951 vs scalar IVW 0.901 / 0.902. Pattern matches Round-2 numbers and is consistent with Wang-Kang 2022 §3. (test_run_v3a/full_run_v3a.log)

4. **Seed aliasing in the Round-2 driver is fixed.** Per-rep seed = `digest::digest2int(scenario|master_seed|rep)`. Adjacent scenarios at the same master seed give different integer seeds (verified for 5 reps in `seed_decoupling_check.txt`). The modular-linear Mersenne-Twister collision that aliased Plei_mult0 with Fsweep_F20 in Round 2 cannot recur. (test_run_v3a/run_tests_v3a.R:120-133)

5. **The `tsls_summary` comparator from Round 2 is correctly renamed `ivw_summary` and documented as scalar IVW (Burgess `MendelianRandomization::mr_ivw`), not 2SLS.** (test_run_v3a/run_tests_v3a.R:100-118, results table prose)

6. **Package tests: 82/82 pass.** Baseline 68 + B's 14 new. No regressions. (`cd rvMR && Rscript -e 'devtools::test()'`)

## 7. Still gap

1. **A did not finish.** 13 of 15 scenarios incomplete (0/18 SignAlt sub-runs, 0/18 confounder-sweep sub-runs, 3/3 F=1 sub-runs). The headline finding -- that scalar IVW COLLAPSES under sign-alternated alpha (the Wang-Kang Fig 6 demo) -- is NOT yet empirically demonstrated by Round 3; only the homogeneous-sign F=0.25/0.5 row exists, and that row already shows scalar IVW under-coverage at the milder ~10% level. The full Wang-Kang-style demonstration requires the SignAlt cells.

2. **A's `results_v3a.md` and `results_v3a.rds` files do not exist at audit time.** The driver writes them only at the end of `run_all_v3a()`. If the background process is interrupted, nothing is salvaged. (Recommend: have A checkpoint to disk after each scenario completes, not only at the end.)

3. **B's "symmetric-part invariance" prose is wrong as a general claim** (§5.1). Should be reworded. The actual code/test changes are fine.

4. **Wang-Kang Table 1 comparison is qualitative (~ "0.31-0.40 band").** B candidly notes that a formal hypothesis test of agreement would need the WK supplement. Acceptable for sanity check, not a formal replication.

5. **mrAR_multi was NOT exercised on the K=25 or K=160 BMI->SBP set.** B explicitly punted this: per-SNP intersection of bounded intervals gives empty joint set in both regimes (consistent with Q-test rejection of homogeneity). A true K-AR analog would call `mrAR_multi` with the full LD-derived `R_xx`. That is a separate validation.

6. **Round-2 generator `R_xy` convention is still in joint-cov layout** (`R_xy[i,j] = cor(b_x_i, b_y_j)`), the transpose of the package convention. Because Round-2's "overlap" mode constructs symmetric R_xy, this is silently fine for the Round-2 results. But the convention mismatch in the generator should still be fixed for clarity. B declined to touch it ("Worker A territory"). A's v3a generator (lines 100-119 of `generate_test_data_v3a.R`) builds R_xy as `diag(rho_xy_diag, K, K)` (symmetric) for the ar/signalt branches and `rho * diag(K) + 0.5*rho*(J - I)` (symmetric) for the overlap branch -- so v3a is silently fine too. Recommend an explicit comment.

## 8. Recommended Round 4 priorities

1. **Finish A's sweep.** Either (a) run unchanged at N=1000 with overnight wall time (~4 hr) and check in `results_v3a.md`/`.rds`, or (b) reduce to N=300 for a smoke pass at all 15 cells in ~70 min. Specifically, the SignAlt F=0.25/0.5/1 cells are the headline-figure cells; they MUST be done before any "Wang-Kang Fig 6 reproduction" claim can be made.

2. **Add per-scenario checkpointing to `run_tests_v3a.R`** so that interrupted runs don't lose state.

3. **Reword B's "symmetric-part invariance" prose** in WORKER_B_REPORT.md and IMPLEMENTATION_NOTES_v3b.md per §5.3 above. The package fix itself does not need to change.

4. **Add one more package regression test**: a transpose-invariance fail test with `se_x != se_y` and asymmetric R_xy, asserting that `mrAR_multi(R_xy)` and `mrAR_multi(t(R_xy))` give DIFFERENT answers (the current first test passes only in the special case). This codifies the correct general behavior.

5. **Optional: run `mrAR_multi` on the BMI->SBP 25-SNP set** with `R_xx` estimated from LD reference. That is the headline K-AR Wang-Kang replication; B's per-SNP work is a useful sanity check but not the published comparison.

---

Supervisor's progress log: `test_run_v3_super/PROGRESS.md` (transcript of polling + audit observations as they occurred).
