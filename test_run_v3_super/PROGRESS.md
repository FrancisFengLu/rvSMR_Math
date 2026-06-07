# Supervisor progress log

Round 3 paired-execution validation cycle. Two workers expected:
- Worker A: sim sweep in `test_run_v3a/`
- Worker B: package fix + real data in `test_run_v3b/`

Each writes a `DONE` file when finished. Polling every 3 min, max 36 min.

## Attempt 0 - initialization
- `test_run_v3a/` exists, empty (created Jun 7 14:55)
- `test_run_v3b/` does NOT yet exist
- Will create polling loop and audit periodically

--- attempt 1 at 14:56:52Z ---
A dir:
  generate_test_data_v3a.R
B dir:
  ls: cannot access '/home/francisfenglu4/projects/rvSMR_Math/test_run_v3b/': No such file or directory
git log:
  ac92fe5 Round 2: complete anchors+pleio+special-cells run; populate VALIDATION_REPORT_v2 with final numbers
  5db73bc Add Round 2 adversarial review
  740be2d Round 2 rvMR validation: fix CRITIQUE issues, replace naive comparator with TSLS, add non-AR DGPs, multi-seed

--- attempt 2 at 14:59:52Z ---
A dir:
  generate_test_data_v3a.R
  results_v3a.md
  results_v3a.rds
  run_tests_v3a.R
  scenarios_v3a.R
B dir:
  ls: cannot access '/home/francisfenglu4/projects/rvSMR_Math/test_run_v3b/': No such file or directory
git log:
  ac92fe5 Round 2: complete anchors+pleio+special-cells run; populate VALIDATION_REPORT_v2 with final numbers
  5db73bc Add Round 2 adversarial review
  740be2d Round 2 rvMR validation: fix CRITIQUE issues, replace naive comparator with TSLS, add non-AR DGPs, multi-seed

--- attempt 3 at 15:02:52Z ---
A dir:
  full_run_v3a.log
  generate_test_data_v3a.R
  run_tests_v3a.R
  scenarios_v3a.R
B dir:
git log:
  ac92fe5 Round 2: complete anchors+pleio+special-cells run; populate VALIDATION_REPORT_v2 with final numbers
  5db73bc Add Round 2 adversarial review
  740be2d Round 2 rvMR validation: fix CRITIQUE issues, replace naive comparator with TSLS, add non-AR DGPs, multi-seed

## ~5 min audit (A) + B dir now appears
- A: full_run_v3a.log shows F=0.25 (lambda=0.75) seed 20260603 currently running; smoke-test outputs from 14:59 deleted, full N=1000 run launched ~15:01
- A generator: `signalt` branch (lines 67-71 of generate_test_data_v3a.R) sets signs = rep(c(1,-1), length.out=K). With K=3 -> (+1,-1,+1). alpha_vec = signs * sqrt(F) * SE_x. AR residual covariance identical to `ar` mode. This correctly stresses scalar IVW: sum(b_x*W*b_y) -> 0 while sum(b_x*W*b_x) -> 0 at low F. AR moment m_k = b_y - beta*b_x is sign-invariant -> AR coverage should hold.
- A driver: hash_seed(scenario_name, master_seed, r) via digest::digest2int -> non-bit-identical seeds across adjacent scenarios (verified in log).
- A scenarios: 15 cells = 3 Fsweep_homog (F=0.25,0.5,1) + 6 SignAlt (F=0.25,0.5,1,2,5,20) + 6 ConfSweep (cs=0.1,0.3,0.5,0.7,1.0,1.5)
- A code does NOT touch rvMR package; load_all only — clean separation.
- B: dir exists at 15:02, still empty. No files yet.

--- attempt 4 at 15:05:52Z ---
A dir:
  full_run_v3a.log
  generate_test_data_v3a.R
  run_tests_v3a.R
  scenarios_v3a.R
  seed_decoupling_check.txt
B dir:
  IMPLEMENTATION_NOTES_v3b.md
  WORKER_B_REPORT.md
  bmi_sbp_results.md
  bmi_sbp_results.rds
  bmi_sbp_sanity.R
  per_snp_160.csv
  per_snp_25.csv
git log:
  ac92fe5 Round 2: complete anchors+pleio+special-cells run; populate VALIDATION_REPORT_v2 with final numbers
  5db73bc Add Round 2 adversarial review
  740be2d Round 2 rvMR validation: fix CRITIQUE issues, replace naive comparator with TSLS, add non-AR DGPs, multi-seed

--- attempt 5 at 15:08:52Z ---
A dir:
  full_run_v3a.log
  generate_test_data_v3a.R
  run_tests_v3a.R
  scenarios_v3a.R
  seed_decoupling_check.txt
B dir:
  DONE
  IMPLEMENTATION_NOTES_v3b.md
  WORKER_B_REPORT.md
  bmi_sbp_results.md
  bmi_sbp_results.rds
  bmi_sbp_sanity.R
  package_diff
  per_snp_160.csv
  per_snp_25.csv
git log:
  6460cc1 Round 3 Worker B: audit R_xy index convention in mrAR_multi, run mr.raps BMI-SBP sanity
  ac92fe5 Round 2: complete anchors+pleio+special-cells run; populate VALIDATION_REPORT_v2 with final numbers
  5db73bc Add Round 2 adversarial review

## ~10 min: B DONE; my independent verification of B's central math claim

### Tests pass
- `Rscript -e 'devtools::test()'` -> PASS 82, FAIL 0, WARN 0. Baseline 68 + B's new 14 assertions in 3 test_that blocks.
- B committed at 6460cc1 (no push needed by supervisor since worker pushed). Package R/mrAR_multi.R + test-mrAR_multi.R modified — within B's scope.

### CRITICAL: B's "symmetric-part invariance" claim is OVERSTATED
B writes "AR(beta; R_xy) == AR(beta; t(R_xy)) for ANY R_xy" and concludes "the bug is mathematically null at the AR layer".
This is FALSE in the general case. My independent check:
- The identity m^T Vb^{-1} m == m^T Vb^{-T} m holds trivially (a scalar equals its transpose).
- This implies AR is invariant to V_xy -> t(V_xy).
- BUT R_xy -> t(R_xy) does NOT correspond to V_xy -> t(V_xy) when Dx != Dy:
    V_xy(R_xy)^T = Dx R_xy^T Dy   while   V_xy(t(R_xy)) = Dy R_xy^T Dx.
- These match only when Dx and Dy are scalar multiples of I (or equal).
- Empirical: with se_x=c(0.05,0.10,0.20), se_y=c(0.10,0.15,0.30), asymmetric R_xy, AR = 48.71 (orig) vs 47.44 (transposed) -- ~2.7% difference. CI intervals also differ.

### Why B's first regression test passes anyway
B's "transpose invariance" test uses se_x = rep(0.1, 3) and se_y = rep(0.05, 3) — both scalar multiples of I. In that special case the symmetric parts ARE equal, so invariance holds. But it's an artifact of the test inputs, not a general property.

### B's third regression test actually catches the real issue
"R_xy: V_xy = Dy R_xy Dx (NOT Dx R_xy Dy)" uses se_x=c(0.20,0.18,0.22), se_y=c(0.05,0.07,0.04) — unequal SEs — and verifies the two conventions give different beta_hat and J_stat. So B did pin the convention against a future swap.

### Bug verdict — corrected
- The package's convention IS the documented one (V_xy = Dy R_xy Dx, i.e. R_xy[i,j] = cor(b_y_i, b_x_j)).
- The Round-2 generator uses the joint-cov layout (R_xy[i,j] = cor(b_x_i, b_y_j)) — the TRANSPOSE.
- For symmetric R_xy (which Round-2 generator's "overlap" mode produces), the difference vanishes regardless of SE heterogeneity.
- For asymmetric R_xy with unequal SEs, the conventions ARE different and give different AR values. If a real-data user supplies an asymmetric R_xy estimated in joint-cov layout, they get a wrong answer (few % off).
- B's roxygen note saying "pass t(R_xy)" if you have joint-cov layout is correct guidance and prevents this.

### Real-data sanity (BMI->SBP) — looks credible
- 25-SNP IVW point 0.324 vs mr.raps 0.354 vs Wang-Kang band ~0.31-0.40 — all consistent.
- CI shape: 23/25 bounded at strong-IV (mean F=33), 72/160 bounded + 31 disconnected + 57 whole_line at weak-IV (mean F=9) — matches Wang-Kang's signature.
- Q-test rejects homogeneity (p=3e-5 / p=0.02) — known pleiotropy. CI narrower than RAPS because no over-dispersion correction; B notes this honestly.
- Caveat: WK Table 1 band is "approximate" (B says so); not a formal hypothesis test.

- A still running.

## ~15 min: A is slow but on track

A's log shows ~300 sec per seed at F=0.25, 2 of 3 seeds done at 15:11. With 15 scenarios x 3 seeds x ~300s estimate = ~3.75 hours total. **A will not finish within the 36-min poll budget.** Per brief, NOT cancelling A; instead will let polling expire naturally and write SUPERVISOR_REPORT on whatever A has produced.

### Preliminary signal from F=0.25 (only completed cell so far)
- AR coverage: 0.964 / 0.955 (slight over-cover; nominal 0.95).
- IVW-of-ratios coverage: 0.982 / 0.971 (over-cover from delta-method buffering at low F).
- **Scalar IVW (ivw_summary) coverage: 0.911 / 0.903 -- UNDER-cover by 4-5% at lambda=0.75.** This is the expected Wang-Kang collapse pattern.
- n_bounded: ~95 / 1000 -- meaning ~90% of AR CIs at F=0.25 are NOT bounded intervals. Canonical AR weak-IV regime.
- This single cell is already on-message; SignAlt cells should crush IVW even harder.

--- attempt 6 at 15:11:52Z ---
A dir:
  full_run_v3a.log
  generate_test_data_v3a.R
  run_tests_v3a.R
  scenarios_v3a.R
  seed_decoupling_check.txt
B dir:
  DONE
  IMPLEMENTATION_NOTES_v3b.md
  WORKER_B_REPORT.md
  bmi_sbp_results.md
  bmi_sbp_results.rds
  bmi_sbp_sanity.R
  package_diff
  per_snp_160.csv
  per_snp_25.csv
git log:
  6460cc1 Round 3 Worker B: audit R_xy index convention in mrAR_multi, run mr.raps BMI-SBP sanity
  ac92fe5 Round 2: complete anchors+pleio+special-cells run; populate VALIDATION_REPORT_v2 with final numbers
  5db73bc Add Round 2 adversarial review

--- attempt 7 at 15:14:52Z ---
A dir:
  full_run_v3a.log
  generate_test_data_v3a.R
  run_tests_v3a.R
  scenarios_v3a.R
  seed_decoupling_check.txt
B dir:
  DONE
  IMPLEMENTATION_NOTES_v3b.md
  WORKER_B_REPORT.md
  bmi_sbp_results.md
  bmi_sbp_results.rds
  bmi_sbp_sanity.R
  package_diff
  per_snp_160.csv
  per_snp_25.csv
git log:
  6460cc1 Round 3 Worker B: audit R_xy index convention in mrAR_multi, run mr.raps BMI-SBP sanity
  ac92fe5 Round 2: complete anchors+pleio+special-cells run; populate VALIDATION_REPORT_v2 with final numbers
  5db73bc Add Round 2 adversarial review

--- attempt 8 at 15:17:52Z ---
A dir:
  full_run_v3a.log
  generate_test_data_v3a.R
  run_tests_v3a.R
  scenarios_v3a.R
  seed_decoupling_check.txt
B dir:
  DONE
  IMPLEMENTATION_NOTES_v3b.md
  WORKER_B_REPORT.md
  bmi_sbp_results.md
  bmi_sbp_results.rds
  bmi_sbp_sanity.R
  package_diff
  per_snp_160.csv
  per_snp_25.csv
git log:
  6460cc1 Round 3 Worker B: audit R_xy index convention in mrAR_multi, run mr.raps BMI-SBP sanity
  ac92fe5 Round 2: complete anchors+pleio+special-cells run; populate VALIDATION_REPORT_v2 with final numbers
  5db73bc Add Round 2 adversarial review

--- attempt 9 at 15:20:52Z ---
A dir:
  full_run_v3a.log
  generate_test_data_v3a.R
  run_tests_v3a.R
  scenarios_v3a.R
  seed_decoupling_check.txt
B dir:
  DONE
  IMPLEMENTATION_NOTES_v3b.md
  WORKER_B_REPORT.md
  bmi_sbp_results.md
  bmi_sbp_results.rds
  bmi_sbp_sanity.R
  package_diff
  per_snp_160.csv
  per_snp_25.csv
git log:
  6460cc1 Round 3 Worker B: audit R_xy index convention in mrAR_multi, run mr.raps BMI-SBP sanity
  ac92fe5 Round 2: complete anchors+pleio+special-cells run; populate VALIDATION_REPORT_v2 with final numbers
  5db73bc Add Round 2 adversarial review

--- attempt 10 at 15:23:52Z ---
A dir:
  full_run_v3a.log
  generate_test_data_v3a.R
  run_tests_v3a.R
  scenarios_v3a.R
  seed_decoupling_check.txt
B dir:
  DONE
  IMPLEMENTATION_NOTES_v3b.md
  WORKER_B_REPORT.md
  bmi_sbp_results.md
  bmi_sbp_results.rds
  bmi_sbp_sanity.R
  package_diff
  per_snp_160.csv
  per_snp_25.csv
git log:
  6460cc1 Round 3 Worker B: audit R_xy index convention in mrAR_multi, run mr.raps BMI-SBP sanity
  ac92fe5 Round 2: complete anchors+pleio+special-cells run; populate VALIDATION_REPORT_v2 with final numbers
  5db73bc Add Round 2 adversarial review

--- attempt 11 at 15:26:52Z ---
A dir:
  full_run_v3a.log
  generate_test_data_v3a.R
  run_tests_v3a.R
  scenarios_v3a.R
  seed_decoupling_check.txt
B dir:
  DONE
  IMPLEMENTATION_NOTES_v3b.md
  WORKER_B_REPORT.md
  bmi_sbp_results.md
  bmi_sbp_results.rds
  bmi_sbp_sanity.R
  package_diff
  per_snp_160.csv
  per_snp_25.csv
git log:
  6460cc1 Round 3 Worker B: audit R_xy index convention in mrAR_multi, run mr.raps BMI-SBP sanity
  ac92fe5 Round 2: complete anchors+pleio+special-cells run; populate VALIDATION_REPORT_v2 with final numbers
  5db73bc Add Round 2 adversarial review

--- attempt 12 at 15:29:52Z ---
A dir:
  full_run_v3a.log
  generate_test_data_v3a.R
  run_tests_v3a.R
  scenarios_v3a.R
  seed_decoupling_check.txt
B dir:
  DONE
  IMPLEMENTATION_NOTES_v3b.md
  WORKER_B_REPORT.md
  bmi_sbp_results.md
  bmi_sbp_results.rds
  bmi_sbp_sanity.R
  package_diff
  per_snp_160.csv
  per_snp_25.csv
git log:
  6460cc1 Round 3 Worker B: audit R_xy index convention in mrAR_multi, run mr.raps BMI-SBP sanity
  ac92fe5 Round 2: complete anchors+pleio+special-cells run; populate VALIDATION_REPORT_v2 with final numbers
  5db73bc Add Round 2 adversarial review
POLL EXITED at 15:32:52Z
