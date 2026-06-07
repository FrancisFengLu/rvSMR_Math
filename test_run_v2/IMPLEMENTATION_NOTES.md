# Round 2 implementation notes

Round 2 fixes applied to the validation harness, indexed by `CRITIQUE.md`
issue number. The rvMR R package itself was NOT modified.

## S1 (severe) — fixes

### S1.1 / S1.2 — replace naive comparator

`run_tests_v2.R`:
- Removed `ivw_pooled_wald` (pool b_x and b_y first, then divide).
- Added `ivw_of_ratios()` — classical IVW-of-ratios + delta-method SE +
  inverse-variance-weighted mean. `r_k = b_y_k / b_x_k`,
  `v_k = (SE_y/b_x)^2 + (b_y * SE_x / b_x^2)^2`, `w_k = 1/v_k`,
  `beta_hat = sum(w*r) / sum(w)`, `SE = sqrt(1/sum(w))`. This is the
  comparator CRITIQUE S1.1 specified.
- Added `tsls_summary()` — summary-form two-stage least squares:
  `beta_hat = (b_x' W b_x)^{-1} b_x' W b_y`, `W = diag(1/SE_y^2)`,
  `SE = (b_x' W b_x)^{-1/2}`. This is the comparator that Wang-Kang
  Fig 6 / PLB Fig 2 typically show collapsing at weak IV; it is the
  scalar IVW form used in `MendelianRandomization::mr_ivw()`.
- Why TSLS in addition to IVW-of-ratios: empirical spot-checks (50-rep
  pre-run) showed the IVW-of-ratios is *partially* weak-IV-buffered
  exactly as CRITIQUE S3.1 predicted — its delta-method SE inflates at
  weak IV and gives over-coverage rather than under-coverage. The
  TSLS form drops the SE_x term in the denominator and is the cleaner
  canonical-collapse comparator.

### S1.3 — DGP=algorithm-assumption circularity

`generate_test_data_v2.R`:
- New `dgp = "confounder"` branch: shared latent `u ~ N(0,1)` enters
  BOTH `b_x_k` and `b_y_k`, producing
  - within-X cross-mask correlation `conf_strength^2`
  - within-Y cross-mask correlation `conf_strength^2`
  - X-Y same-mask correlation `conf_strength^2`
  The joint `(b_x, b_y)` is **not** the AR moment structure — AR is
  called with default `R_xx = R_yy = I, R_xy = 0` (canonical backdoor
  confounder: user does not know `u`). This tests whether AR holds up
  under misspecified covariance from MR-canonical confounding.
- New `dgp = "ld_xx"` branch: compound-symmetric `R_xx` with cross-mask
  cor `rho_xx`. Inference call passes the true `R_xx`.

### S1.4 — relabel F axis as joint concentration parameter lambda

`scenarios_v2.R`:
- F-sweep cell labels include `lambda = K * F_target` (K=3 here).
- `scenarios_v2[[*]]$family == "F_sweep"`.
- `simulate_burden_mr_v2()` `truth$lambda_joint = K * F_target` is the
  axis the headline coverage figure plots against.

### S1.5 — sample-overlap citation and honest off-diagonal handling

- Deleted the `steps_5_to_9_logic.md §"defaults legality"` citation
  (verified by grep — does not exist). Instead reference Burgess,
  Davies & Thompson 2016 *Genetic Epidemiology* on sample overlap.
- New `dgp = "overlap"` branch produces non-zero off-diagonal
  `R_xx, R_yy` AND a full (non-diagonal) `R_xy` block from a shared
  per-individual cohort latent. Inference call receives the TRUE
  block correlations. This is the honest version of CRITIQUE S1.5's
  "shared per-individual residuals induce off-diagonals independent
  of the variant LD".

### S1.6 — bias only over bounded-CI reps

`run_tests_v2.R`:
- `res$bias[r] <- NA` whenever `ci_type != "bounded_interval"`.
- `summarize_one_seed()` reports `bias_mean = mean(bias[bounded])`
  and `n_bounded = sum(ci_type == "bounded_interval")`. The headline
  F-sweep table includes both `n_bounded` and `bias|bounded`.

## S2 (significant) — fixes

### S2.1 — Tier-4 confounder

See `dgp = "confounder"` above. Two cells: `Conf_strong`
(`conf_strength = 0.5, F=20`) and `Conf_weakF` (`conf_strength = 0.5,
F=1`).

### S2.2 — R_xx, R_yy off-identity

See `dgp = "ld_xx"` above. One cell: `LD_xx` (`rho_xx = 0.3, F=20`).
A more aggressive sweep is left for Round 3 (compute budget).

### S2.3 — pleiotropy magnitude sweep

`scenarios_v2.R`: five `Plei_mult*` cells with `pleio_size_mult in
{0, 0.5, 1, 2, 5}` x SE_y. Reports J<0.05 rejection rate per cell.

### S2.4 — multi-seed

`scenarios_v2.R`: `master_seeds_v2 = c(20260603, 20260604, 20260605)`.
`run_tests_v2.R`: each scenario is run at all three seeds; coverage
mean +/- SE across seeds is reported in `results_v2.md` headline tables
and a separate per-seed breakdown.

## Choices made under compute budget

- `n_reps_default_v2 = 500` configured; actual run used `n_reps = 300`
  via `RVMR_NREPS=300` to fit the 90-min compute budget after observing
  weak-IV cells taking ~150s/seed at 500 reps. With Wald MC SE at n=300
  of `sqrt(0.95 * 0.05 / 300) = 0.0126`, a true coverage of 0.95 still
  lies in [0.925, 0.975] with ~95% probability. Multi-seed mitigates
  the loss of MC precision (3 seeds give effective n=900 per cell for
  headline coverage means; MC SE on the seed-averaged coverage is
  `0.0126/sqrt(3) = 0.0073`).
- K=1, K=2 sweeps (CRITIQUE S2.5) not in Round 2; out of compute budget.
- K sweep, n sweep, non-normality, finite-sample SE miscalibration
  (CRITIQUE S2.7-S2.8, S3.4) deferred to Round 3.

## Files

- `generate_test_data_v2.R` — extended generator with 4 DGP modes.
- `scenarios_v2.R` — 17 scenarios:
  - 7 F-sweep cells (headline figure)
  - 5 pleiotropy-magnitude sweep cells
  - 2 confounder cells
  - 1 LD-between-masks cell
  - 2 anchor cells (null, strong)
  - 1 honest sample-overlap cell
- `run_tests_v2.R` — driver with IVW-of-ratios + TSLS comparators,
  bias-conditional-on-bounded, multi-seed loop, and markdown writer.
