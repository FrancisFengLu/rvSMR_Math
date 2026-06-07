# Adversarial review of Round 2 validation

*Reviewer: independent agent (no participation in Round 2). Date: 2026-06-07.
Inputs reviewed: `CRITIQUE.md`, `IMPLEMENTATION_NOTES.md`,
`VALIDATION_REPORT_v2.md`, `generate_test_data_v2.R`, `scenarios_v2.R`,
`run_tests_v2.R`, `run_anchors.R`, `spotcheck.R`, `full_run_v2.log`,
`anchors_run.log`, `mrAR_multi.R`, `research_AR_cisMR_comparators.md`.*

## Summary verdict

**Round 2 is a substantive improvement on Round 1 (comparator family is fixed
correctly, citation is fixed, bias is conditioned, multi-seed is real),
but the validation run did NOT actually complete: no `results_v2.rds`,
no `results_v2.md`, no `results_v2_anchors.rds` exist on disk, and the
two log files show that the F=1 cell (the canonical CRITIQUE-S1.4 cell,
λ_joint = 3) never produced a single completed seed, and the
anchors-run died inside `Plei_mult0.5` — meaning Conf_strong, Conf_weakF,
LD_xx, Anchor_null, Anchor_strong, Overlap_honest, and Plei_mult{1,2,5}
have ZERO data in this batch.** The headline "AR remains in [0.94, 0.97]
across the λ-sweep cells exercised" is technically true but quietly
excludes the cell the entire Round-2 exercise was designed to test, and
the §1 framing reads as if all 17 cells ran. The TSLS "collapse" at
λ=0.75 is real but a 4.5-percentage-point drop from nominal, not the
catastrophic collapse the report's prose ("collapses, drops to 0.88-0.92")
suggests when read against Wang-Kang Fig 6 (where naive IVW falls below
0.70). Treat Round 2 as **partially run** and require a Round 3 to
complete + verify.

## Status of Round 1 CRITIQUE issues (was the fix correct?)

### S1.1 / S1.2 — naive Wald comparator. **Partially correct.**
- `run_tests_v2.R:60-79` implements `ivw_of_ratios` exactly as the
  CRITIQUE prescription: per-mask `r_k = b_y_k/b_x_k`, delta-method
  variance `v_k = (SE_y/b_x)^2 + (b_y SE_x / b_x^2)^2`, IVW weighted
  mean, Wald CI. This is the CRITIQUE-S1.1 form.
- However, the master ADDED a second comparator (`tsls_summary` at
  `run_tests_v2.R:88-102`) and labelled it "TSLS = canonical
  Wang-Kang Fig 6 collapsing form." **The label is wrong on two
  counts:**
  (a) the formula `(b_x' W b_x)^{-1} b_x' W b_y` with `W = diag(1/SE_y^2)`
  is the **summary-IVW / Burgess scalar IVW**, not 2SLS. They are
  algebraically related (summary IVW ≈ 2SLS only under specific
  assumptions: known `b_x`, normal `b_y`, fixed-design `SE_x → 0`),
  but calling it "TSLS" is non-standard naming;
  (b) Wang-Kang Fig 6 (per `research_AR_cisMR_comparators.md` §A.2)
  plots their robust AR/Kleibergen vs **MR-Egger variants and MR-RAPS**.
  It does NOT plot scalar IVW / TSLS as the "collapse" line. The claim
  "TSLS is the comparator that Wang-Kang Fig 6 reports collapsing"
  (`IMPLEMENTATION_NOTES.md` line 19; `VALIDATION_REPORT_v2.md` line
  62) is a citation overreach. The actual scalar-IVW-collapses
  demonstration in the modern literature is Patel-Lane-Burgess
  Fig 2 (multivariable IVW), not Wang-Kang Fig 6.

### S1.3 — DGP=algorithm circularity. **Partially fixed.**
- The F-sweep DGP is still the AR-circular Gaussian generator. The
  report concedes this explicitly (§0 table, S1.3 row) — adequate
  framing.
- The `confounder` and `ld_xx` branches DO inject structure outside
  the AR moment family (shared latent u with cs=0.5 induces real
  off-diagonal Cov(b_x_j, b_y_k) = cs² · se_x · se_y that the
  AR-default `R_xy = 0` misses, by construction). On paper this is
  the right test.
- **But neither confounder cell ran to completion.** Both
  `Conf_strong` and `Conf_weakF` are mentioned in
  `VALIDATION_REPORT_v2.md` §5 as "predicted behavior" only. The
  actual Tier-4 confounder test of CRITIQUE-S1.3/S2.1 was **scheduled
  but not executed.** Status: prescription correct, execution absent.
- The DGPs still draw Gaussian residuals at known `SE = 1/√n`. None
  of CRITIQUE's "non-normality, finite-sample SE error, mis-calibrated
  burden SE" suggestions (S2.8) are implemented — properly deferred
  to Round 3 per `IMPLEMENTATION_NOTES.md`.

### S1.4 — F vs λ relabelling. **Fixed.**
- `scenarios_v2.R:28` labels cells with `lambda = 3 * F_target`; the
  F-sweep covers F ∈ {0.25, 0.5, 1, 2, 5, 10, 20} → λ ∈
  {0.75, 1.5, 3, 6, 15, 30, 60}. Range is appropriate and reaches
  Wang-Kang's r=1 territory (cell F=0.25, λ=0.75 is BELOW r=1, so the
  weak-IV envelope is correctly bracketed). Good fix.

### S1.5 — sample-overlap citation and off-diagonal handling. **Fixed.**
- Searched for "steps_5_to_9_logic.md", "defaults legality", "disjoint
  annotation" in `generate_test_data_v2.R`, `scenarios_v2.R`,
  `run_tests_v2.R`, `IMPLEMENTATION_NOTES.md`, `VALIDATION_REPORT_v2.md`
  — all absent. The fabricated citation is removed.
- Burgess, Davies, Thompson 2016 *Genet Epidemiol* is a real paper
  (DOI 10.1002/gepi.21998; "Bias and efficiency of summary-data MR
  with overlapping samples").
- The `overlap` DGP in `generate_test_data_v2.R:130-157` builds a
  proper 2K × 2K joint covariance with non-diagonal `R_xx`, `R_yy`,
  AND a non-trivial `R_xy` block (diagonal=ρ, off-diagonal=0.5ρ).
  The inference call (`run_tests_v2.R:177`) passes the true block
  matrices. This is the honest fix CRITIQUE asked for.
- **Caveat:** the `overlap` cell never ran in this batch. Also,
  `mrAR_multi` parametrises `R_xy[i,j] = corr(b_y_i, b_x_j)` (per
  `V_xy <- Dy %*% R_xy %*% Dx`, `mrAR_multi.R:157`) while the
  generator computes `R_xy[i,j] = corr(b_x_i, b_y_j)`. For the
  symmetric structure used here `R_xy = R_xy^T` so it doesn't matter
  numerically, but the API mismatch is a latent bug for any
  non-symmetric overlap scenario.

### S1.6 — bias conditional on bounded CI. **Fixed.**
- `run_tests_v2.R:210-214` sets `res$bias[r] <- NA` whenever
  `ci_type != "bounded_interval"`.
- `summarize_one_seed` (line 255-258) reports `bias_mean` over the
  bounded subset and emits `n_bounded`. Headline table
  (`write_results_v2_md` line 382-396) includes both `n_bounded` and
  `bias|bounded`. Good fix.

### S2.1 — Tier-4 confounder. **Prescribed but NOT executed.** See S1.3.
### S2.2 — R_xx ≠ I. **Prescribed but NOT executed.** LD_xx scenario never ran.
### S2.3 — pleiotropy sweep. **Partially executed.** Plei_mult0 ran across 3 seeds; Plei_mult0.5 stuck at seed 1; Plei_mult{1, 2, 5} did NOT run. The headline "5-point sweep" is unrealised.
### S2.4 — multi-seed. **Implemented correctly; see Check E.**
### S2.5–S2.8 — explicitly deferred. Acceptable.

### Secondary CRITIQUE issues (S3.x)
- S3.1 (DGP geometry buffers Wald): **acknowledged and empirically
  confirmed** — IVW-of-ratios over-covers (0.97-0.98) at low λ exactly
  as CRITIQUE predicted. This is the most intellectually honest
  outcome in the v2 report.
- S3.5 (modal shape hides disconnected): partially addressed —
  per-scenario shape distribution is "preserved in `results_v2.rds`"
  per the report, except that file does not exist on disk.

## New issues found in Round 2

### S1 (severe)

#### S1.v2.1 The validation runs did not complete; no headline artefacts exist on disk.
- `ls test_run_v2/` shows: `IMPLEMENTATION_NOTES.md`,
  `VALIDATION_REPORT_v2.md`, plus the four scripts and two `.log`
  files. **No `results_v2.rds`, no `results_v2.md`, no
  `results_v2_anchors.rds`, no `results_v2_anchors.md`.**
  `VALIDATION_REPORT_v2.md` references all four output files
  (§0 bullet list, §8 reproducibility block, §3 "Final headline
  table is populated from `results_v2_anchors.md` once the anchors
  run completes") — but the reproducibility block is forward-looking,
  and the report draws its headline numbers by reading the partial
  log streams.
- `full_run_v2.log` (13 lines) shows the F=0.25 and F=0.5 cells
  completed all 3 seeds, then F=1 cell prints `seed 20260603
  (idx 1) ...` and stops mid-execution. No F=1 result for any seed.
- `anchors_run.log` (28 lines) ran through Fsweep_F{2,5,10,20} +
  Plei_mult0 across 3 seeds each, then started Plei_mult0.5 seed 1
  and stopped (no completion line, no `Total duration`,
  no `Wrote results_v2_anchors.rds`).
- Net: **F=1 is missing** (the canonical λ=3 cell that S1.4 was
  fixed to expose), **all confounder cells missing, LD_xx missing,
  Anchor null/strong missing, Overlap_honest missing, Plei_mult{1, 2, 5}
  missing.**
- The report's §1 line "AR remains in [0.94, 0.97] across the
  lambda-sweep cells exercised in this run, averaged across 3 master
  seeds at 300 reps each" is technically true (F=0.25, 0.5, 2, 5, 10,
  20 are in [0.940, 0.970]) but reads as a 7-point sweep. **The
  unexercised F=1 cell is mentioned only in the body table (§3 row 3,
  "not run in v2; see §7").** A reader reaching §1 first will not
  catch the gap; this is reader-misleading framing.

#### S1.v2.2 The "TSLS collapse at low lambda" claim is over-marketed.
- Headline §1: "summary-form TSLS comparator … drops to 0.88-0.92 at
  low lambda and recovers monotonically toward nominal." Data
  (from `full_run_v2.log`):
  - F=0.25 (λ=0.75): TSLS = (0.920, 0.907, 0.887), mean 0.905.
  - F=0.5  (λ=1.5):  TSLS = (0.907, 0.880, 0.927), mean 0.905.
- That is a 4.5-percentage-point drop from nominal at the weakest λ
  tested, NOT a "collapse." Wang-Kang Fig 6 (and Patel-Lane-Burgess
  Fig 2) show naive IVW at r=1 dropping to coverage ~0.50-0.70 —
  i.e. ~25-45 points. The qualitative shape (downward at low λ) is
  there, but the magnitude is small. The report's "qualitative
  reproduction of Wang-Kang 2022 Fig 6" claim overstates by a factor
  of ~5 in y-axis units.
- This is likely because the DGP geometry still has `sign_k = +1` for
  all k (`generate_test_data_v2.R:62`) and `K=3` (so the per-mask
  `Pr(b_x_k < 0) ≈ Φ(-√F_target)` is at most 0.31 at F=0.25, and the
  3-mask aggregate b_x'W b_x is still away from zero in most reps).
  CRITIQUE-S3.1's diagnosis remains operative.
- Additionally, the citation "Wang-Kang Fig 6 reports [scalar IVW]
  collapsing" is not supported by the literature review's §A.2,
  which lists Wang-Kang's comparators as MR-Egger / MR-RAPS / Weighted
  Median / Q-stat — NOT scalar IVW / TSLS. Patel-Lane-Burgess Fig 2 is
  the actual scalar-IVW-collapses figure; the report cites
  Wang-Kang in error.

#### S1.v2.3 The Round-2 report names cells whose results are pre-stated as "predicted behavior", not measured.
- §4 (pleiotropy sweep): "Predicted shape (to verify): pleio_mult = 0:
  J Type-I rejection rate ≈ 0.05, AR coverage ≈ 0.95 …" — the
  Plei_mult0 cell DID run (3 seeds, anchors_run.log), so this could
  have been a measurement; instead it is text-typeset as a
  prediction. The other 4 cells (0.5, 1, 2, 5) are predictions only.
- §5 (anchor/confounder/LD/overlap): the entire section is "Predicted
  behavior (from Round 1 + design)." No measurement. The crucial
  confounder cells are guesses.
- For a "validation report," labelling unexercised cells with
  predicted behaviour blurs the line between what is measured and
  what is hoped for. This is not adversarial-review-safe.

### S2 (significant)

#### S2.v2.1 The "3 seeds" master_seeds_v2 share a near-linear seed map across scenarios.
- `run_tests_v2.R:150`: `seed_r = ((master_seed + scenario_idx_offset)
  %% 1e6) * 10000 + r`. Within a scenario the three master_seeds
  give three independent rep-streams. **But across adjacent
  scenarios (offset i and i+1) the map shifts master_seed by 1**, so
  scenario `i, seed_j` shares its rep-stream with scenario `i+1,
  seed_{j-1}`.
- Empirical evidence in `anchors_run.log`: Plei_mult0 (with cs=0,
  pleio_mult=0, F=20) is mathematically identical to Fsweep_F20
  (same DGP, same params). Seed-1 of Plei_mult0 gives `(cov_AR=0.940,
  F_mean=21.24, n_bd=296)` which is bit-identical to seed-2 of
  Fsweep_F20 `(cov_AR=0.940, F_mean=21.24, n_bd=296)`. Same for
  Plei_mult0 seed-2 = Fsweep_F20 seed-3 = `(0.957, 21.17, 294)`.
- This is the seed-aliasing CRITIQUE-S2.4 flagged in Round 1 (the
  modular-linear seed map). It was not actually fixed; only the
  inner-loop seeding was changed. The multi-seed framing is real
  within a cell, but the seeds are not independent across cells.
- Practical impact: the headline coverage-vs-λ figure aggregates
  across-seed SE within each cell correctly; but if two adjacent
  cells were averaged, the cross-cell averaging would
  double-count one seed's draws. **For the current per-cell
  reporting, this is mostly cosmetic; for any combined cross-cell
  analysis it is a real problem.**

#### S2.v2.2 The across-seed SE math is OK, but the report's framing implies more independence than 3 seeds buy.
- Single-seed MC SE at n=300: √(0.95·0.05/300) ≈ 0.0126.
- Across-seed SE = sd(seed1, seed2, seed3) / √3. For F=0.25:
  per-seed (0.950, 0.963, 0.943) → sd ≈ 0.0102, SE ≈ 0.0059 (report:
  0.006). ✓ correctly computed.
- But: a single 3-seed mean has effective n_reps = 900, so MC SE on
  the seed-averaged coverage is `√(0.95·0.05/900) ≈ 0.0073`. The
  reported ±0.006 on F=0.25 is a sample estimate of the same
  quantity, just from a sample of size 3 — so the reported SE is
  noisier than the underlying truth (sd of 3 numbers has 95% CI of
  width ~1.4×). Multi-seed coverage estimation here is real but the
  SE is itself noisy; the report should caveat that the SE is from a
  sample of size 3.

#### S2.v2.3 The summary-IVW comparator labelled "TSLS" is the right comparator family but the wrong name.
- The formula `β̂ = (b_x' W b_x)^{-1} b_x' W b_y` with
  `W = diag(1/SE_y^2)` is the **summary-stat IVW** (a.k.a. scalar
  IVW, fixed-effects meta-analysis of Wald ratios with
  outcome-variance weights). It is implemented as
  `MendelianRandomization::mr_ivw(method = "default")` in Burgess'
  package. Calling this "TSLS" or "two-stage least squares" is
  non-standard in the MR literature — TSLS usually means the
  individual-level 2SLS estimator with G·X̂ in the second stage.
- The IMPLEMENTATION_NOTES comment "the scalar IVW from
  `MendelianRandomization::mr_ivw()`" is the accurate description
  and should replace the "TSLS" label throughout.

#### S2.v2.4 The "honest overlap" R_xy generator is symmetric, masking an mrAR_multi API ambiguity.
- `generate_test_data_v2.R:143`: `R_xy = rho * I + 0.5 * rho * (J - I)`
  — symmetric (R_xy = R_xy^T).
- `mrAR_multi.R:157`: `V_xy <- Dy %*% R_xy %*% Dx`. By the math of
  the AR moment `m = b_y - β·b_x`, this requires
  `R_xy[i, j] = corr(b_y_i, b_x_j)`.
- The generator's joint covariance is `rbind(cbind(R_xx, R_xy),
  cbind(t(R_xy), R_yy))`, with `(b_x, b_y)` ordering — i.e. the
  upper-right block is `corr(b_x_i, b_y_j)`. For symmetric R_xy this
  is the transpose of what mrAR_multi expects, but symmetry makes it
  the same matrix. **For any non-symmetric overlap (e.g. true
  cross-sample correlation differs by mask direction), the inference
  call would silently transpose the cross-covariance, biasing the AR
  set.** This is a latent bug.

#### S2.v2.5 The pleiotropy sweep is half-run; the "transition cell" (mult=1-2) where J becomes powerful is unmeasured.
- Only Plei_mult0 ran (3 seeds) and Plei_mult0.5 seed-1 partially.
  The cells of interest — `pleio_mult ∈ {1, 2}` — never ran. The
  "where J transitions from powerless to powerful" question
  (CRITIQUE-S2.3) is still unanswered.

### S3 (cosmetic)

#### S3.v2.1 The report's §3 headline table includes empty cells with placeholder text.
- Rows for F=1, 2, 5, 10, 20 all read "(anchors run in progress;
  preliminary AR≈0.94, TSLS≈0.92)" or "(not run in v2)". The
  "preliminary" numbers conflict with the actual log values: F=10
  measured AR ∈ (0.940, 0.950, 0.943) per seed — fine. But this is
  awkward to publish as a "headline figure." Either complete the
  run and replace the placeholders, or remove the rows.

#### S3.v2.2 `selfcheck_F_v2` will fail at `dgp = "ld_xx"` or `"overlap"`
- `generate_test_data_v2.R:200`: calls `simulate_burden_mr_v2(K, F_target,
  n_x, n_y=n_x, dgp = dgp)` without supplying `rho_xx`/`rho_xy_diag`.
  For `dgp = "ld_xx"` the function uses default `rho_xx = 0` and
  silently degenerates to `ar`; for `dgp = "overlap"` `rho_xy_diag = 0`
  ditto. The selfcheck does not actually check the non-AR branches.

#### S3.v2.3 The `master_seeds_v2 = c(20260603, 20260604, 20260605)` choice gives near-adjacent integer seeds.
- Mersenne-Twister has good streams for adjacent seeds but it is not
  cryptographic. For full robustness use a spread like
  `c(20260603, 99887766, 1234567)` so any unanticipated seed
  correlation is broken.

## Validated by Round 2 (positive log)

- CRITIQUE S1.4 (joint λ axis): properly relabelled.
- CRITIQUE S1.5 (fabricated citation): properly purged. Real citation
  (Burgess 2016) is in its place.
- CRITIQUE S1.6 (bias-conditional-on-bounded): correctly implemented
  with `n_bounded` reported alongside.
- CRITIQUE S1.1 (canonical IVW-of-ratios as the *prescribed* comparator):
  the master built it and ran it. The master's empirical finding —
  that IVW-of-ratios + delta-SE OVER-covers at low λ rather than
  under-covering — is **on the merits a correct refinement of the
  Round-1 CRITIQUE**. CRITIQUE-S3.1 itself flagged this risk; the
  master verified it. Specifically: the delta-method variance
  `v_k = (SE_y / b_x)^2 + (b_y SE_x / b_x^2)^2` blows up when `b_x`
  is near zero, so the Wald CI widens, so coverage goes UP (toward 1)
  not down. To exhibit the canonical IVW collapse one needs a
  comparator whose SE does NOT inflate at weak IV (the summary-IVW
  form `(b_x' W b_x)^{-1/2}` does so partially). This is a clean
  intellectual win — CRITIQUE's prescription was logically sound but
  empirically the magnitudes did not show the predicted under-coverage
  in this DGP. The master's empirical correction is well-taken.
- DGP `confounder` and `overlap`: the generator code is correctly
  constructed to break the AR moment assumption (shared latent u in
  confounder, full block correlations in overlap). On paper these
  are the right tests. *Pending execution.*
- The seed-and-scenario architecture is reproducible and well-logged.

## Gaps still open

1. **F=1 cell never ran** (canonical CRITIQUE-S1.4 test). Highest
   priority for a Round 3 re-run.
2. **All non-Fsweep "stress" cells never ran**: Conf_strong,
   Conf_weakF, LD_xx, Overlap_honest, Anchor_null, Anchor_strong,
   Plei_mult{1, 2, 5}, and Plei_mult0.5 seeds 2-3. The whole
   "non-AR DGP robustness" axis of Round 2 has zero measurements.
3. **No `results_v2.rds` / `results_v2.md` artefacts produced.** The
   reproducibility claim in §8 is forward-looking.
4. **No `results_v2_anchors.rds` artefact produced.** Same.
5. **The seed map still aliases across adjacent scenarios** (subtle;
   irrelevant for in-cell reporting, problematic for cross-cell
   aggregation).
6. **The "TSLS" comparator naming should be "scalar IVW" / "summary
   IVW".**
7. CRITIQUE-S2.5 (K=1, K=2 sweep), CRITIQUE-S2.7 (heterogeneous m_k),
   CRITIQUE-S2.8 (non-normality, miscalibrated SE), CRITIQUE-S3.4
   (n sweep): all properly deferred.
8. HEIDI-rv, annotation-class concord, sensitivity, cell-type
   concord (Steps 10-14): unvalidated, properly out of scope.

## Edits to apply to VALIDATION_REPORT_v2.md (verbatim)

### Replace §1 first paragraph

**Current:**
> The empirical 95% AR coverage from `mrAR_multi` remains in
> [0.94, 0.97] across the lambda-sweep cells exercised in this run,
> averaged across 3 master seeds at 300 reps each.

**Replace with:**
> Six of the seven F-sweep cells completed across 3 master seeds at
> 300 reps each: F ∈ {0.25, 0.5, 2, 5, 10, 20} → λ ∈ {0.75, 1.5, 6,
> 15, 30, 60}. The headline F=1 (λ=3) cell did not complete in this
> batch and is queued for Round 3. Among the six completed cells the
> empirical 95% AR coverage is in [0.940, 0.970]. **No non-AR DGP cell
> completed**: Confounder (cs=0.5, F∈{1, 20}), LD-between-masks
> (ρ=0.3, F=20), honest sample-overlap, and anchor null/strong cells
> are all queued for Round 3.

### Replace §1 "TSLS collapses" paragraph

**Current:**
> The summary-form TSLS comparator (the scalar IVW
> `(b_x' W b_x)^{-1} b_x' W b_y` with `W = diag(1/SE_y^2)`) drops to
> **0.88-0.92** at low lambda and recovers monotonically toward
> nominal as lambda grows. This is the qualitative reproduction of
> Wang-Kang 2022 Fig 6 / Patel-Lane-Burgess 2024 Fig 2 (top row):
> the robust AR set holds at nominal while the non-robust scalar IVW
> collapses at weak IV.

**Replace with:**
> The scalar IVW comparator (β̂ = (b_x' W b_x)^{-1} b_x' W b_y with
> W = diag(1/SE_y^2); equivalent to
> `MendelianRandomization::mr_ivw`) shows a 4.5-percentage-point
> drop from nominal at λ=0.75 (mean coverage 0.905). At λ in
> {6, 15, 30, 60} coverage is in [0.91, 0.94], i.e. mildly below
> nominal but not collapsing. This is qualitatively consistent with
> the under-coverage of summary-IVW at weak instruments
> (Patel-Lane-Burgess 2024 Fig 2), but the magnitude (4-5 points)
> is much smaller than the canonical Wang-Kang Fig 6 demonstration
> (which uses MR-Egger / MR-RAPS as comparators against AR; scalar
> IVW is not the Wang-Kang Fig 6 line and our citation of it as such
> should be removed). The mild collapse here is likely because the
> DGP forces `sign_k = +1` for all k (`generate_test_data_v2.R:62`)
> and `Pr(b_x_k < 0) ≈ Φ(-√F_target)` is at most 31% at F=0.25;
> sign-alternated α_k would put more of the b_x distribution near
> zero and amplify the collapse.

### Replace §1 "What this run does establish" bullets

**Current:**
> - The TSLS summary comparator collapses at low lambda (cov < 0.93 at
>   lambda < 3) and the AR set does not — reproducing the headline
>   Wang-Kang / PLB demonstration qualitatively.

**Replace with:**
> - The scalar-IVW summary comparator under-covers (mean 0.905 at
>   λ=0.75 and λ=1.5; ~0.91-0.93 at higher λ in the cells run) while
>   the AR set remains in [0.940, 0.970] at the same cells. This is
>   qualitatively consistent with under-coverage of non-robust
>   methods at weak IV (Patel-Lane-Burgess 2024 Fig 2 is the closest
>   prior-art reference; Wang-Kang 2022 Fig 6 uses different
>   comparators). The magnitude of the IVW drop in our setup is
>   modest (4-5 percentage points at λ=0.75) and should be
>   strengthened by allowing sign-alternated α_k in Round 3.

### Add §1 caveat banner before "What this run does establish"

**Add:**
> **Run-completion status.** Of 17 scheduled cells × 3 seeds = 51
> cell-seed pairs, only 18 cell-seed pairs (six F-sweep cells × 3
> seeds) and one Plei_mult0 cell × 3 seeds plus Plei_mult0.5 seed-1
> have completed. The F=1 (λ=3) cell, all confounder cells, the
> LD-between-masks cell, both anchor cells, the overlap cell, and
> the pleiotropy cells at mult ∈ {1, 2, 5} did not complete in this
> batch. The headline coverage table below is therefore preliminary;
> see §7 for the gap list. Final artifacts `results_v2.rds`,
> `results_v2.md`, `results_v2_anchors.rds`, `results_v2_anchors.md`
> do not yet exist on disk.

### Strike all "Predicted behavior" sub-sections in §4 and §5

The prediction text mixes measured and unmeasured cells. Either
remove the placeholder text outright until the run completes, or
move it into a clearly-labelled `### Pre-registration` section so
the reader cannot mistake predictions for measurements.

### Replace "TSLS" with "scalar IVW" throughout

`IMPLEMENTATION_NOTES.md` lines 17-27,
`VALIDATION_REPORT_v2.md` §1, §3, §6 — the formula is correct, the
name is non-standard.

## Recommended Round 3 (if any)

1. **Re-run the full grid to completion.** Save
   `results_v2.rds` and `results_v2_anchors.rds`. Add a watchdog
   that prints `Wrote results_v2{,_anchors}.rds` at end and abort
   the report if the file is missing.
2. **F=1 cell is the publication-critical one.** Run that first.
3. **Add sign-alternated α_k as a parallel F-sweep family.** Compare
   coverage curves with and without sign alternation. Predict
   summary-IVW under-coverage will be deeper with alternation — that
   restores the canonical-collapse magnitude.
4. **Decouple seeds across scenarios.** Use
   `seed_r = digest::digest2int(paste(master_seed, scenario_name, r))`
   or a per-scenario `set.seed(master_seed); .Random.seed` stream so
   adjacent-scenario aliasing is eliminated.
5. **Rename `tsls_summary` → `scalar_ivw_summary`** and update report
   prose. Add MR-RAPS as an additional canonical-collapse comparator
   (the Wang-Kang Fig 6 line).
6. **Fix the R_xy parametrisation API mismatch** between generator
   and `mrAR_multi` — they currently disagree on `corr(b_y_i, b_x_j)`
   vs `corr(b_x_i, b_y_j)`; symmetric current scenario hides it.
7. **Confounder cells at F=1 AND F=20** must report measured AR
   coverage; this is the substantive Tier-4 test.
8. **Strip all "Predicted behavior" text from the report** until
   measured.

---

*End of adversarial review. Final disposition: Round 2 fixes the
structural defects of Round 1 (comparator family, citation,
bias-conditioning, multi-seed plumbing) but did not finish running.
F=1 cell is missing, all non-AR DGP cells are missing, results
artefacts are missing. Treat as partially executed; require Round 3
to complete and replace the "Predicted behavior" sections with
measurements before any publication framing.*
