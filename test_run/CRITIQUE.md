# Adversarial review of Phase 1-3 validation

*Reviewer: Claude (adversarial mode). Date: 2026-06-06.
Inputs reviewed: `VALIDATION_REPORT.md`, `audit.md`, `generate_test_data.R`,
`scenarios.R`, `run_tests.R`, `results.md`, `mrAR_multi.R`, `mrAR.R`,
`main.tex` Steps 4-9, `steps_5_to_9_logic.md`, `VALIDATION_PLAN.md`,
`research_AR_cisMR_comparators.md`.*

## Summary verdict

**The "6/6 PASS" headline is not publishable as a weak-instrument robustness
demonstration.** Three of the six PASS cells are tautological (the data
were generated from exactly the model the AR statistic assumes), the naive
Wald comparator is implemented in a form that *cannot* exhibit the
weak-IV failure mode the report claims it exhibits, and the report
narrates over-coverage as "naive Wald collapses" — the opposite of what
the result actually shows. The "headline" coverage figure cannot be drawn
from 3 F-points. Use this run as a *plumbing sanity check*; do not put
"rvSMR-AR is robust at weak IV" into the paper based on these numbers.

## Specific issues found (sorted by severity)

### S1 (severe — blocks publication of these numbers as the robustness demo)

#### S1.1 The "naive Wald comparator" is a non-standard, weak-IV-buffered estimator. The headline claim "AR maintains nominal coverage while naive Wald collapses" is NOT demonstrated by this design.

`run_tests.R:53-65` (`ivw_pooled_wald`) computes
```
pooled_b_x = Σ b_x_k · w_x_k / Σ w_x_k     # IVW pool of exposure betas
pooled_b_y = Σ b_y_k · w_y_k / Σ w_y_k     # IVW pool of outcome betas
β̂_naive   = pooled_b_y / pooled_b_x
```
This is **not** the classical Wald / IVW-of-ratios comparator that
collapses at weak IV in Wang-Kang Fig 6 and PLB Fig 2. The classical
comparator is
```
β̂_IVW = Σ w_k · (b_y_k / b_x_k) / Σ w_k       # per-IV ratio, then meta
```
which inherits weak-IV pathology because each per-IV Wald ratio
`b_y_k / b_x_k` has heavy-tailed Fieller distribution when `b_x_k ≈ 0`.

By pooling `b_x`'s *before* dividing, the implementation in
`run_tests.R` builds a `pooled_b_x` whose denominator is the **sum** of K
near-zero quantities, which is roughly √K times further from zero than
any individual `b_x_k`. At K=3 this gives a hidden 3× boost to first-stage
signal in the comparator. **The "Wald" comparator is itself partially
weak-IV-robust** — it cannot reproduce the canonical Wald collapse, and
so the headline contrast is a straw man with a hidden boost.

#### S1.2 The over-coverage of "naive Wald" (0.985 / 0.988 at F=1 / F=0.5) is described as collapse. It is the opposite — it is over-coverage.

`VALIDATION_REPORT.md` §5 (scenario C) claims:
> "naive Wald is *over*-conservative (0.985) due to inflated delta-method
> SE near the ratio singularity"

This is reading a coverage of 0.985 as a *failure* of the comparator —
but **a confidence interval that contains the truth 98.5% of the time
when nominal is 95% is conservative, not broken**. The canonical
Wang-Kang Fig 6 / PLB Fig 2 demonstration shows naive methods *under-cover*
(drop to e.g. 0.6, 0.4) at weak IV, with intervals that exclude the
truth. We are showing the wrong sign of mis-coverage.

This is partly explained by S1.1 (the comparator is not the standard
Wald) and partly by S3.1 (the DGP avoids the heavy-tail regime). Either
way, the **rhetoric in VALIDATION_REPORT.md §1 and §5** ("naive Wald
collapses") is unsupported by the data and must be deleted.

#### S1.3 The generative model is the AR statistic's own asymptotic regime; coverage near nominal is closer to a tautology than a test.

`generate_test_data.R:84-132` draws
```
b_x_k ~ N(α_k,      SE_x²)
b_y_k ~ N(β·α_k,    SE_y²)        (no overlap)
or jointly Cov(b_x, b_y) = ρ·SE_x·SE_y per mask
```
Then the AR moment `m(β₀) = b_y - β₀·b_x` has *exact* mean zero and
*exactly* the variance matrix `V(β₀)` that `mrAR_multi.R:155-157`
constructs, term for term. There is no model-misspecification slack
between the generator and the statistic. The pivotal statement
`AR ~ χ²_K under H₀` becomes an *algebraic identity* under finite-sample
normal draws (not even asymptotic), modulo numerical inversion.

So: scenarios **A, B, F** (no pleiotropy, normal residuals, R_xx=I)
are essentially **a numerics check on `mrAR_multi`**, not a validation of
"the rvSMR algorithm in real data." Coverage at 0.951 / 0.954 / 0.943 is
the correct outcome of a correct quadratic-form computation. **It does
not establish anything about whether AR will hold up in real summary
statistics.** This must be said explicitly in the report.

#### S1.4 Scenario C "weak IV" doesn't reproduce the canonical regime because F_target=1 is per-mask, joint F ≈ 3.

`generate_test_data.R:71` sets `α_k = sqrt(F_target) * SE_x` so each
mask has `E[F_k] = F_target + 1` (non-central χ²₁). At F_target=1 with
K=3 independent masks, the **joint concentration parameter** (the
quantity that controls Wang-Kang's r) is
```
λ_joint = K · F_target = 3
```
Wang-Kang's r=1 case (their canonical weak-IV cell) is **λ_joint = 1**
— roughly three times weaker than our scenario C in the relevant
identification sense. **The labelling "Weak IV (rvSMR headline regime)
F=1" in `scenarios.R:33` and the report's §3 is misleading.** What we
actually tested is "moderately weak with K=3 independent IVs" —
substantially easier than Wang-Kang r=1 because each of three
independent moments carries 1/3 of the concentration burden.

Worse: with K=3 the AR statistic has 3 df, not 1, so the chi-square
critical value is 7.815 (not 3.84). The CI is correspondingly wider per
mask — this is *baked-in conservatism* from K, not robustness from AR
math. To compare against Wang-Kang the scenario needs to fix the joint
concentration `λ_joint` not the per-mask F.

#### S1.5 The R_xy "harness mismatch" fix is justified by a citation that does not exist.

`VALIDATION_REPORT.md` §6 says:
> "per `steps_5_to_9_logic.md`, masks are *disjoint* annotation classes
> (pLoF / mis:LC / regulatory), so cross-mask LD/sample-overlap
> correlations are structurally zero"

I `grep`'d `steps_5_to_9_logic.md` for "disjoint", "legality",
"cross-mask", and the specific phrase "defaults legality" cited in
`run_tests.R:96-99`. **None of these appear in `steps_5_to_9_logic.md`
or `main.tex`.** The closest statement is `main.tex:204` (a pitfall
warning) which says the *opposite*: if masks share variants, R_xx /
R_yy must be supplied with off-diagonals. The disjointness argument
itself was invented to justify the diag-only fix.

Even granting variant-set disjointness across masks (which is a real
property of STAARpipeline annotation classes), **disjoint variant sets
do NOT imply uncorrelated estimators under sample overlap**. In a one-
sample (overlap) setting, the same individuals contribute genotype data
to *both* `b_x_pLoF` and `b_x_reg` — the score equations from a SAIGE-
QTL or burden regression on different variant subsets *share residuals*
through the per-individual outcome, inducing off-diagonal correlation in
the `b_x` vector independent of any LD between variants. Same for `b_y`.
The diagonal R_xy fix may be a numerical convenience, but the
justification offered does not hold; the validation of Step 9 (sample
overlap) is therefore incomplete. **A proper test needs a generator
that produces non-zero off-diagonals in R_xx, R_yy, and R_xy from a
shared-individuals score, and an inference run that supplies the
corresponding off-diagonals.** Until then, scenario F validates only the
strict block-diagonal case — which is approximately the easy case.

#### S1.6 The bias of +0.61 in scenario C is reported as a meaningful quantity but is not.

`results.md` row C: `bias β̂_AR = +0.6057` while coverage is 0.950.
`main.tex:303` explicitly flags this: "in weak-IV regimes the argmin
can sit far from the bulk of the accepted set." Computing the mean of a
diverging quantity (argmin of a near-flat objective, or fallback to
grid edge when `optimize` doesn't find a minimum on a `whole_line`-shape
56% of the time) over 1000 reps is **statistically meaningless** —
it averages over reps where β̂ is essentially `grid_pad_mult · max|Wald|`
boundary values.

The +0.61 number is being reported in a column labeled "bias" without
any caveat. A reader will read this as "AR estimator is biased high
at low F" — which is wrong; the AR *interval estimator* maintains
coverage. The argmin is not a point estimate one should report.

**Fix:** strip the "bias β̂_AR" column from row C and D, or replace
with "median β̂ over bounded-CI reps only" and footnote that the argmin
is ill-defined in whole-line shape.

### S2 (significant — needs follow-up before final framing)

#### S2.1 Tier 4 (confounder U → MR backdoor) is not tested at all.

The entire point of MR is backdoor confounding from an unobserved U
that affects both X and Y. `generate_test_data.R` line 86-87 has:
```
mu_x = α_vec
mu_y = β·α_vec + pleio_vec
```
There is no `U` term. The only "confounder analog" in the design is
the pleiotropy offset in scenario E, which is **direct horizontal
pleiotropy** (`α_k → b_y_k`), *not* unobserved confounding of `(X, Y)`
through a shared upstream cause. Whether the AR statistic has any
robustness to confounder-induced correlation between `b_x_k` and
`b_y_k` *within* a mask (beyond what `R_xy` captures) is **not tested**.

Real summary statistics from one-sample or overlapping cohorts can
carry exactly this kind of correlation. This is the largest gap.

#### S2.2 R_xx, R_yy never deviate from I_K.

`scenarios.R` never exercises `R_xx ≠ I` or `R_yy ≠ I`. The harness
sources `d$R_xx` and `d$R_yy` from `generate_test_data.R:155-156` which
hard-code `diag(K)`. So we never tested whether `mrAR_multi`'s
inversion of `V(β₀)` survives ill-conditioned `V_yy` from correlated
exposure-side estimates — even though `main.tex:204` and `:243`
explicitly warn that this is a known failure mode.

#### S2.3 Sargan-J 0.807 power at 5×SE pleiotropy is not as impressive as it sounds.

`generate_test_data.R:76` sets default `pleio_size = 5 · mean(SE_y)`.
With `SE_y = 1/√10000 = 0.01`, the offset is `0.05` — **12.5% of
β=0.4** on the outcome scale, applied to one of three masks. That is a
*moderately large* horizontal pleiotropy magnitude, not a subtle one
(Bowden's MR-PRESSO 2018 tests cover offsets down to ~1% of total
signal). 80% power at 5×SE is reasonable but not headline-worthy.

There is **no power sweep** across `pleio_size ∈ {0.5, 1, 2, 5, 10} ·
SE_y`. We don't know where J breaks down. The current cell is far from
the detection threshold; it answers a question reviewers won't ask
("does J work at large pleiotropy?") and skips the question they will
ask ("how small a pleiotropy can J detect?").

#### S2.4 The acceptance band [0.93, 0.97] is fragile single-seed evidence at F.

`results.md` row F shows coverage 0.943. Wald MC SE at n=1000 is
`√(0.95·0.05/1000) ≈ 0.0069`, so a "true" coverage of 0.95 will give
empirical values in `[0.936, 0.964]` with probability 95%. 0.943 is 1
SE below 0.95 — within MC tolerance but it is the *most extreme*
single value of the six and the one closest to the band edge. **With a
different master seed scenario F could drift to 0.93 or below.**

Two specific seed issues:
- `run_tests.R:108`: `seed_r = ((master_seed %% 1e6) * 1000 + r) %% int.max`
  combined with `master_seed + match(nm, names)` from line 203. The
  scenarios are seeded *consecutively* with `master + 1, +2, ..., +6`.
  Reps `r` within a scenario share a tight seed structure, but reps
  across scenarios may share common pseudo-random sub-sequences (the
  `(seed_A + r) → modular_seed` map is nearly linear). This means
  scenarios share more MC correlation than independent draws — the six
  coverage numbers are not 6 independent draws.
- The full 1000×6 run is one master seed. No claim about robustness to
  the seed is supported. Run the harness with 5 master seeds and
  report a range, or it is publishable as a point estimate only.

#### S2.5 K=1 and K=2 are not exercised in the simulation grid.

HANDOVER §2 (referenced in `VALIDATION_REPORT.md`) commits to "K ≥ 3"
in real analyses, but the rvMR package exposes `mrAR` (K=1, closed form)
and `mrAR_multi(K≥1)`. The K=1 cross-check is a numerics anchor (a
single test in `test-mrAR_multi.R:98-115`). **K=2 is the transition
case**: J has df=1, the CI is the disconnect-or-bounded regime, the
math is the cleanest non-trivial multi-IV check. Not testing K=2 is a
gap; real genes often have only 2 viable masks (a gene with no
detected regulatory burden).

#### S2.6 Step 10 (HEIDI-rv), Step 11 (annotation concord), Step 12 (sensitivity), Step 13 (E-value), Step 14 (cell-type Q) are 0% validated.

`audit.md` §3 lists 4 stubs (heidi_rv, annotation_concord, iv_partial_r2,
e_value) plus a missing function (cell-type concordance). The report
correctly flags these as not blocking AR coverage — but the headline
"6/6 PASS, no correctness bugs" is taken across the *whole package*
in §1, which is **overclaim**: the over-id axes and sensitivity
diagnostics are completely unvalidated. The HEIDI-rv `O(1/m)` caveat
in `main.tex:359` is the biggest landmine — a paper that ships the
algorithm without validating that weight choice will get caught by a
methods reviewer.

#### S2.7 Heterogeneous mask sizes (m_k) are not in the generator at all.

`main.tex` Step 10 (HEIDI-rv) and the per-mask burden-construction
discussion both depend on `m_k` (number of variants in mask k). Real
masks have 5-500 variants; the SE of `b_x_k` and `b_y_k` depends on
`m_k` and the per-variant signal. The current generator uses
`SE_x = 1/√n_x` independent of K, m_k, or per-mask MAF — so we cannot
say anything about how the calibration of `(b_x_k, SE_x_k)` from
SAIGE-QTL would survive a real burden where SE varies by 3× across
masks.

#### S2.8 Non-normality, finite-sample SE estimation, and burden-coefficient SE miscalibration are untested.

Real burden coefficients from SAIGE-QTL / Genebass have
heteroscedasticity-robust SE, sometimes miscalibrated when MAC is
small. The generator uses exact normals with known SE plug-in. **Any
real-data deployment will hit at least one of these regimes**, and we
have zero coverage data on them.

### S3 (cosmetic / clarification)

#### S3.1 The DGP cannot reproduce the heavy-tailed Fieller behavior of `b_y/b_x` because b_x is drawn from a Gaussian centered at a positive `α_k`.

`signs <- rep(1, K)` for K≤3 (`generate_test_data.R:70`) forces
`α_k > 0`. So `Pr(b_x_k < 0)` is governed by `Φ(-α_k/SE_x) =
Φ(-√F_target)`. At F_target=1, `Φ(-1) ≈ 0.16` — only 16% of `b_x_k`
draws cross zero, so the per-mask Wald ratio is bimodal-but-bounded for
most reps. Pooling K=3 such draws (S1.1) means
`Pr(pooled_b_x crosses zero)` is much smaller than `Φ(-√(K·F))`. **The
generator's geometry prevents the Wald comparator from ever exhibiting
its true weak-IV pathology.** This is a setup for S1.1 / S1.2.

#### S3.2 `generate_test_data.R:70` silently disables sign alternation for K ≤ 3.

The comment says "For the 'no sign alternation' stability test, force
all positive" but this is the *only* tested K (3). Sign-alternated `α_k`
would put some `b_x_k`'s near zero from below and some from above,
which is the realistic case. The all-positive case is the easy
case for the pooled-Wald comparator (denominator stays away from zero).

#### S3.3 The fallback in `mrAR_multi.R:328-331` can return grid-edge as `β_hat`.

When `optimize` fails or returns non-finite, the fallback is
`grid[which.min(ar_grid_finite)]`. For a flat AR surface (D / C
scenarios) this is effectively a random grid point. The reported "bias"
column is averaging over these grid-edge values. Confirms S1.6.

#### S3.4 Per-scenario sample sizes are all n_x = n_y = 10K.

The PLB 2024 / Yang 2025 simulation grids sweep n. Not sweeping n means
we don't know how coverage behaves at the realistic xQTL n ≈ 1K - 30K
range (a key gap for rare-variant burden tests where n is limited).

#### S3.5 "Modal CI shape" reporting in `VALIDATION_REPORT.md` §4 obscures the disconnected_union rate.

For scenario C, the table shows `whole_line (56.3%)` as modal — but
`disconnected_union (16.8%)` is the *interesting* shape (the rvSMR
"two-bumps" weak-IV signature, `main.tex:228` Eq. 4.6). Reporting only
the modal shape hides this. Show the full distribution in the headline
table.

#### S3.6 "Acceptance band [0.93, 0.97]" cited as "Patel-Lane-Burgess tolerance" in `VALIDATION_PLAN.md:109` — verify the citation.

PLB 2024 §4 reports coverage to 3 decimal places without an explicit
acceptance band. The [0.93, 0.97] band is a reasonable choice (Wald MC
SE at 1000 reps), but presenting it as the PLB tolerance is an
overstatement.

## What the validation *does* correctly establish (positive log)

- `mrAR_multi.R`'s **numerical machinery** (grid + uniroot + classification)
  is internally consistent with its analytical specification. At inputs
  drawn from the assumed Gaussian DGP, coverage matches χ²_K theory.
- The Sargan-J test at K=3 with one offset mask of magnitude 5×SE *does*
  detect pleiotropy at 80% power. Calibration of J ~ χ²_{K-1} appears
  correct at the strong-IV null (Type-I 0.038-0.052 across rows A, B, F).
- The CI-shape classifier produces the four expected shapes at the
  expected rates (whole_line dominates at F<1, bounded dominates at F>10).
- The single-line R_xy fix in `run_tests.R:100-104` is sufficient for
  the diag-only sample-overlap case (coverage moves from 0.899 → 0.943,
  back into MC tolerance).
- The harness is **reproducible** and the master seed is recorded. Re-run
  gives bit-identical results.

In short: this is a **good plumbing / regression test**, suitable as a
testthat-style CI gate. It is **not** a publishable robustness
demonstration.

## What the validation does *not* establish (gap log)

(Roughly ordered by how badly a reviewer would beat the paper for it.)

1. **No real LD-induced or sample-overlap-induced off-diagonals in
   R_xx, R_yy, R_xy.** Test cases R_xx ≠ I_K, R_yy ≠ I_K.
2. **No backdoor-confounding (Tier 4)** — the canonical MR threat model
   is untested.
3. **No combined-hardships scenario** (weak IV + pleiotropy + overlap
   simultaneously). Real genes will be in this corner.
4. **No comparison to a true classical Wald comparator** (per-IV ratio +
   meta-analysis). The current comparator is partially weak-IV-robust.
5. **No F sweep.** The headline "coverage vs F" figure cannot be drawn
   from {0.5, 1, 20} — three points. Wang-Kang Fig 6 has 5-10 F values.
6. **No pleiotropy magnitude sweep.** Don't know where J breaks down.
7. **No K sweep.** K ∈ {1, 2, 3, 5, 10} not tested; only K=3.
8. **No non-normality / finite-sample SE error.** Real data lives there.
9. **No heterogeneous mask sizes m_k.** HEIDI-rv math depends on m;
   ignored.
10. **HEIDI-rv (Step 10), annotation-class concordance (Step 11),
    sensitivity scalars (Step 12-13), cell-type Q (Step 14) all unvalidated.**
11. **No multi-master-seed stability check.** A single seed.
12. **No sub-sampling weak-IV induction** (Patel-Lane-Burgess Fig 5
    tactic).
13. **No selective-reporting / Goodhart simulation** (PLB Fig 3).
14. **No multi-null Type-I test** at β₀ ≠ 0 (Wang-Kang Fig 1 analog).
15. **Other comparators absent** (MR-RAPS, dIVW, MR-Egger, MR-PRESSO,
    RARE, MR-CARV). VALIDATION_PLAN.md §1 lists all of these as the
    expected comparator set; none implemented.

## Specific edits the user should apply to VALIDATION_REPORT.md

(Verbatim replacement text; copy-paste these in.)

### Replace §1 "Headline verdict"

**Current text:**
> rvMR 算法在 6 个测试场景下 6/6 通过覆盖率指标 ... The implemented AR core
> (mrAR, mrAR_multi, wald_burden, utils.R) is correct and behaves as
> predicted by Anderson-Rubin theory across all weak-IV regimes tested.
> **No correctness bugs were uncovered** in the package; one minor
> harness/inference mismatch in the sample-overlap test was diagnosed
> and fixed before the final run (see §6).

**Replace with:**
> Under a Gaussian summary-statistic data-generating process *that matches
> the AR statistic's own asymptotic regime term-for-term*, the
> `mrAR_multi` numerical implementation achieves empirical 95% coverage
> ∈ [0.943, 0.962] across 6 scenarios at K=3, master seed 20260603,
> 1000 reps each. This is a numerics regression test — it confirms that
> the grid + uniroot inversion in `mrAR_multi.R:114-360` correctly inverts
> the χ²_K level set on inputs drawn from the assumed model.
>
> **This run does NOT yet establish that rvSMR-AR is robust at weak IV
> in real summary statistics.** The DGP excludes by construction: (a)
> backdoor confounder-induced correlation between b_x and b_y, (b)
> non-zero R_xx / R_yy off-diagonals from LD or sample overlap, (c)
> heterogeneous per-mask SE calibration, (d) heavy-tail behavior of the
> per-IV Wald ratio at b_x near zero. The naive Wald comparator
> implemented in `run_tests.R:53-65` is the IVW-pooled-then-divided form,
> which is partially weak-IV-buffered and does NOT exhibit the canonical
> Wald collapse from Wang-Kang Fig 6; the reported "naive over-coverage"
> at F=1 (0.985) and F=0.5 (0.988) is therefore consistent with the
> comparator's residual robustness, not with a failure mode. The
> publishable coverage-vs-F figure (`VALIDATION_PLAN.md` §3 Track 1)
> remains to be produced; see §7 for the gap list.

### Replace §5 scenario C interpretation

**Current text:**
> AR's headline regime: coverage **must hold** at F≈1. **PASS** -- AR is
> doing its job ... naive Wald is *over*-conservative (0.985) due to
> inflated delta-method SE near the ratio singularity ... AR is doing
> its job: nominal coverage maintained at F≈1 where standard Wald is
> unreliable.

**Replace with:**
> Coverage 0.950 is on-spec for the AR set at K=3, joint-concentration
> λ ≈ 3 (per-mask F_target=1 with 3 independent masks). Note this is a
> **stronger** identification regime than Wang-Kang's r=1 cell (which
> uses 100 IVs at λ=1, much closer to non-identification). Naive Wald
> over-coverage at 0.985 is consistent with the pooled-Wald comparator
> (`run_tests.R:53-65`) being partially weak-IV-buffered, NOT the
> canonical Wald collapse. The reported `bias β̂_AR = +0.6057` is **not
> a meaningful estimator bias** at this F: 56.3% of CIs are
> `whole_line`, so the AR argmin is fallback grid-edge in those reps
> (see `mrAR_multi.R:328-331`). Strip the bias column for C and D.

### Replace §5 scenario D interpretation

**Current text:**
> F<1: expect many whole_line CIs, but coverage should still ≥ 0.95
> (conservative). PASS

**Replace with:**
> Coverage 0.962 is dominated by the 69.7% `whole_line` rate at F=0.5:
> a `whole_line` set trivially contains β=0.4, so the "coverage"
> measurement reduces to "fraction of non-empty, non-whole-line reps
> that bracket β" + "all whole-line reps count as covered." This is not
> a falsifiable check of AR — it would have passed at F=0.01 too. The
> meaningful check at this F is the **rate of non-empty bounded CIs
> that exclude β** (Type-II error against β); we do not report that.

### Replace §6 last paragraph

**Current text:**
> No bugs found in the rvMR package itself. All deviations from nominal
> were either ... harness mismatch above. The rvMR package was NOT
> modified during this run.

**Add after:**
> **The justification for the R_xy diag-only fix offered above ("per
> steps_5_to_9_logic.md masks are disjoint annotation classes ...
> structurally zero") cites a section that does not exist in
> `steps_5_to_9_logic.md` or `main.tex`. Mask-variant-disjointness does
> NOT imply zero off-diagonal correlations in R_xx / R_yy / R_xy under
> sample overlap (shared per-individual residuals induce off-diagonals
> independent of the variant LD). The diag-only fix is therefore the
> correct fix for the *no-overlap, no-shared-residual* case the
> generator implements, but the validation of Step 9 sample-overlap
> handling is incomplete — a generator that produces realistic
> off-diagonal R_xx / R_yy from a shared-cohort score is required.**

### Add a new §0 caveat banner before §1

> **Caveat.** This report documents a 6-cell numerical regression test
> of `mrAR_multi`, with the data-generating process intentionally
> matched to the AR null distribution. It is not the publishable
> "coverage vs F" figure (VALIDATION_PLAN.md §3 Track 1). Numerical
> targets are met; substantive robustness claims for the paper require
> the extended grid in §7.

## Recommended next experiments (ranked)

1. **Fix the comparator.** Replace `ivw_pooled_wald` with the
   *classical* IVW-of-ratios + delta-SE form `β̂ = Σw_k·r_k / Σw_k`
   where `r_k = b_y_k/b_x_k` and `Var(r_k)` is per-IV delta SE. Show
   the classical comparator collapses (under-covers) at weak IV. This
   alone restores the headline "AR robust, Wald collapses" claim.
2. **Sweep F.** Run F_target ∈ {0.25, 0.5, 1, 2, 5, 10, 20, 50} as a
   single new scenario list, plot coverage vs F with MC bars. This is
   ~1 evening + 2 hours of compute (already in
   `VALIDATION_REPORT.md` §7 step 1).
3. **Sweep pleio_size.** Run pleio_size ∈ {0.5, 1, 2, 3, 5, 10} × SE_y
   at F=20, plot J power as a function of effect magnitude. Locate the
   detection threshold.
4. **Add a confounder (Tier 4).** Add `U_i` to the generator: draw
   `U ~ N(0, σ_U²)`, then `b_x_k` has additional noise from U, and
   `b_y_k = β·α_k + γ·U + ε`. Vary `γ` to control confounding strength.
   Test whether AR maintains coverage under classical MR confounder
   leak. This is the canonical "why MR matters" simulation.
5. **Sweep K.** K ∈ {1, 2, 3, 5, 10}. Verify K=1 reduces to mrAR's
   closed form numerically; verify K=2 J-test calibration; verify K≥5
   doesn't degrade grid inversion.
6. **Multi-seed stability.** Run the existing 6-cell grid with 5
   different master seeds. Report range and confirm 0.943 is not a
   seed-dependent edge of the acceptance band.
7. **Off-diagonal R_xx, R_yy stress test.** Generate b_x with R_xx
   non-trivial (e.g., AR(1) at ρ=0.3); supply correct R_xx to
   mrAR_multi; verify coverage. Then *withhold* R_xx (supply I_K) and
   show the coverage gap — quantifies the cost of misspecifying R_xx.
8. **Implement HEIDI-rv with `CompQuadForm::davies`** and validate the
   eigenvalues-of-V_δ weight choice flagged at `main.tex:359` (the
   "headline trap"). This is the highest-prior-art-overlap reviewer
   risk in the whole package.
9. **Plug in MR-RAPS / dIVW / MR-Egger / MR-PRESSO as additional
   comparators.** VALIDATION_PLAN.md §3 lists them as the expected
   set. Without them the paper has no comparator panel.
10. **Sub-sampling tactic (PLB Fig 5).** Once a real common-variant
    Track 2 substrate is available (e.g., PCSK9 cis-eQTL × LDL),
    sub-sample and produce the heatmap. This is the cheapest tactic
    that produces a real-data-like figure.

---

*End of adversarial review. Final disposition: USE THIS RUN AS REGRESSION
PLUMBING, NOT AS THE HEADLINE ROBUSTNESS DEMO. Rerun items 1-3 above
before drafting any "rvSMR-AR is robust" paragraph for the paper.*
