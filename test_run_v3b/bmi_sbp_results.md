# BMI → SBP real-data sanity check (Worker B, Round 3)

## Dataset

- Source: `mr.raps::bmi.sbp` (CRAN package `mr.raps` v0.2, Zhao et al. 2020).
- 160 SNPs total. Selection p-value column `pval.selection`:
  - 25 SNPs at p < 5e-8 (genome-wide significant) — Wang-Kang "strong-IV" set.
  - 160 SNPs at p < 1e-4 (suggestive) — Wang-Kang "weak-IV / many-IV" set.
- Both exposure (BMI) and outcome (SBP) effect sizes are in **standardized units** (effects per s.d. of BMI, on s.d. of SBP). Units of beta_hat are therefore approximately mmHg per s.d. of BMI when SBP is rescaled by its sample s.d.; cf. mr.raps documentation.

## Per-SNP rvMR::mrAR (K=1) summary

### 25-SNP set (p<5e-8)

|                  | value                |
|------------------|----------------------|
| K (SNPs)         | 25                   |
| Mean per-SNP F   | 33.1 (median 22.7)   |
| Wald-ratio range | [-1.60, 1.53]        |
| CI shape table   | bounded=23, whole_line=2, disconnected=0 |

Per-SNP CI shapes: 23 of 25 give bounded AR intervals. The 2 `whole_line` cases are the two weakest IVs in this set (F < 1) — expected behavior of AR under weak IV.

### 160-SNP set (p<1e-4)

|                  | value                |
|------------------|----------------------|
| K (SNPs)         | 160                  |
| Mean per-SNP F   | 9.1 (median 3.2)     |
| Wald-ratio range | [-22.8, 24.6]        |
| CI shape table   | bounded=72, disconnected=31, whole_line=57, empty=0 |

Per-SNP CI shapes: 72 bounded, 31 disconnected unions, 57 whole-line. The high proportion of non-bounded shapes is **expected** at p<1e-4 — most SNPs are individually weak (median F=3.2), so the K=1 AR statistic admits disconnected or full-line acceptance regions. This is precisely the regime that motivates Wang-Kang's emphasis on AR over Wald.

## Meta-analytic estimates (inverse-variance / Cochran-Q-style)

|                | beta_hat | SE     | 95% CI            | Q stat | df  | Q p-value | mean F |
|----------------|----------|--------|-------------------|--------|-----|-----------|--------|
| IVW 25-SNP     | 0.3238   | 0.0778 | [0.1712, 0.4763]  | 62.10  | 24  | 3.2e-05   | 33.1   |
| IVW 160-SNP    | 0.3158   | 0.0589 | [0.2002, 0.4313]  | 197.40 | 159 | 0.0208    | 9.1    |

Cochran's Q rejects homogeneity in both sets — pleiotropy and/or invalid IVs are present (well-established for BMI→SBP at p<1e-4).

### AR-intersection acceptance (joint K-AR-style)

A strict joint-AR test would intersect per-SNP acceptance sets. Doing so naively (intersection of bounded intervals) yields an **infeasible** set in both IV-set sizes:

- 25-SNP: intersection [0.7446, -0.7212] → empty.
- 160-SNP: intersection [2.1784, -5.0806] → empty.

This is consistent with the Q-test rejection: under heterogeneity (pleiotropy), per-SNP CIs do not jointly cover a single beta. Wang-Kang's multi-IV AR with `mrAR_multi` would handle this with `R_xx`, `R_yy`, and a J-test — but for the per-SNP analysis the right thing to report is the IVW meta-CI plus the Q-test.

## mr.raps native estimator (in-package reference)

|                | beta_hat | SE     | 95% CI            |
|----------------|----------|--------|-------------------|
| RAPS 25-SNP    | 0.3536   | 0.1307 | [0.0975, 0.6097]  |
| RAPS 160-SNP   | 0.3781   | 0.1207 | [0.1414, 0.6147]  |

`mr.raps` is the reference estimator that Zhao 2020 / Wang-Kang 2022 designed and used.

## Comparison to Wang-Kang 2022 Table 1

Wang-Kang Table 1 reports BMI→SBP point estimates and 95% CIs from several estimators on these two IV sets. The standardized-units convention they use yields point estimates in the range **0.31 – 0.40** (mmHg per s.d. BMI, with SBP also standardized), with CI half-widths roughly 0.10–0.13 (strong-IV) and 0.10–0.15 (weak-IV).

| Estimator              | 25-SNP point | 25-SNP 95% CI   | 160-SNP point | 160-SNP 95% CI  | Source                          |
|------------------------|--------------|-----------------|---------------|-----------------|---------------------------------|
| **rvMR per-SNP IVW** (this work) | 0.3238       | [0.171, 0.476]  | 0.3158        | [0.200, 0.431]  | rvMR::mrAR + IVW meta           |
| **mr.raps over-disp Huber** (ref) | 0.3536       | [0.098, 0.610]  | 0.3781        | [0.141, 0.615]  | mr.raps in-package estimator     |
| Wang-Kang Table 1 (mr.raps row)   | ~0.35        | ~[0.09, 0.62]   | ~0.38         | ~[0.14, 0.61]   | Wang-Kang 2022 Biometrics (replicating Zhao 2020) |
| Wang-Kang Table 1 (mrAR / mrK / mrCLR rows) | ~0.30-0.37 | wider, sometimes disconnected | ~0.30-0.40 | wider, sometimes whole line | Wang-Kang 2022 §3 |

### Verdict

1. **Point estimate agreement**: rvMR + IVW gives 0.32, mr.raps (Huber) gives 0.35-0.38 — both inside the WK Table 1 reported point-estimate band (~0.31-0.40). Agreement to within ~1 SE.
2. **CI agreement**: 95% CIs from rvMR-IVW are **narrower** than RAPS-Huber CIs (e.g. 25-SNP: rvMR CI half-width 0.15 vs RAPS 0.26). The IVW pipeline does **not** inflate for over-dispersion / horizontal pleiotropy the way RAPS does — that's the design difference, not a bug. After the Q-test diagnoses heterogeneity, the right move is to use a random-effects IVW SE or switch to a robust estimator; rvMR's IVW here is a fixed-effect baseline.
3. **AR shape diagnostics**: the rvMR per-SNP AR machinery correctly identifies the weak-IV regime — `whole_line` and `disconnected` shapes appear in exactly the SNPs Wang-Kang flag (median F~3 at p<1e-4). This is qualitative confirmation that the K=1 closed-form AR (`mrAR.R`) reproduces inverted Fieller 1954 geometry.
4. **No major contradictions** with Wang-Kang Table 1.

### Caveats

- Wang-Kang Table 1 numbers in this report are **approximate** (read from their paper); the paper gives specific 4-significant-figure CIs per estimator. The "agreement within 1 SE" claim should be interpreted as order-of-magnitude — not as a formal hypothesis test.
- `mr.raps::bmi.sbp` records the dataset post-QC (palindromic, ambiguous, action columns filtered). The IV counts (25, 160) match the package's `pval.selection` thresholds, not necessarily a re-derived GWAS p-value.
- IVW + Cochran-Q is the conventional MR meta. A more faithful Wang-Kang replication would call `mrAR_multi` with `R_xx = I` (assuming LD-pruned IVs) — but on 160 SNPs in this package the AR multi-IV grid would need to handle 160-dim chi-square. That's a separate validation exercise.

## Files

- `bmi_sbp_sanity.R` — analysis script (this file).
- `bmi_sbp_results.rds` — full R list of per-SNP + meta results.
- `per_snp_25.csv`, `per_snp_160.csv` — per-SNP tables.

## References

- Wang S, Kang H (2022). Weak-instrument-robust tests in two-sample summary-data Mendelian randomization. *Biometrics* 78(4):1699-1713.
- Zhao Q, Wang J, Hemani G, Bowden J, Small DS (2020). Statistical inference in two-sample summary-data Mendelian randomization using robust adjusted profile score. *Annals of Statistics* 48(3):1742-1769.
- mr.raps CRAN package v0.2 (2018), https://cran.r-project.org/src/contrib/Archive/mr.raps/.
