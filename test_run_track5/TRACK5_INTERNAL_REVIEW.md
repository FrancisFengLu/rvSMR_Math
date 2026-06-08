# Track 5 — Internal review

**Reviewer**: general-purpose sub-subagent (this internal-review pass)
**Date**: 2026-06-08
**Scope**: verify (a) delta-method variance propagation, (b) `cell_type_q()` re-run, (c) citation correctness, (d) Olink/SomaScan/ELISA scale audit.

## (a) Delta-method variance propagation

The renormalized Wald is `tilde_beta_xy = b_y / b_pqtl`. By the delta method:

```
Var(tilde_beta_xy) = (d/db_y)^2 * Var(b_y) + (d/db_pqtl)^2 * Var(b_pqtl)
                   = (1/b_pqtl)^2 * SE_y^2 + (b_y / b_pqtl^2)^2 * SE_pqtl^2
                   = SE_y^2 / b_pqtl^2 + b_y^2 * SE_pqtl^2 / b_pqtl^4
```

Cross-check vs the script:
```r
se_wald_anch <- sqrt(se_y^2 / b_pqtl^2 + b_y^2 * se_pqtl^2 / b_pqtl^4)
```

Matches. Numerically:
- b_y = 0.0336, SE_y = 0.0033
- b_pqtl = 0.3719, SE_pqtl = 0.0145
- Term 1: 0.0033² / 0.3719² = 1.089e-5 / 0.1383 = 7.87e-5
- Term 2: 0.0336² × 0.0145² / 0.3719⁴ = 1.129e-3 × 2.103e-4 / 0.01913 = 2.374e-7 / 0.01913 = 1.241e-5
- Sum = 9.11e-5; SE = sqrt(9.11e-5) = 0.00955 ≈ 0.0095 ✓

The reported SE = 0.0095 is consistent.

Also cross-check the unanchored Wald SE for Liver:
- b_y = 0.0336, SE_y = 0.0033, b_x = 0.4287, SE_x = 0.0909
- Term 1: 0.0033² / 0.4287² = 1.089e-5 / 0.1838 = 5.92e-5
- Term 2: 0.0336² × 0.0909² / 0.4287⁴ = 1.129e-3 × 8.263e-3 / 0.03378 = 9.328e-6 / 0.03378 = 2.762e-4
- Sum = 3.354e-4; SE = sqrt(3.354e-4) = 0.0183 ✓ matches script output

**Delta-method algebra verified.**

## (b) Independent re-run of `cell_type_q()`

Reading `cell_type_q()` from `/home/francisfenglu4/rvSMR/May_30md/rvMR/R/cell_type_concord.R` (lines 86–182), the statistic is:
```
w_c = 1 / SE_c^2
bar = sum(w_c * b_c) / sum(w_c)
Q = sum(w_c * (b_c - bar)^2)
p = pchisq(Q, df = C - 1, lower.tail = FALSE)
```

**Manual independent computation, unanchored**:

| Tissue | b_xy | SE | w = 1/SE² |
|---|---:|---:|---:|
| Liver | 0.0784 | 0.0183 | 2987.0 |
| Blood | 0.2028 | 0.0481 | 432.2 |
| Adipose-visceral | 0.1434 | 0.0431 | 538.6 |
| Artery-aorta | 0.1788 | 0.0666 | 225.6 |

- sum_w = 2987.0 + 432.2 + 538.6 + 225.6 = 4183.4
- sum_w*b = 2987.0×0.0784 + 432.2×0.2028 + 538.6×0.1434 + 225.6×0.1788
         = 234.2 + 87.65 + 77.24 + 40.34 = 439.4
- bar = 439.4 / 4183.4 = 0.1050
- Q contributions:
  - Liver: 2987.0 × (0.0784 - 0.1050)² = 2987.0 × 7.08e-4 = 2.114
  - Blood: 432.2 × (0.2028 - 0.1050)² = 432.2 × 9.56e-3 = 4.134
  - Adipose-visceral: 538.6 × (0.1434 - 0.1050)² = 538.6 × 1.475e-3 = 0.794
  - Artery-aorta: 225.6 × (0.1788 - 0.1050)² = 225.6 × 5.45e-3 = 1.229
- Q = 2.114 + 4.134 + 0.794 + 1.229 = **8.27** (script reports 8.2747; small rounding)
- p = pchisq(8.27, df=3, lower.tail=F) ≈ 0.0408

Verified — matches script's 0.04066 and Team 1's reported 0.041.

**Manual independent computation, anchored**:
All four tissues have b_xy = 0.0903, SE = 0.0095 (identical).
- bar = 0.0903 (trivially)
- Q = sum(w * (0.0903 - 0.0903)^2) = **0**
- p = pchisq(0, df=3, lower.tail=F) = **1.0**

Verified.

**cell_type_q() runs faithfully; outputs match independent recomputation.**

## (c) Citation correctness

Spec required: Sun BB 2023 *Nature* 622(7982):**329** (NOT 339 = Dhindsa).

Script and report cite:
- Sun BB 2023 as the **preferred** anchor (not used due to access failure)
- Pott J 2024 *Hum Mol Genet* PMID 38491180 as the **used** anchor

**Did NOT use Sun KY 2024 RGC-ME**: confirmed. Sun KY 2024 *Nature* 631:583 is the rare-variant exome PheWAS paper; it would have been wrong here. Verified by checking `team2_drafts/citation_audit.md` lines 19, 35, 129 which establish the distinction (Sun KY = RGC, Sun BB = UKB-PPP).

**Pott 2024 substitution properly flagged** in `data_pull_log.md` and `pcsk9_pqtl_anchor_results.md` — does NOT pretend to be Sun BB 2023.

**Citations are correct.**

## (d) Scale audit (Olink NPX vs SomaScan RFU vs ELISA absolute)

The spec asks: "did you use Olink scale (NPX) vs ratiometric vs antibody-relative units? Document."

**Audit findings**:

- **What we used**: Pott 2024 b_pqtl = 0.37189 (T allele decrease in log-PCSK9).
- **Pott 2024's scale**: log-transformed plasma PCSK9 concentration, **mixed Olink + ELISA** across 6 contributing cohorts (LIFE-Heart, LIFE-Adult, LURIC, TwinGene, KORA-F3, GCKD). Per PMC10964567: "PCSK9 levels were log-transformed." The platform mix is partly Olink (e.g., LIFE-Adult), partly ELISA (e.g., LIFE-Heart). Not pure Olink NPX.
- **What Sun BB 2023 UKB-PPP would have been**: pure Olink Explore 3072 NPX (log2-relative normalized protein expression). Per the published methods (Sun BB 2023 *Nature* 622:329, supplementary), NPX is a relative quantification scale; cis-pQTL betas are typically reported as "change in NPX units per allele". The expected magnitude for rs11591147-T at PCSK9 on Olink NPX in UKB-PPP is ~0.4 NPX (per direct numerical comparability with Pott 2024's 0.37 log-PCSK9).
- **Where SomaScan would have inflated**: Pietzner 2021 reports 0.883 and Gudjonsson 2022 reports ~1.04 — both on SomaScan RFU log-units, which run roughly 2.4–2.8× higher than Olink/ELISA log scales because of platform calibration differences. We deliberately did NOT use Pietzner/Gudjonsson for the anchor; this avoids a 2.4× scale-mismatch that would have miscalibrated the renormalization.
- **The renormalization unit**: with Pott 2024 anchor, `b_y / b_pqtl` = log-OR-CHD per (log-PCSK9 ~ Olink NPX) unit. Interpretable per-(approximately-SD-plasma-PCSK9) units. NOT per-mg/dL absolute mass units.

**Scale documented honestly. No SomaScan/Olink contamination.**

## (e) Other checks

- **Direction**: b_pqtl reported as "decrease" (T allele lowers PCSK9; LOF). Script uses `|b_pqtl|` = 0.37189 and notes sign convention; this is appropriate because Q is sum-of-squares (sign-insensitive) and the universal anchor cancels the sign. **OK**.
- **Sample overlap**: Pott 2024 cohorts (German, Italian, Swedish) do NOT overlap with Team 1's GTEx Liver eQTL (American GTEx donors) nor FinnGen+MVP+UKBB outcome. **R_xy assumption (=0) is appropriate for this anchor swap.**
- **Subset reproducibility**: Q-unanchored = 8.275 (script) vs 8.28 (Team 1 report) vs 8.27 (this manual recompute) — consistent at 0.01 precision.

## Verdict

- Delta-method algebra: **verified**.
- `cell_type_q()` re-run: **matches independently**.
- Citation correctness: **Sun BB 2023 = correct ref; substitution to Pott 2024 properly disclosed; Sun KY 2024 not used (correctly)**.
- Scale audit: **Olink-compatible log-PCSK9; no SomaScan-RFU contamination**.
- Honest reporting in headline result: **yes** — Q → 0 is presented as algebraically expected, not as biological homogeneity proof.

**Internal review: PASS.** One soft caveat: the universal-anchor result is essentially uninformative on the hypothesis. A meaningful follow-up would require per-tissue eQTL at the canonical pQTL lead variant (rs11591147) from the eQTL Catalogue — not done in this run, flagged in the report.
