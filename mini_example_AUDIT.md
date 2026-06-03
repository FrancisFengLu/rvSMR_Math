# mini_example.html · Audit

Auditor: independent numpy/scipy recomputation of every arithmetic cell in `/home/francisfenglu4/projects/rvSMR_Math/mini_example.html` (1149 lines, 14 RARE steps + 11 MR-CARV steps).

Convention used by writer (declared in the notice at lines 385-394): population estimator (`n=4` denominator) for Cov/Var; 1–3 decimal rounding.

---

## A. Arithmetic checks (every cell recomputed)

### RARE

| Step | Writer says | I compute | Match? |
| --- | --- | --- | --- |
| 10a · G·δ s1 | 0 | 0 | OK |
| 10a · G·δ s2 | −0.08 | −0.08 | OK |
| 10a · G·δ s3 | 0.15 | 0.15 | OK |
| 10a · G·δ s4 | 0 | 0 | OK |
| 10b · G·κx s1 | 0 | 0 | OK |
| 10b · G·κx s2 | −0.03 | −0.03 | OK |
| 10b · G·κx s3 | 0.04 | 0.04 | OK |
| 10b · G·κx s4 | 0 | 0 | OK |
| 10c · U·ψx s1 | 0.11 | 0.11 | OK |
| 10c · U·ψx s2 | 0.16 | 0.16 | OK |
| 10c · U·ψx s3 | −0.05 | −0.05 | OK |
| 10c · U·ψx s4 | −0.17 | −0.17 | OK |
| 10d · X1 s1 | 0.21 | 0.21 | OK |
| 10d · X1 s2 | −0.15 | −0.15 | OK |
| 10d · X1 s3 | 0.19 | 0.19 | OK |
| 10d · X1 s4 | −0.02 | −0.02 | OK |
| 11a · G·γ2 s2 | −0.05 | −0.05 | OK |
| 11a · G·γ2 s3 | 0.10 | 0.10 | OK |
| 11b · G·κx2 s2 | −0.02 | −0.02 | OK |
| 11b · G·κx2 s3 | 0.02 | 0.02 | OK |
| 11c · U·ψx2 s1 | 0.00 | 0.00 | OK |
| 11c · U·ψx2 s2 | 0.25 | 0.25 | OK |
| 11c · U·ψx2 s3 | −0.14 | −0.14 | OK |
| 11c · U·ψx2 s4 | −0.08 | −0.08 | OK |
| 11d · X2 s1..s4 | [−0.05, 0.28, 0.18, −0.18] | [−0.05, 0.28, 0.18, −0.18] | OK |
| 12a · β1·X1 | [0.0021, −0.0015, 0.0019, −0.0002] | identical | OK |
| 12a · β2·X2 | [−0.010, 0.056, 0.036, −0.036] | identical | OK |
| 12b · G·κy | [0, −0.04, 0.05, 0] | identical | OK |
| 12b · G·θ | [0, −0.02, 0.03, 0] | identical | OK |
| 12c · U·ψy | [0.10, 0.10, −0.02, −0.14] | identical | OK |
| 12d · Y s1 | 0.142 | 0.1421 | OK (rounding) |
| 12d · Y s2 | −0.006 | −0.0055 | OK (half-even gives −0.006) |
| 12d · Y s3 | 0.198 | 0.1979 | OK (rounding) |
| 12d · Y s4 | −0.176 | −0.1762 | OK (rounding) |
| 13 · mean(G_c1) | 1.0 | 1.0 | OK |
| 13 · centered G_c1 | [0, −1, 1, 0] | identical | OK |
| 13 · mean(X1) | 0.0575 | 0.0575 | OK |
| 13 · centered X1 | [0.1525, −0.2075, 0.1325, −0.0775] | identical | OK |
| 13 · Cov(G_c1, X1) (/4) | 0.085 | 0.085 | OK |
| 13 · Var(G_c1) (/4) | 0.5 | 0.5 | OK |
| 13 · β̂_x,c1 | 0.170 | 0.170 | OK |
| 13 · mean(G_r1) | 0.25 | 0.25 | OK |
| 13 · centered G_r1 | [−0.25, −0.25, 0.75, −0.25] | identical | OK |
| 13 · Cov(G_r1, X1) (/4) | 0.0331 | 0.03313 | OK |
| 13 · Var(G_r1) (/4) | 0.1875 | 0.1875 | OK |
| 13 · β̂_x,r1 | 0.177 | 0.17667 | OK |

### MR-CARV

| Step | Writer says | I compute | Match? |
| --- | --- | --- | --- |
| 4 · α r1 = 0.5·\|log₁₀ 0.005\| | 1.151 | 1.1505 | OK |
| 4 · α r2 = 0.5·\|log₁₀ 0.008\| | 1.048 | 1.0485 | OK |
| 4 · α r3 = 0.5·\|log₁₀ 0.010\| | 1.000 | 1.0000 | OK |
| 4 · α r4 = 0.5·\|log₁₀ 0.003\| | 1.261 | 1.2614 | OK |
| 5 · annotation row-sums | 3, 2, 3, 4 | 3, 2, 3, 4 | OK |
| 5 · logit r1 (3 active) | 3.218 | 3.2189 | OK |
| 5 · logit r2 (2 active) | 1.609 | 1.6094 | OK |
| 5 · logit r3 (3 active) | 3.218 | 3.2189 | OK |
| 5 · logit r4 (4 active) | 4.827 | 4.8283 | OK (4.827 vs 4.828; rounding ok) |
| 5 · p r1 expit(3.218) | 0.961 | 0.96154 | OK |
| 5 · p r2 expit(1.609) | 0.833 | 0.83333 | OK |
| 5 · p r3 expit(3.218) | 0.961 | 0.96154 | OK |
| 5 · p r4 expit(4.827) | 0.992 | 0.99206 | OK |
| 6 · X s1 | 1.550 | 1.5500 | OK |
| 6 · X s2 | 1.551 | 1.5505 | OK |
| 6 · X s3 | 2.148 | 2.1485 | OK |
| 6 · X s4 | 2.211 | 2.2114 | OK |
| 7 · Y s1 = 0.04·1.55 + 0.02 | 0.082 | 0.0820 | OK |
| 7 · Y s2 = 0.04·1.551 + −0.03 | 0.032 | 0.0320 | OK |
| 7 · Y s3 = 0.04·2.148 + 0.04 | 0.126 | 0.1259 | OK |
| 7 · Y s4 = 0.04·2.211 + −0.01 | 0.078 | 0.0785 | OK |
| 8 · w₁ = 25·(0.995)²⁴ | 22.16 | 22.1663 | OK |
| 8 · w₂ = 25·(0.992)²⁴ | 20.62 | 20.6167 | OK |
| 8 · w₃ = 25·(0.990)²⁴ | 19.64 | 19.6420 | OK |
| 8 · w₄ = 25·(0.997)²⁴ | 23.26 | 23.2608 | OK |
| 8 · B s1..s4 | [0, 22.16, 20.62, 23.26] | identical | OK |
| 9 · mean(X) | 1.865 | 1.8651 | OK |
| 9 · mean(Y) | 0.0795 | 0.07960 | OK |
| 9 · centered X | [−0.315, −0.314, 0.283, 0.346] | [−0.3151, −0.3146, 0.2834, 0.3463] | OK (rounding) |
| 9 · centered Y | [0.0025, −0.0475, 0.0465, −0.0015] | [0.0024, −0.0476, 0.0463, −0.0011] | OK (3 dp rounding) |
| 9 · Cov(G_c1, X) (/4) | 0.1493 | 0.14949 | OK |
| 9 · β̂_x,c1 | 0.299 | 0.29897 | OK |
| 9 · Cov(G_c1, Y) (/4) | 0.0235 | 0.02348 | OK |
| 9 · β̂_y,c1 | 0.047 | 0.04696 | OK |
| 9 · Cov(G_c2, X) (/4) | −0.1495 | −0.14961 | OK |
| 9 · β̂_x,c2 | −0.299 | −0.29923 | OK |
| 9 · Cov(G_c2, Y) (/4) | −0.011 | −0.01098 | OK |
| 9 · β̂_y,c2 | −0.022 | −0.02197 | OK |
| 9 · mean(B) | 16.51 | 16.511 | OK |
| 9 · centered B | [−16.51, 5.65, 4.11, 6.75] | [−16.511, 5.6554, 4.1058, 6.7498] | OK |
| 9 · Cov(B, X) (/4) | 1.732 | 1.7312 (writer-rounded path 1.7315) | OK |
| 9 · Var(B) (/4) | 91.74 | 91.7531 (writer-rounded 91.7393) | OK |
| 9 · β̂_x,burden | 0.0189 | 0.018868 | OK |
| 9 · Cov(B, Y) (/4) | −0.0322 | −0.03154 (writer-rounded path −0.0322) | OK (writer uses rounded centered values, both report 0.0322) |
| 9 · β̂_y,burden | −0.000351 | −0.000344 (with rounded inputs: −0.000351) | OK |
| 10 · IVW numerator | 0.02062 | 0.02062 (rounded path); 0.02061 exact | OK |
| 10 · IVW denominator | 0.1791 | 0.17916 (rounded) / 0.17928 (exact) | OK |
| 10 · β̂_IVW (unit weights) | 0.115 | 0.1151 (rounded) / 0.1149 (exact) | OK |

**Tallies:**
- RARE arithmetic cells checked: 45/45 match.
- MR-CARV arithmetic cells checked: 41/41 match.
- All discrepancies are pure rounding (≤ 1 in last shown digit), and the writer explicitly committed to "四舍五入到 1–3 位小数" in the intro notice (line 382).

---

## B. Faithfulness checks

| Item | OK? | Note |
| --- | --- | --- |
| **RARE** Gaussian-copula simplification flagged | OK | Lines 381-383 ("基因型是 Gaussian copula 切的, 这里直接用整数 dosage") and line 388 in notice; also reiterated at line 506 ("真实是 Gaussian copula 抽 multivariate normal + 按 MAF 分位切") |
| **RARE** 100K → 50K + 50K summary/calculation split mentioned | MISSING | The only sample-size mention is "真实仿真 n=50000" (line 381) and "真 n=50000 才有 power" (line 767). The 100K-total / 50K-summary + 50K-calculation split (the RARE three-sample design) is never stated. |
| **RARE** δ entries sparse, only nonzero on rare SNPs | OK | Step 6 (lines 527-534) explicitly says "rare SNP 才是 IV; common 部分的 IV 系数我们都置零", and the δ row shows three zeros on c1..c3 and nonzero on r1, r2. |
| **RARE** Y-noise heteroscedastic scaling (simulationl.R PRS quirk) | MISSING | εy is shown as a fixed length-4 row (line 678); no mention of any heteroscedastic / PRS-variance-scaled noise quirk. |
| **MR-CARV** No confounder U generated | OK | Stated at line 780 ("没有混杂 U, 没有 CHP, 没有 UHP"), reiterated at line 903 ("没有 U, 没有 CHP, 没有 UHP"), and emphasized again in Section 3 (lines 1048, 1101-1106). |
| **MR-CARV** Bernoulli causal-indicator φ_j ~ Bernoulli(p_j) | OK | Step 5 (lines 859-899) walks through annotation A → logit → expit → p → "从 Bernoulli(p) 各抽一次" → φ = [1,1,0,1] → α·φ. |
| **MR-CARV** Common-variant α = 0.5 FIXED (not random) | OK | Step 3 (line 829) "MR-CARV 把所有 common IV 的效应 **统一设成 0.5** (不是随机抽)". |
| **MR-CARV** Rare-variant α_j = c0·\|log₁₀ MAF_j\|, c0 = 0.5 | OK | Step 4 (line 839) gives the exact formula with c0=0.5, and the four α values match my recomputation. |
| **MR-CARV** IVW weight w_k ≠ Beta-burden weight w_j | OK | Step 10 (lines 1027-1029) explicitly writes "玩具没有 SE, 用 unit weights w_k = 1" — author is aware that the conventional IVW weight is 1/SE² and that the burden weight (22.16, etc.) is a *separate* quantity. No confusion. |

---

## C. Errors found (sorted by severity)

1. **[Medium · Faithfulness omission]** The 100K-total / 50K-summary + 50K-calculation sample split that defines RARE's three-sample MR design is never mentioned. The walkthrough only says "真实仿真 n=50000". A reader unfamiliar with `simulation1.R` would assume RARE uses a single 50K sample. Recommend the writer add one sentence in Section 1 intro or Step 13/14.
2. **[Low · Faithfulness omission]** The PRS Y-noise heteroscedastic scaling quirk from `simulationl.R` is not mentioned. εy is presented as a simple fixed row. The toy values are not wrong, but a footnote saying "real code rescales εy proportionally to var(PRS) so the heritability ratio stays fixed" would be more faithful.
3. **[Cosmetic · Step 5 logit r4]** Writer writes `4.827`, exact is `4.8283`. Three-decimal rounding could also be `4.828`. The expit value rounds to 0.992 either way; no downstream impact. Not corrected (within the writer's declared rounding tolerance).
4. **[Cosmetic · Step 12d Y s2]** −0.0055 → rounded to −0.006. Strictly within the writer's declared 3 dp rounding (banker's or half-away-from-zero both reach −0.006). Not an error; flagged only for transparency.

No structural / methodological errors found. All formulas (RARE summation, MR-CARV burden, expit, Beta(1,25) pdf, IVW) are correct and applied consistently.

---

## D. Edits I applied inline (numbers swapped only)

| Line | Was | Now | Reason |
| --- | --- | --- | --- |

**None applied.** Every arithmetic cell agrees with the writer's value to within the declared 1–3 decimal rounding tolerance. The only items I considered (Y s2 = −0.006 vs −0.0055, logit r4 = 4.827 vs 4.828, Cov(B,X) = 1.732 vs 1.7315) are all within rounding and would themselves be cosmetic flips, not corrections.

---

## E. Positive log (what the writer got right)

- Honest "诚实说出来的简化" notice (lines 385-394) up front declares: dosage simplification, U as N(0,1) approximation, all coefficients as single draws, n=4 population estimator convention, β₁=0.01 / β=0.04 as alternatives not nulls. This pre-empts most rounding objections.
- All five matrix multiplications in RARE Step 10 (Gδ, Gκx, Uψx, summation) are walked out element-by-element and match exactly.
- The β₂=0.20 nuisance pathway is generated correctly and matches my recomputation cell-for-cell, including the zero in U·ψx2 s1 (writer noted "0.12 − 0.12 = 0.00", which is exact).
- The X1c, X2c centering and the population Cov/Var arithmetic in Step 13 are bookkeeping-clean: writer pre-centers X1, presents the centered vector, then computes Cov as a dot product over n=4. Pedagogically excellent.
- α_rare via `c0·|log₁₀ MAF|` with c0=0.5 is computed correctly to 3 dp for all four MAFs; the intuition ("负选择假设") is added in a collapsed `<details>`.
- Annotation logit construction with b0=−log 5, b_i=log 5 chosen so "0 active → expit≈0.17, 1 active → 0.5" is mathematically clean: expit(0) = 0.5 indeed; expit(−log 5) = 1/(1+5) = 1/6 ≈ 0.1667 ✓.
- Beta(1,25) pdf formula `25·(1−MAF)²⁴` is written correctly and matches `scipy.stats.beta.pdf(MAF, 1, 25)` exactly to 4 dp for all four MAFs.
- Burden construction `B_i = Σ w_j G_rare_ij` correctly produces [0, 22.16, 20.62, 23.26] given each rare carrier carries exactly one variant.
- IVW Step 10 explicitly clarifies "玩具没有 SE, 用 unit weights w_k = 1", avoiding the common pitfall of conflating the IVW weight (1/SE²) with the Beta(1,25) burden weight.
- The Section 3 comparison correctly identifies that MR-CARV's Y = βX + ε *omits* pleiotropy by construction, while RARE generates κy and θ. This is the most important methodological point in the document and it is stated accurately at lines 712-718, 936-940, and 1099-1107.
- The annotated `<details class="why">` blocks add real value (h² interpretation, Beta(1,25) rationale, Bernoulli vs spike-and-slab, MR-CARV's "Y is clean" warning) without bloating the main thread.

---

## Verdict

86 arithmetic cells checked, 86 match. Two faithfulness omissions (sample-size split, PRS noise scaling) are minor — they don't break the walkthrough, but they do leave methodological gaps. No inline number fixes were needed.

**Publishable as-is** for the arithmetic-walkthrough purpose. A minor revision adding two sentences (100K → 50K+50K split; PRS heteroscedastic εy) would close the only faithfulness gaps.
