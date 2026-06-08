# rvMR stubs — implementation summary

5 stubs filled + 1 new function (cell_type_q) for the rvMR package. All
implementations are pure-R, depend only on `stats` (always), `data.table`
(for output formatting; package already imports), and `CompQuadForm` (Suggests
dependency; one function falls back to a moment-matched scaled chi^2 if
CompQuadForm is not installed).

## 1. `iv_partial_r2(b_x, se_x, n)` — Cinelli–Hazlett IV partial R² and RV

**Source spec**: `main.tex` Step 13 eq (38) / (39).

**Formula**:
- $t = b_x / \mathrm{SE}_x$.
- $R^2_{Z \to X} = t^2 / (t^2 + n - 2)$ — partial coefficient of determination of instrument on exposure.
- $RV = (\sqrt{t^2 + 4} - |t|) / 2$ — minimum partial R² a confounder would need with BOTH instrument and outcome to drive the coefficient to zero.
- $RV_{\alpha=0.05}$ — analog at the 5% significance threshold, using $t_{\rm crit} = qt(0.975, n-2)$ via the Cinelli–Hazlett "excess t" form: $RV_\alpha = 0$ if $|t| \le t_{\rm crit}$, else $(\sqrt{(|t| - t_{\rm crit})^2 + 4} - (|t| - t_{\rm crit})) / 2$.

**Sign convention**: implementation uses $|t|$ (not bare $t$ as the spec literally writes) so that RV is a positive partial-$R^2$ regardless of the sign of $b_x$. With bare $t$ and negative $b_x$, the formula gives RV > 1, which is meaningless.

**Citations**: Cinelli C, Hazlett C (2025). *Biometrika* asaf004 (IV setting, primary). Cinelli C, Hazlett C (2020). *JRSS-B* 82(1):39-67 (OVB framework origin).

**Test coverage**: 7 assertions
- Known-value sanity: t=4, n=50002 → $R^2 = 16/50016$.
- RV matches hand-computed formula.
- RV positive for any finite t (both signs).
- Edge case t=0 → RV=1, $R^2$=0.
- Edge case large t → RV→1/t, $R^2$→1.
- RV_alpha=0 when $|t| < t_{\rm crit}$.
- Input validation (NA, zero/negative SE, n ≤ 2).

## 2. `e_value(b_xy, se_xy)` — VanderWeele–Ding E-value for MR

**Source spec**: `main.tex` Step 13 eq (40).

**Formula**:
- $RR_{\rm approx} = \exp(0.91 \cdot \beta_{\rm std})$ (Swanson–VanderWeele 2020 RR-per-SD conversion).
- $E = RR + \sqrt{RR(RR-1)}$ for $RR \ge 1$; for $RR < 1$ use $RR := 1/RR$ (the formula is symmetric in magnitude of departure from 1).
- $E_{\rm CI}$ uses the bound nearest the null: if the 95% CI on $\beta$ crosses 0, $E_{\rm CI} = 1$ (no displacement); else use the bound with smaller magnitude.

**Citations**: VanderWeele TJ, Ding P (2017). *Ann Intern Med* 167(4):268-274 (E-value origin). Swanson SA, VanderWeele TJ (2020). *Epidemiology* 31(3):e23 (MR-specific RR conversion).

**Test coverage**: 7 assertions
- $\beta = 0$ → $RR = 1$, $E = 1$.
- E monotone in $|\beta|$.
- E symmetric in sign of $\beta$.
- Point E matches VanderWeele–Ding hand calc for $RR = 2$: $E = 2 + \sqrt{2} \approx 3.414$.
- CI-E = 1 when CI crosses 0.
- CI-E < point E when CI excludes 0.
- Input validation.

## 3. `annotation_concord(estimates_list, pqtl_anchor)` — cross-mask Cochran-Q

**Source spec**: `main.tex` Step 11 eq (35)-(36).

**Formula**:
- Per class $k$: if `pqtl_anchor` is non-NULL, renormalize $\tilde\beta^{(k)}_{xy} = \beta^{(k)}_{burden,y} / \beta^{(k)}_{burden \to protein}$, with delta-method variance $\mathrm{Var}(\tilde\beta) = \mathrm{SE}_y^2 / \beta_{prot}^2 + \beta_y^2 \mathrm{SE}_{prot}^2 / \beta_{prot}^4$. If `pqtl_anchor` is NULL, $\tilde\beta = \beta_{xy}$ and $\mathrm{Var} = \mathrm{SE}_{xy}^2$.
- Weights $w_k = 1/\mathrm{Var}(\tilde\beta^{(k)})$, pooled $\bar{\tilde\beta} = \sum w_k \tilde\beta^{(k)} / \sum w_k$.
- Cochran $Q = \sum_k w_k (\tilde\beta^{(k)} - \bar{\tilde\beta})^2 \sim \chi^2_{K-1}$.
- Interpretation tier: "underpowered" if any class SE > 5 × |effect|, else "concordant" if $p > 0.05$, else "discordant_investigate".

**Independence assumption**: classes partition variants; class burden IVs are independent under linkage equilibrium across classes. pQTL anchor and outcome estimators are independent (different traits, different samples).

**Citations**: Cochran WG (1954). *Biometrics* 10(1):101-129. Dhindsa RS et al. (2023). *Nature* 622:339-347 (UKB-PPP rare-variant pQTL anchor source). Sun BB et al. (2023). *Nature* 622:329 (UKB-PPP common-variant). Ferkingstad E et al. (2021). *Nat Genet* 53:1712 (deCODE pQTL).

**Test coverage**: 7 assertions
- Concordant classes (no anchor) → $p > 0.5$.
- One class far off → $p < 0.001$ and interpretation = "discordant_investigate".
- Hand-computed Q (1, 2, 3 with SE=1) → $Q = 2$, $\bar = 2$, df=2.
- pQTL-anchor renormalization: delta-method variance matches by hand.
- Requires pQTL fields when anchor is set.
- Rejects K=1 input.
- Rejects zero pQTL anchor (division by zero).

## 4. `cell_type_q(estimates_list, min_donors)` — cross-cell-type Cochran-Q (NEW)

**Source spec**: `main.tex` Step 12 eq (37). New file `R/cell_type_concord.R` (the spec describes this as a wrapper but I implemented as a standalone function for symmetry with `annotation_concord`).

**Formula**:
- $w_c = 1 / \mathrm{Var}(\hat\beta^{(c)}_{xy})$, $\bar\beta = \sum_c w_c \hat\beta^{(c)} / \sum_c w_c$.
- $Q_{\rm cell} = \sum_c w_c (\hat\beta^{(c)} - \bar\beta)^2 \sim \chi^2_{C-1}$.
- `min_donors` filter drops cell types with `n_donors < min_donors`.

**Citations**: Cochran WG (1954). *Biometrics* 10(1):101-129. Cuomo ASE et al. (2025). *medRxiv* 2025.03.20.25324352 (TenK10K Phase 1 28 PBMC cell types). Zhou W, Cuomo ASE et al. (2024). *medRxiv* 2024.05.15.24307317 (SAIGE-QTL substrate). Ray D et al. (2025). *AJHG* 112(7):1597 (first sc-cis-MR comparator at PBMC resolution).

**Test coverage**: 7 assertions
- Concordant cell types → $p > 0.5$.
- One cell type outlier → $p < 0.001$.
- Hand-computed Q (1, 2, 3, 4 with SE=1) → $Q = 5$, df=3.
- `min_donors` filter drops correctly.
- Rejects < 2 cell types (input or post-filter).
- Input validation.

## 5. `heidi_rv(b_xy_per_variant, se_per_variant, weights)` — within-burden generalized-chi^2 heterogeneity test

**Source spec**: `main.tex` Step 10 eq (32)-(34) — **with applied erratum**.

**Erratum applied**: main.tex specifies $T = \delta^\top V_\delta^+ \delta$ AND Davies weights = eigenvalues of $V_\delta$. These are mathematically inconsistent (Mahalanobis form follows $\chi^2_{\rm rank}$; generalized $\chi^2$ with eigenvalue weights is the law of $\delta^\top\delta$). Empirically verified: literal pairing gives ~99% Type-I at nominal 5%; corrected pairing $T = \delta^\top \delta$ gives 4.8% Type-I (2000-rep null simulation, m=4). Implementation uses the corrected pairing; also reports the Mahalanobis sister statistic and its $\chi^2_{m-1}$ p-value for completeness.

**Formula**:
- Contrast $C$ is the (m-1)×m mean-centering operator: $C_{ij} = \delta_{ij} - 1/m$. Full row rank by construction; annihilates the constant vector.
- $\delta = C \hat b_{xy}^{\rm per-var}$, $V_\delta = C \Sigma_b C^\top$, symmetrized numerically.
- If `weights` is a length-m vector: $\Sigma_b = \mathrm{diag}(\mathrm{SE}^2)$ (independence assumed; per-variant SEs absorb scaling).
- If `weights` is m×m: used directly as $\Sigma_b$ (LD-aware path); symmetry + PSD check enforced.
- Eigendecomposition of $V_\delta$. Drop eigenvalues below `tol = max(|eigs|) * m * sqrt(eps)`. `df_effective` = number retained.
- $T = \delta^\top \delta$, p-value = `CompQuadForm::davies(T, lambda = nonzero eigenvalues)$Qq`.
- Fallback if CompQuadForm unavailable: Satterthwaite moment-match scaled $\chi^2$, with warning.
- Numerical guards: clip p to [0, 1]; fall back to Satterthwaite if Davies `ifault != 0`.

**Sister statistic**: $T_{\rm mahalanobis} = \delta^\top V_\delta^+ \delta$ + `pchisq(., df = df_effective, lower.tail = FALSE)`. Reported under `T_mahalanobis` and `p_mahalanobis` for users who want the standard Mahalanobis reading.

**Citations**: Davies RB (1980). *Appl Stat* 29(3):323-333 (generalized $\chi^2$ tail). Kuonen D (1999). *Biometrika* 86(4):929-935 (saddlepoint sister). Zhu Z et al. (2016). *Nat Genet* 48:481-487 (HEIDI origin).

**Test coverage**: 8 assertions
- Concordant inputs → $p = 1$ within numerical tolerance.
- Outlier variant → $p < 0.001$.
- Mahalanobis sister stat returns sensible $\chi^2$ p-value (large under concordance).
- df_effective = m-1 for non-degenerate covariance.
- Covariance matrix input path (LD-aware) returns finite p.
- **Null calibration (2000-rep null simulation, m=4)**: empirical Type-I at nominal 5% in [0.025, 0.075]; mean p-value in [0.40, 0.60].
- Input validation (NA, m<2, negative weights, zero SE).
- `weights_eig` returned matches eigenvalues of `V_delta`.

## Summary table

| Function | File | LOC added | Tests | Citations |
|---|---|---:|---:|---|
| `iv_partial_r2` | sensitivity.R | ~50 | 7 | Cinelli-Hazlett 2025/2020 |
| `e_value` | sensitivity.R | ~55 | 7 | VanderWeele-Ding 2017; Swanson-VanderWeele 2020 |
| `annotation_concord` | annotation_concord.R | ~135 | 7 | Cochran 1954; Dhindsa 2023 |
| `cell_type_q` | cell_type_concord.R (NEW) | ~190 | 7 | Cochran 1954; Cuomo 2025; Zhou-Cuomo 2024; Ray 2025 |
| `heidi_rv` | heidi_rv.R | ~135 | 8 | Davies 1980; Kuonen 1999; Zhu 2016 |

Total: ~565 LOC + 78 new test assertions across 4 test files.
Baseline: 82 PASS / 0 FAIL → Now: 160 PASS / 0 FAIL.
