# CHECKER_RUBRIC.md — Ground-truth Rubric for rvSMR Algorithm Writeup

*Independent checker pass. Built from HANDOVER_2026-05-27.md §10, briefing_for_wei.md §3–§4, citation_audit_2026-05-27.md, and rvMR/R/*.R sources. Do NOT trust draft against this rubric blindly — also weigh §10 Task A spec exactly: "assumption → derivation → estimand → finite-sample concern" structure per section.*

---

## §1. Notation & Setup

### Required assumptions
- Two-sample summary-stat MR design (default; one-sample / overlap as an extension via $R_{xy}$).
- Per (gene $g$ × mask $k$ × cell type $c$) tuple is the unit of analysis.
- Exposure-side sample size $n_x$; outcome-side sample size $n_y$.
- $R_{xy}$ across-sample correlation block introduced here (default 0 = two-sample independence).
- Linear weighted-sum burden as the instrument variable (NOT SKAT / ACAT quadratic forms — explicit note).
- Per-allele effect scale convention: $\hat\beta_x$, $\hat\beta_y$ are on the **per-allele effect of the weighted burden** $B = \sum_j w_j G_j$, NOT rescaled by allele count / MAF / carrier frequency.

### Required formulas
- Burden definition: $B = \sum_{j \in \mathrm{mask}} w_j G_j$.
- Summary-stat triplet per instrument: $(\hat b_{\mathrm{burden},x}, \mathrm{SE}_x, n_x)$ for exposure; $(\hat b_{\mathrm{burden},y}, \mathrm{SE}_y, n_y)$ for outcome.
- Sample-overlap correlation: $\rho = \operatorname{cor}(\hat\beta_x, \hat\beta_y)$ (scalar K=1) or block matrix $R_{xy} \in \mathbb{R}^{K\times K}$ (K≥2).
- Per-(gene × mask × cell type) indexing: $(g, k, c)$ with $K \ge 3$ masks per gene committed (pLoF, missense:LC, regulatory).

### Required theorem / proposition statements
- None at this level — definitional only. Must state the substrate (SAIGE-QTL on exposure; Genebass / RGC-ME / UKB-PPP on outcome) is exogenous to the math but defines the inputs.

### Required citations
- Pierce BL, Burgess S. 2013. *AJE* 178(7):1177–1184. doi:10.1093/aje/kwt084 — two-sample IV setup.
- Burgess S, Butterworth A, Thompson SG. 2013. *Genetic Epidemiology* 37(7):658–665. doi:10.1002/gepi.21758 — summary-stat MR setup. (NOTE: journal is *Genet Epidemiol*, NOT *Stat Med*.)
- Zhou W, et al. 2022. *Nat Genet* 54(10):1466–1469. doi:10.1038/s41588-022-01178-w — SAIGE-GENE+ (exposure-side substrate, Step 2 burden score test).
- Zhou W, Cuomo ASE et al. 2024. medRxiv 2024.05.15.24307317 — SAIGE-QTL.
- Cuomo et al. 2025. medRxiv 2025.03.20.25324352 — TenK10K Phase 1.

### Common pitfalls / things that would be WRONG
- Writing $\hat\beta_x$ as "per-variant" rather than per-allele burden coefficient.
- Implying SKAT / ACAT can be inverted to a Wald ratio (they're quadratic forms with no signed effect).
- Claiming RareEffect Step 4 BLUP / PEV can substitute for SAIGE-GENE+ Step 2 burden score test (BLUP is Bayesian and breaks the chi-square reference).
- Writing "Sun BB" for RGC-ME (correct is Sun KY).

---

## §2. Burden Construction

### Required assumptions
- **Linear** weighted-sum (instrument linear in genotype) — required for Wald-ratio interpretation.
- Weights $w_j$ are fixed (annotation-driven, not data-driven inside the sample).
- STAAR-style weights default: Beta($\mathrm{MAF}; 1, 25$); explicit warning that STAAR-O is a CCT p-value (need to extract STAAR-B sub-statistic).
- SKAT-style quadratic statistics and ACAT excluded as instruments.
- Cauchy combination test (CCT) is NOT applied across K masks within a gene (would defeat AR over-id); CCT only across orthogonal axes (across-gene / across-trait / across-cell-type) or across MAF-weight nuisance settings within a single mask.

### Required formulas
- $B = \sum_{j \in \mathrm{mask}} w_j G_j$ (Madsen–Browning 2009 weighted-sum statistic).
- Default weights: $w_j \propto \mathrm{Beta}(p_j; 1, 25)$ (Wu et al. 2011 SKAT weighting; $p_j$ = MAF of variant $j$).
- (Optional) Madsen–Browning original weighting: $w_j = 1 / \sqrt{p_j (1-p_j)}$.

### Required theorem / proposition statements
- "Only linear burden statistics preserve Wald-ratio interpretation" — explicit.
- STAAR-B (linear, signed) extracted from STAAR-O (Cauchy-combined p-value).

### Required citations
- Madsen BE, Browning SR. 2009. *PLOS Genet* 5(2):e1000384. doi:10.1371/journal.pgen.1000384 — **primary** for $B = \sum_j w_j G_j$ burden weighted-sum.
- Wu MC, Lee S, Cai T, Li Y, Boehnke M, Lin X. 2011. *AJHG* 89(1):82–93. doi:10.1016/j.ajhg.2011.05.029 — SKAT, source of Beta($\mathrm{MAF}; 1, 25$) weighting.
- Zhou W, et al. 2022. *Nat Genet* 54(10):1466–1469 — SAIGE-GENE+ Step 2 burden score test.
- Li X, et al. 2020. *Nat Genet* 52:969–983. doi:10.1038/s41588-020-0676-4 — STAAR annotation weights. **(Default per user instruction: 2020, NOT 2022 Li Z STAARpipeline.)**
- Liu Y, Chen S, Li Z, et al. ACAT 2019 *AJHG* — only if explicitly named.
- **Open citation decision (b)**: if the draft cites CAST (Morgenthaler-Thilly 2007 *Mut Res* 615:28–56) and/or CMC (Li-Leal 2008 *AJHG* 83:311–321), flag as user-decision item.

### Common pitfalls / things that would be WRONG
- Saying CMC = Morgenthaler-Thilly 2007 (wrong; that's CAST). CMC is Li-Leal 2008 *AJHG* 83:311.
- Calling STAAR-O the burden statistic (it's a CCT p-value).
- Citing "STAAR 2022 *Nat Methods*" for annotation weights (that's the pipeline paper, first author Li Z — not Li X).
- Treating SKAT as the burden — SKAT is quadratic; cannot be a signed IV.
- Asserting Madsen-Browning without a citation (was a previously-flagged gap; the audit landed the fix).

---

## §3. Identification under Additive Linear SMM

### Required assumptions (NON-NEGOTIABLE — all six must appear explicitly)
1. **Relevance**: $\alpha_j \neq 0$ for at least one $j \in \mathrm{mask}$ (i.e., $E[X \mid Z = z]$ depends on $z$).
2. **Exchangeability** (independence): $Z \perp U$ — instrument independent of unmeasured confounders $U$ of $X \to Y$.
3. **Exclusion restriction**: $Z \to Y$ only via $X$ (no direct effect, no pleiotropy).
4. **Monotonicity / no-defier-within-mask** (Imbens-Angrist adaptation): every variant in the mask moves expression in the same sign on average; no defier alleles.
5. **Linearity**: structural equation $E[Y(x) - Y(0)] = \beta \cdot x$ (additive linear SMM).
6. **No Z × X interaction**: Robins 1994 no-interaction structural nested mean model condition; the treatment effect does not depend on the level of the instrument.

### Required formulas (with exact LaTeX-ready form)
- Carrier-frequency-weighted estimand (Delta 2 of briefing §3) — IVW form, the one rvSMR adopts:
$$\beta_{\mathrm{burden}} \;=\; \frac{\sum_j w_j^{\,2}\, p_j(1-p_j)\, \alpha_j \gamma_j}{\sum_j w_j^{\,2}\, p_j(1-p_j)\, \alpha_j^{\,2}} \;=\; \sum_j \pi_j \cdot \frac{\gamma_j}{\alpha_j},$$
with weights summing to 1:
$$\pi_j \;\propto\; w_j^{\,2}\, p_j(1-p_j)\, \alpha_j^{\,2}, \qquad \sum_j \pi_j = 1.$$
- (Optional sibling form, single-pre-specified-burden) First-power weights for the $Z = \sum_j w_j G_j$ alternative:
$$\pi_j^{(1)} \;\propto\; w_j\, p_j(1-p_j)\, \alpha_j.$$
- Population Wald ratio: $\beta_{\mathrm{burden}} = \mathrm{Cov}(Z,Y) / \mathrm{Cov}(Z,X)$, which under the assumptions equals the formula above.

### Required theorem / proposition statements
- **Identification theorem**: "Under (1)–(6), the burden-IV Wald ratio $\hat\beta_y / \hat\beta_x$ identifies the carrier-frequency-weighted average of per-variant local causal slopes $\sum_j \pi_j (\gamma_j / \alpha_j)$ with IVW weights $\pi_j \propto w_j^2 p_j(1-p_j) \alpha_j^2$." (This is the Delta 2 commitment.)
- "Under homoskedasticity, IVW weights $\pi_j \propto w_j^2 p_j(1-p_j) \alpha_j^2$ are the 2SLS-optimal weights" — explicit.
- "The single-pre-specified-burden alternative identifies $\sum_j \pi_j^{(1)} (\gamma_j / \alpha_j)$ with $\pi_j^{(1)} \propto w_j p_j(1-p_j) \alpha_j$; rvSMR adopts the IVW (squared-weight) form because it admits the formal AR weak-IV-robust inference machinery." — explicit.

### Required citations
- Robins JM. 1994. "Correcting for non-compliance in randomized trials using structural nested mean models." *Communications in Statistics — Theory and Methods* 23(8):2379–2412. doi:10.1080/03610929408831393 — additive SMM, no-Z×X-interaction condition. **NOT *Biometrics*.**
- Didelez V, Sheehan N. 2007. "Mendelian randomization as an instrumental variable approach to causal inference." *Statistical Methods in Medical Research* 16(4):309–330. doi:10.1177/0962280206077743 — IV in epidemiology framing.
- Imbens GW, Angrist JD. 1994. "Identification and estimation of local average treatment effects." *Econometrica* 62(2):467–475. doi:10.2307/2951620 — LATE / monotonicity.
- Burgess S, Labrecque JA. 2018. *European Journal of Epidemiology* 33(10):947–952. doi:10.1007/s10654-018-0424-6 — interpretation of MR causal estimates.
- Madsen BE, Browning SR. 2009 (above) — burden weighted-sum primary cite.

### Common pitfalls / things that would be WRONG
- Citing Robins 1994 to *Biometrics* — wrong, it's *Comm Stat Theory Methods* 23(8):2379, doi:10.1080/03610929408831393. **Trap door** flagged by audit.
- Stating the estimand without IVW weights — must show $\pi_j \propto w_j^2 p_j(1-p_j) \alpha_j^2$, not unweighted average.
- Confusing the squared-weight IVW form with the first-power form (the latter is the $Z = \sum_j w_j G_j$ single-burden case; rvSMR adopts the IVW squared-weight form).
- Implying $\pi_j$ are pre-specified annotation weights (they are NOT — they're mechanical from IV identification, contrasting with MR-CARV pre-fixed annotation weights and RARE selection-coefficient prior).
- Failing to state monotonicity adaptation explicitly (no-defier-within-mask is a non-default Imbens-Angrist condition).
- Omitting "no Z × X interaction" — this is the Robins 1994 no-interaction SMM condition; without it the estimand isn't identified.

---

## §4. Inference K = 1 (Closed-form Anderson-Rubin)

### Required assumptions
- All §3 assumptions (identification).
- Single instrument $K = 1$; AR test reduces algebraically to Fieller's 1954 confidence set for a ratio.
- Default $\rho = 0$ (two-sample non-overlap); $\rho \neq 0$ for sample-overlap correction.
- $\chi^2_{1, 1-\alpha}$ asymptotic null (default $\alpha = 0.05$, $c = 3.841$).
- Delta-method variance is **reported for transparency only**, NOT used for CI in the weak-IV regime.

### Required formulas (exact LaTeX-ready)
- **Wald ratio point estimate**: $\hat b_{xy} = \hat b_y / \hat b_x$.
- **Delta-method SE** (first-order): 
$$\operatorname{Var}(\hat b_y / \hat b_x) \;\approx\; \frac{\mathrm{SE}_y^{\,2}}{\hat b_x^{\,2}} + \frac{\hat b_y^{\,2}\, \mathrm{SE}_x^{\,2}}{\hat b_x^{\,4}} - 2\, \frac{\hat b_y}{\hat b_x^{\,3}} \cdot \mathrm{cov}(\hat b_x, \hat b_y),$$
with $\mathrm{cov}(\hat b_x, \hat b_y) = \rho \, \mathrm{SE}_x \, \mathrm{SE}_y$.
- **First-stage F**: $F = (\hat b_x / \mathrm{SE}_x)^2$.
- **AR statistic K=1**:
$$AR(\beta_0) \;=\; \frac{(\hat b_y - \beta_0\, \hat b_x)^2}{\mathrm{SE}_y^{\,2} + \beta_0^{\,2}\, \mathrm{SE}_x^{\,2} - 2\, \beta_0\, \rho\, \mathrm{SE}_x\, \mathrm{SE}_y}.$$
- **AR confidence set** (level set form):
$$\mathcal{C}_{1-\alpha} \;=\; \{\beta_0 : AR(\beta_0) \le \chi^2_{1, 1-\alpha}\}.$$
- **Quadratic inequality** (with $c = \chi^2_{1, 1-\alpha} = 3.841$):
$$\underbrace{(\hat b_x^{\,2} - c\, \mathrm{SE}_x^{\,2})}_{A}\, \beta_0^{\,2} - 2\underbrace{(\hat b_x\, \hat b_y - c\, \rho\, \mathrm{SE}_x\, \mathrm{SE}_y)}_{B}\, \beta_0 + \underbrace{(\hat b_y^{\,2} - c\, \mathrm{SE}_y^{\,2})}_{C} \;\le\; 0.$$
- Discriminant: $\Delta = B^2 - 4AC$ (in the briefing this is written as $\Delta = B^2 - AC$ with a different normalization of $B$; either form is acceptable provided the corresponding $B$ definition is consistent — check the draft uses one form consistently. The R code uses $B = -2(\hat b_x \hat b_y - c\rho \mathrm{SE}_x \mathrm{SE}_y)$ as defined here with $\Delta = B^2 - 4AC$).
- Sign of $A$: $A > 0 \iff F > c$ (i.e., $A$ positive iff first-stage F exceeds the critical value).

### Required theorem / proposition statements — the four-CI-shape table (NON-NEGOTIABLE)
The writeup MUST contain a table with the four CI shapes. Reference table from briefing §4.3:

| Case | Sign of $A = \hat b_x^2 - c \mathrm{SE}_x^2$ | Discriminant $\Delta$ | CI shape |
|------|----------------------------------------------|-----------------------|----------|
| Strong IV, identified | $A > 0$ ($F > c$) | $\Delta > 0$ | **Bounded interval** $[a, b]$ |
| Strong IV, no $\beta_0$ accepted | $A > 0$ ($F > c$) | $\Delta < 0$ | **Empty** (numerical artifact, prob $\le \alpha$ under null) |
| Weak IV | $A < 0$ ($F < c$) | $\Delta > 0$ | **Disconnected union** $(-\infty, a] \cup [b, \infty)$ |
| Weak IV, no information | $A < 0$ ($F < c$) | $\Delta < 0$ | **Whole real line** |

Must explicitly state:
- "Disconnected and whole-line cases are not pathologies of AR — they are the honest geometric answer when data carry too little signal."
- Caveat: at K=1 the AR reduces to inverted Fieller (Fieller 1954 *JRSSB* 16:175); over-identifying power emerges only at K≥2.
- "Coverage stays at the nominal level $1-\alpha$ uniformly over instrument strength — including $F < 1$ — because the test pivots on the reduced-form moment condition, not the first-stage estimator."

### Required citations
- Anderson TW, Rubin H. 1949. *Annals of Mathematical Statistics* 20(1):46–63. doi:10.1214/aoms/1177730090.
- Wang S, Kang H. 2022. *Biometrics* 78(4):1699–1713 — AR for MR.
- Fieller EC. 1954. "Some problems in interval estimation." *JRSSB* 16(2):175–185.
- Lee DS, McCrary J, Moreira MJ, Porter J. 2022. *American Economic Review* 112(10):3260–3290. doi:10.1257/aer.20211063 — tF critical values; rejection of Stock-Yogo $F > 10$ filter.
- Stock JH, Yogo M. 2005 — referenced as the "$F > 10$" filter that rvSMR explicitly does NOT use.
- **Open citation decision (a)**: delta-method variance source — default to flagging as TODO unless the writer picks Burgess-Thompson 2017 textbook / Thomas et al. 2007 / Rothman-Greenland. The Bowden-Vansteelandt 2011 attribution is WRONG (that paper is case-control SMM).

### Common pitfalls / things that would be WRONG
- Combining the AR CI with a Stock-Yogo $F > 10$ filter (do NOT — explicit in roxygen).
- Citing Bowden-Vansteelandt 2011 *Stat Med* 30:678 for delta-method variance (wrong; case-control SMM paper).
- Asserting the AR CI has K=1 over-id power (it doesn't; AR at K=1 ≡ inverted Fieller).
- Writing $A = SE_x^2 (F - c)$ as if exact equality without the carrier "up to positive scalar"; the relation is $A = SE_x^2 \cdot (\hat b_x^2 / SE_x^2 - c) = SE_x^2 (F - c)$, which IS exact for K=1, so the briefing's "up to a positive scalar" is loose phrasing — either form is fine in the draft.
- Forgetting to state $\rho = 0$ default for two-sample.
- Saying the "empty" case is impossible (it can occur with probability $\le \alpha$ under the null; should be flagged as numerical artifact).

---

## §5. Inference K ≥ 2 (Multi-instrument AR + Sargan-J)

### Required assumptions
- All §3 assumptions per instrument.
- $K \ge 2$ instruments stacked (rvSMR commits to $K \ge 3$ masks per gene = pLoF, missense:LC, regulatory).
- Between-instrument correlations $R_{xx}, R_{yy}$ (default $I_K$ = independence under linkage equilibrium between masks).
- Across-sample correlation block $R_{xy}$ (default $\mathbf{0}_{K \times K}$ for two-sample).
- $\chi^2_K$ asymptotic null under $H_0: \beta = \beta_0$.
- $\chi^2_{K-1}$ for the Sargan-J over-identification statistic at the AR argmin.

### Required formulas (exact LaTeX-ready)
- Moment vector: 
$$m(\beta_0) \;=\; \hat{\boldsymbol{b}}_y - \beta_0\, \hat{\boldsymbol{b}}_x \;\in\; \mathbb{R}^K.$$
- Covariance matrix at $\beta_0$:
$$V(\beta_0) \;=\; D_y\, R_{yy}\, D_y + \beta_0^{\,2}\, D_x\, R_{xx}\, D_x - 2\, \beta_0\, D_y\, R_{xy}\, D_x,$$
with $D_x = \mathrm{diag}(\mathrm{SE}_x)$, $D_y = \mathrm{diag}(\mathrm{SE}_y)$.
- AR statistic:
$$AR(\beta_0) \;=\; m(\beta_0)^\top\, V(\beta_0)^{-1}\, m(\beta_0) \;\sim\; \chi^2_K \quad \text{under } H_0: \beta = \beta_0.$$
- AR confidence set: 
$$\mathcal{C}_{1-\alpha} \;=\; \{\beta_0 : AR(\beta_0) \le \chi^2_{K, 1-\alpha}\}.$$
- **Sargan / J over-id statistic**: 
$$J \;=\; AR(\hat\beta_{\mathrm{AR}}) \;=\; \min_{\beta_0}\, AR(\beta_0) \;\sim\; \chi^2_{K-1} \quad \text{under joint validity}.$$
- $\hat\beta_{\mathrm{AR}} = \arg\min_{\beta_0} AR(\beta_0)$ (point estimate).
- $J$ p-value: $1 - F_{\chi^2_{K-1}}(J)$.

### Required theorem / proposition statements
- "Under $H_0: \beta = \beta_0$ and IV validity, $AR(\beta_0) \sim \chi^2_K$ asymptotically."
- "Under joint IV validity (no horizontal pleiotropy), $J = AR(\hat\beta) \sim \chi^2_{K-1}$."
- "The four CI shapes (bounded_interval, disconnected_union, whole_line, empty) generalize from the K=1 closed form to K≥2 via grid + uniroot inversion."
- "Each per-mask first-stage F is strictly weaker than the F of a pooled burden; AR absorbs the hit and reports an honestly-covered wider set rather than burying the weakness." — explicit trade-off.

### Required citations
- Anderson-Rubin 1949 (as above).
- Wang-Kang 2022 (as above).
- Patel A, Lane J, Burgess S. 2024. arXiv:2408.09868 — MVMR-AR extension; primary K≥2 reference.
- Sargan JD. 1958. *Econometrica* 26(3):393–415 — original J-test.
- Hansen LP. 1982. *Econometrica* 50(4):1029–1054 — GMM over-identification generalization.

### Common pitfalls / things that would be WRONG
- Writing $V(\beta_0) = D_y R_{yy} D_y + \beta_0^2 D_x R_{xx} D_x$ without the cross-term $-2\beta_0 D_y R_{xy} D_x$ — incorrect; the cross-term is required.
- Sign error on the cross-term: the cross-term is $-2\beta_0\, D_y R_{xy} D_x$ (NOT $+2\beta_0$); the sign flips with $\beta_0$.
- Stating $J \sim \chi^2_K$ instead of $\chi^2_{K-1}$ at the argmin (df = $K-1$, NOT $K$, because one df is consumed by the point estimate).
- Applying CCT across the K masks (defeats AR over-id; explicit prohibition).
- Calling the AR CI types `bounded` / `disconnected` / `whole_line` / `empty` instead of `bounded_interval` / `disconnected_union` / `whole_line` / `empty` if matching the R code conventions for K≥2 (minor cosmetic — flag if inconsistent with K=1 type labels).
- Writing $V_{xy}$ as $D_x R_{xy} D_y$ — wrong; correct is $V_{xy} = D_y R_{xy} D_x$ (verify against `mrAR_multi.R:157`).

---

## §6. Sample-overlap Correction

### Required assumptions
- $R_{xy}$ is the across-sample correlation block; non-zero for one-sample / overlapping-sample MR.
- The exposure and outcome estimators share individuals → correlation is induced even at independent loci.

### Required formulas
- K=1 cross term in AR denominator: $-2 \beta_0 \rho \mathrm{SE}_x \mathrm{SE}_y$.
- K≥2 cross block: $V(\beta_0)$ contains $-2 \beta_0\, D_y R_{xy} D_x$.
- Sign-flip property: the cross-term changes sign with $\beta_0$ (highlight: as $\beta_0$ moves through 0, the cross-correction sign flips).

### Required theorem / proposition statements
- "Sample-overlap correlation enters AR via the cross-term in $V(\beta_0)$ whose sign flips with $\beta_0$."
- "Default in two-sample non-overlapping MR: $\rho = 0$ (K=1) / $R_{xy} = \mathbf{0}_{K\times K}$ (K≥2)."

### Required citations
- Burgess S, Davies NM, Thompson SG. 2016. *Genet Epidemiol* — sample-overlap correction in MR (as referenced in mrAR.R roxygen).

### Common pitfalls / things that would be WRONG
- Assuming $\rho = 0$ for one-sample analyses.
- Missing the sign-flip discussion (the cross-correction's geometric meaning).

---

## §7. HEIDI-rv Math (Within-burden LOO heterogeneity)

### Required assumptions
- Single burden of size $m$ variants.
- Per-variant burden eQTL effects $(b_{x,j}, \mathrm{SE}_{x,j})$ AND per-variant burden GWAS effects $(b_{y,j}, \mathrm{SE}_{y,j})$ are observable (or equivalently the per-variant score statistics + LD matrix $\Sigma_G$). SAIGE-QTL / Genebass must be re-run in single-variant mode.
- $\Sigma_{\hat b}$ = covariance matrix of the per-variant Wald ratios (or LOO burden Wald ratios).
- Possibility of rank-deficient $V_\delta$ when monomorphic-in-sample variants are present.

### Required formulas
- Per-variant Wald ratio: $\hat b_{xy,j} = b_{y,j} / b_{x,j}$ (or LOO version: $\hat b_{xy}^{(B-j)}$, the Wald ratio with variant $j$ excluded from the burden).
- Deviation vector: $d_j = \hat b_{xy}^{(B-j)} - \hat b_{xy}^{(B)}$.
- Contrast operator: $C$ is full-row-rank $(m-1) \times m$ (any mean-centering contrast).
- $\delta = C\, \hat{\boldsymbol b}_{xy}^{\mathrm{per-var}} \in \mathbb{R}^{m-1}$.
- $V_\delta = C\, \Sigma_{\hat b}\, C^\top$.
- **Test statistic**: 
$$T \;=\; \delta^\top\, V_\delta^{+}\, \delta,$$
where $V_\delta^{+}$ is the Moore-Penrose pseudoinverse.
- **Null distribution**: generalized $\chi^2$ — weighted sum of independent $\chi^2_1$ with weights equal to the **non-zero eigenvalues of $V_\delta$ itself** (not of $V_\delta^+ V_\delta$).
- Tail probability via Davies algorithm: `CompQuadForm::davies(q = T, lambda = eig_nz(V_delta))`.

### Required theorem / proposition statements
- "Under $H_0: \delta = 0$, $T = \delta^\top V_\delta^+ \delta$ follows a generalized $\chi^2$ whose mixing weights are the non-zero eigenvalues of $V_\delta$."
- **NON-NEGOTIABLE proof / justification**: "Using eigenvalues of $V_\delta^+ V_\delta$ (the orthogonal projector onto $\mathrm{range}(V_\delta)$) instead would yield weights $\{1, \ldots, 1, 0, \ldots, 0\}$, collapsing the generalized $\chi^2$ to plain $\chi^2_{\mathrm{rank}(V_\delta)} = \chi^2_{m-1}$ — anti-conservative." (This is the headline correctness claim of HEIDI-rv; the writeup must state it.)
- "Power scales as $\mathcal{O}(1/m)$ — single pleiotropic variant contributes $1/m$ to the burden Wald; HEIDI-rv is a sensitivity diagnostic, not a primary pleiotropy test."

### Required citations
- Zhu Z, Zhang F, Hu H, Bakshi A, Robinson MR, Powell JE, et al. 2016. *Nature Genetics* 48(5):481–487. doi:10.1038/ng.3538 — HEIDI original.
- Davies RB. 1980. "Algorithm AS 155." *Applied Statistics (JRSS Series C)* 29(3):323–333. doi:10.2307/2346911 — generalized $\chi^2$ tail.
- Kuonen D. 1999. *Biometrika* 86(4):929–935. doi:10.1093/biomet/86.4.929 — saddlepoint approximation for quadratic forms.

### Common pitfalls / things that would be WRONG
- Using eigenvalues of $V_\delta^+ V_\delta$ as Davies weights — **wrong**, collapses to $\chi^2_{m-1}$ (anti-conservative). Headline trap.
- Calling HEIDI-rv a "primary pleiotropy test" — wrong; it's a sensitivity diagnostic with $\mathcal{O}(1/m)$ insensitivity.
- Implying STAARpipeline aggregate output is sufficient (it is NOT — need per-variant score statistics).
- Forgetting to drop zero eigenvalues before passing to Davies (preserves effective df under rank deficiency).
- Citing Zhu 2016 to the wrong journal/volume (verified: *Nat Genet* 48:481).

---

## §8. Annotation-class Concordance

### Required assumptions
- Variants partitioned into $K$ functional annotation classes (default 3: LoF / missense:LC / regulatory).
- Within-class burden IV is independent of between-class burden IV under linkage equilibrium between classes (typical for rare-variant aggregates).
- pQTL-anchor mediator-scale normalization preprocessing — required because raw class-specific $b_{\mathrm{burden},y}$ are NOT directly comparable (each class hits the gene with a different effect on protein abundance).
- Sufficient power per class.

### Required formulas
- pQTL-anchor normalization:
$$\tilde\beta^{\mathrm{class}}_{xy} \;=\; \frac{\beta^{\mathrm{class}}_{\mathrm{burden},y}}{\beta^{\mathrm{class}}_{\mathrm{burden} \to \mathrm{protein}}}.$$
- Variance propagation (delta method) on the renormalized class-specific Wald ratio.
- Cochran's Q across $K$ classes:
$$Q \;=\; \sum_{k=1}^{K} w_k\, (\tilde\beta_k - \bar{\tilde\beta})^2 \;\sim\; \chi^2_{K-1},$$
with $w_k = 1 / \mathrm{Var}(\tilde\beta_k)$, $\bar{\tilde\beta} = \sum_k w_k \tilde\beta_k / \sum_k w_k$.
- p-value: $1 - F_{\chi^2_{K-1}}(Q)$.

### Required theorem / proposition statements
- "After pQTL-anchor normalization, class-specific Wald ratios are on a common mediator (protein-abundance) scale and Cochran's Q is the natural homogeneity test with df = $K-1$."
- "Discordance triggers mechanism investigation (dominant-negative, gain-of-function, cell-type-specific regulation), NOT automatic rejection."

### Required citations
- Cochran WG. 1954. "The combination of estimates from different experiments." *Biometrics* 10(1):101–129. doi:10.2307/3001666.
- pQTL anchor data sources:
  - Sun BB, Maranville JC, Peters JE, Stacey D, Staley JR, et al. 2018. *Nature* 558(7708):73–79. doi:10.1038/s41586-018-0175-1 (INTERVAL pQTL).
  - Sun BB, Chiou J, Traylor M, Benner C, Hsu YH, et al. 2023. *Nature* 622(7982):329–338. doi:10.1038/s41586-023-06592-6 (UKB-PPP common-variant).
  - Ferkingstad E, et al. 2021. *Nat Genet* 53(12):1712–1721. doi:10.1038/s41588-021-00978-w (deCODE pQTL).
  - Dhindsa RS, et al. 2023. *Nature* 622(7982):339–347. doi:10.1038/s41586-023-06547-x (UKB-PPP rare-variant — for rvSMR-relevant pQTL anchor).
- **Open citation decision (d)**: sign-concordance reference. If the draft cites Han-Eskin 2011 *AJHG* 88:586 for sign concordance — **WRONG** (that's random-effects meta-analysis). Acceptable alternates: Owen 2009 / Whitlock 2005 (Stouffer-style p-combiners) — FLAG AS TODO unless the writer explicitly justifies.

### Common pitfalls / things that would be WRONG
- Skipping pQTL-anchor normalization (raw class-specific $b_{\mathrm{burden},y}$ are NOT directly comparable).
- Citing Han-Eskin 2011 *AJHG* 88:586 for directional sign concordance — wrong (random-effects meta-analysis, not sign-concordance).
- Citing Sun BB for both INTERVAL and UKB-PPP without distinguishing (different papers, different years).
- Interpreting discordance as automatic rejection (it triggers mechanism investigation).

---

## §9. Cell-type Concordance

### Required assumptions
- Per-(gene × mask) burden IV estimated per cell type $c$.
- Cross-cell-type homogeneity test analogous to §8 but on the cell-type axis.
- Cell-type stratification multiplies $K$ to $3 \cdot c$ in the joint AR formulation.

### Required formulas
- Analog of §8 Cochran's Q across cell types:
$$Q_{\mathrm{cell}} \;=\; \sum_{c} w_c\, (\hat\beta_{xy}^{(c)} - \bar\beta)^2 \;\sim\; \chi^2_{C-1},$$
with $C$ = number of cell types, $w_c = 1/\mathrm{Var}(\hat\beta_{xy}^{(c)})$.
- Joint stratified AR formulation: total instruments = $K \cdot C$ = (masks) × (cell types).

### Required theorem / proposition statements
- "Cross-cell-type homogeneity is the third over-identification axis (alongside HEIDI-rv within-burden and annotation-class)."
- "Cell-type discordance can reflect tissue-specific regulation (true biology) or cell-type-specific pleiotropy (artifact); requires mechanism follow-up."

### Required citations
- Yazar S, et al. 2022. *Science* 376(6589):eabf3041. doi:10.1126/science.abf3041 — OneK1K.
- Cuomo et al. 2025. medRxiv 2025.03.20.25324352 — TenK10K Phase 1 (28 PBMC cell types).
- Ray et al. 2025. *AJHG* 112(7):1597 — sc-cis-MR comparator (first cell-type cis-MR; 14 immune types). **Note: first author is Ray, NOT Ge — Ge 2025 was a confabulation.**

### Common pitfalls / things that would be WRONG
- Saying "Ge et al. 2025" instead of "Ray et al. 2025 *AJHG* 112(7):1597" — confabulation trap.
- Asserting hepatocyte resolution from TenK10K Phase 1 (PBMC-only); hepatocyte awaits Phase 2.
- Overstating cell-type test power — single-cell burden estimates can be noisy.

---

## §10. Sensitivity (Cinelli-Hazlett + Swanson-VanderWeele)

### Required assumptions
- Two-sample summary-stat adaptation of Cinelli-Hazlett (originally one-sample regression).
- "Robustness value" interpretation: minimum partial $R^2$ that an unobserved violator of exclusion would need with **both** the instrument and the outcome to nullify the effect.
- E-value approximation: continuous outcome with standardized $\beta$; uses Swanson-VanderWeele 2020 approximation $RR \approx \exp(0.91 \beta_{\mathrm{std}})$.

### Required formulas (exact LaTeX-ready)
- **IV partial $R^2$** (two-sample summary form):
$$R^2_{Z \to X} \;\approx\; \frac{t^2}{t^2 + n - 2}, \qquad t = \hat b_x / \mathrm{SE}_x, \quad n = n_x.$$
- **Cinelli-Hazlett robustness value**:
$$RV \;=\; \frac{\sqrt{t^2 + 4} - t}{2}.$$
- **Risk-ratio approximation** (Swanson-VanderWeele 2020):
$$RR \;\approx\; \exp(0.91 \cdot \beta_{\mathrm{std}}).$$
- **E-value**:
$$E \;=\; RR + \sqrt{RR \cdot (RR - 1)}.$$
- Report both the point E-value and the CI-bound E-value (the latter at the CI bound nearest the null is more conservative — the "value typically reported").

### Required theorem / proposition statements
- "$RV$ is invariant to scale of the instrument."
- "An unobserved confounder with partial $R^2$ on both the instrument and the outcome below $RV$ cannot nullify the effect."
- "E-value reports the minimum unmeasured-confounder-association strength (RR scale) needed to fully explain away the observed $b_{xy}$."

### Required citations
- **Primary** (IV-specific): Cinelli C, Hazlett C. 2025. "An omitted variable bias framework for sensitivity analysis of instrumental variables." *Biometrika* asaf004. doi:10.1093/biomet/asaf004.
- **Secondary** (OLS origin): Cinelli C, Hazlett C. 2020. "Making Sense of Sensitivity: Extending Omitted Variable Bias." *JRSS-B* 82(1):39–67. doi:10.1111/rssb.12348.
- Swanson SA, VanderWeele TJ. 2020. "E-Values for Mendelian Randomization." *Epidemiology* 31(3):e23–e24. doi:10.1097/EDE.0000000000001164. PMID:31996542. (Letter format — note e-pages.)
- VanderWeele TJ, Ding P. 2017. *Annals of Internal Medicine* 167(4):268–274. doi:10.7326/M16-2607 — E-value origin.
- Wang L, Tchetgen Tchetgen EJ. 2018. *JRSS-B* 80(3):531–550. doi:10.1111/rssb.12262 — bias bounds under invalid IV (for §11 limitations linkage to coherent pleiotropy).

### Common pitfalls / things that would be WRONG
- Citing **only** Cinelli-Hazlett 2020 *JRSS-B* without the 2025 *Biometrika* IV-specific primary — wrong; both required.
- Citing Cinelli-Hazlett 2020 as "*Biometrika*" — wrong, it's *JRSS-B* 82(1):39.
- Using the OLS partial-$R^2$ formula without the two-sample summary-stat adaptation $t^2 / (t^2 + n - 2)$.
- Writing $RV = (t - \sqrt{t^2 + 4})/2$ — sign error; the correct form is $RV = (\sqrt{t^2 + 4} - t)/2$ (positive scalar).
- Citing the E-value to Swanson alone, or to VanderWeele-Ding 2017 instead of Swanson-VanderWeele 2020 *Epidemiology* — the MR-specific E-value is Swanson-VanderWeele 2020.
- Forgetting the "report at CI-bound nearest null" convention (the conservative E-value is the one typically reported).

---

## §11. Limitations (Coherent Pleiotropy & Practical Caveats)

### Required content
- **Coherent pleiotropy**: if every variant in a burden affects expression AND outcome through the *same* unmeasured confounder, none of the three over-id axes can detect it:
  - HEIDI-rv has $\mathcal{O}(1/m)$ insensitivity to single-variant deviations.
  - Annotation-class concordance passes if the confounder is class-shared.
  - Cell-type concordance passes if the confounder is cell-type-shared.
- Mitigation = Cinelli-Hazlett 2025 IV partial-$R^2$/RV + Wang-Tchetgen 2018 bias bounds. Explicit framing: "honest limit, not a fix."
- **Cell-type resolution claim is incomplete until TenK10K Phase 2** (liver, adipose). v1 is PBMC-only. PCSK9, HMGCR are hepatocyte-biased.
- Bridge: pQTL-anchor mediator-scale normalization lets a hepatocyte-relevant exposure be approximated via UKB-PPP / deCODE plasma pQTL even if cis-eQTL is unavailable.
- **Power weak at typical rare-variant MAC** (carrier counts 10–50): per-mask first-stage F often single-digit; joint K=3 F may rescue some genes but not most. AR matters here because Stock-Yogo $F > 10$ would strip nearly all signals.
- **Data-access reality**: TenK10K rare-variant Zenodo files are 214–260 byte placeholders; OneK1K `rv_sign.txt` is filtered to significant only. Real rare-variant exposure inputs not publicly downloadable.

### Required citations
- Wang L, Tchetgen Tchetgen EJ. 2018. *JRSS-B* 80(3):531–550. doi:10.1111/rssb.12262 — partial-identification bounds + multiply-robust estimators under invalid-IV.
- Cinelli-Hazlett 2025 (above).
- Optional: Bowden J et al. MR-Egger / weighted-median pleiotropy literature if the writer extends discussion.

### Common pitfalls / things that would be WRONG
- Claiming any of the three over-id axes detects coherent pleiotropy.
- Asserting hepatocyte coverage in v1 (PBMC-only).
- Citing Stock-Yogo $F > 10$ as the rejection threshold (the whole point is that AR sidesteps it).

---

## § Pre-emptive Issues List (Traps to Flag if the Writer Steps On Them)

### Citation traps (audit-flagged)

1. **Robins 1994 journal**: must be *Comm Stat Theory Methods* 23(8):2379, doi:10.1080/03610929408831393. Any draft saying *Biometrics* is **WRONG**. (A math-foundations agent previously confabulated this.)
2. **Cinelli-Hazlett dual cite**: 2025 *Biometrika* asaf004 is **primary** (IV-specific); 2020 *JRSS-B* 82(1):39 is **secondary** (OLS origin). Citing only one or swapping primary/secondary = flag.
3. **Madsen-Browning 2009**: must appear as the burden weighted-sum citation in §2. *PLOS Genet* 5(2):e1000384.
4. **STAAR**: default to Li X 2020 *Nat Genet* 52:969 for annotation weights (NOT Li Z 2022 *Nat Methods* STAARpipeline). FLAG if writer picked 2022.
5. **RGC-ME first author**: Sun **KY** (Kathie Y. Sun), NOT Sun BB (Benjamin B. Sun). 2024 *Nature* 631(8021):583–592.
6. **CMC vs CAST**: Morgenthaler-Thilly 2007 *Mut Res* 615:28 = **CAST** (not CMC). CMC = Li-Leal 2008 *AJHG* 83:311. FLAG if reversed.
7. **Burgess-Butterworth-Thompson 2013**: *Genet Epidemiol* 37(7):658, NOT *Stat Med*.
8. **Ray 2025 *AJHG* 112(7):1597**: NOT "Ge 2025" (confabulation).
9. **Sun BB for §8 pQTL anchor is fine** (INTERVAL 2018 and UKB-PPP 2023 are real Sun BB papers) — but Sun BB is **NOT** the RGC-ME author. Different person.
10. **SAIGE-QTL**: Zhou W, Cuomo ASE et al. 2024. medRxiv 2024.05.15.24307317.
11. **TenK10K**: Cuomo et al. 2025. medRxiv 2025.03.20.25324352.
12. **Genebass**: Karczewski 2022 *Cell Genomics* 2:100168. doi:10.1016/j.xgen.2022.100168.
13. **Dhindsa 2023 *Nature* 622:339** for rare-variant pQTL (UKB-PPP rare-variant) — distinct from Sun BB 2023 *Nature* 622:329 (UKB-PPP common-variant).

### Mathematical traps

14. **HEIDI-rv weights**: MUST be non-zero eigenvalues of $V_\delta$ itself, NOT eigenvalues of $V_\delta^+ V_\delta$ (the latter projector has only 0/1 eigenvalues and collapses the test to plain $\chi^2_{m-1}$ — anti-conservative). This is a signature correctness claim of HEIDI-rv.
15. **Sargan J df**: $\chi^2_{K-1}$, NOT $\chi^2_K$ (one df consumed by the point estimate at the AR argmin).
16. **K=1 AR df**: $\chi^2_1$ critical value $c = 3.841$ at $\alpha = 0.05$, NOT $\chi^2_K$ for $K > 1$.
17. **Four CI shapes**: bounded / disconnected / whole_line / empty — match briefing §4.3 table. Sign of $A = \hat b_x^2 - c \mathrm{SE}_x^2$ paired with sign of discriminant determines shape. The discriminant convention: with $B = -2(\hat b_x \hat b_y - c \rho \mathrm{SE}_x \mathrm{SE}_y)$, $\Delta = B^2 - 4AC$ (briefing's $\Delta = B^2 - AC$ uses a different normalization but the geometric classification is identical).
18. **AR cross-term sign**: $V(\beta_0)$ contains $\mathbf{-}2\beta_0 D_y R_{xy} D_x$ (negative); ordering is $D_y R_{xy} D_x$ (NOT $D_x R_{xy} D_y$) — sign and ordering both verifiable against `mrAR_multi.R:157`.
19. **AR is weak-IV-robust because**: denominator never inverts $\hat b_x$; pivots on the reduced-form moment, not the first-stage estimator. Must NOT be paired with a Stock-Yogo $F > 10$ filter.
20. **Estimand formula** (briefing Delta 2): $\beta_{\mathrm{burden}} = \sum_j \pi_j (\gamma_j/\alpha_j)$ with $\pi_j \propto w_j^2 p_j(1-p_j) \alpha_j^2$ — IVW (squared-weight) form. Unweighted average = WRONG. First-power form $\pi_j^{(1)} \propto w_j p_j(1-p_j) \alpha_j$ is the **single-pre-specified-burden alternative**, NOT the rvSMR-adopted form.
21. **Sign of robustness value**: $RV = (\sqrt{t^2 + 4} - t)/2$ (positive). Reversed sign = WRONG.
22. **E-value approximation**: $RR \approx \exp(0.91 \beta_{\mathrm{std}})$ (Swanson-VanderWeele 2020). Coefficient 0.91 is the constant — flag any other constant.

### Methodological commitment traps

23. **K ≥ 3 masks per gene** committed. K = 1 sanity-check on the 5 RCT genes is plumbing-only (LoF only on Genebass side); the headline method requires K ≥ 3 = pLoF / missense:LC / regulatory.
24. **Cauchy combination test (CCT) NOT applied across K masks within a gene** — would defeat AR over-id.
25. **Frequentist end-to-end**: SAIGE-GENE+ Step 2 burden **score test** (honest sampling SE), NOT RareEffect Step 4 BLUP / PEV (Bayesian → breaks $\chi^2$ reference).
26. **No-defier-within-mask monotonicity** must appear in §3 assumptions explicitly (non-default Imbens-Angrist adaptation).
27. **No Z × X interaction** (Robins 1994 SMM no-interaction condition) must appear in §3.
28. **Linearity within mask** must appear in §3.
29. **TenK10K v1 is PBMC-only** — hepatocyte resolution awaits Phase 2; this caveat belongs in §11.

### Open citation TODOs (flag if the writer picked one without justification)

30. **Delta-method variance source (a)**: candidates Burgess-Thompson 2017 textbook / Thomas 2007 / Rothman-Greenland; Bowden-Vansteelandt 2011 is WRONG (that's case-control SMM). Flag whatever the writer picks.
31. **CAST vs CMC (b)**: if mentioned, verify attribution (CAST = Morgenthaler-Thilly 2007 *Mut Res* 615:28; CMC = Li-Leal 2008 *AJHG* 83:311). Default expectation: Madsen-Browning 2009 alone covers the linear weighted-sum case; CAST/CMC only needed if explicitly invoked.
32. **STAAR 2020 vs 2022 (c)**: default is 2020 (Li X *Nat Genet* 52:969); if writer picked 2022 (Li Z *Nat Methods* STAARpipeline) — flag for justification.
33. **Sign-concordance reference (d)**: Han-Eskin 2011 is WRONG. If writer picks Owen 2009 or Whitlock 2005 — accept; if not specified — flag as TODO.

### Structural traps

34. Each of §1–§10 must follow the **assumption → derivation → estimand → finite-sample concern** structure (HANDOVER §10 spec). Flag any section missing one of the four phases.
35. **Target length ~30 pages** — flag if dramatically over or under.
36. **§3 formal identification proof** is the work item — it had not been written prior to this writeup. Must be present and load-bearing.
37. **Argos / DFCI** are NOT the compute substrate — rvSMR runs on MGH/Broad. Flag if the writer mentions argos/DFCI in implementation discussion.

---

*End of CHECKER_RUBRIC.md. Diff the writer's draft against this rubric.*
