# Validation strategies in three Mendelian Randomization papers

Research compiled for rvSMR validation planning. Sources: arXiv preprints, PMC, OUP/Biometrics, and the MRBEEX GitHub repository. Where the published Biometrics or BiB main text was not directly accessible, the arXiv preprint or PMC mirror was used; this is flagged in-line.

---

## A. Wang & Kang 2022, Biometrics

### A.1 Paper metadata + URL

- **Title**: "Weak-Instrument Robust Tests in Two-Sample Summary-Data Mendelian Randomization"
- **Authors**: Sheng Wang, Hyunseung Kang (UW-Madison Statistics)
- **Journal**: Biometrics 78(4):1699-1713 (Dec 2022). DOI 10.1111/biom.13524
- **Publisher landing**: https://academic.oup.com/biometrics/article/78/4/1699/7460098
- **Open arXiv**: https://arxiv.org/abs/1909.06950 (v3, 7 Jun 2021)
- **PubMed**: https://pubmed.ncbi.nlm.nih.gov/34213007/
- **bioRxiv**: https://www.biorxiv.org/content/10.1101/769562v1
- Methods proposed: mrAR (Anderson-Rubin), mrK (Kleibergen), mrCLR (conditional likelihood ratio), and a derived point estimator mrLIML.

### A.2 Simulation

**Generative model** (Section 3, Eq. 1/3, structural model):
- Two independent samples, both individual-level data, then summary statistics computed.
- Sample sizes n1 = n2 = 100,000 (matched to the BMI->SBP real-data analysis).
- L = 100 instruments, generated as Binomial(2, p_j) with p_j ~ Unif(0.1, 0.9), mimicking SNP genotypes.
- Bivariate normal random errors with endogeneity parameter rho = 0.1.
- True exposure effect beta varied over a grid; tested H0: beta = 0 and H0: beta = 1.

**Weak-IV regime** (concentration-parameter framing, NOT a simple F-stat):
- IV-exposure relationship gamma is set in a local-to-zero range gamma in [{(r-0.5)/n1}^(1/2), {(r+0.5)/n1}^(1/2)].
- r is varied across {1, 4, 16, 25} and (in supplement) r = 50.
- The authors prove r approximately equals the first-stage F-statistic, so r=1 corresponds to F approx 1 (extremely weak); r=25 to F approx 25 (moderately strong). The empirical replication later achieves F = 25.46 (160 IVs) and F = 58.14 (25 IVs).
- 1,000 simulation replicates per cell.

**Comparators**:
- MR-Egger (robust radial regression with 2nd-order correction; Bowden 2018) = "MR-Egger.r"
- MR-RAPS (profile-likelihood, squared-error loss; Zhao 2020) - equivalent to exact IVW
- Weighted Median (Bowden 2016a)
- Q-statistic (Bowden 2019) - used as the comparator for invalid-IV detection
- Supplement: six MR-Egger variants (Egger.original, Egger.robust, simex, RadialMR1-3)

**Metrics**:
- Type I error (size) under H0: beta = beta0, with beta0 from -2 to 2.
- Power (rejection rate) under fixed nulls H0: beta = 0 and H0: beta = 1.
- Coverage probability of 95% CIs at varying r (their Figure 6 is the canonical weak-IV coverage figure: shows mrAR keeps 95% coverage even when other methods collapse).
- Invalid-IV setting (Sec 3.2): proportion of invalid IVs varied 10-90% of L; direct-effect magnitude alpha_j in {0.003, 0.01, 0.05} representing local, intermediate, and non-local invalidity; rejection rate of mrAR vs Q-stat across beta0 in [-4, 5].
- Correlated-IV setting (Sec 3.3): tests the LD-adjustment procedure including mis-specified correlation matrices.

**Key tables/figures**:
- Figure 1: Type I error vs beta0 at each r. Figure 2: power. Figure 3: rejection rates for invalid-IV detection. Figure 6: coverage probability vs IV strength (the canonical AR weak-IV coverage demo). Supplement Figures S1-S16 extend each panel.

### A.3 Real-data validation

- **Replication of Zhao et al. (2020)**: BMI -> systolic blood pressure (SBP). Effect known a priori to be positive.
- **Datasets**: BMI from GIANT GWAS ("BMI-MAL" for exposure summary stats; "BMI-FEM" for IV pre-screening). SBP from UK Biobank GWAS ("SBP-UKBB"). All publicly available via Zhao 2020.
- **Two IV sets**:
  - 25 SNPs at genome-wide significance (p < 5e-8), default for MR-Base; resulting F = 58.14, Q-test p = 5.58e-8.
  - 160 SNPs at p < 1e-4 (weaker); F = 25.46, Q p = 5.73e-61.
- Each pair of SNPs >=10 Mb apart, LD r^2 < 0.001.
- **Validation logic**: test H0: beta = 0 (no BMI effect on SBP) with progressively expanding instrument sets from the 3 weakest to all 160. Ideal method should reject for all sets. mrCLR rejected even with 3 weakest IVs; mrAR and mrK had smaller p-values than competitors. Table 1 reports 95% CIs across all methods.
- Note: no cis-MR component; this is global polygenic MR.

### A.4 Code/data

- **Data availability**: Both GWAS summary stats are openly available in CRAN: https://cran.r-project.org/package=mr.raps (Zhao 2018).
- **No standalone GitHub for mrAR/mrK/mrCLR was located in the paper or supplement.** The PDF contains no GitHub URL. Methods were implemented by the authors but not packaged in a public repo cited in the paper. Software used for comparators: R packages `MendelianRandomization` (Yavorska & Burgess 2017), `RadialMR` (Bowden 2019), `mr.raps` (Zhao 2020).
- This is the one paper of the three without a clear public code repo. **GitHub: NOT FOUND in the paper.**

### A.5 Lessons for rvSMR

1. **Weak-IV coverage figure** at fixed nominal level, varying IV strength (their r-parameter) is the canonical demonstration. rvSMR should mimic this: plot empirical coverage of 95% CIs as a function of either F or a concentration-style parameter.
2. **Plot Type I error against the null value beta0**, not just beta0 = 0. Many MR methods are size-correct at zero but fail elsewhere. rvSMR-AR should be tested at multiple non-zero nulls.
3. **Empirical "known-sign" replication** of a positive-control trait pair (BMI->SBP for them) is a low-cost real-data sanity check. Equivalent for rvSMR could be a known cis-effect like LDLR->LDL or PCSK9->LDL.
4. **Power-grid simulations under nested IV strengths and effect sizes** with 1,000 replicates is the standard granularity.
5. **Invalid-IV simulations** with magnitudes spanning local-to-zero (alpha ~ 1/sqrt(n)), intermediate, and non-local (alpha ~ n^-1/4) is the standard pleiotropy stress test.

---

## B. Patel, Lane & Burgess 2024, arXiv MVMR weak-IV

### B.1 Paper metadata + URL

- **Title (correct, per arXiv)**: "Weak instruments in multivariable Mendelian randomization: methods and practice." Note: the user-provided title "Anderson-Rubin tests for Mendelian randomization with weak and possibly invalid instruments" does not match the arXiv record at 2408.09868 - the actual title is the methods-and-practice paper, which does cover AR, Kleibergen and a novel adjusted-Kleibergen for MVMR.
- **Authors**: Ashish Patel, James Lane, Stephen Burgess (MRC Biostatistics Unit, Cambridge)
- **arXiv**: https://arxiv.org/abs/2408.09868 (v1, 19 Aug 2024)
- **PDF**: https://arxiv.org/pdf/2408.09868
- Methods covered: MVMR-AR (Anderson-Rubin), Kleibergen CSK, Andrews 2018 linear-combination robust CS (with coverage-distortion cutoff gamma-hat), and a new adjusted-Kleibergen ("Kleibergen-OH") that corrects for overdispersion heterogeneity in outcome associations.

### B.2 Simulation

**Generative model** (Sec 4.1):
- Two-sample summary data, n_X = n_Y = 5,000 individual-level samples per sample; per-SNP summary statistics computed via univariable OLS.
- Baseline: J = 4 instruments, K = 2 exposures. True effects (theta1, theta2) = (1, 0).
- Outcome model Y = theta1*X1 + theta2*X2 + U; X_k = gamma_k' Z + V_k.
- Errors (U, V1, V2) joint normal with Cor(U, V1) = -0.6, Cor(U, V2) = 0.6, Cor(V1, V2) = 0.3 - deliberately strong endogeneity in opposite directions for the two exposures (OLS biased downward for theta1, upward for theta2).
- Instruments correlated: Z ~ multivariate normal, off-diagonal correlation up to 0.4 (a, a' construction).

**Weak-IV regime** - 2-parameter:
- Unconditional strength mu > 0 (overall magnitude).
- Conditional strength xi in [0, 1]: when xi = 0, gamma1 and gamma2 are linearly dependent (no conditional strength); when xi = 1, instruments are perfectly partitioned across exposures.
- gamma1 = [1+xi, 1+xi, 1-xi, 1-xi]' * 0.2 * n_X^{-1/2} * mu; gamma2 = [1-xi, 1-xi, 1+xi, 1+xi]' * n_X^{-1/2} * mu.
- F-statistic range covered: cells where the minimum conditional F crosses 10 (their Figure 3 shows F approx 6, 10, and higher).
- Supplement explores larger J (Figure S2).

**Comparators**:
- Non-robust: multivariable IVW, multivariable GMM, Wald-based CS.
- Robust: Anderson-Rubin CSAR, Kleibergen CSK, Andrews linear-combination CSR (with gamma_min = 0.01), and (later) the new Kleibergen-OH.

**Metrics**:
- Bias in GMM point estimates plotted vs min conditional F-statistic AND vs distortion cutoff gamma-hat (Fig 1).
- Coverage probability of 95% CS for theta_0 = (1, 0) across mu and xi (Fig 2 top row).
- Power: probability the 95% CS excludes the null (0, 0) (Fig 2 bottom row).
- Coverage under selective reporting: only report results with min cond F >= 10 vs no screening (Fig 3) - shows IVW coverage drops from ~90% to <40% under screening. This is a novel publication-bias / selection-bias check.
- Inference-after-selection (Sec 4.5): coverage and CS area when choosing between a "core" 4-IV set and a "full" 8-IV set, using post-selection rules based on conditional F or distortion cutoff (Fig 4).

**Key figures**: Figures 1-5. Fig 5 is the heat-map of 100 robust vs non-robust confidence sets across sub-sample sizes n_X in {2000, 10000, 20000, 50000} in the empirical example.

### B.3 Real-data validation

- **Exposures**: Mean Arterial Pressure (MAP) and Pulse Pressure (PP) - chosen because they are highly correlated (r = 0.56 in UK Biobank), making it intentionally hard to find conditionally strong instruments.
- **Outcome**: Stroke incidence.
- **Datasets**: UK Biobank for MAP & PP associations (n_X = 367,283, mostly European). GIGASTROKE consortium for stroke (n_Y = 727,571 European; Mishra et al. 2022).
- **Instruments**: top 20 uncorrelated variants for each exposure, up to 40 IVs total.
- **Full-sample headline**: cond F = 73.48 (MAP), 113.10 (PP); distortion cutoff converged to gamma_min = 0.05. MAP effect 95% CI Andrews [0.031, 0.049], PP [0.015, 0.027] (per-mmHg log-OR).
- AR CS was empty (signals overdispersion). Kleibergen-OH (new method) gave non-empty CS consistent with Wald/Andrews.
- **Validation move**: re-run with random sub-samples of n_X in {2000, 10000, 20000, 50000} from the same UK Biobank pool - induces realistic weak-IV settings while keeping the truth fixed. Plot heatmaps of 100 CSes per sub-sample size (Fig 5). At n_X = 1000, cond F ~ 2.55, distortion cutoff 0.21.

### B.4 Code/data

- **R package GitHub**: https://github.com/ash-res/mvmr-weakiv/
- Cited directly in the paper text (Sec 1) as "available at github.com/ash-res/mvmr-weakiv/".

### B.5 Lessons for rvSMR

1. **Two-parameter weak-IV stress** is sharper than just sweeping F: separately vary unconditional and conditional strength. For rvSMR-AR with K cell types, an analog is overall rare-variant burden strength vs cell-type-specific independence.
2. **Selective-reporting / Goodhart simulation** (Fig 3): show what happens to coverage if practitioners only report when F >= 10. Highly cited concern.
3. **Plot bias vs both F-stat AND distortion cutoff** - the latter is a novel reliability metric (max coverage loss of the non-robust set) that rvSMR could adopt for cell-type-resolution outputs.
4. **Sub-sampling-based real-data weak-IV induction**: take a large UK Biobank GWAS and downsample to create weak-IV scenarios where the strong-IV result is known. This is a clean validation tactic - rvSMR could downsample a large eQTL or rare-variant burden GWAS.
5. **Heatmap of multiple CS realizations** (their Fig 5) communicates coverage-shape uncertainty better than a single CS per scenario.
6. **Power AND CS-area metrics**: when CS can be unbounded, classical "power" is augmented by CS-area-vs-null reporting.

---

## C. Yang, Lorincz-Comi, Li, Zhu 2025, Briefings in Bioinformatics, cis-MRBEE

### C.1 Paper metadata + URL

- **Title**: "A multivariable cis-Mendelian randomization method robust to weak instrument bias and horizontal pleiotropy bias"
- **Authors**: Yihe Yang, Noah Lorincz-Comi, Gen Li, Xiaofeng Zhu (Case Western Reserve, Population & Quantitative Health Sciences)
- **Journal**: Briefings in Bioinformatics 26(3):bbaf250, 2025.
- **PMC mirror**: https://pmc.ncbi.nlm.nih.gov/articles/PMC12140020/
- **Publisher landing** (gated): https://academic.oup.com/bib/article/26/3/bbaf250/8154932
- **ResearchGate**: https://www.researchgate.net/publication/392441829
- Method: extends MRBEE (bias-corrected MVMR via estimating equations) to the cis-MR setting using SuSiE for variant selection plus a double-penalized minimization (IPOD) that handles horizontal pleiotropy. Companion methods: MRBEE-Mixture, MRBEE-TL (transfer learning).

### C.2 Simulation

(From the PMC mirror and BiB landing page abstract.)

**Cis-window / fine-mapping**:
- Fine-mapping: SuSiE. Variants are flagged as "informative" if they are in a 95% credible set produced by SuSiE.
- Cis-window size not stated in main text (paper uses pre-defined locus boundaries inherited from the real-data ANGPTL3 / CR1 analyses). One PMC reading reported m=200 IVs with 10 exposures; another reported m=100 IVs with 4 exposures - the paper sweeps multiple configurations.

**Generative model**:
- Multi-exposure cis-MR setup. Each exposure has 3 randomly assigned "informative xQTLs" plus 1-2 horizontal-pleiotropy variants per exposure.
- Local heritability of exposures: 0.3 (this is the cis-region heritability, distinct from polygenic models).
- Local heritability of outcome: 0.001, split at a 1:1 ratio between causal-exposure contributions and pleiotropy.
- xQTL sample sizes: 300; 3,000; 30,000 (covers small to medium pQTL/eQTL cohorts).
- Outcome GWAS sample size: 500,000.
- Replication: 500 simulations.

**Comparators**:
- cis-MVIVW (cis-multivariable IVW, the de facto standard).
- PC-GMM (principal-component GMM for correlated IVs).
- MRBEE-IPOD (the global / non-cis MRBEE).
- Adjusted cis-MVIVW (with sparse prediction).
- TGFM (Tissue-of-action Gene-level Fine-Mapping).
- TGVIS (cis-MR via TGFM-style variant selection).

**Metrics**:
- Bias of causal-effect estimate (boxplots vs sample size).
- Statistical power (rejection at p < 0.05).
- Type I error (target 0.05).
- Coverage frequency (target 95%).
- Behavior across exposure-outcome correlation structures.

### C.3 Real-data validation

**Locus 1 - ANGPTL3 region for lipids**:
- Genes in the cis-window: ANGPTL3, APOA1, APOA5, APOC1, APOC3, PCSK9 (a multi-gene region treated as a multivariable cis-MR problem).
- Exposures: pQTLs/eQTLs for these proteins. Reported pQTL dataset size ~69,016 for the ANGPTL3 region.
- 313 variants selected via r^2 < 0.64 clumping from European populations.
- LD reference: 9,680 UK Biobank European individuals, 9.3M variants.
- Outcomes: LDL-C, HDL-C, triglycerides; n = 1,320,016 (GLGC-like meta).
- **Headline finding** (used as biological validation): credible set of {APOA1, APOC1, PCSK9} identified as causal for lipid traits - PCSK9 is the canonical positive control.

**Locus 2 - CR1 for Alzheimer's**:
- Gene: CR1, tested across multiple GTEx tissues and single-cell eQTLs (sceQTLs).
- 51 variants selected.
- Outcomes: AD case-control (n = 103,772); CSF biomarkers (A-beta-42, p-Tau).
- Finding: CR1 expression in specific brain regions is potentially causal for AD.

**Datasets**: pQTLs (likely deCODE/UKB-PPP), GTEx (tissue eQTLs), single-cell eQTL datasets for the CR1 sceQTL piece, UK Biobank for LD reference, GLGC-like meta-analysis for lipids, large AD GWAS for the Alzheimer's outcome. Specific cohort names were not all explicit in the PMC excerpt available.

### C.4 Code/data

- **GitHub**: https://github.com/harryyiheyang/MRBEEX
- The repo packages cis-MRBEE, MRBEE-IPOD, MRBEE-Mixture, and MRBEE-TL together.
- Includes `Tutorial-of-CisMRBEE.pdf` (step-by-step), `MRBEEX_0.1.0.pdf` (package docs), `CisMRBEE_Real_Data.zip` example data via Dropbox.
- SuSiE integration is built into MRBEE-IPOD and MRBEE-Mixture for selecting exposures with non-zero causal effects.

### C.5 Lessons for rvSMR

1. **SuSiE 95%-credible-set filtering** is the canonical fine-mapping step for cis-MR papers in 2024-25. rvSMR's analog of an IV-selection step inside a cis-window should follow this convention so reviewers find it familiar.
2. **Multi-gene cis-window** (treat ANGPTL3 + 5 nearby genes as a multivariable problem) is the modern framing - much more impactful than single-gene cis-MR. rvSMR's cell-type-resolution should be presented analogously.
3. **xQTL sample-size sweep {300, 3K, 30K}** matches realistic protein/expression QTL cohort sizes. rvSMR should sweep at least small to medium burden-test sample sizes.
4. **Lipid panel is the canonical positive control**: PCSK9 -> LDL is the gold-standard hit any cis-MR method must reproduce. APOA1/APOC1 are common secondary validators.
5. **AD/Brain validation**: CR1 -> AD via brain eQTLs is the second-favorite hard real-data benchmark in cis-MR papers.
6. **Comparator set**: cis-IVW, PC-GMM, TGFM/TGVIS are the three families a 2025 cis-MR paper is expected to benchmark against.
7. **500-replicate boxplots** of estimates across sample-size cells is the modern visualization convention.

---

## D. Cross-method synthesis

### Shared validation tactics across all three

- **Coverage of 95% CIs as the headline metric** (not just point-estimate bias). All three papers feature a coverage figure where their robust method holds at 95% while non-robust comparators collapse.
- **F-statistic (or analog) as the x-axis** of the canonical figure. Wang-Kang uses the concentration parameter r ~ F; Patel-Lane-Burgess uses the minimum conditional F; cis-MRBEE uses sample size as a proxy.
- **A nonrobust comparator (IVW / GMM / Wald)** is always shown to demonstrate the weak-IV failure case. This is what makes the robust method look essential.
- **Type I error checked at multiple null values** (not only at zero). Wang-Kang most explicit; Patel-Lane-Burgess via coverage-of-true-vector; cis-MRBEE via type-I-error simulation cells.
- **One headline real-data application** where the truth is at least partially known (BMI->SBP positive; MAP->stroke positive; PCSK9->LDL positive).
- **Public code/reproduction availability** is now standard - Patel-Lane-Burgess and cis-MRBEE have GitHub repos; Wang-Kang is the older paper and only ships data via the `mr.raps` CRAN package.

### Unique to each

- **Wang-Kang**: Invalid-IV stress test with magnitudes spanning local-to-non-local (alpha ~ 1/sqrt(n) up to 1/n^{1/4}) plus the proportion of invalid IVs sweeping 10-90% - this is the single-IV-paper's compensation for not having a multi-IV setup.
- **Patel-Lane-Burgess**: Selective-reporting Goodhart simulation (Fig 3); distortion-cutoff metric (max coverage loss of non-robust set) as a new instrument-strength summary; sub-sampling a real GWAS to manufacture weak-IV scenarios while preserving the truth (Fig 5 heatmap).
- **cis-MRBEE**: Multi-gene cis-window treatment; SuSiE-credible-set IV filtering; sceQTL/single-cell tissue resolution; 500K-N outcome GWAS coupled to small (300-30K) xQTL cohorts.

### Canonical "weak-IV coverage demonstration" standard

The AR literature converges on: **a 2D figure with instrument-strength on x-axis (F-stat or concentration parameter r) and empirical 95% coverage on y-axis**, with one line per method. The proposed AR/Kleibergen method must hold a flat horizontal line at 0.95 across the full F range; competitor methods should drop sharply at low F. This is exactly what Wang-Kang Figure 6 and Patel-Lane-Burgess Figure 2 (top row) both look like. For rvSMR, this is the figure to reproduce.

### Canonical real-data sanity panel for cis-MR

For cis-MR / multivariable cis-MR specifically, the de-facto gene-trait positive-control set is:
- **PCSK9 -> LDL** (gold standard, always present)
- **HMGCR -> LDL** (statin-target validation)
- **IL6R -> CRP** or **IL6R -> CHD** (drug-target inflammation example)
- **APOA1, APOC1, APOC3, APOE** for lipid panels
- **ANGPTL3 -> triglycerides/LDL** (modern cis-MR paper favorite, used by Yang 2025)
- **CR1 -> Alzheimer's** (neurodegeneration positive control, used by Yang 2025)
- **LDLR**, **NPC1L1** (ezetimibe-target)

Yang 2025 (cis-MRBEE) covers PCSK9, APOA1/APOC1, ANGPTL3, CR1 - hitting the lipid and AD canonical pairs. PCSK9-LDL was NOT featured in Wang-Kang or Patel-Lane-Burgess because those are non-cis MR papers - they used the polygenic BMI-SBP / MAP-PP-stroke pairs instead.

---

## E. Top-5 takeaways for rvSMR validation

1. **Reproduce the canonical AR weak-IV coverage figure**: x-axis = some IV-strength summary (F-stat for rvSMR's rare-variant burden, or a rvSMR analog of r), y-axis = empirical 95% CI coverage from >=500 replicates, one line per method (rvSMR-AR, Wald-rvSMR, IVW analog). Show flat 0.95 line for rvSMR-AR while comparators drop at low F. This single figure is what makes both Wang-Kang and Patel-Lane-Burgess persuasive.

2. **Test Type I error at multiple non-zero null values**, not only at the global null. Wang-Kang Figure 1 (size vs beta0) is the template; rvSMR should produce the same with beta0 spanning a realistic effect-size grid for each cell type.

3. **Multi-gene cis-window real-data benchmark**: use the lipid panel (PCSK9, HMGCR, APOA1, APOC1, ANGPTL3, LDLR) plus IL6R-CRP and CR1-AD as the canonical sanity-check set. PCSK9-LDL must be a positive hit. This is now the expected real-data validation in any 2024-25 cis-MR paper (Yang 2025 is explicit about this).

4. **Sub-sample a real large-N rare-variant GWAS to manufacture weak-IV scenarios** (Patel-Lane-Burgess tactic). Take a well-powered burden test (e.g., 500K UK Biobank exomes), random-sub-sample to n in {2K, 10K, 50K}, run rvSMR-AR vs Wald-rvSMR, show coverage of the full-sample truth holds for AR but collapses for Wald. Produces a heatmap-style figure (Fig 5 in their paper).

5. **Include SuSiE-credible-set IV selection in the pipeline** even if rvSMR is rare-variant-based: any modern cis-MR reviewer expects fine-mapping to be in the loop. Either run SuSiE on the burden-defining variants (if applicable) or use a rare-variant fine-mapping analog. Use cis-MRBEE's MRBEEX/SuSiE workflow as the precedent to cite.

**Bonus tactic from Patel-Lane-Burgess**: a "Goodhart simulation" showing that if rvSMR users only reported results with F >= 10, the published coverage would drop dramatically. This is novel, high-impact, and shows the AR method is useful precisely because it removes the need for screening.

---

## F. What we could not find

- **Wang-Kang GitHub repo**: not located. PDF and supplement reference only `mr.raps`, `MendelianRandomization`, `RadialMR` CRAN packages used as comparators. The Wang-Kang methods themselves do not appear to be in a public package. If a reproducibility package exists, it is not advertised in the paper.
- **Cis-MRBEE cis-window size in kb**: not explicitly stated in the accessible PMC excerpt or BiB landing page. The paper uses inherited locus boundaries from the ANGPTL3/CR1 real-data analyses but does not pin down a fixed kb window.
- **Cis-MRBEE pQTL source identity**: the n=69,016 pQTL count suggests deCODE or UKB-PPP, but the specific cohort was not stated in the PMC excerpt we could access.
