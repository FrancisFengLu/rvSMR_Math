# RARE — Validation Strategy Review (for rvSMR planning)

*Compiled 2026-06-01 from the paper landing page, the GitHub repo, and PubMed.*

---

## 1. Paper metadata

- **Full citation**: Cheng Y, Ruan X, Lu X, Yang Y, Wang Y, Yan S, Sun Y, Yan F, Jiang L, Liu T.
  *Accounting for the impact of rare variants on causal inference with RARE: a novel multivariable Mendelian randomization method.*
  **Briefings in Bioinformatics**, 2025, **26(3)**: bbaf214.
  DOI: [10.1093/bib/bbaf214](https://doi.org/10.1093/bib/bbaf214)
- **Published**: 15 May 2025
- **Corresponding authors** (all at the Research Center of Biostatistics and Computational Pharmacy, China Pharmaceutical University, Nanjing):
  - Fangrong Yan
  - Liyun Jiang
  - Tiantian Liu
- **First authors / package maintainers**: Yu Cheng (yucheng.cpu@foxmail.com), Xinjia Ruan (ruan.cpu@foxmail.com)
- **PMID**: 40370099
- **Landing page**: https://academic.oup.com/bib/article/26/3/bbaf214/8131742

Acronym expansion: **R**are variants **A**ccounting for multiple **R**isk factors and shared horizontal pl**E**iotropy. It is an MVMR (multivariable MR) extension, not univariable.

---

## 2. Simulation validation

### Generative model
Individual-level linear additive structural equations (paper Eqs. 7–8, replicated in `simulation/simulationl.R`):

- Primary exposure:   X₁ = G·δ  + G·κₓ + U·ψₓ + ε
- Secondary exposure: X₂ = G·γ₂              + ε
- Outcome:            Y  = β₁·X₁ + β₂·X₂ + G·κᵧ + G·θ + U·ψᵧ + ε

where G is the genotype matrix, U is a latent confounder matrix, κ terms produce **correlated horizontal pleiotropy (CHP)**, and θ produces **uncorrelated horizontal pleiotropy (UHP)**.

### Parameter sweeps (from the simulation script and main text)
| Knob | Values |
|---|---|
| Sample size n | **fixed at 50,000** for both reference and calculation panels (no n-sweep in main text) |
| # common SNPs | 1,000 |
| # rare SNPs | 20 |
| Rare-variant MAF | < 0.01 (common: > 0.05) |
| # latent confounders | 20 |
| Number of exposures K | primary K = 2 (one focal + one secondary); extension supports up to 10 |
| True causal effect β₁ | 0 (null) or 0.01 (signal); β₂ = 0.2 |
| LD r | 0.2 (common SNPs), 0.5 (rare SNPs) |
| Rare-variant heritability sweep | h²_γ1 ∈ {0, 10⁻², 10⁻¹} with h²_δ = (4/5)·h²_γ1 and h²_κx = (1/5)·h²_γ1 |
| CHP / UHP heritability | h²_γ2 = h²_κy = h²_θ = 0.02 (fixed) |
| CHP-UHP correlation ρ_κ | 0.2 – 0.5 |

### Replicates
**300 Monte Carlo replications per cell** (`n_parallel <- 300`, run with 30 cores in the script). The main-text figures appear to summarise 100 of these in box plots; the script saves all 300.

### Metrics
- **Type I error** (rejection when 95 % CI excludes 0 under β₁ = 0)
- **Statistical power** (rejection rate under β₁ = 0.01)
- **Point-estimate bias** (β̂₁ − β₁_true)
- **95 % CI width / coverage** (CIs reported but coverage is *not* explicitly tabulated; credible-interval coverage is also not formally evaluated — see §5)

Compared against: **MV-IVW, MV-Egger, MV-LASSO, GRAPPLE**, each in a "naive" (rare variants ignored) and "standard / RARE-augmented" form.

### Weak-IV sweep
**Not present.** The paper does not vary the F-statistic or examine F < 10 / F < 1 regimes. Authors merely note that the LD-aware formulation "helps mitigate" weak-instrument bias.

### Reference figures
- **Figure 2A**: Type-I error vs LD / pleiotropy conditions
- **Figure 2B**: Statistical power across pleiotropy regimes
- **Figure 2C**: Box plots of β̂₁ over replications (bias visualisation)
- Coverage curve: **none** in main text (CIs shown only as error bars)

---

## 3. Real-data validation

- **Dataset**:
  - GWAS summary statistics from the **Global Lipids Genetics Consortium** (HDL, LDL) and from the GWAS Catalog / fastGWA repositories for the outcomes
  - **UK Biobank** individual-level genotypes for PRS construction (rare-variant arm)
- **Gene panel**: PRS built genome-wide from rare variants in UKBB; no curated drug-target gene panel is used.
- **Exposures**: HDL-C, LDL-C
- **Outcomes**: Type 2 Diabetes (T2D), Coronary Atherosclerosis
- **Positive controls**: well-established lipid → cardiometabolic relationships are used as benchmarks (LDL ↑ → CAD risk ↑; HDL inverse association with CAD). RARE recovers HDL → T2D β = −0.219, p = 0.003 and LDL → T2D β = 0.264, p = 0.007, with a significant rare-variant PRS contribution β_PRS = 0.313, p < 0.001.
- **Negative controls**: **none** reported.
- **Drug-target positive controls (PCSK9-style)**: **none** explicitly reported.
- **Blinded validation against pQTL / external GWAS**: **none**. No held-out replication panel; no pQTL cross-check.
- **Sensitivity analyses**: a comparison of UKBB-genotype-based simulations vs synthetic-genotype simulations is described as producing "comparable results", but this is robustness rather than blinded replication.

---

## 4. Code / data availability

- **GitHub URL**: **https://github.com/Hide-in-lab/RARE** (branch `main`)
- **Install**: `devtools::install_github('Hide-in-lab/RARE@main', force = T)` — requires Rtools / Xcode because the inner Gibbs sampler is in C++.
- **Repo structure**:
  ```
  RARE/
  ├── R/                          # R wrappers; only RcppExports.R is visible
  ├── src/
  │   ├── multi_exposures.cpp     # the Gibbs sampler (functions adapt(), calculateP())
  │   ├── RcppExports.cpp
  │   └── rcpp_hello_world.cpp    # boilerplate
  ├── simulation/
  │   └── simulationl.R           # (sic) the full simulation driver
  ├── man/                        # roxygen docs
  ├── image/Github_RARE.jpg
  ├── DESCRIPTION, NAMESPACE, README.md
  ```
- **Specific file paths**
  - (a) Simulation generator + driver: `simulation/simulationl.R`  (note the typo "simulationl", not "simulation1")
  - (b) Real-data analysis driver: **not present in the public repo**
  - (c) Figure-generating notebook/script: **not present**; figures appear to be produced ad hoc from the CSVs written by `simulationl.R`
  - (d) Gibbs sampler core: `src/multi_exposures.cpp` — function `adapt()` samples Gamma/inverse-Gamma variance hyperparameters (shape = scale = 1) and Normal effects (init var = 0.1); posterior summaries from the second half of iterations (burn-in = first half).
- **Data deposit (Zenodo / Figshare / OSF)**: **none**. UKBB requires application; GWAS summary stats are publicly downloadable but the exact files / hashes are not pinned.
- **License**: not specified in the README (no LICENSE file noted).
- **Reproducibility level**: simulations *are* end-to-end reproducible from `simulationl.R`. Real-data figures are **not** reproducible from the repo as-is (no analysis driver, no preprocessing scripts, no PRS-construction script).

---

## 5. What's specific to rvSMR that RARE does *not* address

1. **Frequentist coverage**: RARE is Bayesian (Gibbs) and reports CIs as posterior 95 % intervals but never tabulates calibration / nominal coverage. rvSMR's commitment to a frequentist **Anderson–Rubin (AR) interval** means we *must* report empirical coverage vs nominal across the parameter grid — RARE provides no template here.
2. **Weak-instrument regime**: RARE has no F-statistic sweep and no F < 10 / F < 1 stress tests. For rvSMR, the AR statistic's main selling point is weak-IV robustness, so we need an explicit weak-IV sweep that RARE skips.
3. **Annotation-class / cell-type concordance**: RARE is trait-level (HDL/LDL → T2D/CAD). It does not stratify instruments by functional annotation, cell type, or tissue. rvSMR's annotation-class concordance axes are not validated anywhere in RARE.
4. **Univariable focus**: RARE is an MVMR method; it requires ≥ 2 exposures and pleiotropy decomposition. rvSMR's univariable framing avoids the secondary-exposure modelling burden but loses RARE's CHP/UHP separation.
5. **Drug-target / RCT-anchored positive controls**: RARE uses only lipid–disease "textbook" positive controls. rvSMR should add PCSK9–LDL, HMGCR–LDL, IL6R–CHD, etc. — gene-level RCT-validated effects RARE never tests.
6. **Negative controls and blinded replication**: RARE has neither. rvSMR can differentiate itself by including negative controls (e.g., implausible exposure–outcome pairs) and an external held-out cohort.
7. **Reproducible figure pipeline**: RARE has no figure-generating notebook. rvSMR should ship a single `make figures` (Snakemake/Quarto) that regenerates every panel from raw simulation CSVs.

---

## 6. Top-3 takeaways for rvSMR validation

1. **Cover the weak-IV gap RARE leaves open**: build an explicit F-statistic sweep (F ∈ {0.5, 1, 2, 5, 10, 30, 100}) and report AR-interval coverage / power across it — this is the single biggest validation differentiator we can claim.
2. **Add drug-target positive controls and matched negative controls**: PCSK9/HMGCR/IL6R for positives and a curated set of biologically implausible exposure–outcome pairs for negatives; RARE's "HDL→CAD recovered, therefore valid" is too weak a bar.
3. **Ship a single-command reproducible figure pipeline with a Zenodo-archived simulation CSV bundle**: RARE has no data deposit and no figure script; meeting that bar is cheap and instantly improves perceived rigor relative to the prior art.

---

### Sources
- Paper landing page: https://academic.oup.com/bib/article/26/3/bbaf214/8131742
- DOI: https://doi.org/10.1093/bib/bbaf214
- GitHub: https://github.com/Hide-in-lab/RARE
- Simulation script: https://raw.githubusercontent.com/Hide-in-lab/RARE/main/simulation/simulationl.R
- Gibbs sampler source: https://github.com/Hide-in-lab/RARE/blob/main/src/multi_exposures.cpp
- PMID: 40370099
