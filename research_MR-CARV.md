# MR-CARV validation strategy — research notes for rvSMR

Compiled 2026-06-01. Sources: Oxford Academic landing page, PubMed record,
medRxiv preprint listing, GitHub repo `yu-zhang-oYo/mr.carv`, Indiana University
ScholarWorks deposit.

---

## 1. Paper metadata

**Full citation**
Zhang Y, Li M, Haas DM, Bairey Merz CN, Li X, Zhao H, Yan Q, et al.
*A novel two-sample Mendelian randomization framework integrating common and
rare variants: application to assess the effect of HDL-C on preeclampsia risk.*
**Briefings in Bioinformatics**, 2026 Jan 7; 27(1):bbaf649.
doi: [10.1093/bib/bbaf649](https://doi.org/10.1093/bib/bbaf649).
Received 2025-09-25, accepted 2025-11-12. ~14 authors total.

**Lead / corresponding authors and affiliations**
- First author: **Yu Zhang** — Dept. of Epidemiology & Biostatistics, Indiana
  University School of Public Health, Bloomington.
- Corresponding: **Xihao Li** — Depts. Biostatistics & Genetics, UNC Chapel Hill.
- Corresponding: **Nianjun Liu** — Dept. Epidemiology & Biostatistics, Indiana
  University.
- Corresponding: **Qi Yan** — Dept. Obstetrics & Gynecology, Columbia
  University.
- Other notable co-authors: Ming Li (IU), David M. Haas (IU School of Medicine),
  C. Noel Bairey Merz (Cedars-Sinai), Hongyu Zhao (Yale).

**Identifiers and URLs**
- **PMID**: 41499219
- Journal landing page:
  https://academic.oup.com/bib/article/27/1/bbaf649/8416446
- medRxiv preprint (2025-08-20, v1/v2):
  https://www.medrxiv.org/content/10.1101/2025.08.20.25334100
- Author-deposit (IU ScholarWorks):
  https://scholarworks.indianapolis.iu.edu/items/49685650-ab6c-44ad-ad9f-5d58d65b6c47

---

## 2. Simulation validation

### Generative model

- **Genotype source**: real haplotypes from the **1000 Genomes Project**.
- **Common variants**: four IVs selected, one independent SNP per chromosome
  from chromosomes 15–18 with MAF > 0.1.
- **Rare variants**: four burden sets, one per chromosome 19–22, each composed
  of 50 SNPs with MAF in 0.001–0.05 (paper text also quotes 0.001–0.5; the
  binding cutoff for "rare" is 0.05).
- **MAF distribution**: empirical MAF from 1000G; weighting within a burden uses
  the Beta(a₁, a₂) family — Beta(1, 25) (upweights rarer alleles, the STAAR
  default) and Beta(1, 1) (equal weights) are both tested.
- **LD**: deliberately suppressed at the IV level by picking one SNP per
  chromosome so that the four common IVs and the four burden IVs are mutually
  independent. A more realistic LD regime (LD pruning rather than per-chrom
  picking) is only explored in Supplementary Section S5.
- **Exposure model**: linear for continuous exposure, logit for binary exposure.
- **Outcome model**: linear (continuous) or logit (binary), giving four
  exposure/outcome type pairs.

### Mixing of common + rare variants

Rare variants in each gene/window are first collapsed into a burden score
weighted by annotation (the STAARpipeline workflow). Burden-level βs and SEs are
then concatenated with common-variant βs and SEs and plugged into a standard MR
estimator (IVW, dIVW, MR-RAPS). P-values across MAF-weight choices /
annotation-mask choices are combined with the **Cauchy combination test (CCT)**.

### Parameter sweeps

| Knob | Values tested |
|---|---|
| Sample size $N_X = N_Y$ | **10,000** (primary); 1,000 and 2,000 (sensitivity) |
| Causal-variant ratio inside burden | ~15% and ~35% |
| Positive-effect ratio | 100%, 80%, 50% (pleiotropy mix proxy) |
| True $\beta$ | 0.04 (cont/cont), 0.06 (cont/bin), 0.2 (bin/cont), 0.4 (bin/bin) |
| Common-variant effect $\alpha_i$ | fixed at **0.5** |
| Rare-variant effect $\alpha_{j\ell}$ | $c_0 \cdot \lvert \log_{10}\text{MAF}\rvert$ with $c_0 = 0.5$ |
| Beta weights | (1,25) and (1,1) |

### Weak-IV regime — **key gap for rvSMR positioning**

- **No F-statistic sweep is reported in the main paper.** Effect sizes are
  fixed; instrument strength is not parameterized. There is no scenario with
  $F < 10$, and certainly nothing approaching $F < 1$.
- This is the central reason the authors do not need (or address) AR-style
  weak-IV-robust inference. They lean on dIVW + MR-RAPS to soak up modest
  weak-IV bias rather than on an Anderson–Rubin construction.

### Comparators benchmarked

Only three baseline MR estimators, each used in two flavors (common-only vs
MR-CARV common+rare):

- IVW
- dIVW (debiased IVW)
- MR-RAPS

**Not included**: MR-Egger, MR-PRESSO, RARE, classical SMR, MR-ROBIN, MR-cML.
The authors note MR-CARV "can be extended" to MR-Egger / multivariable MR but
provide no empirical comparison.

### Replicates

- **Type I error**: 10,000 simulation iterations per cell, $\beta = 0$.
- **Power**: 1,000 simulation iterations per cell.

### Metrics

Type I error rate, power, and effect-estimate accuracy (bias/precision shown as
scatter plots and tables). No explicit nominal-coverage table for 95% CIs in
the main text; coverage is implied through Type-I-rate calibration plus the
power figures.

### Reference figures / tables

- **Figure 1**: Type-I error across the four exposure/outcome combinations.
- **Figure 2**: Power (the headline "+66.3% relative power" figure).
- **Table 1**: Real-data HDL-C → preeclampsia across six estimators (3 common,
  3 MR-CARV).
- **Figure 3**: scatter plot, IVW vs MR-CARV(IVW).
- **Supplementary Section S5**: LD pruning sensitivity.
- **Supplementary Table S2**: LDL-C / TG / TC real-data results.

---

## 3. Real-data validation

### Datasets

- **Exposure (lipid traits)**: published WGS summary statistics from
  **Selvaraj et al. 2022**, $N \approx 66{,}000$. Traits: **HDL-C, LDL-C,
  triglycerides, total cholesterol**.
- **Outcome (preeclampsia)**: **nuMoM2b-HHS** (Nulliparous Pregnancy Outcomes
  Study: Monitoring Mothers-to-be Heart Health Study) — individual-level WGS,
  **486 preeclampsia cases vs 2,821 controls**. Accession: dbGaP
  **phs002808.v1.p1**.

### Gene panel

There is **no curated positive-control gene panel** (e.g., no PCSK9 / HMGCR /
LPA sanity check). Instead, the authors run the genome-wide STAARpipeline
burden test and let it surface significant rare-variant masks.

### Outcomes / traits

Single primary outcome: **preeclampsia (binary)**. Four exposure traits as
above; HDL-C is the headline.

### Positive / negative controls

- **No positive control** such as PCSK9 → LDL or HMGCR → LDL is reported.
- **No negative control** trait or gene is reported.
- Validation is therefore a single "MR-CARV finds something IVW barely misses"
  story rather than a battery of orthogonal sanity checks.

### Specific reported results

- Instruments retained: **63 common variants + 12 rare-variant burdens**
  (9 gene-coding, 1 gene-noncoding, 2 non-gene).
- HDL-C → preeclampsia (headline): MR-CARV(IVW) $\hat\beta = -0.020$,
  SE = 0.0102, **P = 0.047** (protective). Common-only IVW gave P = 0.0659.
  The qualitative claim is: rare variants tip the result over the 0.05 line.
- LDL-C, TG, TC: not significant after MR-CARV (Supplementary Table S2);
  presented as "no false positive" support.

---

## 4. Code / data availability

### GitHub

- **URL**: https://github.com/yu-zhang-oYo/mr.carv
- **Repo type**: a thin R package, NOT a paper-companion analysis repo.
- **Top-level layout**:
  ```
  mr.carv/
  ├── R/                  # function source
  ├── inst/extdata/       # one example workbook (TC.xlsx)
  ├── man/                # roxygen docs
  ├── DESCRIPTION
  ├── NAMESPACE
  ├── README.md
  └── mr.carv.Rproj
  ```
- **Files in `R/`**:
  `Burden_Effect.R`, `Gene_Centric_Coding_Burden.R`,
  `Gene_Centric_Coding_Burden_each.R`, `Gene_Centric_Noncoding_Burden.R`,
  `Gene_Centric_Noncoding_Burden_each.R`, `Individual_Estimate.R`,
  `Window_Burden.R`, `Window_Burden_each.R`, `data.R`,
  `select_LE_Estimate.R`, `select_LE_Estimate_ori.R`.
- **Simulation generator**: **not present in the repo.** No `sim/`,
  `simulation/`, `inst/sim/`, or vignette folder. The simulation code that
  produced Figures 1–2 is not shipped.
- **Real-data analysis driver**: **not present as a runnable script.** The
  README only walks through a toy workflow on `inst/extdata/TC.xlsx`. The
  actual nuMoM2b-HHS / Selvaraj pipeline is not in the repo.
- **Figure-generating notebook**: **none.** No Rmd, Quarto, or Jupyter
  artifact.
- **Toy / example inputs**: yes — `inst/extdata/TC.xlsx` contains four sheets
  (`indv_effect`, `coding`, `noncoding`, `window`) demonstrating the expected
  STAARpipeline-style input format.

### Data deposit

- **No Zenodo / Figshare / OSF DOI** is mentioned anywhere I could reach
  (README, landing page, IU ScholarWorks deposit).
- Exposure data: cited as "use the Selvaraj et al. 2022 summary stats."
- Outcome data: dbGaP **phs002808.v1.p1** (controlled access; nothing
  preprocessed is shipped).

### Reproducibility level

- **Low–medium.** Install path is `devtools::install_github(...)`; no Docker,
  no `renv.lock`, no conda env file, no Snakemake / Nextflow driver, no
  `Makefile`. Heavy upstream dependence on `STAARpipeline`, `SeqArray`,
  `gdsfmt`, `Matrix`, `MendelianRandomization`, `igraph` (all from
  CRAN/Bioconductor).
- A motivated user can rerun the headline analysis if they (a) gain dbGaP
  access to nuMoM2b-HHS, (b) recreate the STAARpipeline preprocessing for the
  Selvaraj exposure side themselves, and (c) write their own loop matching the
  paper's parameter grid. None of those steps are scaffolded in the repo.

---

## 5. What MR-CARV does NOT address that rvSMR will

1. **Weak-IV-robust inference.** MR-CARV's defense against weak instruments
   is to plug stronger burdens into dIVW / MR-RAPS, both of which assume
   "many weak instruments" with bounded $F$. There is **no
   Anderson–Rubin-style construction**, no test that stays valid as
   $F \to 0$, and no simulation that sweeps $F$ at all (let alone $F < 1$).
   This is the cleanest single gap for rvSMR to fill.

2. **Cell-type resolution.** MR-CARV operates at the **trait level** — HDL-C
   serum measurement → preeclampsia. There is no eQTL/sQTL handling, no
   cell-type-specific exposure stratification, and the burden-variant
   annotations are STAAR functional categories (coding / noncoding / window),
   not cell-state-resolved regulatory annotations. rvSMR's per-cell-type SMR
   layer is therefore genuinely new.

3. **Cauchy combination across MAF weights vs across orthogonal axes.**
   MR-CARV applies CCT *across MAF-weighting choices* (Beta(1,25), Beta(1,1))
   and *across annotation masks* (coding / noncoding / window) **within a
   single gene**. rvSMR explicitly does **not** do CCT across masks within a
   gene; instead it reserves CCT for combining genuinely orthogonal axes
   (e.g., across cell types or across independent tissues). The MR-CARV
   choice inflates correlation-induced ties between the combined p-values —
   STAAR-CCT mitigates this but does not eliminate it. rvSMR sidesteps the
   issue by construction.

4. **No positive/negative-control validation harness.** MR-CARV passes
   exactly one real-data analysis; it has no PCSK9-style positive control and
   no negative-control trait/gene. rvSMR can do strictly better with very
   modest additional engineering.

5. **No code for simulations.** Even if a reviewer asks rvSMR to "match
   MR-CARV's simulation," there is no upstream code to clone — the
   simulation has to be re-implemented from the paper's prose. This is an
   opportunity (we control the framing) and a risk (we must read the paper
   carefully to avoid an unfair head-to-head).

---

## 6. Top-3 takeaways for rvSMR validation

1. **Beat MR-CARV on its own simulation grid, then add a weak-IV axis on
   top.** Use the same 1000 Genomes IV construction (4 common + 4 rare
   burdens), same Beta(1,25)/Beta(1,1) weighting, same $\beta$ ladder
   (0.04 / 0.06 / 0.2 / 0.4), same N grid (1k / 2k / 10k), same 10,000-rep
   Type-I-error budget and 1,000-rep power budget — but **add a swept
   F-statistic axis (F = 50, 20, 10, 5, 2, 1, 0.5)**, which is the regime
   MR-CARV does not cover. Show that rvSMR's AR-based test holds nominal
   Type-I rate across that whole axis while dIVW / RAPS / MR-CARV(IVW)
   blow up.

2. **Add the comparators MR-CARV skipped.** MR-Egger, MR-PRESSO, RARE, and
   classical (Pavlides-Wang) SMR are all standard asks from a reviewer.
   Including them costs us little and lets us point directly at MR-CARV's
   missing rows.

3. **Build the positive/negative-control panel MR-CARV lacks.** A small set
   of well-established triples — PCSK9 → LDL → CAD, HMGCR → LDL → CAD,
   LPA → Lp(a) → CAD, plus a negative control such as
   "skin-pigmentation gene → CAD" — gives rvSMR a real-data validation
   story that MR-CARV cannot match, particularly because rvSMR resolves by
   cell type. Pull exposures from a public eQTL resource
   (e.g., OneK1K / eQTL Catalogue) and outcomes from a well-powered GWAS
   (e.g., CARDIoGRAMplusC4D or UKB CAD). This is the single change that
   most decisively differentiates rvSMR from MR-CARV at peer review.

---

## What I could NOT confirm after ~8 fetches

- I could not access the full text body of the published PDF (Oxford
  watermark redirect blocks WebFetch, the IU ScholarWorks PDF download
  returned only metadata, and medRxiv full-text returned 403). Specific
  numbers in §2 (replicate counts, 66.3% headline, β values, α weighting,
  CCT scope) come from the journal landing-page abstract+summary plus the
  README; I have **not** verified them against the body PDF myself. The user
  should spot-check before quoting in the paper.
- I could not confirm whether the body PDF reports a numeric F-statistic
  anywhere (e.g., as a *consequence* of the fixed $\alpha = 0.5$ choice).
  The simulation does not *sweep* F; whether it *reports* F is uncertain.
- I found **no Zenodo/Figshare/OSF DOI** and no archived snapshot of the
  paper-companion analysis code. If one exists, it is not advertised on the
  landing page or in the GitHub README.
