# rvSMR Validation Plan

*Drafted 2026-06-01. Synthesizes `research_RARE.md`, `research_MR-CARV.md`, `research_AR_cisMR_comparators.md` (Wang–Kang 2022, Patel–Lane–Burgess 2024, Yang 2025 cis-MRBEE) against HANDOVER §6 / §10 Task B.*

---

## Executive summary

Three things the comparator literature converges on; rvSMR must do all three:

1. **Coverage-vs-IV-strength figure** with 95% CIs across $F$ (or a concentration analog). The robust method holds a flat line at 0.95; Wald / IVW collapses at $F < 10$. This is the *single* most-cited figure in any weak-IV-robust MR paper (Wang–Kang Fig 6; Patel–Lane–Burgess Fig 2). rvSMR has no choice but to reproduce it.
2. **PCSK9-LDL plus a curated cis-MR panel** as the real-data sanity check (Yang 2025 covers PCSK9, APOA1, APOC1, ANGPTL3, CR1). MR-CARV / RARE both skip this; reviewers ask for it.
3. **Public GitHub + Zenodo + single-command figure pipeline**. Both MR-CARV and RARE ship thin R packages with NO simulation code, NO real-data driver, NO figure notebook. Beating that bar is cheap.

The binding constraint (HANDOVER §6): **real rare-variant exposure summary stats are not publicly downloadable** — TenK10K rare-variant Zenodo files are 214–260-byte placeholders; OneK1K's `rv_sign.txt` is filtered to significant hits only. This blocks the headline real-data analysis but does **not** block simulation, common-variant plumbing, or contacting Wei/Cuomo. The plan below runs four tracks in parallel; outreach unblocks Track 3 on its own schedule.

---

## 1. What comparator validation looks like (1-page cross-ref)

| | RARE (Cheng 2025) | MR-CARV (Zhang 2026) | Wang–Kang 2022 | Patel–Lane–Burgess 2024 | cis-MRBEE (Yang 2025) |
|---|---|---|---|---|---|
| **GitHub** | [Hide-in-lab/RARE](https://github.com/Hide-in-lab/RARE) | [yu-zhang-oYo/mr.carv](https://github.com/yu-zhang-oYo/mr.carv) | not found | [ash-res/mvmr-weakiv](https://github.com/ash-res/mvmr-weakiv) | [harryyiheyang/MRBEEX](https://github.com/harryyiheyang/MRBEEX) |
| **Sim genotype source** | synthetic + UKBB sanity | 1000G real haplotypes | Binomial(2, p) synthetic | mvn synthetic | cis-window 1000G LD |
| **# instruments** | 1000 common + 20 rare | 4 common + 4 burdens | 100 | 4 (baseline) | ~3 per exposure |
| **Sample size** | n = 50 000 fixed | n = 1k / 2k / 10k | n = 100 000 | n = 5 000 (+ subsamples) | n = 300 / 3k / 30k |
| **Weak-IV sweep** | **NONE** | **NONE** | $r$ (= $F$) ∈ {1, 4, 16, 25} | $F_{\min}$ ≈ 6, 10, 30 + $\xi$ | implicit via $n$ |
| **Type-I reps** | 300 | 10 000 | 1 000 | not specified | 500 |
| **Coverage figure** | absent | absent | **Fig 6** (canonical) | **Fig 2** (canonical) | implicit |
| **Comparators** | MV-IVW, MV-Egger, MV-LASSO, GRAPPLE | IVW, dIVW, MR-RAPS | MR-Egger, RAPS, WM, Q | IVW, GMM, AR, Kleibergen, Andrews | cis-MVIVW, PC-GMM, MRBEE-IPOD, TGFM, TGVIS |
| **Real-data positive control** | HDL→T2D (textbook) | HDL→preeclampsia (no PCSK9) | BMI→SBP (Zhao 2020 reprise) | MAP/PP→stroke + UKB subsample | **PCSK9, APOA1, APOC1, ANGPTL3, CR1** ✓ |
| **Negative controls** | none | none | none | none (subsample tactic instead) | none explicit |
| **Sim code shipped** | yes (`simulation/simulationl.R`) | **no** | no (no repo) | yes (R package) | yes (R package) |
| **Figure notebook** | no | no | no | partial | yes (Tutorial pdf) |
| **Zenodo deposit** | no | no | no | no | example data on Dropbox |

Detailed reports: `research_RARE.md` (§1–§6), `research_MR-CARV.md` (§1–§5), `research_AR_cisMR_comparators.md` (§A–§E).

---

## 2. The binding data constraint

| Source | What we need | Public status (verified 2026-05-27, HANDOVER §6) | Path |
|---|---|---|---|
| **TenK10K Phase 1 rare-variant SAIGE-QTL** | per-(gene × mask × cell type) $(\hat b_x, \mathrm{SE}_x)$ for ALL genes, sig + non-sig | Zenodo 17474113 rare-variant zips = **214–260 byte placeholders** | B1: contact Cuomo (Garvan); B2: contact Wei |
| **OneK1K SAIGE-QTL rare-variant** | same | Zenodo 10884040 `rv_sign.txt` (764 KB) = **filtered to significant only** | insufficient for AR (need full pair) |
| **TenK10K Phase 1 common-variant cis-eQTL** | per-cell-type common-variant summary stats | Zenodo 17474113, **real, 14–23 GB / file** | usable now (Track 2) |
| **Genebass burden** | per-gene burden $(\hat b_y, \mathrm{SE}_y)$ for outcomes | GCS requester-pays Hail MT `gs://ukbb-exome-public/500k/results/results.mt`; needs `gsutil -u <project>` | usable with authentication |
| **UKB-PPP rare-variant pQTL** (Dhindsa 2023) | pQTL anchor for §8 annotation-class normalization | public at azphewas.com + Dhindsa browser | usable now |
| **TenK10K individual-level** | for re-running SAIGE-QTL ourselves | EGA EGAS50000001653, **gated** | B3 (high compute, low priority for v1) |

**Implication**: Track 3 (real rare-variant exposure) is on Wei/Cuomo's schedule. Tracks 1, 2, 4, 5 below are unblocked.

---

## 3. Validation tracks (four parallel)

### Track 1 — Simulation harness with a swept $F$ axis (**unblocked, primary**)

**Goal**: produce the canonical coverage-vs-$F$ figure (Wang–Kang Fig 6 / Patel–Lane–Burgess Fig 2 style) showing rvSMR-AR holds at 95% across $F \in \{0.5, 1, 2, 5, 10, 30, 100\}$ while Wald / IVW / dIVW / MR-RAPS collapse below $F < 10$.

**Generative model** — borrow MR-CARV's 1000G-anchored setup, add rvSMR's mask structure + an explicit $F$ knob:

```
genotype  : 1000G real haplotypes; pick one chromosome per mask (suppresses cross-mask LD)
masks     : K = 3 per gene (pLoF, mis:LC, reg); m = 50 variants per mask
MAF       : empirical 1000G MAF in 0.001-0.05; tag rare; STAAR Beta(1, 25) weights
exposure  : X_i = sum over j of alpha_j G_{ij} + U_i + eps_i; calibrate alpha_j to hit target F
outcome   : Y_i = beta X_i + theta * Z_i + U_i + eps_i; theta = pleiotropy mix
sample    : n_x = n_y in {1000, 10000, 50000}; cell-type stratified via subsample factor
overlap   : R_xy fixed at 0 (two-sample); add R_xy = 0.3 sensitivity cell
```

**Swept axes** (the central thing MR-CARV and RARE both skip):

| Axis | Values | Why |
|---|---|---|
| First-stage $F$ | {0.5, 1, 2, 5, 10, 30, 100} | canonical weak-IV figure |
| true $\beta$ | {0, 0.04, 0.1, 0.2, 0.4} | size + power at non-zero nulls (Wang–Kang Fig 1 template) |
| K masks | {1, 2, 3, 5} | identification only at $K \ge 2$; over-id at $K \ge 3$ |
| Pleiotropy fraction | {0%, 10%, 30%, 50%} | annotation-concordance Q stress test |
| Coherent pleiotropy fraction | {0%, 100%} | reveal the HEIDI-rv blind spot (HANDOVER §7 caveat) |
| Sample overlap $R_{xy}$ | {0, 0.3} | sample-overlap correction sanity |
| Cell-type stratification $c$ | {1, 5, 28} | cell-type-Q stress |

**Replicates**: 1 000 per cell for power; 10 000 per cell for Type-I error (MR-CARV's budget).

**Comparators** (cover the rows MR-CARV skipped):
- Wald (delta-method CI) — must collapse at low $F$
- IVW, dIVW (mr.carv reimplementation) — same
- MR-RAPS — same
- MR-Egger — pleiotropy axis
- MR-PRESSO — outlier axis
- RARE (Bayesian) — credible-interval comparison
- classical SMR (Zhu 2016) — common-variant baseline
- **rvSMR-AR (K=1 closed-form, mrAR.R)** — should hold flat
- **rvSMR-AR (K≥2 grid+uniroot, mrAR_multi.R)** — should hold flat + Sargan-J

**Metrics**:
- Empirical 95% CI coverage (the headline plot — flat line for rvSMR-AR)
- Type-I error at $\beta_0 \in \{0, 0.04, 0.1, 0.2, 0.4\}$ (Wang–Kang multi-null)
- Power at true $\beta = 0.04, 0.1, 0.2$
- CI shape distribution: % bounded / disconnected / whole-line / empty across $F$
- Sargan-J Type-I rate under no pleiotropy; power under 10/30/50% pleiotropy
- HEIDI-rv coverage under single-variant vs coherent pleiotropy (highlights the $\mathcal O(1/m)$ caveat honestly)
- Annotation-class Q calibration with / without pQTL-anchor normalization
- $R^2_{Z\to X}$ and $RV$ distributions across $F$

**Decision rule for "rvSMR-AR passes"**: empirical 95% coverage $\in [0.93, 0.97]$ across all $F$ cells (Patel–Lane–Burgess tolerance).

**Deliverables**:
- `inst/sim/generate.R` — generator parameterized by the table above
- `inst/extdata/sim_reference.rds` — fixed reference dataset (regression seed)
- `tests/testthat/test-coverage_vs_F.R` — coverage gate that fails CI if rvSMR-AR drops below 0.93 anywhere
- `figures/coverage_vs_F.qmd` — single-command Quarto regenerating the headline figure
- Zenodo deposit with the full CSV grid + DOI

**Effort estimate**: 2–3 weeks. Generator + 1000G LD setup ≈ 1 week; sweep + comparator wrappers ≈ 1 week; figures ≈ few days.

---

### Track 2 — Common-variant SMR plumbing test (**unblocked**)

**Goal**: end-to-end run on REAL data sizes (14–23 GB Zenodo files) before Track 3 unblocks. Proves the pipe handles real LD / real annotation / real I/O. NOT a methodology test — purely engineering.

**Substrate**:
- Exposure: TenK10K Phase 1 common-variant cis-eQTL (Zenodo 17474113, real public files)
- Outcome: Genebass single-variant outputs (`gs://ukbb-exome-public/500k/results/results.mt`, requester-pays)
- Sanity gene: PCSK9 → LDL-C (canonical positive control; should reproduce direction and rough magnitude from RCT)

**Test**:
- Single end-to-end `mrAR_multi()` call on the PCSK9 common-variant cis-eQTL × Genebass LDL-C burden
- Report $F$ distribution, CI shape distribution, $J$ p-value
- 1-page write-up: "we can do this on real data sizes"

**Effort**: 1 week including GCS authentication setup.

---

### Track 3 — Real rare-variant exposure (**gated on Wei / Cuomo**)

**Headline panel** (HANDOVER §1, §5): 5 RCT-validated genes × 5 outcomes.

| Gene | Outcome | RCT validation |
|---|---|---|
| PCSK9 | LDL-C, CHD | evolocumab / alirocumab (FOURIER / ODYSSEY) |
| ANGPTL3 | LDL-C, TG | evinacumab (ELIPSE) |
| APOC3 | TG | volanesorsen / olezarsen |
| HMGCR | LDL-C, CHD | statins (4S, JUPITER) |
| LPA | Lp(a), CHD | olpasiran / pelacarsen (in-trial) |

For each (gene × outcome): run rvSMR with $K = 3$ masks (pLoF / mis:LC / reg) × $c$ PBMC cell types (TenK10K Phase 1, 28 types). Expected result: bounded AR CI, non-empty, sign matching RCT direction, $J$ p-value > 0.05, annotation-class Q > 0.05.

**Negative-control panel** (rvSMR's differentiator over RARE / MR-CARV): biologically implausible (gene × outcome) pairs from the same data. Candidates: a skin pigmentation gene × CHD; a Y-chromosome gene × osteoporosis; ~5 pairs total. Expected: whole-line CI or wide bounded CI with no significant signal.

**Gating**: requires Wei / Cuomo to populate Zenodo 17474113 OR share internal SAIGE-QTL outputs.

**Outreach drafts** to produce (HANDOVER §12 (c)):
- email to Cuomo (Garvan): narrow scope — the 214-byte placeholder issue + one-paragraph rvSMR pitch
- email to Wei: coordination framing — internal SAIGE-QTL outputs vs re-run on OneK1K/TenK10K

---

### Track 4 — Sub-sampling-induced weak-IV + Goodhart (**unblocked once Track 3 returns ANY real data, but can prototype on common-variant Track 2**)

Two tactics borrowed from Patel–Lane–Burgess 2024:

**(a) Sub-sampling (Fig 5)**: take the largest available rare-variant burden GWAS (Genebass, $N \approx 394$k); sub-sample to $n \in \{2k, 10k, 50k, 100k\}$. For each sub-sample size: 100 random sub-samples → 100 rvSMR-AR CIs. Heatmap of CI shapes. The truth is fixed (= full-sample point estimate); rvSMR-AR should cover the full-sample truth at nominal rate while sub-sampled Wald should collapse.

**(b) Goodhart / selective-reporting (Fig 3)**: run all sub-samples; report coverage two ways — (i) all sub-samples, (ii) only sub-samples where $F > 10$. The (ii) bar shows what publication bias does to coverage if practitioners screen on Stock–Yogo. rvSMR-AR's argument is precisely "you don't need to screen".

Both tactics produce single-pane figures that comparator papers do not have. High impact at peer review.

---

### Track 5 — Annotation-class concordance pilot using publicly available pQTL anchors (**unblocked**)

For each of the 5 RCT genes: pull burden-on-protein from Dhindsa 2023 (UKB-PPP rare-variant, public). Normalize class-specific burden-on-outcome by this anchor (§8 of the algorithm walkthrough). Run Cochran-Q across the three mask classes. Expected: $Q$ p-value > 0.05 (concordance) for the validated genes; report what mask class dominates.

Useful as a sub-figure even before Track 3 unblocks (uses publicly available outcome data + pQTL anchors).

---

## 4. Validation tactics rvSMR can SKIP (defensible non-coverage)

- **Bayesian credible-interval calibration** (RARE-style) — rvSMR is frequentist; coverage is the right metric.
- **MVMR** — rvSMR is univariable (one exposure gene at a time); MVMR is a future extension.
- **Multi-gene cis-window** (cis-MRBEE-style) — rvSMR's unit is the single gene with $K \ge 3$ masks, not a multi-gene region.
- **SuSiE fine-mapping inside the burden** — rare-variant burden has few variants per mask (median ~10–50); fine-mapping is a common-variant convention. Flag this as a *deliberate choice* in the paper; cite cis-MRBEE as the precedent that doesn't apply.

---

## 5. Recommended file structure (extends current `rvMR` R package)

```
rvMR/
├── R/                          (existing — wald_burden, mrAR, mrAR_multi, ...)
├── tests/testthat/             (existing — 68 passing assertions)
├── inst/
│   ├── sim/
│   │   ├── generate.R          NEW: Track 1 generator
│   │   ├── comparators.R       NEW: wrappers for IVW/dIVW/RAPS/Egger/PRESSO/RARE/SMR
│   │   ├── sweep_config.R      NEW: the parameter grid above
│   │   └── README.md           NEW: how to reproduce each figure
│   ├── extdata/
│   │   ├── sim_reference.rds   NEW: fixed reference dataset
│   │   ├── 1000G_LD/           NEW: precomputed per-chromosome LD blocks
│   │   └── pcsk9_smoke.rds     NEW: Track 2 PCSK9 common-variant inputs
│   └── analysis/
│       ├── track2_pcsk9.R      NEW: Track 2 driver
│       ├── track3_rct_panel.R  NEW: Track 3 driver (placeholder until Wei/Cuomo)
│       ├── track4_subsample.R  NEW: Track 4 driver
│       └── track5_pqtl_anchor.R NEW: Track 5 driver
└── figures/
    ├── coverage_vs_F.qmd       NEW: headline figure
    ├── type1_vs_beta0.qmd      NEW: Wang–Kang Fig 1 analog
    ├── ci_shape_dist.qmd       NEW: 4-CI-shape composition vs F
    ├── pleiotropy_J.qmd        NEW: Sargan-J calibration + power
    ├── heidi_rv_oneOverM.qmd   NEW: HEIDI-rv power vs m (honest disclosure)
    ├── annotation_Q.qmd        NEW: with vs without pQTL anchor
    ├── subsample_heatmap.qmd   NEW: Patel-Lane-Burgess Fig 5 analog
    └── goodhart.qmd            NEW: Patel-Lane-Burgess Fig 3 analog
```

Single command `make figures` (or `quarto render figures/`) regenerates every panel. Outputs land in a versioned `outputs/` directory; a `make zenodo` packages outputs + sim CSVs into a Zenodo-ready tarball.

---

## 6. Concrete sequencing (next 8 weeks)

| Week | Track 1 sim | Track 2 plumbing | Track 3 real data | Track 4 sub-sample | Track 5 pQTL anchor |
|---|---|---|---|---|---|
| 1 | 1000G LD blocks; generator skeleton | GCS auth + Genebass cost estimate | draft Cuomo + Wei emails | — | Dhindsa 2023 anchor pull |
| 2 | comparators (IVW/dIVW/RAPS) | PCSK9 cis-eQTL pull | send emails (after Francis review) | — | normalize PCSK9 / ANGPTL3 / APOC3 / HMGCR / LPA |
| 3 | comparators (Egger/PRESSO/RARE/SMR) | PCSK9 × LDL end-to-end | — | — | Cochran-Q per gene |
| 4 | coverage-vs-F sweep (1 000 reps) | Track 2 write-up | follow up | sub-sample setup on Genebass | Track 5 figure |
| 5 | Type-I + power sweep (10 000 reps) | — | — | sub-sample sweep | — |
| 6 | CI-shape composition + Sargan-J | — | (response-dependent) | Goodhart cell | — |
| 7 | HEIDI-rv + annotation-Q figures | — | — | sub-sample heatmap | — |
| 8 | Zenodo bundle + figure pipeline | — | — | — | — |

Tracks 1, 4, 5 run in parallel. Track 2 is sequential plumbing (1 week). Track 3 is on outreach timeline.

---

## 7. Open decisions for Francis

1. **STAARpipeline vs internal mask definition.** If we reuse STAAR's pLoF / mis:LC / reg, we need the variant annotation pipeline to match Genebass and TenK10K masks exactly (HANDOVER §6.2 Wei ask #2). Decide: rely on STAARpipeline output, or write our own crosswalk.
2. **1000G ancestry**: EUR only (matches UKB / Genebass / TenK10K majority) or admixed (more honest but more expensive)? Comparator papers all use EUR.
3. **Negative-control panel composition**: how many pairs? Which genes? Suggest 5 pairs drawn from the same Genebass outcome set so the comparison is apples-to-apples.
4. **Zenodo bundle policy**: deposit raw simulation CSVs (10s of GB potentially) or only the figure-input summaries (~few MB)?
5. **Comparator effort vs novelty**: which of RARE / MR-CARV / cis-MRBEE do we benchmark against in the main figure vs supplement? Suggest IVW + dIVW + RAPS + Wald in main; RARE + MR-CARV + classical SMR + Egger + PRESSO in supplement.

---

## 8. References (URLs and key file paths)

### Papers + GitHub
- RARE: Cheng et al. 2025 *Brief Bioinform* 26(3):bbaf214 — [paper](https://academic.oup.com/bib/article/26/3/bbaf214/8131742) — [GitHub Hide-in-lab/RARE](https://github.com/Hide-in-lab/RARE) — sim driver `simulation/simulationl.R`
- MR-CARV: Zhang et al. 2026 *Brief Bioinform* 27(1):bbaf649 — [paper](https://academic.oup.com/bib/article/27/1/bbaf649/8416446) — [GitHub yu-zhang-oYo/mr.carv](https://github.com/yu-zhang-oYo/mr.carv) — no sim code shipped
- Wang–Kang 2022 *Biometrics* 78(4):1699 — [paper](https://academic.oup.com/biometrics/article/78/4/1699/7460098) — [arXiv 1909.06950](https://arxiv.org/abs/1909.06950) — no GitHub
- Patel–Lane–Burgess 2024 — [arXiv 2408.09868](https://arxiv.org/abs/2408.09868) — [GitHub ash-res/mvmr-weakiv](https://github.com/ash-res/mvmr-weakiv) — title is "Weak instruments in multivariable MR" (NOT "AR tests" as the HANDOVER says; flag this for the paper bibliography)
- cis-MRBEE: Yang et al. 2025 *Brief Bioinform* 26(3):bbaf250 — [paper PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC12140020/) — [GitHub harryyiheyang/MRBEEX](https://github.com/harryyiheyang/MRBEEX)

### Data sources
- TenK10K Phase 1 (common-variant): [Zenodo 17474113](https://zenodo.org/records/17474113) — 14–23 GB / file, real
- TenK10K rare-variant: same Zenodo record but the rare-variant zips are 214–260 byte placeholders (verified)
- OneK1K SAIGE-QTL rare-variant: [Zenodo 10884040](https://zenodo.org/records/10884040) `rv_sign.txt`, 764 KB, filtered-to-significant only
- Genebass: [gs://ukbb-exome-public/500k/results/results.mt](https://app.genebass.org/) — requester-pays GCS Hail MT
- UKB-PPP rare-variant pQTL: [azphewas.com](https://azphewas.com/) + [Dhindsa pQTL browser](https://astrazeneca-cgr-publications.github.io/pqtl-browser/)
- 1000G Phase 3 haplotypes: [1000genomes.org](https://www.internationalgenome.org/data) (for simulation LD reference)

### Cross-references
- Source reports in this repo: `research_RARE.md`, `research_MR-CARV.md`, `research_AR_cisMR_comparators.md`
- Algorithm walkthrough: `main.tex` (each Step's "Code reference" subsection points to the rvMR R function we'd extend)
- Project handover: `/home/francisfenglu4/rvSMR/May_30md/HANDOVER_2026-05-27.md` §6, §10 Task B
