# FinnGen × MVP × UKBB joint meta-analysis — research report for rvSMR validation

*Drafted 2026-06-08 for Francis Lu (Dana-Farber / Sethi lab) — evaluating the public PheWeb at https://public-mvp-ukbb.finngen.fi/ as an outcome-side ($\hat\beta_y$) data source for rare-variant burden MR.*

---

## TL;DR (≤200 words)

**Headline: NOT usable as a rare-variant burden outcome source.** The FinnGen × MVP × UKBB joint meta-analysis is a **single-variant GWAS meta-analysis only**, performed with SAIGE + fixed-effects inverse-variance weighting (IVW) — there are no gene-level burden tests, no pLoF / missense:LC / regulatory masks, no STAARpipeline or SAIGE-GENE+ output. The FinnGen Handbook is explicit: "FinnGen does not run gene-level burden tests such as SKAT-O" and points users to Genebass for burden.

**Key access detail:** Public PheWeb at `https://public-mvp-ukbb.finngen.fi/` (Release 12 = 330 binary phenotypes) is open access, no login. The richer consortium-only version `https://mvp-ukbb.finngen.fi/` (2,003 traits including quantitative lipids LDL-C, ApoB, TC, TG) requires FinnGen account. Per-phenotype TSVs are downloadable from a GCS bucket free of charge; consortium-only version sits in the FinnGen sandbox under `/finngen/library-green/finngen_R12/finngen_R12_analysis_data/meta_analysis/mvp_ukbb/`.

**Fit with VALIDATION_PLAN.md:** Does **not** unblock Track 3 (real rare-variant outcome data). Does **not** supersede Genebass — Genebass remains the only public source of $\hat\beta_{burden,y}$ at scale. The FG-MVP-UKBB resource is **useful for Track 2 (common-variant plumbing)** as a sanity-check on Genebass single-variant outputs at 3× the EUR sample size + meaningful AFR component.

---

## 1. Project metadata

| Field | Value |
|---|---|
| Full name | "FinnGen + MVP + UKBB three-way GWAS meta-analysis" (browser title: "finngen r13 + mvp + ukbb") |
| Cohorts | FinnGen, Million Veteran Program (MVP), UK Biobank (UKBB) |
| Current release | **R12** for the documentation-stable cut (330 binary phenotypes, public); **R13** label appears on the root browser homepage as the current internal cut |
| Public browser | https://public-mvp-ukbb.finngen.fi/ (no login) |
| Consortium browser | https://mvp-ukbb.finngen.fi/ (FinnGen account required) |
| Date of release | R12: 2025 (FinnGen released R12 with the MVP/UKBB meta-analysis as a new addition); R13: 2026 internal |
| File format | TSV, bgzip-compressed, tabix-indexed |
| Citations | **Primary methods/discovery papers**: (1) Pereira, Cho, Gaziano et al. (PMC12799792), "Leveraging large-scale biobanks for therapeutic target discovery" — describes the 2,003-trait FinnGen + MVP + UKBB meta on >1.2M individuals; (2) Ricci et al. 2025 *Orphanet J Rare Dis* — sarcoidosis exemplar using the FG-MVP-UKBB meta browser, DOI 10.1186/s13023-025-04097-1 (PMC12751437). |
| FinnGen flagship | Kurki et al. 2023, *Nature* 613:508–518, DOI 10.1038/s41586-022-05473-8 |
| Senior investigators | J. Michael Gaziano (MVP), Kelly Cho (MVP/VA Boston), Alexandre C. Pereira; FinnGen leadership (Aarno Palotie, Mark Daly, Mari Niemi); MVP leadership through VA Office of Research |

## 2. Sample composition

| Cohort | N (approx.) | Ancestries |
|---|---|---|
| FinnGen R12 | ~520,000 | Finnish (founder-enriched EUR) |
| UK Biobank | ~450,000 (EUR only used in the meta) + small AFR | EUR (primary); small AFR available |
| MVP | ~570,000 | EUR + AFR + HISP (HARE-stratified; the Pereira paper reports MVP "race and ethnicity stratified using HARE") |
| **Total** | **>1.2 million** (~1.5 million reported in some sources) | EUR-dominant with meaningful AFR contribution; no EAS in the trans-cohort meta (confirmed in Ricci 2025) |

**Sarcoidosis paper exemplar sample sizes** (verifies ancestry decomposition):
- FinnGen: 5,411 cases / 492,311 controls (Finnish)
- MVP EUR: 1,509 / 449,523
- MVP AFR: 1,827 / 119,828
- UKBB EUR: 639 / 419,892
- UKBB AFR: 53 / 6,583

**Phenotypes:**
- **Public release**: 330 binary disease endpoints, organized by ICD-10-like FinnGen endpoint codes (E4_*, I9_*, C3_*, etc.). Categories: infectious, neoplasms, blood/immune, endocrine/metabolic, mental, neuro, eye, ear, GI, cardiac, rheumatologic.
- **Consortium-only release (mvp-ukbb.finngen.fi)**: 2,003 harmonized traits including a "Laboratory/Quantitative Measures" category (category index 100) with hemoglobin, cholesterol, ferritin, HbA1c, ALT/AST, glucose, sodium, potassium, creatinine, **triglycerides**, CRP, plus Apolipo_A, Apolipo_B, HDLC, **LDLC**, TotChol (confirmed in Pereira/Cho paper).
- **Lipid biomarkers — public version**: only ICD-coded **disease** endpoints (E4_LIPOPROT "Disorders of lipoprotein metabolism", E4_HYPERCHOL "Pure hypercholesterolaemia"). The continuous quantitative biomarkers (LDL-C, ApoB, TC, TG as labs) are **only in the consortium-only PheWeb**, not the public one.

## 3. Data type & format — THE CRITICAL QUESTION

| Question | Answer |
|---|---|
| Rare-variant burden tests (gene-level $\hat\beta_{burden,y}, \mathrm{SE}_{burden,y}$)? | **NO.** FinnGen Handbook (`docs.finngen.fi/where-to-begin.../im-interested-in-finngen-rare-variant-phenotypes`) states explicitly: "FinnGen does not run gene-level burden tests such as SKAT-O" — and directs users to **Genebass** for burden. The FG-MVP-UKBB meta inherits this: it is single-variant only. |
| Single-SNP GWAS? | **YES — exclusively.** SAIGE mixed-model per cohort, then **fixed-effects inverse-variance weighted** meta with Cochran's Q heterogeneity. |
| Mask definitions (pLoF / mis:LC / reg)? | **N/A — no burden tests.** |
| MAF cutoff for "rare" | Standard GWAS MAF filter (FinnGen typically MAF > 0.001 after imputation); per-variant, not per-mask. |
| File format | TSV bgzip + tabix index. Columns (21): `#CHR, POS, REF, ALT, SNP, rsid`, then per-cohort `beta / SE / pval / AF (overall, cases, controls)`, then meta `N, IVW beta/SE/pval, -log10p, Cochran's Q heterogeneity p`. |
| Bulk download | Per-phenotype TSVs available on GCS. Per the FinnGen "access results" page, summary statistics from GCS are **free of charge**. The exact public bucket path for the FG-MVP-UKBB meta is not in the public docs (sandbox path is `/finngen/library-green/finngen_R12/finngen_R12_analysis_data/meta_analysis/mvp_ukbb/`); could not confirm a public GCS URL for the meta (vs the FinnGen-only sumstats which are at `gs://finngen-public-data-r12/`). **Flag: needs follow-up** to confirm whether the public mvp-ukbb meta TSVs are downloadable outside the sandbox. |
| API for single-gene query | `https://public-mvp-ukbb.finngen.fi/api/gene_phenos/<GENE>` returns top-SNP-per-phenotype JSON. Verified working for PCSK9 and HMGCR (see §7). |

## 4. Access model

| Aspect | Public PheWeb | Consortium PheWeb |
|---|---|---|
| URL | https://public-mvp-ukbb.finngen.fi/ | https://mvp-ukbb.finngen.fi/ |
| Login | None | FinnGen account |
| Application | None | FinnGen consortium membership |
| Cost | Free | Free (membership-gated) |
| IP restriction | None observed | Sandbox-style |
| Terms | FinnGen Data Access policy — summary statistics publicly redistributable with attribution; check FinnGen.fi/en/access_results for exact license. **Flag: confirm publication terms before final manuscript.** |
| API | REST endpoints under `/api/` (phenos, gene_phenos, region, manhattan); could not confirm a stable rate-limited public API but the gene_phenos endpoint returned JSON cleanly in this audit. |

## 5. Comparison to Genebass

| Feature | Genebass (Karczewski 2022) | FG-MVP-UKBB |
|---|---|---|
| Variant resolution | **Gene-level burden + SKAT-O + single-variant** | **Single-variant only** |
| N (effective) | ~394k UKBB exomes | >1.2M trans-cohort (~3× larger) |
| Phenotypes | 4,529 across ICD/quantitative | 330 public binary / 2,003 consortium |
| Ancestry | ~95% EUR | EUR-dominant + meaningful AFR (MVP) + Finnish-enriched (FinnGen); no EAS |
| Mask definitions | pLoF, missense, synonymous; MAF cutoffs 0.001, 0.0001, singleton | N/A |
| File format | Hail MatrixTable, requester-pays GCS | TSV + bgzip + tabix |
| Cost to download | gsutil egress (~$0.12/GB) — small for targeted gene queries | Free |
| Burden $\hat\beta_y$ available | **YES** | **NO** |

**Key takeaway:** Genebass remains the **only public source** of per-gene per-mask $\hat\beta_{burden,y}, \mathrm{SE}_{burden,y}$ at biobank scale. FG-MVP-UKBB does NOT supersede it for our rvSMR use case — but it **complements** Genebass for the **common-variant comparator track** (Track 2 of VALIDATION_PLAN.md) at 3× the sample size with the AFR/Finnish-enriched ancestry mix.

## 6. Fit with rvSMR pipeline

| Question | Answer |
|---|---|
| Can we extract $\hat\beta_{burden,y}$ for our 5 RCT genes × 5 outcomes? | **NO.** No burden tests in this resource. Only top-SNP-per-phenotype hits. |
| Is the mask definition compatible with SAIGE-QTL $\hat\beta_x$? | N/A — no masks. |
| Cross-cohort heterogeneity bias? | The trans-cohort meta uses Cochran's Q per variant; if we were using this as common-variant input, we'd want to filter on Q p > 0.05 or restrict to a single ancestry stratum. Finnish founder enrichment vs UKBB common-allele backbone is a documented source of heterogeneity — flagged but solvable via the per-variant het-p column. |

**Net: this is the wrong granularity for rvSMR's outcome side.** rvSMR's outcome anchor must be gene-level burden; this resource is variant-level.

## 7. Specific URLs for the 5 RCT genes

Browser links (work in any browser, no login):

| Gene | Public browser URL | Verified |
|---|---|---|
| PCSK9 | https://public-mvp-ukbb.finngen.fi/gene/PCSK9 | yes — see below |
| HMGCR | https://public-mvp-ukbb.finngen.fi/gene/HMGCR | yes — see below |
| ANGPTL3 | https://public-mvp-ukbb.finngen.fi/gene/ANGPTL3 | not fetched (API consistent) |
| APOC3 | https://public-mvp-ukbb.finngen.fi/gene/APOC3 | not fetched |
| LPA | https://public-mvp-ukbb.finngen.fi/gene/LPA | not fetched |

**PCSK9 top-SNP-per-phenotype on the public browser** (top variant rs11591147, 1:55039974:G:T, missense p.Arg46Leu):
- E4_LIPOPROT (Disorders of lipoprotein metabolism, ICD): β = -0.473, p = 1.29×10⁻³⁰⁶
- E4_HYPERCHOL (Pure hypercholesterolaemia, ICD): β = -0.405, p = 1.58×10⁻¹⁶⁷
- I9_IHD (Ischaemic heart disease): β = -0.21, p = 6.76×10⁻⁷⁰
- I9_CORATHER (Coronary atherosclerosis): β = -0.226, p = 2.24×10⁻⁶¹
- I9_ANGINA, I9_MI_STRICT, I9_UAP also present
- **LDL-C as a continuous trait: not in public version** (consortium-only)

**HMGCR top-SNP-per-phenotype**:
- E4_LIPOPROT: β = 0.0742, p = 3.39×10⁻¹¹¹ (rs12916, 3' UTR — the classic statin-target proxy)
- E4_HYPERCHOL: β = 0.059, p = 6.92×10⁻⁶¹
- I9_ABAORTANEUR: β = 0.06, p = 6.31×10⁻¹¹
- I9_AORTANEUR: β = 0.0461, p = 5.75×10⁻¹⁰

**Sign and magnitude of PCSK9-rs11591147 → IHD (β ≈ -0.21, OR ≈ 0.81) matches the FOURIER / ODYSSEY RCT direction** — useful as a directional sanity check for common-variant Track 2 even though it is not the rvSMR target unit.

## 8. Verdict for rvSMR validation plan

1. **Does this unblock Track 3 (real rare-variant outcome data) of VALIDATION_PLAN.md?** No. Track 3 requires gene-level burden $\hat\beta_{burden,y}$ with mask stratification; this resource is single-variant GWAS only and the FinnGen Handbook explicitly says FinnGen does not run burden tests.
2. **Is this strictly better than Genebass, or complementary?** Strictly complementary, **not** better. Genebass is the only public source for $\hat\beta_{burden,y}$ at biobank scale; FG-MVP-UKBB is an attractive **common-variant** companion for Track 2 (3× the EUR effective N, plus AFR) but cannot substitute for Genebass on the burden side.

**Items I could not confirm — flagged for follow-up:**
- Whether the public FG-MVP-UKBB meta TSVs are downloadable from a public GCS bucket (vs sandbox-only). The FinnGen "Access results" page promises free GCS downloads for FinnGen sumstats; the analogous bucket for the **3-way meta** is not surfaced in the docs I read. Email finngen-helpdesk@helsinki.fi or check the FinnGen Handbook "Meta-analysis PheWebs" page for the exact public download path.
- Whether the consortium-only quantitative lipid traits (LDL-C, ApoB) are downloadable to non-FinnGen affiliates with a data access application (vs. sandbox-only viewing).
- Exact R13 release date (root browser shows "finngen r13 + mvp + ukbb"); R12 is the documented cut.

## 9. Suggested next step

**This resource is NOT useful for Track 3 burden outcomes — proceed to Genebass for $\hat\beta_{burden,y}$.**

For Track 2 (common-variant plumbing), the FG-MVP-UKBB meta is a strong companion. Concrete 3-bullet action plan if we want to use it as a common-variant outcome benchmark:

1. **Pull per-variant trans-cohort meta TSVs** for E4_LIPOPROT, E4_HYPERCHOL, I9_IHD, I9_CORATHER, I9_MI_STRICT from the public browser (or contact FinnGen for the GCS path). Use these as binary disease-outcome anchors for the 5 RCT genes alongside the Genebass single-variant arm. Note: continuous LDL-C/ApoB will require either Genebass quantitative or, if Francis has FinnGen account access, the consortium PheWeb.
2. **Cross-check PCSK9 rs11591147 → IHD direction and effect-size** against the Genebass single-variant arm to sanity-check Track 2 plumbing end-to-end on a 3× larger sample. Document the Cochran's Q heterogeneity for the Finnish vs UKBB-EUR cohorts — useful for the manuscript's "multi-ancestry robustness" sub-figure.
3. **Do NOT cite this resource as the rvSMR outcome side in the headline panel.** The rvSMR algorithm's outcome anchor is gene-level burden; cite Genebass for the burden outcome and cite FG-MVP-UKBB only where common-variant cross-validation is shown.

If we want to **expand burden outcomes beyond Genebass**, the realistic options remain: (a) AstraZeneca PheWAS portal (azphewas.com) — UKB-PPP/470k exomes, public burden tests; (b) UK Biobank Exome Browser; (c) the in-prep FinnGen rare-variant releases (re-check this resource in 6 months as R13/R14 may add burden output — FinnGen has been hinting at it).

---

## References

- Public PheWeb: https://public-mvp-ukbb.finngen.fi/
- Consortium PheWeb: https://mvp-ukbb.finngen.fi/
- FinnGen Handbook on burden tests: https://docs.finngen.fi/where-to-begin.../im-interested-in-finngen-rare-variant-phenotypes (states "FinnGen does not run gene-level burden tests such as SKAT-O")
- FinnGen Handbook meta-analysis PheWebs index: https://docs.finngen.fi/working-outside-the-sandbox/meta-analysis-phewebs
- FinnGen Handbook file-format spec: https://docs.finngen.fi/finngen-data-specifics/green-library-data-aggregate-data/core-analysis-results-files/ukbb-finngen-meta-analysis-file-formats
- Pereira/Cho/Gaziano "Leveraging large-scale biobanks for therapeutic target discovery": https://pmc.ncbi.nlm.nih.gov/articles/PMC12799792/
- Ricci 2025 sarcoidosis exemplar (DOI 10.1186/s13023-025-04097-1): https://pmc.ncbi.nlm.nih.gov/articles/PMC12751437/
- FinnGen flagship (Kurki 2023, *Nature* 613:508): https://www.nature.com/articles/s41586-022-05473-8
- Cross-reference inside this repo: `/home/francisfenglu4/projects/rvSMR_Math/VALIDATION_PLAN.md` §2 (data constraint table) and §3 Track 2 (common-variant plumbing).
