# FinnGen PheWeb + Open Targets pull log

Date: 2026-06-08
Endpoint host: `https://mvp-ukbb.finngen.fi` (the public FinnGen R12 × MVP × UKBB meta-PheWeb; `public-mvp-ukbb.finngen.fi` redirects here).
- Release: FinnGen R12, ~500 K Finns; MVP three subcohorts (EUR/AFR/AMR); UKBB; matched binary endpoints (330 total).
- ICD-10 derived binary endpoints; lab phenotypes (LOINC codes 30xxxxx) available too.

## API endpoints used

| Purpose | URL pattern | Notes |
|---|---|---|
| Variant autocomplete | `/api/autocomplete?query=<rsid>` | returns rsID → chr-pos-ref-alt map |
| Variant page (all phenotypes) | `/api/variant/<chr>-<pos>-<ref>-<alt>` | returns JSON `results` array indexed by phenocode |
| Manhattan plot (all sig variants for a pheno) | `/api/manhattan/pheno/<pheno>` | gzip-compressed |

## Outcome (b_y) — FinnGen meta-PheWeb pulls

For each lead variant we hit the variant API and extracted the `I9_IHD` row's `beta` and `sebeta` (these are the all-5-cohort inverse-variance-weighted meta estimates). When `I9_IHD` was NA for a variant (i.e., the variant failed the meta-coverage threshold for that endpoint despite being in the global SNP list), we substituted `I9_MI_STRICT` as the canonical CHD-related binary endpoint.

| Gene | Lead variant | chr-pos-ref-alt | Endpoint used | b_y | SE_y | meta_N |
|---|---|---|---|---|---|---|
| PCSK9 | rs11591147 | 1:55039974:G:T | I9_IHD | -0.207 | 0.012 | 5 |
| HMGCR | rs12916 | 5:75360714:T:C | I9_IHD | +0.0141 | 0.00321 | 5 |
| ANGPTL3 | rs10889353 | 1:62652525:A:C | I9_MI_STRICT (I9_IHD NA at this variant) | -0.0123 | 0.00532 | 5 |
| APOC3 | rs964184 (substitute) | 11:116778201:G:C | I9_IHD | -0.057 | 0.00442 | 5 |
| LPA | rs10455872 | 6:160589086:A:G | I9_MI_STRICT (I9_IHD NA at this variant) | +0.290 | 0.00981 | 5 |

### Substitution: APOC3 rs138326449 → rs964184

The canonical APOC3 R19X stop-gain `rs138326449` is in the PheWeb variant list, but `I9_IHD`, `I9_MI_STRICT`, `I9_CORATHER`, and `I9_ANGINA` are all NA at that variant (rs138326449 fg_af ≈ 3e-4, too rare for the meta to call cross-cohort). We substituted `rs964184` (chr11:116778201:G:C, fg_af = 0.855), the canonical APOA5–APOC3 cluster GWAS tag SNP for triglycerides and CHD (Klarin 2018, Do 2013). Sign of the alt-allele effect on TG is negative — same direction as the LoF.

### Substitution: ANGPTL3 and LPA — I9_IHD NA → I9_MI_STRICT

`rs10889353` and `rs10455872` are in the meta SNP list but `I9_IHD` is NA for both (no variant in chr1:62.5-62.8 MB region appears in the I9_IHD manhattan, suggesting an endpoint-level filter excludes these regions). `I9_MI_STRICT` is available and is the appropriate CHD-side endpoint (more specific than wide IHD).

## Exposure (b_x) — Open Targets Platform GraphQL credible-set summary stats

For each lead variant we queried Open Targets Platform GraphQL (`https://api.platform.opentargets.org/api/v4/graphql`) for the credible-set rows on lipid traits.

| Gene | OT Study ID | Trait | b_x | SE_x | p |
|---|---|---|---|---|---|
| PCSK9 | GCST90239659 | LDL cholesterol (Graham GLGC 2021) | -0.485 | 0.0369 | 1.7e-39 |
| HMGCR | GCST90239658 | LDL cholesterol (Graham GLGC 2021) | +0.0701 | 0.00141 | (very small) |
| ANGPTL3 | GCST90662857 | Triglyceride levels | -0.0817 | 0.00441 | 1e-76 |
| APOC3 | GCST004238 | Triglyceride levels (Klarin 2018) | -0.214 | 0.0126 | 2e-64 |
| LPA | GCST90090990 | LDL cholesterol (Lp(a) proxy; no Lp(a) GWAS in OT credible sets for this variant) | +0.208 | 0.0248 | 5e-17 |

### Allele harmonization

Open Targets `variantId` uses `chr_pos_ref_alt` with the same alt-allele convention as the FinnGen meta-PheWeb. Both b_x and b_y are reported per-ALT-allele, so the Wald ratio `b_y / b_x` is harmonized without sign flipping. Spot-check: PCSK9 rs11591147 alt=T (LoF), b_x on LDL = -0.485 (T lowers LDL — correct, since T is the loss-of-function allele); b_y on IHD = -0.207 (T lowers IHD risk — correct, by the protective evolocumab mechanism). Both negative, ratio positive.

## Sources that worked

- FinnGen × MVP × UKBB PheWeb `https://mvp-ukbb.finngen.fi/api/variant/...` — every variant lookup succeeded
- FinnGen autocomplete `https://mvp-ukbb.finngen.fi/api/autocomplete?...` — every rsID lookup succeeded
- Open Targets GraphQL `https://api.platform.opentargets.org/api/v4/graphql` — every credible-set lookup succeeded

## Sources that did not work / were skipped

- `https://public-mvp-ukbb.finngen.fi` — the about-page WebFetch returned empty content; direct `curl` to `mvp-ukbb.finngen.fi` was the working route.
- `https://gwas-api.mrcieu.ac.uk` — deprecated; the new `https://api.opengwas.io` requires a JWT.
- GLGC bulk download (`https://csg.sph.umich.edu/willer/public/glgc-lipids2021/`) — static gzipped files only, no single-variant API; we used Open Targets' per-variant credible-set rows instead (which incorporate GLGC).
- No direct Lp(a) GWAS hit was available for rs10455872 in the Open Targets credible-set rows. We substituted LDL-C as the exposure proxy (Lp(a) particles carry cholesterol; the rs10455872 LDL effect is dominated by Lp(a)). The direction is preserved.

## Raw JSON dumps saved

- `per_gene/PCSK9_pheweb_raw.json`, `HMGCR_pheweb_raw.json`, `ANGPTL3_pheweb_raw.json`, `LPA_pheweb_raw.json` — FinnGen variant pages
- `per_gene/APOC3_pheweb_raw.json` (original rs138326449, all NA for IHD)
- `per_gene/APOC3_pheweb_raw_rs964184.json`, `_rs5128.json`, `_rs2854117.json` — substitute candidates
- `per_gene/PCSK9_opentargets_creds.json`, `HMGCR_*`, `ANGPTL3_*`, `APOC3_*`, `LPA_*` — Open Targets credible-set rows
- `panel_input.json` — harmonized 5-gene panel input to mrAR
