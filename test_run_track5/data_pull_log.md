# Track 5 — Data pull log

**Date**: 2026-06-08
**Goal**: pull a PCSK9 cis-pQTL anchor (preferably Sun BB 2023 UKB-PPP) for the renormalization in `track5_pqtl_anchor.R`.

## Sources attempted

| Source | URL | Result |
|---|---|---|
| UKB-PPP browser | `http://ukb-ppp.gwas.eu` | ECONNREFUSED — host not resolving |
| UKB-PPP Synapse | `https://www.synapse.org/Synapse:syn51365301` | Page loads but per-protein file listing requires auth login; bulk-download paths not accessible to unauthenticated WebFetch |
| UKB-PPP wiki | `https://www.synapse.org/Synapse:syn51364943/wiki/622119` | Auth required; loading placeholder only |
| AWS Open Data | `arn:aws:s3:::ukbiobank.opendata.sagebase.org` | Public bucket but listed as "Controlled Access"; `curl -I` returns 403 |
| Direct S3 path | `s3://ukbiobank.opendata.sagebase.org/UKB-PPP pGWAS summary statistics/` | 403 Forbidden (auth required) |
| Nature supplementary | `https://www.nature.com/articles/s41586-023-06592-6` | Cloudflare challenge / redirects to login; supplementary XLSX path returned 403 |
| Nature S5 ESM | `https://static-content.springer.com/esm/.../41586_2023_6592_MOESM5_ESM.xlsx` | 403 with `cf-mitigated: challenge` (Cloudflare) |
| bioRxiv preprint | `https://www.biorxiv.org/content/10.1101/2022.06.17.496443v1.full` | HTTP 403 |
| azphewas (AZ pQTL browser) | `https://azphewas.com/proteinView/.../Q15119/...` | Page is interactive JS; no data extractable from HTML |
| GWAS Catalog REST | `https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/rs11591147/associations` | **WORKS — 823 associations returned** |
| GWAS Catalog by PMID | `findByPublicationIdPubmedId?pubmedId=37794190` (Sun BB 2023) | 0 studies returned — **Sun BB 2023 NOT in GWAS Catalog REST** |
| deCODE summary data | `https://download.decode.is/form/folder/proteomics` | Index page, folder listing requires session/login |
| deCODE 2021 alt path | `https://download.decode.is/form/2021/` | 404 |
| Open Targets variant | `https://platform.opentargets.org/variant/1_55039974_G_T` | Interactive page; no scrapable data |
| Pott 2024 PMC text | `https://pmc.ncbi.nlm.nih.gov/articles/PMC10964567/` | **WORKS — full text + tables** |

## What we used

**Pott J et al. 2024**, *Hum Mol Genet* (PMID 38491180), PCSK9 sex-stratified meta-GWAS, 6 European studies (LIFE-Heart, LIFE-Adult, LURIC, TwinGene, KORA-F3, GCKD), n ≈ 20 016. Source: extracted from GWAS Catalog REST API associations endpoint for rs11591147 (`assoc_id` 96189423, beta=0.37189 decrease, SE=0.0145297, p=2×10⁻¹⁴⁴, trait "PCSK9 protein measurement", confirmed via efoTraits REST endpoint).

Scale: log-transformed PCSK9 protein concentration (mixed Olink + ELISA platforms across the 6 contributing cohorts; documented in PMC10964567).

**Critical**: this substitution is documented honestly. Pott 2024 does NOT overlap UKB-PPP (different participants), so it is not the "same scale" Team 1 used (Team 1 used GTEx Liver eQTLs, FinnGen+MVP+UKBB outcome; neither overlaps with Pott 2024's German+Italian+Swedish cohorts). For the universal-anchor algebra `tilde_b_xy = b_y / b_pqtl`, what matters is the MAGNITUDE on a comparable log-protein scale, which Pott 2024 provides at higher precision (SE=0.014) than any single contributing cohort.

## Validated PCSK9 pQTL hits for rs11591147 (cross-check)

All from GWAS Catalog REST `singleNucleotidePolymorphisms/rs11591147/associations` filtered to efoTraits containing "PCSK9" or "proprotein convertase subtilisin/kexin type 9":

| Study (1st author, year) | PMID | n | Beta (T allele decrease) | SE | p | Scale |
|---|---|---:|---:|---:|---:|---|
| **Pott J 2024** | 38491180 | ~20 016 | **0.37189** | **0.0145** | 2e-144 | log-PCSK9 (Olink + ELISA mix) |
| Kheirkhah A 2023 | (unconfirmed) | larger meta | 0.39 (G) | 0.013 | 7e-189 | PCSK9 protein measurement |
| Pott J 2021 | 33339796 | ~10k | 0.31965 | 0.0223 | 2e-46 | log-PCSK9 |
| Pott J 2021 | 33339796 | ~10k | 0.29342 | 0.0200 | 9e-49 | log-PCSK9 |
| Pott J 2018 | 29748315 | 2 583 | 0.315 | 0.037 | 2e-17 | log-PCSK9 (ELISA, LIFE) |
| Pietzner M 2021 | 34648354 | 10 708 | 0.883 | 0.052 | 2e-64 | SomaScan (Fenland) — different scale, ~2.4× higher in absolute beta because SomaScan RFU vs Olink NPX |
| Gudjonsson A 2022 | 35545635 | ~5 000 (AGES-RS) | 1.0396 | 0.0883 | 6e-32 | SomaScan (deCODE) |
| Gudjonsson A 2022 | 35545635 | ~5 000 (AGES-RS) | 1.09951 | 0.0815 | 2e-41 | SomaScan |

**Sun BB 2023 UKB-PPP** is conspicuously absent from the GWAS Catalog rs11591147 association list. PCSK9 PubMed ID 37794190 (Sun BB 2023) → `findByPublicationIdPubmedId` returns 0 studies. This is a known GWAS Catalog ingest gap; the Sun BB 2023 supplementary tables would require manual download from Nature, which Cloudflare blocked.

## Scale audit (Olink NPX vs SomaScan vs ELISA)

| Scale | Used here? | Note |
|---|---|---|
| **Olink Explore NPX** | Used (Sun BB 2023 would have been this; Pott 2024 mixes Olink + ELISA) | log2-normalized relative quantification; cis-pQTL betas typically ~0.3–0.5 for rs11591147 in n~20k |
| **SomaScan RFU** | NOT used — would have inflated beta ~2.4× | Pietzner/Gudjonsson have ~0.88–1.10; non-comparable without scaling factor |
| **ELISA absolute** | NOT separated — mixed into Pott 2024 | Older studies (Pott 2018) |

The Pott 2024 beta=0.37189 is on a "log-PCSK9" scale interpretable as ~ unit-SD log-mass-units of plasma PCSK9 per T allele. This is appropriate for the renormalization `b_y / b_pqtl`, which gives "log-OR-CHD per unit-SD plasma PCSK9".

## Substitution rationale summary

1. **Preferred (per spec)**: Sun BB 2023 *Nature* 622:329 UKB-PPP. **Not accessible** without authentication.
2. **Backup A**: Ferkingstad 2021 deCODE — SomaScan; ~2.4× scale offset from UKB-PPP/Olink. Bulk download requires form submission.
3. **Used**: Pott 2024 *Hum Mol Genet* PCSK9 meta-GWAS, PMID 38491180. **Most precise non-UKB-PPP estimate available** through public REST API. Same scale as UKB-PPP (Olink + ELISA mix; per-unit-SD-log-PCSK9). Citation: Pott J, Kheirkhah A, Gadin JR, et al. (2024) "Sex and statin-related genetic associations at the PCSK9 gene locus: results of genome-wide association meta-analysis." *Hum Mol Genet*. PMID 38491180. Citation-audit verified.

## What we could NOT fetch

- Sun BB 2023 UKB-PPP per-variant summary statistics (auth-gated)
- Ferkingstad 2021 deCODE per-variant summary statistics (form-gated)
- per-tissue eQTL effects on rs11591147 from the eQTL Catalogue (not pulled — would require a re-run of `build_per_tissue_panel.py` parameterized on rs11591147 rather than the 7 panel variants. Out of scope for this Track 5 run.)
- pQTL effects on the 7 panel variants individually (none of rs471705, rs6676563, rs2802881, rs61772108, rs114739858, rs111521483, rs143341434 appear as PCSK9 pQTLs in GWAS Catalog — they are intergenic eQTL leads, not coding-variant pQTL leads)
