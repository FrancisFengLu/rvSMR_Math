# Track 5 — Worker report

**Date**: 2026-06-08
**Driver**: `track5_pqtl_anchor.R`
**Output JSON**: `per_tissue_renormalized.json`

## Headline result

| Stage | Q | df | p-value | Interpretation |
|---|---:|---:|---:|---|
| Team 1 unanchored (reproduce) | 8.275 | 3 | **0.0407** | discordant_investigate |
| Track 5 pQTL-anchored | 0.000 | 3 | **1.000** | concordant (by construction) |

p collapsed from 0.041 → 1.000.

## What I did

1. **Pulled** the Team 1 panel and the per-tissue lookup file from `test_run_team1/`.
2. **Found** the PCSK9 cis-pQTL anchor for rs11591147 (canonical PCSK9 LOF coding variant) from the GWAS Catalog REST API. Sun BB 2023 UKB-PPP is NOT in the GWAS Catalog (`findByPublicationIdPubmedId?pubmedId=37794190` returns 0 studies); the Synapse / AWS / Nature paths require authentication that WebFetch could not provide. **Substitution**: **Pott J 2024** *Hum Mol Genet* (PMID 38491180), PCSK9 sex-stratified meta-GWAS, n ≈ 20 016 European, **β = 0.37189, SE = 0.0145, p = 2e-144** on log-PCSK9 (Olink + ELISA mix). Full source-substitution audit in `data_pull_log.md`.
3. **Verified** that none of the 7 Team 1 panel eQTL lead variants (rs471705, rs6676563, rs2802881, rs61772108, rs114739858, rs111521483, rs143341434) are in any PCSK9 pQTL study indexed by GWAS Catalog — they are intergenic eQTL leads ~150 kb–1 Mb from rs11591147. Expected: per-variant matched anchors are not available; **forced to use the spec's "universal anchor" fallback**.
4. **Ran** the renormalization `tilde_b_xy = b_y / b_pqtl` per-tissue with delta-method variance, and ran `cell_type_q()` (from `rvMR`, freshly re-installed from `/home/francisfenglu4/rvSMR/May_30md/rvMR/` source).
5. **Cross-checked** by computing the Cochran-Q manually (Internal Review §b) — matches the rvMR output.

## What I found

The unanchored Q-test reproduces Team 1's reported 0.041 exactly (modulo rounding: 8.275 vs 8.28). The pQTL-anchored Q collapses to **exactly 0** because the universal anchor `b_y(rs471705) / b_pqtl(rs11591147)` is identical across all 4 tissues — the tissue index t doesn't survive the algebra when the same variant is used everywhere.

This is the **algebraically expected** outcome of the universal-anchor fallback the spec describes, **NOT** an independent empirical confirmation of biological homogeneity. The hypothesis ("cross-tissue heterogeneity is a secreted-to-plasma scale artifact") is, in this setup, **confirmed by construction**: the scale variation IS the per-tissue variation in `b_x_tissue`, and the anchor algebraically removes it.

## Interpretation under spec's framework

- **"p increases substantially" (0.041 → 0.3)**: **passed** — actually p → 1.000, the maximal increase. Spec says: "Hypothesis confirmed — heterogeneity was scale variation, not biology. This is a clean positive finding for the rvSMR Step 11/12 design."
- **Caveat**: the collapse is degenerate because the same variant is used across all tissues. The spec's fallback algebra (universal anchor variant) reaches its theoretical maximum effect here. A meaningful test would require per-tissue eQTL leads + per-variant pQTL anchors — not feasible because the panel variants aren't in pQTL databases. See `pcsk9_pqtl_anchor_results.md` §"Honest interpretation".

So: **hypothesis confirmed in the sense the spec defined**, but the confirmation is closer to "the math works out exactly as the universal-anchor algebra predicts" than to "the data refuted real tissue-specific biology". Both readings are coherent.

## Files in `test_run_track5/`

| File | Purpose |
|---|---|
| `track5_pqtl_anchor.R` | Analysis driver |
| `pcsk9_pqtl_anchor_results.md` | Main results document |
| `data_pull_log.md` | Audit of what URLs / sources were attempted and used |
| `pcsk9_pqtl_lookup_raw.json` | Raw per-variant GWAS Catalog PCSK9 pQTL search results (the 7 panel variants — all returned 0 PCSK9 pQTL hits) |
| `per_tissue_renormalized.json` | Per-tissue b_pqtl/SE/tilde_b_xy/SE; primary + secondary analyses |
| `per_tissue_best_variant.json` | Per-tissue F-stat-maximizing panel variant (used in secondary analysis) |
| `TRACK5_INTERNAL_REVIEW.md` | Internal-reviewer sub-pass: algebra check, independent Q-recompute, citation audit, Olink/SomaScan scale audit |
| `WORKER_TRACK5_REPORT.md` | This file |
| `DONE` | sentinel |

## Hard-constraint compliance

- **Did not modify rvMR R package**: confirmed. Reinstalled the package from `/home/francisfenglu4/rvSMR/May_30md/rvMR/` source via `install.packages(..., type="source")` to make `cell_type_q()` (which exists in source but wasn't in the conda-env built copy) available. **No source code modified.**
- **Did not touch `test_run_team1/`, `test_run_v3*/`, `test_run_finngen/`**: confirmed; all reads, no writes.
- **Sun BB 2023 cited correctly** as Sun BB et al. 2023 *Nature* 622(7982):**329–338**. NOT confused with 622:339 (Dhindsa rare-variant) or Sun KY 2024 *Nature* 631:583 (RGC-ME).

## Caveats and what I couldn't do

1. **Could not fetch Sun BB 2023 UKB-PPP per-variant summary stats.** All public access routes (Synapse, AWS S3, Nature supplementary) were auth-gated. Used Pott 2024 as the closest substitute on the Olink-compatible log-PCSK9 scale; documented honestly.
2. **Could not fetch per-variant pQTL anchors for the 7 panel variants.** None are in GWAS Catalog as PCSK9 pQTLs (expected — they are intergenic eQTL leads). The universal-anchor fallback is what the spec calls for in this case; we applied it.
3. **Did not query the eQTL Catalogue for rs11591147 in each tissue.** That would be a meaningful follow-up analysis (anchor on the CANONICAL pQTL variant in each tissue's eQTL data, then compute per-tissue Wald `b_y(rs11591147) / b_pqtl(rs11591147)`). Out of scope for the universal-anchor Track 5 spec; flagged in `pcsk9_pqtl_anchor_results.md`.

## Recommended next steps (if extending)

- Pull `b_x_tissue(rs11591147)` from the eQTL Catalogue across the 9 GTEx tissues (use Team 1's `build_per_tissue_panel.py` pattern, parameterized on rs11591147 instead of rs471705). This gives a TISSUE-VARYING denominator if you anchor as `b_y(rs11591147) / b_x_tissue(rs11591147)`, and a NON-DEGENERATE Q test for tissue heterogeneity AT the canonical pQTL lead variant.
- Pull Sun BB 2023 UKB-PPP supplementary tables via authenticated Synapse access (requires UKB approval) for the gold-standard anchor.
- Run rs562556 (PCSK9 I474V) as a second anchor — both are coding variants with well-characterized pQTL effects on the Olink scale.

---

**Bottom line**: Track 5 ran end-to-end on real data. The unanchored result reproduced Team 1's 0.041. The pQTL-anchored result collapsed Q to 0 (p = 1). The collapse is algebraically expected for the universal-anchor fallback; the hypothesis is confirmed in the sense the spec defines but the test as run is non-discriminating (it would have collapsed Q regardless of the underlying biology). Reported honestly.
