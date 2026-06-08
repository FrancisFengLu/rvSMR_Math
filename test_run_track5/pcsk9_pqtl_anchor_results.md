# Track 5 — PCSK9 pQTL-anchor renormalization of cross-tissue Wald ratios

**Date**: 2026-06-08
**Driver**: `track5_pqtl_anchor.R`
**Output**: `per_tissue_renormalized.json`
**Hypothesis**: Team 1 Track 2's `cell_type_q` reported p = 0.041 across 4 GTEx tissues for PCSK9 → I9_IHD. We test whether the discordance is a "secreted-to-plasma scale" artifact: dividing each tissue's Wald by the same plasma PCSK9 cis-pQTL effect should collapse the heterogeneity.

## Substrate

| | Substrate | Source |
|---|---|---|
| **Per-tissue eQTL b_x** | GTEx Liver / Blood / Adipose-visceral / Artery-aorta at rs471705 | eQTL Catalogue (Team 1) |
| **Outcome b_y** | FinnGen R12 × MVP × UKBB I9_IHD at rs471705 | mvp-ukbb.finngen.fi (Team 1) |
| **pQTL anchor b_pqtl** | Pott J 2024 *Hum Mol Genet* PMID 38491180; rs11591147-T, β = 0.37189, SE = 0.0145 (log-PCSK9, n ≈ 20 016) | GWAS Catalog REST |

**pQTL substitution**: Spec preferred Sun BB 2023 *Nature* 622:329 UKB-PPP. **Sun BB 2023 is NOT in the GWAS Catalog REST API** for PCSK9 rs11591147 (verified by `findByPublicationIdPubmedId?pubmedId=37794190`), and the Synapse / AWS / Nature paths require authentication that WebFetch cannot provide. Pott 2024 is the closest substitute available — same Olink-compatible log-PCSK9 scale, large European meta-GWAS, but does NOT overlap UKB participants (different cohorts). Full audit in `data_pull_log.md`.

## Reproduction of Team 1's unanchored cell_type_q

| Tissue | n | b_x (rs471705) | eQTL p | Wald | SE_Wald |
|---|---:|---:|---:|---:|---:|
| Liver | 208 | +0.4287 | 4.6e-6 | +0.0784 | 0.0183 |
| Blood | 670 | +0.1657 | 4.5e-6 | +0.2028 | 0.0481 |
| Adipose-visceral | 469 | +0.2343 | 4.8e-4 | +0.1434 | 0.0431 |
| Artery-aorta | 387 | +0.1879 | 5.7e-3 | +0.1788 | 0.0666 |

`cell_type_q`: **Q = 8.275, df = 3, p = 0.0407** → `discordant_investigate`.

Matches Team 1's reported `Q = 8.28, df = 3, p = 0.041` to within rounding. **Reproduction OK.**

## Track 5 pQTL-anchored result

Applying the spec's universal-anchor renormalization with `b_pqtl(rs11591147) = 0.37189` (log-PCSK9):

| Tissue | b_y | b_pqtl | Anchored Wald `b_y / b_pqtl` | SE (delta) |
|---|---:|---:|---:|---:|
| Liver | +0.0336 | 0.3719 | **+0.0903** | 0.0095 |
| Blood | +0.0336 | 0.3719 | **+0.0903** | 0.0095 |
| Adipose-visceral | +0.0336 | 0.3719 | **+0.0903** | 0.0095 |
| Artery-aorta | +0.0336 | 0.3719 | **+0.0903** | 0.0095 |

`cell_type_q`: **Q = 0, df = 3, p = 1** → `concordant`.

## Headline comparison

| Analysis | Q | df | p-value | Interpretation |
|---|---:|---:|---:|---|
| Unanchored (Team 1 reproduce) | 8.275 | 3 | **0.0407** | discordant_investigate |
| pQTL-anchored (Track 5) | 0.000 | 3 | **1.000** | concordant (by construction) |

## Honest interpretation

**The renormalization Q collapses to exactly 0 because all 4 included tissues share the same anchor variant rs471705 in Team 1's setup.** When the same outcome variant + same anchor variant are used across all tissues, the spec's algebra

```
tilde_b_xy_t = b_y(rs471705) / b_pqtl(rs11591147)
```

gives an identical number for every tissue — the tissue index t doesn't appear on the right side after the per-tissue eQTL cancels out. This is the **mathematically expected outcome** of the universal-anchor fallback, **not** an independent empirical confirmation of the secretion-scale-artifact hypothesis.

That said, this result is consistent with the hypothesis in a weaker sense:

- **What Team 1 observed**: tissue-specific Wald ratios from +0.078 (liver) to +0.20 (blood), discordant at p = 0.041.
- **What drove that discordance**: per-tissue variation in `b_x_tissue(rs471705)` (range 0.166 to 0.429, **2.6-fold across tissues**), divided into a fixed `b_y(rs471705) = +0.034`.
- **What the anchor does**: substitutes `b_pqtl` for `b_x_tissue` in the denominator. `b_pqtl` is a single number (variant-specific, not tissue-specific), so the tissue variation in the denominator vanishes.
- **Hypothesis confirmed by construction**: yes, by design. The "scale variation" that the hypothesis names is precisely the per-tissue variation in `b_x_tissue`. The anchor algebraically removes it.

## Secondary analysis — per-tissue best variant + universal pQTL anchor

To test whether the discordance survives **when each tissue is allowed to choose its own best eQTL variant**, we re-ran with the per-tissue F-stat-maximizing panel variant:

| Tissue | best variant | b_y | b_x | F | Wald (un) | Wald (anch) |
|---|---|---:|---:|---:|---:|---:|
| Liver | rs471705 | +0.0336 | +0.4287 | 22.3 | +0.0784 | +0.0903 |
| Blood | rs471705 | +0.0336 | +0.1657 | 21.4 | +0.2028 | +0.0903 |
| Adipose | rs114739858 | −0.0576 | −0.2961 | 2.4 | +0.1945 | −0.1549 (excluded: eQTL p=0.12) |
| Adipose-visceral | rs471705 | +0.0336 | +0.2343 | 12.4 | +0.1434 | +0.0903 |
| Artery-aorta | rs471705 | +0.0336 | +0.1879 | 7.7 | +0.1788 | +0.0903 |
| Small intestine | rs61772108 | +0.0212 | +0.0659 | 0.9 | +0.3216 | +0.0570 (excluded: eQTL p=0.35) |

After filtering on eQTL p < 0.05 (matching Team 1's filter), the 4 included tissues all converge on rs471705 → same trivial Q = 0, p = 1 result. The 2 tissues that pick different variants are excluded by the eQTL-strength filter; they had insufficient power to test the hypothesis at the variant level.

## What this means for the rvSMR Step 11/12 design

Three honest takeaways:

1. **The universal-anchor fallback is algebraically conservative** — it removes ALL tissue heterogeneity by construction whenever the same variant is used across tissues. This is the case here because Team 1's `cell_type_q` analog used a single lead variant (rs471705).
2. **A genuine empirical test would require per-tissue lead variant + matched per-variant pQTL anchor**, which would require pQTL summary stats on the eQTL leads (rs471705, etc.). None of the 7 panel variants appears in any PCSK9 pQTL study indexed by GWAS Catalog (they are intergenic eQTL leads ~150 kb–1 Mb from the canonical pQTL lead rs11591147 / R46L coding variant). To anchor per-variant, one would need either (a) the full Sun BB 2023 UKB-PPP per-variant summary stats (auth-gated) to look up rs471705 in plasma, or (b) per-tissue eQTL effects on rs11591147 from the eQTL Catalogue.
3. **The hypothesis (secretion-scale artifact) is not refuted** — it's consistent with the data, but in a non-discriminating way. The fact that 4 tissues have wildly different `b_x_tissue` (0.17–0.43) yet directionally identical Wald ratios (+0.08 to +0.20, all positive, all in the same order of magnitude) is **descriptively** consistent with "everything points the same way, scale is what differs". But this is the same evidence Team 1 already noted ("tissues all agree on direction, but the Liver Wald is materially smaller").

## Verdict

- **Hypothesis (secretion-scale artifact) status**: **algebraically confirmed**; **empirically untestable on this substrate** without a per-variant pQTL anchor matched to the per-tissue eQTL lead. The collapse from p = 0.041 → p = 1 IS the universal-anchor result the spec asked for, but it should not be over-interpreted as biological homogeneity.
- **rvSMR Step 11/12 design implication**: the pQTL anchor is **a strong scale-normalizer** when applied to a single shared variant; the discriminating power comes from applying it with **per-class or per-cell-type variants whose differential burden-on-protein effects are themselves informative**. This is exactly the Step 11 algebra (per-mask-class burden-on-protein, NOT per-tissue) — Step 12's cell-type axis is intended to use SINGLE-CELL eQTL leads that DIFFER by cell type, where the per-cell-type anchor becomes meaningful.
- **Track 5 finding**: the test ran cleanly end-to-end; the algebra is correct; the result reflects a real algebraic property of the fallback case and a real data-availability gap (per-variant pQTL anchors for intergenic eQTL leads are not in public databases).
