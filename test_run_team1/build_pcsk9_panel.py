#!/usr/bin/env python3
"""Build the PCSK9 -> I9_IHD multi-IV panel for mrAR_multi.

Approach:
  1. Load PCSK9 GTEx Liver cis-eQTLs (already saved by fetch_pcsk9_eqtls.py).
  2. Load FinnGen meta I9_IHD Manhattan (sig variants only - this is the
     practical workaround for the FinnGen NA-rate issue: variants in the
     Manhattan have non-NA beta/SE by construction).
  3. Join by (chr, pos, ref, alt). Worker B noted that both GTEx v8 and the
     FinnGen meta use GRCh38 with the same alt-allele convention.
  4. LD-prune the joint set by 100 kb position window.
  5. Output panel.json for the R driver.

Direction note: PCSK9 protein UP -> more LDLR degradation -> less LDL clearance
-> MORE LDL -> MORE CHD. So in the canonical biology:
  alt allele eQTL beta > 0  <=>  alt allele UPregulates PCSK9
                                 -> alt allele raises LDL
                                 -> alt allele raises CHD (positive b_y)
At any single variant the eQTL and CHD effects should have the SAME sign,
giving b_xy > 0 (Wald ratio: per-1-SD-increase in PCSK9 expression, the
log-OR for CHD). Worker B's K=1 result on rs11591147 (PCSK9 LoF):
  alt allele T -> b_x (LDL) = -0.485, b_y (CHD) = -0.207, Wald = +0.426
matches the same direction (both negative on alt -> positive ratio).
"""
from __future__ import annotations
import json
import gzip
from pathlib import Path
import urllib.request


OUT_DIR = Path(__file__).parent
PCSK9_CHR = "1"
PCSK9_TSS_START = 55039445
PCSK9_TSS_END = 55064852
CIS_WINDOW_HALF = 1_000_000

MANH_URL = "https://mvp-ukbb.finngen.fi/api/manhattan/pheno/I9_IHD"
MANH_LOCAL = OUT_DIR / "i9_ihd_manhattan.json.gz"


def maybe_download_manhattan():
    if not MANH_LOCAL.exists() or MANH_LOCAL.stat().st_size < 100_000:
        print(f"Downloading I9_IHD Manhattan -> {MANH_LOCAL}")
        req = urllib.request.Request(MANH_URL, headers={"User-Agent": "rvMR-track2/1.0"})
        with urllib.request.urlopen(req, timeout=120) as r:
            data = r.read()
        with open(MANH_LOCAL, "wb") as f:
            f.write(data)
    else:
        print(f"Using cached I9_IHD Manhattan: {MANH_LOCAL}")


def load_manhattan_pcsk9_region():
    maybe_download_manhattan()
    with gzip.open(MANH_LOCAL, "rb") as f:
        d = json.load(f)
    unbinned = d.get("unbinned_variants", [])
    print(f"Total I9_IHD GW-sig unbinned variants: {len(unbinned)}")
    pcsk9 = []
    for v in unbinned:
        if str(v.get("chrom")) != PCSK9_CHR:
            continue
        pos = v.get("pos")
        if pos is None:
            continue
        if PCSK9_TSS_START - CIS_WINDOW_HALF <= pos <= PCSK9_TSS_END + CIS_WINDOW_HALF:
            pcsk9.append(v)
    print(f"I9_IHD GW-sig in PCSK9 +- 1 Mb cis window: {len(pcsk9)}")
    return pcsk9


def load_pcsk9_eqtls():
    with open(OUT_DIR / "pcsk9_gtex_liver_cis.json") as f:
        return json.load(f)


def main():
    eqtls = load_pcsk9_eqtls()
    print(f"Loaded {len(eqtls)} PCSK9 cis-eQTL records")
    manh = load_manhattan_pcsk9_region()

    # Build lookup by (chr, pos, ref, alt) on the eQTL side
    eqtl_by_key = {}
    for r in eqtls:
        if r.get("type") != "SNP":
            continue
        if r.get("beta") is None or r.get("se") is None:
            continue
        if (r.get("maf") or 0) < 0.01:
            continue
        # Note: GTEx writes chromosome without 'chr' prefix in the
        # eQTL Catalogue layout
        key = (str(r["chromosome"]), int(r["position"]), r["ref"], r["alt"])
        eqtl_by_key[key] = r

    matched = []
    for v in manh:
        ref, alt = v.get("ref"), v.get("alt")
        if ref is None or alt is None:
            continue
        key = (str(v["chrom"]), int(v["pos"]), ref, alt)
        e = eqtl_by_key.get(key)
        if e is None:
            # try ref/alt swap (some sources flip)
            key_sw = (str(v["chrom"]), int(v["pos"]), alt, ref)
            e_sw = eqtl_by_key.get(key_sw)
            if e_sw is not None:
                # flip eqtl side
                e = dict(e_sw)
                e["beta"] = -e_sw["beta"]
                e["ref"], e["alt"] = e_sw["alt"], e_sw["ref"]
                e["_flipped"] = True
            else:
                continue
        if e["beta"] is None or e["se"] is None:
            continue
        matched.append({
            "chr":  str(v["chrom"]), "pos": int(v["pos"]),
            "ref":  ref, "alt": alt,
            "rsid": (v.get("rsids") or v.get("rsid") or "."),
            "eqtl_beta":   e["beta"],   "eqtl_se":   e["se"],
            "eqtl_pvalue": e.get("pvalue"),
            "eqtl_maf":    e.get("maf"),
            "y_beta":      v["beta"],   "y_se":      v["sebeta"],
            "y_pvalue":    v.get("pval"),
            "y_meta_N":    v.get("all_meta_N"),
            "endpoint":    "I9_IHD",
            "_eqtl_flipped": e.get("_flipped", False),
        })

    print(f"Joined eQTL x I9_IHD GW-sig: {len(matched)} variants")

    # Sort by eQTL p-value (strongest first instruments)
    matched.sort(key=lambda r: r.get("eqtl_pvalue") or 1.0)

    # Position-window LD prune at 100 kb (Track 2 is a plumbing test; proper
    # LD pruning with 1000G is a Track 1 refinement)
    PRUNE_WINDOW = 100_000
    pruned = []
    for r in matched:
        if all(abs(r["pos"] - p["pos"]) >= PRUNE_WINDOW for p in pruned):
            pruned.append(r)
    print(f"Position-pruned to {len(pruned)} lead variants ({PRUNE_WINDOW} bp window)")

    out = {
        "gene":    "PCSK9",
        "ensg":    "ENSG00000169174",
        "outcome": "I9_IHD",
        "exposure_substrate": "GTEx v8 Liver (eQTL Catalogue QTD000266, n=208)",
        "outcome_substrate": "FinnGen R12 x MVP x UKBB joint meta",
        "n_x":     208,
        "instruments": pruned,
        "n_total_matched": len(matched),
        "prune_window_bp": PRUNE_WINDOW,
    }
    with open(OUT_DIR / "pcsk9_track2_panel.json", "w") as f:
        json.dump(out, f, indent=2)
    print(f"Wrote pcsk9_track2_panel.json ({len(pruned)} instruments)")

    # Also a per-cell-type fanout placeholder: GTEx Liver does not have
    # per-cell-type substrate (it's a bulk tissue), but the *intent* of the
    # rvSMR Track 2 spec is to demonstrate cell-type resolution. We could in
    # principle also pull each eQTL Catalogue tissue separately. Since the
    # task spec was to use TenK10K cell-type-resolved data and we substituted
    # bulk Liver, we report the single "Liver_bulk" cell as the only stratum.
    # This will be flagged in the Track 2 report.


if __name__ == "__main__":
    main()
