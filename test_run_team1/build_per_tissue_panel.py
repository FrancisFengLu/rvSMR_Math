#!/usr/bin/env python3
"""Build per-tissue PCSK9 panels for cell_type_q() analog.

Substrate sub-substitution: instead of TenK10K Phase 1 28 PBMC cell types
(rare-variant Zenodo placeholders), use GTEx v8 tissues from eQTL Catalogue
as "cell types" for the cell_type_q() test. This is NOT a single-cell
substrate -- it is bulk-tissue resolution -- but it exercises the same
plumbing the cell-type-q() function is designed for.

Tissues chosen (lipid-pathway relevant):
   liver (QTD000266, n=208)     -- PCSK9 native tissue
   blood (QTD000356, n=670)     -- whole blood
   adipose (QTD000116, n=581)
   adipose-visceral (QTD000121, n=469)
   artery-aorta (QTD000131, n=387)
   artery-coronary (QTD000136, n=213)
   artery-tibial (QTD000141, n=584)
   muscle (QTD000281, n=702)
   small intestine (QTD000321, n=174)

For each tissue, look up the lead variant rs471705 (PCSK9 strongest liver
eQTL) and report the per-tissue Wald ratio. Then run cell_type_q() across
the tissues that have non-NA cis-eQTL.
"""
from __future__ import annotations
import json
from pathlib import Path
import pysam


OUT_DIR = Path(__file__).parent
PCSK9_ENSG = "ENSG00000169174"
PCSK9_CHR = "1"
PCSK9_TSS_START = 55039445
PCSK9_TSS_END = 55064852
CIS_WIN = 1_000_000

TISSUES = [
    ("QTD000266", "liver",            208),
    ("QTD000356", "blood",            670),
    ("QTD000116", "adipose",          581),
    ("QTD000121", "adipose_visceral", 469),
    ("QTD000131", "artery_aorta",     387),
    ("QTD000136", "artery_coronary",  213),
    ("QTD000141", "artery_tibial",    584),
    ("QTD000281", "muscle",           702),
    ("QTD000321", "small_intestine",  174),
]

BASE_URL = ("https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/"
            "QTS000015/{ds}/{ds}.all.tsv.gz")

COLS = [
    "molecular_trait_id", "chromosome", "position", "ref", "alt",
    "variant", "ma_samples", "maf", "pvalue", "beta", "se",
    "type", "ac", "an", "r2",
    "molecular_trait_object_id", "gene_id", "median_tpm", "rsid",
]


def fetch_pcsk9(ds: str) -> list[dict]:
    url = BASE_URL.format(ds=ds)
    print(f"  tabix {ds}", flush=True)
    tbx = pysam.TabixFile(url)
    rows = []
    for line in tbx.fetch(PCSK9_CHR, PCSK9_TSS_START - CIS_WIN,
                          PCSK9_TSS_END + CIS_WIN):
        f = line.split("\t")
        if f[0] != PCSK9_ENSG and (len(f) <= 16 or f[16] != PCSK9_ENSG):
            continue
        rec = dict(zip(COLS, f))
        try:
            rec["position"] = int(rec["position"])
        except (ValueError, TypeError):
            continue
        for k in ("maf", "pvalue", "beta", "se"):
            try:
                rec[k] = float(rec[k])
            except (ValueError, TypeError):
                rec[k] = None
        rows.append(rec)
    tbx.close()
    return rows


def main():
    # Load the lead instruments from Track 2 panel
    with open(OUT_DIR / "pcsk9_track2_panel.json") as f:
        panel = json.load(f)
    leads = panel["instruments"]

    # For each lead variant, look up its cis-eQTL across all tissues
    per_tissue = {}
    for ds, label, n in TISSUES:
        rows = fetch_pcsk9(ds)
        print(f"  {label}: {len(rows)} PCSK9 cis-eQTL records")
        key2row = {(str(r["chromosome"]), r["position"], r["ref"], r["alt"]): r
                   for r in rows}
        # Find each lead in this tissue's records
        tis = {"dataset_id": ds, "label": label, "n_donors": n,
               "lead_lookups": []}
        for L in leads:
            key = (str(L["chr"]), L["pos"], L["ref"], L["alt"])
            r = key2row.get(key)
            if r is None:
                # try ref/alt swap
                key_sw = (str(L["chr"]), L["pos"], L["alt"], L["ref"])
                r_sw = key2row.get(key_sw)
                if r_sw is not None:
                    r = dict(r_sw)
                    if r.get("beta") is not None:
                        r["beta"] = -r["beta"]
                    r["_flipped"] = True
                else:
                    tis["lead_lookups"].append({
                        "rsid": L["rsid"], "pos": L["pos"],
                        "found": False, "beta": None, "se": None,
                    })
                    continue
            tis["lead_lookups"].append({
                "rsid": L["rsid"], "pos": L["pos"],
                "found": True, "beta": r.get("beta"), "se": r.get("se"),
                "pvalue": r.get("pvalue"), "maf": r.get("maf"),
                "flipped": r.get("_flipped", False),
            })
        per_tissue[label] = tis

    with open(OUT_DIR / "pcsk9_per_tissue.json", "w") as f:
        json.dump(per_tissue, f, indent=2)
    print(f"\nWrote pcsk9_per_tissue.json ({len(per_tissue)} tissues)")


if __name__ == "__main__":
    main()
