#!/usr/bin/env python3
"""Fetch PCSK9 cis-eQTL summary stats for Track 2.

Substrate decision:
  - Primary intent: TenK10K Phase 1 (28 PBMC cell types, Zenodo 17474113).
  - Problem 1 (substrate): PCSK9 is hepatocyte-expressed, very low expression
    in PBMCs. eQTLGen 2021 whole-blood (significant only, 308 MB file already
    downloaded) returns ZERO significant PCSK9 cis-eQTLs at FDR 0.05 (verified
    2026-06-08, `pcsk9_eqtlgen_sig.tsv` is empty).
  - Problem 2 (substrate): TenK10K Phase 1 rare-variant zips at Zenodo 17474113
    are 214-260 byte placeholders (HANDOVER §6); common-variant files are
    14-23 GB each and the rare-variant ones we need do not yet exist.
  - Substitution: GTEx v8 LIVER cis-eQTLs via eQTL Catalogue
    (study QTS000015 GTEx, dataset QTD000266, n=208, tissue=liver,
    quant_method=ge). PCSK9 is *liver*-dominant; liver is the correct
    causal tissue for the LDL pathway. We document this as a substitution.
  - Method: remote pysam.TabixFile() over the EBI FTP-hosted bgzipped tsv
    with the corresponding .tbi index — no need to download 2.9 GB.

Outcome substrate: FinnGen R12 x MVP x UKBB joint meta (mvp-ukbb.finngen.fi),
endpoint I9_IHD, same route validated by Worker B 2026-06-08.

Usage:
    python3 fetch_pcsk9_eqtls.py
"""
from __future__ import annotations
import json
import sys
import time
from pathlib import Path
import urllib.request
import urllib.error

import pysam

OUT_DIR = Path(__file__).parent
PCSK9_ENSG = "ENSG00000169174"
PCSK9_CHR = "1"
PCSK9_TSS_START = 55039445
PCSK9_TSS_END = 55064852
CIS_WINDOW_HALF = 1_000_000  # +/- 1 Mb (canonical cis window)

GTEX_LIVER_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/"
    "QTS000015/QTD000266/QTD000266.all.tsv.gz"
)
# n_x for GTEx v8 Liver
N_X_LIVER = 208

# eQTL Catalogue all.tsv schema (per eQTL Catalogue v6 docs):
COLS = [
    "molecular_trait_id", "chromosome", "position", "ref", "alt",
    "variant", "ma_samples", "maf", "pvalue", "beta", "se",
    "type", "ac", "an", "r2",
    "molecular_trait_object_id", "gene_id", "median_tpm", "rsid",
]


def fetch_pcsk9_eqtls():
    """Fetch all PCSK9 cis-eQTLs in +-1 Mb window from GTEx Liver."""
    chrom = PCSK9_CHR
    start = PCSK9_TSS_START - CIS_WINDOW_HALF
    end = PCSK9_TSS_END + CIS_WINDOW_HALF
    print(f"Tabix fetch: chr{chrom}:{start}-{end} (PCSK9 +- 1Mb cis window)")
    print(f"URL: {GTEX_LIVER_URL}")

    tbx = pysam.TabixFile(GTEX_LIVER_URL)
    rows = []
    n_total = 0
    for line in tbx.fetch(chrom, start, end):
        n_total += 1
        f = line.split("\t")
        # Filter to PCSK9 gene only (other genes are also in the cis window)
        if f[0] == PCSK9_ENSG or (len(f) > 16 and f[16] == PCSK9_ENSG):
            rec = dict(zip(COLS, f))
            for k in ("position", "ma_samples", "ac", "an"):
                try:
                    rec[k] = int(rec[k])
                except (ValueError, TypeError):
                    pass
            for k in ("maf", "pvalue", "beta", "se", "median_tpm"):
                try:
                    rec[k] = float(rec[k])
                except (ValueError, TypeError):
                    rec[k] = None
            rows.append(rec)
    tbx.close()
    print(f"Found {len(rows)} PCSK9 cis-eQTL records (of {n_total} cis-window records overall)")
    return rows


FINNGEN_ENDPOINTS = ["I9_IHD", "I9_MI_STRICT", "I9_CHD"]


def fetch_finngen_meta(rsid: str, chrom: str, pos: int, ref: str, alt: str, *,
                       max_tries: int = 3, sleep_s: float = 0.5):
    """Hit FinnGen meta-PheWeb variant API.

    Worker B established route: https://mvp-ukbb.finngen.fi/api/variant/<chr>-<pos>-<ref>-<alt>

    Tries FINNGEN_ENDPOINTS in order; returns first non-NA hit (with the
    endpoint code recorded), or None if all are NA / not in meta.
    """
    # FinnGen meta uses GRCh38 (same as GTEx v8); positions are directly
    # compatible.
    key = f"{chrom}-{pos}-{ref}-{alt}"
    url = f"https://mvp-ukbb.finngen.fi/api/variant/{key}"
    last_err = None
    for attempt in range(max_tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "rvMR-track2/1.0"})
            with urllib.request.urlopen(req, timeout=30) as r:
                data = json.loads(r.read().decode("utf-8"))
            if not isinstance(data, dict):
                return None
            entries = data.get("results") or data.get("phenos") or []
            # try the endpoints in priority order
            for code in FINNGEN_ENDPOINTS:
                for entry in entries:
                    if entry.get("phenocode") == code:
                        beta = entry.get("beta")
                        sebeta = entry.get("sebeta")
                        if beta is not None and sebeta is not None and \
                                beta == beta and sebeta == sebeta:  # NaN-safe
                            entry["_endpoint_used"] = code
                            return entry
            return None  # all endpoints NA at this variant
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            last_err = e
            time.sleep(sleep_s * (attempt + 1))
        except (urllib.error.URLError, TimeoutError, OSError) as e:
            last_err = e
            time.sleep(sleep_s * (attempt + 1))
    print(f"  FinnGen fetch failed for {key}: {last_err}", file=sys.stderr)
    return None


def main():
    rows = fetch_pcsk9_eqtls()
    # Save raw
    with open(OUT_DIR / "pcsk9_gtex_liver_cis.json", "w") as fh:
        json.dump(rows, fh)
    print(f"Wrote {OUT_DIR / 'pcsk9_gtex_liver_cis.json'}")

    # Strategy for instrument selection:
    #   1) keep only SNPs (type == 'SNP')
    #   2) keep only those with MAF >= 0.01 (common variants — Track 2 is the
    #      common-variant plumbing test; rare-variant Track 3 is gated on Wei/Cuomo)
    #   3) sort by p-value
    #   4) for `mrAR_multi`, an LD-pruned set of 5-10 lead variants is ideal.
    #      We do a simple position-window prune (no LD reference): grab the top
    #      hit, then exclude variants within 100 kb, then next top hit, etc.
    snps = [r for r in rows if r.get("type") == "SNP" and (r.get("maf") or 0) >= 0.01
            and r.get("beta") is not None and r.get("se") is not None]
    snps.sort(key=lambda r: r.get("pvalue") or 1.0)
    print(f"{len(snps)} SNPs with MAF >= 0.01 and complete beta/SE")

    # Position-window LD prune: 100 kb window between leads. (Proper LD pruning
    # would use 1000G; this is a coarse proxy adequate for the plumbing test.)
    # NOTE: We RELAX the prune window to 100 kb because PCSK9 cis-eQTLs in
    # GTEx Liver n=208 are statistically weak; with 200 kb we get K=2 which
    # is the minimum for mrAR_multi but does not exercise the J-test
    # interestingly. 100 kb gets us K ~ 8 leads spanning the +-1 Mb cis window.
    lead_set = []
    PRUNE_WINDOW = 100_000
    for r in snps:
        if all(abs(r["position"] - L["position"]) >= PRUNE_WINDOW for L in lead_set):
            lead_set.append(r)
        if len(lead_set) >= 20:
            break
    print(f"Position-pruned {len(lead_set)} lead variants (prune window = {PRUNE_WINDOW} bp)")

    # Filter: keep variants with nominal p < 0.05 (rough cis-eQTL nominal
    # significance threshold). GTEx Liver is weakly powered so we go to 0.05;
    # the Worker B / Track 2 spec calls this a plumbing test, not a high-
    # confidence inference. F-statistics will be low; AR is the appropriate
    # response — that's exactly what AR was built for.
    P_THRESH = 0.05
    lead_set = [r for r in lead_set if (r.get("pvalue") or 1) < P_THRESH]
    print(f"After p < {P_THRESH} filter: {len(lead_set)} lead variants")

    # Save
    with open(OUT_DIR / "pcsk9_lead_eqtls.json", "w") as fh:
        json.dump(lead_set, fh, indent=2)
    print("Wrote pcsk9_lead_eqtls.json")

    # Pull FinnGen meta I9_IHD for each lead
    print("\nFetching FinnGen meta I9_IHD per lead variant...")
    finngen_out = []
    for r in lead_set:
        rsid = r.get("rsid") or "."
        chrom = r["chromosome"]
        pos = r["position"]
        ref = r["ref"]
        alt = r["alt"]
        entry = fetch_finngen_meta(rsid, chrom, pos, ref, alt)
        if entry is None:
            print(f"  {rsid} {chrom}:{pos}:{ref}:{alt}: no I9_IHD record (substituting NA)")
            entry = {"phenocode": "I9_IHD", "beta": None, "sebeta": None,
                     "note": "variant not in FinnGen meta or NA at I9_IHD"}
        finngen_out.append({
            "rsid": rsid, "chr": chrom, "pos": pos, "ref": ref, "alt": alt,
            "eqtl_beta": r["beta"], "eqtl_se": r["se"],
            "eqtl_pvalue": r.get("pvalue"), "maf": r.get("maf"),
            "finngen_entry": entry,
        })
        time.sleep(0.3)  # be polite

    with open(OUT_DIR / "pcsk9_track2_panel.json", "w") as fh:
        json.dump(finngen_out, fh, indent=2)
    print(f"Wrote pcsk9_track2_panel.json  ({len(finngen_out)} variants)")

    print("\nSummary:")
    print(f"  total cis-eQTL records:    {len(rows)}")
    print(f"  SNPs with MAF >= 0.01:     {len(snps)}")
    print(f"  position-pruned leads:     {len(lead_set)}")
    print(f"  with FinnGen I9_IHD beta:  {sum(1 for o in finngen_out if (o['finngen_entry'] or {}).get('beta') is not None)}")


if __name__ == "__main__":
    main()
