"""Fetch protein and RNA structures from RCSB into the data lake (never into git).

Writes  <BIO_DATA>/pdb/<ID>.pdb  and  <BIO_DATA>/rna/<ID>.pdb  plus a manifest with
the exact URL, sha256, byte size and fetch time for every file.  A fetch that fails
three times is recorded as ABSENT with the URLs tried.

The entry list spans folds deliberately: all-alpha, all-beta (including barrels),
alpha/beta, alpha+beta and a few large multi-domain chains, plus RNA for contrast.
Fold labels below are the LITERATURE expectation; the fold class actually used in
the separation test is derived by this pipeline from the HELIX/SHEET records
(see expectations.json), not from these comments.

Usage: python fetch_pdb.py
"""
from __future__ import annotations

import json
import sys
import time
import urllib.error
import urllib.request

sys.path.insert(0, "/mnt/disks/disk-socrateai-local-1/wt-topo-bio/topodb_runs/biology")
from bio_common import BIO_DATA, sha256  # noqa: E402

URL = "https://files.rcsb.org/download/{}.pdb"

PROTEIN = {
    # --- expected all-alpha ---
    "256B": "cytochrome b562, four-helix bundle",
    "1BZ1": "haemoglobin-like, all-alpha",
    "2MHR": "myohemerythrin, four-helix bundle",
    "1MBA": "myoglobin (Aplysia), globin fold",
    "1LMB": "lambda repressor, helix-turn-helix",
    "1ENH": "engrailed homeodomain, three helices",
    "2ZTA": "GCN4 leucine zipper, coiled coil",
    "1ROP": "ROP protein, four-helix bundle dimer",
    "1A1Z": "designed alpha-helical",
    "3ICB": "calbindin D9k, EF-hand all-alpha",
    "1PRB": "albumin-binding domain, three-helix",
    "2LZM": "T4 lysozyme, mostly alpha",
    # --- expected all-beta ---
    "1TEN": "tenascin fibronectin type III, beta sandwich",
    "2RH1": "beta2-adrenergic receptor (has a 7TM helical bundle; kept as a stated check)",
    "1BTP": "trypsin, beta barrel protease",
    "2POR": "porin, 16-stranded beta barrel",
    "1QJ8": "OmpX, 8-stranded beta barrel",
    "1BYM": "green fluorescent protein, 11-stranded beta barrel",
    "1EMA": "GFP, beta barrel",
    "1PGB": "protein G B1 domain, alpha+beta",
    "2GB1": "protein G B1 (NMR)",
    "1SHG": "SH3 domain, beta barrel-like",
    "1TTF": "tenth fibronectin type III, beta sandwich",
    "1CBI": "cellular retinol-binding protein, beta barrel (lipocalin-like)",
    "1IFC": "intestinal fatty-acid binding protein, 10-stranded beta barrel",
    "1OPA": "retinol-binding protein, lipocalin beta barrel",
    # --- expected alpha/beta and alpha+beta ---
    "1RNB": "barnase, alpha+beta",
    "3CHY": "CheY, Rossmann-like alpha/beta",
    "1PGA": "protein G, alpha+beta",
    "5CYT": "cytochrome c",
    "1TIM": "triosephosphate isomerase, TIM barrel alpha/beta",
    "4AKE": "adenylate kinase, alpha/beta",
    "1AKE": "adenylate kinase closed form",
    # --- large / multi-domain ---
    "1AON": "GroEL-GroES chaperonin (large, multi-chain)",
    "1HTM": "influenza haemagglutinin HA2 fragment",
    "1GFL": "GFP dimer",
}

RNA = {
    "1EHZ": "yeast tRNA-Phe, 76 nt",
    "2GIS": "SAM-I riboswitch",
    "1GID": "Tetrahymena group I intron P4-P6 domain",
    "1FFK": "50S ribosomal subunit (large RNA, Haloarcula)",
    "1Y26": "adenine riboswitch aptamer",
    "3DIL": "lysine riboswitch",
}


def fetch(pid: str, kind: str, records: list, absent: list) -> None:
    out = BIO_DATA / kind / f"{pid}.pdb"
    url = URL.format(pid)
    if out.exists() and out.stat().st_size > 1000:
        records.append({"id": pid, "kind": kind, "url": url, "path": str(out),
                        "bytes": out.stat().st_size, "sha256": sha256(out), "cached": True})
        return
    tried = []
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=60) as fh:
                data = fh.read()
            if len(data) < 1000:
                raise ValueError(f"suspiciously short: {len(data)} bytes")
            out.write_bytes(data)
            records.append({"id": pid, "kind": kind, "url": url, "path": str(out),
                            "bytes": len(data), "sha256": sha256(out), "cached": False,
                            "fetched_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())})
            return
        except (urllib.error.URLError, urllib.error.HTTPError, ValueError, OSError) as exc:
            tried.append(f"attempt {attempt + 1}: {type(exc).__name__}: {exc}")
            time.sleep(2)
    absent.append({"id": pid, "kind": kind, "urls_tried": [url] * 3, "errors": tried,
                   "status": "ABSENT"})
    print(f"  ABSENT {pid}: {tried[-1]}")


def main() -> int:
    records: list = []
    absent: list = []
    for pid in PROTEIN:
        fetch(pid, "pdb", records, absent)
    for pid in RNA:
        fetch(pid, "rna", records, absent)

    # the three locally verified-genuine files, re-hashed here (byte-identical to RCSB
    # per the earlier validation; we record their own sha256 from this disk)
    local = []
    for pid in ("1CRN", "1UBQ", "4OBE"):
        p = BIO_DATA.parent / "genetics" / "rcsb_reference" / f"{pid}.pdb"
        if p.exists():
            local.append({"id": pid, "kind": "pdb_local_verified", "url": URL.format(pid),
                          "path": str(p), "bytes": p.stat().st_size, "sha256": sha256(p),
                          "note": "verified byte-identical to an RCSB download by the earlier "
                                  "genetics validation run (local_datasets_assessment.json)"})

    man = BIO_DATA / "pdb_manifest.json"
    man.write_text(json.dumps({
        "fetched_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "url_template": URL, "n_ok": len(records), "n_absent": len(absent),
        "entries": records, "local_verified": local, "absent": absent,
        "literature_labels_NOT_used_for_classification": {**PROTEIN, **RNA},
    }, indent=1))
    print(f"fetched/cached {len(records)} + {len(local)} local, ABSENT {len(absent)} -> {man}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
