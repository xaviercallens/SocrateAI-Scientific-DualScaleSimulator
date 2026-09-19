"""
Shared provenance reader for K3 topological data (chi, and where available
b2 and the individual Betti numbers), used by both 04_tadpole_budget.py and
06_chained_rigidity.py so both scripts read the SAME external numbers from
the SAME file, rather than each typing or independently sourcing them.

Poll policy (per ground rules): D-tda/exports.json is preferred if it exists
(poll_exports.sh in this directory was run BEFORE these scripts, up to 20
minutes, to give Track D a chance to publish; its log is poll_exports.log).
If D-tda never gains an exports.json, Track A's A-genus/exports.json is used
instead. As of this poll, D-tda/exports.json DID appear (after the initial
20-minute poll had already completed and A-genus had been used once; see
06's result JSON "provenance" field for the exact history) with schema
{"chi_K3_resolved": 24, "b2_K3_resolved": 22, "b0":1,"b1":0,"b3":0,"b4":1,
"source": "03_resolution_hybrid.py, N=6 (premise-valid), field Z/3", ...} --
an INDEPENDENT computation (GUDHI-based resolution of a Kummer-type
orbifold, per its own "source" string) of chi AND b2 AND the individual
Betti numbers, which is richer than A-genus's chi-only export.

chi_top_from_provenance() (kept for 04's simple use) returns only
(chi, source_file_label, source_key, raw_file_contents).

topology_from_provenance() (used by 06) returns the fuller dict, with
b2/b0/b1/b3/b4 as None when the source file doesn't provide them (04/06
must then say explicitly which are read vs. which are named definitional
inputs -- never silently assumed).
"""
import json
import os

D_TDA_EXPORTS = "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity_v2/D-tda/exports.json"
A_GENUS_EXPORTS = "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity_v2/A-genus/exports.json"


def topology_from_provenance():
    """Returns a dict:
      {chi, b2, b0, b1, b3, b4, source_file, source_keys, raw}
    b2/b0/b1/b3/b4 are None if not present in the source file (read as such,
    never invented). chi is always populated or this raises."""
    if os.path.exists(D_TDA_EXPORTS):
        with open(D_TDA_EXPORTS) as f:
            d = json.load(f)
        chi = None
        chi_key = None
        for key in ("chi_K3_resolved", "chi", "chi_top", "euler_characteristic"):
            if key in d:
                chi, chi_key = int(d[key]), key
                break
        if chi is None:
            raise RuntimeError(f"D-tda/exports.json exists but has no recognized chi key: {list(d.keys())}")
        b2 = None
        b2_key = None
        for key in ("b2_K3_resolved", "b2"):
            if key in d:
                b2, b2_key = int(d[key]), key
                break
        b0 = int(d["b0"]) if "b0" in d else None
        b1 = int(d["b1"]) if "b1" in d else None
        b3 = int(d["b3"]) if "b3" in d else None
        b4 = int(d["b4"]) if "b4" in d else None
        return {
            "chi": chi, "b2": b2, "b0": b0, "b1": b1, "b3": b3, "b4": b4,
            "source_file": "D-tda/exports.json",
            "source_keys": {"chi": chi_key, "b2": b2_key},
            "raw": d,
        }
    if os.path.exists(A_GENUS_EXPORTS):
        with open(A_GENUS_EXPORTS) as f:
            d = json.load(f)
        if "chi_from_genus" not in d:
            raise RuntimeError(f"A-genus/exports.json exists but has no chi_from_genus key: {list(d.keys())}")
        return {
            "chi": int(d["chi_from_genus"]), "b2": None, "b0": None, "b1": None, "b3": None, "b4": None,
            "source_file": "A-genus/exports.json",
            "source_keys": {"chi": "chi_from_genus", "b2": None},
            "raw": d,
        }
    raise RuntimeError(
        "Neither D-tda/exports.json nor A-genus/exports.json exists after polling; "
        "cannot derive chi_top without typing a literal value (forbidden by ground rules)."
    )


def chi_top_from_provenance():
    """Simplified accessor (chi only), for 04_tadpole_budget.py."""
    topo = topology_from_provenance()
    return topo["chi"], topo["source_file"], topo["source_keys"]["chi"], topo["raw"]
