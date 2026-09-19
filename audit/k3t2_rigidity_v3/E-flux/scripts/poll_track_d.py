#!/usr/bin/env python3
"""Poll Track D's exports for the number of T^4/Z2 singular points (declared shared input,
NOT typed here). Polls audit/k3t2_rigidity_v3/D-tda/exports.json for up to POLL_SECONDS
(argument 1, default 1200); if absent, falls back to audit/k3t2_rigidity_v2/D-tda/exports.json
and records the fallback.  Also computes the same count independently (tier B) as the number of
half-lattice fixed points of x -> -x on (R/Z)^n, n=4 (and n=2 for T2/Z2).

Run: cd audit/k3t2_rigidity_v3/E-flux/scripts && <venv python> poll_track_d.py 1200
"""
import sys, time, json, itertools
from fractions import Fraction
from pathlib import Path
from common import *

def fixed_points(n):
    """Points x in (R/Z)^n with x = -x mod 1, enumerated exactly from a candidate grid of denominator 4
    (superset of all solutions: 2x in Z^n forces x in (1/2)Z/Z, grid contains them; the check is exact)."""
    grid = [Fraction(k, 4) for k in range(4)]
    pts = []
    for x in itertools.product(grid, repeat=n):
        if all((2 * xi) % 1 == 0 for xi in x):
            pts.append(x)
    return len(pts)

def find_key(obj, sub):
    """keys containing 'fixed_point' or 'n_singular' (Track D v3 names the count n_singular; v2 num_fixed_points_computed)"""
    if sub == "fixed_point":
        return find_key1(obj, "fixed_point") or find_key1(obj, "n_singular")
    return find_key1(obj, sub)

def find_key1(obj, sub):
    hits = []
    if isinstance(obj, dict):
        for k, v in obj.items():
            if sub in k and isinstance(v, (int, float)):
                hits.append((k, v))
            hits += find_key1(v, sub)
    elif isinstance(obj, list):
        for v in obj:
            hits += find_key1(v, sub)
    return hits

def main():
    secs = int(sys.argv[1]) if len(sys.argv) > 1 else 1200
    v3 = ROOT / "audit" / "k3t2_rigidity_v3" / "D-tda" / "exports.json"
    v2 = ROOT / "audit" / "k3t2_rigidity_v2" / "D-tda" / "exports.json"
    t0 = time.time()
    used = None
    while True:
        if v3.exists():
            try:
                d = json.load(open(v3)); hits = find_key(d, "fixed_point")
                if hits:
                    used = v3; break
            except Exception:
                pass
        if time.time() - t0 > secs:
            break
        time.sleep(15)
    fallback = used is None
    if fallback:
        used = v2
        d = json.load(open(v2)); hits = find_key(d, "fixed_point")
    out = {
        "track_d_export_used": rel(used),
        "v3_export_present_after_polling": not fallback,
        "polled_seconds": round(time.time() - t0),
        "fallback_to_v2": fallback,
        "keys_found": hits,
        "track_d_other_exports_used": {k: json.load(open(used)).get(k) for k in ("chi", "b2")},
        "T4Z2_fixed_points_from_track_D": hits[0][1] if hits else None,
        "T4Z2_fixed_points_computed_here_tierB": fixed_points(4),
        "T2Z2_fixed_points_computed_here_tierB": fixed_points(2),
        "sha256_of_export_used": sha256(used),
        "command": rel_command("poll_track_d.py", str(secs)),
    }
    out["agree_D_vs_here"] = out["T4Z2_fixed_points_from_track_D"] == out["T4Z2_fixed_points_computed_here_tierB"]
    write_json("00_track_d_poll.json", out)
    print(json.dumps(out, indent=1))

if __name__ == "__main__":
    main()
