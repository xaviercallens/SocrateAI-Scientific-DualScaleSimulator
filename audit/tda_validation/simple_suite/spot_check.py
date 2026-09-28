#!/usr/bin/env python3
"""Determinism spot-check: cases run before the committed script (no script_sha256 key)
versus a rerun of one representative per group at the committed script (spot/*.json).
Timing and provenance fields are ignored. Writes spot/spot_check.json."""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
SKIP = {"wall_sec", "wall_sec_total", "peak_rss_mb", "loadavg_start", "started", "git_head", "command", "runtime_sec",
        "pipeline_runtime_sec", "build_sec", "persistence_wall_sec", "topology_sec", "landmark_sec", "status",
        "script_sha256", "code_under_test_sha256", "versions", "count"}


def strip(x):
    if isinstance(x, dict):
        return {k: strip(v) for k, v in x.items() if k not in SKIP}
    if isinstance(x, list):
        return [strip(v) for v in x]
    return x


out = {"pre_commit_case_files": sorted(fn for fn in os.listdir(os.path.join(HERE, "cases"))
                                       if fn.endswith(".json") and "script_sha256" not in json.load(open(os.path.join(HERE, "cases", fn)))),
       "comparisons": {}}
for fn in ("P7_N4000.json", "P6_pipeline_N1000.json"):
    a = strip(json.load(open(os.path.join(HERE, "cases", fn))))
    b = strip(json.load(open(os.path.join(HERE, "spot", fn))))
    out["comparisons"][fn] = {"identical_ignoring_timing": a == b}
a = json.load(open(os.path.join(HERE, "cases", "F3chunk_000.json")))["maps"][:5]
b = json.load(open(os.path.join(HERE, "spot", "F3chunk_000.json")))["maps"]
out["comparisons"]["F3chunk_000.json maps 0-4"] = {"identical": a == b, "n_maps_compared": len(b)}
json.dump(out, open(os.path.join(HERE, "spot", "spot_check.json"), "w"), indent=1)
print(json.dumps(out["comparisons"]), len(out["pre_commit_case_files"]), "pre-commit case files")
