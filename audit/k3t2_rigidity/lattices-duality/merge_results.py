"""Merge all per-item result JSON files in this directory into one results.json."""
import json
import os

DIR = "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity/lattices-duality"

FILES = [
    "01_e8_result.json",
    "02_k3_mukai_gamma_result.json",
    "03_oddz_checks_result.json",
    "04_dual_scale_bound_result.json",
    "05_tadpole_arithmetic_result.json",
    "06_rigidity_a_result.json",
    "07_rigidity_b_result.json",
]

merged = {
    "track": "C",
    "title": "K3 x T2 lattice/T-duality mathematics: exact-arithmetic checks and zero-free-parameter rigidity tests",
    "python": "/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python",
    "worktree": "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2",
    "branch": "loop/k3t2-rigidity",
    "items": {},
    "could_not_do": [],
}

for fname in FILES:
    path = os.path.join(DIR, fname)
    with open(path) as f:
        merged["items"][fname] = json.load(f)

out_path = os.path.join(DIR, "results.json")
with open(out_path, "w") as f:
    json.dump(merged, f, indent=2)

print("wrote", out_path)
