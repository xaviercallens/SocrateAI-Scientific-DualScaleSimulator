#!/usr/bin/env python
"""
Track D, step (4): Betti numbers and Euler characteristic of K3 x T^2, by
the Kunneth formula (field coefficients, so no Tor terms -- exact), applied
to the Betti numbers COMPUTED in 03_resolution_hybrid.py (K3) and
01_controls.py (T^2). Nothing here is typed from memory; every input comes
from a JSON file this track already wrote, and the arithmetic is plain
integer convolution done with Python ints (exact).
"""
import json
import sys

sys.path.insert(0, "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity/tda-gudhi")

with open("/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity/tda-gudhi/03_resolution_hybrid_results.json") as f:
    step3 = json.load(f)
with open("/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity/tda-gudhi/01_controls_results.json") as f:
    step1 = json.load(f)

b_T2 = step1["cubical"]["T2"]["betti"]  # [1,2,1], computed by gudhi PeriodicCubicalComplex

OUT = {"input_b_T2": b_T2, "input_b_T2_source": "01_controls_results.json / cubical.T2.betti"}

for N in ("N=4", "N=6"):
    for field in ("Z3", "Z5"):
        key = f"resolved_K3_betti_{field}"
        r = step3[N][key]
        b_K3 = [r["b0"], r["b1"], r["b2"], r["b3"], r["b4"]]

        # Kunneth over a field: b_k(K3 x T2) = sum_{i+j=k} b_i(K3) * b_j(T2)
        max_dim = (len(b_K3) - 1) + (len(b_T2) - 1)
        b_prod = [0] * (max_dim + 1)
        for i, bi in enumerate(b_K3):
            for j, bj in enumerate(b_T2):
                b_prod[i + j] += bi * bj

        chi_K3 = sum((-1) ** i * bi for i, bi in enumerate(b_K3))
        chi_T2 = sum((-1) ** j * bj for j, bj in enumerate(b_T2))
        chi_prod_from_betti = sum((-1) ** k * b for k, b in enumerate(b_prod))
        chi_prod_multiplicative = chi_K3 * chi_T2  # independent check: chi is multiplicative

        OUT[f"{N}_{field}"] = {
            "input_b_K3": b_K3,
            "input_b_K3_source": f"03_resolution_hybrid_results.json / {N} / {key}",
            "betti_K3xT2_via_kunneth": b_prod,
            "chi_K3": chi_K3,
            "chi_T2": chi_T2,
            "chi_K3xT2_from_product_betti": chi_prod_from_betti,
            "chi_K3xT2_from_multiplicativity_chiK3_times_chiT2": chi_prod_multiplicative,
            "chi_cross_check_matches": chi_prod_from_betti == chi_prod_multiplicative,
        }

OUT_PATH = "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity/tda-gudhi/05_kunneth_results.json"
with open(OUT_PATH, "w") as f:
    json.dump(OUT, f, indent=2)

print(json.dumps(OUT, indent=2))
print("\nWrote", OUT_PATH)
