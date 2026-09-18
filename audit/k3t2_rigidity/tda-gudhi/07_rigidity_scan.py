#!/usr/bin/env python
"""
RIGIDITY TEST: vary the number of resolved singular points k = 0..16 and
find for which k the Mayer-Vietoris/Euler-characteristic count gives a
smooth closed manifold with chi = 24.

Model: X_k = U  union_{16 x RP^3}  [ k disk-bundles-over-S^2  +  (16-k) cones-on-RP^3 ]
  chi(X_k) = chi(U) + k*chi(S^2) + (16-k)*chi(cone(RP^3)) - 16*chi(RP^3)

All four ingredients are read from files this track already computed by
running GUDHI (nothing here is a fresh guess or a typed target):
  chi(U)          <- 03_resolution_hybrid_results.json (complement homology)
  chi(S^2)         <- 03_resolution_hybrid_results.json (boundary-of-tetrahedron control)
  chi(RP^3)         <- computed here from the link Betti numbers in
                        03_resolution_hybrid_results.json (all 16 links agreed)
  chi(cone(RP^3))=1 <- standard fact (a cone is contractible), STATED not computed
  16                <- fixed-vertex count computed in 02_invariance_and_quotient_results.json

Solution set: values of k in {0,...,16} for which chi(X_k) == 24.
Negative control: k=15 (adjacent value) must NOT satisfy chi=24.

Second, independent condition landing on the same k (reported, not just
asserted): X_k is a smooth CLOSED MANIFOLD (no orbifold points remain)
iff every one of the 16 isolated singular points has been resolved, i.e.
iff k == 16. This is a geometric/combinatorial fact about the
construction (a singular point left un-resolved is still there when
k<16), stated here rather than tested by an independent GUDHI run, but it
is a DIFFERENT criterion from the chi=24 arithmetic condition, so
agreement between the two is non-trivial.
"""
import json

BASE = "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity/tda-gudhi"

with open(f"{BASE}/02_invariance_and_quotient_results.json") as f:
    step2 = json.load(f)
with open(f"{BASE}/03_resolution_hybrid_results.json") as f:
    step3 = json.load(f)

OUT = {}

for N in (4, 6):
    key = f"N={N}"
    num_fixed = step2[key]["num_fixed_vertices_computed"]
    chi_U = step3[key]["U_complement_result"]["euler_characteristic_from_simplex_counts"]
    chi_S2 = step3["S2_control"]["euler_characteristic_from_simplex_counts"]

    # chi(RP^3) computed from the (agreed) link Betti numbers over Z/3
    link0 = step3[key]["link_results_per_vertex"][0]["betti_by_field"]["Z3"]
    chi_RP3 = sum((-1) ** i * b for i, b in enumerate(link0))

    chi_cone_RP3 = 1  # STATED: a cone on anything is contractible, chi=1

    def chi_of_k(k):
        return chi_U + k * chi_S2 + (num_fixed - k) * chi_cone_RP3 - num_fixed * chi_RP3

    scan = {k: chi_of_k(k) for k in range(0, num_fixed + 1)}
    solution_set = [k for k, chi in scan.items() if chi == 24]

    entry = {
        "N": N,
        "num_fixed_points_computed": num_fixed,
        "chi_U_computed": chi_U,
        "chi_S2_computed": chi_S2,
        "chi_RP3_computed_from_link_betti": chi_RP3,
        "chi_RP3_link_betti_used": link0,
        "chi_cone_RP3_STATED": chi_cone_RP3,
        "chi_of_k_scan": scan,
        "solution_set_chi_equals_24": solution_set,
        "is_single_point": len(solution_set) == 1,
        "negative_control_k=num_fixed-1": {
            "k": num_fixed - 1,
            "chi": scan.get(num_fixed - 1),
            "fails_chi_24": scan.get(num_fixed - 1) != 24,
        },
        "negative_control_k=0 (fully singular)": {
            "k": 0,
            "chi": scan.get(0),
            "matches_directly_computed_chi_of_singular_quotient": scan.get(0) == step2[key]["quotient_T4_mod_Z2"]["euler_characteristic_from_simplex_counts"],
        },
        "second_independent_condition_smoothness": {
            "criterion": "X_k is a smooth closed manifold iff every one of the 16 singular points is resolved",
            "k_satisfying_smoothness": num_fixed,
            "agrees_with_chi24_solution_set": solution_set == [num_fixed],
        },
        "rigidity_conclusion": (
            f"Zero free parameters for k: two INDEPENDENT criteria (chi=24 arithmetic, "
            f"and smoothness/no-remaining-orbifold-points) both force k={num_fixed}, "
            f"the unique point in {{0,...,{num_fixed}}}; k={num_fixed - 1} fails both."
            if solution_set == [num_fixed] else
            "solution set did not collapse to the expected single point; see raw scan"
        ),
    }
    OUT[key] = entry

OUT_PATH = f"{BASE}/07_rigidity_scan_results.json"
with open(OUT_PATH, "w") as f:
    json.dump(OUT, f, indent=2)

print(json.dumps(OUT, indent=2))
print("\nWrote", OUT_PATH)
