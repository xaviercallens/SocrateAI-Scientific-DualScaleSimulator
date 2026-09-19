#!/usr/bin/env python
"""STEP 0: write every expected result, with its source, BEFORE any computation.
Expected values here come from theory (Kunneth over a field, transfer, Lefschetz,
Kummer construction); the computations in A*/B* scripts do not import this file's
numbers into their chain complexes -- they only compare against them afterwards.
Run: python write_expectations.py   (writes expectations.json next to this file)
"""
import json, os
from math import comb

def kun(a, b):
    """Kunneth over a field: Betti of product from Betti of factors (THEORY, for expectations only)."""
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out

def chi(b):
    return sum((-1) ** k * x for k, x in enumerate(b))

T2 = [1, 2, 1]; S2 = [1, 0, 1]
ORB = [1, 0, 6, 0, 1]           # T^4/Z2 over odd p (transfer)
K3 = [1, 0, 22, 0, 1]           # K3, torsion-free (Kummer; Barth-Hulek-Peters-Van de Ven)
KLEIN_F2 = [1, 2, 1]; KLEIN_F3 = [1, 1, 0]

E = {
  "meta": {
    "written_before_computation": True,
    "date": "2026-09-19",
    "tiers": {"B": "exact arithmetic over Z/p with negative controls", "L": "literature/theory input", "X": "numerics"},
    "note": "Expected values are theory. Kunneth appears ONLY here (expectations), never in the computed chain complexes.",
  },
  "partA": {
    "A1_tori": {
      "expected": {f"T^{n}": [comb(n, k) for k in range(n + 1)] for n in (2, 4, 6)},
      "fields": [2, 3, 5], "grid_N": [2, 3, 4],
      "source": "H_k(T^n; F) = F^binomial(n,k) for every field (exterior algebra on H^1); N-independence of cellular homology",
      "chi": 0,
    },
    "engine_controls": {
      "RP2": {"F2": [1, 1, 1], "F3": [1, 0, 0], "F5": [1, 0, 0]},
      "RP3": {"F2": [1, 1, 1, 1], "F3": [1, 0, 0, 1], "F5": [1, 0, 0, 1]},
      "S3": [1, 0, 0, 1],
      "source": "H_*(RP^n;Z): Z, Z/2, ... (Hatcher Ex. 2.42); field-sensitivity control for the rank engine",
    },
    "A2_orbifold_T4_mod_Z2": {
      "expected_odd_p": {"F3": ORB, "F5": ORB},
      "source_odd_p": "transfer: H(X/G;F) = H(X;F)^G for char F not dividing |G|; v->-v acts on H_k(T^4) by (-1)^k, so invariants = binomial(4,k) for k even, 0 for k odd",
      "expected_F2": "NOT pre-specified (2-torsion from RP^3 links changes it); only constraint: chi = 8",
      "chi": 8,
      "chi_source": "chi(X/G) = (chi(X) + L(sigma))/2 = (0 + 16)/2 = 8; L(sigma)=16 from 16 isolated fixed points of index +1",
      "lefschetz_number_sigma_on_T4": 16,
      "per_degree_action_on_H_k(T4;F3)": "(-1)^k * identity",
      "fixed_cells": "exactly 16 fixed cells (the vertices with v_i in {0,N/2}); no positive-dimensional cell fixed setwise for N even",
    },
    "A3_orbifold_x_T2": {
      "expected": {"F3": kun(ORB, T2), "F5": kun(ORB, T2)},
      "chi": 0,
      "source": "Kunneth over a field (expectation only; the computation is a direct rank computation on the product cell complex)",
    },
    "A4_resolved_K3_and_K3xT2": {
      "K3": {"F2": K3, "F3": K3, "F5": K3, "chi": 24},
      "K3xT2": {"F2": kun(K3, T2), "F3": kun(K3, T2), "F5": kun(K3, T2), "chi": 0},
      "tier_L_inputs": [
        "D(O(-2)) is the mapping cylinder of its circle-bundle projection pi: RP^3=L(2,1) -> S^2 (disc bundle = mapping cylinder of sphere-bundle projection)",
        "pi^*[S^2] generates H^2(RP^3;Z)=Z/2 (Gysin sequence, Euler number +-2)",
      ],
      "local_model_checks": {
        "H_*(Cyl(pi), L; F2)": [0, 0, 1, 0, 1],
        "H_*(Cyl(pi), L; F3)": [0, 0, 1, 0, 1],
        "source": "Lefschetz duality for the oriented 4-manifold D(O(-2)) ~ S^2 with boundary RP^3: H_k(D,dD) = H^{4-k}(D) = Z for k=2,4",
        "wrong_class_control_phi2_eq_0": {"H_*(Cyl,L;F2)": [0, 0, 2, 1, 1], "H_*(Cyl,L;F3)": [0, 0, 1, 0, 1]},
      },
      "cone_reassembly_control": "U union 16 cones on the links must reproduce A2 (same field, same N)",
    },
    "A5_controls": {
      "partial_resolution_X_k": {
        str(k): {"X_k_F3": [1, 0, 6 + k, 0, 1], "X_k_x_T2_F3": kun([1, 0, 6 + k, 0, 1], T2),
                 "chi_X_k": 8 + k, "chi_X_k_x_T2": 0}
        for k in (0, 8, 15, 16)
      },
      "partial_resolution_F2": "not pre-specified except k=16 (K3: (1,0,22,0,1)) and chi(X_k)=8+k",
      "wrong_gluing_class_phi2_eq_0_all16": {
        "F3": K3,
        "F2": "must DIFFER from (1,0,22,0,1); predicted = b(T^4/Z2;F2) + 16 in degree 2 (homologically T^4/Z2 wedge 16 S^2)",
        "meaning": "over odd p the resolution cannot be distinguished from 'orbifold wedge 16 spheres'; only F2 discriminates",
      },
      "S2_instead_of_T2": {
        "orbifold_x_S2_F3": kun(ORB, S2), "K3_x_S2_F3": kun(K3, S2), "K3_x_S2_F2": kun(K3, S2),
      },
      "Klein_instead_of_T2": {
        "Klein": {"F2": KLEIN_F2, "F3": KLEIN_F3},
        "K3_x_Klein": {"F2": kun(K3, KLEIN_F2), "F3": kun(K3, KLEIN_F3)},
        "orbifold_x_Klein_F3": kun(ORB, KLEIN_F3),
        "note": "K3 x Klein over F2 equals K3 x T2 over F2 (1,2,23,44,23,2,1); only F3 separates them",
      },
      "failure_rule": "a route that reports (1,2,23,44,23,2,1) regardless of k FAILS",
    },
  },
  "partB": {
    "tier": "X",
    "criterion_recovered": {
      "definition": "Using the barcode over field F, beta_k(eps) = #H_k bars with birth <= eps < death. The Betti vector is RECOVERED at a sample if there is a window [e1,e2] with e2 >= R*e1, R = 1.5, on which beta_k(eps) equals the expected value for every tested k simultaneously.",
      "R": 1.5,
      "N_min_rule": "smallest N on the tested grid at which ALL 3 seeds (0,1,2) satisfy the criterion",
      "also_reported": "for each k: ratio of persistence of the b_k-th longest bar to the (b_k+1)-th longest bar",
      "no_lowering": "if the criterion fails at the largest feasible N, report failure; R is not lowered after seeing data",
      "rips_reliability": "Rips expanded to dim D gives reliable H_k only for k <= D-1",
      "fields": [2, 3],
    },
    "B1_flat_tori": {
      "T^2 in R^4": [1, 2, 1], "T^3 in R^6": [1, 3, 3, 1], "T^4 in R^8": [1, 4, 6, 4, 1],
      "prior_expectation": "T^2 recovers easily; T^3 with effort; T^4 (needs 5-simplices) likely infeasible within budget",
      "fit": "N_min vs dimension d, fit N_min = A * B^d on the recovered cases",
    },
    "B2_T4_mod_Z2": {"expected_H0_H2_F3": [1, 0, 6], "expected_F2": "not pre-specified (torsion)"},
    "B3_K3_Fermat_quartic": {"expected_H0_H2": [1, 0, 22],
                              "prior_expectation": "not recoverable at feasible N (b2=22 in a 4-manifold with a non-flat metric, 4-dim sampling)",
                              "K3xT2_expected_H0_H2": [1, 2, 23]},
    "B4_null_random_4d": {"expected_H0_H2": [1, 0, 0],
                           "false_positive_rule": "any window with beta_2>=1 of ratio >= 1.5 is a false positive"},
  },
}
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "expectations.json")
with open(out, "w") as f:
    json.dump(E, f, indent=1)
print("wrote", out)
