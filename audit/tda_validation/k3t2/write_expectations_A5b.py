#!/usr/bin/env python
"""Addendum to STEP 0, written AFTER A5 ran and BEFORE A5b is computed (post-A5, pre-A5b).
Motivation: A2/A5 found 5 cyclic 2-primary summands in H_2(T^4/Z2;Z) (from F2 vs F3 Betti).
Theory (tier L): with J = set of CONED singular points (the others resolved), X_J = K3 with the
(-2)-spheres E_i (i in J) contracted, so H_2(X_J;Z) = H_2(K3;Z)/<E_i : i in J> and H_3 = 0.
Its torsion is (saturation of <E_J>)/<E_J> = the subcode of the Kummer code (first-order
Reed-Muller code RM(1,4) on the 16 points = F_2^4, codewords 0, 30 affine hyperplanes, all 16)
supported on J  [Nikulin 1975; Barth-Hulek-Peters-Van de Ven, Compact Complex Surfaces, VIII].
Hence b(X_J;F2) = (1, 0, 6+k+t, t, 1) with k = 16-|J| resolved, t = dim of that subcode.
Points are labelled by bits p_i/(N/2) in F_2^4.
"""
import json, os
E = {"written_before_computing_A5b": True, "tier_of_prediction": "L (Kummer lattice / RM(1,4) code) + exact derivation",
     "cases": {
       "J_all16_coned_k0": {"t": 5, "F2": [1, 0, 11, 5, 1], "note": "already seen in A2/A5 (not a new test)"},
       "J_affine_hyperplane_x0_eq_0_k8": {"t": 1, "F2": [1, 0, 15, 1, 1], "F3": [1, 0, 14, 0, 1]},
       "J_non_hyperplane_8set_k8": {"t": 0, "F2": [1, 0, 14, 0, 1], "F3": [1, 0, 14, 0, 1]},
       "resolved_affine_2plane_k4": {"t": 2, "F2": [1, 0, 12, 2, 1], "F3": [1, 0, 10, 0, 1]},
       "resolved_4_affinely_independent_k4": {"t": 1, "F2": [1, 0, 11, 1, 1], "F3": [1, 0, 10, 0, 1]},
       "resolved_single_point_k1": {"t": 4, "F2": [1, 0, 11, 4, 1], "F3": [1, 0, 7, 0, 1], "why": "RM(1,4) words vanishing at one point: kernel of evaluation, dim 5-1 = 4"},
     },
     "product_rule": "each X_J x T^2(Z_2^2) must equal the F2/F3 Betti vector of X_J tensored with (1,2,1) (field Kunneth, expectation only)"}
json.dump(E, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "expectations_A5b.json"), "w"), indent=1)
print("ok")
