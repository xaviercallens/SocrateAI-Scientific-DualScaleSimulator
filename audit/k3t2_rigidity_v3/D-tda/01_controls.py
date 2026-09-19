#!/usr/bin/env python
"""Track D v3, step 1: known-answer controls (expected values are filled AFTER computing).

Run:  cd audit/k3t2_rigidity_v3/D-tda && /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python 01_controls.py 4 6 8
(arguments = Kuhn grid sizes N for the T^4 control; committed run used 4 6 8)
 (a) T^2, T^4 via gudhi.PeriodicCubicalComplex (Z/3).  (b) T^4 via the Kuhn triangulation, Z/3 and Z/2.
 (c) S^2 = boundary of a tetrahedron.  (d) K3 x T^2 input b(T^2) is taken from (a).
"""
import json, sys, itertools
import gudhi
from lib_common import *

Ns = [int(a) for a in sys.argv[1:]] or [4, 6, 8]
OUT = {"args": Ns}


def cubical(d, cells=3, field=3):
    p = gudhi.PeriodicCubicalComplex(dimensions=[cells] * d, top_dimensional_cells=[0.0] * cells ** d,
                                     periodic_dimensions=[True] * d)
    p.compute_persistence(homology_coeff_field=field, min_persistence=0)
    return list(p.betti_numbers())


OUT["cubical_T2"] = {"betti_Z3": cubical(2), "betti_Z2": cubical(2, field=2)}
OUT["cubical_T4"] = {"betti_Z3": cubical(4), "betti_Z2": cubical(4, field=2)}
OUT["cubical_T2"]["expected_after_computing"] = "C(2,k) = [1,2,1] (standard torus cohomology, tier L)"
OUT["cubical_T4"]["expected_after_computing"] = "C(4,k) = [1,4,6,4,1] (standard torus cohomology, tier L)"
OUT["cubical_T2"]["matches_expected"] = OUT["cubical_T2"]["betti_Z3"] == [1, 2, 1]
OUT["cubical_T4"]["matches_expected"] = OUT["cubical_T4"]["betti_Z3"] == [1, 4, 6, 4, 1]
for N in Ns:
    tops = kuhn_top(N)
    f = fvector(all_faces_by_dim(tops))
    OUT[f"kuhn_T4_N={N}"] = {"betti_Z3": gudhi_betti(tops, 3), "betti_Z2": gudhi_betti(tops, 2),
                             "fvector": f, "chi_from_f": chi_from_f(f)}
S2 = [c for c in itertools.combinations(range(4), 3)]
fs = fvector(all_faces_by_dim(S2))
OUT["S2_tetra_boundary"] = {"betti_Z3": gudhi_betti(S2, 3), "betti_Z2": gudhi_betti(S2, 2), "fvector": fs, "chi_from_f": chi_from_f(fs)}
(HERE / "01_controls_results.json").write_text(json.dumps(OUT, indent=1))
print(json.dumps({k: (v.get("betti_Z3"), v.get("chi_from_f")) for k, v in OUT.items() if isinstance(v, dict)}))
