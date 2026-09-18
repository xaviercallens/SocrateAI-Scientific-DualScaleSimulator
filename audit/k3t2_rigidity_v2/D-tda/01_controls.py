#!/usr/bin/env python
"""
Track D (v2), step 1: known-answer controls, restricted to the valid grid
sizes established in 00_premise_checks.py (N=6, N=8; N=4 fails the
closed-star-disjointness premise and is excluded from every downstream
computation, though its raw torus control below is still reported for
completeness -- it is a control on the UNQUOTIENTED T^4, unaffected by
the quotient premise).

  (a) T^2 via gudhi.PeriodicCubicalComplex -> expect Betti (1,2,1).
      Needed later as an input to the Kunneth computation for K3 x T^2.
  (b) T^4 via gudhi.PeriodicCubicalComplex -> expect Betti (1,4,6,4,1).
  (c) T^4 via the Freudenthal/Kuhn SIMPLICIAL triangulation of Z_N^4, at
      N=4 (reference only), N=6, N=8 -- must reproduce (1,4,6,4,1).

Homology coefficient field: Z/3 (odd prime), consistent with the rest of
this track.
"""
import json
import sys
import time

sys.path.insert(0, "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity_v2/D-tda")
import gudhi
from lib_freudenthal import (
    freudenthal_top_simplices,
    all_faces,
    euler_characteristic_from_faces,
    make_vertex_encoder,
)

FIELD = 3
OUT = {}


def cubical_torus_betti(d, cells_per_dim=3, field=FIELD):
    n_cells = cells_per_dim ** d
    top_cells = [0.0] * n_cells
    pcc = gudhi.PeriodicCubicalComplex(
        dimensions=[cells_per_dim] * d,
        top_dimensional_cells=top_cells,
        periodic_dimensions=[True] * d,
    )
    pcc.compute_persistence(homology_coeff_field=field, min_persistence=0)
    betti = pcc.betti_numbers()
    return betti


t0 = time.time()
b_t2 = cubical_torus_betti(2, cells_per_dim=3)
b_t4 = cubical_torus_betti(4, cells_per_dim=3)
t_cubical = time.time() - t0

OUT["cubical"] = {
    "method": "gudhi.PeriodicCubicalComplex, all dims periodic, 3 cells/dim, filtration=0 everywhere",
    "homology_coeff_field": FIELD,
    "T2": {"betti": b_t2, "expected": [1, 2, 1], "expected_source": "b_k(T^n)=C(n,k), n=2 (standard torus cohomology)"},
    "T4": {"betti": b_t4, "expected": [1, 4, 6, 4, 1], "expected_source": "b_k(T^n)=C(n,k), n=4 (standard torus cohomology)"},
    "runtime_sec": t_cubical,
}

results_simplicial = {}
for N in (4, 6, 8):
    t0 = time.time()
    top, degenerate = freudenthal_top_simplices(N, 4)
    encode, decode = make_vertex_encoder(N, 4)
    st = gudhi.SimplexTree()
    for simp in top:
        st.insert([encode(v) for v in simp], filtration=0.0)
    n_top_expected = (N ** 4) * 24
    st.compute_persistence(homology_coeff_field=FIELD, min_persistence=0, persistence_dim_max=True)
    betti = st.betti_numbers()
    faces = all_faces(top)
    chi, dim_counts = euler_characteristic_from_faces(faces)
    chi_from_betti = sum((-1) ** k * b for k, b in enumerate(betti))
    results_simplicial[str(N)] = {
        "N": N,
        "num_top_4simplices_built": len(top),
        "num_top_4simplices_expected_formula": n_top_expected,
        "degenerate_simplices_dropped": degenerate,
        "num_vertices": st.num_vertices(),
        "simplex_tree_num_simplices": st.num_simplices(),
        "betti": betti,
        "expected_betti": [1, 4, 6, 4, 1],
        "expected_source": "b_k(T^4)=C(4,k), standard torus cohomology",
        "euler_characteristic_from_simplex_counts": chi,
        "euler_characteristic_from_betti_numbers": chi_from_betti,
        "simplex_counts_by_dim": dim_counts,
        "runtime_sec": time.time() - t0,
        "matches_control": betti == [1, 4, 6, 4, 1] and chi == 0 and chi_from_betti == 0,
        "note": "reference only, quotient premise not applicable to raw T^4" if N == 4 else "used downstream (valid N)",
    }

OUT["freudenthal_simplicial_T4"] = results_simplicial

OUT_PATH = "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity_v2/D-tda/01_controls_results.json"
with open(OUT_PATH, "w") as f:
    json.dump(OUT, f, indent=2)

print(json.dumps(OUT, indent=2))
print("\nWrote", OUT_PATH)
