#!/usr/bin/env python
"""
Track D, step (1): known-answer controls.

  (a) T^2 via gudhi.PeriodicCubicalComplex  -> expect Betti (1,2,1)
  (b) T^4 via gudhi.PeriodicCubicalComplex  -> expect Betti (1,4,6,4,1)
  (c) T^4 via the Freudenthal/Kuhn SIMPLICIAL triangulation of Z_N^4,
      at N=4 and N=6 (the load-bearing control: this is the exact same
      construction the T^4/Z_2 quotient in step (2) will be built from,
      so it must reproduce (1,4,6,4,1) at the N we actually use).

All numbers below are computed by GUDHI from the stated construction;
none are hard-coded as targets used anywhere in a comparison other than
the human-readable "expected" field, which is filled from the standard
formula for Betti numbers of a torus, b_k(T^n) = C(n,k), stated as the
source.

Homology coefficient field: Z/3 (odd prime, avoids H^*(T^n;Z) 2-torsion,
which is trivial anyway for a torus, but Z/3 is used throughout this
track for consistency with the quotient/resolution computations where
2-torsion could appear).
"""
import json
import sys
import time

sys.path.insert(0, "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity/tda-gudhi")
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
    """PeriodicCubicalComplex on a d-torus: a grid of `cells_per_dim` cells
    per dimension, all dimensions periodic. Each top-dimensional cell gets
    filtration value 0 so the complex is exactly the d-torus CW structure."""
    n_cells = cells_per_dim ** d
    top_cells = [0.0] * n_cells
    pcc = gudhi.PeriodicCubicalComplex(
        dimensions=[cells_per_dim] * d,
        top_dimensional_cells=top_cells,
        periodic_dimensions=[True] * d,
    )
    pcc.compute_persistence(homology_coeff_field=field, min_persistence=0)
    betti = pcc.persistent_betti_numbers(0, 0) if hasattr(pcc, "persistent_betti_numbers") else None
    # PeriodicCubicalComplex exposes betti_numbers() after compute_persistence
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

# --- (c) Freudenthal simplicial T^4 control, at the N we will quotient ---
results_simplicial = {}
for N in (4, 6):
    t0 = time.time()
    top, degenerate = freudenthal_top_simplices(N, 4)
    encode, decode = make_vertex_encoder(N, 4)
    st = gudhi.SimplexTree()
    for simp in top:
        st.insert([encode(v) for v in simp], filtration=0.0)
    n_top_expected = (N ** 4) * 24  # N^d base cubes * d! permutations, d=4
    st.compute_persistence(homology_coeff_field=FIELD, min_persistence=0, persistence_dim_max=True)
    betti = st.betti_numbers()
    # independent Euler-characteristic check straight from simplex counts
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
    }

OUT["freudenthal_simplicial_T4"] = results_simplicial

OUT_PATH = "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity/tda-gudhi/01_controls_results.json"
with open(OUT_PATH, "w") as f:
    json.dump(OUT, f, indent=2)

print(json.dumps(OUT, indent=2))
print("\nWrote", OUT_PATH)
