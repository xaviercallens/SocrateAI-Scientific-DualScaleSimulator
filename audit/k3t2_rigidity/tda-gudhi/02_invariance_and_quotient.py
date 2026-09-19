#!/usr/bin/env python
"""
Track D, step (2): T^4/Z_2 by an ACTUAL quotient simplicial complex.

For each N in {4, 6}:
  (a) Build the Freudenthal top-4-simplices of Z_N^4 (same construction as
      01_controls.py, already validated there against (1,4,6,4,1)).
  (b) Verify computationally (not assumed) that the involution
          g: v -> -v mod N   (component-wise negation)
      maps the simplex SET to itself: for every top simplex, negate every
      vertex and check the resulting simplex is in the same set.
  (c) Count fixed vertices of g: v with -v = v mod N, i.e. 2v = 0 mod N.
      For N=4 this is v_i in {0,2} for each of 4 coords -> 16 points, for
      N=6 it is v_i in {0,3} -> also 16 points. Report the count computed,
      not assumed.
  (d) Verify no top simplex contains two vertices of a single g-orbit
      {v, g(v)} with v != g(v) (i.e. g acts freely on the SET of top
      4-simplices) -- this is required for the naive quotient (identify
      each simplex with its image) to be a valid simplicial complex
      instead of a self-folded one.
  (e) Build the quotient complex: canonical representative of each vertex
      orbit is min(v, g(v)) in lexicographic order; map every top simplex
      through this canonicalization and insert into a fresh SimplexTree
      (duplicates from {simplex, g(simplex)} collapse automatically since
      GUDHI simplicial trees are sets of simplices).
  (f) Compute Betti numbers of the quotient over Z/3 AND Z/5 (two odd
      primes, cross-check).
  (g) Independent Euler characteristic check from raw simplex counts of
      the quotient SimplexTree (not from the Betti numbers).
  (h) NEGATIVE CONTROL: repeat (e)-(g) with a fixed-point-free action,
      the translation h: v -> v + (2,0,0,0) mod N (order 2, since N is a
      multiple of 2... N=4: translating by 2 twice gives 4=0 mod 4, order
      2; N=6: translating by 2 twice gives 4 != 0 mod 6 -> order 3, so for
      N=6 we instead use shift (3,0,0,0), order 2). This has ZERO fixed
      points, so it demonstrates the pipeline is sensitive to the actual
      group action (a genuine free quotient of T^4 by Z_2, still a
      4-torus, must reproduce Betti (1,4,6,4,1), NOT the singular
      T^4/Z_2 numbers) rather than always outputting some default answer.
"""
import json
import sys
import time
from collections import Counter

sys.path.insert(0, "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity/tda-gudhi")
import gudhi
from lib_freudenthal import (
    freudenthal_top_simplices,
    all_faces,
    euler_characteristic_from_faces,
    negate_mod,
    translate_mod,
)

FIELDS = (3, 5)
D = 4


def check_action_invariance(top, action, N):
    """Return (is_invariant, num_checked)."""
    simplex_set = set(frozenset(s) for s in top)
    for simp in top:
        mapped = frozenset(action(v) for v in simp)
        if mapped not in simplex_set:
            return False, len(top)
    return True, len(top)


def action_acts_freely_on_top_simplices(top, action):
    """No top simplex may contain both v and action(v) for v != action(v)."""
    for simp in top:
        vs = set(simp)
        for v in simp:
            gv = action(v)
            if gv != v and gv in vs:
                return False
    return True


def quotient_betti(top, action, N, fields):
    """Build the quotient complex by canonical-representative collapse and
    compute Betti numbers over each field in `fields`. Also returns raw
    simplex counts for an independent Euler-characteristic check."""
    def canon(v):
        gv = action(v)
        return min(v, gv)

    quot_top = [tuple(canon(v) for v in simp) for simp in top]
    # drop any simplex that degenerated (repeated vertex after canonicalization)
    quot_top_clean = [s for s in quot_top if len(set(s)) == len(s)]
    n_degenerate_after_quotient = len(quot_top) - len(quot_top_clean)

    faces = all_faces(quot_top_clean)
    chi, dim_counts = euler_characteristic_from_faces(faces)

    # vertex encoding for SimplexTree
    all_verts = sorted(set().union(*[set(s) for s in quot_top_clean])) if quot_top_clean else []
    vidx = {v: i for i, v in enumerate(all_verts)}

    results = {}
    for field in fields:
        st = gudhi.SimplexTree()
        for simp in quot_top_clean:
            st.insert([vidx[v] for v in simp], filtration=0.0)
        st.compute_persistence(homology_coeff_field=field, min_persistence=0, persistence_dim_max=True)
        betti = st.betti_numbers()
        chi_from_betti = sum((-1) ** k * b for k, b in enumerate(betti))
        results[f"Z{field}"] = {
            "betti": betti,
            "euler_characteristic_from_betti": chi_from_betti,
        }

    return {
        "num_top_simplices_input": len(top),
        "num_top_simplices_quotient_raw": len(quot_top),
        "num_top_simplices_quotient_after_dedup": len(set(frozenset(s) for s in quot_top_clean if len(set(s)) == len(s))),
        "num_degenerate_after_quotient": n_degenerate_after_quotient,
        "num_vertices_quotient": len(all_verts),
        "euler_characteristic_from_simplex_counts": chi,
        "simplex_counts_by_dim": dim_counts,
        "betti_by_field": results,
    }


OUT = {"D": D}

for N in (4, 6):
    t0 = time.time()
    top, degenerate = freudenthal_top_simplices(N, D)
    assert degenerate == 0, f"N={N}: base Freudenthal triangulation had degenerate simplices"

    neg = lambda v, N=N: negate_mod(v, N)
    is_inv, n_checked = check_action_invariance(top, neg, N)

    fixed_pts = [v for v in {vv for simp in top for vv in simp} if neg(v) == v]
    acts_freely = action_acts_freely_on_top_simplices(top, neg)

    entry = {
        "N": N,
        "negation_action_invariant_on_top_simplices": is_inv,
        "num_top_simplices_checked": n_checked,
        "num_fixed_vertices_computed": len(fixed_pts),
        "fixed_vertices": sorted(fixed_pts),
        "expected_fixed_vertex_count": 16,
        "expected_source": "16 two-torsion points of (Z/N)^4 under x->-x (each coord in {0,N/2}), 2^4=16",
        "action_free_on_top_4simplices": acts_freely,
        "runtime_setup_sec": time.time() - t0,
    }

    if is_inv:
        t1 = time.time()
        q = quotient_betti(top, neg, N, FIELDS)
        q["runtime_sec"] = time.time() - t1
        entry["quotient_T4_mod_Z2"] = q

        # independent chi check: chi(T^4/Z2) should equal (chi(T^4) - F)/2 + F
        # for a quotient by an involution with F fixed points (each fixed
        # point contributes chi=1 alone, each free orbit of 2 points
        # contributes chi=1 to the quotient instead of 2).
        chi_T4 = 0  # computed in 01_controls.py
        F = len(fixed_pts)
        chi_predicted = (chi_T4 - F) // 2 + F
        entry["chi_predicted_from_orbit_counting"] = chi_predicted
        entry["chi_predicted_formula"] = "(chi(T4) - F)/2 + F, F = fixed vertex count, chi(T4)=0 from 01_controls.py"
        entry["chi_predicted_matches_simplex_count"] = (chi_predicted == q["euler_characteristic_from_simplex_counts"])

        # Additional validity check on the canonical-collapse construction:
        # since the action is free on vertices >0-dim... actually on
        # simplices of dim>=1 it is exactly free (no k-simplex, k>=1, can be
        # setwise fixed by negation without a fixed vertex among its
        # vertices being present, given num_fixed=16 and free action on top
        # simplices already verified), so the quotient's simplex COUNT in
        # every dimension >=1 must be EXACTLY HALF of T^4's; only the vertex
        # count deviates (fixed vertices don't pair up). This is a stronger,
        # purely-combinatorial check that no unwanted simplex collisions
        # occurred anywhere in the construction (not just in the top cells).
        t4_dim_counts = {"0": N ** 4, "1": None, "2": None, "3": None, "4": None}
        # recompute T^4's own simplex counts by dimension directly (cheap, same top list)
        t4_faces = all_faces(top)
        _, t4_dim_counts_full = euler_characteristic_from_faces(t4_faces)
        halving_check = {}
        for dim_str, cnt in q["simplex_counts_by_dim"].items():
            dim = int(dim_str)
            t4_cnt = t4_dim_counts_full.get(dim, 0)
            if dim == 0:
                halving_check[dim_str] = {"quotient": cnt, "T4": t4_cnt, "note": "vertices: not expected to halve (16 fixed vertices don't pair up)"}
            else:
                halving_check[dim_str] = {"quotient": cnt, "T4": t4_cnt, "exactly_half": (2 * cnt == t4_cnt)}
        entry["quotient_T4_mod_Z2"]["exact_halving_check_dims_1_to_4"] = halving_check
        entry["quotient_T4_mod_Z2"]["all_dims_1_to_4_exactly_half"] = all(
            v.get("exactly_half") for k, v in halving_check.items() if k != 0
        )
    else:
        entry["quotient_T4_mod_Z2"] = None
        entry["NOTE"] = "negation not invariant at this N; quotient skipped, see could_not_do"

    # --- negative control: fixed-point-free translation action ---
    shift = (2, 0, 0, 0) if N == 4 else (3, 0, 0, 0)
    trans = lambda v, N=N, shift=shift: translate_mod(v, N, shift)
    is_inv_t, _ = check_action_invariance(top, trans, N)
    fixed_pts_t = [v for v in {vv for simp in top for vv in simp} if trans(v) == v]
    acts_freely_t = action_acts_freely_on_top_simplices(top, trans)
    neg_ctrl = {
        "shift": shift,
        "translation_invariant_on_top_simplices": is_inv_t,
        "num_fixed_vertices": len(fixed_pts_t),
        "expected_fixed_vertices": 0,
        "action_free_on_top_4simplices": acts_freely_t,
    }
    if is_inv_t:
        t2 = time.time()
        qt = quotient_betti(top, trans, N, FIELDS)
        qt["runtime_sec"] = time.time() - t2
        neg_ctrl["quotient_result"] = qt
        neg_ctrl["expected_betti"] = [1, 4, 6, 4, 1]
        neg_ctrl["expected_source"] = "free Z2 quotient of T^4 by a translation is again a 4-torus (a smaller-period one), so same Betti numbers as T^4"
        for field in FIELDS:
            b = qt["betti_by_field"][f"Z{field}"]["betti"]
            neg_ctrl[f"matches_expected_Z{field}"] = (b == [1, 4, 6, 4, 1])
    else:
        neg_ctrl["quotient_result"] = None

    entry["negative_control_free_translation_quotient"] = neg_ctrl

    OUT[f"N={N}"] = entry

OUT_PATH = "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity/tda-gudhi/02_invariance_and_quotient_results.json"
with open(OUT_PATH, "w") as f:
    json.dump(OUT, f, indent=2, default=str)

print(json.dumps(OUT, indent=2, default=str))
print("\nWrote", OUT_PATH)
