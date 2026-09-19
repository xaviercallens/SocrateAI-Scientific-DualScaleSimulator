#!/usr/bin/env python
"""
Track D (v2), step 2: T^4/Z_2 by an actual quotient simplicial complex,
restricted to N in {6, 8} (the grid sizes that passed the
closed-star-disjointness premise check in 00_premise_checks.py; N=4 is
EXCLUDED here, not merely footnoted, per that check).

For each valid N:
  (a) Build Freudenthal top-4-simplices of Z_N^4 (validated against
      (1,4,6,4,1) in 01_controls.py at this same N).
  (b) Verify computationally that negation g: v -> -v mod N maps the
      simplex SET to itself.
  (c) Count fixed vertices of g (computed, not assumed -- 16 is expected
      but the number used everywhere downstream is len(fixed_pts), never
      the literal 16).
  (d) Verify g acts freely on the top 4-simplices (no simplex contains
      both v and g(v), v != g(v)).
  (e) Build the quotient complex via canonical-representative collapse.
  (f) Compute Betti numbers of the quotient over Z/3 and Z/5.
  (g) Independent Euler-characteristic check from raw simplex counts.
  (h) NEGATIVE CONTROL: fixed-point-free translation h: v -> v + shift,
      order 2. FIX (v1 defect): the translation control must pass the
      SAME exact-halving check the negation quotient is held to, across
      ALL dimensions INCLUDING dimension 0 (translation has zero fixed
      points, so unlike the negation quotient -- whose vertex count does
      NOT halve because of the 16 fixed vertices -- the translation
      quotient's vertex count MUST also exactly halve; if it does not,
      the control is dropped and flagged rather than reported).
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
    negate_mod,
    translate_mod,
)

FIELDS = (3, 5)
D = 4
BASE = "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity_v2/D-tda"

with open(f"{BASE}/00_premise_checks_results.json") as _f:
    _premise = json.load(_f)
VALID_N = tuple(sorted(
    N for N in (4, 6, 8)
    if _premise[f"N={N}"]["PREMISE_VALID_FOR_MV"]
))  # derived from the premise check, not hardcoded

with open(f"{BASE}/01_controls_results.json") as _f:
    _controls = json.load(_f)


def check_action_invariance(top, action):
    simplex_set = set(frozenset(s) for s in top)
    for simp in top:
        mapped = frozenset(action(v) for v in simp)
        if mapped not in simplex_set:
            return False
    return True


def action_acts_freely_on_top_simplices(top, action):
    for simp in top:
        vs = set(simp)
        for v in simp:
            gv = action(v)
            if gv != v and gv in vs:
                return False
    return True


def quotient_betti(top, action, fields):
    def canon(v):
        gv = action(v)
        return min(v, gv)

    quot_top = [tuple(canon(v) for v in simp) for simp in top]
    quot_top_clean = [s for s in quot_top if len(set(s)) == len(s)]
    n_degenerate = len(quot_top) - len(quot_top_clean)

    faces = all_faces(quot_top_clean)
    chi, dim_counts = euler_characteristic_from_faces(faces)

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
        results[f"Z{field}"] = {"betti": betti, "euler_characteristic_from_betti": chi_from_betti}

    return {
        "num_top_simplices_input": len(top),
        "num_degenerate_after_quotient": n_degenerate,
        "num_vertices_quotient": len(all_verts),
        "euler_characteristic_from_simplex_counts": chi,
        "simplex_counts_by_dim": dim_counts,
        "betti_by_field": results,
    }


_excluded_N = {
    str(N): "fails PREMISE_VALID_FOR_MV (star disjointness and/or nonsingular-link-S3 check), see 00_premise_checks_results.json"
    for N in (4, 6, 8) if N not in VALID_N
}
OUT = {"D": D, "valid_N_used": list(VALID_N), "excluded_N": _excluded_N}

for N in VALID_N:
    t0 = time.time()
    top, degenerate = freudenthal_top_simplices(N, D)
    assert degenerate == 0

    neg = lambda v, N=N: negate_mod(v, N)
    is_inv = check_action_invariance(top, neg)
    fixed_pts = sorted({v for simp in top for v in simp if neg(v) == v})
    acts_freely = action_acts_freely_on_top_simplices(top, neg)

    entry = {
        "N": N,
        "negation_action_invariant_on_top_simplices": is_inv,
        "num_fixed_vertices_computed": len(fixed_pts),
        "fixed_vertices": fixed_pts,
        "action_free_on_top_4simplices": acts_freely,
        "runtime_setup_sec": time.time() - t0,
    }

    assert is_inv and acts_freely, f"N={N}: negation quotient construction invalid, unexpected"

    t1 = time.time()
    q = quotient_betti(top, neg, FIELDS)
    q["runtime_sec"] = time.time() - t1
    entry["quotient_T4_mod_Z2"] = q

    chi_T4 = _controls["freudenthal_simplicial_T4"][str(N)]["euler_characteristic_from_simplex_counts"]
    F = len(fixed_pts)
    chi_predicted = (chi_T4 - F) // 2 + F
    entry["chi_predicted_from_orbit_counting"] = chi_predicted
    entry["chi_predicted_formula"] = "(chi(T4) - F)/2 + F, F = COMPUTED fixed vertex count"
    entry["chi_predicted_matches_simplex_count"] = (chi_predicted == q["euler_characteristic_from_simplex_counts"])

    t4_faces = all_faces(top)
    _, t4_dim_counts_full = euler_characteristic_from_faces(t4_faces)
    halving_check = {}
    for dim, cnt in q["simplex_counts_by_dim"].items():
        dim = int(dim)  # dict keys are ints in-memory but may already be str; normalize
        t4_cnt = t4_dim_counts_full.get(dim, 0)
        if dim == 0:
            halving_check[str(dim)] = {"quotient": cnt, "T4": t4_cnt,
                                        "note": f"vertices: not expected to halve ({F} fixed vertices don't pair up)"}
        else:
            halving_check[str(dim)] = {"quotient": cnt, "T4": t4_cnt, "exactly_half": (2 * cnt == t4_cnt)}
    entry["quotient_T4_mod_Z2"]["exact_halving_check_dims_1_to_4"] = halving_check
    entry["quotient_T4_mod_Z2"]["all_dims_1_to_4_exactly_half"] = all(
        v.get("exactly_half") for k, v in halving_check.items() if int(k) != 0
    )

    # --- negative control: fixed-point-free translation, order 2 ---
    shift = (N // 2, 0, 0, 0)
    trans = lambda v, N=N, shift=shift: translate_mod(v, N, shift)
    is_inv_t = check_action_invariance(top, trans)
    fixed_pts_t = sorted({v for simp in top for v in simp if trans(v) == v})
    acts_freely_t = action_acts_freely_on_top_simplices(top, trans)
    neg_ctrl = {
        "shift": shift,
        "translation_invariant_on_top_simplices": is_inv_t,
        "num_fixed_vertices_computed": len(fixed_pts_t),
        "expected_fixed_vertices": 0,
        "action_free_on_top_4simplices": acts_freely_t,
    }
    include_control = is_inv_t and acts_freely_t and len(fixed_pts_t) == 0
    if include_control:
        t2 = time.time()
        qt = quotient_betti(top, trans, FIELDS)
        qt["runtime_sec"] = time.time() - t2

        # FIX (3): halving check for the translation control, ALL dims
        # including dim 0 (translation has zero fixed points, so unlike the
        # negation quotient the vertex count must also exactly halve).
        t_halving = {}
        for dim_str, cnt in qt["simplex_counts_by_dim"].items():
            dim = int(dim_str)
            t4_cnt = t4_dim_counts_full.get(dim, 0)
            t_halving[dim_str] = {"quotient": cnt, "T4": t4_cnt, "exactly_half": (2 * cnt == t4_cnt)}
        qt["exact_halving_check_all_dims_incl_0"] = t_halving
        all_half = all(v["exactly_half"] for v in t_halving.values())
        qt["all_dims_incl_0_exactly_half"] = all_half

        neg_ctrl["quotient_result"] = qt
        neg_ctrl["passes_same_halving_check_as_negation_quotient"] = all_half
        neg_ctrl["DROPPED"] = not all_half
        if all_half:
            neg_ctrl["expected_betti"] = [1, 4, 6, 4, 1]
            neg_ctrl["expected_source"] = "free Z2 quotient of T^4 by a translation is again a smaller-period 4-torus"
            for field in FIELDS:
                b = qt["betti_by_field"][f"Z{field}"]["betti"]
                neg_ctrl[f"matches_expected_Z{field}"] = (b == [1, 4, 6, 4, 1])
        else:
            neg_ctrl["NOTE"] = "translation control FAILED the halving check at this N; DROPPED per task fix (3), not used as a control"
    else:
        neg_ctrl["quotient_result"] = None
        neg_ctrl["DROPPED"] = True
        neg_ctrl["NOTE"] = "translation action was not invariant/free/fixed-point-free as required; control dropped"

    entry["negative_control_free_translation_quotient"] = neg_ctrl

    OUT[f"N={N}"] = entry

OUT_PATH = "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity_v2/D-tda/02_invariance_and_quotient_results.json"
with open(OUT_PATH, "w") as f:
    json.dump(OUT, f, indent=2, default=str)

print(json.dumps({k: v for k, v in OUT.items() if not isinstance(v, dict) or "N" not in v}, indent=2, default=str))
print("\nWrote", OUT_PATH)
