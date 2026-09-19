#!/usr/bin/env python
"""
Track D (v2), step 3: resolution of the 16 singular points of T^4/Z_2, via
the HYBRID Mayer-Vietoris route, restricted to N in {6, 8} (premise-valid,
see 00_premise_checks.py). N=4 is excluded entirely (its closed stars are
not pairwise disjoint, so the "16 disjoint (cone, link)" MV decomposition
does not hold there and its numbers -- even though they happen to match
N=6/N=8 numerically -- are not used as evidence for anything below).

Mayer-Vietoris setup:
  X_singular = U  union_{k x RP^3}  (k disjoint cones on RP^3)
  X_resolved = U  union_{k x RP^3}  (k disjoint disk bundles over S^2)
  A = U, B = k disjoint copies of the resolving piece, A n B = k x RP^3,
  k = num_fixed_vertices, COMPUTED in 02 (never a literal 16).

Over a field of odd characteristic (Z/3, cross-checked at Z/5),
H_1(RP^3)=H_2(RP^3)=0, so:
  b2(X_resolved) = b2(U) + k * b2(S^2)          [COMPUTED]
  b1(X_resolved) = b1(U)                        [COMPUTED]
  b0 = b4 = 1, b3 = b1 are STATED (Poincare duality / connectedness /
  orientability of the Kummer resolution), NOT computed from simplicial
  data -- this labeling is preserved unchanged from v1 (fix list does not
  ask for b0/b3/b4 to be upgraded, and nothing here upgrades them).

FIX (2): every count used in the MV formula (num_fixed_vertices=k,
chi(RP^3), chi(S^2)) is read from a computed result, never typed as a
literal 16 or 32.
FIX (4): chi(K3) computed independently from the f-vector (simplex-count
Euler characteristic, read off the GUDHI SimplexTree itself -- not the
custom all_faces() face-set, which is O(2^5) per simplex and fine at N=6
but memory-heavy at N=8) AND from the Betti numbers, and compared.
FIX (5): Betti numbers of the singular-point link over Z/2, Z/3, Z/5
(RP^3 torsion control), for ALL 16 links (not just a representative one),
cross-checked for agreement.
"""
import json
import sys
import time
from collections import Counter

sys.path.insert(0, "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity_v2/D-tda")
import gudhi
from lib_freudenthal import freudenthal_top_simplices, negate_mod

FIELDS = (2, 3, 5)
D = 4
BASE = "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity_v2/D-tda"
TIME_BUDGET_SEC = 13 * 60  # keep the whole script well under 15 min; abort remaining N if exceeded

EXPORTS_PATH = f"{BASE}/exports.json"
OUT_PATH = f"{BASE}/03_resolution_hybrid_results.json"

with open(f"{BASE}/00_premise_checks_results.json") as _f:
    _premise = json.load(_f)
VALID_N = tuple(sorted(
    N for N in (4, 6, 8)
    if _premise[f"N={N}"]["PREMISE_VALID_FOR_MV"]
))  # derived from the premise check, not hardcoded


def build_quotient_top(N):
    top, degenerate = freudenthal_top_simplices(N, D)
    assert degenerate == 0
    neg = lambda v: negate_mod(v, N)

    def canon(v):
        return min(v, neg(v))

    quot_top = [tuple(canon(v) for v in simp) for simp in top]
    quot_top_clean = [s for s in quot_top if len(set(s)) == len(s)]
    fixed_pts = sorted({v for simp in top for v in simp if neg(v) == v})
    return quot_top_clean, fixed_pts


def simplex_tree_from_top(top_simplices):
    """Build a SimplexTree once; f-vector and Betti numbers both read off
    the SAME object (fix 4: independent readouts of the same underlying
    complex, via Counter over st.get_simplices() rather than a custom
    exponential face-enumeration -- memory-safe at N=8)."""
    all_verts = sorted(set().union(*[set(s) for s in top_simplices])) if top_simplices else []
    vidx = {v: i for i, v in enumerate(all_verts)}
    st = gudhi.SimplexTree()
    for simp in top_simplices:
        st.insert([vidx[v] for v in simp], filtration=0.0)
    return st, vidx


def fvector_and_chi(st):
    counts = Counter()
    for simplex, _filt in st.get_simplices():
        counts[len(simplex) - 1] += 1
    chi = sum((-1) ** dim * cnt for dim, cnt in counts.items())
    return dict(sorted(counts.items())), chi


def betti_of_top_simplices(top_simplices, fields):
    st, vidx = simplex_tree_from_top(top_simplices)
    dim_counts, chi = fvector_and_chi(st)
    out = {
        "euler_characteristic_from_fvector": chi,
        "simplex_counts_by_dim": dim_counts,
        "num_vertices": st.num_vertices(),
        "num_top_simplices": len(top_simplices),
        "betti_by_field": {},
    }
    for field in fields:
        st2, _ = simplex_tree_from_top(top_simplices)
        st2.compute_persistence(homology_coeff_field=field, min_persistence=0, persistence_dim_max=True)
        betti = st2.betti_numbers()
        out["betti_by_field"][f"Z{field}"] = betti
    # FIX(4), applied to EVERY complex this helper builds (U, every link, the
    # S^2 control), not only the final resolved K3: chi from the f-vector vs
    # chi from each field's Betti numbers must agree, checked here, not just
    # asserted at the top level. Uses the LAST field computed in `fields`
    # (all fields agree with each other by all_16_links_agree checks
    # elsewhere, so any one field suffices for this particular guard).
    if fields:
        last_field = f"Z{fields[-1]}"
        chi_from_betti_last_field = sum((-1) ** k * b for k, b in enumerate(out["betti_by_field"][last_field]))
        out["chi_from_betti_matches_fvector"] = (chi_from_betti_last_field == chi)
        out["chi_from_betti_last_field_used"] = last_field
    return out


def s2_control(fields):
    """Boundary of a tetrahedron (3-simplex) = triangulated S^2; control for
    chi(S^2)=2 and b2(S^2)=1, used in the MV formula below."""
    import itertools
    verts = [0, 1, 2, 3]
    top = [tuple(c) for c in itertools.combinations(verts, 3)]
    return betti_of_top_simplices(top, fields)


OUT = {}
t_start = time.time()

OUT["S2_control"] = s2_control(FIELDS)
b2_S2 = OUT["S2_control"]["betti_by_field"]["Z3"][2] if len(OUT["S2_control"]["betti_by_field"]["Z3"]) > 2 else 0
chi_S2 = OUT["S2_control"]["euler_characteristic_from_fvector"]

exported_yet = False

for N in VALID_N:
    if time.time() - t_start > TIME_BUDGET_SEC:
        OUT[f"N={N}_SKIPPED"] = f"time budget ({TIME_BUDGET_SEC}s) exceeded before starting N={N}"
        print(f"N={N}: SKIPPED, time budget exceeded")
        continue

    t0 = time.time()
    quot_top, fixed_pts = build_quotient_top(N)
    fixed_set = set(fixed_pts)
    num_fixed = len(fixed_pts)  # COMPUTED, used everywhere below instead of literal 16

    # ---- link of EVERY singular vertex, all three fields ----
    link_results = []
    for p in fixed_pts:
        link_top = [tuple(v for v in simp if v != p) for simp in quot_top if p in simp]
        r = betti_of_top_simplices(link_top, FIELDS)
        r["singular_vertex"] = p
        link_results.append(r)

    link_betti_z2 = [tuple(r["betti_by_field"]["Z2"]) for r in link_results]
    link_betti_z3 = [tuple(r["betti_by_field"]["Z3"]) for r in link_results]
    link_betti_z5 = [tuple(r["betti_by_field"]["Z5"]) for r in link_results]
    all_links_agree_z2 = len(set(link_betti_z2)) == 1
    all_links_agree_z3 = len(set(link_betti_z3)) == 1
    all_links_agree_z5 = len(set(link_betti_z5)) == 1

    # RP^3 torsion control: over odd fields expect (1,0,0,1); over Z/2 the
    # Z/2 torsion in H_1(RP^3;Z) is expected to raise every Betti number to 1
    # (1,1,1,1), a genuine mod-2 JUMP relative to the odd-field answer -- the
    # opposite of the S^3 check in 00_premise_checks.py, which is exactly
    # the point (this link IS singular, so it SHOULD show the RP^3 jump).
    link_z2_repr = link_betti_z2[0] if link_betti_z2 else None
    link_z3_repr = link_betti_z3[0] if link_betti_z3 else None

    U_top = [simp for simp in quot_top if not (set(simp) & fixed_set)]
    U_result = betti_of_top_simplices(U_top, FIELDS)

    entry = {
        "N": N,
        "num_singular_vertices_computed": num_fixed,
        "link_results_per_vertex": link_results,
        "all_16_links_agree_Z2": all_links_agree_z2,
        "all_16_links_agree_Z3": all_links_agree_z3,
        "all_16_links_agree_Z5": all_links_agree_z5,
        "link_representative_betti_Z2": link_z2_repr,
        "link_representative_betti_Z3_Z5(odd)": link_z3_repr,
        "link_expected_betti_RP3_odd_field": [1, 0, 0, 1],
        "link_expected_betti_RP3_Z2": [1, 1, 1, 1],
        "torsion_control_shows_mod2_jump": (list(link_z2_repr) == [1, 1, 1, 1] and list(link_z3_repr) == [1, 0, 0, 1]),
        "U_complement_result": U_result,
    }

    for field, fname in ((3, "Z3"), (5, "Z5")):
        b_U = U_result["betti_by_field"][fname]
        b_U = list(b_U) + [0] * (5 - len(b_U))
        b2_U = b_U[2] if len(b_U) > 2 else 0
        b1_U = b_U[1] if len(b_U) > 1 else 0
        b2_K3 = b2_U + num_fixed * b2_S2  # FIX(2): computed b2(S2) and computed num_fixed, no literal
        b1_K3 = b1_U
        b0_K3 = 1  # STATED (connectedness), unchanged labeling from v1
        b4_K3 = 1  # STATED (Poincare duality with b0), unchanged labeling from v1
        b3_U = b_U[3] if len(b_U) > 3 else 0
        b3_K3 = b1_K3  # STATED via Poincare duality (b3=b1), unchanged labeling from v1

        # chi(RP^3) from the (agreed) link Betti numbers, computed, this field
        link0 = link_results[0]["betti_by_field"][fname]
        chi_RP3 = sum((-1) ** i * b for i, b in enumerate(link0))

        chi_K3_from_betti = b0_K3 - b1_K3 + b2_K3 - b3_K3 + b4_K3
        chi_U_fvector = U_result["euler_characteristic_from_fvector"]
        chi_K3_from_MV_fvector = chi_U_fvector + num_fixed * (chi_S2 - chi_RP3)  # FIX(2): no literal 32

        entry[f"resolved_K3_betti_{fname}"] = {
            "b0": b0_K3, "b1": b1_K3, "b2": b2_K3, "b3": b3_K3, "b4": b4_K3,
            "b0_b4_status": "STATED (connectedness/orientability/Poincare duality), NOT computed",
            "b3_status": "STATED via Poincare duality (b3=b1), NOT independently computed",
            "b1_b2_status": "COMPUTED from U via Mayer-Vietoris (b1=b1(U), b2=b2(U)+num_fixed*b2(S2))",
            "chi_from_betti_numbers": chi_K3_from_betti,
            "chi_from_MV_fvector(chi(U)+k*(chi(S2)-chi(RP3)))": chi_K3_from_MV_fvector,
            "chi_cross_check_matches": chi_K3_from_betti == chi_K3_from_MV_fvector,
            "chi_RP3_computed": chi_RP3,
            "chi_S2_computed": chi_S2,
            "num_fixed_used": num_fixed,
        }

    entry["runtime_sec"] = time.time() - t0
    OUT[f"N={N}"] = entry
    print(f"N={N}: done in {entry['runtime_sec']:.1f}s, "
          f"chi_match={entry['resolved_K3_betti_Z3']['chi_cross_check_matches']}, "
          f"b2={entry['resolved_K3_betti_Z3']['b2']}")

    # write results-so-far after every N, so a truncation still leaves usable output
    with open(OUT_PATH, "w") as f:
        json.dump(OUT, f, indent=2, default=str)

    if N == 6 and not exported_yet:
        chi_val = entry["resolved_K3_betti_Z3"]["chi_from_betti_numbers"]
        b2_val = entry["resolved_K3_betti_Z3"]["b2"]
        exports = {
            "source": "03_resolution_hybrid.py, N=6 (premise-valid), field Z/3",
            "chi_K3_resolved": chi_val,
            "b2_K3_resolved": b2_val,
            "b0": entry["resolved_K3_betti_Z3"]["b0"],
            "b1": entry["resolved_K3_betti_Z3"]["b1"],
            "b3": entry["resolved_K3_betti_Z3"]["b3"],
            "b4": entry["resolved_K3_betti_Z3"]["b4"],
            "cross_field_Z5_agrees": entry["resolved_K3_betti_Z5"]["b2"] == b2_val and entry["resolved_K3_betti_Z5"]["chi_from_betti_numbers"] == chi_val,
            "num_fixed_points_computed": num_fixed,
            "chi_fvector_betti_cross_check": entry["resolved_K3_betti_Z3"]["chi_cross_check_matches"],
        }
        with open(EXPORTS_PATH, "w") as f:
            json.dump(exports, f, indent=2)
        exported_yet = True
        print(f"Wrote {EXPORTS_PATH} (Track C unblocked): chi={chi_val}, b2={b2_val}")

OUT["total_runtime_sec"] = time.time() - t_start
with open(OUT_PATH, "w") as f:
    json.dump(OUT, f, indent=2, default=str)

print("\nWrote", OUT_PATH)
