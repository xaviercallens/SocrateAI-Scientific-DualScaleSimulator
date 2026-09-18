#!/usr/bin/env python
"""
Track D, step (3): resolution of the 16 singular points of T^4/Z_2, via the
HYBRID route (ii) from the task: compute with GUDHI the homology of the
link of a singular vertex (must be RP^3) and of the complement U (quotient
complex minus the open stars of the 16 singular vertices), then finish
with an explicit Mayer-Vietoris calculation. Labeled hybrid throughout.

Runs at N=4 (primary; already validated as an invariant, non-degenerate
Freudenthal quotient in 02_invariance_and_quotient.py) and cross-checks
key numbers at N=6.

Mayer-Vietoris setup (all stated explicitly, degree by degree):
  X_singular = U  union_{16 x RP^3}  (16 disjoint cones on RP^3)
  X_resolved = U  union_{16 x RP^3}  (16 disjoint disk bundles over S^2)
  A = U, B = 16 disjoint copies of (resolving piece), A n B = 16 x RP^3.

Over a field (Z/3 here; Z/5 cross-check), H_1(RP^3)=H_2(RP^3)=0, so:
  degree 2: 0 -> H2(U) + H2(16 x S^2) -> H2(X_resolved) -> H1(16xRP^3)=0
            => b2(X_resolved) = b2(U) + 16 * b2(S^2) = b2(U) + 16   [COMPUTED]
  degree 1: H1(16xRP^3)=0 -> H1(U)+H1(16xS^2) -> H1(X_resolved) -> H0(16xRP^3) -> H0(U)+H0(16xS^2)
            the last map is injective (U connected, each S^2-bundle meets a
            distinct RP^3 component, both maps injective on H0 of a
            connected space) so the connecting map H1(X_resolved)->H0(...)
            is zero, giving b1(X_resolved) = b1(U) + 16*b1(S^2) = b1(U)  [COMPUTED]
  b0 = b4 = 1 and b3 = b1 are STATED (not computed here): b0=1 needs
     connectedness of the resolved space (plausible from U connected +
     gluing, not independently verified by a GUDHI run); b4=1 and b3=b1
     invoke closedness + orientability + Poincare duality, standard facts
     about the Kummer resolution, NOT derived from our simplicial data.
     Tier: these four numbers are LABELED, everything else in this file
     is COMPUTED (tier B).
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
    negate_mod,
)

FIELDS = (3, 5)
D = 4


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


def betti_of_top_simplices(top_simplices, fields):
    faces = all_faces(top_simplices)
    chi, dim_counts = euler_characteristic_from_faces(faces)
    all_verts = sorted(set().union(*[set(s) for s in top_simplices])) if top_simplices else []
    vidx = {v: i for i, v in enumerate(all_verts)}
    out = {"euler_characteristic_from_simplex_counts": chi, "simplex_counts_by_dim": dim_counts,
           "num_vertices": len(all_verts), "num_top_simplices": len(top_simplices), "betti_by_field": {}}
    for field in fields:
        st = gudhi.SimplexTree()
        for simp in top_simplices:
            st.insert([vidx[v] for v in simp], filtration=0.0)
        st.compute_persistence(homology_coeff_field=field, min_persistence=0, persistence_dim_max=True)
        betti = st.betti_numbers()
        # pad to dim of complex if trailing zeros were dropped
        out["betti_by_field"][f"Z{field}"] = betti
    return out


def s2_control(fields):
    """Boundary of a 4-simplex = triangulated S^2... wait boundary of a
    TETRAHEDRON (3-simplex, 4 vertices) is S^2. Compute its Betti numbers
    with GUDHI as a control for chi(S^2)=2, used later in the rigidity
    scan and Mayer-Vietoris arithmetic."""
    verts = [0, 1, 2, 3]
    import itertools
    top = [tuple(c) for c in itertools.combinations(verts, 3)]  # the 4 triangular faces
    return betti_of_top_simplices(top, fields)


OUT = {}

t_start = time.time()

# ---- S^2 control (chi(S^2)=2, computed from boundary-of-tetrahedron) ----
OUT["S2_control"] = s2_control(FIELDS)

for N in (4, 6):
    t0 = time.time()
    quot_top, fixed_pts = build_quotient_top(N)
    fixed_set = set(fixed_pts)
    assert len(fixed_pts) == 16

    # ---- link of each singular vertex ----
    link_results = []
    for p in fixed_pts:
        link_top = [tuple(v for v in simp if v != p) for simp in quot_top if p in simp]
        # sanity: each should be a 3-simplex (4 vertices) since quot_top are 4-simplices (5 verts)
        bad = [s for s in link_top if len(s) != 4]
        r = betti_of_top_simplices(link_top, FIELDS)
        r["singular_vertex"] = p
        r["num_bad_facets"] = len(bad)
        link_results.append(r)

    # cross-check all 16 links give the same Betti numbers (RP^3 expected: (1,0,0,1) over Z/3, (1,1,1,1) over Z/2... we use Z/5 not Z/2 here for the primary field pair; add Z/2 separately below for the torsion control requested by the task)
    link_betti_z3 = [tuple(r["betti_by_field"]["Z3"]) for r in link_results]
    link_betti_z5 = [tuple(r["betti_by_field"]["Z5"]) for r in link_results]
    all_links_agree_z3 = len(set(link_betti_z3)) == 1
    all_links_agree_z5 = len(set(link_betti_z5)) == 1

    # torsion control: one representative link over Z/2
    p0 = fixed_pts[0]
    link_top_p0 = [tuple(v for v in simp if v != p0) for simp in quot_top if p0 in simp]
    link_p0_z2 = betti_of_top_simplices(link_top_p0, (2,))

    # ---- complement U: quotient top simplices touching NO singular vertex ----
    U_top = [simp for simp in quot_top if not (set(simp) & fixed_set)]
    U_result = betti_of_top_simplices(U_top, FIELDS)

    entry = {
        "N": N,
        "num_singular_vertices": len(fixed_pts),
        "link_results_per_vertex": link_results,
        "all_16_links_agree_Z3": all_links_agree_z3,
        "all_16_links_agree_Z5": all_links_agree_z5,
        "link_expected_betti_RP3_odd_field": [1, 0, 0, 1],
        "link_expected_source": "RP^3 over a field of odd characteristic: H0=H3=F, H1=H2=0 (standard RP^n homology)",
        "link_representative_vertex0_Z2_torsion_control": link_p0_z2,
        "link_expected_betti_RP3_Z2": [1, 1, 1, 1],
        "link_Z2_expected_source": "RP^3 over F2: all H_k (k=0..3) become rank-1 due to the Z/2 torsion in integral homology of RP^n",
        "U_complement_result": U_result,
    }

    # ---- Mayer-Vietoris derivation of resolved K3 Betti numbers ----
    for field, fname in ((3, "Z3"), (5, "Z5")):
        b_U = U_result["betti_by_field"][fname]
        b_U = list(b_U) + [0] * (5 - len(b_U))
        b2_U = b_U[2] if len(b_U) > 2 else 0
        b1_U = b_U[1] if len(b_U) > 1 else 0
        b2_K3 = b2_U + 16
        b1_K3 = b1_U
        b0_K3 = 1  # STATED, see module docstring
        b4_K3 = 1  # STATED (Poincare duality with b0)
        # b3, partially strengthened via Mayer-Vietoris degree 3 (still not
        # fully computed -- see caveat below):
        #   H3(16xRP3)=F^16 --f--> H3(U)+H3(16xS^2)=H3(U) -> H3(X) -> H2(16xRP3)=0
        # so b3(X) = b3(U) - rank(f), where f is the map induced by inclusion
        # of the 16 boundary RP^3's into U. b3(U) IS computed by GUDHI below.
        # rank(f) is NOT computable from Betti numbers alone (GUDHI gives no
        # chain-level map); it is only ARGUED (not verified) that f is
        # surjective (rank(f)=b3(U)), because the 16 boundary classes are
        # believed to generate H3(U) with the single relation that they
        # collectively bound the compact manifold-with-boundary U. IF that
        # surjectivity holds, b3(X)=b3(U)-b3(U)=0. This is consistent with,
        # but does NOT independently re-derive, the Poincare-duality value
        # b3=b1 used below; it is reported as an argued upper bound only.
        b3_U = b_U[3] if len(b_U) > 3 else 0
        b3_K3_upper_bound_if_f_surjective = 0  # b3(U) - rank(f), rank(f)=b3(U) ASSUMED, not verified
        b3_K3 = b1_K3  # reported number: STATED via Poincare duality (b3=b1), see docstring
        chi_K3_from_betti = b0_K3 - b1_K3 + b2_K3 - b3_K3 + b4_K3
        chi_U = U_result["euler_characteristic_from_simplex_counts"]
        chi_K3_from_MV_simplex_counts = chi_U + 32  # 16*(chi(S^2)-chi(RP^3)) = 16*(2-0)=32, both computed above
        entry[f"resolved_K3_betti_{fname}"] = {
            "b0": b0_K3, "b1": b1_K3, "b2": b2_K3, "b3": b3_K3, "b4": b4_K3,
            "b0_b4_status": "STATED (connectedness/orientability/Poincare duality), NOT computed",
            "b3_status": "STATED via Poincare duality (b3=b1). b3(U)=%d is COMPUTED (GUDHI); a Mayer-Vietoris argument gives b3(K3)=b3(U)-rank(f) for the boundary-inclusion map f, and IF f is surjective (ARGUED, not verified: no chain-level map is available from Betti numbers alone) this evaluates to 0, consistent with but not an independent re-derivation of b3=b1=%d." % (b3_U, b1_K3),
            "b1_b2_status": "COMPUTED from U via Mayer-Vietoris (b1=b1(U), b2=b2(U)+16)",
            "chi_from_betti_numbers": chi_K3_from_betti,
            "chi_from_MV_simplex_counts(chi(U)+32)": chi_K3_from_MV_simplex_counts,
            "chi_cross_check_matches": chi_K3_from_betti == chi_K3_from_MV_simplex_counts,
            "b2_plus_signature_note": "b2 = b2^+ + b2^- ; signs require the intersection form, computed separately in 04_signature.py",
        }

    entry["runtime_sec"] = time.time() - t0
    OUT[f"N={N}"] = entry

OUT["total_runtime_sec"] = time.time() - t_start

OUT_PATH = "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity/tda-gudhi/03_resolution_hybrid_results.json"
with open(OUT_PATH, "w") as f:
    json.dump(OUT, f, indent=2, default=str)

print(json.dumps({k: v for k, v in OUT.items() if k != "N=6"}, indent=2, default=str)[:6000])
print("\n... (N=6 section omitted from stdout for brevity, see file)")
print("\nWrote", OUT_PATH)
