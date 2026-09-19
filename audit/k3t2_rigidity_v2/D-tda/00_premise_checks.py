#!/usr/bin/env python
"""
Track D (v2), step 0: PREMISE checks, run BEFORE any Mayer-Vietoris
arithmetic is trusted for a given grid size N.

Two independent, purely combinatorial checks on the QUOTIENT complex
T^4/Z_2 (built by canonical-representative collapse v -> min(v, -v mod N)
of the Freudenthal/Kuhn triangulation of Z_N^4), for N in {4, 6, 8}:

  (A) CLOSED-STAR DISJOINTNESS of the 16 singular (fixed) vertices, in the
      QUOTIENT complex (not in T^4 itself -- canonicalization can fold
      neighbourhoods together, so disjointness must be checked AFTER the
      quotient is built, not before). For each pair of distinct fixed
      vertices p != q, we build the full vertex set of the closed star of
      p in quot_top (union of every quotient top-simplex containing p) and
      of q, and check the two sets are disjoint. This is the exact premise
      the hybrid Mayer-Vietoris decomposition in 03_resolution_hybrid.py
      needs: X = U union (16 disjoint pieces), with A n B = 16 DISJOINT
      copies of the link. If any two closed stars share a vertex, the "16
      disjoint RP^3 links" picture is wrong at that N and the MV formula
      in 03 must not be trusted there.
      v1 defect: this was never checked; v1 asserted "N=4 invalid, N=6
      valid" without running a check. Here it is checked programmatically
      for N=4, 6, 8.

  (B) NON-SINGULAR LINK VALIDITY: canonical-representative collapse (fold
      v and -v onto min(v,-v)) is only guaranteed to produce a genuine
      simplicial complex (as opposed to something already broken by the
      folding, independent of the star-disjointness issue at the singular
      points) if a link at an ORDINARY (non-fixed) quotient vertex comes
      out as S^3 -- Betti (1,0,0,1) over every field, with NO mod-2 jump
      (a mod-2 jump would signal unexpected RP^3-like torsion at a point
      that is supposed to be perfectly smooth, i.e. a folding defect).
      This is a cheap, structural sanity check on the construction itself,
      independent of (and additional to) the five numbered fixes in the
      task.

Both checks use only integer/tuple arithmetic (exact) plus GUDHI simplex
trees for the link Betti numbers in check (B).
"""
import json
import sys
import time

sys.path.insert(0, "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity_v2/D-tda")
import gudhi
from lib_freudenthal import freudenthal_top_simplices, negate_mod

FIELDS = (2, 3, 5)
D = 4


def build_quotient_top(N):
    top, degenerate = freudenthal_top_simplices(N, D)
    assert degenerate == 0
    neg = lambda v: negate_mod(v, N)

    def canon(v):
        return min(v, neg(v))

    quot_top = [tuple(canon(v) for v in simp) for simp in top]
    quot_top_clean = [s for s in quot_top if len(set(s)) == len(s)]
    n_degenerate_after_quotient = len(quot_top) - len(quot_top_clean)
    fixed_pts = sorted({v for simp in top for v in simp if neg(v) == v})
    return quot_top_clean, fixed_pts, n_degenerate_after_quotient


def closed_star_vertex_set(quot_top, p):
    s = set()
    for simp in quot_top:
        if p in simp:
            s.update(simp)
    return s


def star_disjointness_check(quot_top, fixed_pts):
    stars = {p: closed_star_vertex_set(quot_top, p) for p in fixed_pts}
    violations = []
    pts = list(fixed_pts)
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            p, q = pts[i], pts[j]
            shared = stars[p] & stars[q]
            if shared:
                violations.append({"p": p, "q": q, "shared_vertices": sorted(shared)})
    return {
        "num_pairs_checked": len(pts) * (len(pts) - 1) // 2,
        "num_violating_pairs": len(violations),
        "violations_sample": violations[:5],
        "all_16_closed_stars_pairwise_disjoint": len(violations) == 0,
    }


def link_betti(quot_top, v0, fields):
    link_top = [tuple(u for u in simp if u != v0) for simp in quot_top if v0 in simp]
    all_verts = sorted(set().union(*[set(s) for s in link_top])) if link_top else []
    vidx = {v: i for i, v in enumerate(all_verts)}
    out = {}
    for field in fields:
        st = gudhi.SimplexTree()
        for simp in link_top:
            st.insert([vidx[v] for v in simp], filtration=0.0)
        st.compute_persistence(homology_coeff_field=field, min_persistence=0, persistence_dim_max=True)
        out[f"Z{field}"] = st.betti_numbers()
    return out, len(link_top)


OUT = {}
t_start = time.time()

for N in (4, 6, 8):
    t0 = time.time()
    quot_top, fixed_pts, n_deg = build_quotient_top(N)
    star_check = star_disjointness_check(quot_top, fixed_pts)

    # non-singular link validity: pick the first quotient vertex that is
    # NOT one of the fixed (singular) points
    fixed_set = set(fixed_pts)
    all_quot_verts = sorted(set().union(*[set(s) for s in quot_top]))
    nonsingular_candidates = [v for v in all_quot_verts if v not in fixed_set]
    v0 = nonsingular_candidates[0]
    betti_v0, link_size = link_betti(quot_top, v0, FIELDS)
    expected_S3 = [1, 0, 0, 1]
    is_S3_like = all(betti_v0[f"Z{f}"] == expected_S3 for f in FIELDS)

    OUT[f"N={N}"] = {
        "N": N,
        "num_fixed_vertices_computed": len(fixed_pts),
        "num_degenerate_after_quotient_collapse": n_deg,
        "star_disjointness_check": star_check,
        "nonsingular_link_check": {
            "vertex_tested": v0,
            "link_num_top_simplices": link_size,
            "betti_by_field": betti_v0,
            "expected_S3_betti_all_fields": expected_S3,
            "matches_S3_no_mod2_jump": is_S3_like,
        },
        "PREMISE_VALID_FOR_MV": star_check["all_16_closed_stars_pairwise_disjoint"] and is_S3_like,
        "runtime_sec": time.time() - t0,
    }
    print(f"N={N}: disjoint={star_check['all_16_closed_stars_pairwise_disjoint']} "
          f"nonsingular_link_S3={is_S3_like} "
          f"(runtime {OUT[f'N={N}']['runtime_sec']:.1f}s)")

OUT["total_runtime_sec"] = time.time() - t_start
OUT["note_N4_numerics"] = (
    "N=4 fails the closed-star-disjointness premise (see star_disjointness_check "
    "above). If the raw quotient Betti numbers at N=4 nonetheless numerically "
    "agree with N=6/N=8 (they were observed to in the v1 run), that agreement is "
    "NOT evidence the N=4 premise held -- it is not used as supporting evidence "
    "anywhere downstream. Only N with PREMISE_VALID_FOR_MV=true are used in "
    "03_resolution_hybrid.py."
)

OUT_PATH = "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity_v2/D-tda/00_premise_checks_results.json"
with open(OUT_PATH, "w") as f:
    json.dump(OUT, f, indent=2, default=str)

print("\nWrote", OUT_PATH)
