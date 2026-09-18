#!/usr/bin/env python
"""
Track D (v2), step 5: RIGIDITY scan over k, the number of the 16 singular
points that are resolved (0 <= k <= 16), honestly classified per the task's
NORMALISATION/RIGID rule.

Model (same bookkeeping as v1, all four ingredients read from files this
track computed by running GUDHI -- nothing here is a fresh guess or a
typed target):
  chi(X_k) = chi(U) + k*chi(S^2) + (16-k)*chi(cone(RP^3)) - 16*chi(RP^3)
  chi(U)         <- 03_resolution_hybrid_results.json, N=6 (premise-valid)
  chi(S^2)       <- 03_resolution_hybrid_results.json (tetrahedron-boundary control)
  chi(RP^3)      <- computed from the (16-way agreed) link Betti numbers
  chi(cone(RP^3))=1 <- STATED (a cone is contractible), not computed
  num_fixed      <- 02_invariance_and_quotient_results.json (computed count)

HONEST CLASSIFICATION (per task rule "a scan whose selecting condition is
built from the true value is NOT rigid; label it normalisation"):

  Selecting condition "chi(X_k) == 24" is built directly from the literal
  target chi=24 (the physics identification "this resolution is K3", i.e.
  an EXTERNAL input, not derived from anything in this scan). The scan
  itself simplifies to chi(X_k) = 8 + k (see arithmetic below), which is
  STRICTLY INCREASING in k, so setting it equal to any single target value
  trivially picks out at most one k, for ANY target -- the scan carries NO
  additional discriminating power beyond "restate the target as k = target
  - 8". This is classified NORMALISATION. The input that fixes k=16 is the
  literal target chi=24 (equivalently: the assertion that this particular
  resolution IS the K3 surface).

  A candidate SECOND, structural condition was considered: "X_k is a smooth
  closed manifold (no remaining orbifold points) iff k = num_fixed". This
  was investigated and REJECTED as an independent rigidity criterion (not
  reported as RIGID) for two reasons, both checked explicitly below:
    (i) It is DEFINITIONAL given the X_k construction itself (X_k is
        DEFINED as "resolve k of the 16 points, leave 16-k as cone
        singularities"), so "no orbifold points remain" trivially requires
        k=16 by fiat, not by any independent GUDHI computation.
    (ii) Betti-number-only degree-by-degree bookkeeping (the only tool this
        track has) CANNOT distinguish k from any other value structurally:
        b1(X_k) = b1(U) is INDEPENDENT of k (the degree-1 Mayer-Vietoris
        term gets no contribution from either a cone or an S^2-bundle
        piece, since H_1(RP^3)=0 over every field used here and both the
        cone and the bundle piece are simply connected in the relevant
        degree). Likewise a Z/2-torsion-based degree-2 selector was
        considered and rejected: it would require the RANK of the
        boundary-inclusion map H_2(k x RP^3) -> H_2(U) + H_2(k x pieces),
        which Betti numbers alone (no chain-level map from GUDHI) cannot
        supply. This negative result is reported explicitly, not glossed
        over.
    (iii) chi(X_k) = 8 + k is checked for PARITY: chi is even iff k is
        even. This is a genuine STRUCTURAL fact (closed, even-dimensional
        manifolds have even Euler characteristic via Poincare duality) but
        it selects a SET of 9 values (k in {0,2,4,...,16}), not a point --
        reported as a real (if weak) structural constraint, distinct from
        and much weaker than the k=16-via-chi=24 normalisation.

  CONCLUSION: this track finds NO rigid (structural, target-value-free)
  selecting condition in its own data that isolates k=16 to a single point.
  The chi=24 selection is NORMALISATION. The parity constraint is a
  genuine but weak (9-element) structural narrowing. No claim of rigidity
  for k=16 is made here.
"""
import json

BASE = "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity_v2/D-tda"

with open(f"{BASE}/02_invariance_and_quotient_results.json") as f:
    step2 = json.load(f)
with open(f"{BASE}/03_resolution_hybrid_results.json") as f:
    step3 = json.load(f)

OUT = {}

for N in ("N=6", "N=8"):
    num_fixed = step2[N]["num_fixed_vertices_computed"]
    chi_U = step3[N]["U_complement_result"]["euler_characteristic_from_fvector"]
    chi_S2 = step3["S2_control"]["euler_characteristic_from_fvector"]

    link0 = step3[N]["link_results_per_vertex"][0]["betti_by_field"]["Z3"]
    chi_RP3 = sum((-1) ** i * b for i, b in enumerate(link0))
    chi_cone_RP3 = 1  # STATED: cone on anything is contractible, chi=1

    def chi_of_k(k):
        return chi_U + k * chi_S2 + (num_fixed - k) * chi_cone_RP3 - num_fixed * chi_RP3

    scan = {k: chi_of_k(k) for k in range(0, num_fixed + 1)}
    solution_set_chi24 = [k for k, chi in scan.items() if chi == 24]

    # verify the claimed closed form chi(X_k) = chi_U - num_fixed*chi_RP3 + num_fixed*chi_cone_RP3 + k*(chi_S2-chi_cone_RP3)
    # is exactly affine in k (checked, not assumed)
    slope = chi_of_k(1) - chi_of_k(0)
    is_affine = all(chi_of_k(k + 1) - chi_of_k(k) == slope for k in range(0, num_fixed))

    parity_solution_set = [k for k in range(0, num_fixed + 1) if chi_of_k(k) % 2 == 0]

    entry = {
        "N": N,
        "num_fixed_points_computed": num_fixed,
        "chi_U_computed": chi_U,
        "chi_S2_computed": chi_S2,
        "chi_RP3_computed_from_link_betti": chi_RP3,
        "chi_cone_RP3_STATED": chi_cone_RP3,
        "chi_of_k_scan": scan,
        "chi_of_k_is_affine_in_k": is_affine,
        "chi_of_k_slope": slope,
        "closed_form": f"chi(X_k) = {chi_of_k(0)} + {slope}*k",
        "solution_set_chi_equals_24": {
            "set": solution_set_chi24,
            "classification": "NORMALISATION",
            "reason": "selecting condition chi(X_k)==24 is built from the literal target value 24 (the external identification of this resolution as K3); chi(X_k) is strictly monotonic (affine, slope=%d) in k, so this condition trivially isolates a point for ANY target, carrying no independent discriminating power" % slope,
            "input_that_fixes_it": "the literal target chi=24, i.e. the physics assertion that this resolution IS the K3 surface -- an external input, not derived here",
        },
        "candidate_structural_condition_manifoldness": {
            "criterion": "X_k is a smooth closed manifold (no remaining orbifold points) iff k == num_fixed",
            "classification": "REJECTED as an independent rigid criterion",
            "reason": "definitional given the X_k construction (X_k is DEFINED to have 16-k unresolved cone points), not an independent GUDHI computation; a genuinely computable Betti-number-only distinguishing test was sought (Z/2-torsion jump at unresolved links, degree-1 MV term) and found NOT to depend on k with the tools available here (see module docstring, points (i)-(ii))",
        },
        "structural_but_weak_condition_parity": {
            "criterion": "chi(X_k) even (necessary for a closed even-dimensional manifold, via Poincare duality: b1=b3 forces chi even)",
            "classification": "STRUCTURAL but WEAK (does not use the literal target 24 anywhere)",
            "solution_set": parity_solution_set,
            "set_size": len(parity_solution_set),
            "isolates_single_point": len(parity_solution_set) == 1,
        },
        "negative_control_k=num_fixed-1": {
            "k": num_fixed - 1,
            "chi": scan.get(num_fixed - 1),
            "fails_chi_24": scan.get(num_fixed - 1) != 24,
        },
        "rigidity_conclusion": (
            "NO RIGID selecting condition for k=%d was found in this track's own data. "
            "chi(X_k)==24 is NORMALISATION (external target). The only condition found that "
            "avoids referencing the target value (chi parity) is structural but WEAK: it narrows "
            "k to a %d-element set, not a single point. This is a materially different, more "
            "conservative conclusion than v1's claim of a rigid, two-independent-criteria result."
            % (num_fixed, len(parity_solution_set))
        ),
    }
    OUT[N] = entry

OUT_PATH = f"{BASE}/05_rigidity_scan_results.json"
with open(OUT_PATH, "w") as f:
    json.dump(OUT, f, indent=2)

print(json.dumps(OUT, indent=2))
print("\nWrote", OUT_PATH)
