#!/usr/bin/env python
"""
Master summary: pulls the key computed numbers out of the six result
files this track wrote, and states the two-independent-routes findings
explicitly. Does not recompute anything; every number here is copied
from a file that was itself produced by a GUDHI/sympy run (see the
"source" field on each entry, or the numbered script of the same index).
"""
import json

BASE = "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity/tda-gudhi"

files = {}
for name in [
    "01_controls_results.json",
    "02_invariance_and_quotient_results.json",
    "03_resolution_hybrid_results.json",
    "04_signature_results.json",
    "05_kunneth_results.json",
    "06_pointcloud_results.json",
    "07_rigidity_scan_results.json",
]:
    with open(f"{BASE}/{name}") as f:
        files[name] = json.load(f)

s1, s2, s3, s4, s5, s6, s7 = (files[n] for n in files)

summary = {
    "tier": "B (exact simplicial/cubical arithmetic + a hybrid Mayer-Vietoris step) at best; never 'proved'.",
    "provenance_self_report": [
        "Ran `ls audit/` at session start to confirm the k3t2_rigidity output directory existed; this printed the names of blinded files (PAPER_FACTS.md, physics_claims_ledger.md, etc.) in that listing but none of them were opened, grepped, or read.",
        "A skill named 'k3t2-rigidity-loop' surfaced mid-session (its trigger fired on this task's topic); it was NOT invoked because its described workflow compares results against LeanMaster theorems, which the blindness rule for this run forbids. All computation here used only the mathematical definitions given in the task text.",
        "After committing this track's files, `git status` in the worktree showed three untracked paths belonging to sibling tracks/processes: audit/k3t2_rigidity/targets_sealed.json, audit/k3t2_rigidity/dyons/, audit/k3t2_rigidity/genus-moonshine/. None of these were opened, read, or grepped -- only their names appeared in the status listing, and targets_sealed.json in particular was recognized as likely containing the sealed answer key this blindness exercise is built around, so it was deliberately left untouched.",
    ],
    "route_A_simplicial_TDA": {
        "T2_control": s1["cubical"]["T2"],
        "T4_control_cubical": s1["cubical"]["T4"],
        "T4_control_simplicial_freudenthal": {
            "N=4": s1["freudenthal_simplicial_T4"]["4"]["betti"],
            "N=6": s1["freudenthal_simplicial_T4"]["6"]["betti"],
        },
        "T4_mod_Z2_quotient_singular": {
            "N=4": s2["N=4"]["quotient_T4_mod_Z2"]["betti_by_field"],
            "N=6": s2["N=6"]["quotient_T4_mod_Z2"]["betti_by_field"],
            "chi": {"N=4": s2["N=4"]["quotient_T4_mod_Z2"]["euler_characteristic_from_simplex_counts"],
                     "N=6": s2["N=6"]["quotient_T4_mod_Z2"]["euler_characteristic_from_simplex_counts"]},
        },
        "singular_link_is_RP3": {
            "N=4": s3["N=4"]["link_results_per_vertex"][0]["betti_by_field"],
            "N=4_all_16_agree": s3["N=4"]["all_16_links_agree_Z3"],
            "Z2_torsion_control": s3["N=4"]["link_representative_vertex0_Z2_torsion_control"]["betti_by_field"],
        },
        "resolved_K3_betti": {
            "N=4_Z3": s3["N=4"]["resolved_K3_betti_Z3"],
            "N=6_Z3": s3["N=6"]["resolved_K3_betti_Z3"],
        },
        "K3_x_T2_betti_kunneth": s5["N=4_Z3"]["betti_K3xT2_via_kunneth"],
    },
    "route_B_signature_eigenvalues": {
        "b2_plus": s4["b2_plus_total"],
        "b2_minus": s4["b2_minus_total"],
        "b2_total": s4["b2_total_(should_match_resolved_K3_b2_from_step3)"],
        "signature": s4["signature_b2plus_minus_b2minus"],
    },
    "TWO_INDEPENDENT_ROUTES_TO_rank6_invariant_part_of_H2": {
        "route_1_GUDHI_homology_of_computed_complement_U": s2["N=4"]["quotient_T4_mod_Z2"]["betti_by_field"]["Z3"]["betti"][2],
        "route_2_sympy_exact_rank_of_wedge_pairing_matrix": s4["invariant_H2_T4_basis_size"],
        "agree": s2["N=4"]["quotient_T4_mod_Z2"]["betti_by_field"]["Z3"]["betti"][2] == s4["invariant_H2_T4_basis_size"] == 6,
        "note": "This is the genuinely independent pair: route 1 is GUDHI simplicial homology of the singular quotient's b2 (=6, before resolution); route 2 is the sympy-exact dimension of Lambda^2(Q^4) (=6), the algebraic count of the invariant part of H^2(T^4;Q). Different computations (combinatorial homology vs. exact-arithmetic linear algebra), same rank-6 answer.",
        "caveat": "The subsequent step of adding 16 exceptional classes to reach b2(K3)=22 is a STATED input shared by both routes (route 1 via the Mayer-Vietoris complement-plus-16 argument, route 2 via 16 orthogonal (-2)-classes), so the final '22' is not two fully independent derivations of 22 -- only the rank-6 piece is. What route 2 uniquely adds, and route 1 cannot produce at all, is the SIGNATURE SPLIT (b2+,b2-)=(3,19), sigma=-16.",
    },
    "TWO_INDEPENDENT_ROUTES_TO_T2_TOPOLOGY": {
        "route_1_cubical_and_simplicial_complex": s1["cubical"]["T2"]["betti"],
        "route_2_point_cloud_persistent_homology": {
            "torus_H1_gap_ratio_raw": s6["clifford_torus"]["H1_gap_ratio_after_2nd_bar (life[1]/life[2])"],
            "torus_H2_gap_ratio_raw": s6["clifford_torus"]["H2_gap_ratio_after_1st_bar (life[0]/life[1])"],
            "null_ball_H1_gap_ratio_raw": s6["null_control_uniform_ball"]["H1_gap_ratio_after_2nd_bar (life[1]/life[2])"],
            "null_ball_H2_gap_ratio_raw": s6["null_control_uniform_ball"]["H2_gap_ratio_after_1st_bar (life[0]/life[1])"],
            "threshold_used_for_the_boolean_flags_below": 3.0,
            "threshold_honesty_note": "any threshold in the wide window (1.1, 4.7) separates torus (4.77, 200.5) from ball (1.03, 1.04) identically; the conclusion does not depend on the specific value chosen.",
            "signal_clear": s6["torus_signal_clear"],
            "null_control_clean": s6["null_control_clean_(gap_below_threshold)"],
            "deliberate_substitution": "500 points (not the suggested 1500) with AlphaComplex (not RipsComplex): a Rips expansion to homological dimension 3 in R^4 at 1500 points is not feasible in the time budget; AlphaComplex (exact Delaunay-based) needs no max_edge_length truncation and completed in ~3s at 500 points.",
        },
        "note": "Route 1: exact combinatorial complex, filtration=0, gives (1,2,1) directly. Route 2: 500 points sampled on a Clifford torus in R^4, alpha-complex persistent homology, finds exactly 2 long H1 bars and 1 long H2 bar separated from noise by gap ratios >>3x, with a uniform-ball null control showing no such gap. Independent methods, same topological conclusion (b1=2, b2=1).",
    },
    "RIGIDITY_TESTS": {
        "1_resolution_count_k": {
            "parameter_scanned": "k = number of the 16 singular points resolved, k in {0,...,16}",
            "solution_set_chi_equals_24": {"N=4": s7["N=4"]["solution_set_chi_equals_24"], "N=6": s7["N=6"]["solution_set_chi_equals_24"]},
            "negative_control_k=15_chi": {"N=4": s7["N=4"]["negative_control_k=num_fixed-1"]["chi"], "N=6": s7["N=6"]["negative_control_k=num_fixed-1"]["chi"]},
            "second_independent_criterion": "smoothness (no remaining orbifold points) also forces k=16",
            "conclusion": "Both the chi=24 arithmetic condition and the smoothness condition force k=16 uniquely, out of {0,...,16}: zero free parameters for the resolution count, agreeing at N=4 and N=6.",
        },
        "2_group_action_inserted_as_a_parameter": {
            "condition": "negation x->-x (16 fixed points) vs. free translation (0 fixed points)",
            "negation_result": "singular quotient betti (1,0,6,0,1), chi=8 -- the T^4/Z2 orbifold, matches the fixed-point count exactly (see exact_halving_check in 02_..._results.json)",
            "free_translation_result_N=6": "recovers exactly the T^4 betti (1,4,6,4,1), chi=0, at both Z3 and Z5 -- a genuine smooth free quotient, correctly distinguished from the singular case",
            "free_translation_result_N=4_caveat": "the specific shift (2,0,0,0) has quotient period 2 on one axis, too small for a non-degenerate Freudenthal quotient at N=4 (known discretization limitation, distinct from the negation action's behavior, which IS clean at N=4); N=6 gives the clean instance",
            "conclusion": "The pipeline is demonstrably sensitive to WHICH group action is inserted -- it does not output a fixed answer regardless of input.",
        },
        "3_triangulation_refinement_N_as_a_would_be_parameter": {
            "values_scanned": [4, 6],
            "result": "every Betti number, chi, and the resolved K3 numbers (b2=22, signature=-16) are IDENTICAL at N=4 and N=6",
            "conclusion": "N is an inverse/negative test: a real discretization choice that COULD have mattered, and does not move any reported topological number -- evidence the numbers are properties of the topology, not artifacts of the specific grid size. This is an invariance check, not itself a rigidity constraint.",
        },
    },
    "negative_controls_summary": {
        "free_translation_quotient_N=6_matches_T4_betti": s2["N=6"]["negative_control_free_translation_quotient"]["matches_expected_Z3"],
        "free_translation_quotient_N=4_shift(2,0,0,0)_caveat": "This particular shift has quotient period 2 along one axis, too small for the Freudenthal triangulation to remain non-degenerate under canonical-representative collapse (a DIFFERENT failure mode from the negation action, which is fine at N=4); betti came out as (1,3,69,4,1), chi=64, NOT matching T^4. This is an honestly-reported limitation of the N=4 grid for THIS SPECIFIC shift, not evidence against the main T^4/Z2 quotient (which used the negation action, verified non-degenerate and matching at both N=4 and N=6).",
        "rigidity_scan_k=15_fails": True,
    },
}

with open(f"{BASE}/00_master_summary.json", "w") as f:
    json.dump(summary, f, indent=2, default=str)

print(json.dumps(summary, indent=2, default=str))
