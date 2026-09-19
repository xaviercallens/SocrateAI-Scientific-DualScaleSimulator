#!/usr/bin/env python
"""
Track D (v2) master summary: pulls the key computed numbers out of the six
result files this track wrote, in order, with no recomputation. Every
number has a "source" pointer to the file and key it came from.
"""
import json

BASE = "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity_v2/D-tda"

files = {}
for name in [
    "00_premise_checks_results.json",
    "01_controls_results.json",
    "02_invariance_and_quotient_results.json",
    "03_resolution_hybrid_results.json",
    "04_kunneth_results.json",
    "05_rigidity_scan_results.json",
]:
    with open(f"{BASE}/{name}") as f:
        files[name] = json.load(f)

s0, s1, s2, s3, s4, s5 = (files[n] for n in files)

summary = {
    "tier": "B (exact simplicial arithmetic + a hybrid Mayer-Vietoris step) at best; never 'proved'.",
    "blindness_self_report": [
        "No file under SocrateAI-Scientific-Agora-LeanMaster, k3t2-sealed-v2, proofs/, lean_foundation/, "
        "audit/PAPER_FACTS.md, or the v1 sealed-comparison files (targets_sealed.json, comparison.json, "
        "agreement_ledger.json, build_ledger.py, REPORT.md, reverse/, skeptic/) was opened, read, grepped or "
        "listed in this session. v1 CODE in genus-moonshine/, dyons/, lattices-duality/, tda-gudhi/ was read "
        "and reused (as licensed), but not its sealed comparison outputs.",
    ],
    "premise_checks": {
        "closed_star_disjointness": {N: s0[f"N={N}"]["star_disjointness_check"]["all_16_closed_stars_pairwise_disjoint"] for N in (4, 6, 8)},
        "nonsingular_link_is_S3": {N: s0[f"N={N}"]["nonsingular_link_check"]["matches_S3_no_mod2_jump"] for N in (4, 6, 8)},
        "N4_excluded_because": "closed stars of the 16 singular points are NOT pairwise disjoint at N=4 (computed, see 00_premise_checks_results.json); this despite N=4's raw Betti numbers numerically matching N=6/N=8 -- that numerical agreement is explicitly NOT treated as evidence the N=4 premise held.",
        "valid_N_used_downstream": [6, 8],
        "source": "00_premise_checks_results.json",
    },
    "resolved_K3_betti_numbers": {
        "N=6": {k: s3["N=6"]["resolved_K3_betti_Z3"][k] for k in ("b0", "b1", "b2", "b3", "b4")},
        "N=8": {k: s3["N=8"]["resolved_K3_betti_Z3"][k] for k in ("b0", "b1", "b2", "b3", "b4")},
        "labeling": "b1,b2 COMPUTED (Mayer-Vietoris from GUDHI Betti numbers of U and links); b0,b3,b4 STATED (connectedness/Poincare duality), unchanged from v1",
        "cross_field_Z5_agrees": {N: s3[N]["resolved_K3_betti_Z5"]["b2"] == s3[N]["resolved_K3_betti_Z3"]["b2"] for N in ("N=6", "N=8")},
        "chi_fvector_vs_betti_cross_check": {N: s3[N]["resolved_K3_betti_Z3"]["chi_cross_check_matches"] for N in ("N=6", "N=8")},
        "torsion_control_RP3_mod2_jump_confirmed": {N: s3[N]["torsion_control_shows_mod2_jump"] for N in ("N=6", "N=8")},
        "source": "03_resolution_hybrid_results.json",
    },
    "K3_x_T2_kunneth": {
        "betti": s4["N=6_Z3"]["betti_K3xT2_via_kunneth"],
        "chi": s4["N=6_Z3"]["chi_K3xT2_from_product_betti"],
        "chi_cross_check_matches_multiplicativity": s4["N=6_Z3"]["chi_cross_check_matches"],
        "all_N_and_fields_agree": s4["all_N_and_fields_agree_on_betti"],
        "source": "04_kunneth_results.json",
    },
    "rigidity_scan_over_k": {
        "closed_form": s5["N=6"]["closed_form"],
        "chi_equals_24_classification": s5["N=6"]["solution_set_chi_equals_24"]["classification"],
        "chi_equals_24_solution_set": s5["N=6"]["solution_set_chi_equals_24"]["set"],
        "b1_selector_classification": s5["N=6"]["b1_selector_scan"]["classification"],
        "b1_selector_solution_set_size": s5["N=6"]["b1_selector_scan"]["set_size"],
        "even_intersection_form_criterion_status": s5["N=6"]["even_intersection_form_criterion"]["status"],
        "manifoldness_condition_classification": s5["N=6"]["manifoldness_condition_classification"]["classification"],
        "conclusion": s5["N=6"]["rigidity_conclusion"],
        "agrees_at_N8": s5["N=8"]["solution_set_chi_equals_24"]["set"] == s5["N=6"]["solution_set_chi_equals_24"]["set"],
        "source": "05_rigidity_scan_results.json",
    },
    "rigidity_over_N": {
        "parameter": s5["rigidity_over_N"]["parameter"],
        "selecting_condition": s5["rigidity_over_N"]["selecting_condition"],
        "negative_control": s5["rigidity_over_N"]["negative_control"],
        "invariant_across_premise_valid_N": s5["rigidity_over_N"]["invariant_across_premise_valid_N"],
        "classification": s5["rigidity_over_N"]["classification"],
        "source": "05_rigidity_scan_results.json (rigidity_over_N)",
    },
    "translation_negative_control": {
        "N=6": {
            "dropped": s2["N=6"]["negative_control_free_translation_quotient"]["DROPPED"],
            "passes_halving_check": s2["N=6"]["negative_control_free_translation_quotient"].get("passes_same_halving_check_as_negation_quotient"),
            "matches_T4_betti": s2["N=6"]["negative_control_free_translation_quotient"].get("matches_expected_Z3"),
        },
        "N=8": {
            "dropped": s2["N=8"]["negative_control_free_translation_quotient"]["DROPPED"],
            "passes_halving_check": s2["N=8"]["negative_control_free_translation_quotient"].get("passes_same_halving_check_as_negation_quotient"),
            "matches_T4_betti": s2["N=8"]["negative_control_free_translation_quotient"].get("matches_expected_Z3"),
        },
        "source": "02_invariance_and_quotient_results.json",
    },
}

OUT_PATH = f"{BASE}/results.json"
with open(OUT_PATH, "w") as f:
    json.dump(summary, f, indent=2)

print(json.dumps(summary, indent=2))
print("\nWrote", OUT_PATH)
