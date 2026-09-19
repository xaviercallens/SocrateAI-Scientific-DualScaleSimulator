"""Merge all per-item result JSON files in this directory into one results.json,
and write exports.json (the small cross-track handoff file other tracks may poll).

Run: /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python \
     /mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity_v2/C-lattices/merge_results.py
"""
import json
import os

DIR = "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity_v2/C-lattices"

FILES = [
    "01_e8_result.json",
    "02_k3_mukai_gamma_result.json",
    "03_symbolic_tduality_proofs_result.json",
    "04_tadpole_budget_result.json",
    "06_chained_rigidity_result.json",
]

merged = {
    "track": "C", "version": "v2",
    "title": "K3 x T2 lattice/T-duality mathematics: symbolic exact-arithmetic checks and a "
             "CHAINED (topology -> Noether -> Hodge index -> lattice enumeration) structural "
             "rigidity test, no literal K3 numbers typed into any selecting condition",
    "python": "/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python",
    "worktree": "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2",
    "branch": "loop/k3t2-rigidity",
    "items": {},
    "could_not_do": [
        "Symbolic T-duality involution/spectrum check for the radius map R -> lambda*alpha'/R "
        "(v1 item 'rigidity_b') was not redone in v2: v1's 07_rigidity_b_tduality_radius.py "
        "already used fully symbolic sympy symbols (R,n,w,alpha') with no numeric sampling and "
        "had no literal-answer defect flagged for it, so it was left as a v1 result rather than "
        "duplicated; not re-verified independently in this v2 pass for time-budget reasons.",
        "D-tda/exports.json had NOT appeared during the initial 20-minute poll window (only "
        "A-genus/exports.json had); it appeared later, before this track's final run, and "
        "06_chained_rigidity.py's final result DOES use it (with a chi-value cross-check "
        "against A-genus's independent elliptic-genus computation, both giving chi=24) -- "
        "recorded here only to be explicit that the two-source corroboration was not available "
        "for the entire duration of this track's work, only from partway through.",
        "The Kahler hypothesis used in the chained derivation (step 8 of "
        "06_chained_rigidity.py, b2_plus = 2*h^{2,0}+1) is stated as a named input, not "
        "independently verified computationally in this track (it is a standard K3 fact, but "
        "no script here checks it from more primitive data). Simply-connectedness (b1=0) is, "
        "as of D-tda's later export, READ from an independent GUDHI computation rather than "
        "assumed, and cross-checked against the Betti alternating sum for b2 -- see "
        "06_chained_rigidity_result.json 'derivation_chain' step 5 and '10b_b2_cross_check'.",
        "b2's cross-check in 06_chained_rigidity.py ('10b_b2_internal_consistency_recheck_same_"
        "source') is NOT an independent-computation agreement: b2_read and the alternating-sum "
        "b2 are both derived from the SAME provenance file (D-tda/exports.json's b2_K3_resolved "
        "and its own b0/b1/b3/b4), which already flags its own internal consistency. Only "
        "chi_top is corroborated by two independently-computed sources (Track A elliptic genus, "
        "Track D GUDHI resolution, both giving 24) -- see 'chi_cross_track_corroboration'.",
        "The entry-level symbolic tie-back between the eigenvalue-formula minimizer proof (D1 "
        "in 03_symbolic_tduality_proofs.py) and an explicit symmetric-matrix family (D2) was "
        "done only for d=2; d=3,4 would need symbolic roots of a cubic/quartic characteristic "
        "polynomial, left out for time-budget reasons (see 03's "
        "'implicit_structural_inputs_for_D1_general_d' field for exactly what is and isn't "
        "covered).",
    ],
}

for fname in FILES:
    path = os.path.join(DIR, fname)
    with open(path) as f:
        merged["items"][fname] = json.load(f)

out_path = os.path.join(DIR, "results.json")
with open(out_path, "w") as f:
    json.dump(merged, f, indent=2)
print("wrote", out_path)

chained = merged["items"]["06_chained_rigidity_result.json"]
exports = {
    "track": "C",
    "derived_b2": chained["derived_b2"],
    "derived_signature_p_q": chained["derived_signature_p_q"],
    "selected_mn": chained["classification"]["final_solution_set_mn"],
    "chi_top_used": chained["provenance"]["chi_top_value"],
    "chi_top_source_file": chained["provenance"]["chi_top_read_from_file"],
    "classification": chained["classification"]["label"],
    "e8_cartan_det": None,  # filled below
    "e8_num_roots": None,
}
with open(os.path.join(DIR, "01_e8_result.json")) as f:
    e8 = json.load(f)
exports["e8_cartan_det"] = e8["determinant"]["value"]
exports["e8_num_roots"] = e8["root_enumeration"]["num_roots_found"]

exports_path = os.path.join(DIR, "exports.json")
with open(exports_path, "w") as f:
    json.dump(exports, f, indent=2)
print("wrote", exports_path)
