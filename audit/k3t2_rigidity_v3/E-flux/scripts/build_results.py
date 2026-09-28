#!/usr/bin/env python3
"""E-flux v3: assemble results.json, exports.json and inputs.json from the per-script JSON files in results/.
Never hand-edited.  Run after scripts 01..06 and poll_track_d.py:
  cd audit/k3t2_rigidity_v3/E-flux/scripts && <venv python> build_results.py
"""
import json
from common import *
from inputs_decl import inputs

FILES = ["00_track_d_poll.json", "01_lattice_index.json", "02_tt_worked_examples.json", "03_family41_counts.json", "04_family42_counts.json", "05_second_branch.json", "06_literature_and_moduli.json"]

def main():
    inp = inputs()
    names = {i["name"] for i in inp}
    results, rigidity, lit, could, prov = [], [], [], [], {}
    poll = json.load(open(RES / FILES[0]))
    for f in FILES[1:]:
        d = json.load(open(RES / f))
        for r in d["results"]:
            bad = [x for x in r["shared_inputs"] if x not in names]
            assert not bad, (r["id"], bad)
            results.append(r)
        rigidity += d.get("rigidity", []); lit += d.get("literature_checks", []); could += d.get("could_not_do", [])
        prov[f] = {"script": d["script"], "wall_seconds": d.get("wall_seconds")}
    # T4Z2 count consumed from Track D (declared input value filled from the export)
    for i in inp:
        if i["name"] == "T4Z2_fixed_points_from_track_D":
            i["value"] = f"{poll['T4Z2_fixed_points_from_track_D']} (read from {poll['track_d_export_used']}; also computed here tier B: {poll['T4Z2_fixed_points_computed_here_tierB']})"
    by = {r["id"]: r["computed"] for r in results}
    ex = {"track": "E-flux", "generated_by": "build_results.py", "note": "every value below is read from a results/*.json file produced by a script; none typed",
          "index_of_A6_span_in_unimodular_Gamma33": {"value": by["E01b"]["index_from_transition_matrix_det"], "file": "audit/k3t2_rigidity_v3/E-flux/results/01_lattice_index.json"},
          "alpha_x_sq_residues_mod_8_even_lattice": {"value": by["E01e2"]["residues_alpha_sq_mod_8_even_lattice_U3"], "control_odd": by["E01e2"]["residues_alpha_sq_mod_8_odd_control_I33"], "file": "audit/k3t2_rigidity_v3/E-flux/results/01_lattice_index.json"},
          "alpha_x_sq_values_with_tadpole_family41": {"value": by["E03c"]["even_lattice_U3_alpha_sq_values"], "N_D3": by["E03c"]["N_D3_even"], "file": "audit/k3t2_rigidity_v3/E-flux/results/03_family41_counts.json"},
          "family41_O(Gamma33)_classes": by["E03d"]["total_classes_with_tadpole"], "family41_O(Gamma33)_classes_per_m": by["E03d"]["per_window_(m: (upper,lower))"],
          "second_branch_dim1_allowed_(j,k)": by["E05a"]["allowed_(j,k)_from_scan_up_to_KMAX"], "second_branch_dim1_N_D3": by["E05a"]["N_D3_values"],
          "N_D3_for_TT_examples": {"4.15": by["E02a"]["N_D3_with_half_(2.3)"], "4.31": by["E02b"]["N_D3_with_half"], "5.8": by["E02c"]["N_D3_with_half"]},
          "N_flux_for_N_D3_zero_(2.3)": by["E02d"]["N_flux_for_N_D3=0_with_half_(2.3)"],
          "T4Z2_singular_points_consumed_from_D": poll["T4Z2_fixed_points_from_track_D"], "O7_planes_computed_tierB": poll["T2Z2_fixed_points_computed_here_tierB"],
          "tadpole_total_input_TT_2.3": next(i["value"] for i in inp if i["name"] == "tadpole_total_24_eq2.3"), "tadpole_total_tier": "L (typed input)"}
    out = {"track": "E-flux", "version": "v3", "generated_by": "build_results.py", "provenance": prov, "sources": {"TT_txt_sha256": sha256(TT_TXT), "BLPSSW_txt_sha256": sha256(BL_TXT), "TT_pdf_sha256_in_SHA256SUMS": "ee30a5aa021708a9d340864dcb071651b2e0f366578cedfc81b570f9c2838f18"},
           "results": results, "rigidity": rigidity, "literature_checks": lit, "could_not_do": could + poll_notes(poll) + EXTRA}
    (EDIR / "results.json").write_text(json.dumps(out, indent=1, default=str))
    (EDIR / "exports.json").write_text(json.dumps(ex, indent=1, default=str))
    (EDIR / "inputs.json").write_text(json.dumps(inp, indent=1, default=str))
    print("wrote results.json", len(results), "results;", len(rigidity), "rigidity;", len(lit), "literature checks")

EXTRA = ["orbits of O(Gamma_{3,19}) (full lattice) not computed: orbit classes are for O(Gamma_{3,3}) acting on U^3 only; rank-8 truncation (U^3 + A2(-1)) has Burnside counts mod G only, no O(Gamma) classes",
         "section 4.2: exact verification of (3.16a-c), (3.31), V_flux type and the orbifold test done on a random 400-representative sample per window (seed 7), float classification used for all others (0 disagreements on the sample); the orbifold-free fraction is therefore a sample estimate",
         "first Track D poll (1200 s) matched no key (v3 names it n_singular) and fell back to v2; re-run with 0 s picked up the v3 export, agreeing (16). The committed command is poll_track_d.py 0",
         "S-duality, T-duality/SL(2,Z) of T2 and O(Gamma_{3,19}) are not quotiented in any 'mod G' count; the tadpole 24 is a typed tier-L input, not derived",
         "1/2 vs line 941 decided from three worked examples only; no independent derivation of the 1/2 (would need the orientifold charge computation)"]

def poll_notes(poll):
    return [] if not poll.get("fallback_to_v2") else ["Track D v3 export not present after polling; fell back to v2 export (recorded in 00_track_d_poll.json)"]

if __name__ == "__main__":
    main()
