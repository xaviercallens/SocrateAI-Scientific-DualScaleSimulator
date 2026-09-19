#!/usr/bin/env python
"""Track D v3, step 4: Kunneth for K3 x T^2, scan over the number j of resolved singular points,
rigidity classification, inputs.json, results.json, exports.json (all WRITTEN BY THIS SCRIPT).

Run:  cd audit/k3t2_rigidity_v3/D-tda && /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python 04_kunneth_scan_results.py
(no arguments; reads 00_premise_results.json, 01_controls_results.json, 02_mv_results.json, 03_cup_parity_results.json)
"""
import json
from lib_common import HERE

J = lambda n: json.loads((HERE / n).read_text())
prem, ctl, mv, cup = J("00_premise_results.json"), J("01_controls_results.json"), J("02_mv_results.json"), J("03_cup_parity_results.json")
CMD = lambda s, a="": f"cd audit/k3t2_rigidity_v3/D-tda && /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python {s} {a}".strip()

inputs = [
 {"name": "involution", "value": "v -> -v on Z_N^4, N even; T^4/Z2 as the canonical-representative image of the Kuhn triangulation", "tier": "definition", "why": "defines the orbifold; regularity of the image complex is CHECKED in 00_premise.py (f-vector halving), not assumed"},
 {"name": "kuhn_triangulation", "value": "Freudenthal/Kuhn triangulation of (R/NZ)^4, N in {4,6,8} (odd N and N=3,5,7 used as controls)", "tier": "definition", "why": "the simplicial model; validated by the T^4 Betti control (1,4,6,4,1)"},
 {"name": "resolution_local_model", "value": "each singular point is replaced by a disc bundle over S^2 (Eguchi-Hanson T*S^2) whose boundary is the RP^3 link; homotopy type S^2, Betti (1,0,1); H_0 of the link maps isomorphically to H_0 of its own piece", "tier": "L (FROM MEMORY; standard Kummer construction)", "why": "the resolved space is not triangulated here; Mayer-Vietoris needs the homology of the glued piece. b2, b4, chi of the resolved space are (computed data of Q) + (this input)"},
 {"name": "coefficient_fields", "value": "Z/3 and Z/5 for Mayer-Vietoris (odd, so H_1=H_2=0 for RP^3); Z/2 for the link/Q/U Betti numbers only", "tier": "definition", "why": "with Z/2 the link has H_1,H_2 != 0 and phi_1, phi_2 ranks would be needed (not computed)"},
 {"name": "target_chi_K3", "value": 24, "tier": "L (FROM MEMORY)", "why": "used ONLY in the normalisation column of the j-scan and in 'expected' fields; never in a computation path"},
 {"name": "T2_betti_source", "value": "gudhi.PeriodicCubicalComplex on T^2 (computed in 01_controls.py); not typed", "tier": "B", "why": "input to Kunneth"},
]

# ---------- per-N MV summary ----------
Ns = [n for n in (4, 6, 8) if "Z3" in mv.get(f"N={n}", {})]
per_N = {}
for n in Ns:
    e = mv[f"N={n}"]
    per_N[n] = {"n_singular": e["n_singular"], "chi_Q": e["chi_Q_from_f"], "betti_Q_Z3": e["Z3"]["GUDHI_betti_Q"], "betti_U_Z3": e["betti_U"]["Z3"],
                "rank_phi3": {F: e[F]["rank_phi3_used"] for F in ("Z3", "Z5")},
                "rank_phi3_source": {F: e[F]["rank_phi3_source"] for F in ("Z3", "Z5")},
                "chain_equals_identity": {F: e[F]["chain_equals_identity"] for F in ("Z3", "Z5")},
                "MV_singular_reproduces_GUDHI_Q": {F: e[F]["MV_singular_reproduces_GUDHI_Q"] for F in ("Z3", "Z5")},
                "resolved_betti": {F: e[F]["resolved_betti_b0_b4"] for F in ("Z3", "Z5")},
                "chi_from_betti": {F: e[F]["chi_resolved_from_betti"] for F in ("Z3", "Z5")},
                "chi_from_fvectors": {F: e[F]["chi_resolved_from_fvectors"] for F in ("Z3", "Z5")},
                "chi_agree": {F: e[F]["chi_betti_equals_fvector"] for F in ("Z3", "Z5")}}
allB = {tuple(v["resolved_betti"][F]) for v in per_N.values() for F in ("Z3", "Z5")}
allchi = {v["chi_from_betti"][F] for v in per_N.values() for F in ("Z3", "Z5")}
b = list(next(iter(allB)))
chi_K3 = sum((-1) ** i * x for i, x in enumerate(b))
n_sing = {v["n_singular"] for v in per_N.values()}

# ---------- Kunneth ----------
bT2 = ctl["cubical_T2"]["betti_Z3"]
bprod = [0] * (len(b) + len(bT2) - 1)
for i, x in enumerate(b):
    for j, y in enumerate(bT2): bprod[i + j] += x * y
chiT2 = sum((-1) ** i * x for i, x in enumerate(bT2))
kun = {"b_K3": b, "b_T2": bT2, "b_K3xT2": bprod, "chi_from_product_betti": sum((-1) ** i * x for i, x in enumerate(bprod)),
       "chi_K3_times_chi_T2": chi_K3 * chiT2, "field_coefficients": "Kunneth over a field (no Tor)"}
kun["consistent"] = kun["chi_from_product_betti"] == kun["chi_K3_times_chi_T2"]

# ---------- j-scan (every integer 0..k) at each N, Z/3 ----------
scan = {}
for n in Ns:
    e = mv[f"N={n}"]; k = e["n_singular"]; sc = e["Z3"]["scan_over_j_resolved"]
    rows = []
    for j in range(k + 1):
        bj = sc[str(j)]["betti"]
        rows.append({"j": j, "betti": bj, "chi": sc[str(j)]["chi_from_betti"], "b1": bj[1], "poincare_duality_Z3": bj[1] == bj[3] and bj[0] == bj[4],
                     "chi_equals_target": sc[str(j)]["chi_from_betti"] == inputs[4]["value"],
                     "unresolved_points_with_Z2_link_not_S3": k - j})
    affine = all(rows[j + 1]["chi"] - rows[j]["chi"] == rows[1]["chi"] - rows[0]["chi"] for j in range(k))
    scan[f"N={n}"] = {"rows": rows, "chi_affine_in_j": affine, "slope": rows[1]["chi"] - rows[0]["chi"], "intercept": rows[0]["chi"],
                      "sol_chi_eq_target": [r["j"] for r in rows if r["chi_equals_target"]],
                      "sol_b1_zero": [r["j"] for r in rows if r["b1"] == 0],
                      "sol_poincare_duality_Z3": [r["j"] for r in rows if r["poincare_duality_Z3"]],
                      "sol_no_unresolved_point_with_nonsphere_Z2_link": [r["j"] for r in rows if r["unresolved_points_with_Z2_link_not_S3"] == 0]}
p0 = prem["N=4"]
rig = [
 {"parameter_inserted": "grid size N (with the involution fixed)", "selecting_condition": "PREMISE_VALID_FOR_MV: image complex regular (f-halving), open stars of singular points disjoint, links RP^3 / S^3, no orbit collisions (00_premise.py); does not mention chi, b2 or 24",
  "condition_uses_true_value": False, "solution_set": [n for n in (3, 4, 5, 6, 7, 8) if prem[f"N={n}"]["PREMISE_VALID_FOR_MV"]],
  "classification": "INVARIANCE", "note": "the answer (b0..b4, chi, n_singular) is identical for every valid N and both fields; the premise is a validity check, not a selection of N",
  "negative_control": "odd N=3,5,7 (same parameter): premise fails (regular f-halving false; one fixed vertex, orbit collisions); realizable complexes, just not valid quotients",
  "control_perturbs_same_parameter": True, "shared_inputs": ["involution", "kuhn_triangulation"]},
 {"parameter_inserted": "j = number of resolved singular points, every integer 0..k", "selecting_condition": "chi(X_j) == 24 (target_chi_K3)", "condition_uses_true_value": True,
  "solution_set": scan["N=4"]["sol_chi_eq_target"], "classification": "NORMALISATION", "input_it_depends_on": "target_chi_K3; chi(X_j)=chi(Q)+j(chi(S2)-chi(cone)) is affine, so any target picks a j",
  "negative_control": f"j=k-1 gives chi={scan['N=4']['rows'][-2]['chi']} (computed)", "control_perturbs_same_parameter": True, "shared_inputs": ["target_chi_K3", "resolution_local_model"]},
 {"parameter_inserted": "j", "selecting_condition": "b1(X_j;Z/3)==0", "condition_uses_true_value": False, "solution_set": scan["N=4"]["sol_b1_zero"], "classification": "NON_DISCRIMINATING", "negative_control": "none possible: every j passes", "control_perturbs_same_parameter": True, "shared_inputs": ["resolution_local_model"]},
 {"parameter_inserted": "j", "selecting_condition": "Poincare duality b_n=b_{4-n} over Z/3", "condition_uses_true_value": False, "solution_set": scan["N=4"]["sol_poincare_duality_Z3"], "classification": "NON_DISCRIMINATING", "negative_control": "none possible: every j passes (links are Z/3-homology spheres)", "control_perturbs_same_parameter": True, "shared_inputs": ["resolution_local_model"]},
 {"parameter_inserted": "j", "selecting_condition": "every singular point has a Z/2-homology-sphere link, i.e. no cone point (link Betti mod 2 = (1,1,1,1) computed at each unresolved point)", "condition_uses_true_value": False,
  "solution_set": scan["N=4"]["sol_no_unresolved_point_with_nonsphere_Z2_link"], "classification": "RIGID_GIVEN_DEFINITION",
  "input_it_depends_on": "definition of X_j (k-j points left as cones); resolved pieces assumed manifolds (resolution_local_model). Selects j=k by construction, the value chi=24 is then a computed consequence, not selected",
  "negative_control": "j=k-1: one cone point whose Z/2 link Betti (1,1,1,1) differs from S^3 (computed in 00_premise_results.json)", "control_perturbs_same_parameter": True, "shared_inputs": ["resolution_local_model"]},
 {"parameter_inserted": "none (chi(Q) orbifold formula)", "selecting_condition": "chi(Q)=(chi(T^4)+#fixed)/2 with chi(T^4) and #fixed computed", "condition_uses_true_value": False, "solution_set": "holds at every N with regular quotient", "classification": "VERIFIED_IDENTITY", "negative_control": "odd N: regularity fails", "control_perturbs_same_parameter": False, "shared_inputs": ["involution"]},
]
R = lambda id, q, c, sh, s, a="": {"id": id, "quantity": q, "computed": c, "shared_inputs": sh, "script": s, "command": CMD(s, a)}
results = [
 R("D1", "premise (regularity, open-star disjointness, link homology) for N=3..8", {f"N={n}": {k: prem[f"N={n}"][k] for k in ("n_fixed_points", "quotient_is_regular_f_halving", "open_stars_pairwise_disjoint", "closed_stars_pairwise_disjoint(v2_criterion)", "singular_link_representative", "singular_link_is_RP3_homology", "all_singular_links_agree_per_field", "ordinary_link_is_S3", "chi_Q_from_f", "PREMISE_VALID_FOR_MV")} for n in (3, 4, 5, 6, 7, 8)}, ["involution", "kuhn_triangulation"], "00_premise.py", "3 4 5 6 7 8"),
 R("D2", "Betti numbers b0..b4 and chi of the resolved K3 (Mayer-Vietoris on computed pieces; b3, b4 include the chain-level rank of the H_3 map)", {"betti": b, "chi_from_betti": chi_K3, "all_N_fields_agree": len(allB) == 1 and len(allchi) == 1, "per_N": per_N}, ["involution", "kuhn_triangulation", "resolution_local_model", "coefficient_fields"], "02_mv_resolution.py", "4 6 8 --chain-max-N 8"),
 R("D3", "chi from f-vector (chi(Q)+k(chi(S2)-chi(cone))) equals chi from Betti numbers", {"per_N": {n: v["chi_agree"] for n, v in per_N.items()}, "chi_Q": {n: v["chi_Q"] for n, v in per_N.items()}}, ["resolution_local_model"], "02_mv_resolution.py", "4 6 8 --chain-max-N 8"),
 R("D4", "n_singular (fixed vertices of the involution, computed)", sorted(n_sing), ["involution"], "00_premise.py", "3 4 5 6 7 8"),
 R("D5", "K3 x T^2 Betti numbers and chi by Kunneth", kun, ["resolution_local_model", "T2_betti_source"], "04_kunneth_scan_results.py"),
 R("D6", "j-scan over the number of resolved points (every integer 0..k)", scan, ["resolution_local_model", "target_chi_K3"], "04_kunneth_scan_results.py"),
 R("D7", "mod-2 cup pairing on H^2 (parity attempt), T^4 control and singular quotient; NOT the K3 parity", {k: cup[k] for k in ("N", "T4", "Q", "scope")}, ["involution", "kuhn_triangulation"], "03_cup_parity.py", "4"),
]
could = [
 "Parity of the K3 intersection form: no triangulation of the resolved K3 exists here, so the chain-level cup pairing was computed only on T^4 (control, even, nondegenerate rank 6) and on the singular quotient Q (rank 0 mod 2, q identically 0). This does not decide the K3 parity; no odd-form control complex (e.g. CP^2) was built.",
 "Mayer-Vietoris over Z/2 for the resolved space: needs the ranks of phi_1, phi_2 (link H_1, H_2 nonzero mod 2) and the maps into the disc bundle; not computed. Only Z/3 and Z/5 used.",
 "The homology of the resolved piece (disc bundle over S^2) and its gluing to the RP^3 link are a declared tier-L input, not derived from a triangulation.",
 "Torsion / integral homology of the resolved space not computed (fields only).",
 "No realizable negative control for phi_3 (a wrong rank gives b4 != 1, an impossible closed connected manifold, which is not an object to test).",
]
(HERE / "inputs.json").write_text(json.dumps(inputs, indent=1))
(HERE / "results.json").write_text(json.dumps({"track": "D-tda v3", "inputs_file": "inputs.json", "results": results, "rigidity": rig, "could_not_do": could}, indent=1, default=str))
(HERE / "exports.json").write_text(json.dumps({"chi": chi_K3, "b2": b[2], "n_singular": sorted(n_sing)[0], "betti_b0_b4": b, "chi_K3xT2": kun["chi_from_product_betti"],
   "agree_across_N_4_6_8_and_Z3_Z5": len(allB) == 1 and len(allchi) == 1 and len(n_sing) == 1, "generated_by": "04_kunneth_scan_results.py",
   "status": "computed by exact Mayer-Vietoris on GUDHI Betti numbers; depends on declared input resolution_local_model (tier L)"}, indent=1))
print(json.dumps({"betti": b, "chi": chi_K3, "kun": kun["b_K3xT2"], "chiKxT": kun["chi_from_product_betti"], "n_sing": sorted(n_sing), "scanN4_sol": {k: v for k, v in scan["N=4"].items() if k.startswith("sol")}, "slope": scan["N=4"]["slope"], "intercept": scan["N=4"]["intercept"]}))
