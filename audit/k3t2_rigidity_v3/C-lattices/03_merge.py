"""Track C v3 part 3: merge into results.json + exports.json (script-generated).
Order: cd audit/k3t2_rigidity_v3/C-lattices && PY decl.py && PY 01_lattices.py && PY 02_tduality.py && PY 03_merge.py"""
import json
from decl import HERE
from k3chain import chain
L = json.loads((HERE / "01_lattices.json").read_text()); T = json.loads((HERE / "02_tduality.json").read_text())
ch = chain()
CH_IN = ["K_trivial", "hodge_h0_O", "hodge_h1_O_from_b1", "serre_duality", "noether_formula"]
SIG_IN = CH_IN + ["hirzebruch_signature", "hodge_index", "betti_odd_vanish"]
sel = L["selected"]; cands = L["candidates_mU_nE8_rank_b2"]
tdok = all(all(v.values()) for v in T["by_d"].values()) and T["neg_control_K"]["broken_K_squares_to_identity"] is False
entok = all(v["identity_difference_zero"] for v in T["entry_level_bound"].values())
results = [
 {"id": "C1_chi_top", "quantity": "chi_top(K3) from K=O, b1=0, Serre duality, Noether",
  "computed": f"chi(O)={ch['chi_O']}, chi_top={ch['chi_top']} (verified by exact computation, tier B given tier-L inputs)",
  "shared_inputs": CH_IN, "script": "audit/k3t2_rigidity_v3/C-lattices/k3chain.py",
  "command": "cd audit/k3t2_rigidity_v3/C-lattices && PY 03_merge.py",
  "note": "chi_top = (declared Noether constant 12) x (computed chi(O)=2): not independent of noether_formula"},
 {"id": "C2_signature", "quantity": "(b2+, b2-) from Hirzebruch + Hodge index; cross-check by Hodge diamond",
  "computed": f"tau={ch['tau']} (diamond route {ch['tau_diamond']}), b2={ch['b2']}, (b2+,b2-)=({ch['b2plus']},{ch['b2minus']})",
  "shared_inputs": SIG_IN, "script": "audit/k3t2_rigidity_v3/C-lattices/k3chain.py", "command": "as C1",
  "note": "diamond route shares chi(O)=2, so the two tau routes are not independent of each other's inputs"},
 {"id": "C3_lattice_selection", "quantity": "which mU+n(-E8) has signature (b2+,b2-)",
  "computed": f"rank-{ch['b2']} candidates (m,n)={[ (c['m'],c['n'],tuple(c['signature_ldl'])) for c in cands]}; selected={[(c['m'],c['n']) for c in sel]}; scan 0<=m,n<=30 all integers gives {L['all_mn_0_30_with_signature_equal_hodge']}",
  "shared_inputs": SIG_IN + ["k3_lattice_identification", "U_gram", "e8_root_system"],
  "script": "audit/k3t2_rigidity_v3/C-lattices/01_lattices.py", "command": "cd audit/k3t2_rigidity_v3/C-lattices && PY 01_lattices.py",
  "note": "E8 Gram computed from root system (det 1, 240 roots, even). Controls (11,0),(7,1) are existing even unimodular lattices with other signatures"},
 {"id": "C4_typed_gram_consistency", "quantity": "signature/rank of typed 3U+2(-E8)",
  "computed": f"rank {L['typed_gram_3U_2mE8']['rank']}, signature {L['typed_gram_3U_2mE8']['signature_ldl']}; equals Hodge chain: {tuple(L['typed_gram_3U_2mE8']['signature_ldl'])==(ch['b2plus'],ch['b2minus'])}",
  "shared_inputs": ["typed_gram_3U_2mE8", "U_gram", "e8_root_system"], "script": "audit/k3t2_rigidity_v3/C-lattices/01_lattices.py",
  "command": "as C3", "note": "NOT an independent route: rank and signature follow from the typed Gram"},
 {"id": "C5_mukai_rank", "quantity": "Mukai lattice rank vs chi_top", "computed": f"rank {L['mukai']['rank']}, signature {L['mukai']['signature']}",
  "shared_inputs": ["typed_gram_3U_2mE8", "betti_odd_vanish"], "script": "01_lattices.py", "command": "as C3",
  "note": "consistency only (Gram typed as 4U+2(-E8))"},
 {"id": "C6_tduality_identities", "quantity": "K_i involution, P(A)^T eta P(A)=eta, eta H eta = H(G^-1), d=1..4 generic symbols, plus lambda in R redo",
  "computed": f"all hold: {tdok}; lambda scaling preserves eta and maps H(G)->H(lambda^2 G) for symbolic real lambda; integral lambda among reduced p/q, p,q<=200: {T['lambda_in_R']['integral_lambda_among_reduced_p/q_p,q<=200']}; f=G+1/G invariant under swap: {T['lambda_in_R']['f_invariant_under_swap']}, under lambda scaling: {T['lambda_in_R']['f_invariant_under_lambda_scaling']}",
  "shared_inputs": ["eta_definition"], "script": "audit/k3t2_rigidity_v3/C-lattices/02_tduality.py", "command": "cd audit/k3t2_rigidity_v3/C-lattices && PY 02_tduality.py"},
 {"id": "C7_dual_scale_entry_level", "quantity": "tr G + tr G^-1 - 2d = ||L - L^-T||_F^2 for G=LL^T, L lower-triangular with symbolic entries, d=1..4",
  "computed": f"identity difference zero for d=1..4: {entok}; broken control ||L-L^-1||^2 nonzero for d=2,3,4 (coincides at d=1); G^-1 computed directly from G (LU) for all d; implies bound >=2d, equality iff G=I",
  "shared_inputs": ["eta_definition"], "script": "audit/k3t2_rigidity_v3/C-lattices/02_tduality.py", "command": "as C6"},
]
rig = [
 {"parameter_inserted": "b2+ (equivalently signature) selecting m in mU+n(-E8)", "selecting_condition": "signature (m,m+8n) equals Hodge-index signature from chi_top via Hirzebruch; rank 22 from Noether+Betti",
  "condition_uses_true_value": False, "solution_set": f"{L['all_mn_0_30_with_signature_equal_hodge']} among all integer (m,n) in 0..30; other rank-22 candidates (11,0) sig (11,11), (7,1) sig (7,15) fail",
  "classification": "RIGID_GIVEN_DEFINITION", "input_it_depends_on": "K=O, b1=0, Noether, Hirzebruch, Hodge index, lattice family (declared)",
  "negative_control": "perturb (m,n): every neighbouring existing even unimodular lattice has wrong signature", "control_perturbs_same_parameter": True},
 {"parameter_inserted": "chi_top", "selecting_condition": "chi(O) from K=O, b1=0, Serre duality; Noether with c1^2=0",
  "condition_uses_true_value": False, "solution_set": "chi_top=24 unique given the declared inputs (it is a computation, not a scan)",
  "classification": "RIGID_GIVEN_DEFINITION", "input_it_depends_on": "K_trivial (definition of K3) and noether_formula",
  "negative_control": "none realizable: a surface with chi(O)=3 has K != O, hence is not K3 (not a control); labelled HYPOTHETICAL and not counted",
  "control_perturbs_same_parameter": False},
 {"parameter_inserted": "typed Gram 3U+2(-E8) rank/signature", "selecting_condition": "none (typed)", "condition_uses_true_value": True,
  "solution_set": "n/a", "classification": "CONDITIONAL_ON_INPUT", "input_it_depends_on": "typed_gram_3U_2mE8", "negative_control": "none", "control_perturbs_same_parameter": False},
 {"parameter_inserted": "T-duality / generalized-metric identities d=1..4, lambda in R", "selecting_condition": "universal identities, nothing selected",
  "condition_uses_true_value": False, "solution_set": "all", "classification": "VERIFIED_IDENTITY", "input_it_depends_on": "eta_definition",
  "negative_control": "broken K_i (no diagonal zeroing) fails to square to I; lambda scaling does not preserve f=G+1/G", "control_perturbs_same_parameter": True},
 {"parameter_inserted": "dual-scale bound constant 2d", "selecting_condition": "sum-of-squares identity ||L-L^-T||^2 (equality iff G=I)",
  "condition_uses_true_value": False, "solution_set": "bound tr G + tr G^-1 >= 2d, tight only at G=I", "classification": "VERIFIED_IDENTITY",
  "input_it_depends_on": "none", "negative_control": "||L-L^-1||^2 differs from the true difference for d>=2", "control_perturbs_same_parameter": True},
]
inputs = json.loads((HERE / "inputs.json").read_text())
out = {"track": "C", "version": "v3", "inputs_file": "audit/k3t2_rigidity_v3/C-lattices/inputs.json", "results": results, "rigidity": rig,
       "could_not_do": ["T-duality group O(d,d;Z) generation/discreteness only checked for diagonal lambda scalings (p,q<=200), not the full group",
                        "no independent route to lattice rank/signature beyond the Hodge chain; typed Gram is a declared input",
                        "uniqueness of E8 as even unimodular rank 8 not re-derived (literature)",
                        "no realizable negative control for chi_top: any other chi(O) leaves the K3 class"]}
(HERE / "results.json").write_text(json.dumps(out, indent=2, default=str))
exp = {"chi_top": {"value": ch["chi_top"], "shared_inputs": CH_IN, "path": "audit/k3t2_rigidity_v3/C-lattices/results.json", "label": "RIGID_GIVEN_DEFINITION"},
       "signature": {"value": [ch["b2plus"], ch["b2minus"]], "b2": ch["b2"], "shared_inputs": SIG_IN, "path": "audit/k3t2_rigidity_v3/C-lattices/results.json", "label": "RIGID_GIVEN_DEFINITION"}}
(HERE / "exports.json").write_text(json.dumps(exp, indent=2))
print(json.dumps(exp, indent=1)); print(tdok, entok)
