"""Build agreement_ledger.json for the k3t2-rigidity loop (report stage).

Reads the committed comparator output, the sealed LeanMaster v3.13.1 targets,
the skeptic outputs and the reverse-pass results. Counts are COMPUTED from the
per-row classification below, never typed. The skeptic's per-row judgements
(which rows are disputed, which scripts carry literals) are copied from the
skeptic's audit (workflow input AUDIT, summarised in skeptic/audit_verdicts.json).

For three sealed targets that the comparator never rowed, this script does a
small exact side-check (Hurwitz values, Mukai signature, E8 leading minors).
These are report-writer checks, NOT comparator verdicts, and are NOT counted.

Run:
  /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python \
      audit/k3t2_rigidity/build_ledger.py
"""
import json
import re
from fractions import Fraction as Fr
from pathlib import Path

R = Path(__file__).resolve().parent
PY = "/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python"


def load(p):
    return json.loads((R / p).read_text())


comparison = load("comparison.json")
sealed = load("targets_sealed.json")
skeptic = load("skeptic/audit_verdicts.json")
s2 = load("skeptic/s2_polar_joint_scan_results.json")
p12 = load("reverse/p1_p2_results.json")
p34 = load("reverse/p3_p4_results.json")
p5 = load("reverse/p5_trace_bound_results.json")

sealed_by_id = {r["id"]: r for rows in sealed["tracks"].values() for r in rows}

DYONS_PREP = f"cd audit/k3t2_rigidity/dyons && {PY} theta_forms.py 22"  # QMAX=22 is required (skeptic reproducibility_gap)

# Per-row classification after the skeptic. independence:
#   clean     = no skeptic dispute and no hard-coded literal on the path to the reported value
#   qualified = AGREE verdict kept by the skeptic, but independence/route count overstated or a literal on the path
ROWS = {
    "genus-moonshine_lock_at_identity": dict(
        independence="qualified",
        skeptic_dispute="AGREE stands; factor-2 convention checked (90=2*45, 462=2*231). Qualifier: chi(2A)=8 depends on the tier-L polar-term premise, and H_2A is computed with a literal 8 (run_track_a.py:437), not the solved value.",
        hardcoded_violation="run_track_a.py:312 extract_H(24, ...) uses literal 24 on the mu-term (justified only afterwards by scan (a) at :376); run_track_a.py:437 literal Fr(8,24) for chi(2A).",
        script="genus-moonshine/run_track_a.py -> genus-moonshine/results.json",
        command=f"cd audit/k3t2_rigidity/genus-moonshine && {PY} run_track_a.py"),
    "genus-moonshine_norm_3A": dict(independence=None, skeptic_dispute=None, hardcoded_violation=None,
        script=None, command=None, scope_note="no blind 3A twining"),
    "genus-moonshine_umbral_exact": dict(independence=None, skeptic_dispute=None, hardcoded_violation=None,
        script=None, command=None, scope_note="blind built only Z^(2)=Z_K3"),
    "dyons_c_depends_only_on_D": dict(
        independence="clean", skeptic_dispute=None, hardcoded_violation=None,
        script="genus-moonshine/run_track_a.py; dyons/theta_forms.py (+ dyons/theta_forms_cache.json)",
        command=f"{DYONS_PREP}; cd ../genus-moonshine && {PY} run_track_a.py",
        scope_note="Tracks A and B compute the same theta/phi_{0,1} quantity (skeptic: all 18 shared c(D) identical = common provenance), so this is one blind route against Lean, not two."),
    "dyons_c_first": dict(
        independence="clean", skeptic_dispute=None, hardcoded_violation=None,
        script="genus-moonshine/run_track_a.py; dyons/theta_forms.py",
        command=f"{DYONS_PREP}; cd ../genus-moonshine && {PY} run_track_a.py",
        scope_note="Same common-provenance remark as dyons_c_depends_only_on_D. Normalisation Z_K3 = 2*phi_{0,1} is a definition input."),
    "dyons_p24_values": dict(
        independence="qualified",
        skeptic_dispute="AGREE stands, but 'via DMVV AND via prod(1-p^n)^{-24}' overstates independence: the q^0,y=1 restriction of DMVV collapses analytically to prod(1-p^r)^{-(c0+2c-1)}, and the comparison side is typed with the literal 24 (part1_euler_numbers.py:81). One route checked against a target built from the answer.",
        hardcoded_violation="dyons/part1_euler_numbers.py:81 comb(24+j-1, j) builds the Goettsche table from literal 24.",
        script="dyons/part1_euler_numbers.py -> dyons/part1_euler_numbers_results.json",
        command=f"{DYONS_PREP} && {PY} part1_euler_numbers.py"),
    "lattices-duality_sigK3_eq": dict(
        independence="clean", skeptic_dispute=None, hardcoded_violation=None,
        script="lattices-duality/02_k3_mukai_gamma.py -> 02_k3_mukai_gamma_result.json",
        command=f"cd audit/k3t2_rigidity/lattices-duality && {PY} 02_k3_mukai_gamma.py",
        scope_note="Lean side is pair-addition bookkeeping; blind side computes from the full Gram matrix. The script's after-the-fact 'expected' field says det=1 for K3; the computed value is -1 (recalled expectation wrong, never compared, harmless)."),
    "lattices-duality_sigK3T2_eq": dict(
        independence="clean", skeptic_dispute=None, hardcoded_violation=None,
        script="lattices-duality/02_k3_mukai_gamma.py -> 02_k3_mukai_gamma_result.json",
        command=f"cd audit/k3t2_rigidity/lattices-duality && {PY} 02_k3_mukai_gamma.py"),
    "lattices-duality_thetaShift_isODD": dict(
        independence="clean", skeptic_dispute=None, hardcoded_violation=None,
        script="lattices-duality/03_oddz_checks.py -> 03_oddz_checks_result.json",
        command=f"cd audit/k3t2_rigidity/lattices-duality && {PY} 03_oddz_checks.py",
        scope_note="Blind: 200 samples each at d=2,3. Lean is universal in d and Theta and strictly stronger."),
    "tda-gudhi_k3_euler_characteristic": dict(
        independence="qualified",
        skeptic_dispute="AGREE can stand, but 'four blind routes reaching 24 by different means' is wrong. Track C is arithmetic on typed BETTI_K3=[1,0,22,0,1] (05_tadpole_arithmetic.py:47). Track A Z_K3(tau,0)=24 and Track B c(0)+2c(-1)=24 are the same theta/phi_{0,1} quantity. Track D is invalid at N=4 (Mayer-Vietoris premise fails, U n B = (1,105,0,16)) and valid only at N=6, where it is hybrid with a typed 16. Realistically: one elliptic-genus computation plus one valid GUDHI data point (N=6). 'Stable at N=4,6' withdrawn.",
        hardcoded_violation="lattices-duality/05_tadpole_arithmetic.py:47 typed Betti; tda-gudhi/04_signature.py:102 num_exceptional=16 typed; tda-gudhi/03_resolution_hybrid.py:153,:175 literals 16 and 32; tda-gudhi/07_rigidity_scan.py:45 target literal 24.",
        script="genus-moonshine/run_track_a.py; tda-gudhi/03_resolution_hybrid.py (N=6 column only)",
        command=f"cd audit/k3t2_rigidity/tda-gudhi && {PY} 03_resolution_hybrid.py"),
    "tda-gudhi_kummer_mapper_betti1_eq_376": dict(independence=None, skeptic_dispute=None, hardcoded_violation=None,
        script=None, command=None, scope_note="simulator Mapper constants, not K3xT2 mathematics; out of scope"),
    "tda-gudhi_kummer_mapper_euler_eq_minus_370": dict(independence=None, skeptic_dispute=None, hardcoded_violation=None,
        script=None, command=None, scope_note="simulator Mapper constants, not K3xT2 mathematics; out of scope"),
}

ledger_rows = []
for r in comparison["rows"]:
    tid = r["target_id"]
    c = ROWS[tid]
    st = r["status"]
    if st == "AGREE":
        after = "AGREE_CLEAN" if c["independence"] == "clean" else "AGREE_QUALIFIED"
    else:
        after = st
    ledger_rows.append({
        "target_id": tid,
        "lean_file_line": sealed_by_id[tid]["file_line"],
        "lean_value": r["lean_value"],
        "blind_value": r["blind_value"],
        "status_comparator": st,
        "status_after_skeptic": after,
        "independence": c["independence"],
        "skeptic_dispute": c["skeptic_dispute"],
        "hardcoded_violation": c["hardcoded_violation"],
        "script": c["script"],
        "command": c["command"],
        "comparator_note": r["note"],
        "scope_note": c.get("scope_note"),
    })

# Sealed targets never given a comparison row.
rowed = {r["target_id"] for r in comparison["rows"]}
unrowed = []
for rows in sealed["tracks"].values():
    for t in rows:
        if t["id"] in rowed:
            continue
        unrowed.append({"target_id": t["id"], "lean_file_line": t["file_line"],
                        "status_after_skeptic": "NOT_ATTEMPTED",
                        "statement_head": t["statement_verbatim"].split(":=")[0].strip()[:300]})

# Report-writer side checks (exact; not verdicts, not counted).
side = []
# (1) Hurwitz: Lean list parsed from the sealed statement, blind H(D) from the dyons JSON.
st = sealed_by_id["dyons_hurwitz_values"]["statement_verbatim"]
m = re.search(r"\[([-\d,\s]+)\]\.map h12 = \[([-\d,\s]+)\]", st)
Ds = [int(x) for x in m.group(1).split(",")]
lean12 = [int(x) for x in m.group(2).split(",")]
hz = load("dyons/hurwitz_class_numbers_results.json")
blind12 = [12 * Fr(hz[str(D)]["H"]) for D in Ds]
side.append({"target_id": "dyons_hurwitz_values",
             "lean": dict(zip(map(str, Ds), lean12)),
             "blind_12H": {str(D): str(v) for D, v in zip(Ds, blind12)},
             "match_all": all(Fr(a) == b for a, b in zip(lean12, blind12)),
             "caveat": "D=0 blind value -1/12 is 'by convention, not from the counting method' (per its JSON); the counted cases are D=3..15."})
# (2) Mukai signature.
st = sealed_by_id["lattices-duality_sigMukai_eq"]["statement_verbatim"]
m = re.search(r"⟨(\d+),\s*(\d+)⟩", st)
lean_sig = (int(m.group(1)), int(m.group(2)))
mk = load("lattices-duality/02_k3_mukai_gamma_result.json")["results"]["Mukai_K3_plus_U"]
sa, sb = mk["signature_method_A_congruence_diagonalization"], mk["signature_method_B_sturm_charpoly"]
side.append({"target_id": "lattices-duality_sigMukai_eq", "lean": list(lean_sig),
             "blind_methodA": [sa["p"], sa["q"]], "blind_methodB": [sb["p"], sb["q"]],
             "match_all": (sa["p"], sa["q"]) == lean_sig == (sb["p"], sb["q"])})
# (3) E8 Cartan positive definite <=> all leading principal minors > 0 (Sylvester).
e8 = load("lattices-duality/01_e8_result.json")
minors = [Fr(x) for x in e8["leading_principal_minors"]]
side.append({"target_id": "lattices-duality_cartanE8_posDef", "lean": "cartanE8R.PosDef",
             "blind_leading_minors": [str(x) for x in minors],
             "match_all": all(x > 0 for x in minors),
             "caveat": "Sylvester criterion on the blind integer Cartan matrix; assumes it is the same matrix up to labelling as Lean's cartanE8R (not checked here)."})

# Counts computed from rows.
def count(key, val):
    return sum(1 for x in ledger_rows if x[key] == val)

counts = {
    "sealed_targets_total": len(sealed_by_id),
    "rowed_by_comparator": len(ledger_rows),
    "not_attempted_no_row": len(unrowed),
    "comparator": {s: count("status_comparator", s) for s in ["AGREE", "DISAGREE", "NOT_COMPUTED"]},
    "after_skeptic": {s: count("status_after_skeptic", s)
                      for s in ["AGREE_CLEAN", "AGREE_QUALIFIED", "DISAGREE", "NOT_COMPUTED"]},
    "not_computed_out_of_scope_simulator_constants": sum(
        1 for x in ledger_rows if x["status_after_skeptic"] == "NOT_COMPUTED"
        and "simulator" in (x["scope_note"] or "")),
    "side_checks_matching_not_counted": sum(1 for s in side if s["match_all"]),
}
assert counts["comparator"]["AGREE"] == comparison["agree"]
assert counts["comparator"]["NOT_COMPUTED"] == comparison["not_computed"]

RIGIDITY = [
    ("A", "N (coefficient of the Appell-Lerch mu-term)", "20..28", "{24}", True,
     "Structural condition (j=0 and j=2 y-slices give the same H); N=23,25 give polar terms -3/2,-5/2 vs -2. One linear condition; recovers N = Z_K3(tau,0) from the same elliptic genus, so not independent of the k-scan's 24."),
    ("A", "k (overall factor in Z_K3 = k*phi_{0,1})", "{2,1,3,5/2,2+1/12,2+1/100,0}", "{2}", False,
     "Pinned only by Z(tau,0) == literal 24 (run_track_a.py:114). Integrality allows every integer k."),
    ("A", "chi(2A)", "0..40", "integrality: {8,20,32,...}; + tier-L polar premise -2: {8}", False,
     "Uniqueness needs the tier-L premise fed in as Fr(-2) (:426); downstream H_2A uses literal 8 (:437)."),
    ("B", "kappa (DMVV exponent multiplier)", "n/12, n=0..36", "{1}", False,
     "Target built with literal 24 (part1_euler_numbers.py:81); kappa*24=24 is one linear equation against the answer."),
    ("B", "c(3) or c(4) shifted by +1", "{0,+1} per coefficient", "{0}", True,
     "Deforms a real input; DMZ m=1 identity breaks at 122 coefficients. Only +1 tested; the identity's 9 and 3 are tier L."),
    ("B", "(N,M) polar subtraction", "N 320..328, M 640..656", "{(324,648)}", True,
     "Blind scans were coordinate-wise at literal true values (:205, :250) and on a window that hid a bug (part3_polar.py:141 drops q^{s^2+s}). Skeptic's joint 9x17 grid with corrected A_{2,1}: unique (324,648) over the full window, no D filter. Answer right, blind method wrong."),
    ("C", "(m,n) in m*U + n*(-E8), rank 22", "(11,0),(7,1),(3,2)", "intrinsic: 3 solutions; + signature (3,19): {(3,2)}", False,
     "Only selecting condition is typed signature (3,19) (06_rigidity_a_lattice_enum.py:79), which restates the answer since sig = (m, m+8n)."),
    ("C", "lambda in R -> lambda*alpha'/R", "12 rationals p/q", "{-1,+1}; + positivity: {+1}", False,
     "Involution holds for every lambda; spectrum preservation gives lambda^2=1; +1 needs an added positivity input."),
    ("D", "k (resolved singular points)", "0..16", "{16}", False,
     "chi = 8+k against literal 24 (07_rigidity_scan.py:45); smoothness criterion is definitional; N=4 column rests on an invalid Mayer-Vietoris premise."),
    ("D", "group action (negation vs translation)", "2 actions", "negation -> (1,0,6,0,1); translation -> T^4 at N=6", False,
     "A control, not a scan. N=4 translation control is not a triangulation (1856 vs 1920 edges); N=6 control valid."),
    ("D", "N (triangulation refinement)", "{4,6}", "invariance claimed", False,
     "Invariance check, not a pin; only N=6 is a valid data point."),
]
rigidity = [{"track": t, "parameter": p, "scanned": s, "solution_set": sol,
             "skeptic_genuinely_rigid": g, "skeptic_reason": why}
            for t, p, s, sol, g, why in RIGIDITY]

p1 = p12["P1_goettsche_k6_7_8"]["method_A_direct_partition_convolution"]
reverse = [
    {"id": "P1", "extends": "DualScaleDyons/DMVV.lean p24_values, goettsche (k<=5)",
     "computed": {k: p1[k] for k in ["6", "7", "8"]}, "holds": True,
     "lean_reachable_now": True,
     "caveat": "One route plus a structural derivation, not two numerical routes; exp_correct=24 literal in reverse/p1_p2_dmz516_m4.py:120 (admitted at :91-96). Negative control at exponent 23/25 fails."},
    {"id": "P2", "extends": "DMVV.lean dmz_516_q1/q2, m=4 row (q^1 only)",
     "computed": {k: p12["P2_dmz516_m4"][k]["holds"] for k in ["q^1", "q^2", "q^3", "q^4"]}, "holds": True,
     "lean_reachable_now": False,
     "caveat": "The five RHS coefficients (51,155,93,102,31) are a tier-L ansatz from arXiv:1208.4074 via Lean's docstring. Needs zK3 raised to ellipticGenus N>=10 (q^2)."},
    {"id": "P3", "extends": "polar_part_removes_pole / polar_coefficient_pinned / immortal_exact_m23 (m<=3)",
     "computed": "vanish2 at p24(5)=176256 holds; fails at 176255, 176257", "holds": True,
     "lean_reachable_now": True,
     "caveat": "d1 is identically 0 in every case incl. negative controls; only eval1 discriminates."},
    {"id": "P4", "extends": "same, m=5",
     "computed": "vanish2 at p24(6)=1073720 holds; fails at 1073719, 1073721", "holds": True,
     "lean_reachable_now": False,
     "caveat": "Same d1 caveat. Needs polarDefect re-parametrised from literal K=4 and zK3 raised to N>=12."},
    {"id": "P5", "extends": "blind Track C 04_dual_scale_bound.py (d<=6); NOT a located Lean statement",
     "computed": {d: p5["by_dimension"][d]["random_spd_sample_bound"]["all_f_nonnegative"] for d in ["7", "8", "9", "10"]},
     "holds": True, "lean_reachable_now": None,
     "caveat": "Structurally an algebraic identity (sum (x_i-1)^2/x_i). Note: targets_sealed.json lists DualScaleStream2/DualScale/TraceBound.lean:87 dualScale_eq and :173 dualScale_eq_iff, so the reverse pass's 'no Lean statement located' (it searched StringTheoryFoundation/ and DoubleFieldTheory/) is a search gap, not an absence."},
]

out = {
    "header": "Generated by workflow k3t2-rigidity-loop; computations are tier B at best; physical identifications tier L/C.",
    "leanmaster_tag_sealed": sealed["metadata"]["leanmaster_tag"],
    "leanmaster_commit_comparator": comparison["metadata"]["leanmaster_commit"],
    "reverse_pass_leanmaster_read_at": "7be7626 (per reverse/REVERSE_PASS_FINDINGS.md; differs from the sealed tag)",
    "branch": "loop/k3t2-rigidity",
    "counts": counts,
    "rows": ledger_rows,
    "not_attempted": unrowed,
    "report_writer_side_checks_not_counted": side,
    "rigidity": rigidity,
    "rigidity_genuinely_rigid": sum(1 for x in rigidity if x["skeptic_genuinely_rigid"]),
    "rigidity_total": len(rigidity),
    "reverse": reverse,
    "resolved_discrepancy": {
        "was": "OPEN DISCREPANCY -324 at (n,l)=(4,-5), D=-9 in Track B polar subtraction",
        "now": "code bug at dyons/part3_polar.py:141 (s<=-1 branch of A_{2,1} drops q^{s^2+s}); corrected joint scan unique (324,648) on full window",
        "evidence": {"A21_part3_vs_corrected": s2["A21_part3_vs_corrected_differences_in_window"],
                     "corrected_full_window_joint_solutions": s2["corrected_A21__full_window_no_D_filter"]["joint_solutions"],
                     "part3_full_window_joint_solutions": s2["part3_A21__full_window_no_D_filter"]["joint_solutions"]},
        "fix_applied_to_part3_polar_py": False,
    },
    "blindness": skeptic["blindness"],
    "reproducibility_gap": skeptic["reproducibility_gap"],
    "tda_valid": False,
}
(R / "agreement_ledger.json").write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n")
print(json.dumps(counts, indent=1))
print("side checks:", [(s["target_id"], s["match_all"]) for s in side])
print("rigid:", out["rigidity_genuinely_rigid"], "/", out["rigidity_total"])
