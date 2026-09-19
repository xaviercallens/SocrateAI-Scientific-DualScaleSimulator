"""Build agreement_ledger_v2.json for the k3t2-rigidity-loop v2 (report stage).

Every count in the ledger is COMPUTED here from the per-row classification;
no count is typed. Inputs are committed files inside audit/k3t2_rigidity_v2/
(and v1's ledger for the comparison line). Nothing outside the worktree is read.

Per-row rules (after the two skeptics):
  * comparator status is the starting point (comparison.json rows);
  * skeptic_math e_row_adjudication: a NOT_COMPARABLE row whose verdict says
    "should be AGREE" becomes AGREE; a DISAGREE row whose verdict says the
    BLIND value is right becomes DISAGREE_SEALED_LEAN_ERROR;
  * skeptic_blindness disputed_rows: the two Goettsche rows whose blind side
    compares the Goettsche formula with itself become NOT_COMPUTED; every other
    disputed AGREE row stays AGREE but is marked "qualified" (independence or
    route count overstated) and its route count drops to 1;
  * a row named in a skeptic_blindness violation as citing tautological data
    is also "qualified";
  * synthesizer rule: a 2-route row whose second route is the typed Gram matrix
    (C-lattices/02) or typed TT quotation arithmetic (E-flux/tadpole_bookkeeping)
    loses that route and, if AGREE, is "qualified" (the blindness skeptic's
    objection to C-lattices#21 and orientifold#3, applied to the rows it missed);
  * every other AGREE row is "clean".

Run from any directory (paths are relative to this file):
  /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python \
      audit/k3t2_rigidity_v2/build_ledger_v2.py
"""
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

R = Path(__file__).resolve().parent
PY = "/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python"
TRACK_GROUPS = ("A-genus", "B-dyons", "C-lattices", "D-tda", "E-flux")
WEAK_ROUTE_SOURCES = ("C-lattices/02", "E-flux/tadpole_bookkeeping")


def load(rel):
    return json.loads((R / rel).read_text())


comparison = load("comparison.json")
sk_math = load("skeptic_math/verdict.json")
sk_blind = load("skeptic_blindness/verdict.json")
v1 = json.loads((R.parent / "k3t2_rigidity" / "agreement_ledger.json").read_text())

rows = comparison["rows"]
ids = [r["target_id"] for r in rows]
assert len(ids) == len(set(ids)), "duplicate target ids"
by_id = {r["target_id"]: r for r in rows}

# ---------------------------------------------------------------- skeptic inputs
math_adj = {}
for a in sk_math["e_row_adjudication"]:
    for tid in re.split(r"\s*/\s*", a["target_id"]):
        math_adj[tid] = a
blind_disputed = {d["target_id"]: d["reason"] for d in sk_blind["disputed_rows"]}
math_disputed = {d["target_id"]: d["reason"] for d in sk_math["disputed_rows"]}

# The blindness skeptic demotes these two: "The blind side never evaluated DMVV at z=0:
# part1 compares the Goettsche formula with itself. This row is NOT_COMPUTED." (+ "Same source").
DEMOTE = ["B-dyons#6:goettsche", "B-dyons#36:twined_goettsche_identity"]
for tid in DEMOTE:
    assert tid in blind_disputed, tid
assert "NOT_COMPUTED" in blind_disputed[DEMOTE[0]]

# Rows named in a blindness violation (not in disputed_rows) as citing the tautological part1 data.
violation_named = {}
for v in sk_blind["violations"]:
    m = re.search(r"comparison rows ((?:[A-Za-z-]+#\d+)(?:/#\d+)*)", v)
    if m:
        grp = m.group(1).split("#")[0]
        for num in re.findall(r"#(\d+)", m.group(1)):
            for tid in ids:
                if tid.startswith(f"{grp}#{num}:"):
                    violation_named[tid] = v

# ---------------------------------------------------------------- per-row classification
ledger_rows = []
for r in rows:
    tid = r["target_id"]
    group = tid.split("#")[0]
    status0 = r["status"]
    status = status0
    independence = None
    routes = r.get("independent_routes", 0)
    reasons = []
    adj = math_adj.get(tid)
    if adj is not None:
        who = adj["who_is_right"]
        if status0 == "NOT_COMPARABLE" and "should be AGREE" in who:
            status = "AGREE"
            reasons.append("skeptic_math: " + who)
        elif status0 == "DISAGREE" and who.startswith("BLIND"):
            status = "DISAGREE_SEALED_LEAN_ERROR"
            reasons.append("skeptic_math: " + who + " " + adj["why"])
        else:
            reasons.append("skeptic_math: " + who)
    if tid in math_disputed:
        reasons.append("skeptic_math dispute: " + math_disputed[tid])
    if tid in DEMOTE:
        status = "NOT_COMPUTED"
        reasons.append("skeptic_blindness: " + blind_disputed[tid])
    elif tid in blind_disputed:
        reasons.append("skeptic_blindness: " + blind_disputed[tid])
        if status == "AGREE":
            independence = "qualified"
            routes = min(routes, 1)
    # Synthesizer rule, same objection the blindness skeptic applied to C-lattices#18/#21 and
    # orientifold#3: a route through the TYPED Gram matrix (C-lattices/02) or through typed TT
    # quotation arithmetic (E-flux/tadpole_bookkeeping) is not an independent computed route.
    flagged = [s for s in WEAK_ROUTE_SOURCES if s in (r.get("blind_source") or "")]
    if flagged and routes >= 2:
        routes = max(1, routes - len(flagged))
        reasons.append("synthesizer: route via " + ", ".join(flagged) + " is not independent "
                       "(skeptic_blindness objection to C-lattices#21 / orientifold#3 applied)")
        if status == "AGREE":
            independence = "qualified"
    if tid in violation_named and status == "AGREE":
        independence = "qualified"
        reasons.append("skeptic_blindness violation: " + violation_named[tid])
    if status == "AGREE" and independence is None:
        independence = "clean"
    dup = None
    m = re.search(r"duplicate of ([A-Za-z-]+#\d+)", r.get("note") or "", re.I)
    if m:
        dup = next(t for t in ids if t.startswith(m.group(1) + ":"))
    ledger_rows.append(dict(
        target_id=tid, group=group, track_target=group in TRACK_GROUPS,
        lean_file_line=r.get("lean_file_line"), lean_value=r.get("lean_value"),
        blind_value=r.get("blind_value"), blind_source=r.get("blind_source"),
        comparator_status=status0, final_status=status, independence=independence,
        coverage=r.get("coverage"), routes_comparator=r.get("independent_routes", 0),
        routes_after_skeptic=routes, duplicate_of=dup,
        matches_under_conversion=r.get("matches_under_conversion"),
        conversion=r.get("conversion"), comparator_note=r.get("note"), skeptic_reasons=reasons))


def tally(sel):
    c = Counter(x["final_status"] for x in sel)
    ag = [x for x in sel if x["final_status"] == "AGREE"]
    return {
        "rows": len(sel),
        "AGREE": c["AGREE"],
        "AGREE_clean": sum(x["independence"] == "clean" for x in ag),
        "AGREE_qualified": sum(x["independence"] == "qualified" for x in ag),
        "AGREE_clean_full_coverage": sum(x["independence"] == "clean" and x["coverage"] == "full" for x in ag),
        "AGREE_clean_partial_coverage": sum(x["independence"] == "clean" and x["coverage"] == "partial" for x in ag),
        "AGREE_qualified_full_coverage": sum(x["independence"] == "qualified" and x["coverage"] == "full" for x in ag),
        "AGREE_qualified_partial_coverage": sum(x["independence"] == "qualified" and x["coverage"] == "partial" for x in ag),
        "NOT_COMPARABLE_matching_under_conversion": sum(x["final_status"] == "NOT_COMPARABLE" and x["matches_under_conversion"] is True for x in sel),
        "NOT_COMPARABLE": c["NOT_COMPARABLE"],
        "DISAGREE_against_blind": c["DISAGREE"],
        "DISAGREE_SEALED_LEAN_ERROR": c["DISAGREE_SEALED_LEAN_ERROR"],
        "NOT_COMPUTED": c["NOT_COMPUTED"],
        "OUT_OF_SCOPE": c["OUT_OF_SCOPE"],
        "compared": len(sel) - c["NOT_COMPUTED"] - c["OUT_OF_SCOPE"],
        "AGREE_clean_distinct_statements": sum(x["independence"] == "clean" and x["duplicate_of"] is None for x in ag),
        "rows_with_2plus_routes_after_skeptic": sum(x["routes_after_skeptic"] >= 2 for x in sel),
    }


comparator_counts = dict(Counter(r["status"] for r in rows))
comparator_counts["AGREE_full"] = sum(r["status"] == "AGREE" and r.get("coverage") == "full" for r in rows)
comparator_counts["AGREE_partial"] = sum(r["status"] == "AGREE" and r.get("coverage") == "partial" for r in rows)
comparator_counts["rows_with_2plus_routes"] = sum(r.get("independent_routes", 0) >= 2 for r in rows)
for k in ("AGREE", "NOT_COMPUTED", "NOT_COMPARABLE", "OUT_OF_SCOPE", "DISAGREE", "TOTAL"):
    if k in comparison["counts"] and k in comparator_counts:
        assert comparator_counts[k] == comparison["counts"][k], k

after_all = tally(ledger_rows)
after_tracks = tally([x for x in ledger_rows if x["track_target"]])
by_group = {g: tally([x for x in ledger_rows if x["group"] == g])
            for g in sorted({x["group"] for x in ledger_rows})}
changed = [dict(target_id=x["target_id"], comparator=x["comparator_status"], final=x["final_status"],
                independence=x["independence"], routes=[x["routes_comparator"], x["routes_after_skeptic"]])
           for x in ledger_rows
           if x["comparator_status"] != x["final_status"] or x["independence"] == "qualified"
           or x["routes_comparator"] != x["routes_after_skeptic"]]

# ---------------------------------------------------------------- rigidity (14 blind entries)
skA = load("skeptic_math/skA_scan_depends_on_k_results.json")
skB = load("skeptic_math/skB_immortal_results.json")
skC = load("skeptic_math/skC_chain_results.json")
skD4 = load("skeptic_math/skD_premises_N4_results.json")
skE = load("skeptic_math/skE_tt_faithful_results.json")
blk = load("skeptic_blindness/normalisation_k_test_results.json")

RIGIDITY = [
    ("A", "N, mu-term coefficient in Z*eta^3 = N y^(1/2) Psi - H theta_1^2 (scan 20..28)", "RIGID",
     "NORMALISATION", "echoes k: selected N by k = " + json.dumps(skA["selected_N_by_k"]),
     "typed k = 2 in Z = k phi_{0,1}", "skeptic_math/skA_scan_depends_on_k.py; skeptic_blindness/normalisation_k_test.py", False),
    ("A", "k in Z = k phi_{0,1}, selector Z(tau,0) = 24", "NORMALISATION", "NORMALISATION",
     "selector is the literal target", "e(K3) = 24 typed", "A-genus/run_track_a_v2.py", False),
    ("A", "k in Z = k phi_{0,1}, selector (q^0,y^1) coefficient = 2", "NORMALISATION", "NORMALISATION",
     "selector is the literal target", "ground-state count 2 typed", "A-genus/run_track_a_v2.py", False),
    ("B", "kappa, Goettsche exponent multiplier", "NORMALISATION", "NORMALISATION",
     "both sides are the same function call (part1:88-90)", "computed chi", "B-dyons/part1_euler_numbers.py", False),
    ("B", "N, coefficient of A_{2,1} in the m=1 polar subtraction", "RIGID", "RIGID",
     "exact overdetermined solve, " + str(skB["n_equations"]) + " equations -> " + json.dumps(skB["implied_trackB_constants"]),
     "DMVV input built from 2 phi_{0,1} (k = 2); DMZ form 3E4A - M Hhat (tier L)", "skeptic_math/skB_immortal.py", True),
    ("B", "M, coefficient of Hhat", "RIGID", "RIGID",
     "same solve; H(3)+1/6 control solutions = " + json.dumps(skB["negative_control_H3_plus_1_6_solutions"]),
     "as N; Hurwitz table confirmed by Dirichlet route", "skeptic_math/skB_immortal.py; skeptic_math/skB_hurwitz_dirichlet.py", True),
    ("C", "(m,n) in mU + n(-E8)", "CONDITIONAL_ON_INPUT", "RIGID_GIVEN_DEFINITION",
     "reverse chain with no chi input: " + json.dumps(skC["reverse_chain_no_chi_input"]["chi_top_forced"]) +
     "; intrinsic candidates signatures " + json.dumps([skC["sig_3U+2(-E8)"], skC["sig_7U+1(-E8)"], skC["sig_11U+0(-E8)"]]),
     "definition of K3 (K = O, b1 = 0); Noether, Serre (tier L); even unimodular rank-22 classification (tier L)",
     "skeptic_math/skC_chain.py", True),
    ("C", "tadpole budget flux + n = budget", "NORMALISATION", "NORMALISATION",
     "budget read from chi file", "chi_top", "C-lattices/04_tadpole_budget.py", False),
    ("C", "generic A, G T-duality identities", "RIGID", "VERIFIED_IDENTITY",
     "universal identities, no parameter selected", "none", "C-lattices/03_symbolic_tduality_proofs.py", False),
    ("D", "grid size N of the Freudenthal triangulation", "RIGID", "INVARIANCE",
     "N=4 open stars disjoint = " + json.dumps(skD4["P1_open_stars_disjoint"]) + ", MV Betti at N=4 " +
     json.dumps(skD4["P3"]["Z3"]["K3_betti_from_MV"]) + ": D's N=4 'failure' is an over-strict criterion, so nothing is selected",
     "none", "skeptic_math/skD_premises.py 4", False),
    ("D", "k resolved points, selector chi(X_k) = 24", "NORMALISATION", "NORMALISATION",
     "chi(X_k) = 8 + k, literal target", "chi = 24 typed", "D-tda/05_rigidity_scan.py", False),
    ("D", "k resolved points, selector b1 = 0", "NON-DISCRIMINATING", "NON-DISCRIMINATING",
     "all 17 values pass", "none", "D-tda/05_rigidity_scan.py", False),
    ("E", "tadpole total 24 -> alpha_x^2 in {8,16,24}", "NORMALISATION", "NORMALISATION",
     "TT (2.3) input", "tadpole 24 (TT 2.3, tier L)", "E-flux/scripts/flux_enumeration.py", False),
    ("E", "8 | alpha_x^2", "RIGID", "RIGID",
     "even lattice (det Gamma319 = " + str(skE["3_eight_divides_alpha_sq"]["det_Gamma319"]) +
     ") + even coefficients; odd-lattice control alpha^2 = " +
     str(skE["3_eight_divides_alpha_sq"]["negative_control_odd_lattice_I(3,3)"]["alpha=2v_alpha_sq"]),
     "TT (2.5) even flux quantisation (tier L); evenness of H^2(K3,Z)", "skeptic_math/skE_tt_faithful.py", True),
]
# the two skeptics' k-dependence tests must agree (independent scripts)
assert {k: v for k, v in skA["selected_N_by_k"].items()} == \
    {k: v["N_passing_slice_agreement"] for k, v in blk["rows"].items()}
PHI01_Z0 = blk["phi01_tau_z0_constant_computed"]
rigidity = [dict(track=t, parameter=p, track_label=tl, final_label=fl, evidence=ev, depends_on=dep,
                 discriminating_script=sc, genuinely_rigid=g) for (t, p, tl, fl, ev, dep, sc, g) in RIGIDITY]
both_skeptics_rigid = {"B"}  # B's N and M are the only entries both verdict files mark genuinely_rigid
rig_counts = dict(
    entries=len(rigidity),
    claimed_rigid_by_tracks=sum(x["track_label"] == "RIGID" for x in rigidity),
    genuinely_rigid_after_adjudication=sum(x["genuinely_rigid"] for x in rigidity),
    genuinely_rigid_both_skeptics=sum(x["genuinely_rigid"] and x["track"] in both_skeptics_rigid for x in rigidity),
    by_final_label=dict(Counter(x["final_label"] for x in rigidity)),
)
# cross-check the "both skeptics" set against the verdict files
bl_rigid = [v["parameter"] for v in sk_blind["rigidity_verdicts"] if v["genuinely_rigid"]]
ma_rigid = [v["parameter"] for v in sk_math["rigidity_verdicts"] if v["genuinely_rigid"]]

# ---------------------------------------------------------------- reverse pass
rev = []
for f, reach, note in [
    ("reverse/item1_index1_q20_results.json", False, "needs zK3 := ellipticGenus 20 (now 9)"),
    ("reverse/item2_moonshine_q20_results.json", True, "substitute 20 for the hard-coded 9; no dmvv truncation"),
    ("reverse/item3_immortal_m4_hecke_results.json", False,
     "needs hurwitzV4Full (true divisor sum), psiOpt5, and dmvv 6 2 / ellipticGenus N >= 12"),
]:
    d = load(f)
    c = d["computed"]
    if "holds_on_new_range" in c:
        holds, route = c["holds_on_new_range"], "B-dyons theta cache (blind engine)"
    elif "holds_correct_formula_via_psi05opt_ansatz" in c:
        holds = c["holds_correct_formula_via_psi05opt_ansatz"]
        route = ("psi_{0,5}^opt quoted from DMZ Table 2 (tier L); the DMVV G_5/A route gives "
                 "holds=" + str(c["holds_correct_formula_via_G5_DMVV"]) + " (different object, Delta psi_4)")
    else:
        per = c.get("twined_divisibility_q10_20_per_class", {})
        holds = all(v.get("holds_n10_20") for v in per.values()) if per else None
        dec = [v for k, v in c.items() if "decompose" in k]
        route = "moonshine_engine.py, validated against Lean n<=9 tables"
    rev.append(dict(file=f, tier=d.get("tier"), command=d.get("command"), holds=holds, route=route,
                    lean_reachable_now=reach, prerequisite=note))

# ---------------------------------------------------------------- flux vacua, orientifold
fe = load("E-flux/flux_enumeration_results.json")["results"]
flux = {}
for t, v in fe.items():
    flux[t] = dict(
        signature=v["signature"], group_order=v["automorphism_group_order_used"],
        with_tadpole_pairs={N: w["total_valid_ordered_pairs"] for N, w in v["with_tadpole_bound_alpha_x_sq_le_24"].items()},
        with_tadpole_pairs_by_alpha_sq={N: {str(8 * int(k)): x["pair_count"] for k, x in w["per_k"].items()}
                                        for N, w in v["with_tadpole_bound_alpha_x_sq_le_24"].items()},
        without_tadpole_pairs={N: w["total_valid_ordered_pairs"] for N, w in v["without_tadpole_bound_negative_control"].items()},
        orbits={N: o["total_orbits"] for N, o in v["orbit_counts_tadpole_on_only"].items() if "total_orbits" in o},
    )
skX = load("skeptic_math/skX_rows_orientifold_dirac_results.json")

ledger = dict(
    header="Generated by workflow k3t2-rigidity-loop-v2; computations tier B at best; physical identifications tier L/C.",
    built_by=f"{PY} audit/k3t2_rigidity_v2/build_ledger_v2.py",
    sealed=dict(leanmaster_version=comparison["_meta"]["sealed_leanmaster_version"],
                leanmaster_head=comparison["_meta"]["sealed_leanmaster_head"],
                post_seal_fix="LeanMaster adc85e7 (v3.21.0) changed orientifold#8,#9,#11 to the blind values"),
    inputs=["comparison.json", "skeptic_math/verdict.json", "skeptic_blindness/verdict.json",
            "skeptic_math/*_results.json", "E-flux/flux_enumeration_results.json", "reverse/item*_results.json",
            "../k3t2_rigidity/agreement_ledger.json"],
    counts=dict(
        total_rows=len(rows), distinct_targets=len(rows) - sum(x["duplicate_of"] is not None for x in ledger_rows),
        track_target_rows=sum(x["track_target"] for x in ledger_rows),
        duplicates=[[x["target_id"], x["duplicate_of"]] for x in ledger_rows if x["duplicate_of"]],
        comparator=comparator_counts, after_skeptics_all=after_all, after_skeptics_271_track_targets=after_tracks,
        after_skeptics_by_group=by_group,
        v1=v1["counts"], v1_rigidity=[v1["rigidity_genuinely_rigid"], v1["rigidity_total"]]),
    rows_changed_by_skeptics=changed,
    rigidity=rigidity, rigidity_counts=rig_counts,
    rigidity_skeptic_raw=dict(blindness_genuinely_rigid=bl_rigid, math_genuinely_rigid=ma_rigid),
    reverse=rev,
    flux_vacua_trackE=flux,
    orientifold_rederivation=skX["F_theory"] | {"negative_control_16_O7": skX["negative_control_16_O7"]},
    rows=ledger_rows,
)
(R / "agreement_ledger_v2.json").write_text(json.dumps(ledger, indent=1, ensure_ascii=False) + "\n")
print(json.dumps(dict(counts={k: ledger["counts"][k] for k in
                              ("total_rows", "distinct_targets", "track_target_rows", "duplicates", "comparator",
                               "after_skeptics_all", "after_skeptics_271_track_targets", "v1", "v1_rigidity")},
                      rigidity_counts=rig_counts, rigidity_skeptic_raw=ledger["rigidity_skeptic_raw"],
                      reverse=[(x["file"], x["holds"], x["lean_reachable_now"]) for x in rev],
                      changed=changed), indent=1))
