"""Build agreement_ledger_v3.json for the k3t2-rigidity-loop v3 (report stage).

Every count in the ledger is COMPUTED here from files already committed inside
audit/k3t2_rigidity_v3/ (this directory). Nothing outside the worktree is read,
and no literal physics number is typed into this script -- every number in the
output ledger is read out of a JSON file another script wrote.

Run from the repo root:
    cd audit/k3t2_rigidity_v3 && \
    /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python \
        build_ledger_v3.py

======================================================================
DOCUMENTED ADJUDICATION RULES (ground rule: "applies both skeptics'
disputes by explicit, documented rules")
======================================================================

R1 (comparison, skeptic_math dispute on tadpole_budget and its two siblings
    tadpole_cancellation, d3_tadpole_target_is_24):
    These three AGREE rows have independent_routes=0 and their own note says
    the same literature constant (24, TT eq. 2.3) was typed on both the Lean
    side and the blind side -- not independently derived by either. Per the
    skeptic's recommendation, such rows are relabeled SAME_DECLARED_INPUT
    (a distinct status from AGREE) so a reader scanning AGREE counts does not
    read them as a substantive cross-check.

R2 (comparison, skeptic_blindness dispute on the chi=24 four-row family
    k3_euler_characteristic / k3_euler_eq_24 / euler_K3 / ellipticGenus_z0):
    The dispute is whether "chi_top = (Noether's constant 12) x (computed
    chi(O)=2)" is an independent second route, next to D-tda's GUDHI
    Mayer-Vietoris route. Orchestrator ruling (the skeptic explicitly asked
    for one): a route built by typing a FIXED, universal theorem constant
    (Noether's 12 is the same for every complex surface; nothing tunes it)
    is not the kind of inserted free parameter ground rule 7 warns about,
    which is a CHOSEN normalisation such as Track A's k in Z = k*phi_{0,1}.
    So these four rows keep independent_routes=2 and status AGREE, but are
    flagged reviewed_and_defended=True with the dispute text attached, so
    the reader sees both the objection and the ruling.

R3 (rigidity, skeptic disagreement on the same track/parameter entry):
    When BOTH skeptics reviewed the same rigidity entry (matched by an
    explicit, hand-verified table below, not by fuzzy text matching) and
    they reach a DIFFERENT verdict -- either a different genuinely_rigid
    boolean, or the same boolean but a different proposed final_label --
    the entry is marked status="DISPUTED", genuinely_rigid=False (an
    unresolved dispute does not count as confirmed rigidity), and BOTH
    skeptics' verdicts are preserved verbatim for the reader to weigh.
    This covers: D-tda's j-selector (blind: RIGID_GIVEN_DEFINITION/True,
    math: NORMALISATION/False) and C-lattices' chi_top (blind:
    CONDITIONAL_ON_INPUT, math: VERIFIED_IDENTITY -- both agree "not
    RIGID_GIVEN_DEFINITION" but propose different, not strictly ordered,
    replacement labels).

R4 (rigidity, only one skeptic reviewed the entry):
    That skeptic's verdict (final_label, genuinely_rigid) replaces the
    track's own self-reported label; status="SKEPTIC_CORRECTED" if it
    differs from the track's own label, else "SKEPTIC_CONFIRMED".

R5 (rigidity, neither skeptic reviewed the entry):
    The track's own self-reported label is kept verbatim, status=
    "SELF_REPORTED_UNREVIEWED", and it is NOT counted in
    genuinely_rigid_after_skeptics (self-report alone is not adjudication).

R6 (reverse pass, skeptic disagreement on "stands"):
    When the two skeptics' reverse_verdicts disagree on "stands" for the
    same item (matched by list position -- both files review the same 5
    items in the same order, verified by an assertion below), the item's
    final status is "DISPUTED" and BOTH reasons are kept. This covers only
    reverse_3 (class_number_one): math endorses fully (the classification is
    in fact fully proved for all D by Baker-Heegner-Stark, not just D<=3000);
    blind's substantive objection (5 of the 9 members are consecutive
    discriminants that all pass, so they are not isolated from a neighbour
    the way a RIGID label implies) is a different, unrebutted point about
    the negative control, so it is kept as the operative caveat.
"""
import json
from collections import Counter
from pathlib import Path

R = Path(__file__).resolve().parent
PY = "/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python"
TRACKS = ("A-genus", "B-dyons", "C-lattices", "D-tda", "E-flux", "F-whichk3")

# A rigidity label is not a single linear scale (a VERIFIED_IDENTITY is not
# "weaker" than a CONDITIONAL_ON_INPUT the way a rank order would imply --
# they are different KINDS of claim), so disagreement is never silently
# resolved by picking one; see R3.
RIGID_FAMILY = {"RIGID", "RIGID_GIVEN_DEFINITION"}


def load(rel):
    return json.loads((R / rel).read_text())


def as_list(x):
    """Some tracks (F-whichk3) store results/rigidity as {id: entry} dicts;
    normalise everything to a list of entries carrying their own id."""
    if isinstance(x, dict):
        out = []
        for k, v in x.items():
            if isinstance(v, dict) and "id" not in v:
                v = dict(v, id=k)
            out.append(v)
        return out
    return x


# ---------------------------------------------------------------- load everything
comparison = load("comparison/comparison.json")
chain = load("chain/chain_check_results.json")
sk_blind = load("skeptic_blindness/verdict.json")
sk_math = load("skeptic_math/verdict.json")
reverse = load("reverse/results.json")

track_data = {}
for t in TRACKS:
    d = load(f"{t}/results.json")
    inp = load(f"{t}/inputs.json")
    track_data[t] = dict(
        raw=d,
        results=as_list(d.get("results", [])),
        rigidity=as_list(d.get("rigidity", [])),
        could_not_do=d.get("could_not_do", []),
        inputs=as_list(inp),
    )

# ================================================================== SECTION: comparison
rows = comparison["rows"]
ids = [r["target_id"] for r in rows]
assert len(ids) == len(set(ids)), "duplicate target ids in comparison.json"

TADPOLE_FAMILY = {"tadpole_budget", "tadpole_cancellation", "d3_tadpole_target_is_24"}
assert sk_math["disputed_rows"][0]["target_id"] == "tadpole_budget"
assert "chi_top" in sk_blind["disputed_rows"][0]["target_id"] or \
    "k3_euler_characteristic" in sk_blind["disputed_rows"][0]["target_id"]

# R2 is detected STRUCTURALLY (any row whose shared_inputs contains the fixed
# theorem constant 'noether_formula' AND has >=2 declared routes), not by a
# hand-picked list of target_ids -- otherwise the ruling would silently cover
# only the 4 rows the blindness skeptic happened to name, while 2 more rows
# (k3_second_betti_hodge, rank_K3) share exactly the same Noether-chain
# structure and would be left inconsistently un-ruled.
NOETHER_ROUTE_INPUT = "noether_formula"

comparison_rows = []
for r in rows:
    tid = r["target_id"]
    status0 = r["status"]
    status = status0
    routes = r.get("independent_routes", 0)
    shared = r.get("shared_inputs", [])
    ruling = None
    r2_applies = status0 == "AGREE" and routes >= 2 and NOETHER_ROUTE_INPUT in shared
    if tid in TADPOLE_FAMILY and status0 == "AGREE":
        status = "SAME_DECLARED_INPUT"
        ruling = "R1: " + sk_math["disputed_rows"][0]["reason"]
    elif r2_applies:
        ruling = ("R2 (orchestrator ruling on skeptic_blindness dispute, applied structurally to "
                   "every >=2-route row sharing 'noether_formula', not only the 4 rows the skeptic "
                   "named): a route built by typing a FIXED, universal theorem constant (Noether's "
                   "12, or its Hirzebruch/Hodge-index consequences) x an independently-computed "
                   "quantity counts as an independent route, unlike a chosen/tunable normalisation; "
                   "kept AGREE, routes unchanged, flagged reviewed_and_defended. Dispute text (as "
                   "raised for the chi=24 subset): " + sk_blind["disputed_rows"][0]["reason"])
    comparison_rows.append(dict(
        target_id=tid, sealed_track=r.get("sealed_track"), lean_file_line=r.get("lean_file_line"),
        lean_value=r.get("lean_value"), blind_source=r.get("blind_source"), blind_value=r.get("blind_value"),
        comparator_status=status0, final_status=status, independent_routes=routes,
        reviewed_and_defended=r2_applies,
        shared_inputs=shared, note=r.get("note"), ruling=ruling,
        post_seal_changed=r.get("post_seal_changed", False),
    ))

comp_status_counts = dict(Counter(r["final_status"] for r in comparison_rows))
by_track_tally = {}
for r in comparison_rows:
    g = by_track_tally.setdefault(r["sealed_track"], Counter())
    g[r["final_status"]] += 1
by_track_tally = {k: dict(v) for k, v in by_track_tally.items()}
assert sum(sum(v.values()) for v in by_track_tally.values()) == len(comparison_rows)

comp_counts = dict(
    total_rows=len(comparison_rows),
    comparator_counts=dict(Counter(r["comparator_status"] for r in comparison_rows)),
    final_status_counts=comp_status_counts,
    rows_with_2plus_routes=sum(r["independent_routes"] >= 2 for r in comparison_rows),
    rows_reviewed_and_defended_R2=sum(r["reviewed_and_defended"] for r in comparison_rows),
    rows_with_2plus_routes_if_R2_rejected=sum(
        r["independent_routes"] >= 2 and not r["reviewed_and_defended"] for r in comparison_rows),
    rows_downgraded_same_declared_input=sum(r["final_status"] == "SAME_DECLARED_INPUT" for r in comparison_rows),
    rows_post_seal_changed=sum(r["post_seal_changed"] for r in comparison_rows),
    blind_only_count=len(comparison.get("blind_only", [])),
    by_sealed_track=by_track_tally,
)
for k in ("AGREE", "NOT_COMPUTED"):
    if k in comparison["counts"]:
        assert comp_counts["comparator_counts"].get(k, 0) == comparison["counts"][k], k

# ================================================================== SECTION: declared inputs
declared_inputs = {}
for t in TRACKS:
    consumers = {}
    for res in track_data[t]["results"]:
        for name in (res.get("shared_inputs") or []):
            consumers.setdefault(name, []).append(res.get("id") or res.get("quantity"))
    entries = []
    for item in track_data[t]["inputs"]:
        name = item.get("name")
        entries.append(dict(
            name=name, value=item.get("value"), tier=item.get("tier"), why=item.get("why"),
            depended_on_by=consumers.get(name, []),
        ))
    declared_inputs[t] = entries

# reverse pass has its own declared inputs
rev_inputs = load("reverse/inputs.json")
rev_consumers = {}
for e in reverse["entries"]:
    for name in (e.get("shared_inputs") or []):
        rev_consumers.setdefault(name, []).append(e["id"])
declared_inputs["reverse"] = [
    dict(name=i.get("name"), value=i.get("value"), tier=i.get("tier"), why=i.get("why"),
         depended_on_by=rev_consumers.get(i.get("name"), []))
    for i in as_list(rev_inputs)
]

inputs_counts = {t: len(v) for t, v in declared_inputs.items()}
total_inputs = sum(inputs_counts.values())

# FROM MEMORY tally: an input is "from memory" if its own 'why' field says so
# (ground rule 12 requires tier-L inputs quoted from memory to say so there).
from_memory = []
for t, entries in declared_inputs.items():
    for i in entries:
        if "FROM MEMORY" in str(i.get("why") or ""):
            from_memory.append(dict(track=t, name=i["name"], tier=i.get("tier")))
by_tier = Counter(i.get("tier") for v in declared_inputs.values() for i in v)

# ================================================================== SECTION: full per-track results
# Verbatim (not re-typed) computed values, so sections 6 (flux vacua) and 7
# (which K3) can cite exact numbers straight from the committed JSON.
results_by_track = {}
for t in TRACKS:
    results_by_track[t] = [
        dict(id=res.get("id"), quantity=res.get("quantity"), computed=res.get("computed"),
             shared_inputs=res.get("shared_inputs", []), script=res.get("script"),
             command=res.get("command"), note=res.get("note"))
        for res in track_data[t]["results"]
    ]
results_counts = {t: len(v) for t, v in results_by_track.items()}

# ================================================================== SECTION: could_not_do / not done
could_not_do_by_track = {t: track_data[t]["could_not_do"] for t in TRACKS}
could_not_do_counts = {t: len(v) for t, v in could_not_do_by_track.items()}

# ================================================================== SECTION: rigidity
def flatten_rigidity(track, entries):
    """Expand a per-level classification (A-genus's chi(g) at N in
    {2,3,5,7}, stored either as a real dict or as a "per level: {...}"
    JSON-in-a-string) into one row per level; pass everything else through."""
    out = []
    for e in entries:
        cls = e.get("classification")
        sol = e.get("solution_set")
        if isinstance(cls, str) and cls.startswith("per level:"):
            cls = json.loads(cls.split("per level:", 1)[1].strip())
        if isinstance(sol, str) and sol.strip().startswith("{"):
            try:
                sol = json.loads(sol)
            except json.JSONDecodeError:
                pass
        if isinstance(cls, dict):
            for lvl, c in cls.items():
                out.append(dict(e, classification=c, sub_level=lvl,
                                 solution_set=sol.get(lvl) if isinstance(sol, dict) else sol))
        else:
            out.append(dict(e, sub_level=None))
    return [dict(x, track=track) for x in out]


all_rigidity = []
for t in TRACKS:
    all_rigidity.extend(flatten_rigidity(t, track_data[t]["rigidity"]))

bl_rig = sk_blind["rigidity_verdicts"]
ma_rig = sk_math["rigidity_verdicts"]
assert "D-tda: j" in bl_rig[2]["parameter"] and "Track D: j" in ma_rig[0]["parameter"]
assert "A-genus: chi(g) at level 5" in bl_rig[3]["parameter"] and "Track A: chi(g) at level 5 and 7" in ma_rig[3]["parameter"]
assert "C-lattices: chi_top" in bl_rig[0]["parameter"] and "Track C: chi_top" in ma_rig[1]["parameter"]
assert "U+n(-E8)" in bl_rig[1]["parameter"]
assert "B-dyons: overall factor k" in bl_rig[4]["parameter"]
assert "Track A: k" in ma_rig[2]["parameter"]

rigidity_rows = []
for e in all_rigidity:
    track = e["track"]
    param = e.get("parameter_inserted")
    base_label = e.get("classification")
    bl = ma = None
    if track == "C-lattices" and param and "chi_top" in str(param):
        bl, ma = bl_rig[0], ma_rig[1]
    elif track == "C-lattices" and param and "mU+n(-E8)" in str(param):
        bl, ma = bl_rig[1], None
    elif track == "D-tda" and param == "j" and base_label == "RIGID_GIVEN_DEFINITION":
        bl, ma = bl_rig[2], ma_rig[0]
    elif track == "A-genus" and e.get("sub_level") == "5" and "chi(g)" in str(param):
        bl, ma = bl_rig[3], ma_rig[3]  # blind: "level 5" specifically; math: "level 5 and 7"
    elif track == "A-genus" and e.get("sub_level") == "7" and "chi(g)" in str(param):
        ma = ma_rig[3]  # math's statement explicitly covers level 7 too; blind did not review level 7
    elif track == "A-genus" and param and "overall factor" in str(param):
        ma = ma_rig[2]
    elif track == "B-dyons" and param and "overall factor k" in str(param):
        bl = bl_rig[4]

    reviewers = [x for x in (bl, ma) if x]
    if bl and ma:
        agree = (bl["genuinely_rigid"] == ma["genuinely_rigid"]) and \
                (bl["final_label"].split(" ")[0].split("(")[0].strip() ==
                 ma["final_label"].split(" ")[0].split("(")[0].strip())
        if agree:
            status = "SKEPTIC_CONFIRMED" if base_label in RIGID_FAMILY and bl["genuinely_rigid"] else "SKEPTIC_CORRECTED"
            final_label = bl["final_label"]
            genuinely_rigid = bl["genuinely_rigid"]
        else:
            status = "DISPUTED"
            final_label = f"DISPUTED: blind={bl['final_label']!r} vs math={ma['final_label']!r}"
            genuinely_rigid = False  # R3: unresolved dispute -> not counted as confirmed
    elif reviewers:
        s = reviewers[0]
        final_label = s["final_label"]
        genuinely_rigid = s["genuinely_rigid"]
        status = "SKEPTIC_CONFIRMED" if final_label.split(" ")[0].split("(")[0].strip() == str(base_label) else "SKEPTIC_CORRECTED"
    else:
        status = "SELF_REPORTED_UNREVIEWED"
        final_label = base_label
        genuinely_rigid = None  # not adjudicated; not counted either way

    rigidity_rows.append(dict(
        track=track, parameter_inserted=param, sub_level=e.get("sub_level"),
        base_label_self_reported=base_label, solution_set=e.get("solution_set"),
        negative_control=e.get("negative_control"),
        control_perturbs_same_parameter=e.get("control_perturbs_same_parameter"),
        input_it_depends_on=e.get("input_it_depends_on"),
        status=status, final_label=final_label, genuinely_rigid=genuinely_rigid,
        skeptic_blindness_verdict=bl, skeptic_math_verdict=ma,
    ))

rigidity_counts = dict(
    total_entries=len(rigidity_rows),
    self_claimed_rigid_family=sum(x["base_label_self_reported"] in RIGID_FAMILY for x in rigidity_rows),
    reviewed_by_at_least_one_skeptic=sum(x["status"] != "SELF_REPORTED_UNREVIEWED" for x in rigidity_rows),
    disputed=sum(x["status"] == "DISPUTED" for x in rigidity_rows),
    genuinely_rigid_confirmed=sum(x["genuinely_rigid"] is True for x in rigidity_rows),
    genuinely_rigid_confirmed_by_both_skeptics=sum(
        bool(x["genuinely_rigid"] is True and x["skeptic_blindness_verdict"] and x["skeptic_math_verdict"])
        for x in rigidity_rows),
    self_reported_unreviewed=sum(x["status"] == "SELF_REPORTED_UNREVIEWED" for x in rigidity_rows),
    by_final_status=dict(Counter(x["status"] for x in rigidity_rows)),
)

# ================================================================== SECTION: reverse pass
bl_rev = sk_blind["reverse_verdicts"]
ma_rev = sk_math["reverse_verdicts"]
assert len(bl_rev) == len(ma_rev) == len(reverse["entries"]) == 5
for i, e in enumerate(reverse["entries"]):
    assert e["id"].endswith(bl_rev[i]["item"]) or bl_rev[i]["item"] in e["id"], (e["id"], bl_rev[i]["item"])

reverse_rows = []
for i, e in enumerate(reverse["entries"]):
    b, m = bl_rev[i], ma_rev[i]
    if b["stands"] == m["stands"]:
        final_stands = b["stands"]
        status = "SKEPTIC_AGREE"
    else:
        final_stands = "DISPUTED"
        status = "DISPUTED"
    reverse_rows.append(dict(
        id=e["id"], from_theorem=e.get("from_theorem"), holds=e.get("holds"),
        rigidity_label=e.get("rigidity_label"), status=status, final_stands=final_stands,
        skeptic_blindness_reason=b["reason"], skeptic_math_reason=m["reason"],
    ))

reverse_counts = dict(
    total=len(reverse_rows),
    holds_true=sum(x["holds"] is True for x in reverse_rows),
    holds_false=sum(x["holds"] is False for x in reverse_rows),
    both_skeptics_agree_stands=sum(x["status"] == "SKEPTIC_AGREE" and x["final_stands"] is True for x in reverse_rows),
    disputed=sum(x["status"] == "DISPUTED" for x in reverse_rows),
)

# ================================================================== SECTION: chain
chain_rows = []
for i, link in enumerate(chain["links"]):
    chain_rows.append(dict(index=i + 1, **link) if isinstance(link, dict) else dict(index=i + 1, text=link))
# The same R2 ruling that governs the comparison rows governs one chain link:
# link "C.chi_top (Noether) == D.chi ..." is independent_corroboration=True only
# because a fixed theorem constant (Noether's 12) x a computed number is being
# counted as a genuine cross-method route -- the identical structural call.
r2_chain_links = [l for l in chain_rows if "Noether" in str(l.get("caveat", ""))
                   and l.get("independent_corroboration") is True]
chain_counts = dict(
    n_links=chain.get("n_links", len(chain_rows)),
    all_consistent=chain.get("all_consistent"),
    n_independent_corroborations=chain.get("n_independent_corroborations"),
    n_independent_corroborations_if_R2_rejected=(
        (chain.get("n_independent_corroborations") or 0) - len(r2_chain_links)),
    independent_link_names=chain.get("independent_link_names"),
    r2_dependent_link_names=[l["link"] for l in r2_chain_links],
    n_pointer_links=sum(l.get("method") == "POINTER" for l in chain_rows),
    n_cross_method_links=sum(l.get("method") == "CROSS_METHOD" for l in chain_rows),
)

# ================================================================== SECTION: skeptics summary
skeptics_meta = dict(
    blindness=dict(commit=sk_blind.get("commit"), skeptic=sk_blind.get("skeptic"),
                    n_violations=len(sk_blind.get("violations", [])),
                    n_reruns=len(sk_blind.get("reruns", [])),
                    n_disputed_rows=len(sk_blind.get("disputed_rows", []))),
    math=dict(commit=sk_math.get("commit"), lens=sk_math.get("lens"),
               n_violations=len(sk_math.get("violations", [])),
               n_reruns=len(sk_math.get("reruns", [])),
               n_disputed_rows=len(sk_math.get("disputed_rows", []))),
)

# ================================================================== assemble & write
ledger = dict(
    header="Generated by workflow k3t2-rigidity-loop-v3; computations tier B at best; physical identifications tier L/C.",
    built_by=f"cd audit/k3t2_rigidity_v3 && {PY} build_ledger_v3.py",
    sealed=dict(
        leanmaster_tag=comparison["meta"]["sealed_leanmaster_tag"],
        leanmaster_head=comparison["meta"]["sealed_leanmaster_head"],
        n_sealed_targets=comparison["meta"]["n_sealed_targets"],
    ),
    inputs_read=[
        "comparison/comparison.json", "chain/chain_check_results.json",
        "skeptic_blindness/verdict.json", "skeptic_math/verdict.json", "reverse/results.json",
        "reverse/inputs.json",
    ] + [f"{t}/results.json, {t}/inputs.json" for t in TRACKS],
    rules_applied=["R1", "R2", "R3", "R4", "R5", "R6"],
    comparison=dict(counts=comp_counts, rows=comparison_rows),
    declared_inputs=dict(counts=inputs_counts, total=total_inputs, by_tier=dict(by_tier),
                          from_memory=from_memory, n_from_memory=len(from_memory), by_track=declared_inputs),
    results_by_track=dict(counts=results_counts, by_track=results_by_track),
    could_not_do=dict(counts=could_not_do_counts, total=sum(could_not_do_counts.values()),
                       by_track=could_not_do_by_track),
    rigidity=dict(counts=rigidity_counts, rows=rigidity_rows),
    reverse=dict(counts=reverse_counts, rows=reverse_rows),
    chain=dict(counts=chain_counts, links=chain_rows),
    skeptics=skeptics_meta,
)

(R / "agreement_ledger_v3.json").write_text(json.dumps(ledger, indent=1, ensure_ascii=False) + "\n")
print(json.dumps(dict(
    comparison_counts=comp_counts, inputs_counts=inputs_counts, n_from_memory=len(from_memory),
    results_counts=results_counts, could_not_do_counts=could_not_do_counts,
    rigidity_counts=rigidity_counts, reverse_counts=reverse_counts, chain_counts=chain_counts,
), indent=1))
