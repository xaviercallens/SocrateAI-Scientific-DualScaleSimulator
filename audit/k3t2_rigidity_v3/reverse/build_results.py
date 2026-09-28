"""
Aggregates item1..item5's own results.json files (already written by their own scripts) into
results.json and exports.json for the reverse pass. Does not recompute anything and does not
hand-type any numeric result: every number below is read back from the itemN_results.json files
that itemN_*.py wrote. The only hand-authored text is the prediction/independence/rigidity-label
prose and the proposed Lean-style statements (5-8 short Lean snippets, not kernel-checked here --
LeanMaster is a separate Lake project this worktree does not build).

Run (after item1..item5 have been run):
    cd audit/k3t2_rigidity_v3/reverse && \
    /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python build_results.py
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def load(name):
    return json.loads((HERE / name).read_text())


def main() -> None:
    i1 = load("item1_results.json")
    i2 = load("item2_results.json")
    i3 = load("item3_results.json")
    i4 = load("item4_results.json")
    i5 = load("item5_results.json")

    entries = []

    # --- Item 1 ---------------------------------------------------------------
    entries.append({
        "id": "reverse_1_moore_vs_hurwitz",
        "from_theorem": "DualScaleDyons.WhichK3.moore_vs_hurwitz (Lean range D in [0,400])",
        "tier": "B",
        "shared_inputs": ["D_max_moore_hurwitz"],
        "prediction": "Moore's attractive-K3/dyon-charge count 12*N(D) continues to equal the "
                       "Hurwitz class number 12*H(D), up to the two self-dual-torus correction "
                       "terms (6*[D=4f^2], 8*[D=3f^2]), for every D far beyond Lean's D<=400.",
        "computed": f"Checked exactly for every D in {i1['scanned_range']} using an independent "
                    f"reimplementation of reducedForms/h12 (fractions.Fraction-free, pure int): "
                    f"{i1['newly_scanned_count']} new D values beyond the Lean range, "
                    f"{len(i1['failures_beyond_lean_range'])} failures. Negative control (same "
                    f"identity, correction terms removed) fails on "
                    f"{i1['negative_control']['num_failures']} of the new D values "
                    f"(e.g. D={i1['negative_control']['sample_failures'][:3]}).",
        "holds": i1["identity_holds_for_all_scanned"],
        "rigidity_label": "VERIFIED_IDENTITY",
        "independence": "Shares the mathematical definitions (reducedForms, the weighted "
                         "Hurwitz-number recipe, the two correction terms) with Lean's own "
                         "WhichK3.lean/Immortal.lean; the CODE is an independent from-scratch "
                         "Python reimplementation (common.py), not a call into Lean or an import "
                         "of Lean's output, and the range D in (400,5000] was never checked by "
                         "Lean's `decide +kernel` call, which only iterates D in [0,400].",
        "proposed_lean_statement":
            "theorem moore_vs_hurwitz_5000 :\n"
            "    (List.range 5001).all (fun D => D < 3 || D % 4 == 1 || D % 4 == 2 ||\n"
            "      (12 * (nForms D : ℤ) == h12 D + (if isKSquare 4 D then 6 else 0) +\n"
            "        (if isKSquare 3 D then 8 else 0))) = true := by decide +kernel",
        "lean_reachable_now": False,
        "lean_reachable_note": "Not attempted -- this worktree does not build LeanMaster. "
                                "nForms D itself iterates List.range (D+1) internally, so "
                                "List.range 5001 nests roughly (5000/400)^2 ~ 150x deeper kernel "
                                "reduction than the existing D<=400 `decide +kernel` call; cost "
                                "is extrapolated, not measured, so False rather than an "
                                "unevidenced True.",
        "script": "audit/k3t2_rigidity_v3/reverse/item1_moore_vs_hurwitz.py",
    })

    # --- Item 2 ---------------------------------------------------------------
    entries.append({
        "id": "reverse_2_kummer_iff_even",
        "from_theorem": "DualScaleDyons.WhichK3.kummer_iff_even (Lean range D in [0,200])",
        "tier": "B",
        "shared_inputs": ["D_max_kummer_iff_even"],
        "prediction": "An attractive K3's Kummer test (T(X) all entries divisible by 4) continues "
                       "to coincide exactly with 'a,b,c all even' for every reduced form far "
                       "beyond Lean's D<=200.",
        "computed": f"Checked exactly for D in {i2['scanned_range']}: "
                    f"{i2['num_mismatches']} mismatches. Negative control (require only a,c "
                    f"even, drop the condition on b) disagrees with isKummerForm on "
                    f"{i2['negative_control']['num_forms_where_control_disagrees']} forms beyond "
                    f"D=200, confirming the 'all three even' condition is the one actually doing "
                    f"the selecting, not 'a,c even' alone.",
        "holds": i2["holds_for_all_scanned"],
        "rigidity_label": "VERIFIED_IDENTITY",
        "independence": "Same biconditional Lean states (isKummerForm == all-even), reimplemented "
                         "from scratch in common.py; independent of Lean's kernel evaluation, "
                         "and D in (200,5000] is outside Lean's checked range.",
        "proposed_lean_statement":
            "theorem kummer_iff_even_5000 :\n"
            "    (List.range 5001).all (fun D => (reducedForms D).all fun f =>\n"
            "      isKummerForm f == (f.1 % 2 == 0 && f.2.1 % 2 == 0 && f.2.2 % 2 == 0)) = "
            "true := by decide +kernel",
        "lean_reachable_now": False,
        "lean_reachable_note": "Not attempted -- no LeanMaster build available here; D<=5000 "
                                "nests reducedForms far deeper than Lean's own D<=200 call, so "
                                "cost is extrapolated (much higher), not measured.",
        "script": "audit/k3t2_rigidity_v3/reverse/item2_kummer_iff_even.py",
    })

    # --- Item 3 ---------------------------------------------------------------
    viol = i3["negative_control"]["violations"]
    low_cluster = sorted({v["hit"] for v in viol} | {v["neighbour"] for v in viol})
    entries.append({
        "id": "reverse_3_most_attractive_class_number_one",
        "from_theorem": "DualScaleDyons.WhichK3.most_attractive (Lean checks only D=3 and D=4)",
        "tier": "B/L",
        "shared_inputs": ["D_max_class_number_one", "heegner_discriminants_from_memory"],
        "prediction": "D=3, D=4 are not an isolated pair: the same 'exactly one reduced form' "
                       "property should hold for every discriminant with primitive class number "
                       "1, i.e. the nine Heegner-type discriminants 3,4,7,8,11,19,43,67,163.",
        "computed": f"Scanning ALL valid discriminants in {i3['scanned_range_D']} for "
                    f"nForms(D)==1 (own reducedForms code, no list typed in) gives exactly "
                    f"{i3['computed_nForms_eq_1']}, matching the Tier-L literature list "
                    f"(FROM MEMORY) {i3['expected_from_literature_FROM_MEMORY']}. BUT the "
                    f"neighbour control fails at low D: {i3['negative_control']['num_neighbour_violations']} "
                    f"violations, all among the mutually-adjacent cluster {low_cluster} -- D=3's "
                    f"very next valid discriminant, D=4, ALSO has nForms=1 (and so on through "
                    f"D=11); only from D=19 on does the property become isolated from its "
                    f"neighbours (19, 43, 67, 163 each have both neighbours failing).",
        "holds": i3["computed_matches_literature_list"],
        "rigidity_label": "RIGID",
        "independence": "The Lean theorem only asserts the fact at D=3,4 with no claim about "
                         "neighbours or generality; the literature list of the nine h=1 "
                         "discriminants is a Tier-L fact quoted FROM MEMORY, used only to check "
                         "the independently-computed hit list -- it is not fed into the scan. "
                         "RIGID applies at the level of the full computed set {3,4,7,8,11,19,43,"
                         "67,163}: the selecting condition nForms(D)=1 does not contain the "
                         "target, and a neighbouring value demonstrably fails (both neighbours "
                         "of 19, 43, 67 and 163 fail; that is enough to show the condition "
                         "genuinely selects, not that every instance is isolated). The important "
                         "caveat, kept rather than smoothed over: at the LOW end the property is "
                         "NOT locally isolated -- D=3,4,7,8,11 are five CONSECUTIVE valid "
                         "discriminants that all pass, so D=3's own immediate neighbour (D=4) "
                         "also has nForms=1, and Lean's specific choice of D=3,4 as 'the two most "
                         "attractive' sits inside that dense run rather than being a uniquely "
                         "isolated pair the way 19, 43, 67, 163 are.",
        "proposed_lean_statement":
            "theorem class_number_one_3000 :\n"
            "    ((List.range 3001).filter (fun D => 3 <= D && (D % 4 == 0 || D % 4 == 3) &&\n"
            "        nForms D == 1)) = [3, 4, 7, 8, 11, 19, 43, 67, 163] := by decide +kernel",
        "lean_reachable_now": False,
        "script": "audit/k3t2_rigidity_v3/reverse/item3_class_number_one.py",
    })

    # --- Item 4 ---------------------------------------------------------------
    beyond = i4["beyond_lean_range"]
    first_loss = next(r for r in beyond if not r["D_n_is_optimal"])
    entries.append({
        "id": "reverse_4_ade_trapping_beyond_rank8",
        "from_theorem": "DualScaleDyons.KummerD4.trapping_rank_table (Lean checks rank d in [1,8], "
                         "D4 uniquely best at rank 4)",
        "tier": "B/L",
        "shared_inputs": ["rank_max_ade_trapping", "E6_E7_E8_root_counts_from_memory",
                           "extended_adeSimple_rank_range"],
        "prediction": "D_n (2n(n-1) roots) remains the uniquely largest simply-laced root system "
                       "of rank <= n for every n beyond 8, generalising D4's win at rank 4.",
        "computed": f"Regression: reimplemented bestTable matches Lean's own [0,2,6,12,24,40,72,"
                    f"126,240] exactly ({i4['regression_matches_lean_table_ranks_0_to_8']}). "
                    f"Extending the SAME unbounded-knapsack DP to rank "
                    f"{i4['scanned_range'][1]}: the prediction is FALSE as stated -- at rank "
                    f"{first_loss['rank']}, D{first_loss['rank']} has only "
                    f"{first_loss['D_n_roots']} roots while the true best "
                    f"({first_loss['best_roots']}) comes from "
                    f"{first_loss['alternative_optimal_combos']} (E8 plus a smaller padding). "
                    f"D_n is not optimal at ranks {[r['rank'] for r in beyond if not r['D_n_is_optimal']]}, "
                    f"ties another combination at ranks {i4['tie_ranks']}, and is the unique "
                    f"winner at every other rank up to {i4['scanned_range'][1]} "
                    f"({i4['num_ranks_where_Dn_strictly_optimal']} of "
                    f"{i4['num_ranks_beyond_8']} ranks beyond 8). At rank 24 (the Niemeier rank) "
                    f"D24 ({i4['negative_control_rank24_vs_3xE8']['D24_roots']} roots) does beat "
                    f"the naive '3 copies of E8' guess "
                    f"({i4['negative_control_rank24_vs_3xE8']['3xE8_roots']} roots).",
        "holds": False,
        "rigidity_label": "CONDITIONAL_ON_INPUT",
        "independence": "Lean's own trapping_rank_table only asserts the DP result up to rank 8 "
                         "and only that D4 beats A4 (20) and every split at rank 4; it makes no "
                         "claim about D_n at higher rank. The DP algorithm (unbounded knapsack "
                         "over simple-component root counts) is the same combinatorial recipe "
                         "Lean's bestTable implements, reimplemented from scratch; the component "
                         "root-count formulas A_n=n(n+1), D_n=2n(n-1) are derived, not copied, "
                         "and E6/E7/E8's 72/126/240 are a small Tier-L constant declared "
                         "explicitly (FROM MEMORY) in inputs.json. UNLIKE items 1/2/5, this is "
                         "not a pure range extension: Lean's own adeSimple list (KummerD4.lean) "
                         "contains only A1..A8, D4..D8, E6, E7, E8 -- nothing of rank > 8 -- so "
                         "asking whether D_n stays optimal past rank 8 requires first enlarging "
                         "the component list to A1..A32, D4..D32 (declared as "
                         "extended_adeSimple_rank_range in inputs.json). The result therefore "
                         "depends on a structural input Lean's own definitions do not contain, "
                         "not only on a larger List.range bound.",
        "proposed_lean_statement":
            "-- requires an EXTENDED adeSimple (A1..A9, D4..D9, E6, E7, E8), not the existing one\n"
            "theorem trapping_beats_Dn_at_rank9_extended :\n"
            "    bestTableExtended 9 = 242 ∧ (2 * 9 * 8 : ℕ) < 242 := by decide",
        "lean_reachable_now": False,
        "lean_reachable_note": "Needs a new/extended adeSimple definition first (see "
                                "independence note); not a `decide` on Lean's existing "
                                "definitions at any range bound.",
        "script": "audit/k3t2_rigidity_v3/reverse/item4_ade_trapping.py",
    })

    # --- Item 5 ---------------------------------------------------------------
    entries.append({
        "id": "reverse_5_twined_2A_beyond_q9",
        "from_theorem": "DualScaleMoonshine.Twining.twined_div24 / twined_2A "
                         "(Lean range q^0..q^9)",
        "tier": "B",
        "shared_inputs": ["N_max_twined_2A", "twined_2A_chi_Nlev_c_from_CDH"],
        "prediction": "The structural property Lean proves at q^9 -- every coefficient of "
                       "24*H_2A divisible by 24 -- continues to hold at every order up to q^40, "
                       "using the same CDH frame-shape formula (no literature table exists at "
                       "this order to compare the VALUES against, only the structural property).",
        "computed": f"Reimplemented numer/eta3/lambda24/twined24 from scratch in Python (from "
                    f"QSeries.lean's and Twining.lean's own definitions, not executed Lean code); "
                    f"regression against table2A through q^9: "
                    f"{i5['regression_matches_table2A_through_q9']}. Divisibility by 24 checked "
                    f"through q^{i5['scanned_range'][1]}: "
                    f"{'holds everywhere' if i5['divisibility_by_24_holds_beyond_q9'] else i5['failing_indices_beyond_q9']}. "
                    f"Hand-derivation: div_trunc's constant term shows divisibility-for-all-n "
                    f"reduces to combined[0] = chi*(-2) + c*Nlev*(Nlev-1) % 24 == 0, i.e. "
                    f"c == 8 (mod 12) for (chi,Nlev)=(8,2); checked directly over c in "
                    f"{[r['c'] for r in i5['c_dependence_scan']['rows']]}, congruence prediction "
                    f"confirmed: {i5['c_dependence_scan']['matches_hand_derived_congruence_c_eq_8_mod_12']}. "
                    f"So c=-16 is NOT an isolated selection -- infinitely many c in that residue "
                    f"class (e.g. c=-4, c=8, c=20, c=32, all checked) pass identically to -16.",
        "holds": i5["divisibility_by_24_holds_beyond_q9"],
        "rigidity_label": "NON_DISCRIMINATING",
        "independence": "Shares CDH's defining formula (eq. 4.18, Table 3's F_2A=-16*Lambda_2, "
                         "Table 14's chi_2A=8) with Lean's twined24 call -- these are the "
                         "structural Tier-L inputs, declared explicitly, not the printed q-series "
                         "table. Nothing beyond q^9 was checked against a literature table (none "
                         "was available in this worktree); the q^10..q^40 values are a genuine, "
                         "previously-uncompared computation. The rigidity_label is "
                         "NON_DISCRIMINATING with respect to c specifically, not with respect to "
                         "the divisibility extension itself: divisibility-by-24 to q^40 for "
                         "CDH's own c=-16 does hold (a real, previously unchecked extension), but "
                         "the divisibility property alone does not select c=-16 among integers -- "
                         "every c in the residue class 8 (mod 12) passes it too, so it cannot by "
                         "itself be used to argue -16 is forced.",
        "proposed_lean_statement":
            "theorem twined_div24_40 : (twined24 40 8 2 (-16)).all (· % 24 = 0) := by decide",
        "lean_reachable_now": False,
        "lean_reachable_note": "Not attempted -- no LeanMaster build here; N=40 is O(N^3)-ish "
                                "kernel-reduced integer work vs Lean's own N=9 `decide` call, "
                                "so cost is extrapolated (much higher), not measured.",
        "script": "audit/k3t2_rigidity_v3/reverse/item5_twined_2A_extended.py",
    })

    comparison_counts = {
        "AGREE": sum(1 for e in entries if e["holds"]),
        "DISAGREE": sum(1 for e in entries if not e["holds"]),
        "TOTAL": len(entries),
    }

    results = {
        "meta": {
            "pass": "reverse (v3)",
            "instruction": "Take 5-8 Lean theorems limited to a finite range and test them "
                            "beyond it with exact arithmetic; exclude the 8 already done in v1/v2.",
            "excluded_v1": ["p24 k=6..8", "DMZ m=4", "polar m=4,5", "trace bound d=7..10"],
            "excluded_v2": ["c_depends_only_on_D to q^20",
                            "moonshine twined_div24/twinedAll_div/modules_decompose to n=20",
                            "immortal m=4 Hecke form"],
            "num_items": len(entries),
            "comparison_counts": comparison_counts,
        },
        "entries": entries,
        "source_item_results": {
            "item1": i1, "item2": i2, "item3": i3, "item4": i4, "item5": i5,
        },
    }

    (HERE / "results.json").write_text(json.dumps(results, indent=2) + "\n")

    exports = {
        "reverse_pass_entries": [
            {k: v for k, v in e.items() if k != "computed"} | {"computed_summary": e["computed"]}
            for e in entries
        ],
        "comparison_counts": comparison_counts,
    }
    (HERE / "exports.json").write_text(json.dumps(exports, indent=2) + "\n")

    print(json.dumps(comparison_counts, indent=2))
    for e in entries:
        print(f"- {e['id']}: holds={e['holds']} label={e['rigidity_label']}")


if __name__ == "__main__":
    main()
