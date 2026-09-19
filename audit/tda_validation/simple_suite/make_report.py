#!/usr/bin/env python3
"""Build report.json (narrative + table) from results.json. Every number below is
read from results.json, which simple_suite.py --aggregate builds from cases/*.json.
Tier of every observed number: X (numerics). Expected topology: tier L."""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(HERE, "results.json")))
C = R["cases"]


def f(x, n=3):
    return ("%%.%dg" % n) % x if isinstance(x, (int, float)) else str(x)


def pd(case, path="per_dim_pipeline"):
    return C[case][path]


table = []
for cid in ("P1", "P2", "P3a", "P3b", "P4", "P10a", "P10b"):
    c = C[cid]
    table.append({"id": cid, "expected": c["expected_betti"], "observed_pipeline": c["observed_betti_pipeline"],
                  "pass": c["pass"], "function": c["function"],
                  "dominance_rules": {k: v["rule"] for k, v in c["per_dim_pipeline"].items()},
                  "p_values_vs_matched_poisson": c.get("p_values"),
                  "rips_reference_observed": c["observed_betti_rips_reference"], "rips_reference_pass": c["reference_pass"],
                  "wall_sec": c["wall_sec_total"], "peak_rss_mb": c["peak_rss_mb"]})
for cid in ("P5", "P6"):
    c = C[cid]
    for vn, v in c["variants"].items():
        row = {"id": "%s[%s]" % (cid, vn), "expected": c["expected_betti"], "pass": v["pass"], "parts": {}}
        for part, pr in v["parts"].items():
            row["parts"][part] = {k: (vv["observed"] if isinstance(vv, dict) and "observed" in vv else vv)
                                  for k, vv in pr.items() if k.startswith(("pipeline_Z2", "direct_alpha_Z", "rips_Z", "status", "error", "wall_sec_total", "peak_rss_mb"))}
        if v.get("killed_runs"):
            row["killed_runs"] = [{"status": k["status"], "exit_code": k["exit_code"], "command": k["command"]} for k in v["killed_runs"]]
        table.append(row)
p7 = C["P7"]
table.append({"id": "P7", "expected": [1, 3, 3, 1], "pass": p7["pass"], "function": p7["function"],
              "sweep": [{"N": s["N"], "observed": s["observed"], "pass": s["pass"],
                         "rules": {k: v["rule"] for k, v in (s["per_dim"] or {}).items()},
                         "wall_sec": s["wall_sec_total"], "peak_rss_mb": s["peak_rss_mb"],
                         "simplices_after_expansion": (s["info"] or {}).get("expanded_simplices")} for s in p7["sweep"]],
              "smallest_N_recovered": p7["smallest_N_recovered"],
              "post_hoc_landmarks": [{"n_landmarks": s["n_landmarks"], "pool": s["pool"], "observed": s["observed"], "pass": s["pass"],
                                      "rules": {k: v["rule"] for k, v in s["per_dim"].items()},
                                      "wall_sec": s["wall_sec_total"], "peak_rss_mb": s["peak_rss_mb"]} for s in p7.get("post_hoc_landmarks", [])]})
s42 = C["P8"]["seed42"]
table.append({"id": "P8", "expected": "6 H2 bars, death within 5% of 15", "pass": C["P8"]["pass"], "function": C["P8"]["function"],
              "top8_H2_deaths": [b[1] for b in s42["top8_H2_bars_birth_death"]],
              "death_rel_err_top6": s42["death_rel_err_top6"], "ratio_p6_over_p7_not_gated": s42["ratio_p6_over_p7"],
              "robustness_seeds_passing_of_5": C["P8"]["robustness_pass_count_of_5"],
              "wall_sec": C["P8"]["wall_sec_total"], "peak_rss_mb": C["P8"]["peak_rss_mb"]})
for cid in ("P9a", "P9b", "P9c"):
    c = C[cid]
    table.append({"id": cid, "expected": "no dominant bar (>= 48/50 seeds)", "pass": c["pass"], "function": c["function"],
                  "N": c["N"], "matched_to": c["matched_to"], "seeds_without_dominant_bar": c["n_seeds_no_dominant_bar"],
                  "max_persistence_H1": c["max_persistence_H1"], "max_persistence_H2": c["max_persistence_H2"],
                  "max_ratio_p1_p2_H1": c["ratio_p1_p2_H1"]["max"], "max_ratio_p1_p2_H2": c["ratio_p1_p2_H2"]["max"],
                  "wall_sec": c["wall_sec_total"], "peak_rss_mb": c["peak_rss_mb"]})
F1 = C["F1"]
table.append({"id": "F1", "expected": 7, "pass": F1["pass"],
              "cubical_H0_bars": F1["cubical"]["n_H0_bars_positive_persistence"], "direct_local_minima": F1["direct_local_minima_8nbr"],
              "pipeline_max_nu_b0": F1["pipeline"]["max_nu_b0"], "pipeline_b0_b1_at_top": [F1["pipeline"]["b0_at_max_nu"], F1["pipeline"]["b1_at_max_nu"]],
              "function": [F1["cubical"]["function"], F1["pipeline"]["function"]], "wall_sec": F1["wall_sec_total"], "peak_rss_mb": F1["peak_rss_mb"]})
FS = C["F1shuf"]
table.append({"id": "F1shuf (post hoc, TDA-N3)", "expected": "!= 7", "pass": FS["pass"],
              "cubical_H0_bars": FS["cubical"]["n_H0_bars_positive_persistence"], "direct_local_minima": FS["direct_local_minima_8nbr"],
              "pipeline_max_nu_b0": FS["pipeline"]["max_nu_b0"]})
table.append({"id": "F2", "expected": "binomial(n,k)", "pass": C["F2"]["pass"], "function": C["F2"]["function"],
              "tori": {k: {"infinite_bars": v["infinite_bars_per_dim"], "nonperiodic_control": v["nonperiodic_control_betti"]} for k, v in C["F2"]["tori"].items()}})
F3 = C["F3"]
table.append({"id": "F3", "expected": "uniform p-values (KS p >= 0.01, #{p<0.05} <= 9/100) for b0, b1, chi", "pass": F3["pass"],
              "function": F3["function"],
              "hartlap": {k: {"ks_p": v["hartlap_p"]["ks_p"], "n_below_0.05": v["hartlap_p"]["n_below_0.05"], "gate_pass": v["hartlap_p"]["gate_pass"],
                              "literal_TDA7_le5": v["hartlap_p"]["literal_TDA7_rule_le5"]} for k, v in F3["statistics"].items()},
              "empirical_rank": {k: v["empirical_rank_p"] for k, v in F3["statistics"].items()},
              "topology_diagnostic_betti_full_sky": F3["topology_diagnostic"]["betti_numbers"]})
table.append({"id": "N1", "expected": "no dominant H1/H2", "pass": C["N1"]["pass"], "H1_p1_over_p2": C["N1"]["H1"]["ratio"], "H2_p1_over_p2": C["N1"]["H2"]["ratio"]})
table.append({"id": "N2", "expected": "P3a H1 destroyed", "pass": C["N2"]["pass"], "H1_p2_over_p3": C["N2"]["H1_p2_over_p3"],
              "p_value_second_H1_bar_vs_P9c": C["N2"]["p_value_second_H1_bar_vs_P9c"]})
table.append({"id": "C0", "expected": "collapse_edges preserves the diagram", "pass": C["C0"]["pass"], "bottleneck_by_dim": C["C0"]["bottleneck_by_dim"]})

b1 = F3["statistics"]["b1"]
b1d = b1["post_hoc_bin_diagnostic"]
p6 = C["P6"]["variants"]
p7s = {s["N"]: s for s in p7["sweep"]}
lm = {s["n_landmarks"]: s for s in p7.get("post_hoc_landmarks", [])}
p9c_max = C["P9c"]["max_persistence_H1"]["max"]

def top_equal(parts):
    a = parts["pipeline"]["pipeline_Z2"]["per_dim"]; b = parts["direct"]["direct_alpha_Z/2"]["per_dim"]
    # the essential H0 bar is capped at the largest finite death, which for the direct run also
    # includes H3; compare bar counts, the top-6 H1/H2 persistences, and the finite H0 bars
    return (all(a[k]["n_bars"] == b[k]["n_bars"] for k in ("0", "1", "2"))
            and all(a[k]["top_persistence"] == b[k]["top_persistence"] for k in ("1", "2"))
            and a["0"]["top_persistence"][1:] == b["0"]["top_persistence"][1:])


same_z2 = {"P5 N=8000": top_equal(C["P5"]["variants"]["preregistered"]["parts"]), "P6 N=2000": top_equal(p6["N2000"]["parts"])}

findings_pipeline = [
    {"function": "cosmic_web_tda_scaled.alpha_persistence", "kind": "capability gap (read from source, confirmed by use)",
     "finding": "homology_coeff_field=2 is hard-coded, so the coefficient-field checks (P5, P6 over Z/3) cannot run through it; they ran through gudhi.AlphaComplex directly on the same points. Over Z/2 the pipeline and the direct run have the same bar counts in dims 0-2 and identical top-6 H1/H2 persistences and finite H0 bars (only the cap of the essential H0 bar differs, because the direct run also sees H3): %s." % same_z2,
     "consequence": "a torsion-sensitive result (anything that differs between Z/2 and Z/3) cannot be produced by the pipeline as written"},
    {"function": "cosmic_web_tda_scaled.alpha_persistence / betti_curve / euler_curve", "kind": "capability gap",
     "finding": "by_dim keeps dimensions 0, 1, 2 only; H3 and above are dropped without warning. P7 (beta_3 = 1) had to use gudhi.RipsComplex directly.",
     "consequence": "any 4-D or higher analysis through this function silently loses H3"},
    {"function": "cmb_tda.coarse_stats (with betti_curves_from_topology on the default nu grid [-4, 4])", "kind": "calibration FAILURE (F3)",
     "finding": "for the sublevel b1 curve the Hartlap chi2 p-values of 100 independent Gaussian maps against a 100-map ensemble are not uniform: KS p = %s, %d of 100 below 0.05. Cause (post-hoc diagnostic): the first coarse bin (nu in [-4, -3]) of b1 is exactly 0 in a fraction %s of the 100 ensemble maps (std %s), and the second is a rare-count bin (fraction %s of maps exactly 0, std %s); coarse_stats keeps df = N_BINS = 8 and a 1e-8 ridge, so one degenerate bin shifts the chi2(8) reference and a single loop in the rare-count bin gives a single-bin |z| up to %s. The smallest Hartlap p-values are %s. The empirical rank p-value is closer to uniform for b1 (KS p = %s) but has %d of 100 below 0.05."
               % (f(b1["hartlap_p"]["ks_p"]), b1["hartlap_p"]["n_below_0.05"], f(b1d["coarse_bin_fraction_of_sims_exactly_zero"][0]),
                  f(b1d["coarse_bin_ensemble_std"][0]), f(b1d["coarse_bin_fraction_of_sims_exactly_zero"][1]), f(b1d["coarse_bin_ensemble_std"][1]),
                  f(b1d["test_max_abs_z_single_bin_quantiles"][-1]), [round(x, 4) for x in b1d["test_hartlap_p_sorted_first10"][:3]],
                  f(b1["empirical_rank_p"]["ks_p"]), b1["empirical_rank_p"]["n_below_0.05"]),
     "not_affected": "b0 (KS p = %s) and chi (KS p = %s) pass the pre-registered gate" % (f(F3["statistics"]["b0"]["hartlap_p"]["ks_p"]), f(F3["statistics"]["chi"]["hartlap_p"]["ks_p"])),
     "consequence": "b1 chi2 p-values from the CMB TDA run (same coarse_stats, same nu grid, nside 128 masked) should not be read at face value; drop zero-variance bins (and set df to the number kept), restrict the nu range to where the curve varies, or use the rank p-value with its own calibration"},
    {"function": "cmb_tda.build_topology", "kind": "construction defect with no effect on b0/b1 (predicted in expectations.json before running)",
     "finding": "the full-sky nside-64 complex (V=%d, E=%d, T=%d) has Betti numbers %s instead of (1, 0, 1) for a triangulated S2: every set of 4 mutually 8-adjacent pixels around a pixel corner gives 4 triangles and no tetrahedron (a hollow tetrahedron)."
               % (F3["topology_diagnostic"]["n_vertices"], F3["topology_diagnostic"]["n_edges"], F3["topology_diagnostic"]["n_triangles"], F3["topology_diagnostic"]["betti_numbers"]),
     "consequence": "b0 and b1 are unaffected (F3 topology check b0 = 1, b1 = 0 passes). The curve the pipeline labels 'euler_chi' is b0 - b1, not the Euler characteristic of the complex it builds (b2 is ignored), and persistence spends work on tens of thousands of spurious H2 classes"},
]

findings_other = [
    {"case": "P6 (RP2), pre-registered design", "result": "FAIL (over budget)",
     "cause": "alpha in R4 on the (xy, xz, y^2 - z^2, 2yz) embedding does not fit the 6 GiB / 590 s cap at N = 6000 (the recorded run hit the 590 s timeout, exit 124; an earlier identical attempt, before error capture was added, ended in std::bad_alloc) nor at N = 3000 (std::bad_alloc, peak RSS %s MB). This is a Delaunay-size limit of gudhi/CGAL in R4 for this point set, not pipeline logic."
              % f(p6["N3000"]["parts"]["pipeline"]["peak_rss_mb"], 4),
     "post_hoc": "N = 2000 (first 2000 of the same draw): pipeline Z/2 %s, direct alpha Z/2 %s and Z/3 %s, all PASS (%s simplices, peak RSS %s MB). N = 1000: pipeline Z/2 observed %s, FAIL (H1 rule %s)."
                 % (p6["N2000"]["parts"]["pipeline"]["pipeline_Z2"]["observed"], p6["N2000"]["parts"]["direct"]["direct_alpha_Z/2"]["observed"],
                    p6["N2000"]["parts"]["direct"]["direct_alpha_Z/3"]["observed"], p6["N2000"]["parts"]["direct"]["direct_alpha_n_simplices"],
                    f(p6["N2000"]["parts"]["direct"]["peak_rss_mb"], 4), p6["N1000"]["parts"]["pipeline"]["pipeline_Z2"]["observed"],
                    p6["N1000"]["parts"]["pipeline"]["pipeline_Z2"]["per_dim"]["1"]["rule"])},
    {"case": "P7 (T3 in R6)", "result": "not recovered within budget",
     "sweep": {n: {"observed": s["observed"], "rules": {k: v["rule"] for k, v in (s["per_dim"] or {}).items()}, "wall_sec": s["wall_sec_total"]} for n, s in p7s.items()},
     "reading": "the three H1 bars, three H2 bars and the H3 bar are the three to four longest bars from N = 2000 up, but the 5x rule needs the noise bars to shrink below 1/5 of the capped true bars. At N = 16000 (%s s, %s MB) H1 (%s) and H3 (%s) pass; H0 (%s) and H2 (%s) do not. N = 32000 was not attempted: wall time grew from %s s (N = 8000) to %s s (N = 16000), so N = 32000 is projected far beyond the 590 s cap."
                % (f(p7s[16000]["wall_sec_total"]), f(p7s[16000]["peak_rss_mb"], 4), p7s[16000]["per_dim"]["1"]["rule"], p7s[16000]["per_dim"]["3"]["rule"],
                   p7s[16000]["per_dim"]["0"]["rule"], p7s[16000]["per_dim"]["2"]["rule"], f(p7s[8000]["wall_sec_total"]), f(p7s[16000]["wall_sec_total"])),
     "post_hoc_landmarks": "farthest-point landmarks from a 48000-point pool (not pre-registered): 8000 landmarks give %s (H0 %s); 16000 landmarks pass every 5x rule (H0 %s, H1 %s, H2 %s, H3 %s) but observed beta_0 = %d because that many H0 bars (the second is %s, near the landmark spacing) are >= 0.2 x P_all = %s; so FAIL under the pre-registered observed-beta rule."
                          % (lm[8000]["observed"], lm[8000]["per_dim"]["0"]["rule"], lm[16000]["per_dim"]["0"]["rule"], lm[16000]["per_dim"]["1"]["rule"],
                             lm[16000]["per_dim"]["2"]["rule"], lm[16000]["per_dim"]["3"]["rule"], lm[16000]["observed"][0],
                             f(lm[16000]["per_dim"]["0"]["top_persistence"][1]), f(0.2 * max(lm[16000]["per_dim"][k]["top_persistence"][0] for k in ("1", "2", "3"))))},
    {"case": "N2 (shuffled torus)", "result": "FAIL of the p-value half; the ratio half passes",
     "cause": "shuffling the columns destroys the two dominant H1 bars (p2/p3 = %s) but keeps the torus marginals, so the cloud is structureless but NOT uniform in its bounding box. Its 2nd H1 bar (%s) exceeds the maximum H1 persistence of all 50 uniform-box nulls (P9c max %s), p = %s. The pre-registered uniform null is mismatched for a non-uniform density; this is a null-design lesson (match the marginals or density), not a pipeline defect."
              % (f(C["N2"]["H1_p2_over_p3"]), f(C["N2"]["H1_top5"][1] if "H1_top5" in C["N2"] else float("nan")), f(p9c_max), f(C["N2"]["p_value_second_H1_bar_vs_P9c"]))},
    {"case": "Rips reference paths (P2, P3a, P3b, P5, P10a, P10b)", "result": "FAIL on the pre-registered subsample sizes; P1, P4 and P6 Z/3 pass; P6 Z/2 fails the H1 rule (%s)" % p6["preregistered"]["parts"]["rips"]["rips_Z/2"]["per_dim"]["1"]["rule"],
     "cause": "subsamples of 400-1000 points are too sparse: the longest H0 bar beyond the expected components (nearest-neighbour connection scale) is %s, which exceeds 0.2 x P_all, and in several cases the noise H1 bars are within 5x of the true ones (rules: %s). Pipeline alpha on the full clouds passes the same cases. These are sample-size failures of the reference design, not gudhi or pipeline defects."
              % ({k: f(C[k]["per_dim_reference"]["0"]["top_persistence"][C[k]["expected_betti"][0]]) for k in ("P2", "P3a", "P3b", "P10a", "P10b")},
                 {k: C[k]["per_dim_reference"]["1"]["rule"] for k in ("P3a", "P3b", "P10b")})},
]

margins = {
    "P3a_H1": pd("P3a")["1"]["rule"], "P10a_H1_zero_rule": pd("P10a")["1"]["rule"],
    "P5_Z3_H1": C["P5"]["variants"]["preregistered"]["parts"]["direct"]["direct_alpha_Z/3"]["per_dim"]["1"]["rule"],
    "P6_N2000_Z2_H1": p6["N2000"]["parts"]["pipeline"]["pipeline_Z2"]["per_dim"]["1"]["rule"],
    "P8_ratio_p6_over_p7_not_gated": s42["ratio_p6_over_p7"],
    "note": "these cases pass with modest margins against the 5x or 0.2 thresholds",
}

loads = []
for fn in os.listdir(os.path.join(HERE, "cases")):
    d = json.load(open(os.path.join(HERE, "cases", fn)))
    if isinstance(d.get("loadavg_start"), list):
        loads.append(d["loadavg_start"][0])

report = {
    "title": "TDA known-answer suite on simple spaces (LeanFlow TEST_USE_CASES section 2)",
    "generated_by": "make_report.py from results.json",
    "tiers": "expected topology: L (Hatcher 2002; Kunneth; construction). Every observed number: X (numerics). Nothing here is proved.",
    "hardware": "Intel Xeon @ 2.20GHz, 8 cores, 29 GB RAM, shared; 1-min load average at case start ranged %s-%s over the recorded runs (per-case loadavg_start in results.json). Every process: prlimit --as=6442450944, timeout <= 590 s, OMP_NUM_THREADS=1. Wall times are under contention and are not benchmarks." % (f(min(loads)), f(max(loads))),
    "code_under_test_sha256": R["code_under_test_sha256"],
    "summary": R["summary"],
    "table": table,
    "pipeline_function_findings": findings_pipeline,
    "other_failures": findings_other,
    "narrow_margins": margins,
    "p_value_note": "P1, P2, P3a p-values compare the weakest expected bar to the MAX bar of each of 50 matched uniform-box nulls; the minimum attainable p is 1/51 = 0.0196, and all three reached it: %s, %s, %s."
                    % (C["P1"]["p_values"], C["P2"]["p_values"], C["P3a"]["p_values"]),
    "deviations_from_expectations": [
        "evaluate() first lacked the pre-registered clause 'if the case expects no bar in any dim >= 1, use p_1 < 5 p_2'; it was added before any affected case was aggregated and P6 Rips was rerun (only P6 over Z/3 expects no bar in dims >= 1).",
        "P6 post-hoc N = 1000, 2000 (and a failed N = 3000) after N = 6000 exceeded the budget; the pre-registered result stays FAIL.",
        "P7 extended past the pre-registered N list (8000, 16000), following the task's 'N as large as fits the budget'; post-hoc landmark runs P7lm are labelled and do not change the P7 result.",
        "F1shuf (shuffled-site control) added post hoc because TEST_USE_CASES gained TDA-N3 during this run.",
        "F3 per-bin diagnostic added post hoc after the b1 KS failure; it does not change the gate.",
        "The P5 direct run frees the AlphaComplex object before persistence (memory); the complex is unchanged.",
        "All cheap and medium cases were rerun at the committed script (9da5954); outputs matched the earlier runs exactly apart from timing fields. P6 N-variants, P7 and F3 chunks were run with the pre-commit working copy (same code path for those cases).",
    ],
    "killed_or_crashed_runs": R.get("killed_or_crashed_runs", []),
}
txt = json.dumps(report, indent=1, default=float)
# rule strings carry the threshold; make it unambiguous that it is a requirement, not a claim
txt = txt.replace(" >= 5", " (required >= 5)").replace(" < 0.2", " (required < 0.2)").replace(" < 5\"", " (required < 5)\"")
json.loads(txt)
open(os.path.join(HERE, "report.json"), "w").write(txt)
print(json.dumps(R["summary"]))
