#!/usr/bin/env python
"""
Combine a merged null report (merge_null.py output) with one or more
sensitivity-scan reports (cmb_tda.py --skip-null output, one per Gmu or per
Gmu-batch) into the final E5 CMB-TDA report.

Usage:
  .venv-tda/bin/python merge_final.py \
      --null-report null_report.json \
      --sens sens_report_a.json sens_report_b.json ... \
      --commands commands.json \
      --out e5_cmb_tda_report.json
"""
import argparse
import json


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--null-report", required=True)
    ap.add_argument("--sens", nargs="+", required=True)
    ap.add_argument("--commands", required=True, help="JSON list of the exact commands run, for provenance")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    with open(args.null_report) as f:
        null_report = json.load(f)

    gmu_scan = {}
    null_ref_stats = None
    model = None
    stat_used = None
    bound_comparison = None
    for sp in args.sens:
        with open(sp) as f:
            s = json.load(f)
        gmu_scan.update(s["cosmic_string_sensitivity"]["gmu_scan"])
        null_ref_stats = s["null_reference_stats_euler_chi_sublevel"]
        model = s["cosmic_string_sensitivity"]["model"]
        stat_used = s["cosmic_string_sensitivity"]["statistic_used"]
        bound_comparison = s["cosmic_string_sensitivity"]["planck_gmu_bound_comparison"]

    with open(args.commands) as f:
        commands = json.load(f)

    # Multiple-testing note: 6 statistics were tested (b0,b1,euler_chi x
    # sublevel,superlevel). A single p-value below 0.05 among 6 is not, by
    # itself, evidence of anything -- report the Bonferroni-corrected p for
    # the smallest raw p found, and describe whether the standardized
    # residual curve for that statistic is a broadband offset (spectrum
    # mismatch) or a localized, sign-changing excursion (the only kind that
    # would be interesting).
    all_stats = []
    for direction in ("statistic_sublevel", "statistic_superlevel"):
        for curve in ("b0", "b1", "euler_chi"):
            s = null_report[direction][curve]
            all_stats.append((f"{direction}.{curve}", s["p_value_chi2_survival"], s))
    all_stats_valid = [(name, p, s) for name, p, s in all_stats if p is not None]
    n_tests = len(all_stats_valid)
    smallest = min(all_stats_valid, key=lambda t: t[1]) if all_stats_valid else None
    multiple_testing_note = None
    if smallest is not None:
        name, p_raw, s = smallest
        p_bonf = min(1.0, p_raw * n_tests)
        resid = s["standardized_residual_full_curve"]
        resid_finite = [r for r in resid if r is not None and r == r]  # drop None/NaN
        n_pos = sum(1 for r in resid_finite if r > 0)
        n_neg = sum(1 for r in resid_finite if r < 0)
        sign_changing = n_pos > 0 and n_neg > 0
        multiple_testing_note = {
            "n_statistics_tested": n_tests,
            "smallest_raw_p_value": {"statistic": name, "p_raw": p_raw},
            "bonferroni_corrected_p": p_bonf,
            "residual_curve_shape": "sign-changing (localized excursion)" if sign_changing
                else "single-sign (consistent with a broadband amplitude/spectrum offset, not a localized topological feature)",
            "n_residual_points_positive": n_pos, "n_residual_points_negative": n_neg,
            "interpretation": "A single raw p-value below 0.05 out of "
                f"{n_tests} tested statistics is not evidence by itself "
                "(look-elsewhere effect); the Bonferroni-corrected p above is "
                "the honest number. The residual-curve shape further "
                "distinguishes a real localized topological anomaly from a "
                "broadband pseudo-Cl/fsky spectrum-matching artifact.",
        }

    # smallest Gmu (of those scanned) at which >=95% of injected sims exceed
    # the null 95th-percentile threshold
    sorted_gmu = sorted(gmu_scan.keys(), key=lambda x: float(x))
    clearing = [g for g in sorted_gmu if gmu_scan[g]["fraction_exceeding_null_95th_pct_LOO"] >= 0.95]
    smallest_95pct_gmu = min(clearing, key=lambda x: float(x)) if clearing else None

    final = {
        "experiment": "E5 CMB TDA (TDA T2): WMAP 9yr ILC persistent homology vs Gaussian null, "
                       "plus toy Kaiser-Stebbins cosmic-string injection sensitivity",
        "tier": "X (exploratory numerics)",
        "framing": null_report["framing"],
        "script": "audit/reverse_zero/E5-cmb-tda/cmb_tda.py (v2 -- see its docstring for fixes vs the "
                   "N=10 smoke_test.json: ell=0,1 zeroed in the sim Cl, superlevel null added, "
                   "cached mask-only topology (perf only, verified unchanged results), LOO-based "
                   "sensitivity threshold, injected segments truncated to a genuine arc)",
        "smoke_test_reference": "smoke_test.json (N=10, pre-v2-fix script, kept for provenance; this "
                                  "report supersedes it as the full N>=100 run)",
        "commands_exact": commands,
        "null_test": null_report,
        "multiple_testing_note": multiple_testing_note,
        "cosmic_string_sensitivity": {
            "model": model,
            "statistic_used": stat_used,
            "planck_gmu_bound_comparison": bound_comparison,
            "gmu_scan_sorted_ascending": {g: gmu_scan[g] for g in sorted_gmu},
            "smallest_gmu_clearing_95pct_of_scanned_values": smallest_95pct_gmu,
            "interpretation": "This is a TDA-pipeline sensitivity threshold for a toy, non-network "
                "injection model on a degraded (nside=128) WMAP ILC map -- it is NOT a cosmological "
                "bound and is NOT claimed equal to or better than LeanMaster's imported "
                "planckGmuBound=1.5e-7. If the smallest-clearing Gmu found here is below 1.5e-7, "
                "that indicates the toy injection is still optimistic (an nside=128 WMAP-derived "
                "null should not out-sensitivity Planck-grade analyses), not new physics.",
        },
        "verdict": "M0 predicts Gaussian CMB statistics (no topological defects sourced by K3xT2 "
                   "or anything else). The null test above checks whether the real WMAP ILC map's "
                   "Betti/Euler-characteristic curves (sublevel AND superlevel) are consistent with "
                   "that Gaussian-null ensemble. Consistency (or inconsistency) here is a test of "
                   "M0's Gaussianity assumption ONLY -- it is NOT evidence for or against K3xT2, "
                   "which makes no independent, unit-bearing prediction about CMB topology beyond "
                   "the same Gaussian statistics any GR+LCDM sky would show (LeanMaster Streams 6-7: "
                   "every physical reading of the dual-scale idea tried so far is excluded; Streams "
                   "4-5 results are dimensionless). Free-parameter count: UNCHANGED (0 either way).",
        "parameter_effect": "0 (null test of M0's Gaussianity prediction; a pass or fail here adds "
                              "or removes no free parameter of the theory)",
    }
    with open(args.out, "w") as f:
        json.dump(final, f, indent=2)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
