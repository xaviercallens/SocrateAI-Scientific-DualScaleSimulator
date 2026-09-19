#!/usr/bin/env python3
"""
EXPERIMENT E2 -- is the PTA sector GR? (no-data branch)

Ground truth checked THIS session (commands below, both re-run and
recorded verbatim in the JSON output):
  1. `find data -iname '*hellings*' -o -iname '*angular*'` over the whole
     worktree data/ tree -> empty.
  2. `ls data/real/pulsar_timing/` -> only nanograv_kde_freespectrum.zip
     (a ceffyl free-spectrum KDE, a function of FREQUENCY, not angle).
This matches and does not override the data manifest's own
"nanograv_15yr_hd_angular_correlation": "ABSENT" entry, which already
records 3+ fetch attempts and a URL trail (NANOGrav 15yr
stochastic-analysis repo figure_1 data, PTArcade, Zenodo searches) and
notes that a binned Gamma(theta) table exists only as an intermediate
product of the enterprise_extensions pipeline, computed from a ~GB raw
dataset in a ~18 min run -- reconstructing it here would violate the
evidence-bound ground rule (we did not run that pipeline).

Because no angular-correlation dataset is in hand, we do NOT fit
Gamma(theta) = A*(HD(theta) + c4*P4(theta)) to anything. That fit is only
meaningful once a real binned correlation (with covariance) exists.

Instead this script does the one thing that is legitimate on a fixed,
already-defined formula with NO data: it uses
scripts/param_loop_sim.compute_pta_observable (imported, not
reimplemented) to state a *procurement spec* -- exactly how sensitive a
future angular table would need to be to move c4_pta_product off zero at
several thresholds. This is tier X/B (exact arithmetic on a fixed,
pre-existing formula), not a fit, and not a bound on c4 itself.
"""
import json
import math
import os
import subprocess
import sys

REPO_ROOT = "/mnt/disks/disk-socrateai-local-1/dualscale-wt-reverse"
sys.path.insert(0, os.path.join(REPO_ROOT, "scripts"))
from param_loop_sim import compute_pta_observable  # noqa: E402

OUT_DIR = os.path.dirname(os.path.abspath(__file__))


def grep_for_angular_data():
    checks = {}
    find_cmd = ["find", "data", "-iname", "*hellings*", "-o", "-iname", "*angular*"]
    r = subprocess.run(find_cmd, cwd=REPO_ROOT, stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE, text=True, check=False)
    checks["find_hellings_or_angular_in_data"] = {
        "command": "cd " + REPO_ROOT + " && " + " ".join(find_cmd),
        "stdout": r.stdout.strip(),
        "stderr": r.stderr.strip(),
    }
    ls_cmd = ["ls", "-la", "data/real/pulsar_timing/"]
    r2 = subprocess.run(ls_cmd, cwd=REPO_ROOT, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, text=True, check=False)
    checks["ls_pulsar_timing_dir"] = {
        "command": "cd " + REPO_ROOT + " && " + " ".join(ls_cmd),
        "stdout": r2.stdout.strip(),
        "stderr": r2.stderr.strip(),
    }
    return checks


def sensitivity_scan():
    """
    For a grid of c4_pta_product values (with pta_suppression=1.0 so
    c4_c0_ratio == c4_pta_product directly, per compute_pta_observable's
    own definition c4_pta_product = c4_c0_ratio * pta_suppression),
    report max_deviation_from_hd across the harness's fixed bins, and the
    fractional deviation relative to the HD curve's own dynamic range.
    This is EXACT arithmetic on the existing formula -- no fit, no data.
    """
    grid = [1e-4, 1e-3, 1e-2, 3e-2, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 1.0]
    rows = []
    # HD curve dynamic range (reference), computed once at c4=0
    ref = compute_pta_observable(pta_suppression=1.0, c4_c0_ratio=0.0)
    hd = ref["hd_curve"]
    hd_range = max(hd) - min(hd)
    for c4 in grid:
        res = compute_pta_observable(pta_suppression=1.0, c4_c0_ratio=c4)
        frac_of_hd_range = res["max_deviation_from_hd"] / hd_range
        rows.append({
            "c4_pta_product": c4,
            "max_deviation_from_hd": res["max_deviation_from_hd"],
            "hd_curve_dynamic_range": hd_range,
            "fractional_perturbation_of_hd_range": frac_of_hd_range,
        })
    # Invert (linearly -- compute_pta_observable is linear in c4_c0_ratio
    # at fixed pta_suppression, verified below) for round-number thresholds
    # of fractional perturbation: 1%, 5%, 10%.
    # Linearity check: max_deviation_from_hd(c4) / c4 should be constant.
    ratios = [r["max_deviation_from_hd"] / r["c4_pta_product"] for r in rows if r["c4_pta_product"] > 0]
    linearity_spread = (max(ratios) - min(ratios)) / (sum(ratios) / len(ratios))
    slope = sum(ratios) / len(ratios)
    thresholds = {}
    for frac in (0.01, 0.05, 0.10):
        target_dev = frac * hd_range
        c4_needed = target_dev / slope
        thresholds[f"c4_for_{int(frac*100)}pct_of_HD_dynamic_range"] = c4_needed
    return {
        "grid_rows": rows,
        "linearity_check": {
            "note": "compute_pta_observable's gamma_theta is exactly linear "
                     "in c4_c0_ratio at fixed pta_suppression (l4_response "
                     "does not depend on c4); slope = max_deviation_from_hd / c4",
            "slope_max_dev_per_unit_c4": slope,
            "fractional_spread_of_slope_across_grid": linearity_spread,
        },
        "procurement_thresholds": thresholds,
    }


def main():
    out = {
        "experiment": "E2_pta_sector_GR_test",
        "branch": "NO_DATA",
        "ground_rule_evidence": grep_for_angular_data(),
        "manifest_absent_entry_quoted": (
            "nanograv_15yr_hd_angular_correlation status=ABSENT: "
            "'3+ sources checked (nanograv/15yr_stochastic_analysis "
            "figure_1, PTArcade, Zenodo searches). Only posterior chains "
            "and a pipeline notebook exist that COMPUTES the binned "
            "correlation from a ~GB raw dataset (~18min enterprise_extensions "
            "run); not reconstructed per ground rules.'"
        ),
        "preregistration_line_quoted": (
            "PRE_REGISTRATION.md:65 -- 'c4_pta_product (l=4 term in the PTA "
            "angular correlation) | bound it; a measured correlation "
            "consistent with pure Hellings-Downs would exclude large values "
            "| NANOGrav, IPTA, then SKA'"
        ),
        "what_would_decide_it": {
            "required_product": "A binned angular two-point pulsar-timing "
                "correlation table: bin centers in angular separation theta "
                "(degrees), the correlation estimate rho(theta) or Gamma(theta) "
                "per bin, AND the bin-to-bin covariance matrix (HD-curve "
                "estimators from a fixed pulsar array are correlated across "
                "angular bins; diagonal errors alone would understate chi2 "
                "and bias the c4 interval).",
            "named_candidates": [
                "NANOGrav 15yr stochastic-analysis release: the exact table "
                "underlying its published Hellings-Downs figure (rho vs "
                "angle with uncertainties/covariance), if released as data "
                "rather than only as a plotted figure -- checked, not found "
                "as of this session (see manifest URL trail).",
                "IPTA DR3 (or DR2) combined angular-correlation data product, "
                "if/when released with covariance.",
                "Running the enterprise_extensions binned-correlation "
                "pipeline on the full NANOGrav 15yr raw dataset ourselves "
                "(~1 GB, ~18 min per the manifest) -- NOT done here because "
                "it is a multi-step derived computation, not a fetch of an "
                "existing public product; flagged as a concrete follow-up.",
            ],
            "sufficiency_note": "The NANOGrav/IPTA table would be sufficient "
                "ALONE only if published with its covariance; a figure-only "
                "release (points with error bars, no covariance) is usable "
                "for a rough diagonal-only chi2 but not for a defensible "
                "68/95% interval on c4_pta_product.",
        },
        "sensitivity_scan": sensitivity_scan(),
        "parameter_effect": (
            "none -- free-parameter count stays at 2 (mu_sym, "
            "c4_pta_product). No angular-correlation dataset exists to fit, "
            "so c4_pta_product is neither constrained nor removed by data "
            "this round. Setting c4_pta_product = 0 in the M0 hypothesis is "
            "a HYPOTHESIS CHANGE (the PTA sector set to its GR/pure-"
            "Hellings-Downs value), tier L for GR itself, NOT a derivation "
            "from K3xT2 and NOT a data-driven reduction. LeanMaster Streams "
            "6-8 (E2-E4) confirm mu_sym and c4_pta_product are unchanged by "
            "every result reviewed there."
        ),
    }
    out_path = os.path.join(OUT_DIR, "e2_pta_result.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"wrote {out_path}")
    print(json.dumps(out["sensitivity_scan"]["procurement_thresholds"], indent=2))


if __name__ == "__main__":
    main()
