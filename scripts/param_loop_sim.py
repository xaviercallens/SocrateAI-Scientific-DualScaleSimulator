#!/usr/bin/env python3
"""
param_loop_sim.py -- REQ-LOOP-01: Parametrized experiment harness for the
Hypothesis -> Experimentation -> TDA -> reduce-parameter loop.

Takes a 6-vector of theory parameters (a_pot, b_pot, mu_sym, lambda_sym,
pta_suppression, c4_c0_ratio) and computes, from a single call into
workshopcosmo.py, three observable blocks:

  1. dark-energy observables  (quintessence ODE -> w(z), w0/wa CPL fit,
     H(z)/H0, dimensionless D_M*H0 "BAO-style" distance ratios, and a
     shape-only distance modulus, all on a z-grid and at six fixed z's)
  2. screening observable     (run_symmetron_screening_simulation)
  3. PTA angular-correlation observable Gamma(theta) on 15 bins
     (Hellings-Downs + l=4 hexadecapole term), replicating exactly the
     formula used in run_nanograv_hexadecapole_simulation but with
     pta_suppression wired in as a real function argument instead of the
     hard-coded local constant it is in that function today.

Evidence trail: every number in the output JSON comes from this script's
own deterministic computation (seed 42, no randomness is actually used --
the ODE/BVP solves are deterministic -- the seed is set defensively per
the loop ground rules) or from a direct call into workshopcosmo.py.

CLI:
    param_loop_sim.py --params '{"a_pot":1.0,...}' --out file.json
    param_loop_sim.py --selftest
"""

import argparse
import copy
import json
import math
import os
import subprocess
import sys
import time
from typing import Any, Dict, List

import numpy as np

np.random.seed(42)  # ground rule: seed every random call (seed 42); no RNG is actually
                     # invoked by the deterministic ODE/BVP solvers used below.

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, REPO_ROOT)
import workshopcosmo as wc  # noqa: E402

# -----------------------------------------------------------------------
# Defaults: exactly the hard-coded values workshopcosmo.py used before this
# harness existed (a_pot=1.0, b_pot=0.01 from compute_potential;
# mu_sym=1.0, lambda_sym=1.0 from run_symmetron_screening_simulation;
# pta_suppression=0.005, c4_c0_ratio=16.07 from run_nanograv_hexadecapole_simulation).
# -----------------------------------------------------------------------
DEFAULT_PARAMS: Dict[str, float] = {
    "a_pot": 1.0,
    "b_pot": 0.01,
    "mu_sym": 1.0,
    "lambda_sym": 1.0,
    "pta_suppression": 0.005,
    "c4_c0_ratio": 16.07,
}

Z_POINTS: List[float] = [0.3, 0.5, 0.7, 1.0, 1.5, 2.3]
Z_GRID = np.linspace(0.0, 2.3, 40)

# Fixed simulation resolution for the dark-energy ODE, held constant across
# parameter points for a fair comparison. NOTE (reported in "problems"):
# every dark-energy distance/modulus observable below is defined relative to
# a0 = a(t_max), i.e. "today" is a convention tied to t_max, which is a
# harness setting, not one of the six theory parameters.
QUINT_T_MAX = 70.0
QUINT_NUM_POINTS = 500

PTA_N_BINS = 15


# =============================================================================
# Block 1: dark-energy observables
# =============================================================================
def compute_dark_energy_observables(a_pot: float, b_pot: float) -> Dict[str, Any]:
    quint = wc.run_quintessence_simulation(
        t_max=QUINT_T_MAX, num_points=QUINT_NUM_POINTS, a_pot=a_pot, b_pot=b_pot
    )
    # Regression check against the pre-existing workshopcosmo code path (unchanged
    # by this harness). Its (0.2 <= a <= 1.0) mask assumes a is normalized to 1
    # "today"; here a grows unbounded (a0 = a(t_max) ~ 1e4-1e5), so that mask
    # actually selects an EARLY window of the trajectory, not a late-time one.
    # Kept only as an existing-behavior regression signal; see "problems".
    legacy_obs = wc.run_observables_analysis(quint)

    a_arr = np.asarray(quint["a"], dtype=np.float64)
    w_arr = np.asarray(quint["w_phi"], dtype=np.float64)
    H_arr = np.asarray(quint["H"], dtype=np.float64)

    a0 = float(a_arr[-1])  # convention: "today" = state at t_max
    z_arr = a0 / a_arr - 1.0

    order = np.argsort(z_arr)
    z_sorted = z_arr[order]
    w_sorted = w_arr[order]
    H_sorted = H_arr[order]
    a_norm_sorted = a_arr[order] / a0

    valid = z_sorted >= 0.0
    z_v = z_sorted[valid]
    w_v = w_sorted[valid]
    H_v = H_sorted[valid]
    a_norm_v = a_norm_sorted[valid]

    if len(z_v) < 5 or z_v.max() < 2.3:
        raise RuntimeError(
            f"insufficient late-time sampling: {len(z_v)} points with z>=0, "
            f"max z={z_v.max() if len(z_v) else float('nan')} (need >=2.3)"
        )

    w_of_z_grid = np.interp(Z_GRID, z_v, w_v)
    w_of_z_points = np.interp(Z_POINTS, z_v, w_v)

    # Own late-time CPL fit: w(a) = w0 + wa*(1-a_norm), restricted to the
    # actual "today" window z in [0, 2.3], unlike the legacy mask above.
    fit_mask = z_v <= 2.3
    x_feat = 1.0 - a_norm_v[fit_mask]
    w_feat = w_v[fit_mask]
    poly = np.polyfit(x_feat, w_feat, 1)
    wa_cpl_latetime = float(poly[0])
    w0_cpl_latetime = float(poly[1])

    # H(z)/H0, H0 := H at z=0 (interpolated)
    H0 = float(np.interp(0.0, z_v, H_v))
    Hz_H0_v = H_v / H0
    Hz_H0_grid = np.interp(Z_GRID, z_v, Hz_H0_v)
    Hz_H0_points = np.interp(Z_POINTS, z_v, Hz_H0_v)

    # Dimensionless comoving distance H0*D_M(z) = integral_0^z dz'/(H(z')/H0)
    inv_Hz_H0 = 1.0 / Hz_H0_v
    seg = 0.5 * (inv_Hz_H0[1:] + inv_Hz_H0[:-1]) * np.diff(z_v)
    cum = np.concatenate(([0.0], np.cumsum(seg)))
    DM_H0_grid = np.interp(Z_GRID, z_v, cum)
    DM_H0_points = np.interp(Z_POINTS, z_v, cum)

    # Shape-only distance modulus: mu_shape(z) = 5*log10(D_L*H0), D_L = (1+z)*D_M.
    # The additive absolute offset (5*log10(c/H0) in physical units) is
    # UNDETERMINED here because this model has no physical calibration of H0;
    # only the z-shape is meaningful for a later BAO/SN comparison.
    DL_H0_grid = (1.0 + Z_GRID) * DM_H0_grid
    DL_H0_points = (1.0 + np.asarray(Z_POINTS)) * DM_H0_points
    mu_shape_grid = 5.0 * np.log10(np.clip(DL_H0_grid, 1e-12, None))
    mu_shape_points = 5.0 * np.log10(np.clip(DL_H0_points, 1e-12, None))

    return {
        "success": bool(quint["success"]),
        "a0_today_convention": a0,
        "z_grid": Z_GRID.tolist(),
        "w_of_z_grid": w_of_z_grid.tolist(),
        "z_points": list(Z_POINTS),
        "w_of_z_points": w_of_z_points.tolist(),
        "w0_cpl_latetime": w0_cpl_latetime,
        "wa_cpl_latetime": wa_cpl_latetime,
        "H_of_z_over_H0_grid": Hz_H0_grid.tolist(),
        "H_of_z_over_H0_points": Hz_H0_points.tolist(),
        "D_M_times_H0_grid": DM_H0_grid.tolist(),
        "D_M_times_H0_points": DM_H0_points.tolist(),
        "distance_modulus_shape_grid": mu_shape_grid.tolist(),
        "distance_modulus_shape_points": mu_shape_points.tolist(),
        "legacy_w0_fit": legacy_obs["w0_fit"],
        "legacy_wa_fit": legacy_obs["wa_fit"],
        "legacy_late_time_w_current": legacy_obs["late_time_w_current"],
    }


# =============================================================================
# Block 2: screening observable
# =============================================================================
def compute_screening_observable(mu_sym: float, lambda_sym: float) -> Dict[str, Any]:
    res = wc.run_symmetron_screening_simulation(mu_sym=mu_sym, lambda_sym=lambda_sym)
    ssf = float(res["screening_suppression_factor"])
    pcr = float(res["phi_center_ratio"])
    # workshopcosmo's run_symmetron_screening_simulation always reports
    # success=True even when the BVP/relaxation fallback has actually
    # diverged to NaN/Inf (observed empirically at mu_sym=10.0, lambda_sym=1.0
    # -- see "problems"). Detect that here rather than trusting its flag.
    numerically_stable = bool(np.isfinite(ssf) and np.isfinite(pcr))
    return {
        "success": bool(res["success"]),
        "numerically_stable": numerically_stable,
        "screening_suppression_factor": ssf,
        "phi_center_ratio": pcr,
        "is_screened": bool(res["is_screened"]),
    }


# =============================================================================
# Block 3: PTA observable
# =============================================================================
def pta_suppression_grep_evidence() -> str:
    """
    Ground-rule requirement: if pta_suppression appears nowhere in any
    computation, report that fact verbatim via grep rather than inventing a
    role for it. Here it DOES appear (workshopcosmo.py line ~764), but only
    as a hard-coded local constant inside run_nanograv_hexadecapole_simulation,
    never as a function parameter -- this harness is what first exposes it as
    a real, externally-settable argument.
    """
    try:
        out = subprocess.run(
            ["grep", "-n", "pta_suppression", os.path.join(REPO_ROOT, "workshopcosmo.py")],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False,
        )
        return out.stdout.strip()
    except Exception as exc:  # pragma: no cover
        return f"ABSENT (grep failed: {exc})"


def compute_pta_observable(pta_suppression: float, c4_c0_ratio: float) -> Dict[str, Any]:
    """
    Gamma(theta) on PTA_N_BINS angular bins = Hellings-Downs curve + l=4
    hexadecapole term, reproducing EXACTLY the formula found at
    workshopcosmo.py:run_nanograv_hexadecapole_simulation (lines ~752-766):
        hd_curve  = 1.5*x*log(x) - 0.25*x + 0.5,  x = (1-cos theta)/2
        l4        = (35 cos^4 - 30 cos^2 + 3) / 8            (Legendre P4)
        mod_curve = hd_curve + c4_c0_ratio * pta_suppression * l4
    The only change from that function is (a) pta_suppression is a real
    argument instead of the hard-coded local 0.005, and (b) the +/-0.15
    "cosmic variance envelope" and is_hidden flag from that function are an
    uncited, invented bound and are deliberately NOT reproduced here.
    """
    theta = np.linspace(0.01, math.pi, PTA_N_BINS)
    cos_theta = np.cos(theta)
    x = (1.0 - cos_theta) / 2.0
    hd_curve = 1.5 * x * np.log(x + 1e-12) - 0.25 * x + 0.5
    l4_response = (35.0 * cos_theta ** 4 - 30.0 * cos_theta ** 2 + 3.0) / 8.0
    gamma_theta = hd_curve + c4_c0_ratio * pta_suppression * l4_response
    return {
        "n_bins": PTA_N_BINS,
        "theta_rad": theta.tolist(),
        "hd_curve": hd_curve.tolist(),
        "l4_response": l4_response.tolist(),
        "gamma_theta": gamma_theta.tolist(),
        "pta_suppression": pta_suppression,
        "c4_c0_ratio": c4_c0_ratio,
        "c4_pta_product": float(c4_c0_ratio * pta_suppression),
        "max_deviation_from_hd": float(np.max(np.abs(gamma_theta - hd_curve))),
        "grep_evidence_pta_suppression_workshopcosmo": pta_suppression_grep_evidence(),
    }


# =============================================================================
# Full point evaluation
# =============================================================================
def evaluate_point(params: Dict[str, float]) -> Dict[str, Any]:
    merged = dict(DEFAULT_PARAMS)
    merged.update(params)
    t0 = time.time()
    dark_energy = compute_dark_energy_observables(merged["a_pot"], merged["b_pot"])
    screening = compute_screening_observable(merged["mu_sym"], merged["lambda_sym"])
    pta = compute_pta_observable(merged["pta_suppression"], merged["c4_c0_ratio"])
    elapsed = time.time() - t0
    return {
        "params": merged,
        "dark_energy": dark_energy,
        "screening": screening,
        "pta": pta,
        "wall_time_seconds": elapsed,
    }


# =============================================================================
# Self-test
# =============================================================================
def _scalar_diff(a: Any, b: Any) -> float:
    return float(np.max(np.abs(np.asarray(a, dtype=np.float64) - np.asarray(b, dtype=np.float64))))


def _compare_observable(baseline_val: Any, varied_val: Any, rtol: float = 1e-6, atol: float = 1e-9) -> str:
    """
    Classifies a (baseline, x10-varied) observable pair as "changed",
    "unchanged", or "nonfinite". Uses a relative+absolute tolerance (not a
    bare 1e-9 absolute one) because some observables here span many orders
    of magnitude (e.g. screening_suppression_factor ~1e-4-1e-14), where a
    real, reproducible change can sit well under a fixed absolute epsilon.
    NaN/Inf in either side is reported as its own category rather than
    silently comparing equal (NaN == NaN is False in a naive np.abs diff,
    but a NaN vs NaN diff of 0 would otherwise be misread as "unchanged" --
    this happened during development for mu_sym x10, where the symmetron
    BVP fallback silently diverged to NaN; see "problems").
    """
    a = np.asarray(baseline_val, dtype=np.float64)
    b = np.asarray(varied_val, dtype=np.float64)
    if not (np.all(np.isfinite(a)) and np.all(np.isfinite(b))):
        return "nonfinite"
    diff = float(np.max(np.abs(a - b)))
    scale = float(np.max(np.abs(a)))
    threshold = atol + rtol * scale
    return "changed" if diff > threshold else "unchanged"


OBSERVABLE_BLOCKS = {
    "dark_energy.w0_cpl_latetime": lambda r: r["dark_energy"]["w0_cpl_latetime"],
    "dark_energy.wa_cpl_latetime": lambda r: r["dark_energy"]["wa_cpl_latetime"],
    "dark_energy.w_of_z_grid": lambda r: r["dark_energy"]["w_of_z_grid"],
    "dark_energy.H_of_z_over_H0_grid": lambda r: r["dark_energy"]["H_of_z_over_H0_grid"],
    "dark_energy.D_M_times_H0_grid": lambda r: r["dark_energy"]["D_M_times_H0_grid"],
    "dark_energy.distance_modulus_shape_grid": lambda r: r["dark_energy"]["distance_modulus_shape_grid"],
    "screening.screening_suppression_factor": lambda r: r["screening"]["screening_suppression_factor"],
    "screening.phi_center_ratio": lambda r: r["screening"]["phi_center_ratio"],
    "pta.gamma_theta": lambda r: r["pta"]["gamma_theta"],
}


def run_selftest() -> Dict[str, Any]:
    report: Dict[str, Any] = {"checks": []}
    ok_all = True

    # --- Check A: defaults reproduce current workshopcosmo w0/wa -----------
    print("[selftest] Check A: defaults reproduce current workshopcosmo w0/wa ...")
    quint_direct = wc.run_quintessence_simulation(t_max=QUINT_T_MAX, num_points=QUINT_NUM_POINTS)
    obs_direct = wc.run_observables_analysis(quint_direct)
    de = compute_dark_energy_observables(DEFAULT_PARAMS["a_pot"], DEFAULT_PARAMS["b_pot"])
    match_w0 = de["legacy_w0_fit"] == obs_direct["w0_fit"]
    match_wa = de["legacy_wa_fit"] == obs_direct["wa_fit"]
    checkA = {
        "name": "defaults_reproduce_workshopcosmo_w0_wa",
        "direct_w0_fit": obs_direct["w0_fit"],
        "direct_wa_fit": obs_direct["wa_fit"],
        "harness_legacy_w0_fit": de["legacy_w0_fit"],
        "harness_legacy_wa_fit": de["legacy_wa_fit"],
        "pass": bool(match_w0 and match_wa),
    }
    print(f"    direct w0={obs_direct['w0_fit']!r} wa={obs_direct['wa_fit']!r}")
    print(f"    harness w0={de['legacy_w0_fit']!r} wa={de['legacy_wa_fit']!r}")
    print(f"    PASS={checkA['pass']}")
    report["checks"].append(checkA)
    ok_all = ok_all and checkA["pass"]

    # --- Check B: x10 sensitivity per parameter -----------------------------
    print("[selftest] Check B: x10 parameter sensitivity per observable block ...")
    baseline = evaluate_point(DEFAULT_PARAMS)
    baseline_vals = {name: fn(baseline) for name, fn in OBSERVABLE_BLOCKS.items()}
    sensitivity: Dict[str, Dict[str, Any]] = {}
    for pname in DEFAULT_PARAMS:
        varied_params = dict(DEFAULT_PARAMS)
        varied_params[pname] = DEFAULT_PARAMS[pname] * 10.0
        result = evaluate_point(varied_params)
        changed, unchanged, nonfinite = [], [], []
        for oname, fn in OBSERVABLE_BLOCKS.items():
            verdict = _compare_observable(baseline_vals[oname], fn(result))
            {"changed": changed, "unchanged": unchanged, "nonfinite": nonfinite}[verdict].append(oname)
        sensitivity[pname] = {"changed": changed, "unchanged": unchanged, "nonfinite": nonfinite}
        print(f"    {pname} x10 -> changed: {changed}")
        print(f"    {pname} x10 -> UNCHANGED (structurally unobservable via this param): {unchanged}")
        if nonfinite:
            print(f"    {pname} x10 -> NONFINITE (solver diverged to NaN/Inf -- numerical instability, NOT evidence of insensitivity): {nonfinite}")
    report["sensitivity_x10"] = sensitivity

    # --- Check C: negative control, c4_c0_ratio/pta_suppression degeneracy -
    print("[selftest] Check C: negative control -- c4_c0_ratio*pta_suppression product degeneracy ...")
    base_pta = compute_pta_observable(DEFAULT_PARAMS["pta_suppression"], DEFAULT_PARAMS["c4_c0_ratio"])
    scaled_pta = compute_pta_observable(DEFAULT_PARAMS["pta_suppression"] / 10.0, DEFAULT_PARAMS["c4_c0_ratio"] * 10.0)
    max_diff = _scalar_diff(base_pta["gamma_theta"], scaled_pta["gamma_theta"])
    checkC = {
        "name": "c4_c0_ratio_pta_suppression_product_degeneracy",
        "base_pta_suppression": DEFAULT_PARAMS["pta_suppression"],
        "base_c4_c0_ratio": DEFAULT_PARAMS["c4_c0_ratio"],
        "scaled_pta_suppression": DEFAULT_PARAMS["pta_suppression"] / 10.0,
        "scaled_c4_c0_ratio": DEFAULT_PARAMS["c4_c0_ratio"] * 10.0,
        "max_abs_diff_gamma_theta": max_diff,
        "pass": bool(max_diff < 1e-12),
    }
    print(f"    max|Gamma(theta) diff| = {max_diff:.3e}  (threshold 1e-12)  PASS={checkC['pass']}")
    print("    -> c4_c0_ratio and pta_suppression enter Gamma(theta) ONLY as the product")
    print("       c4_c0_ratio*pta_suppression (see workshopcosmo.py:766); they are exactly")
    print("       degenerate in every observable this harness computes. Per ground rule (c)")
    print("       this combination counts as a single free parameter, not two: 6 -> 5.")
    report["checks"].append(checkC)
    ok_all = ok_all and checkC["pass"]

    report["all_pass"] = ok_all
    return report


# =============================================================================
# CLI
# =============================================================================
def main() -> int:
    parser = argparse.ArgumentParser(description="Parametrized dark-energy/screening/PTA experiment harness")
    parser.add_argument("--params", type=str, default=None, help='JSON dict, e.g. \'{"a_pot":1.0,...}\'')
    parser.add_argument("--out", type=str, default=None, help="Output JSON file path")
    parser.add_argument("--selftest", action="store_true", help="Run self-test suite instead of a single point")
    args = parser.parse_args()

    if args.selftest:
        report = run_selftest()
        if args.out:
            with open(args.out, "w", encoding="utf-8") as f:
                json.dump(report, f, indent=2)
        print(f"[selftest] ALL_PASS={report['all_pass']}")
        return 0 if report["all_pass"] else 1

    params: Dict[str, float] = {}
    if args.params:
        params = json.loads(args.params)
        for k in params:
            if k not in DEFAULT_PARAMS:
                raise SystemExit(f"Unknown parameter '{k}'. Valid keys: {sorted(DEFAULT_PARAMS)}")

    result = evaluate_point(params)
    if result["wall_time_seconds"] >= 20.0:
        print(f"WARNING: point took {result['wall_time_seconds']:.2f}s (>= 20s budget)", file=sys.stderr)

    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(result, f, indent=2)
        print(f"Wrote {args.out}")
    else:
        print(json.dumps(result, indent=2))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
