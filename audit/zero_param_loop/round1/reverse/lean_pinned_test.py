#!/usr/bin/env python3
"""
REVERSE LOOP secondary test (tier-C, NOT counted as a reduction per ground
rules -- forward round-1 explicitly flagged this NOT ADOPTED): nested-model
test of the candidate identification a_pot,b_pot -> LeanMaster
dark_energy_w0:=-1, dark_energy_wa:=0 (DualScaleValidation/Observables.lean:
97-98, verified present at this exact line this round via grep, same as
round 1).

Step 1: is (w0_cpl_latetime, wa_cpl_latetime) = (-1, 0) even REACHABLE by
the model's (a_pot,b_pot) image? Multi-start local search minimizing
(w0+1)^2 + wa^2 over log10(a_pot) in [-2,2], log10(b_pot) in [-4,0] (the
SAME bounds round 1's sweep used). If the achieved minimum residual is not
close to zero, the pair is unreachable within the searched domain and the
identification is falsified outright, cleanly, without needing a chi2
comparison at all.

Step 2 (only if reachable): chi2 at the pinned point vs the free best-fit
chi2 (703.8819, forward-reported / reduced_chi2.py-confirmed), with 2 fewer
fitted parameters (w0,wa exactly fixed rather than 2 free scan directions).

Writes: lean_pinned_report.json
"""
import json
import os
import sys

import numpy as np
from scipy.optimize import minimize

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, REPO_ROOT)
sys.path.insert(0, os.path.join(REPO_ROOT, "scripts"))
sys.path.insert(0, HERE)

import param_loop_sim as pls  # noqa: E402
pls.Z_GRID = np.linspace(0.0, 2.45, 80)
DESI_Z_UNIQUE = [0.295, 0.510, 0.706, 0.930, 1.317, 1.491, 2.330]
pls.Z_POINTS = list(DESI_Z_UNIQUE)

import reduced_model as rm  # noqa: E402
from reduced_chi2 import (  # noqa: E402
    model_base_bao_from_de, mu_shape_from_de, chi2_bao_given_base, chi2_sn_given_mu,
)

LEAN_CONST_VERIFIED = {
    "name": "dark_energy_w0, dark_energy_wa",
    "file_line": "DualScaleValidation/Observables.lean:97-98",
    "value": "-1, 0",
    "grep_check_command": "grep -n 'dark_energy_w0\\|dark_energy_wa' "
                           "/home/callensxavier_gmail_com/SocrateAI-Scientific-Agora-LeanMaster/DualScaleValidation/Observables.lean",
}


def w0_wa_at(a_pot, b_pot):
    de = pls.compute_dark_energy_observables(a_pot, b_pot)
    if not de["success"]:
        return None, None
    return de["w0_cpl_latetime"], de["wa_cpl_latetime"]


def residual(log10_ab):
    a_pot = 10.0 ** log10_ab[0]
    b_pot = 10.0 ** log10_ab[1]
    w0, wa = w0_wa_at(a_pot, b_pot)
    if w0 is None or not (np.isfinite(w0) and np.isfinite(wa)):
        return 1e6
    return (w0 - (-1.0)) ** 2 + (wa - 0.0) ** 2


def run():
    starts = [(-2, -4), (-2, 0), (0, -2), (2, -4), (2, 0), (0, 0), (-1, -1), (1, -3)]
    best = None
    trials = []
    for s in starts:
        res = minimize(residual, x0=np.array(s, dtype=float), method="Nelder-Mead",
                        options={"xatol": 1e-6, "fatol": 1e-12, "maxiter": 2000})
        a_pot, b_pot = 10.0 ** res.x[0], 10.0 ** res.x[1]
        # clip to searched bounds for reporting (Nelder-Mead is unconstrained)
        in_bounds = bool(-2.0 <= res.x[0] <= 2.0 and -4.0 <= res.x[1] <= 0.0)
        w0, wa = w0_wa_at(a_pot, b_pot)
        trial = {"start": s, "a_pot": float(a_pot), "b_pot": float(b_pot),
                 "in_searched_bounds": in_bounds,
                 "w0_achieved": w0, "wa_achieved": wa, "residual": float(res.fun)}
        trials.append(trial)
        if best is None or (trial["in_searched_bounds"] and trial["residual"] < best["residual"]) or \
           (not best["in_searched_bounds"] and trial["in_searched_bounds"]):
            best = trial

    dist_to_target = float(np.sqrt(best["residual"])) if np.isfinite(best["residual"]) else float("nan")
    reachable = bool(best["in_searched_bounds"] and dist_to_target < 0.02)

    report = {
        "lean_constant": LEAN_CONST_VERIFIED,
        "target": {"w0": -1.0, "wa": 0.0},
        "search_bounds_log10": {"a_pot": [-2.0, 2.0], "b_pot": [-4.0, 0.0]},
        "multi_start_trials": trials,
        "best_trial": best,
        "distance_to_target_in_w0_wa_plane": dist_to_target,
        "reachable_within_searched_bounds": reachable,
    }

    if reachable:
        de = pls.compute_dark_energy_observables(best["a_pot"], best["b_pot"])
        base = model_base_bao_from_de(de)
        bao = chi2_bao_given_base(base)
        mu_model = mu_shape_from_de(de)
        sn = chi2_sn_given_mu(mu_model)
        chi2_pinned = bao["chi2"] + sn["chi2"]
        chi2_free_bestfit = 703.8819057544594
        delta_chi2 = chi2_pinned - chi2_free_bestfit
        report["pinned_point_chi2"] = {
            "a_pot": best["a_pot"], "b_pot": best["b_pot"],
            "w0_achieved": best["w0_achieved"], "wa_achieved": best["wa_achieved"],
            "chi2_bao": bao["chi2"], "chi2_sn": sn["chi2"], "chi2_total": chi2_pinned,
            "chi2_free_bestfit_for_comparison": chi2_free_bestfit,
            "delta_chi2_pinned_minus_free": delta_chi2,
            "n_fewer_fitted_params": 2,
            "verdict": ("Lean-pinned (w0=-1,wa=0) point costs Delta-chi2=%.3f for 2 fewer fitted params "
                        "(threshold for '>2 per param' = 4)" % delta_chi2),
            "supported_by_data": bool(delta_chi2 <= 4.0),
        }
    else:
        report["conclusion"] = (
            f"(w0,wa)=(-1,0) is UNREACHABLE within the searched (a_pot,b_pot) bounds "
            f"(closest achieved distance {dist_to_target:.4f} at a_pot={best['a_pot']:.4g}, "
            f"b_pot={best['b_pot']:.4g}, giving w0={best['w0_achieved']:.4f}, wa={best['wa_achieved']:.4f}). "
            "This falsifies the identification outright within this toy model's bounded parameter "
            "domain -- the Lean-quoted (w0,wa)=(-1,0) point is not in the model's observable image "
            "for reasonable a_pot,b_pot, so it cannot be adopted as a derived replacement value this round."
        )

    with open(os.path.join(HERE, "lean_pinned_report.json"), "w") as f:
        json.dump(report, f, indent=2)
    print(json.dumps(report, indent=2))
    return report


if __name__ == "__main__":
    run()
