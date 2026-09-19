#!/usr/bin/env python3
"""
REVERSE LOOP step 3 (EXPERIMENT): refit the reduced (5-param) model to the
SAME real data (DESI 2024 BAO + Pantheon+ SH0ES) and compare chi2/AIC/BIC
with the full model (round-1's own sweep.csv/chi2_report.json) and with
flat LCDM (recomputed HERE, independently, not imported from round 1's
chi2_and_jacobian.py, as a consistency cross-check).

KEY STRUCTURAL FACT, verified below rather than assumed: the tier-B
reduction (pta_suppression,c4_c0_ratio -> c4_pta_product) touches ONLY the
PTA sector (workshopcosmo.py's run_nanograv_hexadecapole_simulation path).
It does not touch compute_dark_energy_observables(a_pot,b_pot) at all, and
mu_sym/lambda_sym only touch the screening sector. Neither BAO nor SN chi2
uses the PTA or screening blocks. Therefore the reduced model's chi2(a_pot,
b_pot) is LITERALLY THE SAME FUNCTION as the full model's -- not
approximately, not "expected to be similar." This script proves that with
a byte-level check at forward's own best-fit point, rather than re-deriving
a new best fit from a different, non-comparable grid (round-1's grid vs a
fresh grid would let an optimizer beat a grid-min and manufacture a
negative-looking delta-chi2 for a nested model, which is impossible and
would just be a grid-vs-optimum artifact -- see advisor guidance recorded
in this round's transcript).

Writes: reduced_chi2_report.json
"""
import json
import math
import os
import sys

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
ROUND1_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, REPO_ROOT)
sys.path.insert(0, os.path.join(REPO_ROOT, "scripts"))
sys.path.insert(0, HERE)

import param_loop_sim as pls  # noqa: E402
# CRITICAL: reproduce round-1's param_sweep.py _worker_init grid extension
# EXACTLY. Forward's sweep.csv/chi2_report.json were computed with
# pls.Z_GRID extended to linspace(0,2.45,80) and pls.Z_POINTS set to DESI's
# 7 z's (covering DESI z=2.330 and Pantheon+ max z~2.26 WITHOUT boundary
# clipping). Without this override, the module-level default Z_GRID
# (linspace(0,2.3,40)) clips np.interp at z=2.33 (DESI's highest z) and
# near z~2.26-2.3 (many Pantheon+ SNe), corrupting chi2_bao/chi2_sn -- this
# was caught empirically: an earlier run of this script (without the
# override) gave chi2_total=709.0044 instead of the byte-identical
# 703.8819 the equivalence check requires, entirely from this clipping.
pls.Z_GRID = np.linspace(0.0, 2.45, 80)

DESI_ROWS = [
    (0.295, 7.92512927, "DV_over_rs"),
    (0.510, 13.62003080, "DM_over_rs"),
    (0.510, 20.98334647, "DH_over_rs"),
    (0.706, 16.84645313, "DM_over_rs"),
    (0.706, 20.07872919, "DH_over_rs"),
    (0.930, 21.70841761, "DM_over_rs"),
    (0.930, 17.87612922, "DH_over_rs"),
    (1.317, 27.78720817, "DM_over_rs"),
    (1.317, 13.82372285, "DH_over_rs"),
    (1.491, 26.07217182, "DV_over_rs"),
    (2.330, 39.70838281, "DM_over_rs"),
    (2.330, 8.52256583, "DH_over_rs"),
]
DESI_Z_UNIQUE = sorted(set(z for z, _, _ in DESI_ROWS))
pls.Z_POINTS = list(DESI_Z_UNIQUE)  # match round-1 worker_init exactly

import reduced_model as rm  # noqa: E402

DESI_DATA_VEC = np.array([v for _, v, _ in DESI_ROWS])
DESI_COV = np.loadtxt(os.path.join(REPO_ROOT, "data", "real", "dark_energy", "desi_2024_bao_all_cov.txt"))
assert DESI_COV.shape == (12, 12), DESI_COV.shape
DESI_COV_INV = np.linalg.inv(DESI_COV)


def chi2_bao_given_base(base_vec):
    if not np.all(np.isfinite(base_vec)):
        return {"chi2": float("nan"), "s_fit": float("nan"), "ok": False}
    denom = float(base_vec @ DESI_COV_INV @ base_vec)
    if denom <= 0 or not np.isfinite(denom):
        return {"chi2": float("nan"), "s_fit": float("nan"), "ok": False}
    u_star = float(base_vec @ DESI_COV_INV @ DESI_DATA_VEC) / denom
    if u_star <= 0:
        return {"chi2": float("nan"), "s_fit": float("nan"), "ok": False}
    resid = DESI_DATA_VEC - u_star * base_vec
    return {"chi2": float(resid @ DESI_COV_INV @ resid), "s_fit": float(1.0 / u_star), "ok": True}


def load_pantheon():
    path = os.path.join(REPO_ROOT, "data", "real", "dark_energy", "pantheon_plus_sh0es.dat")
    df = pd.read_csv(path, sep=r"\s+", engine="python")
    df = df[(df["zHD"] > 0.01) & (df["zHD"] <= 2.4) & np.isfinite(df["MU_SH0ES"]) & (df["MU_SH0ES_ERR_DIAG"] > 0)]
    return df["zHD"].to_numpy(), df["MU_SH0ES"].to_numpy(), df["MU_SH0ES_ERR_DIAG"].to_numpy()


SN_Z, SN_MU, SN_ERR = load_pantheon()
SN_W = 1.0 / SN_ERR ** 2


def chi2_sn_given_mu(mu_model_at_SNz):
    if not np.all(np.isfinite(mu_model_at_SNz)):
        return {"chi2": float("nan"), "offset_fit": float("nan"), "n_sn": int(len(SN_Z)), "ok": False}
    resid0 = SN_MU - mu_model_at_SNz
    offset = float(np.sum(SN_W * resid0) / np.sum(SN_W))
    resid = resid0 - offset
    return {"chi2": float(np.sum(SN_W * resid ** 2)), "offset_fit": offset, "n_sn": int(len(SN_Z)), "ok": True}


def model_base_bao_from_de(de):
    # Use the harness's OWN "_points" arrays (evaluated by np.interp directly
    # against the raw high-resolution ODE trajectory, since pls.Z_POINTS was
    # overridden to DESI_Z_UNIQUE above) rather than re-interpolating the
    # already-downsampled 80-point z_grid. An earlier version of this
    # function used the grid re-interpolation and got chi2_bao=16.494
    # instead of forward's reported 16.672 -- a double-interpolation
    # discrepancy of ~0.18 in chi2_total, caught by the byte-level
    # equivalence check below (see param_loop_sim.py:116-138 for the
    # exact _points computation this now matches).
    DM_H0_points = de["D_M_times_H0_points"]
    Hz_H0_points = de["H_of_z_over_H0_points"]
    xM = {z: float(DM_H0_points[j]) for j, z in enumerate(DESI_Z_UNIQUE)}
    xH = {z: 1.0 / float(Hz_H0_points[j]) for j, z in enumerate(DESI_Z_UNIQUE)}
    base = []
    for z, _, q in DESI_ROWS:
        if q == "DM_over_rs":
            base.append(xM[z])
        elif q == "DH_over_rs":
            base.append(xH[z])
        elif q == "DV_over_rs":
            base.append((z * xM[z] ** 2 * xH[z]) ** (1.0 / 3.0))
    return np.array(base, dtype=float)


def mu_shape_from_de(de):
    z_grid = np.array(de["z_grid"])
    DM_H0_grid = np.array(de["D_M_times_H0_grid"])
    DM_H0_q = np.interp(SN_Z, z_grid, DM_H0_grid)
    DL_H0_q = (1.0 + SN_Z) * DM_H0_q
    return 5.0 * np.log10(np.clip(DL_H0_q, 1e-12, None))


# ---------------------------------------------------------------------------
# 1) Byte-level equivalence proof: reduced model at forward's best-fit
#    (a_pot,b_pot), any c4_pta_product, must give IDENTICAL dark-energy
#    observables and hence identical chi2 to forward's reported best fit.
# ---------------------------------------------------------------------------
FORWARD_BEST = {
    "idx": 282,
    "a_pot": 3.4714217520576933,
    "b_pot": 0.0004276423734589,
    "chi2_bao_reported": 16.672122654385387,
    "chi2_sn_reported": 687.209783100074,
    "chi2_total_reported": 703.8819057544594,
    "w0_cpl_latetime_reported": -0.694151106969398,
    "wa_cpl_latetime_reported": 0.954542950540572,
}


def verify_equivalence_and_compute_reduced_chi2_at_bestfit():
    r = rm.evaluate_reduced_point({"a_pot": FORWARD_BEST["a_pot"], "b_pot": FORWARD_BEST["b_pot"]})
    de = r["dark_energy"]
    w0_diff = abs(de["w0_cpl_latetime"] - FORWARD_BEST["w0_cpl_latetime_reported"])
    wa_diff = abs(de["wa_cpl_latetime"] - FORWARD_BEST["wa_cpl_latetime_reported"])
    base = model_base_bao_from_de(de)
    bao = chi2_bao_given_base(base)
    mu_model = mu_shape_from_de(de)
    sn = chi2_sn_given_mu(mu_model)
    chi2_total = bao["chi2"] + sn["chi2"]
    return {
        "reduced_model_params_used": r["reduced_params"],
        "w0_cpl_latetime": de["w0_cpl_latetime"],
        "wa_cpl_latetime": de["wa_cpl_latetime"],
        "abs_diff_w0_vs_forward_reported": w0_diff,
        "abs_diff_wa_vs_forward_reported": wa_diff,
        "chi2_bao": bao["chi2"], "s_fit_bao": bao["s_fit"],
        "chi2_sn": sn["chi2"], "offset_fit_sn": sn["offset_fit"], "n_sn": sn["n_sn"],
        "chi2_total": chi2_total,
        "abs_diff_chi2_total_vs_forward_reported": abs(chi2_total - FORWARD_BEST["chi2_total_reported"]),
        "byte_level_equivalence_confirmed": bool(w0_diff < 1e-9 and wa_diff < 1e-9 and
                                                   abs(chi2_total - FORWARD_BEST["chi2_total_reported"]) < 1e-6),
    }


# ---------------------------------------------------------------------------
# 2) Independent LCDM recompute (own implementation, same public data files,
#    NOT calling round-1's chi2_and_jacobian.py) -- consistency check
#    against round-1's reported chi2_total=697.3402981903151.
# ---------------------------------------------------------------------------
def lcdm_bao_base(Om):
    zmax = 2.4
    zg = np.linspace(0.0, zmax, 4000)
    Ez = np.sqrt(Om * (1.0 + zg) ** 3 + (1.0 - Om))
    inv_E = 1.0 / Ez
    DC = np.concatenate(([0.0], np.cumsum(0.5 * (inv_E[1:] + inv_E[:-1]) * np.diff(zg))))
    xM = {z: float(np.interp(z, zg, DC)) for z in DESI_Z_UNIQUE}
    xH = {z: float(np.interp(z, zg, inv_E)) for z in DESI_Z_UNIQUE}
    base = []
    for z, _, q in DESI_ROWS:
        if q == "DM_over_rs":
            base.append(xM[z])
        elif q == "DH_over_rs":
            base.append(xH[z])
        elif q == "DV_over_rs":
            base.append((z * xM[z] ** 2 * xH[z]) ** (1.0 / 3.0))
    return np.array(base, dtype=float), zg, DC


def lcdm_bao_best_fit_independent():
    def neg_ll(Om):
        base, *_ = lcdm_bao_base(Om)
        return chi2_bao_given_base(base)["chi2"]
    res = minimize_scalar(neg_ll, bounds=(0.02, 0.98), method="bounded", options={"xatol": 1e-8})
    base, *_ = lcdm_bao_base(res.x)
    r = chi2_bao_given_base(base)
    r["Om_fit"] = float(res.x)
    return r


def lcdm_mu_shape(Om):
    zmax = float(np.max(SN_Z)) + 0.05
    zg = np.linspace(0.0, zmax, 4000)
    Ez = np.sqrt(Om * (1.0 + zg) ** 3 + (1.0 - Om))
    inv_E = 1.0 / Ez
    DC = np.concatenate(([0.0], np.cumsum(0.5 * (inv_E[1:] + inv_E[:-1]) * np.diff(zg))))
    DC_q = np.interp(SN_Z, zg, DC)
    DL_q = (1.0 + SN_Z) * DC_q
    return 5.0 * np.log10(np.clip(DL_q, 1e-12, None))


def lcdm_sn_best_fit_independent():
    def neg_ll(Om):
        return chi2_sn_given_mu(lcdm_mu_shape(Om))["chi2"]
    res = minimize_scalar(neg_ll, bounds=(0.02, 0.98), method="bounded")
    r = chi2_sn_given_mu(lcdm_mu_shape(res.x))
    r["Om_fit"] = float(res.x)
    return r


def run():
    eq = verify_equivalence_and_compute_reduced_chi2_at_bestfit()
    lcdm_bao = lcdm_bao_best_fit_independent()
    lcdm_sn = lcdm_sn_best_fit_independent()
    chi2_lcdm_total = lcdm_bao["chi2"] + lcdm_sn["chi2"]

    forward_reported_lcdm_total = 697.3402981903151  # chi2_report.json lcdm_baseline.chi2_total
    lcdm_consistency_abs_diff = abs(chi2_lcdm_total - forward_reported_lcdm_total)

    chi2_full = FORWARD_BEST["chi2_total_reported"]
    chi2_reduced = eq["chi2_total"]
    delta_chi2_reduced_minus_full = chi2_reduced - chi2_full
    n_removed_params = 1  # (pta_suppression, c4_c0_ratio) -> c4_pta_product
    fit_degraded = bool(delta_chi2_reduced_minus_full > 2.0 * n_removed_params)

    n_bao, n_sn = 12, eq["n_sn"]
    N_data = n_bao + n_sn
    ln_N = math.log(N_data)

    # k_fit convention (PRIMARY, per advisor guidance): count only params
    # actually driving the fitted observables + nuisances, identically
    # defined for all three models so differences carry the data content.
    k_scan_full, k_scan_reduced, k_scan_lcdm = 2, 2, 1   # (a_pot,b_pot) vs (Om)
    k_nuis = 2  # s (BAO) + offset (SN), same for all three
    k_fit_full = k_scan_full + k_nuis
    k_fit_reduced = k_scan_reduced + k_nuis
    k_fit_lcdm = k_scan_lcdm + k_nuis

    aic = lambda chi2, k: chi2 + 2 * k
    bic = lambda chi2, k: chi2 + k * ln_N

    aic_full = aic(chi2_full, k_fit_full)
    aic_reduced = aic(chi2_reduced, k_fit_reduced)
    aic_lcdm = aic(chi2_lcdm_total, k_fit_lcdm)
    bic_full = bic(chi2_full, k_fit_full)
    bic_reduced = bic(chi2_reduced, k_fit_reduced)
    bic_lcdm = bic(chi2_lcdm_total, k_fit_lcdm)

    # Nominal-k bookkeeping (SECONDARY, explicitly labeled contentless):
    # counts the THEORY's total nominal free-parameter list (6 full / 5
    # reduced / 1 LCDM), not the fit-relevant subset.
    k_nominal_full, k_nominal_reduced, k_nominal_lcdm = 6, 5, 1
    aic_nominal_full = aic(chi2_full, k_nominal_full)
    aic_nominal_reduced = aic(chi2_reduced, k_nominal_reduced)
    aic_nominal_lcdm = aic(chi2_lcdm_total, k_nominal_lcdm)

    report = {
        "note": "PTA chi2 remains ABSENT this round for the same reason as round 1 (no real angular Gamma(theta) dataset); only BAO+SN enter chi2 here, same as round 1.",
        "equivalence_check": eq,
        "chi2_full_forward_reported": chi2_full,
        "chi2_reduced": chi2_reduced,
        "chi2_lcdm_independent_recompute": {
            "bao": lcdm_bao, "sn": lcdm_sn, "chi2_total": chi2_lcdm_total,
            "consistency_vs_forward_reported_697.3402981903151": {
                "forward_reported": forward_reported_lcdm_total,
                "this_round_independent": chi2_lcdm_total,
                "abs_diff": lcdm_consistency_abs_diff,
                "consistent": bool(lcdm_consistency_abs_diff < 0.01),
            },
        },
        "delta_chi2_reduced_minus_full": delta_chi2_reduced_minus_full,
        "n_removed_params_this_step": n_removed_params,
        "fit_degraded_rule": "delta_chi2(reduced-full) > 2 per removed parameter",
        "fit_degraded": fit_degraded,
        "why_delta_chi2_is_exactly_zero": "The removed direction (pta_suppression vs c4_c0_ratio split) lives entirely in the PTA sector, and pta_chi2 is ABSENT (no real angular-correlation dataset this round, same as round 1) -- so removing it cannot change a chi2 that never depended on it. This delta_chi2=0 has ZERO discriminating power about whether the reduction is a good idea; it only shows the reduction did not break the (unrelated) dark-energy fit.",
        "aic_bic_k_fit_convention_PRIMARY": {
            "definition": "k_fit = n_scanned_theory_params_relevant_to_chi2 + n_nuisance (2: s_BAO, offset_SN), identical nuisance convention for all three models so only real data content survives in the deltas.",
            "k_fit_full": k_fit_full, "k_fit_reduced": k_fit_reduced, "k_fit_lcdm": k_fit_lcdm,
            "aic_full": aic_full, "aic_reduced": aic_reduced, "aic_lcdm": aic_lcdm,
            "bic_full": bic_full, "bic_reduced": bic_reduced, "bic_lcdm": bic_lcdm,
            "delta_aic_reduced_minus_full": aic_reduced - aic_full,
            "delta_aic_reduced_minus_lcdm": aic_reduced - aic_lcdm,
            "delta_bic_reduced_minus_full": bic_reduced - bic_full,
            "delta_bic_reduced_minus_lcdm": bic_reduced - bic_lcdm,
            "interpretation": "delta_aic_reduced_minus_full ~ 0 (expected: chi2 identical, k_fit identical -- the reduction is invisible to this data, not validated by it). delta_aic_reduced_minus_lcdm > 0 means LCDM is STILL preferred over the reduced model by ~the same margin as round 1 found for the full model (this reduction does not change that verdict either way)."
        },
        "aic_nominal_k_bookkeeping_SECONDARY_NOT_DATA_CONTENT": {
            "definition": "k = total nominal theory free-parameter count (6/5/1), NOT fit-relevant. Decrementing this by construction gives ANY reduction a -2 'AIC improvement' with zero data involvement -- reported only as bookkeeping, never as evidence the reduction is supported by data.",
            "k_nominal_full": k_nominal_full, "k_nominal_reduced": k_nominal_reduced, "k_nominal_lcdm": k_nominal_lcdm,
            "aic_nominal_full": aic_nominal_full, "aic_nominal_reduced": aic_nominal_reduced, "aic_nominal_lcdm": aic_nominal_lcdm,
            "delta_aic_nominal_reduced_minus_full": aic_nominal_reduced - aic_nominal_full,
        },
        "dof_total": (n_bao - k_scan_full - 1) + (n_sn - k_scan_full - 1),
        "n_data_points_used": N_data,
    }
    with open(os.path.join(HERE, "reduced_chi2_report.json"), "w") as f:
        json.dump(report, f, indent=2)
    print(json.dumps({k: v for k, v in report.items() if k not in ("equivalence_check",)}, indent=2))
    print("equivalence_check.byte_level_equivalence_confirmed =", eq["byte_level_equivalence_confirmed"])
    return report


if __name__ == "__main__":
    run()
