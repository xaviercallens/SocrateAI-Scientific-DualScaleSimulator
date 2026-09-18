#!/usr/bin/env python3
"""
ROUND 2 pre-registered DECISIVE EXPERIMENT (from the harness's staged
plan): flat LCDM with Omega_Lambda = 0.68885, H0 = c / hubbleRadius_m
FROZEN to the LeanMaster constants (DarkEnergyScale.lean:84,
SelfDualCutoff.lean:100; both kernel-fixed Lean `def`s, tier L as physics,
identification with this fit tier C) vs LCDM with Omega_m FITTED, on the
SAME DESI 2024 BAO (12x12 covariance) + Pantheon+ SH0ES (1701 SNe,
diagonal errors) data used throughout this loop. Also: CPL (w0, wa
fitted) on the same data, and negative control Omega_Lambda = 0.5.

IMPORTANT STRUCTURAL CAVEAT (flagged per advisor review before running):
in this harness both chi2_bao_given_base (a fitted overall scale u_star,
degenerate with r_d/H0) and chi2_sn_given_mu (a fitted additive offset,
degenerate with H0/M_B) already marginalize over H0. So "freezing H0"
changes NOTHING in this fit -- H0 never appears outside those two fitted
nuisances. The ONLY real shape parameter being frozen here is Omega_m
(equivalently Omega_Lambda). This is 1 dof, not 2. The task's
pre-registered rule (reject if delta_chi2 > 11.8, stated for "2
parameters") is applied AS GIVEN, but the correct threshold for the
actual 1-dof comparison being made (3-sigma, k=1) is 9.0. Both are
reported; this script does not silently launder the mismatch.

A genuinely independent H0-in-km/s/Mpc test (e.g. freezing the SN
absolute magnitude M_B to a fetched SH0ES calibration value) is NOT run
here: no SH0ES M_B fiducial was fetched and verified in this session, so
that specific frozen-offset variant is reported ABSENT rather than
invented, per ground rules.

Writes: decisive_experiment_report.json
"""
import json
import math
import os
import sys

import numpy as np
import pandas as pd
from scipy.optimize import minimize, minimize_scalar

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, REPO_ROOT)

# --- Lean constants (quoted verbatim, pinned v3.9.0 per round-1 verification) ---
LEAN_OMEGA_LAMBDA = 0.68885          # DarkEnergyScale.lean:84
LEAN_HUBBLE_RADIUS_M = 1.3672e26     # SelfDualCutoff.lean:100 (= c/H0)
C_LIGHT_M_S = 299792458.0
MPC_IN_M = 3.0856775814913673e22
LEAN_H0_KM_S_MPC = (C_LIGHT_M_S / LEAN_HUBBLE_RADIUS_M) * MPC_IN_M / 1000.0
OMEGA_M_FROZEN = 1.0 - LEAN_OMEGA_LAMBDA

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


def comoving_distance(zg, Ez):
    inv_E = 1.0 / Ez
    return np.concatenate(([0.0], np.cumsum(0.5 * (inv_E[1:] + inv_E[:-1]) * np.diff(zg)))), inv_E


def lcdm_Ez(zg, Om):
    return np.sqrt(Om * (1.0 + zg) ** 3 + (1.0 - Om))


def cpl_Ez(zg, Om, w0, wa):
    de_term = (1.0 + zg) ** (3.0 * (1.0 + w0 + wa)) * np.exp(-3.0 * wa * zg / (1.0 + zg))
    return np.sqrt(np.clip(Om * (1.0 + zg) ** 3 + (1.0 - Om) * de_term, 1e-300, None))


def bao_base_from_Ez(Ez_func, params):
    zmax = 2.4
    zg = np.linspace(0.0, zmax, 4000)
    Ez = Ez_func(zg, *params)
    DC, inv_E = comoving_distance(zg, Ez)
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
    return np.array(base, dtype=float)


def sn_mu_shape_from_Ez(Ez_func, params):
    zmax = float(np.max(SN_Z)) + 0.05
    zg = np.linspace(0.0, zmax, 4000)
    Ez = Ez_func(zg, *params)
    DC, _ = comoving_distance(zg, Ez)
    DC_q = np.interp(SN_Z, zg, DC)
    DL_q = (1.0 + SN_Z) * DC_q
    return 5.0 * np.log10(np.clip(DL_q, 1e-12, None))


def eval_point(Ez_func, params):
    bao = chi2_bao_given_base(bao_base_from_Ez(Ez_func, params))
    sn = chi2_sn_given_mu(sn_mu_shape_from_Ez(Ez_func, params))
    return {"bao": bao, "sn": sn, "chi2_total": bao["chi2"] + sn["chi2"]}


def run():
    # 1) LCDM Om FITTED (baseline the frozen theory is compared to)
    def neg_ll_om(Om):
        return eval_point(lcdm_Ez, (Om,))["chi2_total"]
    res_fit = minimize_scalar(neg_ll_om, bounds=(0.02, 0.98), method="bounded", options={"xatol": 1e-10})
    Om_fit = float(res_fit.x)
    lcdm_fitted = eval_point(lcdm_Ez, (Om_fit,))
    lcdm_fitted["Om_fit"] = Om_fit

    # 2) LCDM Om, H0 FROZEN to LeanMaster constants
    lcdm_frozen = eval_point(lcdm_Ez, (OMEGA_M_FROZEN,))
    lcdm_frozen["Om_frozen"] = OMEGA_M_FROZEN
    lcdm_frozen["H0_implied_km_s_Mpc"] = LEAN_H0_KM_S_MPC
    lcdm_frozen["note_H0_role"] = ("H0 does not independently enter this dimensionless-ratio chi2 (it is degenerate "
                                    "with the always-fitted BAO scale u_star and SN offset), so 'freezing H0' here "
                                    "changes nothing beyond freezing Om; see module docstring.")

    # 3) negative control: Omega_Lambda = 0.5 (Om = 0.5), must be worse than fitted
    lcdm_neg_control = eval_point(lcdm_Ez, (0.5,))
    lcdm_neg_control["Om"] = 0.5

    # 4) CPL fit (w0, wa fitted; Om fixed at frozen Lean value, matching what the harness's own
    #    dark-energy sector varies -- w(z), not Om -- so this answers "does the fetched data ask
    #    for a CPL departure from LCDM" using the SAME Om baseline as (2), not a re-fit of Om too)
    def neg_ll_cpl(x):
        w0, wa = x
        return eval_point(cpl_Ez, (OMEGA_M_FROZEN, w0, wa))["chi2_total"]
    res_cpl = minimize(neg_ll_cpl, x0=[-1.0, 0.0], method="Nelder-Mead",
                        options={"xatol": 1e-8, "fatol": 1e-8, "maxiter": 2000})
    w0_fit, wa_fit = float(res_cpl.x[0]), float(res_cpl.x[1])
    cpl_fitted = eval_point(cpl_Ez, (OMEGA_M_FROZEN, w0_fit, wa_fit))
    cpl_fitted["w0_fit"] = w0_fit
    cpl_fitted["wa_fit"] = wa_fit
    cpl_fitted["Om_held_at_frozen_lean_value"] = OMEGA_M_FROZEN

    # 5) also fit (Om, w0, wa) jointly for a fully-free CPL comparison
    def neg_ll_cpl_free_om(x):
        Om, w0, wa = x
        if not (0.02 < Om < 0.98):
            return 1e12
        return eval_point(cpl_Ez, (Om, w0, wa))["chi2_total"]
    res_cpl_free = minimize(neg_ll_cpl_free_om, x0=[Om_fit, -1.0, 0.0], method="Nelder-Mead",
                             options={"xatol": 1e-8, "fatol": 1e-8, "maxiter": 4000})
    Om_cpl, w0_cpl, wa_cpl = [float(v) for v in res_cpl_free.x]
    cpl_free_om = eval_point(cpl_Ez, (Om_cpl, w0_cpl, wa_cpl))
    cpl_free_om.update({"Om_fit": Om_cpl, "w0_fit": w0_cpl, "wa_fit": wa_cpl})

    delta_chi2_frozen_minus_fitted = lcdm_frozen["chi2_total"] - lcdm_fitted["chi2_total"]
    preregistered_threshold = 11.8   # as stated in the task, "2 parameters" 3-sigma
    correct_threshold_1dof_3sigma = 9.0  # 3.0**2, since only Om is actually frozen here (see docstring)

    reject_frozen_preregistered_rule = bool(delta_chi2_frozen_minus_fitted > preregistered_threshold)
    reject_frozen_correct_1dof_rule = bool(delta_chi2_frozen_minus_fitted > correct_threshold_1dof_3sigma)

    neg_control_ok = bool(lcdm_neg_control["chi2_total"] > lcdm_fitted["chi2_total"])

    report = {
        "lean_constants_used": {
            "omegaLambda": {"value": LEAN_OMEGA_LAMBDA, "file_line": "DarkEnergyScale.lean:84",
                             "tier_as_lean_def": "kernel-fixed", "tier_as_physics": "L (Planck 2018 literature value)"},
            "hubbleRadius_m": {"value": LEAN_HUBBLE_RADIUS_M, "file_line": "SelfDualCutoff.lean:100",
                                "tier_as_lean_def": "kernel-fixed", "tier_as_physics": "L"},
            "H0_implied_km_s_Mpc": LEAN_H0_KM_S_MPC,
            "Omega_m_frozen": OMEGA_M_FROZEN,
            "identification_with_this_fit_tier": "C (the claim that THIS toy model's Omega_m equals the Lean value is asserted, not derived)",
        },
        "lcdm_om_fitted": lcdm_fitted,
        "lcdm_frozen_lean": lcdm_frozen,
        "lcdm_negative_control_om_0.5": lcdm_neg_control,
        "cpl_fit_om_frozen_w0wa_free": cpl_fitted,
        "cpl_fit_all_free": cpl_free_om,
        "delta_chi2_frozen_minus_fitted": delta_chi2_frozen_minus_fitted,
        "preregistered_threshold_as_stated_in_task": preregistered_threshold,
        "reject_frozen_theory_per_preregistered_rule": reject_frozen_preregistered_rule,
        "actual_dof_removed_by_freezing": 1,
        "correct_1dof_3sigma_threshold": correct_threshold_1dof_3sigma,
        "reject_frozen_theory_per_correct_1dof_rule": reject_frozen_correct_1dof_rule,
        "dof_caveat": "H0 is degenerate with the always-fitted BAO scale and SN offset in this dimensionless-ratio "
                       "setup, so freezing it removes 0 real dof; only Om is actually frozen (1 dof), not 2. The "
                       "task's pre-registered 11.8 (stated for '2 parameters') is reported alongside the correct "
                       "1-dof threshold (9.0) rather than silently substituted.",
        "negative_control_omega_0.5_worse_than_fitted": neg_control_ok,
        "sh0es_frozen_offset_variant": "ABSENT: no SH0ES M_B fiducial was fetched and verified in this session; "
                                        "inventing one would violate the ground rule against invented experimental "
                                        "bounds, so this specific sub-variant (SN offset frozen to a Cepheid-"
                                        "calibrated absolute magnitude, expected to show the Hubble tension as a "
                                        "chi2 penalty) is not run.",
        "cpl_preference_in_fetched_data": {
            "delta_chi2_cpl_vs_lcdm_frozen_om": cpl_fitted["chi2_total"] - lcdm_frozen["chi2_total"],
            "delta_chi2_cpl_free_om_vs_lcdm_fitted": cpl_free_om["chi2_total"] - lcdm_fitted["chi2_total"],
            "interpretation": "negative delta favors CPL over LCDM at that same Om baseline; AIC/BIC penalty for "
                               "the 2 (or 3) extra fitted parameters is not subtracted here, only raw chi2.",
        },
        "n_bao": 12, "n_sn": int(lcdm_fitted["sn"]["n_sn"]), "n_data": 12 + int(lcdm_fitted["sn"]["n_sn"]),
    }
    with open(os.path.join(HERE, "decisive_experiment_report.json"), "w") as f:
        json.dump(report, f, indent=2)
    print(json.dumps({k: v for k, v in report.items()}, indent=2))
    return report


if __name__ == "__main__":
    run()
