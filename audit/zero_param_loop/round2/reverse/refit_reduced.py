#!/usr/bin/env python3
"""
REVERSE LOOP (round 2), step 3: EXPERIMENT. Refit the reduced model to the
same real data (DESI 2024 BAO 12x12-covariance + Pantheon+ SH0ES SNe) and
compare chi2/AIC with the full model and with LCDM.

This is an INDEPENDENT re-implementation (own DESI/Pantheon loading, own
chi2_bao/chi2_sn, own LCDM Ez), not an import of
../decisive_experiment.py, so its two overlapping numbers
(chi2_reduced == frozen-LCDM chi2_total, chi2_lcdm == Om-fitted chi2_total)
serve as a genuine cross-check, not a restated copy. chi2_full is
explicitly REUSED (not recomputed) from ../sweep_chi2.csv, since
reproducing the full quintessence-ODE optimizer is out of scope here; this
is stated, not hidden.

mu_sym and c4_pta_product contribute ZERO to chi2_reduced: no dataset in
data/real/ constrains either (fifth_force_screening and PTA angular
Gamma(theta) are both ABSENT, data/real/MANIFEST.json). chi2_reduced is
therefore the chi2 of a model with 0 parameters fit to this data, not 2.
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

HERE = os.path.dirname(os.path.abspath(__file__))
ROUND2_DIR = os.path.abspath(os.path.join(HERE, ".."))
REPO_ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))

# --- Lean constants (same pins as ../decisive_experiment.py and reduced_model.py) ---
LEAN_OMEGA_LAMBDA = 0.68885           # DarkEnergyScale.lean:84
LEAN_HUBBLE_RADIUS_M = 1.3672e26      # SelfDualCutoff.lean:100
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
assert DESI_COV.shape == (12, 12)
DESI_COV_INV = np.linalg.inv(DESI_COV)


def load_pantheon():
    path = os.path.join(REPO_ROOT, "data", "real", "dark_energy", "pantheon_plus_sh0es.dat")
    df = pd.read_csv(path, sep=r"\s+", engine="python")
    df = df[(df["zHD"] > 0.01) & (df["zHD"] <= 2.4) & np.isfinite(df["MU_SH0ES"]) & (df["MU_SH0ES_ERR_DIAG"] > 0)]
    return df["zHD"].to_numpy(), df["MU_SH0ES"].to_numpy(), df["MU_SH0ES_ERR_DIAG"].to_numpy()


SN_Z, SN_MU, SN_ERR = load_pantheon()
SN_W = 1.0 / SN_ERR ** 2


def chi2_bao_given_base(base_vec):
    denom = float(base_vec @ DESI_COV_INV @ base_vec)
    u_star = float(base_vec @ DESI_COV_INV @ DESI_DATA_VEC) / denom
    resid = DESI_DATA_VEC - u_star * base_vec
    return {"chi2": float(resid @ DESI_COV_INV @ resid), "s_fit": float(1.0 / u_star)}


def chi2_sn_given_mu(mu_model_at_SNz):
    resid0 = SN_MU - mu_model_at_SNz
    offset = float(np.sum(SN_W * resid0) / np.sum(SN_W))
    resid = resid0 - offset
    return {"chi2": float(np.sum(SN_W * resid ** 2)), "offset_fit": offset, "n_sn": int(len(SN_Z))}


def lcdm_Ez(zg, Om):
    return np.sqrt(Om * (1.0 + zg) ** 3 + (1.0 - Om))


def comoving_distance(zg, Ez):
    inv_E = 1.0 / Ez
    return np.concatenate(([0.0], np.cumsum(0.5 * (inv_E[1:] + inv_E[:-1]) * np.diff(zg)))), inv_E


def bao_base(Om):
    zg = np.linspace(0.0, 2.4, 4000)
    Ez = lcdm_Ez(zg, Om)
    DC, inv_E = comoving_distance(zg, Ez)
    xM = {z: float(np.interp(z, zg, DC)) for z in DESI_Z_UNIQUE}
    xH = {z: float(np.interp(z, zg, inv_E)) for z in DESI_Z_UNIQUE}
    base = []
    for z, _, q in DESI_ROWS:
        if q == "DM_over_rs":
            base.append(xM[z])
        elif q == "DH_over_rs":
            base.append(xH[z])
        else:
            base.append((z * xM[z] ** 2 * xH[z]) ** (1.0 / 3.0))
    return np.array(base)


def sn_mu_shape(Om):
    zmax = float(np.max(SN_Z)) + 0.05
    zg = np.linspace(0.0, zmax, 4000)
    Ez = lcdm_Ez(zg, Om)
    DC, _ = comoving_distance(zg, Ez)
    DC_q = np.interp(SN_Z, zg, DC)
    DL_q = (1.0 + SN_Z) * DC_q
    return 5.0 * np.log10(np.clip(DL_q, 1e-12, None))


def eval_lcdm(Om):
    bao = chi2_bao_given_base(bao_base(Om))
    sn = chi2_sn_given_mu(sn_mu_shape(Om))
    return {"bao": bao, "sn": sn, "chi2_total": bao["chi2"] + sn["chi2"]}


def run():
    # chi2_reduced: frozen Om (0 free params fit to this data)
    reduced = eval_lcdm(OMEGA_M_FROZEN)
    reduced["Om_frozen"] = OMEGA_M_FROZEN
    reduced["H0_implied_km_s_Mpc"] = LEAN_H0_KM_S_MPC
    chi2_reduced = reduced["chi2_total"]

    # chi2_lcdm: Om fitted (1 free param)
    res = minimize_scalar(lambda Om: eval_lcdm(Om)["chi2_total"], bounds=(0.02, 0.98),
                           method="bounded", options={"xatol": 1e-10})
    Om_fit = float(res.x)
    lcdm_fitted = eval_lcdm(Om_fit)
    lcdm_fitted["Om_fit"] = Om_fit
    chi2_lcdm = lcdm_fitted["chi2_total"]

    # chi2_full: REUSED from the forward round's Sobol sweep of (a_pot,b_pot)
    sweep_chi2_path = os.path.join(ROUND2_DIR, "sweep_chi2.csv")
    full_source = "ABSENT"
    chi2_full = None
    if os.path.exists(sweep_chi2_path):
        df_full = pd.read_csv(sweep_chi2_path)
        col = "chi2_total" if "chi2_total" in df_full.columns else None
        if col:
            df_full_valid = df_full[np.isfinite(df_full[col])]
            chi2_full = float(df_full_valid[col].min())
            full_source = f"{sweep_chi2_path} (min of {len(df_full_valid)} finite points, reused not recomputed)"

    # cross-check against forward round's decisive_experiment_report.json
    cross_check = {}
    dec_path = os.path.join(ROUND2_DIR, "decisive_experiment_report.json")
    if os.path.exists(dec_path):
        with open(dec_path) as f:
            dec = json.load(f)
        forward_frozen = dec["lcdm_frozen_lean"]["chi2_total"]
        forward_fitted = dec["lcdm_om_fitted"]["chi2_total"]
        cross_check = {
            "forward_lcdm_frozen_chi2_total": forward_frozen,
            "this_script_chi2_reduced": chi2_reduced,
            "diff_frozen": chi2_reduced - forward_frozen,
            "forward_lcdm_fitted_chi2_total": forward_fitted,
            "this_script_chi2_lcdm": chi2_lcdm,
            "diff_fitted": chi2_lcdm - forward_fitted,
        }

    # AIC bookkeeping: k=0 reduced, k=1 fitted-Om LCDM, k=2 full quintessence
    # (a_pot,b_pot only -- mu_sym/lambda_sym/c4_pta_product are not fit to
    # this data by ANY of the three models, so excluded from k throughout)
    k_reduced, k_lcdm, k_full = 0, 1, 2
    aic_reduced = chi2_reduced + 2 * k_reduced
    aic_lcdm = chi2_lcdm + 2 * k_lcdm
    aic_full = (chi2_full + 2 * k_full) if chi2_full is not None else None

    delta_aic_vs_lcdm = aic_reduced - aic_lcdm
    delta_aic_vs_full = (aic_reduced - aic_full) if aic_full is not None else None

    delta_chi2_reduced_minus_full = (chi2_reduced - chi2_full) if chi2_full is not None else None
    n_removed_chi2_relevant = 2  # a_pot, b_pot only
    n_removed_conservative = 3   # + lambda_sym, even though it doesn't touch chi2 either
    fit_degraded_2removed = bool(delta_chi2_reduced_minus_full is not None and
                                  delta_chi2_reduced_minus_full > 2 * n_removed_chi2_relevant)
    fit_degraded_3removed = bool(delta_chi2_reduced_minus_full is not None and
                                  delta_chi2_reduced_minus_full > 2 * n_removed_conservative)

    report = {
        "chi2_reduced": chi2_reduced,
        "chi2_reduced_detail": reduced,
        "chi2_reduced_data_facing_free_params": 0,
        "chi2_reduced_caveat": "mu_sym and c4_pta_product are free in the model but contribute 0 to this "
                                "number: no fifth-force or angular-PTA dataset exists in data/real/ to fit "
                                "them against. This is the chi2 of a 0-parameter model, not a 2-parameter fit.",
        "chi2_lcdm": chi2_lcdm,
        "chi2_lcdm_detail": lcdm_fitted,
        "chi2_full": chi2_full,
        "chi2_full_source": full_source,
        "cross_check_vs_forward_decisive_experiment": cross_check,
        "aic": {"k_reduced": k_reduced, "k_lcdm": k_lcdm, "k_full": k_full,
                "aic_reduced": aic_reduced, "aic_lcdm": aic_lcdm, "aic_full": aic_full,
                "sign_convention": "delta_aic_vs_X = AIC_reduced - AIC_X; negative favors the reduced model"},
        "delta_aic_vs_lcdm": delta_aic_vs_lcdm,
        "delta_aic_vs_full": delta_aic_vs_full,
        "delta_chi2_reduced_minus_full": delta_chi2_reduced_minus_full,
        "fit_degraded_threshold_2removed_params": 2 * n_removed_chi2_relevant,
        "fit_degraded_2removed": fit_degraded_2removed,
        "fit_degraded_threshold_3removed_params_conservative": 2 * n_removed_conservative,
        "fit_degraded_3removed_conservative": fit_degraded_3removed,
        "fit_degraded": fit_degraded_2removed,
        "n_bao": 12, "n_sn": int(reduced["sn"]["n_sn"]),
    }
    with open(os.path.join(HERE, "refit_reduced_report.json"), "w") as f:
        json.dump(report, f, indent=2)
    print(json.dumps({k: v for k, v in report.items() if k not in ("chi2_reduced_detail", "chi2_lcdm_detail")}, indent=2))
    return report


if __name__ == "__main__":
    run()
