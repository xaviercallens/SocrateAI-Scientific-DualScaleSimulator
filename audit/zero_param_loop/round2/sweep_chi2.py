#!/usr/bin/env python3
"""
ROUND 2: chi2(BAO+SN) for every numerically-stable row of sweep.csv (the
5-param model, real DESI 2024 BAO 12x12 covariance + Pantheon+ SH0ES 1590
SNe, same machinery as decisive_experiment.py / round1's chi2_and_
jacobian.py), to find this round's best-fit point honestly from the
sweep grid (not a continuous re-optimization, so it cannot manufacture an
artificially low chi2 the model doesn't actually reach elsewhere).

Writes: sweep_chi2.csv (sweep.csv + chi2_bao, chi2_sn, chi2_total columns)
and sweep_chi2_report.json (best-fit row, comparison to LCDM baselines
from decisive_experiment.py).
"""
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, REPO_ROOT)

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
DESI_COV_INV = np.linalg.inv(DESI_COV)


def chi2_bao_given_base(base_vec):
    if not np.all(np.isfinite(base_vec)):
        return float("nan"), float("nan")
    denom = float(base_vec @ DESI_COV_INV @ base_vec)
    if denom <= 0 or not np.isfinite(denom):
        return float("nan"), float("nan")
    u_star = float(base_vec @ DESI_COV_INV @ DESI_DATA_VEC) / denom
    if u_star <= 0:
        return float("nan"), float("nan")
    resid = DESI_DATA_VEC - u_star * base_vec
    return float(resid @ DESI_COV_INV @ resid), float(1.0 / u_star)


def load_pantheon():
    path = os.path.join(REPO_ROOT, "data", "real", "dark_energy", "pantheon_plus_sh0es.dat")
    df = pd.read_csv(path, sep=r"\s+", engine="python")
    df = df[(df["zHD"] > 0.01) & (df["zHD"] <= 2.4) & np.isfinite(df["MU_SH0ES"]) & (df["MU_SH0ES_ERR_DIAG"] > 0)]
    return df["zHD"].to_numpy(), df["MU_SH0ES"].to_numpy(), df["MU_SH0ES_ERR_DIAG"].to_numpy()


SN_Z, SN_MU, SN_ERR = load_pantheon()
SN_W = 1.0 / SN_ERR ** 2


def chi2_sn_given_mu(mu_model):
    if not np.all(np.isfinite(mu_model)):
        return float("nan"), float("nan")
    resid0 = SN_MU - mu_model
    offset = float(np.sum(SN_W * resid0) / np.sum(SN_W))
    resid = resid0 - offset
    return float(np.sum(SN_W * resid ** 2)), offset


def row_chi2(row):
    xM = {z: row[f"DM_H0_z{z}"] for z in DESI_Z_UNIQUE}
    xH = {z: 1.0 / row[f"HzH0_z{z}"] if row[f"HzH0_z{z}"] not in (0, None) and np.isfinite(row[f"HzH0_z{z}"]) else float("nan")
          for z in DESI_Z_UNIQUE}
    base_vec = []
    for z, _, q in DESI_ROWS:
        if q == "DM_over_rs":
            base_vec.append(xM[z])
        elif q == "DH_over_rs":
            base_vec.append(xH[z])
        elif q == "DV_over_rs":
            base_vec.append((z * xM[z] ** 2 * xH[z]) ** (1.0 / 3.0) if np.isfinite(xM[z]) and np.isfinite(xH[z]) else float("nan"))
    base_vec = np.array(base_vec, dtype=float)
    chi2_bao, s_fit = chi2_bao_given_base(base_vec)

    z_grid = np.array(json.loads(row["z_grid_json"]))
    DM_H0_grid = np.array(json.loads(row["DM_H0_grid_json"]))
    DM_H0_q = np.interp(SN_Z, z_grid, DM_H0_grid)
    DL_H0_q = (1.0 + SN_Z) * DM_H0_q
    mu_model = 5.0 * np.log10(np.clip(DL_H0_q, 1e-12, None))
    chi2_sn, offset = chi2_sn_given_mu(mu_model)

    return chi2_bao, s_fit, chi2_sn, offset


def run():
    df = pd.read_csv(os.path.join(HERE, "sweep.csv"))
    stable = df[(df["error"].isna()) & (df["numerically_stable"] == True)].copy()  # noqa: E712
    chi2_bao_list, s_fit_list, chi2_sn_list, offset_list = [], [], [], []
    for _, row in stable.iterrows():
        cb, sf, cs, off = row_chi2(row)
        chi2_bao_list.append(cb); s_fit_list.append(sf); chi2_sn_list.append(cs); offset_list.append(off)
    stable["chi2_bao"] = chi2_bao_list
    stable["s_fit_bao"] = s_fit_list
    stable["chi2_sn"] = chi2_sn_list
    stable["offset_sn"] = offset_list
    stable["chi2_total"] = stable["chi2_bao"] + stable["chi2_sn"]

    finite = stable[np.isfinite(stable["chi2_total"])].copy()
    print(f"[sweep_chi2] {len(stable)} numerically_stable rows; {len(finite)} with finite chi2_total")
    best_idx = finite["chi2_total"].idxmin()
    best_row = finite.loc[best_idx]

    out_csv = os.path.join(HERE, "sweep_chi2.csv")
    stable.to_csv(out_csv, index=False)

    report = {
        "n_stable_rows": int(len(stable)),
        "n_finite_chi2_rows": int(len(finite)),
        "best_fit_row": {
            "idx": int(best_row["idx"]),
            "a_pot": float(best_row["a_pot"]), "b_pot": float(best_row["b_pot"]),
            "mu_sym": float(best_row["mu_sym"]), "lambda_sym": float(best_row["lambda_sym"]),
            "c4_pta_product": float(best_row["c4_pta_product"]),
            "chi2_bao": float(best_row["chi2_bao"]), "chi2_sn": float(best_row["chi2_sn"]),
            "chi2_total": float(best_row["chi2_total"]),
        },
        "note": "This is a GRID minimum over the 513-point Sobol sweep, not a continuous re-optimization; "
                "it is a lower bound on how good this sweep's best point looks, used only to pick a Jacobian "
                "probe point, not as a claimed global best fit. Compare to decisive_experiment.py's independent "
                "LCDM chi2_total=702.68 (Om fitted) / 702.91 (Om,H0 frozen to Lean).",
    }
    with open(os.path.join(HERE, "sweep_chi2_report.json"), "w") as f:
        json.dump(report, f, indent=2)
    print(json.dumps(report, indent=2))
    return report


if __name__ == "__main__":
    run()
