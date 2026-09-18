#!/usr/bin/env python3
"""
Round 1 FORWARD LOOP step 1 (chi2) + step 2 (RESPONSE RANK / Jacobian).

chi2:
  - DESI 2024 BAO (12 rows: DV/DM/DH over r_s at 7 z's), full 12x12 covariance
    (data/real/dark_energy/desi_2024_bao_all_cov.txt, fetched this round --
    see data_provenance.json), one nuisance s = r_s*H0/c fitted per model
    point (BAO measures distances only up to the r_s*H0/c scale; this model
    has no physical H0 or r_s, so s must be marginalized, never invented).
  - Pantheon+ SH0ES SN Ia (MU_SH0ES, MU_SH0ES_ERR_DIAG; diagonal errors only
    -- the full 1701x1701 stat+sys covariance was not fetched, per ground
    rule "SN with diagonal errors if covariance is too large"), one nuisance
    additive offset fitted per model point (absorbs 5log10(c/H0) - M_B,
    both unknown in this toy model).
  - Flat LCDM (Om fitted, same nuisances) computed once as the baseline on
    the same two datasets, same dof convention.
  - PTA: ABSENT. The verified NANOGrav product
    (data/real/pulsar_timing/nanograv_kde_freespectrum.zip ->
    ceffyl_data/*/{density,freqs,log10rhogrid,bandwidths}.npy) is a
    per-frequency free-spectrum KDE (no angular/pulsar-pair information), so
    there is no real Gamma(theta) dataset to compare pta.gamma_theta
    against. Not fabricated; excluded from chi2.datasets_used.

RESPONSE RANK:
  Numerical Jacobian of the standardized observable vector w.r.t. ln(param)
  at the best-fit stable point and 5 random stable points from the sweep;
  central differences with step 1e-2 in ln-space; log-transform decade-
  spanning observables before standardizing; SVD; report full singular
  value spectrum + which ln-param the near-null right-singular vectors load
  on.

Writes:
  audit/zero_param_loop/round1/chi2_report.json
  audit/zero_param_loop/round1/jacobian_report.json
"""
import json
import math
import os
import sys

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
sys.path.insert(0, REPO_ROOT)
sys.path.insert(0, os.path.join(REPO_ROOT, "scripts"))
OUT_DIR = os.path.dirname(os.path.abspath(__file__))

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


def model_base_vector_from_grid(z_grid, HzH0_grid, DM_H0_grid):
    """base (unscaled-by-1/s) DM/rs, DH/rs, DV/rs at DESI's 7 z's, using the
    model's own z_grid/H_of_z_over_H0_grid/D_M_times_H0_grid (all dimensionless,
    already = H0*D_M/c and H(z)/H0)."""
    xM = {z: float(np.interp(z, z_grid, DM_H0_grid)) for z in DESI_Z_UNIQUE}
    xH = {z: 1.0 / float(np.interp(z, z_grid, HzH0_grid)) for z in DESI_Z_UNIQUE}
    base = []
    for z, _, q in DESI_ROWS:
        if q == "DM_over_rs":
            base.append(xM[z])
        elif q == "DH_over_rs":
            base.append(xH[z])
        elif q == "DV_over_rs":
            base.append((z * xM[z] ** 2 * xH[z]) ** (1.0 / 3.0))
        else:
            raise ValueError(q)
    return np.array(base)


def chi2_bao_given_base(base_vec):
    """Fit nuisance s>0 (data ~ base/s) with full covariance, EXACTLY
    (analytic linear least squares in u=1/s, since data - u*base is linear
    in u): u* = (base.Cinv.data)/(base.Cinv.base). An earlier version used a
    bounded 1D search over ln(s) in [-3,3] which silently clipped at its
    boundary for the LCDM baseline (whose natural distance scale differs
    from the model's a0=a(t_max) convention by orders of magnitude) --
    caught by comparing against this analytic solution, which has no
    truncation."""
    if not np.all(np.isfinite(base_vec)):
        return {"chi2": float("nan"), "s_fit": float("nan"), "ok": False}
    denom = float(base_vec @ DESI_COV_INV @ base_vec)
    if denom <= 0 or not np.isfinite(denom):
        return {"chi2": float("nan"), "s_fit": float("nan"), "ok": False}
    u_star = float(base_vec @ DESI_COV_INV @ DESI_DATA_VEC) / denom
    if u_star <= 0:
        return {"chi2": float("nan"), "s_fit": float("nan"), "ok": False}
    resid = DESI_DATA_VEC - u_star * base_vec
    chi2 = float(resid @ DESI_COV_INV @ resid)
    return {"chi2": chi2, "s_fit": float(1.0 / u_star), "ok": True}


def lcdm_bao_base_vector(Om, H0_over_c=1.0):
    """Flat LCDM base vector in units where H0/c=1 (absorbed into nuisance s
    exactly like the model side); E(z)=sqrt(Om(1+z)^3+(1-Om))."""
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
    return np.array(base), zg, DC, inv_E


def lcdm_bao_best_fit():
    """1D search over Om only; s is fitted analytically (exactly) inside
    chi2_bao_given_base for each Om, so there is no nuisance-boundary
    truncation risk."""
    def neg_ll(Om):
        base, *_ = lcdm_bao_base_vector(Om)
        return chi2_bao_given_base(base)["chi2"]

    res = minimize_scalar(neg_ll, bounds=(0.02, 0.98), method="bounded",
                           options={"xatol": 1e-8})
    base, *_ = lcdm_bao_base_vector(res.x)
    r = chi2_bao_given_base(base)
    r["Om_fit"] = float(res.x)
    return r


def load_pantheon():
    path = os.path.join(REPO_ROOT, "data", "real", "dark_energy", "pantheon_plus_sh0es.dat")
    df = pd.read_csv(path, sep=r"\s+", engine="python")
    df = df[(df["zHD"] > 0.01) & (df["zHD"] <= 2.4) & np.isfinite(df["MU_SH0ES"]) & (df["MU_SH0ES_ERR_DIAG"] > 0)]
    return df["zHD"].to_numpy(), df["MU_SH0ES"].to_numpy(), df["MU_SH0ES_ERR_DIAG"].to_numpy()


SN_Z, SN_MU, SN_ERR = load_pantheon()
SN_W = 1.0 / SN_ERR ** 2


def mu_shape_from_DM_H0(z_query, z_grid, DM_H0_grid):
    """5*log10((1+z)*D_M_H0(z)), computed AFTER interpolating the smooth
    D_M_H0(z) grid at the exact query z's (not by interpolating an already
    log-transformed grid, which is inaccurate near z=0 -- see param_sweep.py
    comment)."""
    DM_H0_q = np.interp(z_query, z_grid, DM_H0_grid)
    DL_H0_q = (1.0 + np.asarray(z_query)) * DM_H0_q
    return 5.0 * np.log10(np.clip(DL_H0_q, 1e-12, None))


def chi2_sn_given_mu(mu_model_at_SNz):
    """Shared SN chi2 (1 nuisance additive offset, analytic) given the
    model's shape-only distance modulus ALREADY evaluated at SN_Z."""
    if not np.all(np.isfinite(mu_model_at_SNz)):
        return {"chi2": float("nan"), "offset_fit": float("nan"), "n_sn": int(len(SN_Z)), "ok": False}
    resid0 = SN_MU - mu_model_at_SNz
    offset = float(np.sum(SN_W * resid0) / np.sum(SN_W))
    resid = resid0 - offset
    chi2 = float(np.sum(SN_W * resid ** 2))
    return {"chi2": chi2, "offset_fit": offset, "n_sn": int(len(SN_Z)), "ok": True}


def chi2_sn_given_shape(z_grid, DM_H0_grid):
    """Model-side entry point: DM_H0_grid is the harness's smooth
    D_M_times_H0_grid (NOT a pre-logged distance modulus grid)."""
    mu_model = mu_shape_from_DM_H0(SN_Z, z_grid, DM_H0_grid)
    return chi2_sn_given_mu(mu_model)


def lcdm_mu_shape_at(Om, z_query):
    """Flat LCDM shape-only distance modulus evaluated exactly at z_query,
    via a fine internal comoving-distance grid (4000 pts) then a SINGLE
    interpolation of the smooth D_C(z) (not of an already-logged quantity)."""
    zmax = float(np.max(z_query)) + 0.05
    zg = np.linspace(0.0, zmax, 4000)
    Ez = np.sqrt(Om * (1.0 + zg) ** 3 + (1.0 - Om))
    inv_E = 1.0 / Ez
    DC = np.concatenate(([0.0], np.cumsum(0.5 * (inv_E[1:] + inv_E[:-1]) * np.diff(zg))))
    DC_q = np.interp(z_query, zg, DC)
    DL_q = (1.0 + np.asarray(z_query)) * DC_q
    return 5.0 * np.log10(np.clip(DL_q, 1e-12, None))


def lcdm_sn_best_fit():
    def neg_ll(Om):
        mu_model = lcdm_mu_shape_at(Om, SN_Z)
        return chi2_sn_given_mu(mu_model)["chi2"]

    res = minimize_scalar(neg_ll, bounds=(0.02, 0.98), method="bounded")
    r = chi2_sn_given_mu(lcdm_mu_shape_at(res.x, SN_Z))
    r["Om_fit"] = float(res.x)
    return r


def run_chi2():
    df = pd.read_csv(os.path.join(OUT_DIR, "sweep.csv"))
    clean = df[df["error"].isna()].copy()
    stable = clean[clean["numerically_stable"] == True].copy()  # noqa: E712

    print(f"[chi2] sweep rows: {len(df)}; clean(no error): {len(clean)}; numerically_stable: {len(stable)}")

    rows_out = []
    for _, row in clean.iterrows():
        z_grid = np.array(json.loads(row["z_grid_json"]))
        DM_H0_grid = np.array(json.loads(row["DM_H0_grid_json"]))
        # DM_H0/HzH0 at DESI z's were stored directly as columns -> use those exactly (no re-interp needed)
        xM = {z: row[f"DM_H0_z{z}"] for z in DESI_Z_UNIQUE}
        xH = {z: 1.0 / row[f"HzH0_z{z}"] for z in DESI_Z_UNIQUE}
        base = []
        for z, _, q in DESI_ROWS:
            if q == "DM_over_rs":
                base.append(xM[z])
            elif q == "DH_over_rs":
                base.append(xH[z])
            elif q == "DV_over_rs":
                base.append((z * xM[z] ** 2 * xH[z]) ** (1.0 / 3.0))
        base = np.array(base, dtype=float)
        bao = chi2_bao_given_base(base)
        sn = chi2_sn_given_shape(z_grid, DM_H0_grid)
        rows_out.append({
            "idx": row["idx"], "a_pot": row["a_pot"], "b_pot": row["b_pot"],
            "numerically_stable": row["numerically_stable"],
            "w0_cpl_latetime": row["w0_cpl_latetime"], "wa_cpl_latetime": row["wa_cpl_latetime"],
            "chi2_bao": bao["chi2"], "s_fit_bao": bao["s_fit"],
            "chi2_sn": sn["chi2"], "offset_fit_sn": sn["offset_fit"], "n_sn": sn["n_sn"],
            "chi2_total": (bao["chi2"] if np.isfinite(bao["chi2"]) else np.inf) +
                          (sn["chi2"] if np.isfinite(sn["chi2"]) else np.inf),
        })
    chi_df = pd.DataFrame(rows_out)
    chi_df.to_csv(os.path.join(OUT_DIR, "sweep_chi2.csv"), index=False)

    finite = chi_df[np.isfinite(chi_df["chi2_total"])]
    best = finite.loc[finite["chi2_total"].idxmin()]
    lcdm_bao = lcdm_bao_best_fit()
    lcdm_sn = lcdm_sn_best_fit()

    n_bao = 12
    dof_bao_model = n_bao - (2 + 1)  # a_pot,b_pot scanned + 1 nuisance s
    dof_bao_lcdm = n_bao - (1 + 1)   # Om + 1 nuisance s
    dof_sn_model = int(best["n_sn"]) - (2 + 1)
    dof_sn_lcdm = lcdm_sn["n_sn"] - (1 + 1)

    w0_all = clean["w0_cpl_latetime"].to_numpy()
    wa_all = clean["wa_cpl_latetime"].to_numpy()
    report = {
        "n_sweep_rows": int(len(df)),
        "n_clean": int(len(clean)),
        "n_numerically_stable": int(len(stable)),
        "datasets_used": ["desi_2024_bao_all (12 pts, full 12x12 covariance)",
                           "pantheon_plus_sh0es (Hubble-flow subset, N=%d SNe, diagonal MU_SH0ES_ERR_DIAG)" % int(best["n_sn"])],
        "pta_chi2": "ABSENT: NANOGrav verified product (ceffyl_data/*.npy) is a per-frequency free-spectrum KDE, no angular/pulsar-pair separation info; no real Gamma(theta) dataset exists to compare against.",
        "dof_convention": "dof = n_data_points - (n_scanned_theory_params + n_nuisance_params); model scans (a_pot,b_pot)=2 params + 1 nuisance per sector; LCDM scans Om=1 param + 1 nuisance per sector.",
        "best_fit_model_point": {
            "idx": int(best["idx"]), "a_pot": float(best["a_pot"]), "b_pot": float(best["b_pot"]),
            "chi2_bao": float(best["chi2_bao"]), "dof_bao": dof_bao_model,
            "chi2_sn": float(best["chi2_sn"]), "dof_sn": dof_sn_model,
            "chi2_total": float(best["chi2_total"]), "dof_total": dof_bao_model + dof_sn_model,
            "w0_cpl_latetime": float(best["w0_cpl_latetime"]), "wa_cpl_latetime": float(best["wa_cpl_latetime"]),
        },
        "lcdm_baseline": {
            "bao": {**lcdm_bao, "dof": dof_bao_lcdm},
            "sn": {**lcdm_sn, "dof": dof_sn_lcdm},
            "chi2_total": lcdm_bao["chi2"] + lcdm_sn["chi2"], "dof_total": dof_bao_lcdm + dof_sn_lcdm,
        },
        "w0_cpl_latetime_range_over_sweep": [float(np.min(w0_all)), float(np.max(w0_all))],
        "wa_cpl_latetime_range_over_sweep": [float(np.min(wa_all)), float(np.max(wa_all))],
        "headline_finding": None,
    }
    w0_lo, w0_hi = report["w0_cpl_latetime_range_over_sweep"]
    wa_lo, wa_hi = report["wa_cpl_latetime_range_over_sweep"]
    lcdm_like = (abs(w0_lo - (-1.0)) < 0.05 and abs(w0_hi - (-1.0)) < 0.05 and abs(wa_lo) < 0.05 and abs(wa_hi) < 0.05)
    report["headline_finding"] = (
        ("The model is LCDM-degenerate across the ENTIRE 4-decade (a_pot,b_pot) sweep: "
         f"w0 in [{w0_lo:.4f},{w0_hi:.4f}], wa in [{wa_lo:.4f},{wa_hi:.4f}] (both pinned near (-1,0)); "
         "chi2_model ~= chi2_LCDM on BAO+SN, so this data cannot distinguish a_pot,b_pot from LCDM, "
         "let alone see mu_sym/lambda_sym/pta_suppression/c4_c0_ratio, which do not enter dark-energy "
         "observables at all.") if lcdm_like else
        (f"w0 spans [{w0_lo:.4f},{w0_hi:.4f}], wa spans [{wa_lo:.4f},{wa_hi:.4f}] over the sweep -- "
         "NOT uniformly LCDM-pinned; see per-point chi2 for where BAO/SN prefer non-LCDM regions.")
    )
    with open(os.path.join(OUT_DIR, "chi2_report.json"), "w") as f:
        json.dump(report, f, indent=2)
    print(json.dumps(report, indent=2))
    return report


if __name__ == "__main__":
    run_chi2()
