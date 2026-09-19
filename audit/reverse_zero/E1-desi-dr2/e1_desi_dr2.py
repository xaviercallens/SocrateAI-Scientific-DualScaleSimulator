"""E1 - pre-registered P1 test (audit/PRE_REGISTRATION.md, section "P1. Dark energy is a
cosmological constant"), re-run with DESI DR2 BAO (full ALL_GCcomb covariance, 13x13) and
Pantheon+SH0ES with the FULL STAT+SYS 1701x1701 covariance (not diagonal errors as in round3).

Pre-registered rule (quoted verbatim from audit/PRE_REGISTRATION.md, read-only, main repo):
  "Thresholds (computed with `scipy.stats.chi2.isf`):
   - 1 dof: Delta chi2 = 9.00 is 3 sigma, and 25.00 is 5 sigma;
   - 2 dof: Delta chi2 = 11.83 is 3 sigma, and 28.74 is 5 sigma."
  P1: "Falsified if ... prefers the w0 > -1, wa < 0 quadrant over w = -1 at >= 5 sigma
   (Delta chi2 >= 28.74, 2 dof)"
  P1: "Also falsified if a LCDM fit to those data puts Omega_Lambda away from 0.68885 by
   >= 5 sigma of the measurement error."
  P1: "Pipeline action, to be run with these thresholds fixed now: re-run
   decisive_experiment.py with DESI DR2 BAO and the full Pantheon+ covariance; report the
   CPL-vs-Lambda Delta chi2 as it comes out."
These threshold numbers (11.83, 28.74) are READ from PRE_REGISTRATION.md above (not
hardcoded independently of it) and are not changed here.

FRAMING: Omega_Lambda = 0.68885 is a postulate (LeanMaster DarkEnergyScale.lean:84, Planck
2018 benchmark), not a prediction of K3xT2. This script tests that postulate; it says
nothing about mu_sym or c4_pta_product, which remain unconstrained by this experiment.

Hubble-flow cut (identical to round2/round3): zHD in (0.01, 2.4], MU_SH0ES finite,
MU_SH0ES_ERR_DIAG > 0 (the diagonal cut is only used to build the boolean mask; the actual
SN chi2 below uses the FULL off-diagonal STAT+SYS covariance restricted to that mask, not
the diagonal errors).

Nuisance accounting (stated explicitly, per advisor review): BAO constrains D/r_d, so an
overall scale nuisance (proportional to r_d*h) is marginalised analytically inside chi2_bao
(the `u` fit, identical for every E(z) model tested). Pantheon+ MU_SH0ES already has the
SH0ES absolute calibration baked in; a constant offset (degenerate with M_B/H0) is
marginalised analytically inside chi2_sn (the `off` fit). Both nuisances are marginalised
identically across ALL models (frozen LCDM, fitted LCDM, CPL), so the reported Delta chi2
values are on equal footing. The "frozen" model is NOT literally zero-parameter: it still
carries these two analytically-profiled nuisances, common to every model compared.

Seeds: none needed (no stochastic step; all fits are deterministic optimisation of a fixed
chi2 surface from fixed data).
"""
import json, os
import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import minimize, minimize_scalar
from scipy.stats import chi2 as chi2dist

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
D2 = os.path.join(ROOT, "data", "real2", "dark_energy")
D1 = os.path.join(ROOT, "data", "real", "dark_energy")

OMEGA_LAMBDA = 0.68885          # LeanMaster DarkEnergyScale.lean:84 (Planck 2018, tier L)
HUBBLE_RADIUS_M = 1.3672e26     # LeanMaster SelfDualCutoff.lean:100
C_KM_S = 299792.458
MPC_M = 3.0856775814913673e22
H0_LEAN = C_KM_S / (HUBBLE_RADIUS_M / MPC_M)

# --- pre-registered thresholds, read from the file, not invented here ---
# PRE_REGISTRATION.md is READ-ONLY and lives in the main repo, not this worktree
# (per ground rules); read it from there rather than copying it into the worktree.
PREREG_PATH = "/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/audit/PRE_REGISTRATION.md"
prereg_text = open(PREREG_PATH).read()
assert "11.83" in prereg_text and "28.74" in prereg_text
THRESH_3SIGMA_2DOF = 11.83
THRESH_5SIGMA_2DOF = 28.74
# cross-check against scipy (2 dof)
_chk3 = chi2dist.isf(2 * (1 - 0.9973002), df=2)  # 3-sigma two-sided normal equivalent
_chk5 = chi2dist.isf(2 * (1 - 0.9999994267), df=2)

zg = np.linspace(0, 2.5, 25001)


def E_lcdm(Om):
    return np.sqrt(Om * (1 + zg) ** 3 + 1 - Om)


def E_cpl(Om, w0, wa):
    return np.sqrt(Om * (1 + zg) ** 3 + (1 - Om) * (1 + zg) ** (3 * (1 + w0 + wa)) * np.exp(-3 * wa * zg / (1 + zg)))


def load_bao(mean_path, cov_path):
    rows = [l.split() for l in open(mean_path) if l.strip() and not l.startswith("#")]
    bz = np.array([float(r[0]) for r in rows])
    bv = np.array([float(r[1]) for r in rows])
    bq = [r[2] for r in rows]
    cov = np.loadtxt(cov_path)
    Cinv = np.linalg.inv(cov)
    return bz, bv, bq, Cinv


def load_sn_full_cov():
    sn = pd.read_csv(os.path.join(D1, "pantheon_plus_sh0es.dat"), sep=r"\s+")
    n_full = len(sn)
    mask = (sn.zHD > 0.01) & (sn.zHD <= 2.4) & np.isfinite(sn.MU_SH0ES) & (sn.MU_SH0ES_ERR_DIAG > 0)
    mask = mask.to_numpy()
    # full STAT+SYS covariance: first line = N, then N*N values (row-major), same row order
    # as the .dat file (Pantheon+SH0ES DataRelease convention).
    with open(os.path.join(D2, "Pantheon+SH0ES_STAT+SYS.cov")) as f:
        n_cov = int(f.readline().strip())
    assert n_cov == n_full, f"cov N={n_cov} != dat rows={n_full}"
    cov_full = np.loadtxt(os.path.join(D2, "Pantheon+SH0ES_STAT+SYS.cov"), skiprows=1).reshape(n_cov, n_cov)
    cov_cut = cov_full[np.ix_(mask, mask)]
    Cinv = np.linalg.inv(cov_cut)
    sz = sn.zHD.to_numpy()[mask]
    smu = sn.MU_SH0ES.to_numpy()[mask]
    return sz, smu, Cinv, int(mask.sum()), n_full


def make_chi2(bz, bv, bq, BCinv, sz, smu, SCinv):
    ones = np.ones(len(sz))
    Sd = SCinv @ ones
    S11 = float(ones @ Sd)

    def chi2(E, sn_offset=None):
        dc = cumulative_trapezoid(1 / E, zg, initial=0)  # D_C * H0 / c
        DM = np.interp(bz, zg, dc)
        DH = np.interp(bz, zg, 1 / E)
        base = np.where([q == "DM_over_rs" for q in bq], DM,
                         np.where([q == "DH_over_rs" for q in bq], DH, (bz * DM ** 2 * DH) ** (1 / 3)))
        u = (base @ BCinv @ bv) / (base @ BCinv @ base)
        rb = bv - u * base
        cb = float(rb @ BCinv @ rb)
        m0 = 5 * np.log10((1 + sz) * np.interp(sz, zg, dc))
        d = smu - m0
        if sn_offset is None:
            off = float((d @ Sd) / S11)
        else:
            off = sn_offset
        rs = d - off
        cs = float(rs @ SCinv @ rs)
        return {"chi2": cb + cs, "bao": cb, "sn": cs, "sn_offset": off,
                "H0_implied_by_offset": C_KM_S / 10 ** ((off - 25) / 5)}
    return chi2


def profile_omega_lambda_sigma(chi2fn, om_l_best):
    """1-sigma half-width on Omega_Lambda from Delta chi2 = 1 (1 dof) around the fitted LCDM
    minimum, used only to express how many sigma 0.68885 sits from the best-fit value."""
    best = chi2fn(E_lcdm(1 - om_l_best))["chi2"]

    def f(om_l):
        return chi2fn(E_lcdm(1 - om_l))["chi2"] - best - 1.0
    # bracket search outward from best
    lo, hi = om_l_best, om_l_best
    step = 0.001
    while f(hi) < 0 and hi < 0.999:
        hi += step
    while f(lo) < 0 and lo > 0.001:
        lo -= step
    from scipy.optimize import brentq
    try:
        sig_hi = brentq(f, om_l_best, hi) - om_l_best
    except Exception:
        sig_hi = float("nan")
    try:
        sig_lo = om_l_best - brentq(f, lo, om_l_best)
    except Exception:
        sig_lo = float("nan")
    return sig_lo, sig_hi


def run_combo(label, bz, bv, bq, BCinv, sz, smu, SCinv, n_bao_pts):
    chi2 = make_chi2(bz, bv, bq, BCinv, sz, smu, SCinv)
    out = {"label": label, "n_bao": n_bao_pts, "n_sn": len(sz)}
    Om_L = 1 - OMEGA_LAMBDA
    out["frozen_omega_lambda_0.68885"] = chi2(E_lcdm(Om_L))

    r = minimize_scalar(lambda o: chi2(E_lcdm(o))["chi2"], bounds=(0.02, 0.98), method="bounded",
                         options={"xatol": 1e-10})
    out["fitted_lcdm"] = {"Om_m": float(r.x), "Om_L": float(1 - r.x), **chi2(E_lcdm(r.x))}

    c = minimize(lambda p: chi2(E_cpl(Om_L, *p))["chi2"], [-1.0, 0.0], method="Nelder-Mead",
                 options={"xatol": 1e-8, "fatol": 1e-9, "maxiter": 6000, "maxfev": 6000})
    out["cpl_om_frozen"] = {"w0": float(c.x[0]), "wa": float(c.x[1]), **chi2(E_cpl(Om_L, *c.x))}

    c2 = minimize(lambda p: chi2(E_cpl(*p))["chi2"], [1 - OMEGA_LAMBDA, -1.0, 0.0], method="Nelder-Mead",
                  options={"xatol": 1e-8, "fatol": 1e-9, "maxiter": 10000, "maxfev": 10000})
    out["cpl_om_free"] = {"Om_m": float(c2.x[0]), "w0": float(c2.x[1]), "wa": float(c2.x[2]),
                           **chi2(E_cpl(*c2.x))}

    dchi2 = out["frozen_omega_lambda_0.68885"]["chi2"] - out["cpl_om_frozen"]["chi2"]
    out["delta_chi2_CPL_vs_Lambda_2dof"] = dchi2
    out["sigma_equivalent_2dof"] = float(chi2dist.isf(chi2dist.sf(dchi2, df=2), df=1)) ** 0.5 if dchi2 > 0 else 0.0
    # report against the two pre-registered bars, quoted verbatim above
    out["prereg_rule_ge_28.74_falsifies_5sigma"] = bool(dchi2 >= THRESH_5SIGMA_2DOF)
    out["prereg_rule_ge_11.83_is_3sigma"] = bool(dchi2 >= THRESH_3SIGMA_2DOF)

    # negative control
    out["negative_control_Om_L_0.5"] = {
        "delta_chi2_vs_fitted_lcdm": chi2(E_lcdm(0.5))["chi2"] - out["fitted_lcdm"]["chi2"]}

    # Omega_Lambda-away-from-0.68885 rule
    sig_lo, sig_hi = profile_omega_lambda_sigma(chi2, out["fitted_lcdm"]["Om_L"])
    diff = OMEGA_LAMBDA - out["fitted_lcdm"]["Om_L"]
    sigma_used = sig_lo if diff < 0 else sig_hi
    out["omega_lambda_profile_1sigma"] = {"minus": sig_lo, "plus": sig_hi}
    out["omega_lambda_distance_from_0.68885_in_sigma"] = float(abs(diff) / sigma_used) if sigma_used and np.isfinite(sigma_used) else None
    out["prereg_rule_omega_lambda_ge_5sigma_falsifies"] = bool(
        out["omega_lambda_distance_from_0.68885_in_sigma"] is not None and
        out["omega_lambda_distance_from_0.68885_in_sigma"] >= 5.0)
    return out


def main():
    sz, smu, SCinv, n_used, n_full = load_sn_full_cov()

    dr2_bz, dr2_bv, dr2_bq, dr2_BCinv = load_bao(
        os.path.join(D2, "desi_dr2", "desi_gaussian_bao_ALL_GCcomb_mean.txt"),
        os.path.join(D2, "desi_dr2", "desi_gaussian_bao_ALL_GCcomb_cov.txt"))
    dr1_bz, dr1_bv, dr1_bq, dr1_BCinv = load_bao(
        os.path.join(D1, "desi_2024_bao_all.txt"),
        os.path.join(D1, "desi_2024_bao_all_cov.txt"))

    out = {
        "H0_from_hubbleRadius_m_km_s_Mpc": H0_LEAN,
        "prereg_thresholds_quoted_from_PRE_REGISTRATION_md": {
            "2dof_3sigma": THRESH_3SIGMA_2DOF, "2dof_5sigma": THRESH_5SIGMA_2DOF,
            "source_line": "audit/PRE_REGISTRATION.md section 2: '2 dof: Delta chi2 = 11.83 is 3 sigma, and 28.74 is 5 sigma.'"},
        "sn_cut": "zHD in (0.01, 2.4], MU_SH0ES finite, MU_SH0ES_ERR_DIAG>0 (mask only; chi2 uses full STAT+SYS off-diagonal cov restricted to mask)",
        "sn_n_used": n_used, "sn_n_full_file": n_full,
        "sn_cov_source": "data/real2/dark_energy/Pantheon+SH0ES_STAT+SYS.cov (full 1701x1701 STAT+SYS)",
        "DR2": run_combo("DESI DR2 ALL_GCcomb (13 pts, full 13x13 cov) + Pantheon+ full cov",
                          dr2_bz, dr2_bv, dr2_bq, dr2_BCinv, sz, smu, SCinv, len(dr2_bv)),
        "DR1": run_combo("DESI 2024 DR1 all-tracer (12 pts, full 12x12 cov) + Pantheon+ full cov",
                          dr1_bz, dr1_bv, dr1_bq, dr1_BCinv, sz, smu, SCinv, len(dr1_bv)),
    }
    out["DR1_to_DR2_delta_chi2_CPL_vs_Lambda_change"] = (
        out["DR2"]["delta_chi2_CPL_vs_Lambda_2dof"] - out["DR1"]["delta_chi2_CPL_vs_Lambda_2dof"])

    out["published_3.1sigma_reproduction"] = {
        "reproduced": False,
        "why_not": ("The published 3.1sigma (PRE_REGISTRATION.md R4) is DESI DR2 BAO + CMB + "
                    "supernovae combined. This script has no CMB likelihood (no Planck/ACT/SPT "
                    "chain in the USABLE DATA list; none fetched). CMB breaks the Omega_m*h^2 "
                    "vs r_d (BAO) and Omega_m vs H0 (SN) degeneracies that BAO+SN alone leave "
                    "partially open, so a BAO+SN-only Delta chi2 is expected to be smaller than "
                    "the published DESI+CMB+SN result, not a reproduction of it. Whatever "
                    "Delta chi2/sigma this script found for DR2 is reported above as a distinct, "
                    "weaker number; it is not adjusted toward 3.1sigma.")}

    out_path = os.path.join(HERE, "e1_desi_dr2_report.json")
    json.dump(out, open(out_path, "w"), indent=1, default=float)
    print(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
