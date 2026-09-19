#!/usr/bin/env python3
"""Skeptic (statistics lens) check of E1 (DESI DR2 + Pantheon+ full cov) and E3 (M0 vs fitted LCDM).
Imports E1's own chi2 machinery (audit/reverse_zero/E1-desi-dr2/e1_desi_dr2.py) and checks:
 (1) pre-registered thresholds recomputed with scipy (2-sided normal equivalents);
 (2) CPL minimum robustness: 6 Nelder-Mead starts + coarse (w0,wa) grid with Om_m profiled;
 (3) sigma conversion of Delta chi2 (2 dof -> two-sided normal);
 (4) SN-cut consistency E1 (1590, calibrators kept) vs E3 (1580, IS_CALIBRATOR==0): Delta chi2(CPL vs LCDM)
     and Delta chi2(frozen Om_L=0.68885 vs fitted LCDM, 1 dof) on the E3 cut -- cross-validates E3's 0.746
     from E1's independent code path;
 (5) diagonal-vs-full SN covariance on DR2 AND on DR1 (R3's setup: DR1 + 1590 SNe diagonal) to test E1's
     claim that R3's 2.0 sigma (Delta chi2 = 6.09) was inflated by neglecting SN covariance.
No RNG. Command (from worktree root):
  /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python audit/reverse_zero/skeptic_statistics/check_e1_e3_chi2.py
"""
import json, os, sys
import numpy as np, pandas as pd
from scipy.optimize import minimize, minimize_scalar
from scipy.stats import chi2 as C, norm
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "E1-desi-dr2"))
import e1_desi_dr2 as e1  # noqa: E402

NM = dict(method="Nelder-Mead", options={"xatol": 1e-8, "fatol": 1e-9, "maxiter": 20000, "maxfev": 20000})


def fits(chi):
    f = lambda p: chi(e1.E_cpl(*p))["chi2"]
    r = minimize(f, [0.30, -0.9, -0.2], **NM)
    l = minimize_scalar(lambda o: chi(e1.E_lcdm(o))["chi2"], bounds=(0.02, 0.98), method="bounded", options={"xatol": 1e-10})
    frozen = chi(e1.E_lcdm(1 - e1.OMEGA_LAMBDA))["chi2"]
    return {"cpl": [float(x) for x in r.x], "cpl_chi2": float(r.fun), "lcdm_Om": float(l.x), "lcdm_chi2": float(l.fun),
            "dchi2_cpl_vs_lcdm_2dof": float(l.fun - r.fun), "sigma_2dof": float(norm.isf(C.sf(l.fun - r.fun, 2) / 2)),
            "dchi2_frozen_vs_fitted_lcdm_1dof": float(frozen - l.fun)}


out = {"thresholds": {"2dof_3sig": float(C.isf(2 * norm.sf(3), 2)), "2dof_5sig": float(C.isf(2 * norm.sf(5), 2)),
                      "1dof_3sig": float(C.isf(2 * norm.sf(3), 1)), "1dof_5sig": float(C.isf(2 * norm.sf(5), 1))}}
D2, D1 = e1.D2, e1.D1
sn = pd.read_csv(os.path.join(D1, "pantheon_plus_sh0es.dat"), sep=r"\s+")
cov = np.loadtxt(os.path.join(D2, "Pantheon+SH0ES_STAT+SYS.cov"), skiprows=1).reshape(len(sn), len(sn))
mE1 = ((sn.zHD > 0.01) & (sn.zHD <= 2.4) & np.isfinite(sn.MU_SH0ES) & (sn.MU_SH0ES_ERR_DIAG > 0)).to_numpy()
mE3 = mE1 & (sn.IS_CALIBRATOR == 0).to_numpy()
out["n_sn_E1_cut"] = int(mE1.sum()); out["n_sn_E3_cut"] = int(mE3.sum()); out["n_calibrators_in_E1_cut"] = int((mE1 & ~mE3).sum())
bao = {"DR2": e1.load_bao(os.path.join(D2, "desi_dr2", "desi_gaussian_bao_ALL_GCcomb_mean.txt"), os.path.join(D2, "desi_dr2", "desi_gaussian_bao_ALL_GCcomb_cov.txt")),
       "DR1": e1.load_bao(os.path.join(D1, "desi_2024_bao_all.txt"), os.path.join(D1, "desi_2024_bao_all_cov.txt"))}
z, mu = sn.zHD.to_numpy(), sn.MU_SH0ES.to_numpy()
for tag, (bz, bv, bq, BCi) in bao.items():
    for cut_name, m in (("E1cut_1590", mE1), ("E3cut_1580", mE3)):
        full = e1.make_chi2(bz, bv, bq, BCi, z[m], mu[m], np.linalg.inv(cov[np.ix_(m, m)]))
        diag = e1.make_chi2(bz, bv, bq, BCi, z[m], mu[m], np.diag(1 / sn.MU_SH0ES_ERR_DIAG.to_numpy()[m] ** 2))
        out[f"{tag}_{cut_name}_fullcov"] = fits(full)
        out[f"{tag}_{cut_name}_diagcov"] = fits(diag)
# CPL multistart + grid on the E1 headline configuration
bz, bv, bq, BCi = bao["DR2"]
chi = e1.make_chi2(bz, bv, bq, BCi, z[mE1], mu[mE1], np.linalg.inv(cov[np.ix_(mE1, mE1)]))
f = lambda p: chi(e1.E_cpl(*p))["chi2"]
starts = [[0.30, -1, 0], [0.28, -0.7, -1], [0.33, -1.2, 0.5], [0.25, -0.5, -2], [0.35, -0.9, -0.5], [0.31, -0.8, 0.3]]
out["DR2_E1cut_cpl_multistart"] = [[*map(float, minimize(f, s, **NM).x), float(minimize(f, s, **NM).fun)] for s in starts]
best = (1e99, None)
for w0 in np.linspace(-1.3, -0.5, 17):
    for wa in np.linspace(-2.5, 1.0, 15):
        r = minimize_scalar(lambda o: f([o, w0, wa]), bounds=(0.15, 0.5), method="bounded")
        if r.fun < best[0]:
            best = (float(r.fun), [float(r.x), float(w0), float(wa)])
out["DR2_E1cut_cpl_grid_best"] = {"chi2": best[0], "Om_w0_wa": best[1], "note": "coarse grid; must be >= multistart minimum"}
out["committed_E1_DR2_dchi2"] = json.load(open(os.path.join(HERE, "..", "E1-desi-dr2", "e1_desi_dr2_report.json")))["DR2"]["delta_chi2_CPL_vs_Lambda_2dof"]
out["committed_E3_dchi2"] = json.load(open(os.path.join(HERE, "..", "E3-M0", "m0_result.json")))["delta_chi2_M0_minus_LCDM_fitted"]
json.dump(out, open(os.path.join(HERE, "check_e1_e3_chi2.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
