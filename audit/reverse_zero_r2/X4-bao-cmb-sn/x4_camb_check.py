"""X4 CAMB check. (1) builds data/camb_rd_table.npz (CAMB r_drag on an (omega_b, omega_m) grid) and validates
the spline against direct CAMB at 40 seeded random points (seed 9100000). (2) Checks the published prior
table against a recomputation in CAMB at the parameters implied by the table: solve (Omega_m, h) so the
paper's own formulas reproduce (R, l_A) at omega_b = 0.02235976, then run CAMB at that point.
Usage: cd audit/reverse_zero_r2/X4-bao-cmb-sn && <venv python> x4_camb_check.py
Writes x4_camb_check.json."""
import json, pathlib, numpy as np, camb
from scipy.optimize import fsolve
import x4_core as X
HERE = pathlib.Path(__file__).resolve().parent
OUT = {}
omnu = 0.06 / 93.14


def camb_rd(wb, wm):
    h = 0.75
    p = camb.CAMBparams()
    p.set_cosmology(H0=100 * h, ombh2=wb, omch2=wm - wb - omnu, mnu=0.06, nnu=3.046, num_massive_neutrinos=1,
                    TCMB=X.TCMB, omk=0)
    return float(camb.get_background(p).get_derived_params()["rdrag"])


tab = HERE / "data" / "camb_rd_table.npz"
wbg = np.linspace(0.0185, 0.0275, 10)
wmg = np.exp(np.linspace(np.log(0.04), np.log(0.42), 60))
if not tab.exists():
    rd = np.array([[camb_rd(a, b) for b in wmg] for a in wbg])
    np.savez(tab, wb=wbg, wm=wmg, rd=rd)
X._RD.clear()
rng = np.random.default_rng(9100000)
errs = []
for _ in range(40):
    wb = rng.uniform(0.0195, 0.0265); wm = float(np.exp(rng.uniform(np.log(0.05), np.log(0.40))))
    errs.append(X.rd_camb_interp(wb, wm) / camb_rd(wb, wm) - 1)
OUT["rd_spline_vs_camb"] = {"seed": 9100000, "n": 40, "max_abs_rel_err": float(np.max(np.abs(errs))),
                            "grid": {"wb": [float(wbg[0]), float(wbg[-1]), len(wbg)], "wm": [float(wmg[0]), float(wmg[-1]), len(wmg)]}}

# --- table vs CAMB
wb = float(X.CMB_D[2])
f = lambda x: np.array([X.cmb_paper(x[0], x[1], wb)[k] for k in ("R", "lA")]) - X.CMB_D[:2]
sol = fsolve(f, [0.316, 0.673], xtol=1e-12)
Om, h = map(float, sol)
pp = X.cmb_paper(Om, h, wb)
cc = X.camb_exact(Om, h, wb)
sigR = 0.0046; sigA = 0.0895
OUT["table_implied_point"] = {"Omega_m": Om, "h": h, "omega_m": Om * h * h, "omega_b": wb,
                              "paper_formulas": pp, "camb_theta_star": cc,
                              "table_means": {"R": float(X.CMB_D[0]), "lA": float(X.CMB_D[1])},
                              "camb_minus_table_R_over_sigmaR": (cc["R"] - X.CMB_D[0]) / sigR,
                              "camb_minus_table_lA_over_sigmalA": (cc["lA"] - X.CMB_D[1]) / sigA,
                              "chi2_cmb_paper_formulas": X.chi2_cmb_from(pp["R"], pp["lA"], wb),
                              "chi2_cmb_camb": X.chi2_cmb_from(cc["R"], cc["lA"], wb),
                              "note": "paper-formula lA reproduces the table mean by construction of the solve; the CAMB theta*-based lA at the same parameters is offset, which is why the primary fit uses the paper's own definition"}
# --- component comparison at that point: paper r_s vs CAMB r_*
Or = X.Or_paper(Om, h); c_H0 = X.C / (100 * h)
zs = cc["zstar"]
E2 = X.E2_fn(Om, Or)
OUT["components_at_camb_zstar"] = {"zstar_camb": zs, "zstar_HS_code": pp["zstar"],
    "DM_paper_formula_Mpc": float(c_H0 * X.comoving_over_c_H0(E2, zs)), "DM_camb_Mpc": cc["DMstar"],
    "rs_paper_formula_Mpc": float(c_H0 * X.rs_over_c_H0(Or, Om, 1 - Om - Or, -1, 0, wb, zs)), "rs_camb_Mpc": cc["rstar"]}

# --- Planck-like point quoted in the paper's (commented) TeX line 78: H0=67.36, Om h^2=0.1430, ob h^2=0.02237
h2 = 0.6736; wm = 0.1430; Om2 = wm / h2 ** 2
OUT["paper_line78_point"] = {"H0": 67.36, "omega_m": wm, "omega_b": 0.02237,
                             "paper_formulas": X.cmb_paper(Om2, h2, 0.02237), "camb_theta_star": X.camb_exact(Om2, h2, 0.02237),
                             "note": "TeX line 78 is commented out and gives the TT,TE,EE+lowE+lensing values; used only as a second reference point"}
# --- covariance cross-check: inverse of published inverse covariance vs table sigma and correlations
C_pub = np.linalg.inv(X.CMB_INVC)
sig = np.sqrt(np.diag(C_pub)); corr = C_pub / np.outer(sig, sig)
C_tab = np.outer(X.CMB_SIG_TAB, X.CMB_SIG_TAB) * X.CMB_CORR_TAB
OUT["covariance_check"] = {"sigma_from_published_invcov": sig.tolist(), "sigma_table": X.CMB_SIG_TAB.tolist(),
                           "corr_from_published_invcov": corr.tolist(), "corr_table": X.CMB_CORR_TAB.tolist(),
                           "chi2_diff_at_1sigma_R_offset_published_vs_table_cov": float(
                               np.array([0.0046, 0, 0]) @ X.CMB_INVC @ np.array([0.0046, 0, 0]) - np.array([0.0046, 0, 0]) @ np.linalg.inv(C_tab) @ np.array([0.0046, 0, 0]))}
json.dump(OUT, open(HERE / "x4_camb_check.json", "w"), indent=1)
print(json.dumps(OUT, indent=1))
