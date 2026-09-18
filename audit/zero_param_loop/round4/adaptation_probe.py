"""Round-4 report-writer probe (exploratory, tier X unless stated). Writes adaptation_probe.json.
Grounds two numbers the round-4 brainstorm quoted without a file:
 (1) DATA-FIRST idea 1: self-dual vacuum rho_DE(a) with x = a/a*, four discrete variants, zero fitted
     shape parameters. chi2 against DESI 2024 BAO (12x12 cov) + Pantheon+ (1590 SNe, diagonal errors),
     same pipeline as round3/skeptic/chi2_independent.py (lines copied, not imported, because that file
     runs its whole experiment at import time). Variant B (mirror of A) is the negative control.
 (2) DERIVE target A: Legendre coefficients of the Hellings-Downs curve used in scripts/param_loop_sim.py
     (hd = 1.5 x ln x - 0.25 x + 0.5, x = (1-cos)/2), compared with the closed form
     a_l = (3/2)(2l+1)(l-2)!/(l+2)! for l>=2 and a_0 = a_1 = 0 (known-answer control).
Deterministic: no random calls (seed 42 set anyway)."""
import json, os, math
import numpy as np, pandas as pd
from scipy.integrate import cumulative_trapezoid, quad
from scipy.optimize import minimize
np.random.seed(42)
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(os.path.join(HERE, *[".."] * 3))
D = os.path.join(ROOT, "data", "real", "dark_energy")
OMEGA_LAMBDA = 0.68885          # LeanMaster DarkEnergyScale.lean:84 (tier L)
C_KM_S = 299792.458
rows = [l.split() for l in open(os.path.join(D, "desi_2024_bao_all.txt")) if l.strip() and not l.startswith("#")]
bz = np.array([float(r[0]) for r in rows]); bv = np.array([float(r[1]) for r in rows]); bq = [r[2] for r in rows]
Cinv = np.linalg.inv(np.loadtxt(os.path.join(D, "desi_2024_bao_all_cov.txt")))
sn = pd.read_csv(os.path.join(D, "pantheon_plus_sh0es.dat"), sep=r"\s+")
sn = sn[(sn.zHD > 0.01) & (sn.zHD <= 2.4) & np.isfinite(sn.MU_SH0ES) & (sn.MU_SH0ES_ERR_DIAG > 0)]
sz, smu, sw = sn.zHD.to_numpy(), sn.MU_SH0ES.to_numpy(), 1 / sn.MU_SH0ES_ERR_DIAG.to_numpy() ** 2
zg = np.linspace(0, 2.5, 25001); ag = 1 / (1 + zg)
def chi2(E):
    dc = cumulative_trapezoid(1 / E, zg, initial=0)
    DM = np.interp(bz, zg, dc); DH = np.interp(bz, zg, 1 / E)
    base = np.where([q == "DM_over_rs" for q in bq], DM, np.where([q == "DH_over_rs" for q in bq], DH, (bz * DM ** 2 * DH) ** (1 / 3)))
    u = (base @ Cinv @ bv) / (base @ Cinv @ base); r = bv - u * base; cb = float(r @ Cinv @ r)
    m0 = 5 * np.log10((1 + sz) * np.interp(sz, zg, dc))
    off = float(np.sum(sw * (smu - m0)) / np.sum(sw)); cs = float(np.sum(sw * (smu - m0 - off) ** 2))
    return {"chi2": cb + cs, "bao": cb, "sn": cs, "H0_implied_by_offset": C_KM_S / 10 ** ((off - 25) / 5)}
Om = 1 - OMEGA_LAMBDA
a_star_eq = (Om / OMEGA_LAMBDA) ** (1 / 3)
def f_A(x): return 2 * x / (1 + x ** 2)          # peaks at x=1, symmetric under x -> 1/x
def f_B(x): return (x + 1 / x) / 2               # mirror (1/f_A): negative control
def dlnf_dlnx(kind, x): return (1 - x ** 2) / (1 + x ** 2) if kind == "A" else (x ** 2 - 1) / (x ** 2 + 1)
def variant(kind, a_star):
    f = f_A if kind == "A" else f_B
    E = np.sqrt(Om * (1 + zg) ** 3 + OMEGA_LAMBDA * f(ag / a_star) / f(1 / a_star))
    x0 = 1 / a_star; eps = 1e-6
    w = lambda a: -1 - dlnf_dlnx(kind, a / a_star) / 3            # w = -1 - (1/3) dln rho/dln a
    w0 = w(1.0); wa = -(w(1 + eps) - w(1 - eps)) / (2 * eps)       # CPL: wa = -dw/da at a=1
    return {"kind": kind, "a_star": a_star, "w0": w0, "wa": wa, **chi2(E)}
out = {"a_star_matter_lambda_equality": a_star_eq}
lcdm = chi2(np.sqrt(Om * (1 + zg) ** 3 + OMEGA_LAMBDA)); out["frozen_lcdm"] = lcdm
out["variants"] = {"A_astar_eq": variant("A", a_star_eq), "B_astar_eq": variant("B", a_star_eq),
                   "C_A_astar_1": variant("A", 1.0), "D_B_astar_1": variant("B", 1.0)}
for v in out["variants"].values(): v["delta_chi2_vs_frozen_lcdm"] = v["chi2"] - lcdm["chi2"]
def E_cpl(w0, wa): return np.sqrt(Om * (1 + zg) ** 3 + OMEGA_LAMBDA * (1 + zg) ** (3 * (1 + w0 + wa)) * np.exp(-3 * wa * zg / (1 + zg)))
c = minimize(lambda p: chi2(E_cpl(*p))["chi2"], [-1, 0], method="Nelder-Mead", options={"xatol": 1e-7, "fatol": 1e-8, "maxiter": 4000})
out["cpl_Om_frozen"] = {"w0": float(c.x[0]), "wa": float(c.x[1]), **chi2(E_cpl(*c.x))}
out["cpl_Om_frozen"]["delta_chi2_vs_frozen_lcdm"] = out["cpl_Om_frozen"]["chi2"] - lcdm["chi2"]
out["fraction_of_cpl_gain_recovered_by_A"] = out["variants"]["A_astar_eq"]["delta_chi2_vs_frozen_lcdm"] / out["cpl_Om_frozen"]["delta_chi2_vs_frozen_lcdm"]
out["look_elsewhere_note"] = "4 discrete variants tried, A chosen after seeing the data; threshold 11.8 pre-registered in round 2."
# (2) Hellings-Downs Legendre coefficients
def hd(t):
    x = (1 - t) / 2
    return 1.5 * x * math.log(x) - 0.25 * x + 0.5 if x > 0 else 0.5
def P(l, t): return float(np.polynomial.legendre.Legendre.basis(l)(t))
hdc = {}
for l in range(0, 7):
    num = (2 * l + 1) / 2 * quad(lambda t: hd(t) * P(l, t), -1, 1, limit=200, epsabs=1e-13, epsrel=1e-12)[0]
    closed = 0.0 if l < 2 else 1.5 * (2 * l + 1) * math.factorial(l - 2) / math.factorial(l + 2)
    hdc[str(l)] = {"numeric": num, "closed_form": closed, "abs_diff": abs(num - closed)}
out["hd_legendre_coefficients"] = hdc
out["hd_l4_closed_form_rational"] = "3/80"
out["model_l4_coefficient_nominal"] = 3 / 80 + 0.08035
out["note_pta"] = ("param_loop_sim adds c4_pta_product*P4 on top of hd, so the model's total l=4 coefficient is "
                   "3/80 + c4_pta_product; c4_pta_product = 0 is the GR (massless transverse-traceless, isotropic) value. "
                   "The integral is tier B here (quadrature vs closed form); identifying the model's PTA sector with GR is tier C.")
json.dump(out, open(os.path.join(HERE, "adaptation_probe.json"), "w"), indent=1, default=float)
print(json.dumps({k: out[k] for k in ("a_star_matter_lambda_equality", "frozen_lcdm", "variants", "cpl_Om_frozen", "fraction_of_cpl_gain_recovered_by_A", "hd_legendre_coefficients")}, indent=1, default=float))
