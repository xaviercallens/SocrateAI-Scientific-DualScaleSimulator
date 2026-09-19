"""Round-3 skeptic: independent re-implementation (no import of round2/decisive_experiment)
of the pre-registered decisive experiment, plus the frozen-SN-offset variant, plus harness
probes of mu_sym and c4_pta_product (x10, x0.1) at the nominal point.
Data: data/real/dark_energy/{desi_2024_bao_all.txt, desi_2024_bao_all_cov.txt,
pantheon_plus_sh0es.dat} (sha256 verified in earlier rounds). Writes chi2_independent.json."""
import json, os, sys
import numpy as np, pandas as pd
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import minimize, minimize_scalar
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(os.path.join(HERE, *[".."] * 4))
D = os.path.join(ROOT, "data", "real", "dark_energy")
OMEGA_LAMBDA = 0.68885          # LeanMaster DarkEnergyScale.lean:84 (Planck 2018, tier L)
HUBBLE_RADIUS_M = 1.3672e26     # LeanMaster SelfDualCutoff.lean:100
C_KM_S = 299792.458; MPC_M = 3.0856775814913673e22
H0_LEAN = C_KM_S / (HUBBLE_RADIUS_M / MPC_M)
rows = [l.split() for l in open(os.path.join(D, "desi_2024_bao_all.txt")) if l.strip() and not l.startswith("#")]
bz = np.array([float(r[0]) for r in rows]); bv = np.array([float(r[1]) for r in rows]); bq = [r[2] for r in rows]
Cinv = np.linalg.inv(np.loadtxt(os.path.join(D, "desi_2024_bao_all_cov.txt")))
sn = pd.read_csv(os.path.join(D, "pantheon_plus_sh0es.dat"), sep=r"\s+")
sn = sn[(sn.zHD > 0.01) & (sn.zHD <= 2.4) & np.isfinite(sn.MU_SH0ES) & (sn.MU_SH0ES_ERR_DIAG > 0)]
sz, smu, sw = sn.zHD.to_numpy(), sn.MU_SH0ES.to_numpy(), 1 / sn.MU_SH0ES_ERR_DIAG.to_numpy() ** 2
zg = np.linspace(0, 2.5, 25001)
def E_lcdm(Om): return np.sqrt(Om * (1 + zg) ** 3 + 1 - Om)
def E_cpl(Om, w0, wa): return np.sqrt(Om * (1 + zg) ** 3 + (1 - Om) * (1 + zg) ** (3 * (1 + w0 + wa)) * np.exp(-3 * wa * zg / (1 + zg)))
def chi2(E, sn_offset=None):
    dc = cumulative_trapezoid(1 / E, zg, initial=0)          # D_C * H0 / c
    DM = np.interp(bz, zg, dc); DH = np.interp(bz, zg, 1 / E)
    base = np.where([q == "DM_over_rs" for q in bq], DM, np.where([q == "DH_over_rs" for q in bq], DH, (bz * DM ** 2 * DH) ** (1 / 3)))
    u = (base @ Cinv @ bv) / (base @ Cinv @ base); r = bv - u * base; cb = float(r @ Cinv @ r)
    m0 = 5 * np.log10((1 + sz) * np.interp(sz, zg, dc))
    off = float(np.sum(sw * (smu - m0)) / np.sum(sw)) if sn_offset is None else sn_offset
    cs = float(np.sum(sw * (smu - m0 - off) ** 2))
    return {"chi2": cb + cs, "bao": cb, "sn": cs, "sn_offset": off, "H0_implied_by_offset": C_KM_S / 10 ** ((off - 25) / 5)}
out = {"H0_from_hubbleRadius_m": H0_LEAN, "n_bao": int(len(bv)), "n_sn": int(len(sz))}
Om_L = 1 - OMEGA_LAMBDA
out["frozen_lean"] = chi2(E_lcdm(Om_L))
r = minimize_scalar(lambda o: chi2(E_lcdm(o))["chi2"], bounds=(0.05, 0.95), method="bounded", options={"xatol": 1e-9})
out["fitted_lcdm"] = {"Om": float(r.x), **chi2(E_lcdm(r.x))}
out["delta_chi2_frozen_minus_fitted"] = out["frozen_lean"]["chi2"] - out["fitted_lcdm"]["chi2"]
c = minimize(lambda p: chi2(E_cpl(Om_L, *p))["chi2"], [-1, 0], method="Nelder-Mead", options={"xatol": 1e-7, "fatol": 1e-8, "maxiter": 4000})
out["cpl_Om_frozen"] = {"w0": float(c.x[0]), "wa": float(c.x[1]), **chi2(E_cpl(Om_L, *c.x))}
c2 = minimize(lambda p: chi2(E_cpl(*p))["chi2"], [0.3, -1, 0], method="Nelder-Mead", options={"xatol": 1e-7, "fatol": 1e-8, "maxiter": 8000})
out["cpl_Om_free"] = {"Om": float(c2.x[0]), "w0": float(c2.x[1]), "wa": float(c2.x[2]), **chi2(E_cpl(*c2.x))}
out["delta_chi2_frozen_minus_cpl_Om_frozen"] = out["frozen_lean"]["chi2"] - out["cpl_Om_frozen"]["chi2"]
out["controls"] = {f"Om_L={1-o:.2f}": chi2(E_lcdm(o))["chi2"] - out["fitted_lcdm"]["chi2"] for o in (0.5, 0.28, 0.35)}
# frozen SN offset: MU_SH0ES carries the SH0ES Cepheid M_B calibration inside the data file;
# freezing H0 = c/hubbleRadius_m fixes offset = 5 log10(c/H0 [Mpc]) + 25.
off_frozen = 5 * np.log10(C_KM_S / H0_LEAN) + 25
out["frozen_lean_sn_offset_frozen"] = {"offset": off_frozen, **chi2(E_lcdm(Om_L), sn_offset=off_frozen)}
out["delta_chi2_offset_frozen_minus_offset_fitted"] = out["frozen_lean_sn_offset_frozen"]["chi2"] - out["frozen_lean"]["chi2"]
out["threshold_preregistered"] = 11.8
out["reject_frozen_offset_fitted"] = bool(out["delta_chi2_frozen_minus_fitted"] > 11.8)
out["reject_frozen_offset_frozen"] = bool(out["delta_chi2_offset_frozen_minus_offset_fitted"] > 11.8)
# harness probes at nominal point
sys.path.insert(0, os.path.join(ROOT, "scripts")); import param_loop_sim as pls
def probe(**kw):
    p = dict(pls.DEFAULT_PARAMS); p.update(kw); e = pls.evaluate_point(p)
    return {"ssf": e["screening"]["screening_suppression_factor"], "pcr": e["screening"]["phi_center_ratio"],
            "stable": e["screening"]["numerically_stable"], "pta_max_dev": e["pta"]["max_deviation_from_hd"],
            "w0_cpl": e["dark_energy"]["w0_cpl_latetime"]}
out["probes"] = {"nominal": probe(), "mu_sym_x10": probe(mu_sym=10.0), "mu_sym_x0.1": probe(mu_sym=0.1),
                 "c4_pta_product_x10": probe(c4_c0_ratio=160.7), "c4_pta_product_x0.1": probe(c4_c0_ratio=1.607)}
out["probes_note"] = "chi2 (BAO+SN) is computed from frozen LCDM E(z), which takes neither mu_sym nor c4_pta_product; no screening or PTA-angular dataset exists (ABSENT), so chi2 is identical at every probe."
json.dump(out, open(os.path.join(HERE, "chi2_independent.json"), "w"), indent=1, default=float)
print(json.dumps(out, indent=1, default=float))
