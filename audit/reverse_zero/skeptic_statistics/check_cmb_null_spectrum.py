#!/usr/bin/env python3
"""Skeptic (statistics lens) check of E5 CMB TDA: do the Gaussian null sims share the
DATA's effective smoothing / small-scale spectrum, and does fixing any mismatch change the
6 null-test p-values?

Step A (spectrum check): data masked pseudo-Cl vs mean masked pseudo-Cl of N_CHECK sims drawn
  exactly as cmb_tda.run_null draws them (synfast(Cl_pseudo/fsky, lmax=2*nside), same mask,
  same mean-subtraction), up to lmax=3*nside-1. Band ratios + z (sim scatter), and the
  spectral moment <l(l+1)> that sets peak/loop densities.
Step B (recalibrated null): Cl_corr iterated (<=4 times) as Cl_corr * R(l), R = data/sim pseudo-Cl ratio
  (running mean over 11 ell), extended to lmax=3*nside-1, until every band has |z|<3 (gate); only then
  N_NULL fresh sims through the SAME betti_curves_from_topology + coarse_stats code of
  cmb_tda.py; report the 6 p-values next to the committed ones.
Seeds: step A 424242+i; recheck 525252+1000*iter+i; step B null 626262+i. All numpy legacy seeding,
  as in cmb_tda.py.
Command (from worktree root):
  prlimit --as=8589934592 -- /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python \
    audit/reverse_zero/skeptic_statistics/check_cmb_null_spectrum.py --n-check 20 --n-null 100
"""
import argparse, json, sys, time
from pathlib import Path
import numpy as np, healpy as hp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "E5-cmb-tda"))
import cmb_tda as C  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--n-check", type=int, default=20)
ap.add_argument("--n-null", type=int, default=100)
ap.add_argument("--out", default=str(HERE / "check_cmb_null_spectrum.json"))
a = ap.parse_args()

m, mask, _ = C.load_data_map_and_mask()
cl, fsky = C.estimate_cl(m, mask, C.NSIDE_WORK)
L = 3 * C.NSIDE_WORK - 1
ell = np.arange(L + 1)
BANDS = [(2, 10), (11, 30), (31, 60), (61, 100), (101, 150), (151, 200), (201, 256), (257, 383)]


def pcl(x):
    y = np.where(mask > 0, x - x[mask > 0].mean(), 0.0)
    return hp.anafast(y, lmax=L)


def mom(c):
    w = (2 * ell + 1) * c
    return float((w * ell * (ell + 1)).sum() / w.sum())


def compare(d, sims):
    mu = sims.mean(0)
    out = []
    for lo, hi in BANDS:
        sl = slice(lo, hi + 1)
        bs = sims[:, sl].sum(1)
        out.append({"ell": [lo, hi], "data_over_sim": float(d[sl].sum() / mu[sl].sum()),
                    "z": float((d[sl].sum() - bs.mean()) / bs.std())})
    return out, mom(mu), float(np.std([mom(s) for s in sims]))


t0 = time.time()
d = pcl(m)
simsA = []
for i in range(a.n_check):
    np.random.seed(424242 + i)
    simsA.append(pcl(hp.synfast(cl, nside=C.NSIDE_WORK, new=True)))
simsA = np.array(simsA)
bandsA, momA, momA_sd = compare(d, simsA)
res = {"fsky": fsky, "lmax_pipeline_cl": len(cl) - 1, "lmax_check": L, "n_check": a.n_check,
       "stepA_pipeline_null": {"bands": bandsA, "mean_l(l+1)_data": mom(d), "mean_l(l+1)_sims": momA,
                               "mean_l(l+1)_sims_std": momA_sd}}

# Step B: recalibrate iteratively (v2: a single multiplicative step overshot ell>256 by ~4x because the
# pipeline sims carry no input power there; iterate until every band is within |z|<3 or MAX_ITER reached).
MAX_ITER = 4
k = 11
cl_corr = np.zeros(L + 1)
cl_corr[:len(cl)] = cl
cl_corr[len(cl):] = cl[-1]
sims_prev = simsA
iters = []
for it in range(MAX_ITER):
    ratio = d / np.where(sims_prev.mean(0) > 0, sims_prev.mean(0), np.nan)
    ratio[:2] = 1.0
    ratio = np.nan_to_num(ratio, nan=1.0)
    ratio_s = np.convolve(np.pad(ratio, k // 2, mode="edge"), np.ones(k) / k, mode="valid")
    cl_corr = cl_corr * ratio_s
    cl_corr[:2] = 0.0
    simsR = []
    for i in range(a.n_check):
        np.random.seed(525252 + 1000 * it + i)
        simsR.append(pcl(hp.synfast(cl_corr, nside=C.NSIDE_WORK, new=True)))
    simsR = np.array(simsR)
    bandsR, momR, momR_sd = compare(d, simsR)
    iters.append({"iter": it, "bands": bandsR, "mean_l(l+1)_sims": momR, "mean_l(l+1)_sims_std": momR_sd})
    sims_prev = simsR
    print(f"recal iter {it}: max|z|={max(abs(b['z']) for b in bandsR):.2f}", flush=True)
    if max(abs(b["z"]) for b in bandsR) < 3.0:
        break
gate_ok = max(abs(b["z"]) for b in iters[-1]["bands"]) < 3.0
res["stepB_recalibrated_spectrum_check"] = {"iterations": iters, "gate_all_bands_abs_z_lt_3": bool(gate_ok),
    "seeds": "525252+1000*iter+i",
    "method": "Cl_corr_{k+1} = Cl_corr_k * running-mean(11) of data/sim pseudo-Cl ratio; start = pipeline Cl extended flat to lmax=383"}
if not gate_ok:
    res["stepB_recalibrated_null_pvalues"] = "NOT RUN: recalibrated spectrum failed its own |z|<3 band gate"
    json.dump(res, open(a.out, "w"), indent=1)
    print(json.dumps(res, indent=1))
    sys.exit(0)
print(json.dumps(res, indent=1), flush=True)

nu = np.linspace(C.NU_MIN, C.NU_MAX, C.NU_STEP_GRID)
unm, E, T = C.build_topology(mask, C.NSIDE_WORK)
dat = {True: C.betti_curves_from_topology(m, unm, E, T, nu, sublevel=True),
       False: C.betti_curves_from_topology(m, unm, E, T, nu, sublevel=False)}
curves = {True: ([], [], []), False: ([], [], [])}
for i in range(a.n_null):
    np.random.seed(626262 + i)
    s = hp.synfast(cl_corr, nside=C.NSIDE_WORK, new=True)
    for sub in (True, False):
        b0, b1, ch = C.betti_curves_from_topology(s, unm, E, T, nu, sublevel=sub)
        curves[sub][0].append(b0); curves[sub][1].append(b1); curves[sub][2].append(ch)
    if (i + 1) % 10 == 0:
        print(f"null {i+1}/{a.n_null} t={time.time()-t0:.0f}s", flush=True)

committed = json.load(open(HERE.parent / "E5-cmb-tda" / "null_report.json"))
pv = {}
for sub, name in ((True, "sublevel"), (False, "superlevel")):
    for j, stat in enumerate(("b0", "b1", "euler_chi")):
        st = C.coarse_stats(np.array(curves[sub][j]), dat[sub][j])
        old = committed[f"statistic_{name}"][stat]
        resid = np.array(st["standardized_residual_full_curve"], dtype=float)
        pv[f"{name}_{stat}"] = {"recalibrated_p_chi2": st["p_value_chi2_survival"], "recalibrated_rank_p": st["empirical_rank_p"],
                                "committed_p_chi2": old["p_value_chi2_survival"], "committed_rank_p": old["empirical_rank_p"],
                                "recalibrated_max_abs_resid": float(np.nanmax(np.abs(resid)))}
res["stepB_recalibrated_null_pvalues"] = pv
res["n_null_recalibrated"] = a.n_null
res["bonferroni_min_recal_p_chi2_x6"] = float(min(1.0, 6 * min(v["recalibrated_p_chi2"] for v in pv.values())))
res["runtime_sec"] = time.time() - t0
json.dump(res, open(a.out, "w"), indent=1)
print(json.dumps(pv, indent=1))
