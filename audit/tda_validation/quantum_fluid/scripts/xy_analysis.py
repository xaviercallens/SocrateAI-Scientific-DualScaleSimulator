"""Test 1 analysis: helicity-modulus cross-check (E1a), logistic-regression TDA
crossover (E1b), training-window control (E1c), vortex/Betti identity (E1d),
site-shuffle negative control (E1e), alpha-path pairing statistic (E1f).
Inputs: DATA_ROOT/xy/L*/T*.npz (xy_mc.py) and DATA_ROOT/xy_tda/L*/T*.npz
(xy_tda_features.py). Everything as pre-stated in expectations.json.
Command: prlimit --as=8589934592 -- .venv-tda/bin/python xy_analysis.py
Output: ../results/xy_analysis.json
"""
import glob
import json
import os
import sys

import numpy as np
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qf_common as qc  # noqa: E402

TBKT = 0.89289


def crossing(T, y, level):
    """first T where y - level changes sign from + to -, linear interpolation"""
    d = np.asarray(y) - level
    for k in range(len(T) - 1):
        if d[k] > 0 and d[k + 1] <= 0:
            return float(T[k] + (T[k + 1] - T[k]) * d[k] / (d[k] - d[k + 1]))
    return None


def crossing_up(T, y, level):
    d = np.asarray(y) - level
    for k in range(len(T) - 1):
        if d[k] < 0 and d[k + 1] >= 0:
            return float(T[k] + (T[k + 1] - T[k]) * (-d[k]) / (d[k + 1] - d[k]))
    return None


def helicity(L):
    rows = []
    for fn in sorted(glob.glob(os.path.join(qc.DATA_ROOT, "xy", f"L{L}", "T*.npz"))):
        d = np.load(fn)
        rows.append((float(d["T"]), float(d["upsilon"]), float(d["upsilon_err"]), float(d["vortex_density"]),
                     float(d["vortex_density_err"]), float(d["energy"]), float(d["tau_int_energy_meas_units"]), int(d["nmeas"])))
    rows.sort()
    a = np.array(rows)
    T, U, Ue = a[:, 0], a[:, 1], a[:, 2]
    tc = crossing(T, U - 2 * T / np.pi, 0.0)
    # bootstrap error of the crossing from the jackknife errors (Gaussian resampling, seed 0)
    rng = np.random.default_rng(0)
    bs = [crossing(T, U + rng.normal(0, Ue) - 2 * T / np.pi, 0.0) for _ in range(500)]
    bs = np.array([b for b in bs if b is not None])
    return {"L": L, "T": T.tolist(), "upsilon": U.tolist(), "upsilon_err": Ue.tolist(),
            "vortex_density": a[:, 3].tolist(), "vortex_density_err": a[:, 4].tolist(), "energy": a[:, 5].tolist(),
            "tau_int_energy_meas_units_max": float(a[:, 6].max()), "nmeas": int(a[0, 7]),
            "T_Ups_crossing": tc, "T_Ups_crossing_err": float(bs.std()) if bs.size else None}


def load_feats(L):
    out = []
    for fn in sorted(glob.glob(os.path.join(qc.DATA_ROOT, "xy_tda", f"L{L}", "T*.npz"))):
        out.append(dict(np.load(fn)))
    out.sort(key=lambda d: float(d["T"]))
    return out


def classifier_crossing(feats, L, low, high, key0="b0", key1="b1", seed=0):
    X, y, Ts = [], [], []
    for d in feats:
        T = float(d["T"])
        F = np.hstack([d[key0], d[key1]]).astype(float) / (L * L)
        for row in F:
            X.append(row); Ts.append(T)
    X = np.array(X); Ts = np.array(Ts)
    lowm = (Ts >= low[0] - 1e-9) & (Ts <= low[1] + 1e-9)
    highm = (Ts >= high[0] - 1e-9) & (Ts <= high[1] + 1e-9)
    idx = np.where(lowm | highm)[0]
    rng = np.random.default_rng(seed)
    rng.shuffle(idx)
    tr = idx[: len(idx) // 2]; te = idx[len(idx) // 2:]
    lab = highm.astype(int)
    clf = make_pipeline(StandardScaler(), LogisticRegression(C=1.0, max_iter=5000))
    clf.fit(X[tr], lab[tr])
    acc = float(clf.score(X[te], lab[te]))
    P = clf.predict_proba(X)[:, 1]
    Tu = np.unique(Ts)
    Pm = np.array([P[Ts == t].mean() for t in Tu])
    tc = crossing_up(Tu, Pm, 0.5)
    # bootstrap over configurations within each T (seed 1) for the crossing uncertainty
    rng2 = np.random.default_rng(1)
    bs = []
    for _ in range(200):
        Pb = np.array([rng2.choice(P[Ts == t], size=(Ts == t).sum()).mean() for t in Tu])
        c = crossing_up(Tu, Pb, 0.5)
        if c is not None:
            bs.append(c)
    return {"window_low": low, "window_high": high, "test_accuracy": acc, "T": Tu.tolist(), "mean_P_high": Pm.tolist(),
            "T_TDA": tc, "T_TDA_boot_std": float(np.std(bs)) if bs else None}


def main():
    res = {"T_BKT_ref": TBKT}
    hel = {}
    for L in (32, 64, 128):
        if glob.glob(os.path.join(qc.DATA_ROOT, "xy", f"L{L}", "T*.npz")):
            hel[L] = helicity(L)
    res["E1a_helicity"] = {str(k): v for k, v in hel.items()}
    Ls = [L for L in hel if hel[L]["T_Ups_crossing"] is not None]
    if len(Ls) >= 2:
        xs = np.array([1 / np.log(L) ** 2 for L in Ls]); ys = np.array([hel[L]["T_Ups_crossing"] for L in Ls])
        A = np.c_[np.ones_like(xs), xs]
        coef, *_ = np.linalg.lstsq(A, ys, rcond=None)
        res["E1a_extrapolation"] = {"L": Ls, "T_Ups": ys.tolist(), "fit": "T_Ups(L) = T_inf + c/(ln L)^2",
                                    "T_inf": float(coef[0]), "c": float(coef[1]),
                                    "residuals": (ys - A @ coef).tolist()}
        ok = 0.863 <= coef[0] <= 0.923 and (hel[max(Ls)]["T_Ups_crossing"] <= hel[min(Ls)]["T_Ups_crossing"])
        res["E1a_pass"] = bool(ok and 128 in Ls and 32 in Ls)
    # TDA
    tda = {}
    for L in (32, 64, 128):
        feats = load_feats(L)
        n_mc = len(glob.glob(os.path.join(qc.DATA_ROOT, "xy", f"L{L}", "T*.npz")))
        if not feats or len(feats) < n_mc:
            print(f"L={L}: TDA features incomplete ({len(feats)}/{n_mc} temperatures), skipped")
            continue
        r = {"n_T": len(feats), "n_cfg_per_T": int(len(feats[0]["b0"]))}
        r["window_A"] = classifier_crossing(feats, L, (0.40, 0.60), (1.30, 1.60))
        r["window_B"] = classifier_crossing(feats, L, (0.40, 0.75), (1.30, 1.60))
        r["window_C"] = classifier_crossing(feats, L, (0.40, 0.60), (1.10, 1.60))
        # E1d correlation N_v vs b0 at -3pi/4 (index 4 of the 33-point grid)
        corr = {}
        for d in feats:
            T = float(d["T"])
            if any(abs(T - t) < 1e-6 for t in (0.70, 0.75, 0.80, 0.85, 0.90)):
                nv = d["nv"]; b = d["b0"][:, 4]
                rho = stats.spearmanr(nv, b).correlation if nv.std() > 0 and b.std() > 0 else None
                corr[f"{T:.3f}"] = {"spearman_rho": None if rho is None else float(rho), "mean_nv": float(nv.mean()),
                                    "mean_b0_m3pi4": float(b.mean())}
        r["E1d_corr"] = corr
        # E1e shuffle control
        sh = {}
        for d in feats:
            if len(d["b0_shuf"]) == 0:
                continue
            T = float(d["T"])
            n = len(d["b0_shuf"])
            real = d["b0"][:n, 16] / L ** 2; shuf = d["b0_shuf"][:, 16] / L ** 2
            diff = shuf - real
            se = diff.std(ddof=1) / np.sqrt(n)
            sh[f"{T:.3f}"] = {"real_b0_theta_le_0_per_site": float(real.mean()), "shuffled": float(shuf.mean()),
                              "diff_over_se": float(diff.mean() / se) if se > 0 else None,
                              "real_b1_theta_le_0_per_site": float(d["b1"][:n, 16].mean() / L ** 2),
                              "shuffled_b1": float(d["b1_shuf"][:, 16].mean() / L ** 2)}
        r["E1e_shuffle"] = sh
        if any(len(d["b0_shuf"]) for d in feats):
            fs = [d for d in feats if len(d["b0_shuf"])]
            r["E1e_classifier_on_shuffled_windowA"] = classifier_crossing(fs, L, (0.40, 0.60), (1.30, 1.60), "b0_shuf", "b1_shuf")
        shdir = os.path.join(qc.DATA_ROOT, "xy_tda", f"L{L}_shufonly")
        if os.path.isdir(shdir) and len(glob.glob(os.path.join(shdir, "T*.npz"))) == len(feats):
            fs = sorted([dict(np.load(f)) for f in glob.glob(os.path.join(shdir, "T*.npz"))], key=lambda d: float(d["T"]))
            r["E1e_classifier_on_shuffled_windowA"] = classifier_crossing(fs, L, (0.40, 0.60), (1.30, 1.60), "b0_shuf", "b1_shuf")
            r["E1e_classifier_on_shuffled_note"] = "EXTENSION: shuffle-only run (xy_tda_features.py --shuffle-only 50)"
        # E1f alpha pairing
        al = {}
        for d in feats:
            fsd = d["fshort"]; fr = d["fshort_rand"]
            ok = np.isfinite(fsd) & np.isfinite(fr)
            if ok.sum() >= 5:
                al[f"{float(d['T']):.3f}"] = {"n_cfg": int(ok.sum()), "f_short": float(fsd[ok].mean()),
                                              "f_short_rand": float(fr[ok].mean()),
                                              "S_pair": float(fsd[ok].mean() / fr[ok].mean()) if fr[ok].mean() > 0 else None}
        r["E1f_alpha_pairing"] = al
        tda[str(L)] = r
    res["tda"] = tda
    # verdicts
    v = {}
    if "32" in tda and "64" in tda:
        t32 = tda["32"]["window_A"]["T_TDA"]; t64 = tda["64"]["window_A"]["T_TDA"]
        br = all(t is not None and 0.80 <= t <= 1.05 for t in (t32, t64))
        if "128" in tda:
            t128 = tda["128"]["window_A"]["T_TDA"]
            br = br and t128 is not None and 0.80 <= t128 <= 1.05
        v["E1b_pass"] = bool(br and t64 <= t32 + 0.02)
    if "64" in tda:
        r = tda["64"]
        tA, tB, tC = r["window_A"]["T_TDA"], r["window_B"]["T_TDA"], r["window_C"]["T_TDA"]
        v["E1c_shift_B_minus_A"] = None if None in (tA, tB) else tB - tA
        v["E1c_shift_C_minus_A"] = None if None in (tA, tC) else tC - tA
        v["E1c_pass"] = bool(None not in (tA, tB, tC) and abs(tB - tA) <= 0.04 and abs(tC - tA) <= 0.04)
        rhos = [r["E1d_corr"].get(k, {}).get("spearman_rho") for k in ("0.700", "0.750", "0.800")]
        rhos = [x for x in rhos if x is not None]
        v["E1d_median_rho"] = float(np.median(rhos)) if rhos else None
        v["E1d_pass"] = bool(rhos and np.median(rhos) > 0.3)
        sh = {k: x for k, x in r["E1e_shuffle"].items() if float(k) <= 1.0 + 1e-9}
        v["E1e_min_diff_over_se_T_le_1"] = float(min(x["diff_over_se"] for x in sh.values())) if sh else None
        v["E1e_pass"] = bool(sh and all(x["diff_over_se"] is not None and x["diff_over_se"] > 5 for x in sh.values()))
        al = r["E1f_alpha_pairing"]
        s085 = al.get("0.850", {}).get("S_pair"); s150 = al.get("1.500", {}).get("S_pair")
        v["E1f_S_pair_0.85"] = s085; v["E1f_S_pair_1.50"] = s150
        v["E1f_pass"] = bool(s085 is not None and s150 is not None and s085 > 2 and s150 < 1.5)
    res["verdicts"] = v
    with open(os.path.join(qc.RESULTS, "xy_analysis.json"), "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps({"E1a": res.get("E1a_extrapolation"), "E1a_pass": res.get("E1a_pass"),
                      "Tups": {k: (h["T_Ups_crossing"], h["T_Ups_crossing_err"]) for k, h in res["E1a_helicity"].items()},
                      "TDA": {L: {w: (tda[L][w]["T_TDA"], tda[L][w]["T_TDA_boot_std"], tda[L][w]["test_accuracy"]) for w in ("window_A", "window_B", "window_C")} for L in tda},
                      "verdicts": v}, indent=1))


if __name__ == "__main__":
    main()
