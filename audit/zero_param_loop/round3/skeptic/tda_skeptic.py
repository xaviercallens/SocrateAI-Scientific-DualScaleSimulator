"""Round-3 skeptic TDA check (independent of round3/tda.py's reporting logic).
1. gudhi really imported: version printed into the report.
2. Known-answer circle at seed 42 (same generator recipe) and seed 44: what does
   st.betti_numbers() (the field round3/tda.py reports as 'betti_numbers') say, vs the
   persistent H1 bar? -> tests whether the reported Betti numbers can see a circle at all.
3. Model cloud (round3/sweep.csv stable rows, same columns/standardisation/PCA as tda.py):
   b0 at final scale with sparse=0.2 vs sparse=None, and the H0 death values of the
   largest gaps -> is b0=2 a sparse artefact or a real gap?
4. Null vs real discrimination: 10 Poisson box-matched nulls (seeds 100..109) and
   10 real re-subsamples (seeds 200..209); H1 max persistence and bottleneck(H1) to a
   fixed real subsample. Real is distinguishable only if the two distributions separate.
Writes tda_skeptic_report.json next to this file. Seeds fixed."""
import json, math, os, sys
import numpy as np, pandas as pd
import gudhi
from scipy.spatial.distance import pdist
HERE = os.path.dirname(os.path.abspath(__file__))
R3 = os.path.abspath(os.path.join(HERE, ".."))
ROOT = os.path.abspath(os.path.join(R3, "..", "..", ".."))
sys.path.insert(0, os.path.join(R3, "..", "round2"))
import tda_gudhi as tg  # reused only for rescale_unit_diam, pca_reduce, comoving distance

def persist(X, sparse=0.2):
    X2, _ = tg.rescale_unit_diam(X)
    d = pdist(X2); mel = 2.0 * float(np.median(d))
    st = gudhi.RipsComplex(points=X2, max_edge_length=mel, sparse=sparse).create_simplex_tree(max_dimension=2)
    diag = st.persistence(homology_coeff_field=2, min_persistence=0.0)
    bars = {k: np.array([[b, dth] for dm, (b, dth) in diag if dm == k]) for k in (0, 1, 2)}
    return st, bars, mel

def h1_stats(bars):
    h = bars[1]
    if len(h) == 0: return {"n": 0, "max_pers": 0.0, "ratio": None}
    p = np.sort(h[:, 1] - h[:, 0])[::-1]
    return {"n": int(len(p)), "max_pers": float(p[0]), "ratio": float(p[0] / p[1]) if len(p) > 1 and p[1] > 0 else None,
            "top_bar": h[np.argmax(h[:, 1] - h[:, 0])].tolist()}

def circle(seed, n=300, s=0.05):
    r = np.random.RandomState(seed); th = r.uniform(0, 2 * math.pi, n); rr = 1 + r.normal(0, s, n)
    return np.column_stack([rr * np.cos(th), rr * np.sin(th)])

out = {"gudhi_version": gudhi.__version__, "gudhi_file": gudhi.__file__}
# 2. circle
out["circle"] = {}
for seed in (42, 44):
    st, bars, mel = persist(circle(seed))
    hs = h1_stats(bars); mid = 0.5 * (hs["top_bar"][0] + hs["top_bar"][1])
    out["circle"][seed] = {"final_scale_betti_numbers_as_reported_by_tda_py": st.betti_numbers(),
                           "max_edge_length": mel, "h1": hs,
                           "persistent_betti_at_mid_of_top_H1_bar": st.persistent_betti_numbers(mid, mid)}
# 3. model cloud
df = pd.read_csv(os.path.join(R3, "sweep.csv"))
cols = ["log10_ssf", "log10_pcr"] + [f"gamma_theta_{i}" for i in range(15)] + ["pta_max_deviation_from_hd"]
X = df[df["numerically_stable"] == True][cols].to_numpy(float); X = X[np.all(np.isfinite(X), 1)]
rng = np.random.RandomState(42); X = X[rng.choice(len(X), 300, replace=False)]
Xz = (X - X.mean(0)) / np.where(X.std(0) < 1e-12, 1, X.std(0))
Xp, k, frac = tg.pca_reduce(Xz, max_dim=8, var_target=0.95)
out["model"] = {"pca_k": k}
for sp in (0.2, None):
    st, bars, mel = persist(Xp, sparse=sp)
    h0 = bars[0]; fin = np.sort(h0[np.isfinite(h0[:, 1]), 1])[::-1]
    out["model"][f"sparse={sp}"] = {"final_scale_betti": st.betti_numbers(), "max_edge_length": mel,
                                    "n_infinite_H0": int(np.sum(~np.isfinite(h0[:, 1]))),
                                    "largest_finite_H0_deaths": fin[:5].tolist(), "h1": h1_stats(bars)}
# 4. null vs real
path = os.path.join(ROOT, "data", "real", "dark_energy", "pantheon_plus_sh0es.dat")
p = pd.read_csv(path, sep=r"\s+", engine="python")
p = p[(p["zHD"] > 0.005) & (p["zHD"] <= 2.4) & np.isfinite(p["RA"]) & np.isfinite(p["DEC"])]
ra, dec = np.radians(p["RA"].to_numpy()), np.radians(p["DEC"].to_numpy())
rc = tg.comoving_distance_flat_lcdm(p["zHD"].to_numpy())
R = np.column_stack([rc * np.cos(dec) * np.cos(ra), rc * np.cos(dec) * np.sin(ra), rc * np.sin(dec)])
def sub(seed): return R[np.random.RandomState(seed).choice(len(R), 300, replace=False)]
ref = sub(42); _, bref, _ = persist(ref)
def summarize(clouds):
    rows = []
    for Y in clouds:
        _, b, _ = persist(Y)
        rows.append({"h1_max_pers": h1_stats(b)["max_pers"], "h1_n": h1_stats(b)["n"],
                     "bottleneck_H1_to_ref": float(gudhi.bottleneck_distance(b[1], bref[1])),
                     "bottleneck_H0_to_ref": float(gudhi.bottleneck_distance(b[0][np.isfinite(b[0][:, 1])], bref[0][np.isfinite(bref[0][:, 1])]))})
    return rows
reals = summarize([sub(s) for s in range(200, 210)])
nulls = []
for s in range(100, 110):
    Y = sub(200 + s - 100); lo, hi = Y.min(0), Y.max(0)
    nulls.append(np.random.RandomState(s).uniform(lo, hi, size=Y.shape))
nulls = summarize(nulls)
def rng_of(rows, key): v = [r[key] for r in rows]; return [float(min(v)), float(np.median(v)), float(max(v))]
out["null_vs_real"] = {"ref": "real subsample seed 42 (N=300)",
    "real_resubsamples_min_med_max": {k: rng_of(reals, k) for k in reals[0]},
    "poisson_nulls_min_med_max": {k: rng_of(nulls, k) for k in nulls[0]}}
for k in ("h1_max_pers", "bottleneck_H0_to_ref", "bottleneck_H1_to_ref"):
    a, b = rng_of(reals, k), rng_of(nulls, k)
    out["null_vs_real"][f"separated_on_{k}"] = bool(a[2] < b[0] or b[2] < a[0])
json.dump(out, open(os.path.join(HERE, "tda_skeptic_report.json"), "w"), indent=1, default=str)
print(json.dumps(out, indent=1, default=str))
