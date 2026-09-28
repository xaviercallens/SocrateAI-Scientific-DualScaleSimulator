"""POST HOC (not pre-stated): (a) offset-1 vs offset-2 mean contact in the observed E. coli and GM12878 matrices
(the linear null sets offset-1 = E(2) for every case, the observed matrices keep measured offset-1 except Caulobacter);
(b) exercise the pipeline's betti_curve and euler_curve on the circle control and the E. coli alpha diagram.
Command: python posthoc_check2.py"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tda_common import PIPE, alpha_persistence, betti_curve, euler_curve, classical_mds, dump
import test2_hic as T2

HERE = os.path.dirname(os.path.abspath(__file__))
cases, _ = T2.build_cases()
out = {"NOTE": "POST HOC; does not change any verdict"}
for k in ("ecoli_full", "gm12878_chr1q", "caulo_full"):
    C = cases[k][0]
    e1, e2 = float(np.mean(np.diagonal(C, 1))), float(np.mean(np.diagonal(C, 2)))
    out[f"{k}_E1_E2_ratio"] = {"E1": e1, "E2": e2, "E1_over_E2": e1 / e2}


def curves(xyz, label):
    res, by_dim = alpha_persistence(np.ascontiguousarray(xyz, float), float("inf"), label, None)
    deaths = [d for dim in (0, 1, 2) for (_, d) in by_dim[dim] if np.isfinite(d)]
    r_trunc = float(max(deaths))
    grid = np.linspace(0, r_trunc, 21)
    chi, b0, b1, b2 = euler_curve(by_dim, grid, r_trunc)
    b1b = betti_curve(by_dim, 1, grid, r_trunc)
    return {"r_grid": grid.tolist(), "beta0": b0.tolist(), "beta1": b1.tolist(), "beta2": b2.tolist(), "euler": chi.tolist(),
            "betti_curve_dim1_equals_euler_component": bool(np.array_equal(b1, b1b)),
            "max_beta1": float(b1.max())}


out["circle_n200_curves"] = curves(PIPE.build_circle_cloud(200, noise_sigma=0.05, seed=0), "circle")
X3, _ = classical_mds(T2.to_dist(cases["ecoli_full"][0]), 3)
out["ecoli_mds3_curves"] = curves(X3, "ecoli")
dump(out, os.path.join(HERE, "posthoc_check2.json"))
print({k: v for k, v in out.items() if "E1" in k})
print("circle beta1", out["circle_n200_curves"]["beta1"], "\necoli beta1", out["ecoli_mds3_curves"]["beta1"])
