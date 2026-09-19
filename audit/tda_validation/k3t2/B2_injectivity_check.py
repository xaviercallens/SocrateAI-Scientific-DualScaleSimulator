#!/usr/bin/env python
"""
B2 embedding check (tier X): F(th) = ((cos th_i)_i, (sin th_i sin th_j)_{i<=j}) in R^14.
(1) invariance: |F(th) - F(-th)| over random th (should be 0 up to float);
(2) injectivity on T^4/+-1: over random pairs and over near-collisions found by nearest neighbours
    in F-space, report min |F(a)-F(b)| / d_Q(a,b) and the largest d_Q among pairs with small F-distance.
Seed 20260919. Writes B2_injectivity.json.
Run: prlimit --as=8589934592 -- <venv-python> B2_injectivity_check.py
"""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.spatial import cKDTree
from lib_ph import orbifold_embed, quotient_dist

rng = np.random.default_rng(20260919)
M = 400000
th = rng.uniform(0, 2 * np.pi, size=(M, 4))
F = orbifold_embed(th)
inv = float(np.abs(F - orbifold_embed(-th)).max())
# random pairs
a = rng.integers(M, size=2_000_000); b = rng.integers(M, size=2_000_000)
dq = quotient_dist(th[a], th[b]); dF = np.linalg.norm(F[a] - F[b], axis=1)
m = dq > 1e-9
ratio = dF[m] / dq[m]
# near-collisions: 5 nearest neighbours in F-space
tree = cKDTree(F)
dd, ii = tree.query(F, k=6)
dFn = dd[:, 1:].ravel(); j = ii[:, 1:].ravel(); i = np.repeat(np.arange(M), 5)
dqn = quotient_dist(th[i], th[j])
out = {"tier": "X", "seed": 20260919, "n_points": M, "embedding": "R^14: cos th_i (4), sin th_i sin th_j i<=j (10)",
       "invariance_max_abs": inv,
       "random_pairs": {"n": int(m.sum()), "min_dF_over_dQ": float(ratio.min()),
                        "quantiles_dF_over_dQ": [float(x) for x in np.quantile(ratio, [0.001, 0.01, 0.5])]},
       "nn_pairs": {"n": int(len(dFn)), "max_dQ_among_pairs_with_dF_below_0.05": float(dqn[dFn < 0.05].max()) if (dFn < 0.05).any() else None,
                    "max_dQ_over_all_5NN_pairs": float(dqn.max()), "max_dF_5NN": float(dFn.max())},
       "note": "injective on the quotient by construction argument (cos fixes |th_i|, sin_i sin_j fix relative signs); "
               "measured: min dF/dQ over 2e6 random pairs is reported above (bounded away from 0 on this sample); near a singular point "
               "the sin_i sin_j block is quadratic in th, so for pairs BOTH within r of that point dF/dQ can shrink like r "
               "(not probed below r=0.3 here)."}
# ratio near singular points vs generic
sing = np.array([[x * np.pi for x in bits] for bits in np.ndindex(2, 2, 2, 2)])
out["min_dF_over_dQ_pairs_within_0.3_of_a_singular_point"] = None
near = np.min(np.stack([quotient_dist(th[a], s[None, :]) for s in sing]), axis=0) < 0.3
if (near & m).any():
    out["min_dF_over_dQ_pairs_within_0.3_of_a_singular_point"] = float((dF[near & m] / dq[near & m]).min())
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "B2_injectivity.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
