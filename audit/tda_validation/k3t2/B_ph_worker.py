#!/usr/bin/env python
"""
Part B worker (tier X): one persistence computation, spec given as a JSON string argument.
  spec keys: space in {T2,T3,T4,orbifold,K3,K3xT2,null4}, N, seed, method in {rips,alpha,sparse},
             tau (max edge length / max alpha radius), sparse (eps for sparse Rips),
             expected (list), fields (list), fps (bool: farthest-point subsample of N from a pool of 10N),
             save_diagram_dims (list), null_D, null_diam
Prints one JSON line. Called by B_driver.py under prlimit + timeout.
"""
import json, os, resource, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import gudhi
from lib_ph import (sample_fermat_quadric, sample_torus, orbifold_embed, sample_fermat_k3, projector_embed, sample_null_cube,
                    farthest_point_subsample, recovery, bars_from_st)

spec = json.loads(sys.argv[1])
rng = np.random.default_rng(spec["seed"])
N = spec["N"]
t0 = time.time()
info = {}
pool = 10 * N if spec.get("fps") else N


def gen(M):
    sp = spec["space"]
    if sp in ("T2", "T3", "T4"):
        return sample_torus(int(sp[1]), M, rng)[0]
    if sp == "orbifold":
        th = rng.uniform(0, 2 * np.pi, size=(M, 4))
        return orbifold_embed(th)
    if sp in ("K3", "K3xT2"):
        Z, resid = sample_fermat_k3((M + 3) // 4, rng)
        info["fermat_residual_max"] = resid
        X = projector_embed(Z)[:M]
        if sp == "K3xT2":
            T, _ = sample_torus(2, M, rng)
            X = np.concatenate([X, spec.get("torus_scale", 0.5) * T], axis=1)
        return X
    if sp == "quadric":
        Z, resid = sample_fermat_quadric((M + 1) // 2, rng)
        info["fermat_residual_max"] = resid
        return projector_embed(Z)[:M]
    if sp == "null4":
        return sample_null_cube(M, spec["null_D"], spec["null_diam"], rng)
    if sp == "null4m":
        # density-matched null: uniform 4-cube in R^16, rescaled so that its median 10-NN distance equals
        # that of a Fermat-K3 projector sample of the same size (independent draw from the same rng stream)
        from scipy.spatial import cKDTree
        Zk, _ = sample_fermat_k3((M + 3) // 4, rng)
        K = projector_embed(Zk)[:M]
        Y = sample_null_cube(M, 16, 2.0, rng)
        mk = float(np.median(cKDTree(K).query(K, k=11)[0][:, 10]))
        my = float(np.median(cKDTree(Y).query(Y, k=11)[0][:, 10]))
        info["null_scale_factor"] = mk / my
        info["median_10NN_K3"] = mk
        return Y * (mk / my)
    raise ValueError(sp)


X = gen(pool)
if spec.get("fps"):
    X, cov = farthest_point_subsample(X, N, rng)
    info["fps_covering_radius_of_pool"] = cov
info["ambient_dim"] = int(X.shape[1])
info["diameter_est"] = float(np.linalg.norm(X - X[0], axis=1).max() * 1.0)
tau = spec["tau"]
maxk = len(spec["expected"]) - 1
tb = time.time()
if spec["method"] == "alpha":
    st = gudhi.AlphaComplex(points=X).create_simplex_tree(max_alpha_square=tau ** 2)
    # convert filtration to radius
    for s, f in list(st.get_simplices()):
        st.assign_filtration(s, float(np.sqrt(max(f, 0.0))))
    st.make_filtration_non_decreasing()
elif spec["method"] == "rips":
    st = gudhi.RipsComplex(points=X, max_edge_length=tau).create_simplex_tree(max_dimension=1)
    info["n_edges_before_collapse"] = st.num_simplices() - st.num_vertices()
    st.collapse_edges(nb_iterations=spec.get("collapse_iter", 3))
    info["n_edges_after_collapse"] = st.num_simplices() - st.num_vertices()
    st.expansion(maxk + 1)
elif spec["method"] == "witness":
    # strong witness complex: landmarks = farthest-point subsample of L points of the N-point sample (witnesses = all N).
    # GUDHI filtration is a squared relaxation; converted to its square root (length units) as for alpha.
    nL = spec.get("landmarks") or N // spec["witness_ratio"]
    Lm, lcov = farthest_point_subsample(X, nL, rng)
    info["n_landmarks"] = nL; info["landmark_covering_radius"] = lcov
    st = gudhi.EuclideanStrongWitnessComplex(witnesses=X, landmarks=Lm).create_simplex_tree(
        max_alpha_square=tau ** 2, limit_dimension=maxk + 1)
    for s, f in list(st.get_simplices()):
        st.assign_filtration(s, float(np.sqrt(max(f, 0.0))))
    st.make_filtration_non_decreasing()
elif spec["method"] == "sparse":
    st = gudhi.RipsComplex(points=X, max_edge_length=tau, sparse=spec["sparse"]).create_simplex_tree(max_dimension=maxk + 1)
else:
    raise ValueError(spec["method"])
info["n_simplices"] = st.num_simplices()
info["build_seconds"] = round(time.time() - tb, 2)
res = {"spec": spec, "info": info, "by_field": {}}
for p in spec.get("fields", [3]):
    tp = time.time()
    st.compute_persistence(homology_coeff_field=p, persistence_dim_max=False)
    bars = bars_from_st(st, maxk)
    rec = recovery(bars, spec["expected"], tau)
    rec["seconds"] = round(time.time() - tp, 2)
    rec["n_bars"] = {k: len(v) for k, v in bars.items()}
    from lib_ph import betti_curve
    grid = [round(tau * i / 20, 4) for i in range(1, 20)]
    rec["beta_curve"] = {str(e): [int(betti_curve(bars[k], np.array([e]))[0]) for k in range(maxk + 1)] for e in grid}
    # top 12 bars per dim by persistence (for the report)
    rec["top_bars"] = {k: sorted(v, key=lambda x: -(min(x[1], tau) - x[0]))[:12] for k, v in bars.items()}
    for k in spec.get("save_diagram_dims", []):
        rec[f"diagram_H{k}"] = [(b, d if np.isfinite(d) else None) for b, d in bars[k]]
    res["by_field"][f"F{p}"] = rec
res["seconds_total"] = round(time.time() - t0, 2)
res["maxrss_MB"] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1)
print(json.dumps(res, default=lambda o: None if (isinstance(o, float) and not np.isfinite(o)) else str(o)))
