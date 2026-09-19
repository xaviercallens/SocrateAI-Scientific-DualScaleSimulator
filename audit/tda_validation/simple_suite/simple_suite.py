#!/usr/bin/env python3
"""
Known-answer TDA suite on simple spaces (LeanFlow TEST_USE_CASES section 2).

Expected answers, sources, seeds and pass rules are in expectations.json
(committed alone, before this script was run). This script imports the
pipeline functions under test (does not reimplement them) and runs ONE case
per process:

    python simple_suite.py --case P1
    python simple_suite.py --case P7 --n 500
    python simple_suite.py --case F3chunk --start 0 --count 50
    python simple_suite.py --aggregate          # -> results.json

Each case writes cases/<case>.json with its observed bars, pass flags,
function used, wall time, peak RSS (ru_maxrss of this process), load average
and the exact command. Tier of every observed number: X (numerics).
See run_suite.sh for the wrapper (timeout + prlimit --as=6442450944).
"""
import argparse
import hashlib
import importlib.util
import json
import math
import os
import resource
import subprocess
import sys
import time

import numpy as np
import gudhi
import gudhi.subsampling

HERE = os.path.dirname(os.path.abspath(__file__))
WT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
CASE_DIR = os.path.join(HERE, "cases")
EXPECT = json.load(open(os.path.join(HERE, "expectations.json")))

CW_REL = "audit/reverse_zero/E5-cosmic-web-tda-scaled/cosmic_web_tda_scaled.py"
CMB_REL = "audit/reverse_zero/E5-cmb-tda/cmb_tda.py"


def load_module(name, rel):
    """Import a pipeline file under a non-__main__ name so its main() never runs."""
    spec = importlib.util.spec_from_file_location(name, os.path.join(WT, rel))
    mod = importlib.util.module_from_spec(spec)
    assert mod.__name__ != "__main__"
    spec.loader.exec_module(mod)
    return mod


def sha256(rel):
    return hashlib.sha256(open(os.path.join(WT, rel), "rb").read()).hexdigest()


CW = load_module("cosmic_web_tda_scaled_under_test", CW_REL)
CMB = load_module("cmb_tda_under_test", CMB_REL)

L_REF = 1.2          # Rips reference threshold (expectations.json)
L_T3 = 1.6           # P7 threshold


# ----------------------------------------------------------------- helpers
def git_head():
    try:
        return subprocess.check_output(["git", "-C", WT, "rev-parse", "HEAD"], text=True).strip()
    except Exception as e:  # report, do not hide
        return "unavailable: %r" % e


def peak_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def alpha_by_dim_from_pipeline(xyz, max_alpha_sq, label):
    """Pipeline path: cosmic_web_tda_scaled.alpha_persistence (Z/2, dims 0-2, radius units)."""
    res, by_dim = CW.alpha_persistence(xyz, max_alpha_sq, label, HERE)
    bars = {k: np.array(v, dtype=float).reshape(-1, 2) for k, v in by_dim.items()}
    return bars, {"n_points": res["n_points"], "pipeline_runtime_sec": res["runtime_sec"],
                  "betti_numbers_at_truncation": res["betti_numbers_at_truncation"]}


def intervals(st, maxdim, sqrt=False):
    out = {}
    for k in range(maxdim + 1):
        iv = np.array(st.persistence_intervals_in_dimension(k), dtype=float).reshape(-1, 2)
        if sqrt:
            iv = np.where(np.isfinite(iv), np.sqrt(np.maximum(iv, 0.0)), np.inf)
        out[k] = iv
    return out


def rips_tree(X, L, maxdim, collapse=True):
    st = gudhi.RipsComplex(points=X, max_edge_length=L).create_simplex_tree(max_dimension=1)
    n_graph = st.num_simplices()
    n_iter = 0
    if collapse:
        while True:
            n0 = st.num_simplices()
            st.collapse_edges()
            n_iter += 1
            if st.num_simplices() == n0 or n_iter > 100:
                break
    n_coll = st.num_simplices()
    st.expansion(maxdim + 1)
    return st, {"L": L, "n_points": int(X.shape[0]), "graph_simplices": int(n_graph),
                "after_collapse_simplices": int(n_coll), "collapse_iterations": n_iter,
                "expanded_simplices": int(st.num_simplices()), "expansion_dim": maxdim + 1}


def rips_diagram(X, L, maxdim, fields=(2,), collapse=True):
    st, info = rips_tree(X, L, maxdim, collapse)
    out = {}
    for f in fields:
        st.compute_persistence(homology_coeff_field=f)
        out[f] = intervals(st, maxdim)
    return out, info


def evaluate(bars, expected, cap=None, gate_dims=None, top=6):
    """Apply the pre-registered dominant/zero rules (expectations.json 'definitions')."""
    dims = sorted(bars.keys())
    if cap is None:
        fin = [bars[k][np.isfinite(bars[k][:, 1]), 1] for k in dims if bars[k].size]
        fin = np.concatenate(fin) if fin else np.array([0.0])
        cap = float(fin.max()) if fin.size else 0.0
    pers = {}
    capped = {}
    for k in dims:
        b = bars[k]
        if b.size == 0:
            pers[k] = np.zeros(0)
            capped[k] = np.zeros((0, 2))
            continue
        d = np.where(np.isfinite(b[:, 1]), np.minimum(b[:, 1], cap), cap)
        p = d - b[:, 0]
        o = np.argsort(p)[::-1]
        pers[k] = p[o]
        capped[k] = np.c_[b[o, 0], d[o]]
    p_all = max([pers[k][0] for k in dims if k >= 1 and pers[k].size] + [0.0])
    # expectations.json zero_rule clause: if the case expects no bar in any dim >= 1
    # (e.g. RP2 over Z/3), the no-dominant-bar rule p_1 < 5 p_2 applies in dims >= 1 and
    # H0 is judged against its own capped essential bar.
    acyclic = all(e == 0 for e in expected[1:])
    if gate_dims is None:
        gate_dims = list(range(len(expected)))
    per_dim = {}
    ok = True
    for k in dims:
        p = pers[k]
        pj = lambda j: float(p[j - 1]) if j <= p.size else 0.0
        exp_k = expected[k] if k < len(expected) else 0
        gated = k in gate_dims
        if acyclic and k >= 1:
            ratio = pj(1) / pj(2) if pj(2) > 0 else (float("inf") if pj(1) > 0 else 0.0)
            rule = ratio < 5.0
            obs = 0 if rule else 1
            rule_txt = "acyclic clause: p_1 / p_2 = %.4g < 5" % ratio
            dim_pass = bool(rule and obs == exp_k)
            if gated:
                ok = ok and dim_pass
            per_dim[str(k)] = {"expected": exp_k, "observed_beta": obs, "rule": rule_txt, "rule_pass": bool(rule),
                               "pass": dim_pass, "gated": gated, "n_bars": int(p.size),
                               "top_persistence": [float(x) for x in p[:top]],
                               "top_bars_birth_death_capped": [[float(a), float(c)] for a, c in capped[k][:top]],
                               "n_infinite": int(np.sum(~np.isfinite(bars[k][:, 1]))) if bars[k].size else 0}
            continue
        ref = pj(1) if (acyclic and k == 0) else p_all
        obs = int(np.sum(p >= 0.2 * ref)) if ref > 0 else int(p.size)
        if exp_k >= 1:
            ratio = pj(exp_k) / pj(exp_k + 1) if pj(exp_k + 1) > 0 else float("inf")
            rule = ratio >= 5.0
            rule_txt = "p_%d / p_%d = %.4g >= 5" % (exp_k, exp_k + 1, ratio)
        else:
            ratio = pj(1) / p_all if p_all > 0 else 0.0
            rule = pj(1) < 0.2 * p_all
            rule_txt = "p_1 / P_all = %.4g < 0.2" % ratio
        dim_pass = bool(rule and obs == exp_k)
        if gated:
            ok = ok and dim_pass
        per_dim[str(k)] = {
            "expected": exp_k, "observed_beta": obs, "rule": rule_txt, "rule_pass": bool(rule),
            "pass": dim_pass, "gated": gated, "n_bars": int(p.size),
            "top_persistence": [float(x) for x in p[:top]],
            "top_bars_birth_death_capped": [[float(a), float(c)] for a, c in capped[k][:top]],
            "n_infinite": int(np.sum(~np.isfinite(bars[k][:, 1]))) if bars[k].size else 0,
        }
    return {"cap": cap, "P_all": float(p_all), "per_dim": per_dim, "pass": bool(ok),
            "observed_betti_vector": [per_dim[str(k)]["observed_beta"] for k in dims]}


def ratio12(bars, k):
    b = bars.get(k, np.zeros((0, 2)))
    if b.size == 0:
        return {"p1": 0.0, "p2": 0.0, "ratio": 0.0}
    p = np.sort(b[:, 1] - b[:, 0])[::-1]
    p1 = float(p[0]); p2 = float(p[1]) if p.size > 1 else 0.0
    return {"p1": p1, "p2": p2, "ratio": (p1 / p2) if p2 > 0 else float("inf")}


def sub(X, n, seed):
    rng = np.random.default_rng(seed)
    return X[rng.choice(X.shape[0], size=min(n, X.shape[0]), replace=False)]


# ----------------------------------------------------------------- clouds
def sphere_pts(n, rng):
    v = rng.normal(size=(n, 3))
    return v / np.linalg.norm(v, axis=1, keepdims=True)


def cloud(case):
    if case == "P1":
        return CW.build_circle_cloud(n=2000, noise_sigma=0.05, seed=11)
    if case == "P2":
        rng = np.random.default_rng(12)
        return sphere_pts(5000, rng) + rng.normal(0, 0.02, (5000, 3))
    if case == "P3a":
        rng = np.random.default_rng(13)
        R, r, N = 2.5, 1.0, 8000
        u_list, v_list = [], []
        while sum(len(x) for x in u_list) < N:
            u = rng.uniform(0, 2 * np.pi, 4 * N); v = rng.uniform(0, 2 * np.pi, 4 * N)
            keep = rng.uniform(0, 1, 4 * N) < (R + r * np.cos(v)) / (R + r)
            u_list.append(u[keep]); v_list.append(v[keep])
        u = np.concatenate(u_list)[:N]; v = np.concatenate(v_list)[:N]
        X = np.c_[(R + r * np.cos(v)) * np.cos(u), (R + r * np.cos(v)) * np.sin(u), r * np.sin(v)]
        return X + rng.normal(0, 0.02, X.shape)
    if case == "P3b":
        rng = np.random.default_rng(14)
        a = rng.uniform(0, 2 * np.pi, 8000); b = rng.uniform(0, 2 * np.pi, 8000)
        return np.c_[np.cos(a), np.sin(a), np.cos(b), np.sin(b)]
    if case == "P4":
        rng = np.random.default_rng(15)
        N = 2000
        side = np.where(rng.uniform(size=N) < 0.5, -1.0, 1.0)
        t = rng.uniform(0, 2 * np.pi, N)
        X = np.c_[side + np.cos(t), np.sin(t), np.zeros(N)]
        return X + rng.normal(0, 0.05, X.shape)
    if case == "P5":
        rng = np.random.default_rng(16)
        R, r, N = 2.0, 1.0, 8000
        u = rng.uniform(0, 2 * np.pi, N); v = rng.uniform(0, 2 * np.pi, N)
        return np.c_[(R + r * np.cos(v)) * np.cos(u), (R + r * np.cos(v)) * np.sin(u),
                     r * np.sin(v) * np.cos(u / 2), r * np.sin(v) * np.sin(u / 2)]
    if case in ("P6", "P6v"):
        rng = np.random.default_rng(17)
        s = sphere_pts(6000, rng)
        if case == "P6" and ARGS.n:
            # POST-HOC deviation (see report.json): N=6000 exceeded the 6 GiB cap in the
            # pipeline alpha_persistence (std::bad_alloc); use the first n of the same iid draw.
            s = s[:ARGS.n]
        x, y, z = s[:, 0], s[:, 1], s[:, 2]
        if case == "P6":
            return np.c_[x * y, x * z, y ** 2 - z ** 2, 2 * y * z]
        q = math.sqrt(2.0)
        return np.c_[x ** 2, y ** 2, z ** 2, q * x * y, q * x * z, q * y * z]
    if case == "P10a":
        rng = np.random.default_rng(18)
        X = np.r_[1.0 * sphere_pts(1500, rng), 2.5 * sphere_pts(4500, rng)]
        return X + rng.normal(0, 0.02, X.shape)
    if case == "P10b":
        rng = np.random.default_rng(19)
        S = sphere_pts(3000, rng)
        t = rng.uniform(0, 2 * np.pi, 1000)
        C = np.c_[4.0 + np.cos(t), np.sin(t), np.zeros(1000)]
        X = np.r_[S, C]
        return X + rng.normal(0, 0.02, X.shape)
    if case == "N1":
        return np.random.default_rng(31).normal(size=(2000, 3))
    if case == "N2":
        X = cloud("P3a").copy()
        for c in range(3):
            X[:, c] = np.random.default_rng(32 + c).permutation(X[:, c])
        return X
    raise KeyError(case)


EXPECTED = {"P1": [1, 1, 0], "P2": [1, 0, 1], "P3a": [1, 2, 1], "P3b": [1, 2, 1], "P4": [1, 2, 0],
            "P10a": [2, 0, 2], "P10b": [2, 1, 1]}
REF_N = {"P1": (400, 101), "P2": (400, 102), "P3a": (1000, 103), "P3b": (1000, 104), "P4": (400, 105),
         "P10a": (600, 108), "P10b": (600, 109)}


# ----------------------------------------------------------------- cases
def case_manifold(case):
    X = cloud(case)
    exp = EXPECTED[case]
    t0 = time.time()
    bars, info = alpha_by_dim_from_pipeline(X, float("inf"), case)
    t_alpha = time.time() - t0
    ev = evaluate(bars, exp)
    # pipeline top_bars on H1 (with the cap used above) for the record
    by_dim_lists = {k: bars[k].tolist() for k in bars}
    tb = CW.top_bars(by_dim_lists, 1, k=4, r_trunc=ev["cap"]).tolist()
    t1 = time.time()
    n, seed = REF_N[case]
    rd, rinfo = rips_diagram(sub(X, n, seed), L_REF, 2, fields=(2,))
    ev_r = evaluate(rd[2], exp, cap=L_REF)
    t_rips = time.time() - t1
    bbox = [X.min(axis=0).tolist(), X.max(axis=0).tolist()]
    return {
        "expected_betti": exp,
        "primary": {"function": "cosmic_web_tda_scaled.alpha_persistence(max_alpha_sq=inf) [pipeline, Z/2]",
                    "info": info, "wall_sec": t_alpha, "evaluation": ev,
                    "pipeline_top_bars_H1": tb},
        "reference": {"function": "gudhi.RipsComplex + collapse_edges + expansion(3), Z/2", "info": rinfo,
                      "subsample_seed": seed, "wall_sec": t_rips, "evaluation": ev_r},
        "pass": bool(ev["pass"]),
        "reference_pass": bool(ev_r["pass"]),
        "bbox": bbox,
    }


def case_coeff(case):
    """P5 (Klein bottle) and P6 (RP2): pipeline (Z/2), direct alpha (Z/2, Z/3), Rips reference (Z/2, Z/3)."""
    exp = {"P5": {2: [1, 2, 1], 3: [1, 1, 0]}, "P6": {2: [1, 1, 1], 3: [1, 0, 0]}}[case]
    X = cloud(case)
    out = {"expected_betti": {"Z/2": exp[2], "Z/3": exp[3]}}
    part = ARGS.part
    if part in ("pipeline", "all"):
        t0 = time.time()
        bars, info = alpha_by_dim_from_pipeline(X, float("inf"), case)
        ev = evaluate(bars, exp[2])
        out["pipeline_Z2"] = {"function": "cosmic_web_tda_scaled.alpha_persistence (Z/2 hard-coded)",
                              "info": info, "wall_sec": time.time() - t0, "evaluation": ev}
    if part in ("direct", "all"):
        t0 = time.time()
        ac = gudhi.AlphaComplex(points=X)
        st = ac.create_simplex_tree()
        del ac  # free the CGAL triangulation before persistence (the pipeline run peaked at 5.3 GB)
        nsimp = st.num_simplices()
        tb = time.time() - t0
        d = {}
        for f in (2, 3):
            t1 = time.time()
            st.compute_persistence(homology_coeff_field=f)
            iv = intervals(st, 3, sqrt=True)
            ev = evaluate(iv, exp[f], gate_dims=[0, 1, 2])
            d["Z/%d" % f] = {"evaluation": ev, "persistence_wall_sec": time.time() - t1,
                             "betti_numbers_end_of_filtration": st.betti_numbers()}
        out["direct_alpha"] = {"function": "gudhi.AlphaComplex -> one SimplexTree -> compute_persistence(field 2, then 3)",
                               "n_simplices": int(nsimp), "build_sec": tb, "fields": d}
    if part in ("rips", "all"):
        t0 = time.time()
        Xr = cloud("P6v") if case == "P6" else X
        n, seed = {"P5": (1000, 106), "P6": (600, 107)}[case]
        rd, rinfo = rips_diagram(sub(Xr, n, seed), L_REF, 2, fields=(2, 3))
        out["rips_reference"] = {"function": "gudhi.RipsComplex + collapse + expansion(3), fields 2 and 3 on one tree",
                                 "embedding": "R6 Veronese" if case == "P6" else "R4 (same as alpha)",
                                 "info": rinfo, "subsample_seed": seed, "wall_sec": time.time() - t0,
                                 "fields": {"Z/%d" % f: {"evaluation": evaluate(rd[f], exp[f], cap=L_REF)} for f in (2, 3)}}
    return out


def case_p7(n):
    rng = np.random.default_rng(7)
    A = rng.uniform(0, 2 * np.pi, (n, 3))
    X = np.c_[np.cos(A[:, 0]), np.sin(A[:, 0]), np.cos(A[:, 1]), np.sin(A[:, 1]), np.cos(A[:, 2]), np.sin(A[:, 2])]
    t0 = time.time()
    rd, info = rips_diagram(X, L_T3, 3, fields=(2,))
    ev = evaluate(rd[2], [1, 3, 3, 1], cap=L_T3)
    return {"N": n, "expected_betti": [1, 3, 3, 1],
            "function": "gudhi.RipsComplex(L=1.6) + collapse_edges + expansion(4), Z/2 (pipeline alpha_persistence cannot report H3)",
            "info": info, "wall_sec": time.time() - t0, "evaluation": ev, "pass": bool(ev["pass"])}


def case_p7lm(n_land, n_pool=48000):
    """POST-HOC (not pre-registered): farthest-point (maxmin) landmarks from a larger uniform
    pool, to test whether even spacing recovers (1,3,3,1) within budget where random N did not."""
    rng = np.random.default_rng(77)
    A = rng.uniform(0, 2 * np.pi, (n_pool, 3))
    X = np.c_[np.cos(A[:, 0]), np.sin(A[:, 0]), np.cos(A[:, 1]), np.sin(A[:, 1]), np.cos(A[:, 2]), np.sin(A[:, 2])]
    t0 = time.time()
    L = gudhi.subsampling.choose_n_farthest_points(points=X, nb_points=n_land, starting_point=0)
    L = np.asarray(L)
    t_land = time.time() - t0
    rd, info = rips_diagram(L, L_T3, 3, fields=(2,))
    ev = evaluate(rd[2], [1, 3, 3, 1], cap=L_T3)
    return {"post_hoc": True, "n_landmarks": n_land, "pool": n_pool, "pool_seed": 77, "landmark_sec": t_land,
            "function": "gudhi.subsampling.choose_n_farthest_points + RipsComplex(L=1.6) + collapse + expansion(4), Z/2",
            "info": info, "wall_sec": time.time() - t0, "evaluation": ev, "pass": bool(ev["pass"])}


def case_p8():
    out = {}
    for seed in (42, 0, 1, 2, 3, 4):
        t0 = time.time()
        xyz, centers = CW.build_planted_void_cloud(n_fill=15000, box=100.0, void_radius=15.0, n_voids=6, seed=seed)
        res, bd = CW.alpha_persistence(xyz, 37.5 ** 2, "void_seed%d" % seed, HERE)
        tb = CW.top_bars(bd, 2, k=8, r_trunc=37.5)
        deaths = tb[:, 1]; pers = tb[:, 1] - tb[:, 0]
        within = np.abs(deaths - 15.0) / 15.0 < 0.05
        rule = bool(len(centers) == 6 and within[:6].all() and (len(within) < 7 or not within[6]))
        within20 = int(np.sum(np.abs(deaths - 15.0) / 15.0 < 0.20))
        gap = float(pers[5] - pers[6]) if len(pers) > 6 else None
        out[str(seed)] = {
            "n_points": int(xyz.shape[0]), "n_voids_planted": int(len(centers)),
            "top8_H2_bars_birth_death": tb.tolist(), "persistence": pers.tolist(),
            "death_rel_err_top6": (np.abs(deaths[:6] - 15.0) / 15.0).tolist(),
            "tda5_rule_pass": rule, "ratio_p6_over_p7": float(pers[5] / pers[6]) if len(pers) > 6 else None,
            "pipeline_gate_20pct_count": within20, "pipeline_gate_gap": gap,
            "pipeline_gate_pass": bool(within20 == len(centers) and (gap is None or gap > 1.0)),
            "wall_sec": time.time() - t0,
        }
    return {"function": "cosmic_web_tda_scaled.build_planted_void_cloud + alpha_persistence(37.5^2) + top_bars(2, k=8, r_trunc=37.5)",
            "seeds": out, "pass": out["42"]["tda5_rule_pass"],
            "robustness_pass_count": int(sum(v["tda5_rule_pass"] for k, v in out.items() if k != "42"))}


def case_p9(which):
    src, n, s0 = {"P9a": ("P1", 2000, 5000), "P9b": ("P2", 5000, 6000), "P9c": ("P3a", 8000, 7000)}[which]
    X = cloud(src)
    lo, hi = X.min(axis=0), X.max(axis=0)
    rows = []
    for s in range(s0, s0 + 50):
        P = lo + (hi - lo) * np.random.default_rng(s).random((n, 3))
        bars, _ = alpha_by_dim_from_pipeline(P, float("inf"), "poisson%d" % s)
        r1, r2 = ratio12(bars, 1), ratio12(bars, 2)
        rows.append({"seed": s, "H1": r1, "H2": r2, "no_dominant": bool(r1["ratio"] < 5 and r2["ratio"] < 5)})
    m1 = np.array([r["H1"]["p1"] for r in rows]); m2 = np.array([r["H2"]["p1"] for r in rows])
    q = lambda a: {"min": float(a.min()), "p5": float(np.percentile(a, 5)), "median": float(np.median(a)),
                   "p95": float(np.percentile(a, 95)), "max": float(a.max())}
    nd = int(sum(r["no_dominant"] for r in rows))
    return {"matched_to": src, "N": n, "box_lo": lo.tolist(), "box_hi": hi.tolist(), "seeds": [s0, s0 + 49],
            "function": "cosmic_web_tda_scaled.alpha_persistence(max_alpha_sq=inf)",
            "max_persistence_H1": q(m1), "max_persistence_H2": q(m2),
            "ratio_p1_p2_H1": q(np.array([r["H1"]["ratio"] for r in rows])),
            "ratio_p1_p2_H2": q(np.array([r["H2"]["ratio"] for r in rows])),
            "n_seeds_no_dominant_bar": nd, "pass": bool(nd >= 48),
            "null_max_H1": m1.tolist(), "null_max_H2": m2.tolist(), "per_seed": rows}


def case_neg(case):
    X = cloud(case)
    t0 = time.time()
    bars, info = alpha_by_dim_from_pipeline(X, float("inf"), case)
    r1, r2 = ratio12(bars, 1), ratio12(bars, 2)
    p1 = np.sort(bars[1][:, 1] - bars[1][:, 0])[::-1]
    out = {"function": "cosmic_web_tda_scaled.alpha_persistence(max_alpha_sq=inf)", "info": info,
           "wall_sec": time.time() - t0, "H1": r1, "H2": r2, "H1_top5": p1[:5].tolist(),
           "H2_top3": np.sort(bars[2][:, 1] - bars[2][:, 0])[::-1][:3].tolist()}
    if case == "N1":
        out["pass"] = bool(r1["ratio"] < 5 and r2["ratio"] < 5)
    else:
        p2, p3 = float(p1[1]), float(p1[2])
        out["H1_p2_over_p3"] = p2 / p3
        out["H1_second_bar_persistence"] = p2
        out["pass_ratio_part"] = bool(p2 / p3 < 5)  # the p-value part is applied in --aggregate
    return out


def case_c0():
    X = sub(cloud("P3a"), 300, 40)
    a, ia = rips_diagram(X, L_REF, 2, collapse=False)
    b, ib = rips_diagram(X, L_REF, 2, collapse=True)
    bd = {}
    for k in (0, 1, 2):
        A = np.where(np.isfinite(a[2][k]), a[2][k], L_REF); B = np.where(np.isfinite(b[2][k]), b[2][k], L_REF)
        bd[str(k)] = float(gudhi.bottleneck_distance(A, B))
    return {"function": "gudhi.RipsComplex with vs without SimplexTree.collapse_edges", "no_collapse": ia,
            "collapse": ib, "bottleneck_by_dim": bd, "pass": bool(max(bd.values()) <= 1e-12)}


F1_CENTRES = None


def case_f1_field():
    global F1_CENTRES
    n, sig = 256, 10.0
    rng = np.random.default_rng(21)
    C = []
    while len(C) < 7:
        c = rng.uniform(30, 226, 2)
        if all(np.linalg.norm(c - d) >= 50 for d in C):
            C.append(c)
    C = np.array(C)
    F1_CENTRES = C
    ii, jj = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    f = np.zeros((n, n))
    for c in C:
        f -= np.exp(-((ii - c[0]) ** 2 + (jj - c[1]) ** 2) / (2 * sig ** 2))
    return f


def lower_star_counts(f):
    n = f.shape[0]
    t0 = time.time()
    cc = gudhi.CubicalComplex(top_dimensional_cells=f)
    cc.compute_persistence()
    h0 = np.array(cc.persistence_intervals_in_dimension(0)).reshape(-1, 2)
    n_cub = int(np.sum((h0[:, 1] - h0[:, 0] > 0) | ~np.isfinite(h0[:, 1])))
    t_cub = time.time() - t0
    # direct strict 8-neighbour local minima
    pad = np.pad(f, 1, constant_values=np.inf)
    ismin = np.ones_like(f, dtype=bool)
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di == 0 and dj == 0:
                continue
            ismin &= f < pad[1 + di:1 + di + n, 1 + dj:1 + dj + n]
    n_direct = int(ismin.sum())
    # pipeline lower-star path on a Freudenthal triangulation of the grid
    t1 = time.time()
    idx = np.arange(n * n).reshape(n, n)
    e = np.r_[np.c_[idx[:, :-1].ravel(), idx[:, 1:].ravel()],
              np.c_[idx[:-1, :].ravel(), idx[1:, :].ravel()],
              np.c_[idx[:-1, :-1].ravel(), idx[1:, 1:].ravel()]]
    tr = np.r_[np.c_[idx[:-1, :-1].ravel(), idx[:-1, 1:].ravel(), idx[1:, 1:].ravel()],
               np.c_[idx[:-1, :-1].ravel(), idx[1:, :-1].ravel(), idx[1:, 1:].ravel()]]
    temp = f.ravel()
    sd = temp.std()
    nu = np.linspace(temp.min() / sd, temp.max() / sd, 4001)
    b0, b1, chi = CMB.betti_curves_from_topology(temp, np.arange(n * n), e.astype(np.int64), tr.astype(np.int64), nu, sublevel=True)
    t_pipe = time.time() - t1
    return {"grid": n,
            "cubical": {"function": "gudhi.CubicalComplex(top_dimensional_cells=f)", "n_H0_bars_positive_persistence": n_cub,
                        "n_H0_bars_total": int(h0.shape[0]), "wall_sec": t_cub},
            "direct_local_minima_8nbr": n_direct,
            "pipeline": {"function": "cmb_tda.betti_curves_from_topology (Freudenthal grid triangulation, sublevel)",
                         "n_edges": int(e.shape[0]), "n_triangles": int(tr.shape[0]), "nu_points": int(nu.size),
                         "max_nu_b0": int(b0.max()), "b0_at_max_nu": int(b0[-1]), "b1_at_max_nu": int(b1[-1]),
                         "max_nu_b1": int(b1.max()), "wall_sec": t_pipe},
            "n_minimum_values_listed": min(n_direct, 20),
            "minimum_values_first20": sorted(f[ismin].tolist())[:20]}


def case_f1():
    f = case_f1_field()
    r = lower_star_counts(f)
    ok = bool(r["cubical"]["n_H0_bars_positive_persistence"] == 7 and r["pipeline"]["max_nu_b0"] == 7
              and r["direct_local_minima_8nbr"] == 7 and r["pipeline"]["b0_at_max_nu"] == 1 and r["pipeline"]["b1_at_max_nu"] == 0)
    return {"centres": F1_CENTRES.tolist(), "sigma_px": 10.0, **r, "pass": ok}


def case_f1shuf():
    """POST-HOC control added after TEST_USE_CASES gained TDA-N3 (2026-09-19, during this run):
    the F1 field with its sites shuffled must NOT give 7 minima."""
    base = case_f1_field()
    n = base.shape[0]
    f = np.random.default_rng(23).permutation(base.ravel()).reshape(n, n)
    r = lower_star_counts(f)
    ok = bool(r["pipeline"]["max_nu_b0"] != 7 and r["cubical"]["n_H0_bars_positive_persistence"] != 7)
    return {**r, "post_hoc": True, "pass": ok,
            "pass_rule": "pipeline max_nu b0 != 7 and CubicalComplex count != 7 on the shuffled field"}


def case_f2():
    from math import comb
    out = {}
    allok = True
    for n, m in ((1, 64), (2, 24), (3, 12), (4, 6)):
        rng = np.random.default_rng(22 + n)
        field = 1.0 + 1e-3 * rng.normal(size=(m,) * n)
        t0 = time.time()
        pc = gudhi.PeriodicCubicalComplex(top_dimensional_cells=field, periodic_dimensions=[True] * n)
        pc.compute_persistence()
        inf_counts = [int(np.sum(~np.isfinite(np.array(pc.persistence_intervals_in_dimension(k)).reshape(-1, 2)[:, 1])))
                      for k in range(n + 1)]
        betti = pc.betti_numbers()
        cc = gudhi.CubicalComplex(top_dimensional_cells=field)
        cc.compute_persistence()
        neg = cc.betti_numbers()
        exp = [comb(n, k) for k in range(n + 1)]
        ok = bool(inf_counts == exp and list(betti)[:n + 1] == exp and neg[0] == 1 and all(x == 0 for x in neg[1:]))
        allok = allok and ok
        out["T%d" % n] = {"grid": [m] * n, "expected": exp, "infinite_bars_per_dim": inf_counts,
                          "betti_numbers": list(betti), "nonperiodic_control_betti": list(neg),
                          "pass": ok, "wall_sec": time.time() - t0}
    return {"function": "gudhi.PeriodicCubicalComplex / gudhi.CubicalComplex", "tori": out, "pass": allok}


def f3_cl():
    lmax = 3 * 64 - 1
    ell = np.arange(lmax + 1)
    sb = math.radians(3.0) / math.sqrt(8 * math.log(2))
    cl = np.zeros(lmax + 1)
    cl[2:] = np.exp(-ell[2:] * (ell[2:] + 1) * sb ** 2) / (ell[2:] * (ell[2:] + 1))
    return cl


def f3_topology():
    import healpy as hp  # noqa
    mask = np.ones(12 * 64 ** 2, dtype=np.uint8)
    return CMB.build_topology(mask, 64)


def case_f3chunk(start, count):
    import healpy as hp
    t0 = time.time()
    unm, e, tr = f3_topology()
    t_topo = time.time() - t0
    cl = f3_cl()
    nu = np.linspace(CMB.NU_MIN, CMB.NU_MAX, CMB.NU_STEP_GRID)
    maps = []
    for i in range(start, start + count):
        np.random.seed(30000 + i)
        m = hp.synfast(cl, nside=64, new=True)
        b0, b1, chi = CMB.betti_curves_from_topology(m, unm, e, tr, nu, sublevel=True)
        maps.append({"i": i, "seed": 30000 + i, "b0": b0.tolist(), "b1": b1.tolist(), "chi": chi.tolist()})
    return {"function": "cmb_tda.build_topology(full-sky mask, 64) + betti_curves_from_topology(sublevel)",
            "start": start, "count": count, "topology_sec": t_topo, "wall_sec": time.time() - t0,
            "n_vertices": int(unm.size), "n_edges": int(e.shape[0]), "n_triangles": int(tr.shape[0]), "maps": maps}


def case_f3topo():
    t0 = time.time()
    unm, e, tr = f3_topology()
    st = gudhi.SimplexTree()
    for i in range(unm.size):
        st.insert([i], 0.0)
    for a, b in e:
        st.insert([int(a), int(b)], 0.0)
    for a, b, c in tr:
        st.insert([int(a), int(b), int(c)], 0.0)
    st.compute_persistence(persistence_dim_max=True)
    betti = st.betti_numbers()
    V, E, T = unm.size, e.shape[0], tr.shape[0]
    return {"function": "cmb_tda.build_topology(full-sky mask, nside 64) -> gudhi.SimplexTree, constant filtration",
            "n_vertices": int(V), "n_edges": int(E), "n_triangles": int(T),
            "euler_char_V_minus_E_plus_T": int(V - E + T), "betti_numbers": list(betti),
            "expected_if_triangulated_S2": [1, 0, 1],
            "pass_b0_b1": bool(betti[0] == 1 and (len(betti) < 2 or betti[1] == 0)),
            "b2_equals_1": bool(len(betti) > 2 and betti[2] == 1), "wall_sec": time.time() - t0}


# ----------------------------------------------------------------- aggregate
def aggregate():
    from scipy.stats import kstest
    cases = {}
    for fn in sorted(os.listdir(CASE_DIR)):
        if fn.endswith(".json"):
            cases[fn[:-5]] = json.load(open(os.path.join(CASE_DIR, fn)))
    res = {"generated_by": "simple_suite.py --aggregate", "git_head_at_aggregate": git_head(),
           "code_under_test_sha256": {CW_REL: sha256(CW_REL), CMB_REL: sha256(CMB_REL)},
           "tier_observed": "X", "tier_expected": "L", "cases": {}}

    def meta(c):
        return {k: c.get(k) for k in ("command", "wall_sec_total", "peak_rss_mb", "loadavg_start", "started")}

    null_of = {"P1": "P9a", "P2": "P9b", "P3a": "P9c", "N2": "P9c"}

    def pval(obs, arr):
        arr = np.asarray(arr)
        return float((1 + np.sum(arr >= obs)) / (1 + arr.size))

    for cid in ("P1", "P2", "P3a", "P3b", "P4", "P10a", "P10b"):
        if cid not in cases:
            continue
        c = cases[cid]
        r = {"expected_betti": c["expected_betti"],
             "observed_betti_pipeline": c["primary"]["evaluation"]["observed_betti_vector"],
             "observed_betti_rips_reference": c["reference"]["evaluation"]["observed_betti_vector"],
             "pass": c["pass"], "reference_pass": c["reference_pass"],
             "function": c["primary"]["function"], "reference_function": c["reference"]["function"],
             "per_dim_pipeline": c["primary"]["evaluation"]["per_dim"],
             "per_dim_reference": c["reference"]["evaluation"]["per_dim"], **meta(c)}
        if cid in null_of and null_of[cid] in cases:
            nl = cases[null_of[cid]]
            pd = c["primary"]["evaluation"]["per_dim"]
            pv = {}
            if c["expected_betti"][1] >= 1:
                b = c["expected_betti"][1]
                pv["H1_bar_%d" % b] = pval(pd["1"]["top_persistence"][b - 1], nl["null_max_H1"])
            if c["expected_betti"][2] >= 1:
                pv["H2_bar_1"] = pval(pd["2"]["top_persistence"][0], nl["null_max_H2"])
            r["null"] = null_of[cid]
            r["p_values"] = pv
            r["p_value_gate_pass"] = bool(all(v <= 0.05 for v in pv.values()))
            r["pass"] = bool(r["pass"] and r["p_value_gate_pass"])
        res["cases"][cid] = r
    for cid in ("P5", "P6"):
        parts = {k: v for k, v in cases.items() if k.startswith(cid + "_")}
        if not parts:
            continue
        variants = {}
        for k, v in parts.items():
            part, _, nsuf = k[len(cid) + 1:].partition("_")
            variants.setdefault(nsuf or "preregistered", {})[part] = v
        r = {"expected_betti": next(iter(parts.values())).get("expected_betti"), "variants": {}}
        for vname, vparts in sorted(variants.items()):
            vr = {"parts": {}}
            flags = []
            for part, v in vparts.items():
                pr = {kk: v.get(kk) for kk in ("status", "error", "command", "wall_sec_total", "peak_rss_mb", "loadavg_start")}
                if v.get("status") == "error":
                    flags.append(False)
                if "pipeline_Z2" in v:
                    ev = v["pipeline_Z2"]["evaluation"]
                    pr["pipeline_Z2"] = {"observed": ev["observed_betti_vector"], "pass": ev["pass"], "per_dim": ev["per_dim"],
                                         "n_points": v["pipeline_Z2"]["info"]["n_points"], "function": v["pipeline_Z2"]["function"]}
                    flags.append(ev["pass"])
                if "direct_alpha" in v:
                    pr["direct_alpha_n_simplices"] = v["direct_alpha"]["n_simplices"]
                    pr["function"] = v["direct_alpha"]["function"]
                    for f in ("Z/2", "Z/3"):
                        ev = v["direct_alpha"]["fields"][f]["evaluation"]
                        pr["direct_alpha_" + f] = {"observed": ev["observed_betti_vector"], "pass": ev["pass"], "per_dim": ev["per_dim"]}
                        flags.append(ev["pass"])
                if "rips_reference" in v:
                    pr["rips_info"] = v["rips_reference"]["info"]
                    pr["function"] = v["rips_reference"]["function"]
                    pr["embedding"] = v["rips_reference"]["embedding"]
                    for f in ("Z/2", "Z/3"):
                        ev = v["rips_reference"]["fields"][f]["evaluation"]
                        pr["rips_" + f] = {"observed": ev["observed_betti_vector"], "pass": ev["pass"], "per_dim": ev["per_dim"]}
                vr["parts"][part] = pr
            has_alpha = any(("pipeline_Z2" in v or "direct_alpha" in v or v.get("status") == "error") for v in vparts.values())
            vr["pass"] = bool(has_alpha and flags and all(flags))
            vr["pass_rule"] = "every alpha run in this variant (pipeline Z/2, direct Z/2 and Z/3) passes; an error counts as FAIL; Rips reference reported, not gated"
            r["variants"][vname] = vr
        r["pass"] = r["variants"].get("preregistered", {}).get("pass", False)
        r["pass_note"] = "pass = pre-registered design only; post-hoc variants (other N) are reported under variants and never change it"
        res["cases"][cid] = r
    p7 = {k: v for k, v in cases.items() if k.startswith("P7_N")}
    if p7:
        rows = sorted(p7.values(), key=lambda v: v["N"])
        passed = [v["N"] for v in rows if v.get("pass")]
        res["cases"]["P7"] = {"expected_betti": [1, 3, 3, 1],
                              "sweep": [{"N": v["N"], "status": v.get("status", "ok"), "observed": v.get("evaluation", {}).get("observed_betti_vector"),
                                         "pass": v.get("pass"), "per_dim": v.get("evaluation", {}).get("per_dim"),
                                         "info": v.get("info"), "command": v.get("command"),
                                         "wall_sec_total": v.get("wall_sec_total"), "peak_rss_mb": v.get("peak_rss_mb")} for v in rows],
                              "smallest_N_recovered": min(passed) if passed else None,
                              "function": rows[0].get("function"), "pass": bool(passed)}
    if "P8" in cases:
        c = cases["P8"]
        res["cases"]["P8"] = {"expected": "6 H2 bars with death within 5% of 15; bar 7 not", "pass": c["pass"],
                              "seed42": c["seeds"]["42"], "robustness_pass_count_of_5": c["robustness_pass_count"],
                              "robustness": {k: {kk: v[kk] for kk in ("tda5_rule_pass", "death_rel_err_top6", "ratio_p6_over_p7", "pipeline_gate_pass")}
                                             for k, v in c["seeds"].items() if k != "42"},
                              "function": c["function"], **meta(c)}
    for cid in ("P9a", "P9b", "P9c"):
        if cid in cases:
            c = cases[cid]
            res["cases"][cid] = {k: c[k] for k in ("matched_to", "N", "box_lo", "box_hi", "seeds", "function", "max_persistence_H1",
                                                  "max_persistence_H2", "ratio_p1_p2_H1", "ratio_p1_p2_H2",
                                                  "n_seeds_no_dominant_bar", "pass")}
            res["cases"][cid].update(meta(c))
    if "N1" in cases:
        c = cases["N1"]
        res["cases"]["N1"] = {"H1": c["H1"], "H2": c["H2"], "pass": c["pass"], "function": c["function"], **meta(c)}
    if "N2" in cases and "P9c" in cases:
        c = cases["N2"]
        p = pval(c["H1_second_bar_persistence"], cases["P9c"]["null_max_H1"])
        res["cases"]["N2"] = {"H1_top5": c["H1_top5"], "H1_p2_over_p3": c["H1_p2_over_p3"], "H2_top3": c["H2_top3"],
                              "p_value_second_H1_bar_vs_P9c": p,
                              "pass": bool(c["pass_ratio_part"] and p > 0.05), "function": c["function"], **meta(c)}
    killed = [v for k, v in cases.items() if k.startswith("killed_")]
    res["killed_or_crashed_runs"] = killed
    if "P6" in res["cases"]:
        kp6 = [kv for kv in killed if "--case P6" in kv.get("args", "") and "--n" not in kv.get("args", "")]
        pre = res["cases"]["P6"]["variants"].setdefault("preregistered", {"parts": {}, "pass": False})
        pre["killed_runs"] = kp6
        if kp6:
            pre["pass"] = False
            res["cases"]["P6"]["pass"] = False
    lm = {k: v for k, v in cases.items() if k.startswith("P7lm_N")}
    if lm and "P7" in res["cases"]:
        res["cases"]["P7"]["post_hoc_landmarks"] = [
            {"n_landmarks": v["n_landmarks"], "pool": v["pool"], "status": v.get("status"),
             "observed": v["evaluation"]["observed_betti_vector"], "pass": v["pass"],
             "per_dim": v["evaluation"]["per_dim"], "info": v["info"], "command": v["command"],
             "wall_sec_total": v["wall_sec_total"], "peak_rss_mb": v["peak_rss_mb"]}
            for v in sorted(lm.values(), key=lambda v: v["n_landmarks"])]
    for kv in killed:
        if "--case P7" in kv.get("args", ""):
            nn = int(kv["args"].split("--n")[1].split()[0])
            if "P7_N%d" % nn not in p7:
                res["cases"].setdefault("P7", {"sweep": []})["sweep"].append({"N": nn, "status": kv["status"], "exit_code": kv["exit_code"],
                                                                             "pass": False, "command": kv["command"]})
    for cid in ("C0", "F1", "F2", "F1shuf"):
        if cid in cases:
            c = cases[cid]
            res["cases"][cid] = {k: v for k, v in c.items() if k not in ("command", "wall_sec_total", "peak_rss_mb", "loadavg_start", "started")}
            res["cases"][cid].update(meta(c))
    chunks = [v for k, v in cases.items() if k.startswith("F3chunk_")]
    if chunks:
        maps = sorted([m for c in chunks for m in c["maps"]], key=lambda m: m["i"])
        idx = [m["i"] for m in maps]
        f3 = {"n_maps": len(maps), "indices_complete_0_199": idx == list(range(200))}
        if f3["indices_complete_0_199"]:
            ens = [m for m in maps if m["i"] < 100]
            test = [m for m in maps if m["i"] >= 100]
            stats = {}
            ok = True
            for key in ("b0", "b1", "chi"):
                S = np.array([m[key] for m in ens], dtype=float)
                ph, pr = [], []
                for m in test:
                    cs = CMB.coarse_stats(S, np.array(m[key], dtype=float))
                    ph.append(cs["p_value_chi2_survival"]); pr.append(cs["empirical_rank_p"])
                ph = np.array(ph); pr = np.array(pr)
                ksh = kstest(ph, "uniform"); ksr = kstest(pr, "uniform")
                n05 = int(np.sum(ph < 0.05))
                g = bool(ksh.pvalue >= 0.01 and n05 <= 9)
                ok = ok and g
                # POST-HOC diagnostic (added after the b1 KS failure): per-bin ensemble spread of the
                # coarse curves and the largest single-bin contribution to each test chi2.
                Sc = np.array([CMB.coarsen(c, CMB.N_BINS) for c in S])
                sd_bins = Sc.std(axis=0, ddof=1)
                frac_zero = [float(np.mean(Sc[:, j] == 0)) for j in range(CMB.N_BINS)]
                mu = Sc.mean(axis=0)
                zmax = []
                for m in test:
                    tc = CMB.coarsen(np.array(m[key], dtype=float), CMB.N_BINS)
                    z = np.abs(tc - mu) / np.where(sd_bins > 0, sd_bins, np.nan)
                    zmax.append(float(np.nanmax(z)) if np.isfinite(z).any() else float("nan"))
                diag = {"coarse_bin_ensemble_std": sd_bins.tolist(), "coarse_bin_ensemble_mean": mu.tolist(),
                        "coarse_bin_fraction_of_sims_exactly_zero": frac_zero,
                        "n_bins_with_std_below_0.05": int(np.sum(sd_bins < 0.05)),
                        "test_max_abs_z_single_bin_quantiles": np.nanpercentile(zmax, [50, 90, 99, 100]).tolist(),
                        "test_hartlap_p_sorted_first10": np.sort(ph)[:10].tolist()}
                stats[key] = {"post_hoc_bin_diagnostic": diag, "hartlap_p": {"ks_stat": float(ksh.statistic), "ks_p": float(ksh.pvalue), "n_below_0.05": n05,
                                            "literal_TDA7_rule_le5": bool(n05 <= 5), "gate_pass": g,
                                            "deciles": np.percentile(ph, np.arange(0, 101, 10)).tolist()},
                              "empirical_rank_p": {"ks_stat": float(ksr.statistic), "ks_p": float(ksr.pvalue),
                                                   "n_below_0.05": int(np.sum(pr < 0.05))}}
            f3["statistics"] = stats
            f3["pass"] = ok
        f3["function"] = chunks[0]["function"] + " + cmb_tda.coarse_stats"
        f3["chunks"] = [{k: c.get(k) for k in ("start", "count", "command", "wall_sec_total", "peak_rss_mb", "loadavg_start")} for c in chunks]
        if "F3topo" in cases:
            f3["topology_diagnostic"] = {k: v for k, v in cases["F3topo"].items() if k != "command"}
        res["cases"]["F3"] = f3
    order = ["P1", "P2", "P3a", "P3b", "P4", "P5", "P6", "P7", "P8", "P9a", "P9b", "P9c", "P10a", "P10b",
             "F1", "F2", "F3", "N1", "N2", "C0", "F1shuf"]
    res["summary"] = [{"id": k, "pass": res["cases"][k].get("pass")} for k in order if k in res["cases"]]
    res["missing"] = [k for k in order if k not in res["cases"]]
    json.dump(res, open(os.path.join(HERE, "results.json"), "w"), indent=1, default=float)
    print(json.dumps(res["summary"]), "missing:", res["missing"])


# ----------------------------------------------------------------- main
if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--case")
    ap.add_argument("--n", type=int)
    ap.add_argument("--part", default="all", choices=["all", "pipeline", "direct", "rips"])
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--count", type=int, default=50)
    ap.add_argument("--aggregate", action="store_true")
    ARGS = ap.parse_args()
    if ARGS.aggregate:
        aggregate()
        sys.exit(0)
    os.makedirs(CASE_DIR, exist_ok=True)
    t0 = time.time()
    la = os.getloadavg()
    c = ARGS.case
    names = {"P5": "P5_%s" % ARGS.part, "P6": "P6_%s%s" % (ARGS.part, "_N%d" % ARGS.n if ARGS.n else ""),
             "P7": "P7_N%s" % ARGS.n, "P7lm": "P7lm_N%s" % ARGS.n, "F3chunk": "F3chunk_%03d" % ARGS.start}
    name = names.get(c, c)
    dispatch = {"P5": lambda: case_coeff("P5"), "P6": lambda: case_coeff("P6"), "P7": lambda: case_p7(ARGS.n), "P7lm": lambda: case_p7lm(ARGS.n),
                "P8": case_p8, "C0": case_c0, "F1": case_f1, "F1shuf": case_f1shuf, "F2": case_f2, "F3topo": case_f3topo,
                "F3chunk": lambda: case_f3chunk(ARGS.start, ARGS.count)}
    for k in EXPECTED:
        dispatch[k] = (lambda kk: (lambda: case_manifold(kk)))(k)
    for k in ("P9a", "P9b", "P9c"):
        dispatch[k] = (lambda kk: (lambda: case_p9(kk)))(k)
    for k in ("N1", "N2"):
        dispatch[k] = (lambda kk: (lambda: case_neg(kk)))(k)
    if c not in dispatch:
        raise SystemExit("unknown case %r" % c)
    try:
        out = dispatch[c]()
        out.setdefault("status", "ok")
    except (Exception, MemoryError) as e:  # fail closed: record, never pass
        out = {"status": "error", "error": repr(e), "pass": False}
        if c == "P7":
            out["N"] = ARGS.n
    out["command"] = (os.environ.get("SUITE_WRAPPER", "") + " " + sys.executable + " " + " ".join(sys.argv)).strip()
    out["started"] = time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(t0))
    out["loadavg_start"] = la
    out["wall_sec_total"] = time.time() - t0
    out["peak_rss_mb"] = peak_mb()
    out["git_head"] = git_head()
    out["script_sha256"] = hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest()
    out["code_under_test_sha256"] = {CW_REL: sha256(CW_REL), CMB_REL: sha256(CMB_REL)}
    out["versions"] = {"python": sys.version.split()[0], "gudhi": gudhi.__version__, "numpy": np.__version__}
    json.dump(out, open(os.path.join(CASE_DIR, name + ".json"), "w"), indent=1, default=float)
    print(name, "pass=", out.get("pass"), "wall=%.1fs" % out["wall_sec_total"], "rss=%.0fMB" % out["peak_rss_mb"])
