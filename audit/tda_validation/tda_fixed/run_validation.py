#!/usr/bin/env python3
"""Acceptance runs for the two TDA fixes, one step per process.

Every gate and every threshold was declared in fix_expectations.json
(commit 1bc7ebc) BEFORE this script was written and before any number below
was computed.  Nothing here is "proved"; tier of every observed number is X.

Steps (each writes parts/<name>.json; --assemble builds validation_results.json):

  topology    D2 gates 1-3: combinatorial V/E/F at nside 8,16,32,64; gudhi
              Betti at 32 and 64; disk-mask sanity at 32
  smooth      D2 gate 4: old vs fixed complex b0/b1 on one smooth map
  f1          D2 gate 5: the F1 grid case through the fixed filtration code
  p2          P2 regression through the ORIGINAL alpha path + the wrapper
  d1old       D1 columns A/B/C on the committed F3 curves (old topology)
  f3curves    recompute F3 curves with the FIXED topology (--start/--count)
  d1new       D1 column D: the full fixed pipeline
  negctrl     D1 negative control: 20 maps at a different beam
  e5          recompute the committed E5 CMB p-values with the fixed statistic

Run under:  timeout 590 prlimit --as=8589934592 -- <python> run_validation.py --step X
"""
import argparse
import glob
import hashlib
import json
import math
import os
import resource
import subprocess
import sys
import time

import numpy as np
from scipy.stats import chi2 as _scipy_chi2


def _chi2sf(x, df):
    return np.asarray(_scipy_chi2.sf(x, df=df), dtype=float)

HERE = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.join(HERE, "..", "simple_suite")
WT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
PARTS = os.path.join(HERE, "parts")
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..")))

from tda_fixed import cmb_topology as ct          # noqa: E402
from tda_fixed import stats as fs                 # noqa: E402
from tda_fixed import pointcloud as pc            # noqa: E402

CMB_REL = "audit/reverse_zero/E5-cmb-tda/cmb_tda.py"
CW_REL = "audit/reverse_zero/E5-cosmic-web-tda-scaled/cosmic_web_tda_scaled.py"

NU_MIN, NU_MAX, NU_N = -4.0, 4.0, 41
N_BINS = 8
F3_SEED_BASE = 30000          # suite expectations.json case F3
NEG_SEED_BASE = 40000         # declared in fix_expectations.json
NEG_FWHM_DEG = 6.0            # declared: control beam (F3 uses 3 deg)
N_NEG = 20


def sha256(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def load_original():
    """Import the ORIGINAL pipeline file under a non-__main__ name."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("cmb_tda_original", os.path.join(WT, CMB_REL))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def cl_for(fwhm_deg, nside=64):
    """The F3 spectrum: C_l = exp(-l(l+1) sigma_b^2)/(l(l+1)), 2 <= l <= 3*nside-1."""
    lmax = 3 * nside - 1
    ell = np.arange(lmax + 1)
    sb = math.radians(fwhm_deg) / math.sqrt(8 * math.log(2))
    cl = np.zeros(lmax + 1)
    cl[2:] = np.exp(-ell[2:] * (ell[2:] + 1) * sb ** 2) / (ell[2:] * (ell[2:] + 1))
    return cl


def nu_grid():
    return np.linspace(NU_MIN, NU_MAX, NU_N)


def write(name, obj):
    os.makedirs(PARTS, exist_ok=True)
    obj["_peak_rss_mb"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
    obj["_command"] = " ".join([sys.executable] + sys.argv)
    obj["_generated"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    json.dump(obj, open(os.path.join(PARTS, name + ".json"), "w"), indent=1, default=float)
    print("wrote", name, "rss=%.0fMB" % obj["_peak_rss_mb"])


# --------------------------------------------------------------- D2 gates
def step_topology():
    import healpy as hp
    out = {"gate": "D2 gates 1-3", "levels": {}, "tier": "X"}
    for nside in (8, 16, 32, 64):
        t0 = time.time()
        unm, e, tr, info = ct.build_topology_fixed(np.ones(hp.nside2npix(nside), np.uint8),
                                                   nside, return_info=True)
        V, E, F = info["n_vertices"], info["n_edges"], info["n_triangles"]
        rec = {k: info[k] for k in (
            "npix", "n_missing_neighbour_entries", "missing_neighbour_entries_by_direction",
            "n_missing_at_edge_sharing_positions", "n_quad_corners", "n_tri_corners",
            "predicted_n_quad_corners_V_minus_6", "predicted_n_tri_corners",
            "predicted_full_sky_edges_2V_plus_Q", "predicted_full_sky_faces_2Q_plus_T3",
            "degenerate_quad_rows", "degenerate_tri_rows",
            "dedup_factor_quads", "dedup_factor_tris")}
        rec.update({
            "n_vertices": V, "n_edges": E, "n_triangles": F,
            "euler_char_V_minus_E_plus_F": V - E + F,
            "gate_1_combinatorial_pass": bool(
                V - E + F == 2
                and info["n_quad_corners"] == info["predicted_n_quad_corners_V_minus_6"]
                and info["n_tri_corners"] == 8
                and E == info["predicted_full_sky_edges_2V_plus_Q"]
                and F == info["predicted_full_sky_faces_2Q_plus_T3"]),
            "build_sec": time.time() - t0,
        })
        if nside in (32, 64):
            t1 = time.time()
            b = ct.complex_betti(V, e, tr)
            rec["gudhi_betti"] = b
            rec["gate_2_betti_1_0_1_pass"] = bool(b == [1, 0, 1])
            rec["betti_sec"] = time.time() - t1
        out["levels"][str(nside)] = rec

    # --- the OLD complex on the same full sky, for the before column
    orig = load_original()
    for nside in (32, 64):
        t0 = time.time()
        unm_o, e_o, tr_o = orig.build_topology(np.ones(hp.nside2npix(nside), np.uint8), nside)
        V, E, F = int(unm_o.size), int(e_o.shape[0]), int(tr_o.shape[0])
        b = ct.complex_betti(V, e_o, tr_o)
        out["levels"][str(nside)]["OLD_build_topology"] = {
            "n_vertices": V, "n_edges": E, "n_triangles": F,
            "euler_char_V_minus_E_plus_F": V - E + F, "gudhi_betti": b,
            "wall_sec": time.time() - t0,
            "note": "original cmb_tda.build_topology, unchanged, for the before column"}

    # --- gate 3: disk mask at nside 32 (all pixels within 30 deg of the north pole)
    nside = 32
    vec = np.array(hp.pix2vec(nside, np.arange(hp.nside2npix(nside))))
    mask = (vec[2] >= math.cos(math.radians(30.0))).astype(np.uint8)
    unm, e, tr, info = ct.build_topology_fixed(mask, nside, return_info=True)
    b = ct.complex_betti(int(unm.size), e, tr)
    V, E, F = int(unm.size), int(e.shape[0]), int(tr.shape[0])
    unm_o, e_o, tr_o = orig.build_topology(mask, nside)
    b_o = ct.complex_betti(int(unm_o.size), e_o, tr_o)
    out["gate_3_disk_mask"] = {
        "mask": "nside 32, pixels with z >= cos(30 deg) (a polar cap)",
        "n_unmasked": V, "n_edges": E, "n_triangles": F,
        "euler_char_V_minus_E_plus_F": V - E + F,
        "expected_betti_disc": [1, 0, 0], "gudhi_betti_fixed": b,
        "gate_3_pass": bool(b == [1, 0, 0]),
        "OLD_build_topology_betti": b_o,
        "OLD_n_edges": int(e_o.shape[0]), "OLD_n_triangles": int(tr_o.shape[0]),
        "OLD_euler_char": int(unm_o.size) - int(e_o.shape[0]) + int(tr_o.shape[0]),
    }
    write("topology", out)


def step_smooth():
    """D2 gate 4: b0/b1 curves on one smooth map, old complex vs fixed complex."""
    import healpy as hp
    orig = load_original()
    nside = 64
    cl = cl_for(3.0, nside)
    np.random.seed(F3_SEED_BASE)
    m = hp.synfast(cl, nside=nside, new=True)
    mask = np.ones(hp.nside2npix(nside), np.uint8)
    nu = nu_grid()

    t0 = time.time()
    unm_o, e_o, tr_o = orig.build_topology(mask, nside)
    b0_o, b1_o, chi_o = orig.betti_curves_from_topology(m, unm_o, e_o, tr_o, nu, sublevel=True)
    t_old = time.time() - t0

    t1 = time.time()
    unm_n, e_n, tr_n = ct.build_topology_fixed(mask, nside)
    r = ct.betti_curves_from_topology(m, unm_n, e_n, tr_n, nu, sublevel=True)
    t_new = time.time() - t1

    b0_n, b1_n = r["b0"], r["b1"]
    d0 = np.abs(b0_n - b0_o).max()
    d1 = np.abs(b1_n - b1_o).max()
    tol0 = 0.15 * b0_o.max()
    tol1 = 0.15 * b1_o.max()
    step = nu[1] - nu[0]
    out = {
        "gate": "D2 gate 4: smooth-field b0/b1 regression",
        "map": "hp.synfast(cl(FWHM 3 deg), nside=64, new=True), np.random.seed(%d)" % F3_SEED_BASE,
        "nu_grid": nu.tolist(),
        "old": {"b0": b0_o.tolist(), "b1": b1_o.tolist(), "euler_chi_as_returned_b0_minus_b1": chi_o.tolist(),
                "n_edges": int(e_o.shape[0]), "n_triangles": int(tr_o.shape[0]), "wall_sec": t_old},
        "fixed": {"b0": b0_n.tolist(), "b1": b1_n.tolist(),
                  "b0_minus_b1": r["b0_minus_b1"].tolist(),
                  "euler_char_true": r["euler_char_true"].tolist(),
                  "b2_implied": r["b2_implied"].tolist(),
                  "n_edges": int(e_n.shape[0]), "n_triangles": int(tr_n.shape[0]), "wall_sec": t_new},
        "max_abs_delta_b0": float(d0), "max_abs_delta_b1": float(d1),
        "max_b0_old": float(b0_o.max()), "max_b1_old": float(b1_o.max()),
        "declared_tolerance_0.15_of_peak": {"b0": float(tol0), "b1": float(tol1)},
        "rms_delta_b0": float(np.sqrt(np.mean((b0_n - b0_o) ** 2))),
        "rms_delta_b1": float(np.sqrt(np.mean((b1_n - b1_o) ** 2))),
        "peak_nu_b0_old": float(nu[int(np.argmax(b0_o))]), "peak_nu_b0_fixed": float(nu[int(np.argmax(b0_n))]),
        "peak_nu_b1_old": float(nu[int(np.argmax(b1_o))]), "peak_nu_b1_fixed": float(nu[int(np.argmax(b1_n))]),
        "grid_step": float(step),
        "gate_4_pass": bool(d0 <= tol0 and d1 <= tol1
                            and abs(nu[int(np.argmax(b0_o))] - nu[int(np.argmax(b0_n))]) <= step + 1e-9
                            and abs(nu[int(np.argmax(b1_o))] - nu[int(np.argmax(b1_n))]) <= step + 1e-9),
        "euler_check_at_final_nu": {
            "old_euler_chi_key_b0_minus_b1": int(chi_o[-1]),
            "fixed_b0_minus_b1": int(r["b0_minus_b1"][-1]),
            "fixed_euler_char_true": int(r["euler_char_true"][-1]),
            "note": "the full sphere has chi = 2 and b0-b1 = 1; the original returned 1 under the name euler_chi"},
        "tier": "X",
    }
    write("smooth", out)


def step_f1():
    """D2 gate 5: the F1 grid case must be untouched by the topology fix."""
    orig = load_original()
    n, sig = 256, 10.0
    rng = np.random.default_rng(21)
    C = []
    while len(C) < 7:
        c = rng.uniform(30, 226, 2)
        if all(np.linalg.norm(c - d) >= 50 for d in C):
            C.append(c)
    C = np.array(C)
    ii, jj = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    f = np.zeros((n, n))
    for c in C:
        f -= np.exp(-((ii - c[0]) ** 2 + (jj - c[1]) ** 2) / (2 * sig ** 2))
    idx = np.arange(n * n).reshape(n, n)
    e = np.r_[np.c_[idx[:, :-1].ravel(), idx[:, 1:].ravel()],
              np.c_[idx[:-1, :].ravel(), idx[1:, :].ravel()],
              np.c_[idx[:-1, :-1].ravel(), idx[1:, 1:].ravel()]].astype(np.int64)
    tr = np.r_[np.c_[idx[:-1, :-1].ravel(), idx[:-1, 1:].ravel(), idx[1:, 1:].ravel()],
               np.c_[idx[:-1, :-1].ravel(), idx[1:, :-1].ravel(), idx[1:, 1:].ravel()]].astype(np.int64)
    temp = f.ravel()
    nu = np.linspace(temp.min() / temp.std(), temp.max() / temp.std(), 4001)
    v = np.arange(n * n)

    t0 = time.time()
    b0_o, b1_o, chi_o = orig.betti_curves_from_topology(temp, v, e, tr, nu, sublevel=True)
    t_old = time.time() - t0
    t1 = time.time()
    r = ct.betti_curves_from_topology(temp, v, e, tr, nu, sublevel=True)
    t_new = time.time() - t1

    committed = json.load(open(os.path.join(SUITE, "cases", "F1.json")))["pipeline"]
    out = {
        "gate": "D2 gate 5: F1 (256x256 grid, Freudenthal, 7 Gaussian wells, seed 21)",
        "committed_suite_result": {k: committed[k] for k in
                                   ("max_nu_b0", "b0_at_max_nu", "b1_at_max_nu", "max_nu_b1",
                                    "n_edges", "n_triangles", "nu_points")},
        "original_function_rerun": {"max_nu_b0": int(b0_o.max()), "b0_at_max_nu": int(b0_o[-1]),
                                    "b1_at_max_nu": int(b1_o[-1]), "max_nu_b1": int(b1_o.max()),
                                    "wall_sec": t_old},
        "fixed_function": {"max_nu_b0": int(r["b0"].max()), "b0_at_max_nu": int(r["b0"][-1]),
                           "b1_at_max_nu": int(r["b1"][-1]), "max_nu_b1": int(r["b1"].max()),
                           "wall_sec": t_new},
        "curves_identical_b0": bool(np.array_equal(b0_o, r["b0"])),
        "curves_identical_b1": bool(np.array_equal(b1_o, r["b1"])),
        "euler_char_true_at_final_nu": int(r["euler_char_true"][-1]),
        "b0_minus_b1_at_final_nu": int(r["b0_minus_b1"][-1]),
        "euler_note": "a 256x256 disc has chi = 1 and b0-b1 = 1; both agree here because b2 = 0",
        "gate_5_pass": bool(np.array_equal(b0_o, r["b0"]) and np.array_equal(b1_o, r["b1"])
                            and int(r["b0"].max()) == 7 and int(r["b0"][-1]) == 1
                            and int(r["b1"][-1]) == 0
                            and int(r["b0"].max()) == committed["max_nu_b0"]),
        "tier": "X",
    }
    write("f1", out)


def step_p2():
    """P2 regression: the alpha path is not implicated; check it is undisturbed."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("cw_original", os.path.join(WT, CW_REL))
    cw = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cw)

    X = pc.sphere_cloud_P2()
    t0 = time.time()
    res_o, by_dim_o = cw.alpha_persistence(X, float("inf"), "P2_regression", HERE)
    t_o = time.time() - t0
    t1 = time.time()
    res_n, by_dim_n = pc.alpha_persistence(X, float("inf"), "P2_regression")
    t_n = time.time() - t1

    same = all(np.allclose(np.sort(np.array(by_dim_o[k]).reshape(-1, 2), axis=0),
                           np.sort(np.array(by_dim_n[k]).reshape(-1, 2), axis=0),
                           equal_nan=True)
               for k in (0, 1, 2) if by_dim_o[k])
    # dominant/zero evaluation, same rules as the suite
    def dominant(bars, cap):
        b = np.array(bars, float).reshape(-1, 2)
        d = np.where(np.isfinite(b[:, 1]), np.minimum(b[:, 1], cap), cap)
        return np.sort(d - b[:, 0])[::-1]
    fin = np.concatenate([np.array(by_dim_o[k], float).reshape(-1, 2)[:, 1] for k in (0, 1, 2) if by_dim_o[k]])
    cap = float(fin[np.isfinite(fin)].max())
    p = {k: dominant(by_dim_o[k], cap) for k in (0, 1, 2)}
    P_all = max(p[1][0] if p[1].size else 0.0, p[2][0] if p[2].size else 0.0)
    obs = [int(np.sum(p[k] >= 0.2 * P_all)) for k in (0, 1, 2)]
    committed = json.load(open(os.path.join(SUITE, "cases", "P2.json")))
    out = {
        "gate": "P2 regression (sphere S2 through the alpha path)",
        "numbers_produced_by": "cosmic_web_tda_scaled.alpha_persistence (the ORIGINAL pipeline function)",
        "n_points": int(X.shape[0]),
        "expected_betti": [1, 0, 1],
        "observed_betti_vector": obs,
        "betti_numbers_at_truncation_original": res_o["betti_numbers_at_truncation"],
        "betti_numbers_at_truncation_wrapper": res_n["betti_numbers_at_truncation"],
        "top_persistence_H0": p[0][:3].tolist(), "top_persistence_H1": p[1][:3].tolist(),
        "top_persistence_H2": p[2][:3].tolist(),
        "cap": cap, "P_all": float(P_all),
        "wrapper_bars_identical_to_original": bool(same),
        "committed_suite_observed_betti_vector": committed["primary"]["evaluation"]["observed_betti_vector"],
        "matches_committed": bool(obs == committed["primary"]["evaluation"]["observed_betti_vector"]),
        "wall_sec_original": t_o, "wall_sec_wrapper": t_n,
        "pass": bool(obs == [1, 0, 1] and same),
        "tier": "X",
    }
    write("p2", out)


# --------------------------------------------------------------- D1 gates
def load_committed_f3_curves():
    maps = []
    for fn in sorted(glob.glob(os.path.join(SUITE, "cases", "F3chunk_*.json"))):
        maps += json.load(open(fn))["maps"]
    maps = sorted(maps, key=lambda m: m["i"])
    assert [m["i"] for m in maps] == list(range(200)), "committed F3 curves incomplete"
    return maps


def calibration(ens, test, label, mode, orig=None, thresholds=None):
    """KS-against-uniform and #{p<0.05} for one column.

    mode: 'A' original coarse_stats | 'B' minimal fix | 'C' full fix
    ens/test: dicts key -> (n, L) arrays
    """
    from scipy.stats import kstest
    out = {"column": label, "mode": mode, "keys": {}}
    for key in ("b0", "b1", "chi"):
        S, T = ens[key], test[key]
        ph, pr, extra = [], [], []
        for row in T:
            if mode == "A":
                cs = orig.coarse_stats(S, row, n_bins=N_BINS)
                ph.append(cs["p_value_chi2_survival"]); pr.append(cs["empirical_rank_p"])
                extra.append({"df": cs["n_bins"]})
            else:
                cs = fs.coarse_stats_fixed(S, row, n_bins=N_BINS,
                                           drop_atomic=(mode == "C"),
                                           full_curve_diagnostics=False)
                ph.append(cs["p_value_chi2_survival"]); pr.append(cs["empirical_rank_p"])
                extra.append({"df": cs["df"], "kept": cs["n_bins_kept"],
                              "dropped": cs["dropped_bin_indices"],
                              "reasons": cs["dropped_bin_reasons"],
                              "cond": cs.get("kept_block_condition_number")})
        ph = np.array(ph, float); pr = np.array(pr, float)
        ksh = kstest(ph, "uniform"); ksr = kstest(pr, "uniform")
        n05 = int(np.sum(ph < 0.05))
        rec = {
            "n_test_maps": int(len(ph)),
            "chi2_p": {"ks_p_vs_uniform": float(ksh.pvalue), "ks_stat": float(ksh.statistic),
                       "n_below_0.05": n05, "min_p": float(ph.min()), "median_p": float(np.median(ph)),
                       "deciles": np.percentile(ph, np.arange(0, 101, 10)).tolist()},
            "rank_p": {"ks_p_vs_uniform": float(ksr.pvalue), "ks_stat": float(ksr.statistic),
                       "n_below_0.05": int(np.sum(pr < 0.05)), "min_p": float(pr.min()),
                       "median_p": float(np.median(pr))},
            "gate_ks_gt_0.05_and_n05_le_9": bool(ksh.pvalue > 0.05 and n05 <= 9),
        }
        if mode != "A":
            rec["bin_selection"] = {"df": extra[0]["df"], "n_bins_kept": extra[0]["kept"],
                                    "dropped_bin_indices": extra[0]["dropped"],
                                    "dropped_bin_reasons": extra[0]["reasons"],
                                    "kept_block_condition_number": extra[0]["cond"],
                                    "identical_for_all_test_maps": bool(
                                        all(e["dropped"] == extra[0]["dropped"] for e in extra))}
        else:
            rec["bin_selection"] = {"df": extra[0]["df"], "note": "original: df = n_bins always"}
        out["keys"][key] = rec

    if thresholds:
        sens = {}
        for th in thresholds:
            sens[str(th)] = {}
            for key in ("b0", "b1", "chi"):
                S, T = ens[key], test[key]
                ps = []
                for row in T:
                    cs = fs.coarse_stats_fixed(S, row, n_bins=N_BINS, drop_atomic=True,
                                               atomic_modal_mass_max=th, with_rank_p=False,
                                               full_curve_diagnostics=False)
                    ps.append(cs["p_value_chi2_survival"])
                ps = np.array(ps, float)
                k = kstest(ps, "uniform")
                kept = fs.select_bins(np.array([fs.coarsen(c, N_BINS) for c in S]), th)[0]
                sens[str(th)][key] = {"ks_p": float(k.pvalue), "n_below_0.05": int(np.sum(ps < 0.05)),
                                      "n_bins_kept": len(kept), "kept_bin_indices": kept,
                                      "gate_pass": bool(k.pvalue > 0.05 and np.sum(ps < 0.05) <= 9)}
        out["atomic_threshold_sensitivity"] = sens
    return out


def _split(maps):
    ens = {k: np.array([m[k] for m in maps if m["i"] < 100], float) for k in ("b0", "b1", "chi")}
    test = {k: np.array([m[k] for m in maps if m["i"] >= 100], float) for k in ("b0", "b1", "chi")}
    return ens, test


def step_d1old():
    orig = load_original()
    maps = load_committed_f3_curves()
    ens, test = _split(maps)
    out = {"case": "F3 (nside 64), curves as committed by the suite (OLD build_topology)",
           "curves_from": "audit/tda_validation/simple_suite/cases/F3chunk_*.json",
           "tier": "X", "columns": {}}
    out["columns"]["A_before"] = calibration(ens, test, "A: before (original coarse_stats)", "A", orig=orig)
    out["columns"]["B_minimal"] = calibration(ens, test, "B: dead-bin drop + df=kept + eigen conditioning", "B")
    out["columns"]["C_full"] = calibration(ens, test, "C: B + atomic-bin drop (the declared fix)", "C",
                                           thresholds=[0.3, 0.4, 0.5, 0.6, 0.7, 0.9, 1.0])
    # reproduce the suite's own before numbers as a determinism check
    committed = json.load(open(os.path.join(SUITE, "results.json")))["cases"]["F3"]["statistics"]
    out["reproduces_committed_before_numbers"] = {
        k: {"committed_ks_p": committed[k]["hartlap_p"]["ks_p"],
            "recomputed_ks_p": out["columns"]["A_before"]["keys"][k]["chi2_p"]["ks_p_vs_uniform"],
            "match": bool(abs(committed[k]["hartlap_p"]["ks_p"]
                              - out["columns"]["A_before"]["keys"][k]["chi2_p"]["ks_p_vs_uniform"]) < 1e-9)}
        for k in ("b0", "b1", "chi")}
    write("d1old", out)


def step_f3curves(start, count, fwhm=3.0, seed_base=F3_SEED_BASE, tag="f3curves"):
    """Recompute F3-style curves on the FIXED complex."""
    import healpy as hp
    t0 = time.time()
    nside = 64
    unm, e, tr, info = ct.build_topology_fixed(np.ones(hp.nside2npix(nside), np.uint8),
                                               nside, return_info=True)
    t_topo = time.time() - t0
    cl = cl_for(fwhm, nside)
    nu = nu_grid()
    maps = []
    for i in range(start, start + count):
        np.random.seed(seed_base + i)
        m = hp.synfast(cl, nside=nside, new=True)
        r = ct.betti_curves_from_topology(m, unm, e, tr, nu, sublevel=True)
        maps.append({"i": i, "seed": seed_base + i,
                     "b0": r["b0"].tolist(), "b1": r["b1"].tolist(),
                     "chi": r["b0_minus_b1"].tolist(),
                     "euler_char_true": r["euler_char_true"].tolist()})
    write("%s_%03d" % (tag, start), {
        "function": "tda_fixed.cmb_topology.build_topology_fixed + betti_curves_from_topology",
        "nside": nside, "fwhm_deg": fwhm, "seed_base": seed_base,
        "start": start, "count": count, "topology_sec": t_topo,
        "n_vertices": info["n_vertices"], "n_edges": info["n_edges"], "n_triangles": info["n_triangles"],
        "euler_char_V_minus_E_plus_F": info["euler_char_V_minus_E_plus_F"],
        "wall_sec": time.time() - t0, "maps": maps, "tier": "X",
        "chi_key_note": "'chi' is b0-b1 (the quantity the suite called chi, kept for a like-for-like "
                        "comparison with column A); 'euler_char_true' is the true Euler characteristic"})


def _load_parts(tag):
    maps = []
    for fn in sorted(glob.glob(os.path.join(PARTS, "%s_*.json" % tag))):
        maps += json.load(open(fn))["maps"]
    return sorted(maps, key=lambda m: m["i"])


def step_d1new():
    maps = _load_parts("f3curves")
    assert [m["i"] for m in maps] == list(range(200)), \
        "fixed-topology F3 curves incomplete: %d maps" % len(maps)
    ens, test = _split(maps)
    out = {"case": "F3 (nside 64), curves recomputed with the FIXED topology",
           "curves_from": "parts/f3curves_*.json", "tier": "X", "columns": {}}
    out["columns"]["D_full_fixed_pipeline"] = calibration(
        ens, test, "D: fixed topology + fixed statistic", "C",
        thresholds=[0.3, 0.4, 0.5, 0.6, 0.7, 0.9, 1.0])
    orig = load_original()
    out["columns"]["D0_fixed_topology_old_stats"] = calibration(
        ens, test, "D0: fixed topology + ORIGINAL coarse_stats (isolates the topology change)",
        "A", orig=orig)
    # the true Euler curve as its own statistic, only available after the fix
    ens_e = {"b0": np.array([m["euler_char_true"] for m in maps if m["i"] < 100], float)}
    test_e = {"b0": np.array([m["euler_char_true"] for m in maps if m["i"] >= 100], float)}
    ens_e["b1"] = ens_e["chi"] = ens_e["b0"]
    test_e["b1"] = test_e["chi"] = test_e["b0"]
    ce = calibration(ens_e, test_e, "E: TRUE Euler characteristic curve, fixed statistic", "C")
    out["columns"]["E_true_euler_curve"] = {"column": ce["column"], "mode": "C",
                                            "keys": {"euler_char_true": ce["keys"]["b0"]}}
    write("d1new", out)


def step_negctrl():
    """D1 negative control: maps from a DIFFERENT law must be rejected."""
    from scipy.stats import kstest
    ens_maps = [m for m in _load_parts("f3curves") if m["i"] < 100]
    assert len(ens_maps) == 100, len(ens_maps)
    ctrl = _load_parts("negctrl")
    assert len(ctrl) == N_NEG, "negative-control curves incomplete: %d" % len(ctrl)
    out = {"gate": "D1 negative control",
           "construction": "20 maps, same C_l form with Gaussian beam FWHM %.1f deg (F3 uses 3.0), "
                           "np.random.seed(%d + i), i = 0..19, nside 64, FIXED complex, scored "
                           "against the same 100-map 3-deg ensemble" % (NEG_FWHM_DEG, NEG_SEED_BASE),
           "declared_gate": "at least 15 of 20 control maps have chi2 p < 0.05 under the FIXED statistic",
           "tier": "X", "keys": {}}
    for key in ("b0", "b1", "chi"):
        S = np.array([m[key] for m in ens_maps], float)
        ph, pr = [], []
        for m in ctrl:
            cs = fs.coarse_stats_fixed(S, np.array(m[key], float), n_bins=N_BINS,
                                       drop_atomic=True, full_curve_diagnostics=False)
            ph.append(cs["p_value_chi2_survival"]); pr.append(cs["empirical_rank_p"])
        ph = np.array(ph, float); pr = np.array(pr, float)
        out["keys"][key] = {
            "n_control_maps": len(ph),
            "chi2_p_values": ph.tolist(), "rank_p_values": pr.tolist(),
            "median_chi2_p": float(np.median(ph)), "n_chi2_p_below_0.05": int(np.sum(ph < 0.05)),
            "n_rank_p_below_0.05": int(np.sum(pr < 0.05)),
            "ks_p_vs_uniform": float(kstest(ph, "uniform").pvalue),
            "gate_pass_15_of_20": bool(np.sum(ph < 0.05) >= 15),
        }
    write("negctrl", out)


def step_e5():
    """Recompute the committed E5 CMB p-values with the fixed statistic."""
    orig = load_original()
    npz = os.path.join(WT, "audit/reverse_zero/E5-cmb-tda/chunks/null_merged.npz")
    c = np.load(npz)
    out = {"what": "audit/reverse_zero/E5-cmb-tda committed run, statistics recomputed",
           "ensemble": npz, "n_sims": int(c["sim_b0"].shape[0]),
           "nside_of_that_run": 128, "mask": "WMAP KQ85 (ud_graded, >= 0.5)",
           "SCOPE": "ONLY the statistics defect (D1) is recomputed. The curves themselves were "
                    "produced by the DEFECTIVE build_topology at nside 128; recomputing them means "
                    "regenerating 100 synfast sims x 2 filtration directions at nside 128 and is "
                    "NOT done here. This is not a revalidation of E5.",
           "tier": "X", "statistics": {}}
    for direction, suffix in (("sublevel", ""), ("superlevel", "_super")):
        for key, arr, dat in (("b0", "sim_b0" + suffix, "b0_d" + suffix),
                              ("b1", "sim_b1" + suffix, "b1_d" + suffix),
                              ("euler_chi", "sim_chi" + suffix, "chi_d" + suffix)):
            S = np.asarray(c[arr], float)
            d = np.asarray(c[dat], float)
            old = orig.coarse_stats(S, d, n_bins=N_BINS)
            new = fs.coarse_stats_fixed(S, d, n_bins=N_BINS, drop_atomic=True,
                                        full_curve_diagnostics=False)
            minimal = fs.coarse_stats_fixed(S, d, n_bins=N_BINS, drop_atomic=False,
                                            full_curve_diagnostics=False)
            out["statistics"]["%s_%s" % (direction, key)] = {
                "OLD_chi2": old["data_chi2_hartlap"], "OLD_p_chi2": old["p_value_chi2_survival"],
                "OLD_rank_p": old["empirical_rank_p"], "OLD_df": old["n_bins"],
                "MINIMAL_p_chi2": minimal["p_value_chi2_survival"], "MINIMAL_df": minimal["df"],
                "NEW_chi2": new["data_chi2_hartlap"], "NEW_p_chi2": new["p_value_chi2_survival"],
                "NEW_rank_p": new["empirical_rank_p"], "NEW_df": new["df"],
                "NEW_n_bins_kept": new["n_bins_kept"],
                "NEW_dropped_bin_indices": new["dropped_bin_indices"],
                "NEW_dropped_bin_reasons": new["dropped_bin_reasons"],
                "NEW_test_differs_in_dropped_bin": new["test_differs_in_dropped_bin"],
                "delta_p_chi2": (None if old["p_value_chi2_survival"] is None
                                 else new["p_value_chi2_survival"] - old["p_value_chi2_survival"]),
                "crosses_0.05": bool((old["p_value_chi2_survival"] < 0.05)
                                     != (new["p_value_chi2_survival"] < 0.05)),
            }
    write("e5", out)


def step_calibcheck():
    """Is #{p < 0.05} = 11 or 12 of 100 evidence that the rank p-value is
    NOT exchangeable, or is it the spread induced by 100 test maps sharing
    ONE 100-member ensemble?

    Discriminating check on synthetic multivariate normal data, where the
    statistic's assumptions hold exactly by construction: draw an ensemble of
    100 and 100 test vectors from the SAME normal law, score with the same
    coarse_stats_fixed code path, and record #{rank p < 0.05} per ensemble.
    Repeat over many INDEPENDENT ensembles.  If the pooled rate is 5/101 and
    the per-ensemble count has a standard deviation well above the binomial
    one, then the count is correlated across test maps and 11-12 is ordinary.
    If the pooled rate itself exceeds 5/101, the implementation is not
    exchangeable and that is a defect.
    """
    rng = np.random.default_rng(777)
    p, n, n_test, nrep = 6, 100, 100, 100
    A = rng.normal(size=(p, p))
    Sig = A @ A.T + np.eye(p)
    L = np.linalg.cholesky(Sig)
    tau = fs.EIGEN_TRUNCATION_TAU

    def loo_scores(P):
        """All m leave-one-out Mahalanobis scores of an (m, p) pool, by
        rank-one downdate of the pool's sum and second-moment matrix."""
        m = P.shape[0]
        S1 = P.sum(axis=0)
        Q = P.T @ P
        nn = m - 1
        out = np.empty(m)
        for i in range(m):
            xi = P[i]
            mu = (S1 - xi) / nn
            C = (Q - np.outer(xi, xi) - nn * np.outer(mu, mu)) / (nn - 1)
            cinv, k, _ = fs.conditioned_inverse(C, tau)
            h = (nn - k - 2) / (nn - 1) if nn > k + 2 else 1.0
            d = xi - mu
            out[i] = float(d @ (cinv * h) @ d)
        return out

    counts = []
    pooled = []
    counts_chi2 = []
    pooled_chi2 = []
    t0 = time.time()
    for _ in range(nrep):
        E = rng.normal(size=(n, p)) @ L.T
        T = rng.normal(size=(n_test, p)) @ L.T
        # chi2 branch: one ensemble mean/cov for all test vectors (as in the gate)
        mu = E.mean(axis=0)
        cinv, k, _ = fs.conditioned_inverse(np.cov(E, rowvar=False), tau)
        h = (n - k - 2) / (n - 1)
        D = T - mu
        chi2v = np.einsum("ij,jk,ik->i", D, cinv * h, D)
        pc2 = _chi2sf(chi2v, k)
        counts_chi2.append(int((pc2 < 0.05).sum()))
        pooled_chi2.append(pc2)
        ps = np.empty(n_test)
        for j in range(n_test):
            sc = loo_scores(np.vstack([E, T[j][None, :]]))
            ps[j] = (np.sum(sc[:-1] >= sc[-1]) + 1) / (n + 1)
        counts.append(int((ps < 0.05).sum()))
        pooled.append(ps)
    counts = np.array(counts)
    pooled = np.concatenate(pooled)
    counts_chi2 = np.array(counts_chi2)
    pooled_chi2 = np.concatenate(pooled_chi2)
    target = float(np.sum(np.arange(1, n + 2) / (n + 1) < 0.05) / (n + 1))
    binom_sd = math.sqrt(n_test * target * (1 - target))
    out = {
        "question": "is #{rank p < 0.05} = 11-12 of 100 a broken exchangeability, or the "
                    "spread from 100 test vectors sharing ONE 100-member ensemble?",
        "construction": "multivariate normal, p = %d, ensemble n = %d, %d test vectors per "
                        "ensemble, %d INDEPENDENT ensembles, np.random.default_rng(777); "
                        "same fs.conditioned_inverse / Hartlap / pooled-exchangeable code path"
                        % (p, n, n_test, nrep),
        "exact_uniform_target_P_rank_p_below_0.05": target,
        "pooled_rate_over_all_%d_test_vectors" % (nrep * n_test): float((pooled < 0.05).mean()),
        "per_ensemble_count_mean": float(counts.mean()),
        "per_ensemble_count_sd": float(counts.std(ddof=1)),
        "binomial_sd_if_independent": binom_sd,
        "variance_inflation_factor": float((counts.std(ddof=1) / binom_sd) ** 2),
        "per_ensemble_count_quantiles_5_50_90_95_99": np.percentile(counts, [5, 50, 90, 95, 99]).tolist(),
        "per_ensemble_count_max": int(counts.max()),
        "empirical_P_count_ge_11": float((counts >= 11).mean()),
        "empirical_P_count_ge_12": float((counts >= 12).mean()),
        "empirical_P_count_ge_13": float((counts >= 13).mean()),
        "counts": counts.tolist(),
        "chi2_p_branch": {
            "note": "the branch the acceptance gate actually uses: ONE ensemble mean/cov per "
                    "ensemble, chi2 survival with df = retained rank, Gaussian data so the "
                    "reference distribution is exact",
            "pooled_rate_P_chi2_p_below_0.05": float((pooled_chi2 < 0.05).mean()),
            "per_ensemble_count_mean": float(counts_chi2.mean()),
            "per_ensemble_count_sd": float(counts_chi2.std(ddof=1)),
            "binomial_sd_if_independent": math.sqrt(n_test * 0.05 * 0.95),
            "variance_inflation_factor": float((counts_chi2.std(ddof=1) / math.sqrt(n_test * 0.05 * 0.95)) ** 2),
            "per_ensemble_count_quantiles_5_50_90_95_99": np.percentile(counts_chi2, [5, 50, 90, 95, 99]).tolist(),
            "per_ensemble_count_max": int(counts_chi2.max()),
            "empirical_P_count_ge_10": float((counts_chi2 >= 10).mean()),
            "empirical_P_count_ge_11": float((counts_chi2 >= 11).mean()),
            "empirical_P_count_ge_12": float((counts_chi2 >= 12).mean()),
            "counts": counts_chi2.tolist(),
            "implication_for_the_gate": "this is the false-failure rate of the declared "
                                        "#{p<0.05} <= 9 clause for a PERFECTLY calibrated "
                                        "statistic, caused by the 100 test maps sharing one "
                                        "ensemble; it is a property of the F3 design, not of the fix",
        },
        "wall_sec": time.time() - t0,
        "tier": "X",
    }
    out["verdict"] = (
        "pooled rate matches the exact uniform target to within Monte-Carlo error, so the "
        "pooled-exchangeable rank p is calibrated; the per-ensemble count is over-dispersed "
        "relative to binomial by the stated factor because the 100 test vectors share one "
        "ensemble"
        if abs(out["pooled_rate_over_all_%d_test_vectors" % (nrep * n_test)] - target) < 3 * math.sqrt(
            target * (1 - target) / (nrep * n_test)) * out["variance_inflation_factor"] ** 0.5
        else "pooled rate differs from the exact uniform target: the implementation is NOT "
             "exchangeable and this is a defect")
    write("calibcheck", out)


def step_gate4full():
    """D2 gate 4 over all 200 maps instead of one.

    The committed suite curves (cases/F3chunk_*.json, OLD complex) and the
    recomputed ones (parts/f3curves_*.json, FIXED complex) use the SAME seeds
    30000+i, i = 0..199, the same C_ell and the same nside, so they are the
    same 200 realisations through the two complexes.  No simulation is
    re-run here; this is arithmetic on committed numbers.
    """
    old = {m["i"]: m for m in load_committed_f3_curves()}
    new = {m["i"]: m for m in _load_parts("f3curves")}
    common = sorted(set(old) & set(new))
    assert len(common) == 200, len(common)
    rows = []
    for i in common:
        assert old[i]["seed"] == new[i]["seed"], (i, old[i]["seed"], new[i]["seed"])
        b0o = np.array(old[i]["b0"], float); b0n = np.array(new[i]["b0"], float)
        b1o = np.array(old[i]["b1"], float); b1n = np.array(new[i]["b1"], float)
        nu = nu_grid()
        rows.append({
            "i": i, "seed": old[i]["seed"],
            "max_abs_delta_b0": float(np.abs(b0n - b0o).max()),
            "max_abs_delta_b1": float(np.abs(b1n - b1o).max()),
            "max_b0_old": float(b0o.max()), "max_b1_old": float(b1o.max()),
            "frac_b0": float(np.abs(b0n - b0o).max() / b0o.max()) if b0o.max() else 0.0,
            "frac_b1": float(np.abs(b1n - b1o).max() / b1o.max()) if b1o.max() else 0.0,
            "peak_nu_shift_b0": float(nu[int(np.argmax(b0n))] - nu[int(np.argmax(b0o))]),
            "peak_nu_shift_b1": float(nu[int(np.argmax(b1n))] - nu[int(np.argmax(b1o))]),
        })
    f0 = np.array([r["frac_b0"] for r in rows])
    f1 = np.array([r["frac_b1"] for r in rows])
    s0 = np.array([abs(r["peak_nu_shift_b0"]) for r in rows])
    s1 = np.array([abs(r["peak_nu_shift_b1"]) for r in rows])
    step = float(nu_grid()[1] - nu_grid()[0])
    out = {
        "gate": "D2 gate 4, all 200 maps (the pre-registered version ran on 1 map, seed 30000)",
        "n_maps": len(rows),
        "same_seeds_verified": True,
        "declared_tolerance": "max |delta| <= 0.15 * peak of the old curve, and the peak nu within "
                              "one grid step (%.2f sigma)" % step,
        "frac_delta_b0": {"max": float(f0.max()), "mean": float(f0.mean()),
                          "p50": float(np.percentile(f0, 50)), "p95": float(np.percentile(f0, 95)),
                          "n_over_0.15": int((f0 > 0.15).sum())},
        "frac_delta_b1": {"max": float(f1.max()), "mean": float(f1.mean()),
                          "p50": float(np.percentile(f1, 50)), "p95": float(np.percentile(f1, 95)),
                          "n_over_0.15": int((f1 > 0.15).sum())},
        "peak_shift_b0": {"max_abs": float(s0.max()), "n_over_one_grid_step": int((s0 > step + 1e-9).sum())},
        "peak_shift_b1": {"max_abs": float(s1.max()), "n_over_one_grid_step": int((s1 > step + 1e-9).sum())},
        "gate_4_all_maps_pass": bool((f0 <= 0.15).all() and (f1 <= 0.15).all()
                                     and (s0 <= step + 1e-9).all() and (s1 <= step + 1e-9).all()),
        "per_map": rows,
        "tier": "X",
    }
    write("gate4full", out)


def step_guard():
    """Run the regression guard BOTH ways and record the outcome.

    A guard that has never been seen to fail is not a guard.  Mode 1 runs the
    tests against the fixed library (must be all green); mode 2 rebinds the
    three fixed entry points to the ORIGINAL implementations via
    tests/conftest.py and must produce failures.
    """
    out = {"gate": "regression guard: the tests must FAIL on the originals and PASS on the fix",
           "tests": "audit/tda_validation/tda_fixed/tests", "tier": "X", "modes": {}}
    for mode, env in (("against_fixed_library", {}),
                      ("against_original_implementations", {"TDA_GUARD_AGAINST_ORIGINAL": "1"})):
        e = dict(os.environ); e.update(env); e["OMP_NUM_THREADS"] = "1"
        t0 = time.time()
        r = subprocess.run([sys.executable, "-m", "pytest", os.path.join(HERE, "tests"),
                            "-q", "--no-header", "-rf"],
                           capture_output=True, text=True, env=e, cwd=WT)
        tail = r.stdout.strip().splitlines()
        failed = sorted({ln.split("::", 1)[1].split(" ")[0]
                         for ln in tail if ln.startswith("FAILED") and "::" in ln})
        summary = [ln for ln in tail if " passed" in ln or " failed" in ln][-1:] or ["(no summary)"]
        out["modes"][mode] = {
            "command": "%s OMP_NUM_THREADS=1 %s -m pytest audit/tda_validation/tda_fixed/tests -q"
                       % (" ".join("%s=%s" % kv for kv in env.items()), sys.executable),
            "exit_code": r.returncode, "summary_line": summary[-1],
            "n_failed": len(failed), "failed_tests": failed, "wall_sec": time.time() - t0}
    a = out["modes"]["against_fixed_library"]
    b = out["modes"]["against_original_implementations"]
    out["guard_pass"] = bool(a["exit_code"] == 0 and b["exit_code"] != 0 and b["n_failed"] > 0)
    out["verdict"] = ("the guard is green on the fixed library and fails on %d tests when the "
                      "original implementations are substituted" % b["n_failed"])
    write("guard", out)


# --------------------------------------------------------------- assemble
def step_assemble():
    parts = {}
    for fn in sorted(glob.glob(os.path.join(PARTS, "*.json"))):
        parts[os.path.basename(fn)[:-5]] = json.load(open(fn))
    f3c = {k: v for k, v in parts.items() if k.startswith("f3curves_") or k.startswith("negctrl_")}
    for k in list(f3c):
        parts.pop(k)

    def head():
        try:
            return subprocess.check_output(["git", "-C", WT, "rev-parse", "HEAD"], text=True).strip()
        except Exception as e:
            return "unavailable: %r" % e

    res = {
        "title": "Before/after for the two TDA defects found by the known-answer suite",
        "generated_by": "audit/tda_validation/tda_fixed/run_validation.py --step assemble",
        "pre_registration": "audit/tda_validation/tda_fixed/fix_expectations.json (commit 1bc7ebc, "
                            "written before any number here)",
        "git_head_at_assemble": head(),
        "tier_of_every_observed_number": "X (numerics)",
        "never_claimed": "no result here is 'proved'; the acceptances are numerical checks against "
                         "pre-declared gates, each with a negative control",
        "code_under_test_sha256": {
            CMB_REL + " (now, with the defect comment blocks)": sha256(os.path.join(WT, CMB_REL)),
            CMB_REL + " (at the suite run, from the committed case JSONs)":
                json.load(open(os.path.join(SUITE, "cases", "F3topo.json")))["code_under_test_sha256"][CMB_REL],
            CW_REL + " (now)": sha256(os.path.join(WT, CW_REL)),
        },
        "fixed_library_sha256": {f: sha256(os.path.join(HERE, f))
                                 for f in sorted(os.listdir(HERE)) if f.endswith(".py")},
        "versions": {"python": sys.version.split()[0], "numpy": np.__version__},
        "parts": parts,
        "f3_curve_chunks": {k: {kk: v.get(kk) for kk in
                                ("start", "count", "nside", "fwhm_deg", "seed_base", "n_edges",
                                 "n_triangles", "euler_char_V_minus_E_plus_F", "wall_sec", "_command")}
                            for k, v in sorted(f3c.items())},
    }
    try:
        import gudhi, healpy
        res["versions"].update({"gudhi": gudhi.__version__, "healpy": healpy.__version__})
    except Exception as e:
        res["versions"]["import_error"] = repr(e)

    # ---- the summary table
    t = {}
    d1o = parts.get("d1old", {}).get("columns", {})
    d1n = parts.get("d1new", {}).get("columns", {})
    for key in ("b0", "b1", "chi"):
        row = {}
        for col, src in (("A_before", d1o.get("A_before")), ("B_minimal", d1o.get("B_minimal")),
                         ("C_full", d1o.get("C_full")),
                         ("D0_fixed_topology_old_stats", d1n.get("D0_fixed_topology_old_stats")),
                         ("D_full_fixed_pipeline", d1n.get("D_full_fixed_pipeline"))):
            if src and key in src.get("keys", {}):
                k = src["keys"][key]
                row[col] = {"ks_p": k["chi2_p"]["ks_p_vs_uniform"],
                            "n_below_0.05": k["chi2_p"]["n_below_0.05"],
                            "df": k["bin_selection"].get("df"),
                            "rank_ks_p": k["rank_p"]["ks_p_vs_uniform"],
                            "rank_n_below_0.05": k["rank_p"]["n_below_0.05"],
                            "gate_pass": k["gate_ks_gt_0.05_and_n05_le_9"]}
        t[key] = row
    res["D1_summary_table"] = t
    res["D1_acceptance"] = {
        "gate": "KS of the 100 chi2 p-values against U(0,1) has p > 0.05 AND #{p<0.05} <= 9, "
                "for b0, b1 and the Euler curve",
        "column_C_pass": {k: t[k].get("C_full", {}).get("gate_pass") for k in t},
        "column_D_pass": {k: t[k].get("D_full_fixed_pipeline", {}).get("gate_pass") for k in t},
        "binomial_context": "P(X >= 10 | p = 0.05, n = 100) = 0.0282; P(X >= 6) = 0.384",
        "known_caveat": "the 100 test p-values share one 100-sim ensemble, so they are correlated; "
                        "KS against uniform assumes independence and is anti-conservative here. "
                        "This is a property of the F3 design (commit 033f017), not of the fix.",
    }
    cal = parts.get("calibcheck", {})
    cal2 = cal.get("chi2_p_branch", {})
    neg = parts.get("negctrl", {})
    guard = parts.get("guard", {})
    res["headline_findings"] = {
        "1_D1_acceptance_met_on_the_declared_case": (
            "column C (fixed statistic on the suite's committed F3 curves, which is the case the "
            "acceptance names) passes for b0, b1 and the Euler curve: KS p = %.4g / %.4g / %.4g, "
            "#{p<0.05} = %s / %s / %s"
            % tuple([t[k]["C_full"]["ks_p"] for k in ("b0", "b1", "chi")]
                    + [t[k]["C_full"]["n_below_0.05"] for k in ("b0", "b1", "chi")])
            if all("C_full" in t[k] for k in t) else "column C missing"),
        "2_the_minimal_enumerated_fix_is_not_sufficient": (
            "column B (dead-bin drop + df = kept bins + eigen conditioning, i.e. exactly the three "
            "steps enumerated in the defect report) leaves b1 at #{p<0.05} = %s > 9, KS p = %.4g. "
            "The atomic-bin rule of column C is what closes it. Reported as its own column rather "
            "than merged into a single 'after'."
            % (t["b1"]["B_minimal"]["n_below_0.05"], t["b1"]["B_minimal"]["ks_p"])
            if "B_minimal" in t.get("b1", {}) else "column B missing"),
        "3_full_fixed_pipeline_b1_still_exceeds_the_count_clause": (
            "column D (fixed topology AND fixed statistic) gives b1 KS p = %.4g (passes) but "
            "#{p<0.05} = %s (the gate allows 9), so the gate as declared is NOT met for b1 on the "
            "recomputed curves. The threshold is not moved: --step calibcheck measures, on exactly "
            "Gaussian data with the same code path, P(#{p<0.05} >= 11) = %.3g for a PERFECTLY "
            "calibrated statistic, because the 100 test maps share one 100-member ensemble."
            % (t["b1"]["D_full_fixed_pipeline"]["ks_p"],
               t["b1"]["D_full_fixed_pipeline"]["n_below_0.05"],
               cal2.get("empirical_P_count_ge_11", float("nan")))
            if "D_full_fixed_pipeline" in t.get("b1", {}) else "column D missing"),
        "4_residual_limitation_of_the_chi2_branch": (
            "even after the fix the chi2 branch is mildly anti-conservative: on exactly Gaussian "
            "data it rejects at %.4f instead of 0.05, because with an ESTIMATED covariance the "
            "exact small-sample reference is Hotelling's T^2 (an F distribution), not chi2, and the "
            "Hartlap factor corrects the mean of the inverse covariance rather than the whole "
            "distribution. The pooled-exchangeable rank p-value has no such residual (%.4f against "
            "an exact %.4f) and is the one to prefer."
            % (cal2.get("pooled_rate_P_chi2_p_below_0.05", float("nan")),
               cal.get("pooled_rate_over_all_10000_test_vectors", float("nan")),
               cal.get("exact_uniform_target_P_rank_p_below_0.05", float("nan")))),
        "5_the_gate_itself_has_a_false_failure_rate": (
            "the declared '#{p<0.05} <= 9 of 100' clause fails with probability %.3g per curve for a "
            "perfectly calibrated statistic under the F3 design (variance inflation %.2f from the "
            "shared ensemble). This is a property of the F3 design (commit 033f017), not of the fix."
            % (cal2.get("empirical_P_count_ge_10", float("nan")),
               cal2.get("variance_inflation_factor", float("nan")))),
        "6_negative_controls_do_reject": (
            "20 maps drawn at beam FWHM 6 deg against the 3 deg ensemble: #{chi2 p < 0.05} = %s/20, "
            "%s/20, %s/20 for b0, b1, Euler (median p %.3g, %.3g, %.3g). The statistic can fail."
            % tuple([neg.get("keys", {}).get(k, {}).get("n_chi2_p_below_0.05") for k in ("b0", "b1", "chi")]
                    + [neg.get("keys", {}).get(k, {}).get("median_chi2_p", float("nan")) for k in ("b0", "b1", "chi")])
            if neg.get("keys") else "negative control missing"),
        "7_regression_guard_fails_on_the_originals": guard.get("verdict", "guard not run"),
        "8_gate_4_is_n_equals_1_as_pre_registered_and_the_200_map_extension_is_not_uniform": (
            "the declared gate 4 compares the old and fixed complexes on ONE map (seed 30000) and "
            "passes. Extending it to all 200 F3 realisations (same seeds through both complexes, "
            "arithmetic on committed curves): the b0 difference exceeds the 15%% tolerance for "
            "%s of 200 maps (max %.3f, median %.3f) and b1 for %s of 200 (max %.3f); the peak nu "
            "moves by more than one grid step for %s maps in b0 and %s in b1 (max %.1f sigma). "
            "The 1-skeleton genuinely changed - one diagonal per quad instead of two - so exact "
            "agreement was never expected; and the peak-location clause is brittle because argmax "
            "on a nearly flat maximum jumps for small changes. Reported rather than absorbed into "
            "the n=1 gate."
            % (parts["gate4full"]["frac_delta_b0"]["n_over_0.15"],
               parts["gate4full"]["frac_delta_b0"]["max"],
               parts["gate4full"]["frac_delta_b0"]["p50"],
               parts["gate4full"]["frac_delta_b1"]["n_over_0.15"],
               parts["gate4full"]["frac_delta_b1"]["max"],
               parts["gate4full"]["peak_shift_b0"]["n_over_one_grid_step"],
               parts["gate4full"]["peak_shift_b1"]["n_over_one_grid_step"],
               parts["gate4full"]["peak_shift_b0"]["max_abs"])
            if "gate4full" in parts else "gate4full not run"),
    }
    res["regression_guard"] = guard
    res["e5_impact"] = parts.get("e5", {}).get("statistics", {})
    res["calibration_context"] = {k: v for k, v in cal.items()
                                  if k not in ("counts", "_command") and not k.startswith("_")}

    top = parts.get("topology", {})
    res["D2_acceptance"] = {
        "gate_1_combinatorial": {n: v.get("gate_1_combinatorial_pass")
                                 for n, v in top.get("levels", {}).items()},
        "gate_2_betti_1_0_1": {n: v.get("gate_2_betti_1_0_1_pass")
                               for n, v in top.get("levels", {}).items() if "gudhi_betti" in v},
        "gate_3_disk_mask": top.get("gate_3_disk_mask", {}).get("gate_3_pass"),
        "gate_4_smooth_field": parts.get("smooth", {}).get("gate_4_pass"),
        "gate_4_note": "the pre-registered gate 4 is n = 1 (seed 30000) and it passes: "
                       "max |delta b0| = 8 against a tolerance of 15.9, max |delta b1| = 6 "
                       "against 16.2, both peaks unmoved. The post-hoc 200-map extension "
                       "(step gate4full, arithmetic on committed curves, no simulation) is "
                       "stronger and does NOT hold for every realisation - see "
                       "gate_4_all_200_maps.",
        "gate_4_all_200_maps": {k: parts.get("gate4full", {}).get(k) for k in
                                ("n_maps", "frac_delta_b0", "frac_delta_b1", "peak_shift_b0",
                                 "peak_shift_b1", "gate_4_all_maps_pass")},
        "gate_5_f1_untouched": parts.get("f1", {}).get("gate_5_pass"),
        "before": {n: v.get("OLD_build_topology") for n, v in top.get("levels", {}).items()
                   if "OLD_build_topology" in v},
    }
    json.dump(res, open(os.path.join(HERE, "validation_results.json"), "w"), indent=1, default=float)
    print(json.dumps({"D1_C": res["D1_acceptance"]["column_C_pass"],
                      "D1_D": res["D1_acceptance"]["column_D_pass"],
                      "D2": {k: v for k, v in res["D2_acceptance"].items() if k != "before"}}, indent=1))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--step", required=True)
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--count", type=int, default=50)
    ARGS = ap.parse_args()
    if ARGS.step == "f3curves":
        step_f3curves(ARGS.start, ARGS.count)
    elif ARGS.step == "negctrl_curves":
        step_f3curves(ARGS.start, ARGS.count, fwhm=NEG_FWHM_DEG, seed_base=NEG_SEED_BASE, tag="negctrl")
    else:
        {"topology": step_topology, "smooth": step_smooth, "f1": step_f1, "p2": step_p2,
         "d1old": step_d1old, "d1new": step_d1new, "negctrl": step_negctrl, "calibcheck": step_calibcheck, "guard": step_guard, "gate4full": step_gate4full,
         "e5": step_e5, "assemble": step_assemble}[ARGS.step]()
