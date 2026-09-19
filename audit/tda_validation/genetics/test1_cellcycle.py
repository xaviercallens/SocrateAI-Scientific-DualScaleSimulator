"""Test 1: cell-cycle loop in single-cell RNA-seq, pre-stated in expectations.json.

Usage:
  python test1_cellcycle.py --dataset u2os  --n-null 1000 --n-rips-null 100 --n-neg 20 --n-perm 10000
  python test1_cellcycle.py --dataset mesc  --n-null 1000 --n-rips-null 100 --n-neg 20 --n-perm 10000 --n-tri 200
Seeds: null seeds 0..n_null-1; negative-control draws 0..n_neg-1; label
permutation seed 0; Rips subsample seed 0; within-phase control seeds 0..n_tri-1.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")  # machine is oversubscribed; 1 BLAS thread is ~10x faster here

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tda_common import (DATA_ROOT, alpha_h1_summary, rips_h1_summary, pval_ge, circ_corr,
                        pca3, mad_low, dump, pipeline_sha256, sha256, resumable_seed_loop, Budget)

HERE = os.path.dirname(os.path.abspath(__file__))


def load_u2os():
    d = os.path.join(DATA_ROOT, "u2os")
    counts = pd.read_csv(os.path.join(d, "GSE146773_Counts.csv.gz"), index_col=0)
    fucci = pd.read_csv(os.path.join(d, "GSE146773_fucci_coords.csv.gz")).set_index("cell")
    common = counts.index.intersection(fucci.index)
    fucci = fucci.loc[common]
    fucci = fucci[np.isfinite(fucci["polar_coord_phi"].astype(float))]
    counts = counts.loc[fucci.index]
    genes = [g for g in counts.columns if not g.startswith("ERCC")]
    ercc = [g for g in counts.columns if g.startswith("ERCC")]
    X = counts[genes].values.astype(float)
    info = {"n_cells_counts": int(len(counts.index)), "n_cells_fucci_matched": int(len(fucci)),
            "n_genes": len(genes), "n_ercc_cols": len(ercc)}
    lab = {"phi": fucci["polar_coord_phi"].values.astype(float),
           "time": fucci["fucci_time_hrs"].values.astype(float)}
    cc = pd.read_csv(os.path.join(DATA_ROOT, "genelists", "Homo_sapiens.csv"))
    return X, np.array(genes), None, lab, cc, info


def load_mesc():
    d = os.path.join(DATA_ROOT, "buettner")
    mats, labs = [], []
    for ph in ("G1", "S", "G2M"):
        t = pd.read_csv(os.path.join(d, f"{ph}_singlecells_counts.txt"), sep="\t")
        ids = t["EnsemblGeneID"].values
        m = t.iloc[:, 4:].values.astype(float).T
        mats.append(m)
        labs += [ph] * m.shape[0]
    X_all = np.vstack(mats)
    is_ercc = np.array([str(g).startswith("ERCC") for g in ids])
    X = X_all[:, ~is_ercc]
    E = X_all[:, is_ercc]
    info = {"n_cells": int(X.shape[0]), "n_genes": int((~is_ercc).sum()), "n_ercc": int(is_ercc.sum())}
    cc = pd.read_csv(os.path.join(DATA_ROOT, "genelists", "Mus_musculus.csv"))
    return X, ids[~is_ercc], E, {"phase": np.array(labs)}, cc, info


def qc(X, E):
    lib = np.log10(X.sum(1) + 1)
    det = (X > 0).sum(1)
    lo_lib, _ = mad_low(lib)
    lo_det, _ = mad_low(det.astype(float))
    keep = (lib >= lo_lib) & (det >= lo_det)
    thr = {"log10_lib_min": lo_lib, "detected_min": lo_det}
    if E is not None:
        frac = E.sum(1) / (E.sum(1) + X.sum(1))
        _, hi = mad_low(frac)
        keep &= frac <= hi
        thr["ercc_frac_max"] = hi
    return keep, thr


def zscore(M):
    sd = M.std(0)
    sd[sd == 0] = 1.0
    return (M - M.mean(0)) / sd


def loop_stats(Z):
    X3, varfrac = pca3(Z)
    s = alpha_h1_summary(X3)
    s["pca_var_frac_top3"] = varfrac[:3].tolist()
    return s, X3


def gene_perm(Z, rng):
    Zp = Z.copy()
    for j in range(Z.shape[1]):
        Zp[:, j] = Z[rng.permutation(Z.shape[0]), j]
    return Zp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", choices=["u2os", "mesc"], required=True)
    ap.add_argument("--n-null", type=int, default=1000)
    ap.add_argument("--n-rips-null", type=int, default=100)
    ap.add_argument("--n-neg", type=int, default=20)
    ap.add_argument("--n-perm", type=int, default=10000)
    ap.add_argument("--n-tri", type=int, default=200)
    ap.add_argument("--budget-sec", type=float, default=540.0,
                    help="wall-clock budget; null loops are cached and the run is resumed by re-invoking")
    a = ap.parse_args()
    cache_dir = os.path.join(DATA_ROOT, "cache")
    os.makedirs(cache_dir, exist_ok=True)
    tag = f"test1_{a.dataset}"

    def loop(name, n, fn):
        return resumable_seed_loop(os.path.join(cache_dir, f"{tag}_{name}_n{n}.json"), n, fn, a.budget_sec)

    if a.dataset == "u2os":
        X, genes, E, lab, cc, info = load_u2os()
    else:
        X, genes, E, lab, cc, info = load_mesc()
    keep, thr = qc(X, E)
    X = X[keep]
    lab = {k: v[keep] for k, v in lab.items()}
    info.update({"qc_thresholds": thr, "n_cells_after_qc": int(keep.sum())})

    L = np.log2(X / X.sum(1, keepdims=True) * 1e6 + 1)
    det_frac = (X > 0).mean(0)
    cc_ids = set(cc["geneID"].astype(str))
    is_cc = np.array([str(g) in cc_ids for g in genes])
    cc_sel = is_cc & (det_frac >= 0.05)
    info.update({"cc_list_size": len(cc_ids), "cc_genes_in_matrix": int(is_cc.sum()),
                 "cc_genes_used": int(cc_sel.sum())})
    Z = zscore(L[:, cc_sel])
    obs, X3 = loop_stats(Z)

    # label angle
    if a.dataset == "u2os":
        phi = lab["phi"]
        phi_t = 2 * np.pi * lab["time"] / np.max(lab["time"])
    else:
        amap = {"G1": 0.0, "S": 2 * np.pi / 3, "G2M": 4 * np.pi / 3}
        phi = np.array([amap[p] for p in lab["phase"]])
        phi_t = None
    theta = np.arctan2(X3[:, 1], X3[:, 0])
    rho = circ_corr(theta, phi)
    rng = np.random.RandomState(0)
    perm_rho = np.array([abs(circ_corr(theta, phi[rng.permutation(len(phi))])) for _ in range(a.n_perm)])
    pos = {"rho_cc": rho, "abs_rho_cc": abs(rho), "p_perm": pval_ge(perm_rho, abs(rho)), "n_perm": a.n_perm}
    if phi_t is not None:
        pos["rho_cc_fucci_time"] = circ_corr(theta, phi_t)
    if a.dataset == "mesc":
        pos["circular_mean_theta_by_phase_deg"] = {
            p: float(np.degrees(np.angle(np.mean(np.exp(1j * theta[lab["phase"] == p]))))) for p in ("G1", "S", "G2M")}

    # gene-wise permutation null (alpha path)
    def _null(sd):
        s, _ = loop_stats(gene_perm(Z, np.random.RandomState(sd)))
        return [s["S"], s["P1_over_P2"]]
    nv = np.array(loop("null", a.n_null, _null))
    nullS, nullR = nv[:, 0], nv[:, 1]
    null = {"n_null": a.n_null, "p_S": pval_ge(nullS, obs["S"]),
            "null_S_quantiles_50_95_99_max": np.quantile(nullS, [0.5, 0.95, 0.99, 1.0]).tolist(),
            "null_frac_P1_over_P2_ge_2": float(np.mean(np.array(nullR) >= 2.0))}

    # negative control: random expressed non-cc genes
    pool = np.where((~is_cc) & (det_frac >= 0.20))[0]
    neg = []
    for sd in range(a.n_neg):
        idx = np.random.RandomState(sd).choice(pool, size=int(cc_sel.sum()), replace=False)
        Zn = zscore(L[:, idx])
        s, Xn = loop_stats(Zn)
        th_n = np.arctan2(Xn[:, 1], Xn[:, 0])
        rec = {"seed": sd, "P1_over_P2": s["P1_over_P2"], "S": s["S"], "abs_rho_cc": abs(circ_corr(th_n, phi)),
               "top_h1_bars_birth_death": s["top_h1_bars_birth_death"][:3]}
        if sd == 0:
            nS = np.array(loop("neg0null", a.n_null, lambda k: loop_stats(gene_perm(Zn, np.random.RandomState(k)))[0]["S"]))
            rec["p_S_vs_own_null"] = pval_ge(nS, s["S"])
            rec["own_null_n"] = a.n_null
            primary = rec
        neg.append(rec)
    negc = {"pool_size": int(pool.size), "n_draws": a.n_neg, "primary_draw0": primary,
            "draws": neg,
            "frac_draws_P1_over_P2_ge_2": float(np.mean([r["P1_over_P2"] >= 2 for r in neg])),
            "median_abs_rho_cc": float(np.median([r["abs_rho_cc"] for r in neg]))}
    neg_dominant_and_sig = bool(primary["P1_over_P2"] >= 2.0 and primary["p_S_vs_own_null"] <= 0.01)

    # Rips cross-check on full-dimensional z-scored cc genes
    rs = np.random.RandomState(0)
    sub = rs.choice(Z.shape[0], size=min(400, Z.shape[0]), replace=False)
    rips = rips_h1_summary(points=Z[sub])
    rnull = loop("ripsnull", a.n_rips_null, lambda sd: rips_h1_summary(points=gene_perm(Z, np.random.RandomState(sd))[sub])["S_rips"])
    rips_out = {k: v for k, v in rips.items() if k != "all_h1_persistence_sorted"}
    rips_out.update({"n_sub": int(sub.size), "sub_seed": 0, "n_null": a.n_rips_null,
                     "p_S_rips": pval_ge(rnull, rips["S_rips"]),
                     "null_S_rips_quantiles_50_95_max": np.quantile(rnull, [0.5, 0.95, 1.0]).tolist()})
    # also Rips on the 3-D PCA embedding of the same subsample (to separate embedding from complex type)
    rips3 = rips_h1_summary(points=X3[sub])
    rips_out["rips_on_pca3_same_subsample"] = {k: v for k, v in rips3.items() if k != "all_h1_persistence_sorted"}

    out = {"dataset": a.dataset, "args": vars(a), "pipeline_file_sha256": pipeline_sha256(),
           "info": info, "observed_alpha_pca3": obs, "null_gene_permutation": null,
           "positive_check_circular_correlation": pos, "negative_control_random_noncc_genes": negc,
           "rips_crosscheck_full_dim": rips_out}

    if a.dataset == "mesc":
        # within-phase gene-wise permutation control (keeps the phase triangle)
        def _tri(sd):
            r = np.random.RandomState(sd)
            Zt = Z.copy()
            for p in ("G1", "S", "G2M"):
                ix = np.where(lab["phase"] == p)[0]
                for j in range(Z.shape[1]):
                    Zt[ix, j] = Z[ix[r.permutation(ix.size)], j]
            return loop_stats(Zt)[0]["S"]
        tri = loop("tri", a.n_tri, _tri)
        out["continuity_diagnostic"] = {
            "dominant_birth_over_death": obs["dominant_birth_over_death"],
            "within_phase_perm_n": a.n_tri, "p_tri": pval_ge(tri, obs["S"]),
            "within_phase_S_quantiles_50_95_max": np.quantile(tri, [0.5, 0.95, 1.0]).tolist(),
            "continuous_loop_rule": "birth/death <= 0.5 and p_tri <= 0.05",
            "continuous": bool(obs["dominant_birth_over_death"] is not None and obs["dominant_birth_over_death"] <= 0.5
                               and pval_ge(tri, obs["S"]) <= 0.05)}

    crit = {"i_dominant_P1_over_P2_ge_2": bool(obs["P1_over_P2"] >= 2.0),
            "ii_S_significant_p_le_0.01": bool(null["p_S"] <= 0.01),
            "iii_negative_control_not_dominant_and_significant": not neg_dominant_and_sig,
            "iv_abs_rho_ge_0.3_and_p_le_0.01": bool(abs(rho) >= 0.3 and pos["p_perm"] <= 0.01)}
    out["criteria"] = crit
    out["verdict"] = "PASS" if all(crit.values()) else "FAIL"
    dump(out, os.path.join(HERE, f"test1_{a.dataset}_results.json"))
    print(out["verdict"], crit)
    print("obs", obs["top_h1_bars_birth_death"][:3], "P1/P2", obs["P1_over_P2"], "S", obs["S"], "p", null["p_S"])
    print("pos", pos)
    print("neg0", primary)
    print("rips", rips_out["top_h1_bars_birth_death"][:3], rips_out["P1_over_P2"], rips_out["p_S_rips"])
    if "continuity_diagnostic" in out:
        print("cont", out["continuity_diagnostic"])


if __name__ == "__main__":
    try:
        main()
    except Budget as e:
        print("INCOMPLETE (budget reached, re-invoke with the same arguments to resume):", e)
        sys.exit(3)
