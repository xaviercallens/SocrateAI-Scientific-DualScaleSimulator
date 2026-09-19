"""Test 2: circular bacterial chromosome from Hi-C (pre-stated in expectations.json).

Usage (resumable; re-invoke with the same arguments until it prints a verdict):
  python test2_hic.py --n-null 500 --n-rips-null 50 --budget-sec 540 [--cases caulo_full,caulo_cut_ter,...]
  (null sizes reduced from the pre-stated 1000/200 for compute; recorded in expectations_addendum_ecoli.json)
  --cases runs only the listed cases (to fill the caches in parallel processes); the final
  invocation without --cases assembles all cases from the caches and writes the JSON.
Test 2b (E. coli, expectations_addendum_ecoli.json) is evaluated in the same script.
Seeds: linear-null seeds 0..n_null-1 (alpha path), 0..n_rips_null-1 (Rips path).

Deviation from expectations.json, recorded in the output: the authors' iteratively
corrected Caulobacter matrix has its main diagonal AND first off-diagonal set to 0
for every bin (masked by processing, not measured zeros). Treating those as zero
contacts would make adjacent bins maximally distant, which is not the intent of the
pre-stated 'zero contacts -> half the smallest positive contact' rule (written for
sporadic zeros). Masked first-offset entries (linear offset 1) are imputed with the
mean of the offset-2 diagonal E(2); the same rule is applied to every null matrix.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import argparse
import glob
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tda_common import (DATA_ROOT, alpha_h1_summary, rips_h1_summary, pval_ge, classical_mds, dump,
                        pipeline_sha256, resumable_seed_loop, Budget)

HERE = os.path.dirname(os.path.abspath(__file__))
ALPHA_EXP = 1.0 / 3.0


def load_caulo():
    p = glob.glob(os.path.join(DATA_ROOT, "caulobacter", "GSM1120445_*after_normalization.txt.gz"))[0]
    return np.loadtxt(p), p


def load_gm():
    p = os.path.join(DATA_ROOT, "gm12878", "gm12878_insitu_combined_chr1_145000000_249250621_250kb_KR_observed.npy")
    M = np.load(p)
    M = np.nan_to_num(M, nan=0.0)  # absent sparse records = zero contacts; KR-NaN bins give all-zero rows
    return M, p


def load_ecoli():
    import pandas as pd
    p = glob.glob(os.path.join(DATA_ROOT, "ecoli", "GSM2870407_*.txt.gz"))[0]
    df = pd.read_csv(p, sep="\t", index_col=0)
    df = df.loc[:, [c for c in df.columns if not str(c).startswith("Unnamed")]]  # trailing tab -> empty column
    M = df.values.astype(float)
    assert M.shape[0] == M.shape[1] and np.isfinite(M).all(), M.shape
    n5 = M.shape[0]
    nb = (n5 + 1) // 2
    B = np.zeros((nb, nb))
    idx = np.arange(n5) // 2
    np.add.at(B, (idx[:, None], idx[None, :]), M)
    np.fill_diagonal(B, 0.0)
    rs = B.sum(1)
    med = np.median(rs); mad = 1.4826 * np.median(np.abs(rs - med))
    keep = rs >= med - 3 * mad
    B = B[np.ix_(keep, keep)]
    it_done, dev = 0, None
    for it in range(200):
        r = B.sum(1); r = r / r.mean()
        dev = float(np.max(np.abs(r - 1)))
        if dev < 1e-6:
            break
        B = B / np.outer(r, r)
        it_done = it + 1
    info = {"raw_shape": list(M.shape), "n_bins_10kb": int(nb), "n_bins_dropped_low_coverage": int((~keep).sum()),
            "ice_iterations": it_done, "ice_final_max_dev": dev}
    return B, np.where(keep)[0], info, p


def clean(C, masked_offset1):
    """Symmetrise-check, drop zero-sum bins, impute masked offset-1 (if flagged),
    replace remaining off-diagonal zeros by half the smallest positive contact."""
    info = {"max_abs_asym": float(np.abs(C - C.T).max())}
    keep = C.sum(1) > 0
    info["n_bins_in"] = int(C.shape[0])
    info["n_zero_sum_bins_dropped"] = int((~keep).sum())
    C = C[np.ix_(keep, keep)].copy()
    n = C.shape[0]
    iu = np.arange(n - 1)
    if masked_offset1:
        e2 = float(np.mean(np.diagonal(C, 2)))
        C[iu, iu + 1] = e2
        C[iu + 1, iu] = e2
        info["offset1_imputed_with_E2"] = e2
    off = ~np.eye(n, dtype=bool)
    pos = C[off & (C > 0)]
    fill = 0.5 * float(pos.min())
    nz = off & (C <= 0)
    info["n_offdiag_zero_replaced"] = int(nz.sum() // 2)
    C[nz] = fill
    np.fill_diagonal(C, 0.0)
    return C, info


def to_dist(C, a=ALPHA_EXP):
    D = np.zeros_like(C)
    off = ~np.eye(C.shape[0], dtype=bool)
    D[off] = C[off] ** (-a)
    return D


def analyse(C, label, with_rips=True):
    D = to_dist(C)
    X3, w = classical_mds(D, 3)
    s = alpha_h1_summary(X3, label)
    s["mds_top3_eig"] = w[:3].tolist()
    s["mds_frac_negative_eig_mass"] = float(np.sum(np.abs(w[w < 0])) / np.sum(np.abs(w)))
    out = {"alpha_path_mds3": s}
    if with_rips:
        r = rips_h1_summary(distance_matrix=D)
        r.pop("all_h1_persistence_sorted")
        out["rips_full_matrix"] = r
    # sensitivity: alpha exponent 1 on the alpha path
    X3b, _ = classical_mds(to_dist(C, 1.0), 3)
    sb = alpha_h1_summary(X3b, label + "_a1")
    out["alpha_path_mds3_exponent1_sensitivity"] = {k: sb[k] for k in ("P1", "P2", "P1_over_P2", "S", "top_h1_bars_birth_death")}
    return out


def linear_null_factory(C, circular):
    """E(s) for s <= N/2 from the data (linear and circular offsets coincide there),
    power law fitted on N/8 <= s <= N/2 extrapolated beyond; residuals relative to
    E(circular-or-linear offset) permuted across pairs."""
    n = C.shape[0]
    I, J = np.triu_indices(n, 2)          # offset-1 is imputed, not data; excluded from residuals
    s_lin = J - I
    s_eff = np.minimum(s_lin, n - s_lin) if circular else s_lin
    half = n // 2
    E = np.zeros(n)
    for s in range(1, half + 1):
        E[s] = np.mean(np.diagonal(C, s))
    ss = np.arange(max(2, n // 8), half + 1)
    slope, icpt = np.polyfit(np.log(ss), np.log(E[ss]), 1)
    for s in range(half + 1, n):
        E[s] = np.exp(icpt + slope * np.log(s))
    resid = np.log(C[I, J]) - np.log(E[s_eff])
    e2 = float(np.mean(np.diagonal(C, 2)))

    def make(seed):
        rng = np.random.RandomState(seed)
        Cn = np.zeros_like(C)
        v = E[s_lin] * np.exp(resid[rng.permutation(resid.size)])
        Cn[I, J] = v
        Cn[J, I] = v
        iu = np.arange(n - 1)
        Cn[iu, iu + 1] = e2
        Cn[iu + 1, iu] = e2
        return Cn
    return make, {"powerlaw_slope_fit_N8_to_N2": float(slope), "E_offsets_2_5_10_half": [float(E[2]), float(E[5]), float(E[10]), float(E[half])]}


def run_case(C, label, circular, a, cache_dir):
    obs = analyse(C, label)
    make, nullinfo = linear_null_factory(C, circular)

    def f_alpha(sd):
        X3, _ = classical_mds(to_dist(make(sd)), 3)
        s = alpha_h1_summary(X3)
        return [s["S"], s["P1_over_P2"]]

    def f_rips(sd):
        r = rips_h1_summary(distance_matrix=to_dist(make(sd)))
        return [r["S_rips"], r["P1_over_P2"]]

    na = np.array(resumable_seed_loop(os.path.join(cache_dir, f"test2_{label}_alphanull_n{a.n_null}.json"), a.n_null, f_alpha, a.budget_sec))
    nr = np.array(resumable_seed_loop(os.path.join(cache_dir, f"test2_{label}_ripsnull_n{a.n_rips_null}.json"), a.n_rips_null, f_rips, a.budget_sec))
    sA = obs["alpha_path_mds3"]
    sR = obs["rips_full_matrix"]
    obs["linear_null"] = dict(nullinfo, **{
        "n_null_alpha": a.n_null, "p_S_alpha": pval_ge(na[:, 0], sA["S"]),
        "null_S_alpha_q50_95_99_max": np.quantile(na[:, 0], [0.5, 0.95, 0.99, 1]).tolist(),
        "null_frac_alpha_P1_over_P2_ge_2": float(np.mean(na[:, 1] >= 2)),
        "n_null_rips": a.n_rips_null, "p_S_rips": pval_ge(nr[:, 0], sR["S_rips"]),
        "null_S_rips_q50_95_max": np.quantile(nr[:, 0], [0.5, 0.95, 1]).tolist(),
        "null_frac_rips_P1_over_P2_ge_2": float(np.mean(nr[:, 1] >= 2))})
    obs["dominant_and_significant_alpha"] = bool(sA["P1_over_P2"] >= 2.0 and obs["linear_null"]["p_S_alpha"] <= 0.01)
    obs["dominant_and_significant_rips"] = bool(sR["P1_over_P2"] >= 2.0 and obs["linear_null"]["p_S_rips"] <= 0.01)
    return obs


ALL_CASES = ["caulo_full", "caulo_cut_ter", "caulo_cut_ori", "gm12878_chr1q", "ecoli_full", "ecoli_cut_ter"]


def build_cases():
    Craw, pC = load_caulo()
    C, info_c = clean(Craw, masked_offset1=True)
    n = C.shape[0]
    assert n == 405, n  # bin indices of the cut controls assume no dropped bins
    order_ter = np.r_[np.arange(222, n), np.arange(0, 182)]   # ter cut: bins 182-221 removed, chain 222..404,0..181
    order_ori = np.arange(20, 385)                              # ori cut: bins 0-19 and 385-404 removed
    Graw, pG = load_gm()
    G, info_g = clean(Graw, masked_offset1=False)
    E, kept, info_e, pE = load_ecoli()
    E, info_e2 = clean(E, masked_offset1=False)
    info_e.update(info_e2)
    chain = [k for k, b in enumerate(kept) if b >= 182] + [k for k, b in enumerate(kept) if b < 136]
    info_e["ter_cut_original_bins_removed"] = "136-181"
    cases = {
        "caulo_full": (C, True), "caulo_cut_ter": (C[np.ix_(order_ter, order_ter)], False),
        "caulo_cut_ori": (C[np.ix_(order_ori, order_ori)], False), "gm12878_chr1q": (G, False),
        "ecoli_full": (E, True), "ecoli_cut_ter": (E[np.ix_(chain, chain)], False)}
    meta = {"caulobacter": {"file": pC, "clean_info": info_c}, "gm12878_chr1q": {"file": pG, "clean_info": info_g},
            "ecoli": {"file": pE, "clean_info": info_e}}
    return cases, meta


def crit_block(full, cut, ext, path):
    k = "alpha_path_mds3" if path == "alpha" else "rips_full_matrix"
    pk = "p_S_alpha" if path == "alpha" else "p_S_rips"
    dk = "dominant_and_significant_" + path
    return {"i_P1_over_P2_ge_2": bool(full[k]["P1_over_P2"] >= 2.0),
            "ii_S_p_le_0.01_vs_linear_null": bool(full["linear_null"][pk] <= 0.01),
            "iii_ter_cut_not_dominant_and_significant": not cut[dk],
            "iv_gm12878_linear_not_dominant_and_significant": not ext[dk]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-null", type=int, default=500)
    ap.add_argument("--n-rips-null", type=int, default=50)
    ap.add_argument("--budget-sec", type=float, default=540.0)
    ap.add_argument("--cases", default=",".join(ALL_CASES))
    a = ap.parse_args()
    cache_dir = os.path.join(DATA_ROOT, "cache")
    os.makedirs(cache_dir, exist_ok=True)
    cases, meta = build_cases()
    sel = a.cases.split(",")
    out = {}
    for name in sel:
        Cm, circ = cases[name]
        out[name] = run_case(Cm, name, circ, a, cache_dir)
        print(name, "done")
    if set(sel) != set(ALL_CASES):
        return
    res = {"args": vars(a), "pipeline_file_sha256": pipeline_sha256(),
           "distance": "d_ij = C_ij^(-1/3); alpha path = classical MDS to 3-D -> alpha_persistence",
           "deviations_from_expectations": ["Caulobacter masked offset-1 diagonal imputed with E(2) (see module docstring)",
                                            "null sizes 500 (alpha) / 50 (Rips) instead of 1000 / 200 (compute); recorded in expectations_addendum_ecoli.json before the E. coli computation"],
           "inputs": meta, "cases": out}
    res["test2_caulobacter"] = {"criteria_alpha_binding": crit_block(out["caulo_full"], out["caulo_cut_ter"], out["gm12878_chr1q"], "alpha"),
                                "criteria_rips_crosscheck": crit_block(out["caulo_full"], out["caulo_cut_ter"], out["gm12878_chr1q"], "rips")}
    res["test2b_ecoli"] = {"criteria_alpha_binding": crit_block(out["ecoli_full"], out["ecoli_cut_ter"], out["gm12878_chr1q"], "alpha"),
                           "criteria_rips_crosscheck": crit_block(out["ecoli_full"], out["ecoli_cut_ter"], out["gm12878_chr1q"], "rips")}
    for t in ("test2_caulobacter", "test2b_ecoli"):
        res[t]["verdict"] = "PASS" if all(res[t]["criteria_alpha_binding"].values()) else "FAIL"
        res[t]["verdict_rips_crosscheck"] = "PASS" if all(res[t]["criteria_rips_crosscheck"].values()) else "FAIL"
        print(t, res[t])
    dump(res, os.path.join(HERE, "test2_hic_results.json"))
    for k, v in out.items():
        s = v["alpha_path_mds3"]; r = v["rips_full_matrix"]; ln = v["linear_null"]
        print(k, "alpha top", np.round(s["top_h1_bars_birth_death"][:2], 4).tolist(), "P1/P2 %.2f S %.3f p %.4f" % (s["P1_over_P2"], s["S"], ln["p_S_alpha"]),
              "| rips top", np.round(r["top_h1_bars_birth_death"][:2], 4).tolist(), "P1/P2 %.2f p %.4f" % (r["P1_over_P2"], ln["p_S_rips"]))


if __name__ == "__main__":
    try:
        main()
    except Budget as e:
        print("INCOMPLETE (budget reached, re-invoke with the same arguments to resume):", e)
        sys.exit(3)
