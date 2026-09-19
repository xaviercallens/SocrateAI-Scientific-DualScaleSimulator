"""POST-HOC diagnostics (NOT pre-stated; do not change any verdict).
Command: python posthoc_diagnostics.py
(A) U2OS/mESC: is the cell-cycle loop hidden by the embedding (3-D alpha filling a 2-D loop through the
    noise dimension) or by density? Alpha on PC1-PC2 (points padded with z=0 so the pipeline's 3-D function
    is still the one used), PCs whitened, subsamples of 250 cells (seeds 0..9), Rips on PC1-PC2.
(B) Test-2 machinery known-answer control: synthetic noiseless circular / linear contact matrices
    C_ij = s^-1 with s = min(|i-j|, N-|i-j|) (ring) or |i-j| (chain), N = 405, through the SAME
    to_dist -> classical_mds -> alpha path and the Rips path. Also the MDS negative-eigenvalue mass and the
    genomic location of the top Rips H1 birth edge on the real Caulobacter matrix."""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys
import numpy as np
import gudhi
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tda_common import alpha_h1_summary, rips_h1_summary, pca3, dump
import test1_cellcycle as T1
import test2_hic as T2

HERE = os.path.dirname(os.path.abspath(__file__))
KEYS = ("P1", "P2", "P1_over_P2", "S", "dominant_birth_over_death", "top_h1_bars_birth_death")


def slim(d, keys=KEYS):
    return {k: d[k] for k in keys if k in d}


def cc_matrix(ds):
    X, genes, E, lab, cc, info = T1.load_u2os() if ds == "u2os" else T1.load_mesc()
    keep, _ = T1.qc(X, E)
    X = X[keep]
    L = np.log2(X / X.sum(1, keepdims=True) * 1e6 + 1)
    det = (X > 0).mean(0)
    ids = set(cc["geneID"].astype(str))
    sel = np.array([str(g) in ids for g in genes]) & (det >= 0.05)
    return T1.zscore(L[:, sel])


def diag_cells(ds):
    Z = cc_matrix(ds)
    X3, _ = pca3(Z)
    out = {"n_cells": int(Z.shape[0])}
    out["alpha_pca3_as_prestated"] = slim(alpha_h1_summary(X3))
    X2 = np.column_stack([X3[:, 0], X3[:, 1], np.zeros(len(X3))])
    out["alpha_pc12_z0"] = slim(alpha_h1_summary(X2))
    W = X3 / X3.std(0)
    out["alpha_pca3_whitened"] = slim(alpha_h1_summary(W))
    W2 = np.column_stack([W[:, 0], W[:, 1], np.zeros(len(W))])
    out["alpha_pc12_whitened_z0"] = slim(alpha_h1_summary(W2))
    subs = []
    for sd in range(10):
        ix = np.random.RandomState(sd).choice(len(X3), size=min(250, len(X3)), replace=False)
        s3 = alpha_h1_summary(X3[ix]); s2 = alpha_h1_summary(W2[ix])
        subs.append({"seed": sd, "pca3_P1_over_P2": s3["P1_over_P2"], "pca3_b_over_d": s3["dominant_birth_over_death"],
                     "pc12w_P1_over_P2": s2["P1_over_P2"], "pc12w_b_over_d": s2["dominant_birth_over_death"]})
    out["subsample_250_seeds0_9"] = subs
    r = rips_h1_summary(points=X3[:, :2]); r.pop("all_h1_persistence_sorted")
    out["rips_pc12"] = r
    # gene-permutation null for the whitened PC1-2 alpha statistic, 200 seeds (post hoc)
    nullR, nullS = [], []
    for sd in range(200):
        Xn, _ = pca3(T1.gene_perm(Z, np.random.RandomState(sd)))
        Wn = Xn / Xn.std(0)
        s = alpha_h1_summary(np.column_stack([Wn[:, 0], Wn[:, 1], np.zeros(len(Wn))]))
        nullR.append(s["P1_over_P2"]); nullS.append(s["S"])
    o = out["alpha_pc12_whitened_z0"]
    out["alpha_pc12_whitened_z0_null200"] = {
        "p_S": float((1 + np.sum(np.array(nullS) >= o["S"])) / 201),
        "p_P1_over_P2": float((1 + np.sum(np.array(nullR) >= o["P1_over_P2"])) / 201),
        "null_S_q50_95_max": np.quantile(nullS, [0.5, 0.95, 1]).tolist(),
        "null_P1_over_P2_q50_95_max": np.quantile(nullR, [0.5, 0.95, 1]).tolist()}
    return out


def synth(kind, n=405):
    i = np.arange(n)
    s = np.abs(i[:, None] - i[None, :]).astype(float)
    if kind == "ring":
        s = np.minimum(s, n - s)
    C = np.zeros((n, n))
    off = s > 0
    C[off] = s[off] ** -1.0
    return C


def diag_hic():
    out = {}
    for kind in ("ring", "chain"):
        C = synth(kind)
        D = T2.to_dist(C)
        X3, w = T2.classical_mds(D, 3)
        a = alpha_h1_summary(X3)
        r = rips_h1_summary(distance_matrix=D); r.pop("all_h1_persistence_sorted")
        out[f"synthetic_{kind}_noiseless"] = {"alpha_mds3": slim(a), "rips": slim(r, ("P1", "P2", "P1_over_P2", "S_rips", "top_h1_bars_birth_death")),
                                             "mds_frac_negative_eig_mass": float(np.sum(np.abs(w[w < 0])) / np.sum(np.abs(w)))}
    Craw, _ = T2.load_caulo()
    C, _ = T2.clean(Craw, masked_offset1=True)
    D = T2.to_dist(C)
    X3, w = T2.classical_mds(D, 3)
    out["caulo_mds_frac_negative_eig_mass"] = float(np.sum(np.abs(w[w < 0])) / np.sum(np.abs(w)))
    out["caulo_mds_top6_eig"] = w[:6].tolist()
    # genomic location of the top Rips H1 bars (birth/death edges)
    rc = gudhi.RipsComplex(distance_matrix=D)
    st = rc.create_simplex_tree(max_dimension=2)
    st.compute_persistence(homology_coeff_field=2)
    gens = st.flag_persistence_generators()
    h1 = gens[1][0] if len(gens[1]) else np.empty((0, 4), int)
    rows = []
    for (a1, b1, a2, b2) in h1:
        birth, death = D[a1, b1], D[a2, b2]
        rows.append((death - birth, birth, death, int(a1), int(b1), int(a2), int(b2)))
    rows.sort(reverse=True)
    out["caulo_rips_top5_h1_generators_bins"] = [
        {"pers": r[0], "birth": r[1], "death": r[2], "birth_edge_bins": [r[3], r[4]], "death_edge_bins": [r[5], r[6]],
         "birth_edge_circular_offset": int(min(abs(r[3] - r[4]), 405 - abs(r[3] - r[4])))} for r in rows[:5]]
    # inter-arm vs along-arm contacts: mean contact at pairs (i, N-i) (arm-mirror) vs pairs at offset 20
    n = C.shape[0]
    mirror = np.mean([C[i, n - 1 - i] for i in range(20, 180)])
    off20 = float(np.mean(np.diagonal(C, 20)))
    far = float(np.mean(np.diagonal(C, 150)))
    out["caulo_contact_mirror_pairs_i_Nminus1minusi_20_180"] = float(mirror)
    out["caulo_contact_mean_offset20"] = off20
    out["caulo_contact_mean_offset150"] = far
    return out


if __name__ == "__main__":
    res = {"NOTE": "POST HOC, not pre-stated in expectations.json; does not change any verdict.",
           "hic": diag_hic()}
    print(res["hic"])
    for ds in ("u2os", "mesc"):
        res[ds] = diag_cells(ds)
        print(ds, {k: v for k, v in res[ds].items() if k != "subsample_250_seeds0_9"})
    dump(res, os.path.join(HERE, "posthoc_diagnostics.json"))
