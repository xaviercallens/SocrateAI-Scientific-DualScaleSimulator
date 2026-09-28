"""Single-cell RNA-seq: Rips on the cell-cell correlation metric.

Three datasets:
  u2os    GSE146773 U2OS FUCCI  (Mahdessian et al. 2021 Nature)   -- cell-cycle loop expected
  mesc    E-MTAB-2805 mESC      (Buettner et al. 2015 Nat Biotech) -- cell-cycle loop expected
  pbmc3k  10x Genomics PBMC 3k                                     -- BREADTH only, no loop expected

METHOD: d = sqrt(2 (1 - rho_spearman)) between cells over the selected genes, fed to
Rips DIRECTLY (no PCA, no MDS, no embedding).  This differs from the earlier
genetics validation, which ran alpha on a 3-D PCA and FAILED both cell-cycle
datasets (P1/P2 = 1.53 and 1.61; gene-permutation p = 0.37 and 0.08).  The present
setting was pre-stated in expectations.json before running, and whatever it gives
is reported -- a second failure is the finding, and nothing is tuned.

NULL: gene-wise permutation -- each gene's values are independently permuted across
cells, which destroys the cell-cell covariance structure while preserving every
gene's marginal distribution; the metric and the diagram are recomputed each time.
NEGATIVE CONTROL: the same number of randomly chosen expressed non-cell-cycle genes.

Pre-stated criteria (expectations.json): (i) dominance >= 2, (ii) p <= 0.01.

Usage:  python run_scrna.py --dataset u2os --n-null 200 [--budget-sec 480]
        python run_scrna.py --assemble
Seeds: null seeds 0..n-1; cell subsample seed 20260920.
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

sys.path.insert(0, "/mnt/disks/disk-socrateai-local-1/wt-topo-bio/topodb_runs/biology")
from bio_common import (BIO_DATA, GEN_DATA, RESULTS, Block, S_stat, dominance, finite,  # noqa: E402
                        max_pers, rank_p, rips_collapsed_from_distance, sha256, top_bars)

SEED = 20260920
N_CELLS_MAX = 400
SCRIPT = "topodb_runs/biology/run_scrna.py"
CACHE = RESULTS / "scrna_cache"
DATASETS = ["u2os", "mesc", "pbmc3k"]


# ---------------------------------------------------------------------- loading
def mad_low(v):
    m = np.median(v)
    s = 1.4826 * np.median(np.abs(v - m))
    return m - 3 * s, m + 3 * s


def load_u2os():
    d = GEN_DATA / "u2os"
    counts = pd.read_csv(d / "GSE146773_Counts.csv.gz", index_col=0)
    fucci = pd.read_csv(d / "GSE146773_fucci_coords.csv.gz").set_index("cell")
    fucci = fucci.loc[counts.index.intersection(fucci.index)]
    fucci = fucci[np.isfinite(fucci["polar_coord_phi"].astype(float))]
    counts = counts.loc[fucci.index]
    genes = [g for g in counts.columns if not g.startswith("ERCC")]
    ercc = [g for g in counts.columns if g.startswith("ERCC")]
    cc = pd.read_csv(GEN_DATA / "genelists" / "Homo_sapiens.csv")
    return (counts[genes].values.astype(float), np.array(genes),
            counts[ercc].values.astype(float) if ercc else None,
            {"phi": fucci["polar_coord_phi"].values.astype(float)}, cc,
            str(d / "GSE146773_Counts.csv.gz"),
            "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146773",
            "Mahdessian et al. 2021 Nature 590:649 (U2OS FUCCI scRNA-seq)")


def load_mesc():
    d = GEN_DATA / "buettner"
    mats, labs, ids = [], [], None
    for ph in ("G1", "S", "G2M"):
        t = pd.read_csv(d / f"{ph}_singlecells_counts.txt", sep="\t")
        ids = t["EnsemblGeneID"].values
        m = t.iloc[:, 4:].values.astype(float).T
        mats.append(m)
        labs += [ph] * m.shape[0]
    X_all = np.vstack(mats)
    is_e = np.array([str(g).startswith("ERCC") for g in ids])
    cc = pd.read_csv(GEN_DATA / "genelists" / "Mus_musculus.csv")
    return (X_all[:, ~is_e], ids[~is_e], X_all[:, is_e], {"phase": np.array(labs)}, cc,
            str(d / "G1_singlecells_counts.txt"),
            "https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-2805/",
            "Buettner et al. 2015 Nat Biotechnol 33:155 (FACS-sorted mESC)")


def load_pbmc3k():
    from scipy.io import mmread
    d = BIO_DATA / "scrna" / "filtered_gene_bc_matrices" / "hg19"
    # kept sparse and gene-filtered BEFORE densifying: the full 32738 x 2700 float64
    # array is ~707 MB, which is wasteful under --as=8G beside other running jobs
    S = mmread(str(d / "matrix.mtx")).tocsr()                 # genes x cells, sparse
    keep = np.asarray((S > 0).sum(axis=1)).ravel() >= 0.05 * S.shape[1]
    g = pd.read_csv(d / "genes.tsv", sep="\t", header=None)
    M = np.asarray(S[keep].todense(), dtype=float).T          # cells x kept genes
    return (M, g[0].values[keep], None, {}, None,
            str(BIO_DATA / "scrna" / "pbmc3k_filtered_gene_bc_matrices.tar.gz"),
            "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/"
            "pbmc3k_filtered_gene_bc_matrices.tar.gz",
            "10x Genomics PBMC 3k (the dataset scanpy.datasets.pbmc3k downloads)")


LOADERS = {"u2os": load_u2os, "mesc": load_mesc, "pbmc3k": load_pbmc3k}


# --------------------------------------------------------------------- pipeline
def qc(X, E):
    lib = np.log10(X.sum(1) + 1)
    det = (X > 0).sum(1).astype(float)
    keep = (lib >= mad_low(lib)[0]) & (det >= mad_low(det)[0])
    thr = {"log10_lib_min": float(mad_low(lib)[0]), "detected_min": float(mad_low(det)[0])}
    if E is not None and E.size:
        frac = E.sum(1) / (E.sum(1) + X.sum(1))
        keep &= frac <= mad_low(frac)[1]
        thr["ercc_frac_max"] = float(mad_low(frac)[1])
    return keep, thr


def zscore(M):
    sd = M.std(0)
    sd[sd == 0] = 1.0
    return (M - M.mean(0)) / sd


def spearman_dist(Z):
    """d = sqrt(2 (1 - rho)) between cells: a proper metric on the correlation sphere."""
    R = np.apply_along_axis(rankdata, 1, Z)
    R = R - R.mean(1, keepdims=True)
    n = np.linalg.norm(R, axis=1, keepdims=True)
    n[n == 0] = 1.0
    rho = (R / n) @ (R / n).T
    np.clip(rho, -1.0, 1.0, out=rho)
    D = np.sqrt(np.maximum(2.0 * (1.0 - rho), 0.0))
    np.fill_diagonal(D, 0.0)
    return D


def summarise(Z):
    dg = rips_collapsed_from_distance(spearman_dist(Z), max_hom_dim=1)
    return dg, {"h1_dominance_P1_over_P2": dominance(dg[1]), "h1_S": S_stat(dg[1]),
                "h1_max_persistence": max_pers(dg[1]),
                "h1_finite_bar_count": float(len(finite(dg[1])))}


def gene_perm(Z, rng):
    Zp = Z.copy()
    for j in range(Z.shape[1]):
        Zp[:, j] = Z[rng.permutation(Z.shape[0]), j]
    return Zp


def prepare(name):
    X, genes, E, lab, cc, path, url, citation = LOADERS[name]()
    keep, thr = qc(X, E)
    X = X[keep]
    lab = {k: v[keep] for k, v in lab.items()}
    L = np.log2(X / np.maximum(X.sum(1, keepdims=True), 1) * 1e6 + 1)
    det = (X > 0).mean(0)
    if cc is not None:
        cc_ids = set(cc["geneID"].astype(str))
        is_cc = np.array([str(g) in cc_ids for g in genes])
        sel = is_cc & (det >= 0.05)
        gene_set = "cell-cycle genes (tinyatlas list) detected in >= 5% of cells"
    else:
        expressed = det >= 0.05
        v = L.var(0)
        v[~expressed] = -np.inf
        sel = np.zeros(len(genes), dtype=bool)
        sel[np.argsort(v)[-500:]] = True
        gene_set = "500 most variable genes detected in >= 5% of cells (no cell-cycle list used)"
    neg_pool = np.where((det >= 0.05) & ~sel)[0]
    rng = np.random.default_rng(SEED)
    idx = np.arange(L.shape[0])
    if L.shape[0] > N_CELLS_MAX:
        idx = np.sort(rng.choice(L.shape[0], N_CELLS_MAX, replace=False))
    info = {"n_cells_loaded": int(len(keep)), "n_cells_after_qc": int(keep.sum()),
            "n_cells_used": int(len(idx)), "qc_thresholds": thr, "n_genes_selected": int(sel.sum()),
            "gene_set": gene_set, "subsample_seed": SEED if L.shape[0] > N_CELLS_MAX else None,
            "path": path, "url": url, "citation": citation}
    return zscore(L[np.ix_(idx, np.where(sel)[0])]), L[idx], neg_pool, int(sel.sum()), lab, idx, info


def run_dataset(name, n_null, n_neg, budget):
    Z, Lsub, neg_pool, n_sel, lab, idx, info = prepare(name)
    t = time.time()
    dg, obs = summarise(Z)
    wall = round(time.time() - t, 2)
    print(f"  {name}: n_cells={Z.shape[0]} n_genes={Z.shape[1]} dom={obs['h1_dominance_P1_over_P2']:.3f} "
          f"S={obs['h1_S']:.4f} ({wall}s)")

    CACHE.mkdir(parents=True, exist_ok=True)
    cp = CACHE / f"{name}_genenull_n{n_null}.json"
    done = json.loads(cp.read_text()) if cp.exists() else {}
    t0 = time.time()
    for sd in range(n_null):
        if str(sd) in done:
            continue
        if time.time() - t0 > budget:
            break
        _, s = summarise(gene_perm(Z, np.random.default_rng(sd)))
        done[str(sd)] = [s["h1_S"], s["h1_dominance_P1_over_P2"]]
        if sd % 10 == 0:
            cp.write_text(json.dumps(done))
    cp.write_text(json.dumps(done))
    print(f"    gene-permutation null {len(done)}/{n_null} ({round(time.time() - t0, 1)}s)")
    if len(done) < n_null:
        return None

    # negative control: random non-selected expressed genes, same count
    cn = CACHE / f"{name}_neg_n{n_neg}.json"
    negd = json.loads(cn.read_text()) if cn.exists() else {}
    t1 = time.time()
    for sd in range(n_neg):
        if str(sd) in negd:
            continue
        if time.time() - t1 > budget / 3:
            break
        rng = np.random.default_rng(1000 + sd)
        pick = rng.choice(neg_pool, min(n_sel, len(neg_pool)), replace=False)
        _, s = summarise(zscore(Lsub[:, pick]))
        negd[str(sd)] = [s["h1_S"], s["h1_dominance_P1_over_P2"]]
        cn.write_text(json.dumps(negd))
    cn.write_text(json.dumps(negd))

    nv = np.array([done[str(i)] for i in range(n_null)])
    ng = np.array([negd[str(i)] for i in sorted(map(int, negd))]) if negd else np.empty((0, 2))
    out = {"dataset": name, "info": info, "observed": obs, "wall_obs_sec": wall,
           "n_null": n_null, "p_S": float(rank_p(obs["h1_S"], nv[:, 0])),
           "p_dominance": float(rank_p(obs["h1_dominance_P1_over_P2"], nv[:, 1])),
           "null_S_q50_q95_q99_max": [float(np.quantile(nv[:, 0], q)) for q in (.5, .95, .99, 1)],
           "null_frac_dominance_ge_2": float(np.mean(nv[:, 1] >= 2)),
           "negative_control": {"n_draws": int(len(ng)),
                                "median_dominance": float(np.median(ng[:, 1])) if len(ng) else None,
                                "frac_dominance_ge_2": float(np.mean(ng[:, 1] >= 2)) if len(ng) else None,
                                "median_S": float(np.median(ng[:, 0])) if len(ng) else None},
           "top_h1_bars": [[round(b, 5), round(d, 5)] for b, d in top_bars(dg[1], 5)],
           "diagram_h1": [[float(b), float(d)] for b, d in dg[1]]}
    (RESULTS / f"scrna_{name}.json").write_text(json.dumps(out, indent=1))
    print(f"    p_S={out['p_S']:.4f} p_dom={out['p_dominance']:.4f} "
          f"neg_median_dom={out['negative_control']['median_dominance']}")
    return out


def observed_only(name):
    """Record the observed diagram with NO p-value, when the null has not finished.

    The database refuses a p_value without its null model, which is the right rule;
    this writes what was actually measured and says plainly that no null completed,
    rather than leaving a downloaded, hashed dataset out of the record entirely.
    """
    Z, Lsub, neg_pool, n_sel, lab, idx, info = prepare(name)
    t = time.time()
    dg, obs = summarise(Z)
    wall = round(time.time() - t, 2)
    cp = CACHE / f"{name}_genenull_n200.json"
    n_done = len(json.loads(cp.read_text())) if cp.exists() else 0
    out = {"dataset": name, "info": info, "observed": obs, "wall_obs_sec": wall,
           "n_null": 0, "null_incomplete_draws_reached": n_done,
           "no_pvalue_reason": f"the gene-permutation null did not complete: {n_done} of 200 "
                               f"draws at roughly 13 s each on a {Z.shape[0]}-cell correlation "
                               f"metric. No p-value is stored, because a p-value without its "
                               f"null is not storable.",
           "negative_control": {"n_draws": 0, "median_dominance": None,
                                "frac_dominance_ge_2": None, "median_S": None},
           "top_h1_bars": [[round(b, 5), round(d, 5)] for b, d in top_bars(dg[1], 5)],
           "diagram_h1": [[float(b), float(d)] for b, d in dg[1]]}
    (RESULTS / f"scrna_{name}_observed_only.json").write_text(json.dumps(out, indent=1))
    print(f"  {name} OBSERVED ONLY: n_cells={Z.shape[0]} n_genes={Z.shape[1]} "
          f"dom={obs['h1_dominance_P1_over_P2']:.3f} S={obs['h1_S']:.4f} "
          f"(null {n_done}/200, no p-value stored)")
    return out


def assemble():
    blk = Block("scrna", SCRIPT, "see per-run command field")
    summary = {}
    for name in DATASETS:
        f = RESULTS / f"scrna_{name}.json"
        if not f.exists():
            f = RESULTS / f"scrna_{name}_observed_only.json"
        if not f.exists():
            summary[name] = {"status": "ABSENT: not computed"}
            continue
        r = json.loads(f.read_text())
        has_null = r["n_null"] > 0
        info, obs = r["info"], r["observed"]
        expects_loop = name != "pbmc3k"
        ds_id = f"biology/scrna_{name}"
        blk.dataset(id=ds_id, domain="biology",
                    title=f"scRNA-seq {name}: {info['n_cells_used']} cells x "
                          f"{info['n_genes_selected']} genes",
                    source=info["url"], provenance="experiment",
                    n_objects=info["n_cells_used"], ambient_dim=info["n_genes_selected"],
                    units="d = sqrt(2(1 - Spearman rho)) between cells, dimensionless",
                    sha256=sha256(info["path"]), local_path=info["path"],
                    notes=f"{info['citation']}. genes: {info['gene_set']}. QC "
                          f"{json.dumps(info['qc_thresholds'])}; {info['n_cells_after_qc']} cells "
                          f"pass QC, {info['n_cells_used']} used"
                          + (f" (subsample seed {info['subsample_seed']})"
                             if info["subsample_seed"] else "")
                          + (". BREADTH dataset: discrete immune cell types, no cell-cycle loop "
                             "is expected and none is claimed." if not expects_loop else ""))
        dom = obs["h1_dominance_P1_over_P2"]
        ok = bool(has_null and dom >= 2 and r["p_S"] <= 0.01)
        blk.run(dataset_id=ds_id, method="rips", coeff_field=2, max_dim=1,
                params={"metric": "sqrt(2(1-spearman rho)) between cells, fed directly to Rips "
                                  "(NO PCA, NO MDS)", "max_hom_dim": 1, "max_edge_length": None,
                        "edge_collapse": True, "expansion_dim": 2,
                        "n_cells": info["n_cells_used"], "n_genes": info["n_genes_selected"],
                        "cell_subsample_seed": info["subsample_seed"],
                        "criteria": "dominance >= 2 AND p_S <= 0.01"},
                preprocessing=f"CPM log2(x+1); {info['gene_set']}; z-scored per gene",
                seed=f"cells {SEED}; null seeds 0..{r['n_null'] - 1}", tier="X",
                wall_sec=r["wall_obs_sec"],
                diagrams={1: [(b, d) for b, d in r["diagram_h1"]]},
                betti={1: 1 if ok else 0},
                expected_betti=({1: 1} if expects_loop else None),
                stats=([{"name": "h1_dominance_P1_over_P2", "value": dom,
                         "null_model": "gene-wise permutation across cells", "n_null": r["n_null"],
                         "p_value": r["p_dominance"], "p_method": "rank"},
                        {"name": "h1_S_longest_over_total", "value": obs["h1_S"],
                         "null_model": "gene-wise permutation across cells", "n_null": r["n_null"],
                         "p_value": r["p_S"], "p_method": "rank"},
                        {"name": "null_frac_dominance_ge_2",
                         "value": r["null_frac_dominance_ge_2"]}] if has_null else
                       [{"name": "h1_dominance_P1_over_P2", "value": dom},
                        {"name": "h1_S_longest_over_total", "value": obs["h1_S"]},
                        {"name": "null_draws_reached_before_stopping",
                         "value": float(r.get("null_incomplete_draws_reached", 0))}])
                      + [{"name": "h1_max_persistence", "value": obs["h1_max_persistence"]},
                         {"name": "h1_finite_bar_count", "value": obs["h1_finite_bar_count"]}]
                      + ([{"name": "negative_control_median_dominance",
                           "value": r["negative_control"]["median_dominance"]}]
                         if r["negative_control"]["median_dominance"] is not None else []),
                controls=[{"kind": "shuffle",
                           "description": "gene-wise permutation destroys cell-cell covariance and "
                                          "keeps every gene's marginal",
                           "passed": (bool(r["null_frac_dominance_ge_2"] < 0.1) if has_null
                                      else None),
                           "detail": (f"{r['null_frac_dominance_ge_2']:.3f} of {r['n_null']} "
                                      f"permutations reach dominance >= 2; null S q95 = "
                                      f"{r['null_S_q50_q95_q99_max'][1]:.4f}") if has_null
                                     else f"NOT RUN TO COMPLETION: {r.get('no_pvalue_reason')}"},
                          {"kind": "negative",
                           "description": "the same number of random expressed genes outside the "
                                          "selected set",
                           "passed": (None if r["negative_control"]["median_dominance"] is None
                                      else bool(r["negative_control"]["frac_dominance_ge_2"] < 0.5)),
                           "detail": f"{r['negative_control']['n_draws']} draws, median dominance "
                                     f"{r['negative_control']['median_dominance']}"}]
                         + ([{"kind": "known_answer",
                              "description": "published cell-cycle loop should appear as a dominant H1",
                              "passed": (bool(ok) if has_null else None),
                              "detail": (f"dominance {dom:.3f}, p_S {r['p_S']:.4f}" if has_null
                                         else f"dominance {dom:.3f}; NOT RUN TO COMPLETION: "
                                              f"{r.get('no_pvalue_reason')}")}]
                            if expects_loop else []),
                findings=[{"claim": f"{name}: Rips on the cell-cell Spearman metric over "
                                    f"{info['n_genes_selected']} genes gives H1 dominance "
                                    f"{dom:.2f} and S {obs['h1_S']:.4f}"
                                    + (f" (p {r['p_dominance']:.4f} and {r['p_S']:.4f} against "
                                       f"{r['n_null']} gene-wise permutations)." if has_null
                                       else f". NO p-value: {r.get('no_pvalue_reason')}")
                                    + (f" Pre-stated criteria (dominance >= 2 and p <= 0.01): "
                                       f"{'MET' if ok else 'NOT MET'}." if expects_loop else
                                       " No loop was expected and none is claimed."),
                           "verdict": ("recovered" if (expects_loop and ok) else
                                       "failed" if expects_loop else "null"),
                           "tier": "X",
                           "caveat": ("This is a different pipeline from the earlier genetics "
                                      "validation (Rips on the correlation metric in full dimension, "
                                      "not alpha on a 3-D PCA), pre-stated before running. Nothing "
                                      "was tuned after seeing the result."
                                      if expects_loop else
                                      "breadth dataset; the gene set is the 500 most variable genes, "
                                      "not a cell-cycle list"),
                           "reference": "topodb_runs/biology/expectations.json"}])
        summary[name] = {"dominance": round(dom, 4), "S": round(obs["h1_S"], 5),
                         "p_S": (round(r["p_S"], 5) if has_null else None),
                         "p_dom": (round(r["p_dominance"], 5) if has_null else None),
                         "null_complete": has_null,
                         "null_draws_reached": r.get("null_incomplete_draws_reached", r["n_null"]),
                         "n_cells": info["n_cells_used"], "n_genes": info["n_genes_selected"],
                         "expects_loop": expects_loop, "criteria_met": bool(ok),
                         "neg_median_dom": r["negative_control"]["median_dominance"]}
    (RESULTS / "scrna_summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))
    blk.write()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default="")
    ap.add_argument("--n-null", type=int, default=200)
    ap.add_argument("--n-neg", type=int, default=10)
    ap.add_argument("--budget-sec", type=float, default=440)
    ap.add_argument("--assemble", action="store_true")
    ap.add_argument("--observed-only", action="store_true")
    a = ap.parse_args()
    if a.assemble:
        assemble()
        return 0
    if a.observed_only:
        for name in (a.dataset.split(",") if a.dataset else DATASETS):
            observed_only(name.strip())
        return 0
    for name in (a.dataset.split(",") if a.dataset else DATASETS):
        run_dataset(name.strip(), a.n_null, a.n_neg, a.budget_sec)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
