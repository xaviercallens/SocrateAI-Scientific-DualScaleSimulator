"""Hi-C: circular bacterial chromosomes vs a linear human chromosome arm.

Known-answer axis: a circular genome should give H1 = 1 (one dominant loop); a
linear chromosome arm should give H1 = 0.

BINDING METHOD RULE (lesson (a) of the earlier genetics validation): the
contact-derived distance d = C^(-1/3) is NOT Euclidean -- the earlier run measured
41% negative eigenvalue mass in its classical MDS of the Caulobacter matrix -- so
it is fed to Rips DIRECTLY.  No MDS, no alpha.  bio_common has no MDS function.

Cases:
  caulo_full / caulo_cut_ter / caulo_cut_ori   Caulobacter crescentus GSM1120445
  ecoli_full / ecoli_cut_ter                   E. coli GSM2870407
  bsub_full  / bsub_cut_ter                    B. subtilis PY79 GSM1671399  (new here)
  gm12878_chr1q                                human GM12878 chr1q, LINEAR control

Null: a linear null built from the dataset's own distance-decay curve E(s) with
linear (non-wrapping) offsets and the observed log-residuals permuted across pairs;
a circular genome's wrap contacts are thereby destroyed while the decay is kept.

Usage (one case per process, foreground, then assemble):
  python run_hic.py --cases caulo_full --n-null 200
  python run_hic.py --assemble
Seeds: null seeds 0..n_null-1.
"""
from __future__ import annotations

import argparse
import glob
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, "/mnt/disks/disk-socrateai-local-1/wt-topo-bio/topodb_runs/biology")
from bio_common import (BIO_DATA, GEN_DATA, RESULTS, Block, S_stat, betti_from_bars,  # noqa: E402
                        dominance, finite, max_pers, rank_p, rips_collapsed_from_distance,
                        sha256, top_bars)

SEED = 20260920
ALPHA_EXP = 1.0 / 3.0
SCRIPT = "topodb_runs/biology/run_hic.py"
CACHE = RESULTS / "hic_cache"

ALL_CASES = ["caulo_full", "caulo_cut_ter", "caulo_cut_ori", "ecoli_full", "ecoli_cut_ter",
             "bsub_full", "bsub_cut_ter", "gm12878_chr1q"]


# ------------------------------------------------------------------- loading
def load_caulo():
    p = glob.glob(str(GEN_DATA / "caulobacter" / "GSM1120445_*after_normalization.txt.gz"))[0]
    return np.loadtxt(p), p, "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1120445"


def load_bsub():
    p = str(BIO_DATA / "hic" / "GSM1671399_bsubtilis_PY79_HindIII_HiC.matrix.txt.gz")
    return np.loadtxt(p), p, "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1671nnn/GSM1671399/suppl/GSM1671399_01_Rudnerlab_HindIII_HiC_PY79.matrix.txt.gz"


def load_gm():
    p = str(GEN_DATA / "gm12878" / "gm12878_insitu_combined_chr1_145000000_249250621_250kb_KR_observed.npy")
    return np.nan_to_num(np.load(p), nan=0.0), p, "Rao et al. 2014 GSE63525 GM12878 in-situ combined, chr1:145-249 Mb, 250 kb, KR"


def load_ecoli():
    import pandas as pd
    p = glob.glob(str(GEN_DATA / "ecoli" / "GSM2870407_*.txt.gz"))[0]
    df = pd.read_csv(p, sep="\t", index_col=0)
    df = df.loc[:, [c for c in df.columns if not str(c).startswith("Unnamed")]]
    M = df.values.astype(float)
    nb = (M.shape[0] + 1) // 2                     # 5 kb -> 10 kb
    B = np.zeros((nb, nb))
    idx = np.arange(M.shape[0]) // 2
    np.add.at(B, (idx[:, None], idx[None, :]), M)
    np.fill_diagonal(B, 0.0)
    rs = B.sum(1)
    med = np.median(rs)
    mad = 1.4826 * np.median(np.abs(rs - med))
    keep = rs >= med - 3 * mad
    B = B[np.ix_(keep, keep)]
    for _ in range(200):                            # ICE
        r = B.sum(1)
        r = r / r.mean()
        if float(np.max(np.abs(r - 1))) < 1e-6:
            break
        B = B / np.outer(r, r)
    return B, p, "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2870407"


def clean(C, masked_offset1):
    """Drop zero-sum bins; impute the processing-masked offset-1 diagonal with E(2);
    replace remaining off-diagonal zeros with half the smallest positive contact."""
    info = {"max_abs_asym": float(np.abs(C - C.T).max()), "n_bins_in": int(C.shape[0])}
    keep = C.sum(1) > 0
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
    fill = 0.5 * float(C[off & (C > 0)].min())
    nz = off & (C <= 0)
    info["n_offdiag_zero_replaced"] = int(nz.sum() // 2)
    C[nz] = fill
    np.fill_diagonal(C, 0.0)
    info["n_bins_out"] = int(n)
    return C, info


def coarsen(C, target=210):
    """Sum adjacent bins until at most `target` bins remain.

    Declared compute choice: Rips on the full 404-464 bin matrices costs ~4 s per
    diagram, so a 200-null test would need hours per case; at ~200 bins it costs
    ~0.1 s.  The SAME coarsening is applied to the observed matrix and to every null
    matrix, so the comparison is like for like, and the observed statistics are also
    reported at full resolution as a sensitivity control (stat
    `h1_dominance_full_resolution`).  200 bins around a chromosome is ample to
    resolve a single ring.
    """
    n = C.shape[0]
    k = int(np.ceil(n / target))
    if k <= 1:
        return C, 1
    m = (n // k) * k
    C = C[:m, :m]
    idx = np.arange(m) // k
    B = np.zeros((m // k, m // k))
    np.add.at(B, (idx[:, None], idx[None, :]), C)
    np.fill_diagonal(B, 0.0)
    return B, k


def to_dist(C, a=ALPHA_EXP):
    D = np.zeros_like(C)
    off = ~np.eye(C.shape[0], dtype=bool)
    D[off] = C[off] ** (-a)
    return D


def linear_null_factory(C, circular):
    n = C.shape[0]
    I, J = np.triu_indices(n, 2)
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
    return make, {"powerlaw_slope_fit_N8_to_N2": float(slope),
                  "E_offsets_2_5_10_half": [float(E[2]), float(E[5]), float(E[10]), float(E[half])]}


# ------------------------------------------------------------------- cases
def build_case(name):
    """Returns (contact matrix, circular?, provenance dict)."""
    if name.startswith("caulo"):
        raw, p, url = load_caulo()
        C, info = clean(raw, masked_offset1=True)
        meta = {"organism": "Caulobacter crescentus NA1000 swarmer", "bin_kb": 10,
                "citation": "Le, Imakaev, Mirny, Laub 2013 Science (GEO GSM1120445)"}
        if name == "caulo_cut_ter":
            order = np.r_[np.arange(222, C.shape[0]), np.arange(0, 182)]
            C = C[np.ix_(order, order)]
            info["cut"] = "ter cut: bins 182-221 removed, chain 222..end,0..181"
            return C, False, dict(meta, path=p, url=url, clean=info)
        if name == "caulo_cut_ori":
            order = np.arange(20, 385)
            C = C[np.ix_(order, order)]
            info["cut"] = "ori cut: bins 0-19 and 385-end removed"
            return C, False, dict(meta, path=p, url=url, clean=info)
        return C, True, dict(meta, path=p, url=url, clean=info)
    if name.startswith("bsub"):
        raw, p, url = load_bsub()
        C, info = clean(raw, masked_offset1=True)
        meta = {"organism": "Bacillus subtilis PY79 wild type", "bin_kb": 10,
                "citation": "Wang, Le, Hughes, Kohler, Rudner 2015 Genes Dev (GEO GSM1671399)"}
        if name == "bsub_cut_ter":
            n = C.shape[0]
            lo, hi = n // 2 - 20, n // 2 + 20      # 40 bins around the antipode of ori
            order = np.r_[np.arange(hi, n), np.arange(0, lo)]
            C = C[np.ix_(order, order)]
            info["cut"] = f"ter-region cut: bins {lo}-{hi - 1} removed, chain {hi}..end,0..{lo - 1}"
            return C, False, dict(meta, path=p, url=url, clean=info)
        return C, True, dict(meta, path=p, url=url, clean=info)
    if name.startswith("ecoli"):
        raw, p, url = load_ecoli()
        C, info = clean(raw, masked_offset1=False)
        meta = {"organism": "Escherichia coli wild type, minimal medium 30C", "bin_kb": 10,
                "citation": "Lioy et al. 2018 Cell (GEO GSM2870407)"}
        if name == "ecoli_cut_ter":
            n = C.shape[0]
            order = np.r_[np.arange(182, n), np.arange(0, 136)]
            C = C[np.ix_(order, order)]
            info["cut"] = "ter cut: bins 136-181 removed"
            return C, False, dict(meta, path=p, url=url, clean=info)
        return C, True, dict(meta, path=p, url=url, clean=info)
    if name == "gm12878_chr1q":
        raw, p, url = load_gm()
        C, info = clean(raw, masked_offset1=False)
        return C, False, {"organism": "Homo sapiens GM12878, chr1q (LINEAR control)",
                          "bin_kb": 250, "citation": "Rao et al. 2014 Cell (GEO GSE63525)",
                          "path": p, "url": url, "clean": info}
    raise ValueError(name)


def analyse(C):
    D = to_dist(C)
    dg = rips_collapsed_from_distance(D, max_hom_dim=1)
    return dg, {"h1_dominance_P1_over_P2": dominance(dg[1]), "h1_S": S_stat(dg[1]),
                "h1_max_persistence": max_pers(dg[1]),
                "h1_finite_bar_count": float(len(finite(dg[1]))),
                "n_bins": float(C.shape[0])}


def run_case(name, n_null, budget_sec):
    C_full, circular, meta = build_case(name)
    t = time.time()
    _, obs_full = analyse(C_full)
    wall_full = round(time.time() - t, 2)
    C, kfac = coarsen(C_full)
    meta["coarsening"] = {"factor": kfac, "n_bins_full": int(C_full.shape[0]),
                          "n_bins_used": int(C.shape[0]),
                          "effective_bin_kb": meta["bin_kb"] * kfac}
    meta["full_resolution_observed"] = obs_full
    meta["wall_full_resolution_sec"] = wall_full
    t = time.time()
    dg, obs = analyse(C)
    wall_obs = round(time.time() - t, 2)
    make, nullinfo = linear_null_factory(C, circular)
    CACHE.mkdir(parents=True, exist_ok=True)
    # the cache key carries the bin count: a cache written at another resolution is
    # a different statistic and must never be reused (this bit us once)
    cpath = CACHE / f"{name}_nbins{C.shape[0]}_null_n{n_null}.json"
    done = json.loads(cpath.read_text()) if cpath.exists() else {}
    t0 = time.time()
    for sd in range(n_null):
        if str(sd) in done:
            continue
        if time.time() - t0 > budget_sec:
            break
        _, s = analyse(make(sd))
        done[str(sd)] = [s["h1_S"], s["h1_dominance_P1_over_P2"]]
        if sd % 25 == 0:
            cpath.write_text(json.dumps(done))
    cpath.write_text(json.dumps(done))
    n_have = len(done)
    print(f"  {name}: {n_have}/{n_null} nulls cached ({round(time.time() - t0, 1)}s)")
    if n_have < n_null:
        return None
    nv = np.array([done[str(i)] for i in range(n_null)])
    out = {"case": name, "circular": circular, "meta": meta, "observed": obs,
           "wall_obs_sec": wall_obs, "n_null": n_null,
           "p_S_vs_linear_null": float(rank_p(obs["h1_S"], nv[:, 0])),
           "p_dominance_vs_linear_null": float(rank_p(obs["h1_dominance_P1_over_P2"], nv[:, 1])),
           "null_S_q50_q95_q99_max": [float(np.quantile(nv[:, 0], q)) for q in (.5, .95, .99, 1)],
           "null_frac_dominance_ge_2": float(np.mean(nv[:, 1] >= 2)),
           "null_info": nullinfo,
           "top_h1_bars": [[round(b, 5), round(d, 5)] for b, d in top_bars(dg[1], 5)],
           "b0_infinite": betti_from_bars(dg[0]),
           "diagram_h1": [[float(b), float(d)] for b, d in dg[1]],
           "diagram_h0": [[float(b), (None if not np.isfinite(d) else float(d))] for b, d in dg[0]]}
    (RESULTS / f"hic_{name}.json").write_text(json.dumps(out, indent=1))
    print(f"  {name}: dom={obs['h1_dominance_P1_over_P2']:.3f} S={obs['h1_S']:.4f} "
          f"p_S={out['p_S_vs_linear_null']:.4f} p_dom={out['p_dominance_vs_linear_null']:.4f}")
    return out


def assemble():
    blk = Block("hic", SCRIPT, "see per-run command field")
    summary = {}
    for name in ALL_CASES:
        f = RESULTS / f"hic_{name}.json"
        if not f.exists():
            summary[name] = {"status": "ABSENT: not computed"}
            continue
        r = json.loads(f.read_text())
        meta, obs = r["meta"], r["observed"]
        ds_id = f"biology/hic_{name}"
        blk.dataset(id=ds_id, domain="biology",
                    title=f"Hi-C contact map: {meta['organism']}"
                          + (f" [{meta['clean'].get('cut')}]" if meta["clean"].get("cut") else ""),
                    source=meta["url"], provenance="experiment",
                    n_objects=int(obs["n_bins"]), ambient_dim=int(obs["n_bins"]),
                    units=f"{meta['bin_kb']} kb bins; distance = contact^(-1/3), dimensionless",
                    sha256=sha256(meta["path"]), local_path=meta["path"],
                    notes=f"{meta['citation']}. cleaning: {json.dumps(meta['clean'])}. "
                          f"circular chromosome: {r['circular']}. coarsening: "
                          f"{json.dumps(meta['coarsening'])}. The contact-derived metric is "
                          f"NOT Euclidean and is fed to Rips directly; no MDS is computed "
                          f"anywhere in this block (lesson (a)).")
        dom = obs["h1_dominance_P1_over_P2"]
        p_S = r["p_S_vs_linear_null"]
        expect_b1 = 1 if r["circular"] else 0
        got_b1 = 1 if (dom >= 2 and p_S <= 0.01) else 0
        blk.run(dataset_id=ds_id, method="rips", coeff_field=2, max_dim=1,
                params={"metric": "contact^(-1/3) fed directly to Rips (NO embedding)",
                        "exponent": ALPHA_EXP, "max_hom_dim": 1, "max_edge_length": None,
                        "edge_collapse": True, "expansion_dim": 2, "n_bins": obs["n_bins"],
                        "bin_kb": meta["bin_kb"], "circular_chromosome": r["circular"],
                        "coarsening": meta["coarsening"],
                        "criterion_b1_eq_1": "dominance >= 2 AND p_S <= 0.01 vs the linear null"},
                preprocessing=json.dumps(meta["clean"]), seed=f"null seeds 0..{r['n_null'] - 1}",
                tier="X", wall_sec=r["wall_obs_sec"],
                diagrams={0: [(b, d) for b, d in r["diagram_h0"]] if False else
                          [(b, float("inf") if d is None else d) for b, d in r["diagram_h0"]],
                          1: [(b, d) for b, d in r["diagram_h1"]]},
                betti={1: got_b1}, expected_betti={1: expect_b1},
                stats=[{"name": "h1_dominance_P1_over_P2", "value": dom,
                        "null_model": f"linear null: E(s) decay of the same matrix with "
                                      f"non-wrapping offsets, log-residuals permuted across pairs",
                        "n_null": r["n_null"], "p_value": r["p_dominance_vs_linear_null"],
                        "p_method": "rank"},
                       {"name": "h1_S_longest_over_total", "value": obs["h1_S"],
                        "null_model": "same linear null", "n_null": r["n_null"],
                        "p_value": p_S, "p_method": "rank"},
                       {"name": "h1_max_persistence", "value": obs["h1_max_persistence"]},
                       {"name": "h1_finite_bar_count", "value": obs["h1_finite_bar_count"]},
                       {"name": "null_frac_dominance_ge_2", "value": r["null_frac_dominance_ge_2"]},
                       {"name": "h1_dominance_full_resolution",
                        "value": meta["full_resolution_observed"]["h1_dominance_P1_over_P2"]},
                       {"name": "h1_S_full_resolution",
                        "value": meta["full_resolution_observed"]["h1_S"]}],
                controls=[{"kind": "null_calibration",
                           "description": "linear null destroys the circular wrap while keeping the "
                                          "contact-distance decay of the same dataset",
                           "passed": bool(r["null_frac_dominance_ge_2"] < 0.1),
                           "detail": f"{r['null_frac_dominance_ge_2']:.3f} of {r['n_null']} null "
                                     f"matrices reach dominance >= 2"},
                          {"kind": "injection",
                           "description": "resolution sensitivity: the same statistic recomputed on "
                                          "the FULL-resolution matrix, which has no null attached "
                                          "because it costs ~40x more per diagram",
                           "passed": bool((meta["full_resolution_observed"]["h1_dominance_P1_over_P2"] >= 2)
                                          == (dom >= 2)),
                           "detail": f"coarsened ({meta['coarsening']['n_bins_used']} bins, factor "
                                     f"{meta['coarsening']['factor']}) dominance {dom:.3f} vs full "
                                     f"({meta['coarsening']['n_bins_full']} bins) "
                                     f"{meta['full_resolution_observed']['h1_dominance_P1_over_P2']:.3f}"},
                          {"kind": "known_answer",
                           "description": f"{'circular' if r['circular'] else 'linear'} genome -> "
                                          f"b1 = {expect_b1}",
                           "passed": bool(got_b1 == expect_b1),
                           "detail": f"dominance {dom:.3f}, p_S {p_S:.4f} -> b1 = {got_b1}"}],
                findings=[{"claim": f"{name}: Rips on the contact metric gives H1 dominance "
                                    f"{dom:.2f} (p {r['p_dominance_vs_linear_null']:.4f}) and "
                                    f"S = {obs['h1_S']:.4f} (p {p_S:.4f}) against a linear null of "
                                    f"{r['n_null']}. Expected b1 = {expect_b1} "
                                    f"({'circular' if r['circular'] else 'linear/cut'}); observed "
                                    f"b1 = {got_b1}.",
                           "verdict": ("recovered" if got_b1 == expect_b1 == 1 else
                                       "null" if got_b1 == expect_b1 == 0 else "failed"),
                           "tier": "X",
                           "caveat": "the contact-to-distance exponent 1/3 is a modelling choice, "
                                     "not a measurement; no MDS is used, so this is not the false-loop "
                                     "route of the earlier validation",
                           "reference": "topodb_runs/biology/expectations.json"}])
        summary[name] = {"circular": r["circular"], "dominance": round(dom, 4),
                         "S": round(obs["h1_S"], 5), "p_S": round(p_S, 5),
                         "p_dom": round(r["p_dominance_vs_linear_null"], 5),
                         "expected_b1": expect_b1, "observed_b1": got_b1,
                         "n_bins": int(obs["n_bins"]), "n_null": r["n_null"],
                         "dominance_full_res": round(
                             meta["full_resolution_observed"]["h1_dominance_P1_over_P2"], 4),
                         "n_bins_full": meta["coarsening"]["n_bins_full"]}
    (RESULTS / "hic_summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))
    blk.write()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", default="")
    ap.add_argument("--n-null", type=int, default=200)
    ap.add_argument("--budget-sec", type=float, default=480)
    ap.add_argument("--assemble", action="store_true")
    a = ap.parse_args()
    if a.assemble:
        assemble()
        return 0
    for name in (a.cases.split(",") if a.cases else ALL_CASES):
        run_case(name.strip(), a.n_null, a.budget_sec)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
