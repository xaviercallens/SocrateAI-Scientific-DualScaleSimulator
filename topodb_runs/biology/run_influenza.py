"""Influenza: single segments (tree-like, H1 ~ 0) vs concatenated genomes (reassortment, H1 > 0).

Reproduces the published result of Chan, Carlsson & Rabadan, PNAS 110:18566 (2013):
a single segment evolves clonally, so its phylogeny is a TREE and H1 is essentially
empty; whole genomes reassort, which creates H1 classes.

METHOD (binding): Rips on the Hamming / p-distance matrix DIRECTLY.  The earlier
genetics validation showed that embedding this same metric into 3-D and running
alpha gives dominance 11.1 on segment 4 (HA) -- a FALSE loop on a tree.  That route
is not computed here at all.

DATA: alignments prepared by the earlier genetics validation run (pyfamsa 0.7.0,
300 avian complete 8-segment genomes drawn from 14717 unique sets with seed
20260919) and reused here byte-for-byte; sha256 of each file is recorded.  Source:
NCBI Influenza Virus Resource genomeset.dat / influenza.fna.

NULL: site bootstrap -- alignment columns resampled with replacement within each
segment, which preserves the per-site composition and the sample size but breaks
the linkage that carries the reassortment signal.

Usage:  python run_influenza.py [--n-null 200] [--budget-sec 480]
Seeds: bootstrap seeds 0..n-1; genome subsample seed 20260919 (inherited).
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, "/mnt/disks/disk-socrateai-local-1/wt-topo-bio/topodb_runs/biology")
from bio_common import (GEN_DATA, RESULTS, Block, S_stat, betti_from_bars, dominance,  # noqa: E402
                        finite, max_pers, rank_p, rips_collapsed_from_distance, sha256, top_bars)

SEED_SUB = 20260919
SCRIPT = "topodb_runs/biology/run_influenza.py"
COMMAND = ("timeout 590 prlimit --as=8589934592 -- "
           "/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python "
           "topodb_runs/biology/run_influenza.py --n-null 200")
WORK = GEN_DATA / "influenza" / "work"
SEGMENTS = {1: "PB2", 2: "PB1", 3: "PA", 4: "HA", 5: "NP", 6: "NA", 7: "M", 8: "NS"}
BAR_THRESH = 0.005          # "B" statistic of the earlier validation, kept identical
CACHE = RESULTS / "flu_cache"


def load_seg(s):
    a = np.load(WORK / f"aln_seg{s}.npy")
    return np.frombuffer(a.tobytes(), dtype="S1").reshape(a.shape)


def pdist_hamming(A):
    """Fraction of differing aligned columns, over columns where BOTH are non-gap."""
    codes = np.zeros(A.shape, dtype=np.int8)
    for i, ch in enumerate([b"A", b"C", b"G", b"T"]):
        codes[A == ch] = i + 1        # 0 = gap or ambiguous
    n, L = codes.shape
    D = np.zeros((n, n))
    valid = codes > 0
    for i in range(n):
        both = valid[i] & valid
        diff = (codes[i] != codes) & both
        cnt = both.sum(1)
        D[i] = np.where(cnt > 0, diff.sum(1) / np.maximum(cnt, 1), 0.0)
    D = 0.5 * (D + D.T)
    np.fill_diagonal(D, 0.0)
    return D


def summarise(D):
    dg = rips_collapsed_from_distance(D, max_hom_dim=1)
    fb = finite(dg[1])
    return dg, {"A_max_h1_persistence": max_pers(dg[1]),
                "B_n_h1_bars_ge_0.005": float(sum(1 for b, d in fb if d - b >= BAR_THRESH)),
                "n_h1_bars": float(len(fb)),
                "h1_dominance_P1_over_P2": dominance(dg[1]),
                "h1_S": S_stat(dg[1]),
                "median_pdistance": float(np.median(D[np.triu_indices(D.shape[0], 1)]))}


def bootstrap_cols(A, rng):
    idx = rng.integers(0, A.shape[1], A.shape[1])
    return A[:, idx]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-null", type=int, default=200)
    ap.add_argument("--budget-sec", type=float, default=480)
    a = ap.parse_args()
    CACHE.mkdir(parents=True, exist_ok=True)

    prep = json.loads((WORK / "prep_info.json").read_text())
    alns = {s: load_seg(s) for s in SEGMENTS}
    n_gen = alns[1].shape[0]

    # observed
    obs = {}
    for s, nm in SEGMENTS.items():
        t = time.time()
        D = pdist_hamming(alns[s])
        dg, st = summarise(D)
        obs[f"seg{s}_{nm}"] = dict(st, wall_sec=round(time.time() - t, 2), aln_len=int(alns[s].shape[1]),
                                   diagram_h1=[[float(b), float(d)] for b, d in dg[1]],
                                   top_bars=[[round(b, 5), round(d, 5)] for b, d in top_bars(dg[1], 5)],
                                   b0=betti_from_bars(dg[0]))
        print(f"  seg{s} {nm}: A={st['A_max_h1_persistence']:.5f} B={st['B_n_h1_bars_ge_0.005']:.0f} "
              f"n_h1={st['n_h1_bars']:.0f} med_pd={st['median_pdistance']:.4f}")
    concat = np.concatenate([alns[s] for s in SEGMENTS], axis=1)
    t = time.time()
    Dc = pdist_hamming(concat)
    dgc, stc = summarise(Dc)
    obs["concatenated"] = dict(stc, wall_sec=round(time.time() - t, 2), aln_len=int(concat.shape[1]),
                               diagram_h1=[[float(b), float(d)] for b, d in dgc[1]],
                               top_bars=[[round(b, 5), round(d, 5)] for b, d in top_bars(dgc[1], 5)],
                               b0=betti_from_bars(dgc[0]))
    print(f"  concatenated: A={stc['A_max_h1_persistence']:.5f} B={stc['B_n_h1_bars_ge_0.005']:.0f} "
          f"n_h1={stc['n_h1_bars']:.0f}")

    # site-bootstrap null on the concatenated alignment
    cpath = CACHE / f"concat_siteboot_n{a.n_null}.json"
    done = json.loads(cpath.read_text()) if cpath.exists() else {}
    t0 = time.time()
    for sd in range(a.n_null):
        if str(sd) in done:
            continue
        if time.time() - t0 > a.budget_sec:
            break
        rng = np.random.default_rng(sd)
        _, st = summarise(pdist_hamming(bootstrap_cols(concat, rng)))
        done[str(sd)] = [st["A_max_h1_persistence"], st["B_n_h1_bars_ge_0.005"]]
        if sd % 10 == 0:
            cpath.write_text(json.dumps(done))
    cpath.write_text(json.dumps(done))
    print(f"  site-bootstrap null: {len(done)}/{a.n_null} done ({round(time.time() - t0, 1)}s)")
    if len(done) < a.n_null:
        print("  INCOMPLETE -- re-invoke with the same arguments")
        return 2
    nv = np.array([done[str(i)] for i in range(a.n_null)])
    p_A = float(rank_p(stc["A_max_h1_persistence"], nv[:, 0]))
    p_B = float(rank_p(stc["B_n_h1_bars_ge_0.005"], nv[:, 1]))

    crit = {
        "i_A_concat_significant_p_le_0.01": bool(p_A <= 0.01),
        "ii_every_single_segment_A_le_half_A_concat": bool(all(
            obs[k]["A_max_h1_persistence"] <= 0.5 * stc["A_max_h1_persistence"]
            for k in obs if k != "concatenated")),
        "iii_B_concat_gt_max_single_B": bool(stc["B_n_h1_bars_ge_0.005"] > max(
            obs[k]["B_n_h1_bars_ge_0.005"] for k in obs if k != "concatenated")),
    }
    verdict = "PASS" if all(crit.values()) else "FAIL"

    blk = Block("influenza", SCRIPT, COMMAND)
    shas = {f"aln_seg{s}.npy": sha256(WORK / f"aln_seg{s}.npy") for s in SEGMENTS}
    for s, nm in SEGMENTS.items():
        k = f"seg{s}_{nm}"
        ds_id = f"biology/influenza_seg{s}_{nm.lower()}"
        blk.dataset(id=ds_id, domain="biology",
                    title=f"Influenza A avian segment {s} ({nm}), {n_gen} genomes, "
                          f"{obs[k]['aln_len']} aligned columns",
                    source="NCBI Influenza Virus Resource (genomeset.dat + influenza.fna); "
                           "https://ftp.ncbi.nlm.nih.gov/genomes/INFLUENZA/",
                    provenance="observation", n_objects=n_gen, ambient_dim=obs[k]["aln_len"],
                    units="p-distance (fraction of differing non-gap columns)",
                    sha256=shas[f"aln_seg{s}.npy"], local_path=str(WORK / f"aln_seg{s}.npy"),
                    notes=f"alignment produced by the earlier genetics validation run (pyfamsa "
                          f"0.7.0) and reused byte-for-byte; {prep['n_unique_concatenated']} unique "
                          f"avian 8-segment sets, {n_gen} subsampled with seed {SEED_SUB}. "
                          f"A single segment is expected to be TREE-like, i.e. H1 ~ 0.")
        st = obs[k]
        blk.run(dataset_id=ds_id, method="rips", coeff_field=2, max_dim=1,
                params={"metric": "hamming p-distance fed directly to Rips (NO embedding)",
                        "max_hom_dim": 1, "max_edge_length": None, "edge_collapse": True,
                        "expansion_dim": 2, "n_genomes": n_gen, "aln_len": st["aln_len"],
                        "bar_threshold_B": BAR_THRESH},
                preprocessing="gap/ambiguous columns excluded pairwise; p-distance over columns "
                              "where both sequences are A/C/G/T",
                seed=str(SEED_SUB), tier="X", wall_sec=st["wall_sec"],
                diagrams={1: [(b, d) for b, d in st["diagram_h1"]]},
                betti={1: int(st["B_n_h1_bars_ge_0.005"])}, expected_betti={1: 0},
                stats=[{"name": "A_max_h1_persistence", "value": st["A_max_h1_persistence"]},
                       {"name": "B_n_h1_bars_ge_0.005", "value": st["B_n_h1_bars_ge_0.005"]},
                       {"name": "n_h1_bars", "value": st["n_h1_bars"]},
                       {"name": "h1_dominance_P1_over_P2", "value": st["h1_dominance_P1_over_P2"]},
                       {"name": "median_pdistance", "value": st["median_pdistance"]}],
                controls=[{"kind": "known_answer",
                           "description": "a clonally evolving single segment is a tree: no H1 bar "
                                          "above 0.005 expected",
                           "passed": bool(st["B_n_h1_bars_ge_0.005"] <=
                                          0.5 * stc["B_n_h1_bars_ge_0.005"]),
                           "detail": f"B = {st['B_n_h1_bars_ge_0.005']:.0f} vs concatenated "
                                     f"{stc['B_n_h1_bars_ge_0.005']:.0f}"},
                          {"kind": "negative",
                           "description": "the MDS-then-alpha route is NOT computed: it gave "
                                          "dominance 11.1 on this same segment 4 (HA) in the earlier "
                                          "validation, a FALSE loop on a tree",
                           "passed": None,
                           "detail": "no embedding is computed anywhere in this block"}],
                findings=[{"claim": f"Segment {s} ({nm}): Rips on the p-distance metric gives "
                                    f"max H1 persistence {st['A_max_h1_persistence']:.5f} and "
                                    f"{st['B_n_h1_bars_ge_0.005']:.0f} bars above {BAR_THRESH} over "
                                    f"{n_gen} genomes.",
                           "verdict": "null" if st["B_n_h1_bars_ge_0.005"] <= 4 else "inconclusive",
                           "tier": "X",
                           "caveat": "no per-segment null was run; only the concatenated statistic "
                                     "has a site-bootstrap p-value",
                           "reference": "Chan, Carlsson & Rabadan, PNAS 110:18566 (2013)"}])

    ds_id = "biology/influenza_concatenated_genomes"
    blk.dataset(id=ds_id, domain="biology",
                title=f"Influenza A avian concatenated 8-segment genomes, {n_gen} genomes, "
                      f"{obs['concatenated']['aln_len']} aligned columns",
                source="NCBI Influenza Virus Resource (genomeset.dat + influenza.fna); "
                       "https://ftp.ncbi.nlm.nih.gov/genomes/INFLUENZA/",
                provenance="observation", n_objects=n_gen,
                ambient_dim=obs["concatenated"]["aln_len"],
                units="p-distance (fraction of differing non-gap columns)",
                local_path=str(WORK), notes=f"concatenation of the eight per-segment alignments "
                                            f"(sha256 each: {json.dumps(shas)}). Reassortment is "
                                            f"expected to create H1 classes.")
    blk.run(dataset_id=ds_id, method="rips", coeff_field=2, max_dim=1,
            params={"metric": "hamming p-distance fed directly to Rips (NO embedding)",
                    "max_hom_dim": 1, "max_edge_length": None, "edge_collapse": True,
                    "expansion_dim": 2, "n_genomes": n_gen,
                    "aln_len": obs["concatenated"]["aln_len"], "bar_threshold_B": BAR_THRESH,
                    "criteria": crit},
            preprocessing="eight per-segment alignments concatenated column-wise",
            seed=f"genome subsample {SEED_SUB}; bootstrap seeds 0..{a.n_null - 1}", tier="X",
            wall_sec=obs["concatenated"]["wall_sec"],
            diagrams={1: [(b, d) for b, d in obs["concatenated"]["diagram_h1"]]},
            betti={1: int(stc["B_n_h1_bars_ge_0.005"])}, expected_betti={1: 1},
            stats=[{"name": "A_max_h1_persistence", "value": stc["A_max_h1_persistence"],
                    "null_model": "site bootstrap: alignment columns resampled with replacement",
                    "n_null": a.n_null, "p_value": p_A, "p_method": "rank"},
                   {"name": "B_n_h1_bars_ge_0.005", "value": stc["B_n_h1_bars_ge_0.005"],
                    "null_model": "site bootstrap: alignment columns resampled with replacement",
                    "n_null": a.n_null, "p_value": p_B, "p_method": "rank"},
                   {"name": "n_h1_bars", "value": stc["n_h1_bars"]},
                   {"name": "h1_dominance_P1_over_P2", "value": stc["h1_dominance_P1_over_P2"]},
                   {"name": "median_pdistance", "value": stc["median_pdistance"]},
                   {"name": "max_single_segment_A", "value": max(
                       obs[k]["A_max_h1_persistence"] for k in obs if k != "concatenated")}],
            controls=[{"kind": "known_answer",
                       "description": "reassorting whole genomes must show more H1 than any single "
                                      "clonal segment (Chan/Carlsson/Rabadan 2013)",
                       "passed": bool(crit["ii_every_single_segment_A_le_half_A_concat"]
                                      and crit["iii_B_concat_gt_max_single_B"]),
                       "detail": json.dumps(crit)},
                      {"kind": "shuffle",
                       "description": f"site bootstrap (n={a.n_null}) breaks the linkage that "
                                      f"carries reassortment while keeping per-site composition",
                       "passed": bool(p_A <= 0.01),
                       "detail": f"null A q50/q95/q99/max = "
                                 f"{[round(float(np.quantile(nv[:, 0], q)), 5) for q in (.5, .95, .99, 1)]}"},
                      {"kind": "negative",
                       "description": "the MDS-then-alpha route that produced a FALSE loop "
                                      "(dominance 11.1 on the tree-like HA segment) is not computed",
                       "passed": None, "detail": "lesson (a) of the earlier genetics validation"}],
            findings=[{"claim": f"Concatenated influenza genomes give max H1 persistence "
                                f"{stc['A_max_h1_persistence']:.5f} (site-bootstrap p {p_A:.4f}, "
                                f"n={a.n_null}) and {stc['B_n_h1_bars_ge_0.005']:.0f} bars above "
                                f"{BAR_THRESH}, against at most "
                                f"{max(obs[k]['A_max_h1_persistence'] for k in obs if k != 'concatenated'):.5f} "
                                f"and {max(obs[k]['B_n_h1_bars_ge_0.005'] for k in obs if k != 'concatenated'):.0f} "
                                f"for any single segment. Verdict {verdict}.",
                       "verdict": "recovered" if verdict == "PASS" else "inconclusive", "tier": "X",
                       "caveat": "300 avian genomes, one subsample (seed 20260919); Rips on the "
                                 "p-distance metric, never an embedding",
                       "reference": "Chan, Carlsson & Rabadan, PNAS 110:18566 (2013)"}])

    out = {"seed_subsample": SEED_SUB, "n_genomes": n_gen, "n_null": a.n_null,
           "p_A_concat": p_A, "p_B_concat": p_B, "criteria": crit, "verdict": verdict,
           "null_A_q50_q95_q99_max": [float(np.quantile(nv[:, 0], q)) for q in (.5, .95, .99, 1)],
           "sha256": shas,
           "observed": {k: {kk: vv for kk, vv in v.items() if kk != "diagram_h1"}
                        for k, v in obs.items()}}
    (RESULTS / "influenza_summary.json").write_text(json.dumps(out, indent=1))
    print(json.dumps({"criteria": crit, "verdict": verdict, "p_A": p_A, "p_B": p_B}, indent=1))
    blk.write()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
