"""POST-HOC addendum: a null for influenza reassortment that is actually a null.

Why this exists.  The pre-stated criterion (i) in expectations.json used a site
bootstrap of the concatenated alignment, and it gave p = 0.756.  That is not a weak
signal; it is an ill-posed test.  Resampling alignment COLUMNS leaves every genome's
row intact, so the pattern of which genome carries which segment -- the reassortment
itself -- survives the resampling untouched.  A site bootstrap measures how stable
A is under site resampling; it cannot ask whether reassortment happened.  The
original run and its FAIL verdict are kept in TopoDB exactly as computed; this adds
the test that the question needed, and it is labelled post_hoc everywhere.

The null used here is the one the earlier genetics validation used: a POOLED
SINGLE-SEGMENT site bootstrap.  Each single segment is clonal, so its genealogy is
a tree and its H1 is the H1 of tree-like data at this sample size.  Bootstrapping
sites within each segment and pooling the resulting A values across all eight gives
the distribution of A for NON-reassorting data, which is the right comparison for
A_concatenated.

Its limitation is stated rather than hidden: the single segments are shorter than
the concatenation (865-2369 vs 14008 columns), so this null differs from the
observed statistic in alignment length as well as in reassortment.  A longer
non-reassorting alignment of the same organisms does not exist to draw from.

Usage:  python run_influenza_null2.py [--n-per-segment 25]
Seeds: bootstrap seeds 100000 + 1000*segment + draw.
"""
from __future__ import annotations

import argparse
import json
import sys
import time

import numpy as np

sys.path.insert(0, "/mnt/disks/disk-socrateai-local-1/wt-topo-bio/topodb_runs/biology")
from bio_common import RESULTS, Block, rank_p, sha256  # noqa: E402
from run_influenza import (BAR_THRESH, SEGMENTS, WORK, bootstrap_cols, load_seg,  # noqa: E402
                           pdist_hamming, summarise)

SCRIPT = "topodb_runs/biology/run_influenza_null2.py"
COMMAND = ("timeout 1200 prlimit --as=8589934592 -- "
           "/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python "
           "topodb_runs/biology/run_influenza_null2.py --n-per-segment 25")
CACHE = RESULTS / "flu_cache"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-per-segment", type=int, default=25)
    ap.add_argument("--budget-sec", type=float, default=1000)
    a = ap.parse_args()
    CACHE.mkdir(parents=True, exist_ok=True)

    prev = json.loads((RESULTS / "influenza_summary.json").read_text())
    obs = prev["observed"]
    A_obs = obs["concatenated"]["A_max_h1_persistence"]
    B_obs = obs["concatenated"]["B_n_h1_bars_ge_0.005"]

    cpath = CACHE / f"pooled_single_segment_null_{a.n_per_segment}per.json"
    done = json.loads(cpath.read_text()) if cpath.exists() else {}
    t0 = time.time()
    for s in SEGMENTS:
        A = load_seg(s)
        for d in range(a.n_per_segment):
            key = f"{s}_{d}"
            if key in done:
                continue
            if time.time() - t0 > a.budget_sec:
                break
            seed = 100000 + 1000 * s + d
            _, st = summarise(pdist_hamming(bootstrap_cols(A, np.random.default_rng(seed))))
            done[key] = [st["A_max_h1_persistence"], st["B_n_h1_bars_ge_0.005"]]
            cpath.write_text(json.dumps(done))
        print(f"  segment {s}: {sum(1 for k in done if k.startswith(f'{s}_'))}/"
              f"{a.n_per_segment} ({round(time.time() - t0, 1)}s)")
    n_total = 8 * a.n_per_segment
    if len(done) < n_total:
        print(f"  INCOMPLETE {len(done)}/{n_total} -- re-invoke")
        return 2

    nv = np.array(list(done.values()))
    p_A = float(rank_p(A_obs, nv[:, 0]))
    p_B = float(rank_p(B_obs, nv[:, 1]))
    crit = {"i_A_concat_p_le_0.01_vs_pooled_single_segment_bootstrap": bool(p_A <= 0.01),
            "ii_every_single_segment_A_le_half_A_concat":
                prev["criteria"]["ii_every_single_segment_A_le_half_A_concat"],
            "iii_B_concat_gt_max_single_B":
                prev["criteria"]["iii_B_concat_gt_max_single_B"]}
    verdict = "PASS" if all(crit.values()) else "FAIL"

    blk = Block("influenza_null2", SCRIPT, COMMAND)
    blk.run(dataset_id="biology/influenza_concatenated_genomes", method="rips", coeff_field=2,
            max_dim=1,
            params={"test": "POST-HOC: pooled single-segment site-bootstrap null for A_concatenated",
                    "metric": "hamming p-distance fed directly to Rips (NO embedding)",
                    "max_hom_dim": 1, "edge_collapse": True, "expansion_dim": 2,
                    "n_per_segment": a.n_per_segment, "n_null_total": n_total,
                    "prespecified": False,
                    "reason": "the pre-stated site bootstrap of the concatenated alignment cannot "
                              "be a null for reassortment: resampling columns leaves each genome's "
                              "row, and therefore its segment combination, intact",
                    "bar_threshold_B": BAR_THRESH},
            preprocessing="sites resampled with replacement WITHIN each single segment; the "
                          "resulting A and B values pooled across all eight segments",
            seed="100000 + 1000*segment + draw", tier="X",
            stats=[{"name": "A_max_h1_persistence_concat", "value": A_obs,
                    "null_model": "pooled single-segment site bootstrap (clonal, non-reassorting "
                                  "data at this sample size)", "n_null": n_total,
                    "p_value": p_A, "p_method": "rank",
                    "multiplicity": "post hoc; not pre-registered and not corrected"},
                   {"name": "B_n_h1_bars_ge_0.005_concat", "value": B_obs,
                    "null_model": "same pooled single-segment site bootstrap", "n_null": n_total,
                    "p_value": p_B, "p_method": "rank"},
                   {"name": "null_A_q50", "value": float(np.quantile(nv[:, 0], 0.5))},
                   {"name": "null_A_q95", "value": float(np.quantile(nv[:, 0], 0.95))},
                   {"name": "null_A_max", "value": float(nv[:, 0].max())},
                   {"name": "null_B_q95", "value": float(np.quantile(nv[:, 1], 0.95))}],
            controls=[{"kind": "null_calibration",
                       "description": "this null is drawn from CLONAL data (single segments are "
                                      "trees), so it is the distribution of A when no reassortment "
                                      "has occurred -- unlike the pre-stated site bootstrap of the "
                                      "concatenation, which preserves reassortment",
                       "passed": bool(p_A <= 0.01),
                       "detail": f"observed A {A_obs:.5f} vs null q50 "
                                 f"{np.quantile(nv[:, 0], 0.5):.5f}, q95 "
                                 f"{np.quantile(nv[:, 0], 0.95):.5f}, max {nv[:, 0].max():.5f}"},
                      {"kind": "negative",
                       "description": "DISCLOSED LIMITATION: the null's alignments are shorter "
                                      "(865-2369 columns) than the concatenation (14008), so this "
                                      "null differs from the observed statistic in length as well "
                                      "as in reassortment. No long non-reassorting alignment of "
                                      "these organisms exists to draw from.",
                       "passed": None,
                       "detail": "length confound stated, not corrected"}],
            findings=[{"claim": f"POST-HOC, with a null that can actually be one: A_concatenated = "
                                f"{A_obs:.5f} against a pooled single-segment site-bootstrap null "
                                f"of {n_total} draws (median {np.quantile(nv[:, 0], 0.5):.5f}, max "
                                f"{nv[:, 0].max():.5f}) gives p = {p_A:.5f}. With criteria (ii) and "
                                f"(iii) already met, the published contrast of Chan, Carlsson & "
                                f"Rabadan (2013) is {'reproduced' if verdict == 'PASS' else 'not reproduced'} "
                                f"on this criterion set: {verdict}.",
                       "verdict": "recovered" if verdict == "PASS" else "inconclusive", "tier": "X",
                       "caveat": "POST HOC: this null was chosen after the pre-stated one turned "
                                 "out to be ill-posed. The pre-stated run and its FAIL verdict "
                                 "remain in the database unchanged. The null's alignments are "
                                 "shorter than the concatenation, which is a stated confound.",
                       "reference": "Chan, Carlsson & Rabadan, PNAS 110:18566 (2013); the same null "
                                    "the 2026-09-19 genetics validation used (its p_A was 0.0025 "
                                    "over 400 draws)"}])
    out = {"post_hoc": True, "A_obs": A_obs, "B_obs": B_obs, "p_A": p_A, "p_B": p_B,
           "n_null_total": n_total, "criteria": crit, "verdict": verdict,
           "null_A_q50_q95_q99_max": [float(np.quantile(nv[:, 0], q)) for q in (.5, .95, .99, 1)],
           "null_B_q50_q95_max": [float(np.quantile(nv[:, 1], q)) for q in (.5, .95, 1)],
           "earlier_validation_p_A": 0.0024937655860349127}
    (RESULTS / "influenza_null2_summary.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))
    blk.write()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
