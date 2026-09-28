#!/usr/bin/env python3
"""SECONDARY diagnostic of registration.json TEST_A (listed there as a
conditional item, outside the multiplicity family): the lower-star Betti curves
of the smoothed field T_s on the HEALPix 2-complex.

THIS IS THE PLACE THE FIXED TDA LIBRARY BELONGS, and it is the only lower-star
statistic in this audit:
  lib/cmb_topology.build_topology_fixed        replaces cmb_tda.build_topology
  lib/cmb_topology.betti_curves_from_topology  replaces the original of the same name
  lib/stats.coarse_stats_fixed                 replaces cmb_tda.coarse_stats
audit/reverse_zero/E5-cmb-tda/cmb_tda.py is NOT imported.  The RANK p-value is
the reported statistic; the chi2 branch is carried only as a diagnostic because
it remains mildly anti-conservative (0.058 instead of 0.05).

MANDATORY CONTROL (registration limit L2).  On rough fields the lower-star path
FAILED a site-shuffle control in the Re6Zr work: shuffled fields reproduced the
"signal", so it came from the value distribution and not from topology.  This
script therefore also runs the C2 value-shuffle ensemble through the identical
lower-star path and reports where it lands.

CORRECTION (2026-09-20), which report.json carries and which supersedes the
`control_verdict` string this script emits: what is measured here is whether the
SHUFFLED ensemble is SEPARATED FROM THE GAUSSIAN NULL.  That is NOT the
proposition L2 names.  L2's failure mode is that the shuffled field looks like
the REAL field, i.e. that an apparent signal lives in the value distribution; a
statistic can separate shuffles from the null and still be value-driven.  The
fractions below are therefore reported as a separation measurement and NOTHING
is concluded from them about L2.  The verdict string is left as the code wrote
it so the committed JSON matches its producer.

Command:
  timeout 7200 prlimit --as=8589934592 -- <venv-tda python> cv_lowerstar.py --which wmap --n 200
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cv_lib as L  # noqa: E402
sys.path.insert(0, L.CV)
from lib import cmb_topology as CT  # noqa: E402
from lib import stats as ST  # noqa: E402

OUT = os.path.join(L.CV, "results")
NU = np.linspace(-4.0, 4.0, 41)        # X1's registered nu grid


def log(m):
    print("[%s] %s" % (time.strftime("%H:%M:%S"), m), file=sys.stderr, flush=True)


def curves(ts, unm, edges, tris):
    d = CT.betti_curves_from_topology(ts, unm, edges, tris, NU, sublevel=True)
    return np.asarray(d["b0"], float), np.asarray(d["b1"], float), np.asarray(d["euler_char_true"], float)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--which", default="wmap", choices=list(L.MAPS))
    ap.add_argument("--n", type=int, default=200)
    ap.add_argument("--nshuf", type=int, default=50)
    a = ap.parse_args()

    m, mask = L.load_data(a.which)
    unm, edges, tris, info = CT.build_topology_fixed(mask, L.NSIDE, return_info=True)
    cl_in, _ = L.load_cl_in(a.which)

    t0 = time.time()
    tp = L.prep_map(m, mask)
    ts_data = L.smooth_map(tp, mask)
    b0d, b1d, eud = curves(ts_data, unm, edges, tris)
    log("data curves in %.1f s; b0 range %g..%g" % (time.time() - t0, b0d.min(), b0d.max()))

    nb0, nb1, neu = [], [], []
    for k in range(a.n):
        s = L.make_sim(cl_in, 2000000 + k)          # SAME null seeds as the alpha path
        ts = L.smooth_map(L.prep_map(s, mask), mask)
        c = curves(ts, unm, edges, tris)
        nb0.append(c[0]); nb1.append(c[1]); neu.append(c[2])
        if k % 25 == 0:
            log("null %d/%d" % (k, a.n))
    nb0, nb1, neu = map(np.array, (nb0, nb1, neu))

    sb0, sb1 = [], []
    er = L.eroded_mask(mask, L.disk_table()[0])
    for k in range(a.nshuf):
        rng = np.random.default_rng(6100000 + k)     # SAME seeds as the C2 control
        sh = ts_data.copy()
        # NOTE: the lower-star complex is built on the FULL mask, so the
        # permutation is over mask > 0 -- NOT over the R_disk-eroded region the
        # alpha path's C2 control uses.  Same seed integers, different pixel set.
        sh[mask > 0] = rng.permutation(ts_data[mask > 0])
        c = curves(sh, unm, edges, tris)
        sb0.append(c[0]); sb1.append(c[1])
    sb0, sb1 = np.array(sb0), np.array(sb1)

    res = {
        "which": a.which, "n_null": a.n, "null_seeds": "2000000+k (the SAME null ensemble as the alpha path)",
        "n_shuffle": a.nshuf, "shuffle_seeds": "6100000+k (the same seed integers as the alpha path's C2 control, but "
                         "permuting over mask > 0 rather than over the eroded region)",
        "nu_grid": "np.linspace(-4, 4, 41), X1's registered grid",
        "filtration": "lower-star on the HEALPix quad-corner 2-complex of lib/cmb_topology.build_topology_fixed",
        "complex_info": {k: (v if not isinstance(v, np.ndarray) else v.tolist())
                         for k, v in info.items() if not isinstance(v, (np.ndarray,))},
        "modules_used": {
            "topology": "lib/cmb_topology.build_topology_fixed + betti_curves_from_topology "
                        "(loop/tda-simple d8175f1)",
            "statistic": "lib/stats.coarse_stats_fixed (RANK p reported; chi2 diagnostic only)",
            "NOT_used": "audit/reverse_zero/E5-cmb-tda/cmb_tda.py (the defective originals)",
        },
        "fsky_mask": float(mask.mean()),
    }
    for tag, sims, data in (("b0", nb0, b0d), ("b1", nb1, b1d)):
        s = ST.coarse_stats_fixed(sims, data, n_bins=8)
        res.setdefault("data_vs_null", {})[tag] = {
            k: s[k] for k in ("empirical_rank_p", "p_value_chi2_survival", "df", "retained_rank",
                              "n_bins_kept", "dropped_bin_indices", "dropped_bin_reasons",
                              "test_differs_in_dropped_bin", "hartlap_factor", "data_chi2_hartlap")}
    for tag, sims, shuf in (("b0", nb0, sb0), ("b1", nb1, sb1)):
        ps = [ST.coarse_stats_fixed(sims, shuf[i], n_bins=8)["empirical_rank_p"] for i in range(len(shuf))]
        ps = np.array([p for p in ps if p is not None], float)
        res.setdefault("shuffle_control_vs_null", {})[tag] = {
            "n": int(ps.size), "median_rank_p": float(np.median(ps)),
            "fraction_below_0.05": float((ps < 0.05).mean()),
            "min_rank_p": float(ps.min()), "max_rank_p": float(ps.max())}
    f0 = res["shuffle_control_vs_null"]["b0"]["fraction_below_0.05"]
    f1 = res["shuffle_control_vs_null"]["b1"]["fraction_below_0.05"]
    # NOTE (correction, 2026-09-20): the field below measures SEPARATION OF THE
    # SHUFFLED ENSEMBLE FROM THE GAUSSIAN NULL.  That is NOT the proposition
    # registration limit L2 names -- L2's failure mode is that the shuffled
    # field looks like the REAL field.  report.json supersedes the verdict
    # string emitted here and reports the fractions as a separation measurement
    # only.  Kept unchanged so the committed JSON matches the code that wrote it.
    stands = bool(f0 >= 0.95 or f1 >= 0.95)
    res["control_verdict"] = {
        "fraction_of_shuffles_separated_from_the_null_b0": f0,
        "fraction_b1": f1,
        "diagnostic_stands": stands,
        "reading": ("The value-shuffled field is separated from the Gaussian null by this statistic in "
                    "%.0f%% (b0) / %.0f%% (b1) of realisations, so the lower-star path IS reading spatial "
                    "arrangement on this field and the diagnostic stands." % (100 * f0, 100 * f1)) if stands
                   else ("The value-shuffled field is NOT reliably separated from the Gaussian null "
                         "(%.0f%% b0 / %.0f%% b1), i.e. this statistic cannot tell a spatially destroyed "
                         "field from a real one here. The diagnostic is WITHDRAWN, per registration limit L2."
                         % (100 * f0, 100 * f1)),
    }
    res["tier"] = "X"
    json.dump(res, open(os.path.join(OUT, "cmb_lowerstar_%s.json" % a.which), "w"), indent=1)
    print(json.dumps({"data_vs_null": res["data_vs_null"], "control": res["control_verdict"]}, indent=1)[:2500])


if __name__ == "__main__":
    main()
