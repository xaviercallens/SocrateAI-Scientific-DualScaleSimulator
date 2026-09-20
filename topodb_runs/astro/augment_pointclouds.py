#!/usr/bin/env python3
"""Add the separation statistics the rank p-value cannot express.

WHY THIS EXISTS.  With n_null = 20 the smallest attainable TWO-SIDED rank
p-value is 2/21 = 0.0952.  That is above the pre-declared Bonferroni/6
threshold of 0.00833, so no statistic can ever pass it, however far the data
lies from the null.  On every DESI sample measured here the data sits OUTSIDE
the entire range of the 20 null draws -- e.g. BGS_BRIGHT-21.5 NGC has 114 long
H2 bars against 208.1 +/- 10.9 in the randoms -- and the p-value at its floor
says nothing about that.

Reporting only the floored p-value would therefore misrepresent a large,
unambiguous separation as a null result.  Reporting a small p-value that 20
draws cannot support would be worse.  The honest record is both:

  * the rank p AT ITS FLOOR, already stored, with the floor in `p_method`;
  * `null_separation_sigma` = (data - mean(null)) / std(null), stored WITHOUT a
    p-value, because a z-score from 20 draws has no calibrated null here and
    TopoDB correctly refuses a p-value without one;
  * `data_outside_null_range`, 1 or 0.

None of this changes what the comparison means: the DESI randoms are
UNCLUSTERED, so a large separation was expected in advance and is a measurement
of the pipeline's sensitivity to 3-D clustering, not of any anomaly.

Tier X.
"""
from __future__ import annotations

import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RESULTS, topodb  # noqa: E402

STATS = ["n_h1_bars_over_5mpc", "n_h2_bars_over_5mpc", "total_persistence_h1",
         "total_persistence_h2", "max_persistence_h2", "h0_death_iqr_over_median"]


def main():
    n = 0
    for f in sorted(glob.glob(os.path.join(RESULTS, "pointclouds_*.json"))):
        d = json.load(open(f))
        for ds, v in d.items():
            if not isinstance(v, dict) or "null" not in v or not v["null"]:
                continue
            rid = v["run_id"]
            outside = {}
            with topodb() as db:
                for s in STATS:
                    nv = np.array([x[s] for x in v["null"]], float)
                    obs = float(v["data"][s])
                    sd = float(nv.std(ddof=1))
                    z = (obs - float(nv.mean())) / sd if sd > 0 else float("nan")
                    out = bool(obs < nv.min() or obs > nv.max())
                    outside[s] = out
                    if np.isfinite(z):
                        db.add_statistic(rid, f"null_separation_sigma__{s}", float(z))
                    db.add_statistic(rid, f"data_outside_null_range__{s}", float(out))
                k = sum(outside.values())
                db.add_control(rid, "null_calibration",
                               "rank-p resolution: with n_null=20 the two-sided floor is 2/21=0.0952, "
                               "ABOVE the pre-declared Bonferroni/6 threshold 0.00833, so no statistic "
                               "can pass it however far the data lies from the null; the separation is "
                               "therefore also recorded in null-sigma units, without a p-value",
                               passed=None,
                               detail=(f"{k}/{len(STATS)} statistics lie outside the FULL range of the "
                                       f"20 null draws; "
                                       + ", ".join(f"{s}: z={((float(v['data'][s]) - np.mean([x[s] for x in v['null']])) / max(np.std([x[s] for x in v['null']], ddof=1), 1e-12)):+.1f}"
                                                   for s in STATS)))
            print(f"{ds} (run {rid}): {sum(outside.values())}/{len(STATS)} outside the null range")
            n += 1
    print(f"augmented {n} runs")
    return 0


if __name__ == "__main__":
    sys.exit(main())
