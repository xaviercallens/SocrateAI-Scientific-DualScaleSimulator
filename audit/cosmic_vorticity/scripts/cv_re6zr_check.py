#!/usr/bin/env python3
"""KNOWN-ANSWER CHECK, run before any CV number is interpreted.

LeanFlow CLAUDE.md rule 4 ("Known answers before new answers"): the statistic
this project transfers to the sky must first reproduce, on the SAME saved point
clouds, the numbers the Re6Zr run published.  Two things are checked:

  (1) cv_lib.alpha_h0_deaths + cv_lib.spread_stats reproduce
      audit/tda_validation/quantum_fluid/results/stm_vortex_tda.json's
      h0_pooled.iqr_over_median for all 11 fields, from the saved vortex-core
      point clouds (DATA_ROOT/stm/minima_*kOe_460mK.npz).

  (2) TRUNCATION INERTNESS.  The registered CV statistic keeps Re6Zr's
      max_alpha_square = (3a)^2.  Truncation removes the upper tail of the H0
      deaths and so can bias IQR/median.  S1 is therefore recomputed with NO
      truncation (max_alpha_square = inf) and the two values compared.  If they
      differ materially on the source data, the transferred truncation is not
      inert and that is disclosed rather than assumed.

Command:
  timeout 900 prlimit --as=8589934592 -- <venv-tda python> cv_re6zr_check.py
"""
import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cv_lib as L  # noqa: E402

STM = "/mnt/disks/disk-socrateai-local-1/dualscale-data-r3/tda_validation/quantum_fluid/stm"
PUB = ("/mnt/disks/disk-socrateai-local-1/dualscale-wt-tdaval/audit/tda_validation/"
       "quantum_fluid/results/stm_vortex_tda.json")
PHI0 = 2.067833848e-15      # stm_vortex_tda.py, verbatim


def main():
    pub = json.load(open(PUB))
    rows = []
    for fn in sorted(glob.glob(os.path.join(STM, "minima_*kOe_460mK.npz")),
                     key=lambda p: float(os.path.basename(p).split("_")[1][:-3])):
        H = float(os.path.basename(fn).split("_")[1][:-3])
        a = 1.075 * np.sqrt(PHI0 / (H / 10.0)) * 1e9        # a_tri(B) in nm, verbatim
        d = np.load(fn)
        trunc, full = [], []
        for k in sorted(d.files):
            pts = d[k]
            trunc.append(L.alpha_h0_deaths(pts, (3 * a) ** 2)[0])
            full.append(L.alpha_h0_deaths(pts, float("inf"))[0])
        st_t = L.spread_stats(np.concatenate(trunc))
        st_f = L.spread_stats(np.concatenate(full))
        key = "%gkOe" % H
        ref = pub["fields"][key]["h0_pooled"]
        rows.append({
            "field": key, "a_tri_nm": float(a), "n_images": len(d.files),
            "published_iqr_over_median": ref["iqr_over_median"],
            "recomputed_truncated_3a": st_t["iqr_over_median"],
            "abs_diff_vs_published": abs(st_t["iqr_over_median"] - ref["iqr_over_median"]),
            "published_n_bars": ref["n"], "recomputed_n_bars_truncated": st_t["n"],
            "recomputed_untruncated": st_f["iqr_over_median"],
            "n_bars_untruncated": st_f["n"],
            "truncation_shift_in_S1": st_f["iqr_over_median"] - st_t["iqr_over_median"],
            "published_median_nm": ref["median"], "recomputed_median_nm": st_t["median"],
        })
        print("%-7s published %.6f  recomputed %.6f  diff %.2e | untruncated %.6f (shift %.2e) "
              "| bars %d vs %d" % (key, ref["iqr_over_median"], st_t["iqr_over_median"],
                                   rows[-1]["abs_diff_vs_published"], st_f["iqr_over_median"],
                                   rows[-1]["truncation_shift_in_S1"], ref["n"], st_t["n"]), flush=True)
    md = max(r["abs_diff_vs_published"] for r in rows)
    ms = max(abs(r["truncation_shift_in_S1"]) for r in rows)
    out = {
        "purpose": "known-answer check of the transferred S1 statistic + truncation-inertness test",
        "source_point_clouds": STM,
        "published_reference": PUB,
        "statistic": "cv_lib.spread_stats(cv_lib.alpha_h0_deaths(pts, (3*a_tri)^2))['iqr_over_median']",
        "rows": rows,
        "max_abs_diff_vs_published": md,
        "reproduces_published": bool(md < 1e-12),
        "max_abs_truncation_shift_in_S1": ms,
        "truncation_inert_on_source_data": bool(ms < 1e-3),
        "tier": "X",
    }
    json.dump(out, open(os.path.join(L.CV, "results", "re6zr_known_answer_check.json"), "w"), indent=1)
    print("max |diff vs published| = %.3e  reproduces = %s" % (md, out["reproduces_published"]))
    print("max |truncation shift|  = %.3e  inert = %s" % (ms, out["truncation_inert_on_source_data"]))


if __name__ == "__main__":
    main()
