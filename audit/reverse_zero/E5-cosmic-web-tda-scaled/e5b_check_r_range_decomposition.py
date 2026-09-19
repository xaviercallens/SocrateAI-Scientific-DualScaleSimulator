#!/usr/bin/env python3
"""
Post-hoc advisor-review check on e5b_lognormal_and_poisson_null.py's
final (third) full run. Result folded into
e5b_lognormal_and_poisson_null_report.json's
post_hoc_advisor_review_checks.r_range_decomposition_of_the_betti_excess.

WHY: the main run's xi(r) VALIDATION GATE only covers r in [20,60]
Mpc/h (XI_FIT_RANGE_MPC_H) -- that is the ONLY range in which the
lognormal mocks are shown to reproduce the data's own two-point
function. The Betti-curve comparison (STAT_R_RANGE) runs over (2,40)
Mpc/h, so HALF of that range, r in [2,20), is claimed as "diagnostic"
by the gate's own stated rule ("this gate MUST pass for the Betti-curve
comparison to be treated as diagnostic of anything") without actually
being validated there. Separately, the lognormal field is Gaussian-
smoothed at R_SMOOTH=1 grid cell=2.0 Mpc/h (the fix for the sigma_G^2
blowup documented elsewhere in this report) -- the smoothing window is
at half power around r~5 Mpc/h, so the mocks are close to POISSON
(uniform-within-cell) below that scale by construction, while real
galaxies have genuine small-scale clustering there. Both point at r<20
Mpc/h as a place a Betti-curve "excess" could be a MOCK ARTIFACT rather
than a detection.

TEST: split the SAME pre-stated L2 statistic (leave-one-out rank test,
identical definition and code to the main run's) into r in [2,20) and
r in [20,40] Mpc/h -- the gate's own validated/unvalidated boundary --
using the mock Betti-curve stacks already saved by the main run
(e5b_mock_betti_stacks.npz), not a new run (same mocks, same seeds).

RESULT (2026-09-19): the ENTIRE excess is in r in [2,20). At r in
[20,40] (the gate-validated range), rank=22/22 (the real data is LESS
extreme than every single leave-one-out mock) and the margin is
undefined (median mock-to-mock L2 in that range is exactly 0 for all
three homology dimensions -- i.e. essentially no Betti-curve variation
at all, mock or real, above r=20 Mpc/h at this N and truncation). At r
in [2,20), the margins reproduce the full-range numbers exactly (33.8x,
17.9x, 10.2x for H0/H1/H2) -- confirming ALL of the reported "excess"
lives below the range the xi(r) gate actually validates, and inside/
near the mock's own Gaussian smoothing scale.

CONCLUSION: this is a LIMITATION of the lognormal-mock null below ~20
Mpc/h (a range the mocks are not shown to be valid in, and are smoothed
to be nearly featureless in by construction), NOT evidence that real
cosmic-web topology exceeds a properly-resolved null -- the properly-
resolved (gate-validated) range shows NO excess at all. Superseded any
earlier framing of this result as a "detection of structure beyond the
two-point function."

Command:
  /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python \
    audit/reverse_zero/E5-cosmic-web-tda-scaled/e5b_check_r_range_decomposition.py
"""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
NPZ_PATH = os.path.join(HERE, "e5b_mock_betti_stacks.npz")
GATE_VALIDATED_RANGE = (20.0, 40.0)  # intersection of XI_FIT_RANGE_MPC_H=(20,60) and STAT_R_RANGE=(2,40)
UNVALIDATED_RANGE = (2.0, 20.0)


def l2(a, b, r, r_range):
    mask = (r >= r_range[0]) & (r <= r_range[1])
    return float(np.sqrt(np.mean((a[mask] - b[mask]) ** 2)))


def run():
    z = np.load(NPZ_PATH)
    r = z["r_grid"]
    out = {}
    for key, rkey, label in [("b0", "b0_real", "H0"), ("b1", "b1_real", "H1"), ("b2", "b2_real", "H2")]:
        stack, real = z[key], z[rkey]
        n = stack.shape[0]
        entry = {}
        for name, r_range in [("unvalidated_2_20", UNVALIDATED_RANGE), ("gate_validated_20_40", GATE_VALIDATED_RANGE)]:
            loo = np.array([l2(stack[i], np.delete(stack, i, axis=0).mean(axis=0), r, r_range) for i in range(n)])
            real_l2 = l2(real, stack.mean(axis=0), r, r_range)
            rank = int(np.sum(loo >= real_l2))
            median_loo = float(np.median(loo))
            margin = float(real_l2 / median_loo) if median_loo > 0 else None
            entry[name] = {"r_range": list(r_range), "real_l2": real_l2, "rank": rank, "n_mocks": n,
                           "median_loo_l2": median_loo, "margin": margin,
                           "p_value": (rank + 1) / (n + 1)}
            print(label, name, entry[name])
        out[label] = entry
    with open(os.path.join(HERE, "e5b_r_range_decomposition.json"), "w") as f:
        json.dump(out, f, indent=2)
    return out


if __name__ == "__main__":
    run()
