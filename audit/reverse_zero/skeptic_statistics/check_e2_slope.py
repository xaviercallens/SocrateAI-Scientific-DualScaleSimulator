#!/usr/bin/env python3
"""Skeptic check of E2's sensitivity-scan slope: E2_RESULT.md:78 and the returned numbers block say
slope = 6.171; recompute from scripts/param_loop_sim.compute_pta_observable and compare with the committed
e2_pta_result.json. No RNG. Command (worktree root):
  /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python audit/reverse_zero/skeptic_statistics/check_e2_slope.py
"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from param_loop_sim import compute_pta_observable  # noqa: E402
ref = compute_pta_observable(1.0, 0.0); hd = np.array(ref["hd_curve"]); th = np.degrees(ref["theta_rad"])
res = compute_pta_observable(1.0, 0.1)
dev = np.abs(np.array(res["gamma_theta"]) - hd)
committed = json.load(open(os.path.join(HERE, "..", "e2_pta_and_tda_extension", "e2_pta_result.json")))
md = open(os.path.join(HERE, "..", "e2_pta_and_tda_extension", "E2_RESULT.md")).read().splitlines()
out = {"slope_recomputed": res["max_deviation_from_hd"] / 0.1,
       "slope_committed_json": committed["sensitivity_scan"]["linearity_check"]["slope_max_dev_per_unit_c4"],
       "E2_RESULT_md_line78": md[77] if len(md) > 77 else None,
       "hd_dynamic_range": float(hd.max() - hd.min()),
       "argmax_deviation_theta_deg": [float(th[i]) for i in np.where(dev >= dev.max() - 1e-12)[0]],
       "threshold_1pct_recomputed": 0.01 * float(hd.max() - hd.min()) / (res["max_deviation_from_hd"] / 0.1)}
json.dump(out, open(os.path.join(HERE, "check_e2_slope.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
