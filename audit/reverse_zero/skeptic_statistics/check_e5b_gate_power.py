#!/usr/bin/env python3
"""Skeptic (statistics lens) check of E5b's xi(r) validation gate: does the gate have power?
Reads the committed gate arrays from e5b_lognormal_and_poisson_null_report.json (lognormal_gate) and
re-evaluates the SAME pre-registered statistic (RMS over r in [20,60] Mpc/h of
(xi_mock - xi_data)/sigma_total, sigma_total = sqrt(sigma_mock^2 (1+1/n) + sigma_shot^2), PASS iff < 3.0,
GATE_Z_THRESHOLD read from the e5b script source) for counterfactual mocks:
  - xi_mock = 0 (an unclustered Poisson catalogue),
  - xi_mock = s * committed mock-mean xi, s in a grid (bias rescaled by sqrt(s)).
If a Poisson catalogue passes, the gate cannot validate the lognormal recipe's two-point function.
No RNG. Command (from worktree root):
  /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python audit/reverse_zero/skeptic_statistics/check_e5b_gate_power.py
"""
import json, os, re
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
E5 = os.path.join(HERE, "..", "E5-cosmic-web-tda-scaled")
src = open(os.path.join(E5, "e5b_lognormal_and_poisson_null.py")).read()
thr = float(re.search(r"^GATE_Z_THRESHOLD\s*=\s*([0-9.]+)", src, re.M).group(1))
g = json.load(open(os.path.join(E5, "e5b_lognormal_and_poisson_null_report.json")))["lognormal_gate"]
xd = np.array(g["xi_data_fit_range"]); xm = np.array(g["xi_mock_mean_fit_range"])
sm = np.array(g["sigma_mock_field_scatter_fit_range"]); sd = np.array(g["sigma_data_shotnoise_fit_range"])
n = g["n_mocks_used_for_gate"]
st = np.sqrt(sm ** 2 * (1 + 1 / n) + sd ** 2)
rms = lambda x: float(np.sqrt(np.mean(((x - xd) / st) ** 2)))
scales = [0, 0.25, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5]
out = {"gate_threshold_from_source": thr, "committed_rms_z": g["rms_z"], "recomputed_rms_z": rms(xm),
       "rms_z_poisson_xi0": rms(np.zeros_like(xd)), "poisson_passes_gate": rms(np.zeros_like(xd)) < thr,
       "rms_z_vs_scale_of_mock_xi": {str(s): rms(s * xm) for s in scales},
       "passing_scale_range": [s for s in scales if rms(s * xm) < thr],
       "median_sigma_total": float(np.median(st)), "median_sigma_mock_field": float(np.median(sm)),
       "median_sigma_shot": float(np.median(sd)), "max_xi_data_fit_range": float(xd.max()),
       "chi2_diag_committed_mock": float(np.sum(((xm - xd) / st) ** 2)), "chi2_diag_poisson": float(np.sum((xd / st) ** 2)),
       "n_bins": int(xd.size),
       "note": "RMS-z ignores bin-to-bin correlation; the counterfactuals use the committed sigma arrays unchanged."}
json.dump(out, open(os.path.join(HERE, "check_e5b_gate_power.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
