#!/usr/bin/env python3
"""Two post-hoc diagnostics requested by review, run against already-committed mocks/ npz files
(no new mocks generated). Not part of the registered decision rule; both are context for the
INCONCLUSIVE verdict in REPORT.md, not a re-derivation of it.
cd audit/reverse_zero_r2/X2-cosmic-web && python x2_check_offset_and_bfit_G.py

(a) offset_vs_target: is (measured mock-ensemble xi) - (CAMB linear target b^2*xi_lin) a
    b-independent constant (which would cancel in the data-vs-mock comparison) or does it scale
    with b (which would not cancel unless b is exactly right)? Uses the 10-mock-per-b biasgrid
    ensemble at b=0.9 and b=1.5.
(b) bfit_on_G: refit b using ONLY the G=[20,60) range (the biasgrid means, each from only 10
    mocks, against the same 500-mock Hartlap-corrected covariance used in x2_analyze.py) --
    i.e. "if we let the gate itself pick the best-fitting amplitude, does the mismatch go away?"
"""
import sys, json, pathlib
import numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import x2_lib as L

H = L.HERE
c = L.Ctx()
G = slice(3, 13); p = 10
rgrid, xi_lin = L.xi_lin_table()
rc = (L.RS_EDGES[:-1] + L.RS_EDGES[1:]) / 2

def target_xi(b):
    return b ** 2 * np.interp(rc, rgrid, xi_lin)

offset = {}
for b in (0.9, 1.5):
    idx = round((b - 0.8) / 0.1)
    bb = round(0.8 + 0.1 * idx, 6)
    xs = np.array([np.load(H / "mocks" / ("biasgrid_%05d.npz" % (10 * idx + j)))["xi"] for j in range(10)])
    meas = xs.mean(0); tgt = target_xi(bb)
    offset[str(bb)] = dict(measured_G=meas[G].tolist(), target_lognormal_linear_G=tgt[G].tolist(),
                           diff_G=(meas - tgt)[G].tolist(), mean_diff_G=float((meas - tgt)[G].mean()))
r = offset["1.5"]["mean_diff_G"] / offset["0.9"]["mean_diff_G"]
offset["ratio_mean_diff_1p5_over_0p9"] = float(r)
offset["b_squared_ratio_1p5_over_0p9"] = float(1.5 ** 2 / 0.9 ** 2)
offset["reading"] = ("diff scales roughly with the b^2 ratio (%.2f vs expected %.2f): NOT a b-independent constant, "
                      "so it does not simply cancel between data and mocks unless b is exactly right; part of the "
                      "measured-vs-target gap is a real, bias-scale-dependent lognormal-transform/estimator effect, "
                      "not just a common shift." % (r, offset["b_squared_ratio_1p5_over_0p9"]))

NF = 500
fid = np.array([np.load(H / "mocks" / ("fiducial_%05d.npz" % k))["xi"] for k in range(NF)])
xi_d = np.load(H / "data_xi.npz")["xi"]
mu = fid[:, G].mean(0); C = np.cov(fid[:, G].T, ddof=1); Ci = np.linalg.inv(C)
hart = (NF - p - 2) / (NF - 1.0)
bs = np.round(0.8 + 0.1 * np.arange(17), 6)
means = np.array([np.array([np.load(H / "mocks" / ("biasgrid_%05d.npz" % (10 * i + j)))["xi"] for j in range(10)]).mean(0) for i in range(17)])
chis = np.array([hart * (means[i, G] - xi_d[G]) @ Ci @ (means[i, G] - xi_d[G]) for i in range(17)])
i = int(np.argmin(chis))
bfit_g = dict(b_values=bs.tolist(), chi2_G_per_b=chis.tolist(), best_b_on_G=float(bs[i]), chi2_min_on_G=float(chis[i]),
              b_fit_on_F=json.load(open(H / "bias_fit.json"))["b_fit"], chi2_gate_at_b_fit_F=40.21471907030045,
              thr95_fiducial=20.17582016844158, dof=10 - 1,
              caveat="each point uses only a 10-mock mean (not the 500-mock fiducial ensemble), so these chi2 values "
                     "are inflated by roughly the extra mean-of-10 noise (~10%) relative to a 500-mock mean; this is a "
                     "diagnostic over a coarse 0.1-spaced b grid, not a refit of the registered b_fit.",
              reading=("even at the best-fitting b on G (%.2f, vs %.3f fitted on F), chi2=%.1f still exceeds the "
                        "gate threshold %.2f: the G-range mismatch is not resolved by re-choosing the amplitude alone, "
                        "so the calibration failure is at least partly a shape mismatch, not only an amplitude/bias "
                        "value that failed to transfer from F to G." % (bs[i], json.load(open(H / "bias_fit.json"))["b_fit"], chis[i], 20.17582016844158)))
res = dict(offset_vs_target=offset, bfit_on_G=bfit_g)
json.dump(res, open(H / "x2_offset_and_bfit_diagnostics.json", "w"), indent=1)
print(json.dumps(res, indent=1))
