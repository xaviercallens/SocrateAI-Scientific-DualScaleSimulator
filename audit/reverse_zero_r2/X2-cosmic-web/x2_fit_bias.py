#!/usr/bin/env python3
"""Bias fit on F=[8,20) Mpc/h only (bins 0..2). cd audit/reverse_zero_r2/X2-cosmic-web && python x2_fit_bias.py
Reads mocks/biasgrid_*.npz (17 b values x 10 mocks), the data xi from x2_lib.Ctx. Writes bias_fit.json and data_xi.npz."""
import sys, json, pathlib
import numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import x2_lib as L
c = L.Ctx()
D = c.cgrid(c.data_xyz)
xi_d = c.xi_from_grid(D)
np.savez(L.HERE / "data_xi.npz", xi=xi_d)
bs = np.round(0.8 + 0.1 * np.arange(17), 6)
X = np.array([[np.load(L.HERE / "mocks" / ("biasgrid_%05d.npz" % (10 * i + j)))["xi"][:3] for j in range(10)] for i in range(17)])
mean = X.mean(1); var = X.var(1, ddof=1)
grid = np.arange(0.8, 2.4001, 0.001)
chi = np.zeros_like(grid); mods = []
for k in range(3):
    pm = np.polyfit(bs, mean[:, k], 2); pv = np.polyfit(bs, np.log(var[:, k]), 2)
    m = np.polyval(pm, grid); v = np.exp(np.polyval(pv, grid))
    chi += (xi_d[k] - m) ** 2 / v; mods.append(pm.tolist())
i = int(np.argmin(chi)); bfit = float(grid[i])
# 1-sigma range from delta chi2 = 1
inside = grid[chi <= chi[i] + 1.0]
res = dict(b_fit=bfit, chi2_min_F=float(chi[i]), dof=3 - 1, b_range_dchi2_1=[float(inside.min()), float(inside.max())],
           at_grid_edge=bool(i in (0, len(grid) - 1)), data_xi_F=xi_d[:3].tolist(), mock_mean_by_b=mean.tolist(),
           mock_var_by_b=var.tolist(), b_values=bs.tolist(), poly_deg2_mean_coeffs=mods,
           method="per-bin quadratic-in-b fit to 10-mock means; variance = exp(quadratic fit to log per-b sample variance); diagonal chi2 on F only",
           seeds="6000000+100*i+j", command="cd audit/reverse_zero_r2/X2-cosmic-web && python x2_fit_bias.py")
json.dump(res, open(L.HERE / "bias_fit.json", "w"), indent=1)
print(json.dumps({k: res[k] for k in ("b_fit", "chi2_min_F", "b_range_dchi2_1", "at_grid_edge", "data_xi_F")}))
