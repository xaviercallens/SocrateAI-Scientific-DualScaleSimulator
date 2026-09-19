#!/usr/bin/env python3
"""Cross-check the gridded FFT xi(r) estimator (x2_lib.Ctx.xi_from_grid) against a direct
Landy-Szalay pair count with cKDTree, on the real data and its registered randoms (5xN_D).
cd audit/reverse_zero_r2/X2-cosmic-web && python x2_paircount_check.py
Full pair counting at N~193k/~968k randoms is too slow for a check script, so this compares the
two estimators on the SAME random sub-sample of n=6000 data points (and alpha-matched randoms
drawn from the shared randoms grid), which the two estimators are free to disagree on due to
NGP/grid binning -- agreement within ~10% validates the grid method's basic normalisation."""
import sys, json, pathlib
import numpy as np
from scipy.spatial import cKDTree
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import x2_lib as L

c = L.Ctx()
rng = np.random.default_rng(1)
n = 6000
di = rng.choice(c.N, n, replace=False)
D = c.data_xyz[di]
NRsub = n * 20
rand_full = None
# regenerate a fresh random subsample directly (cheap) rather than reusing the 5x grid (keeps this script standalone)
rxyz = c.draw_shell_points(np.bincount(np.clip(np.searchsorted(c.edges, np.sqrt((D ** 2).sum(1)), "right") - 1, 0, c.nsh - 1), minlength=c.nsh),
                            np.random.default_rng(2), mult=20)
alpha = n / len(rxyz)
tD, tR = cKDTree(D), cKDTree(rxyz)
def paircount(t1, t2, edges):
    return np.diff(t1.count_neighbors(t2, edges)).astype(float)
# count_neighbors(t,t,edges) counts each unordered pair twice (i,j) and (j,i); edges[0]=8>0 excludes i==j.
DD = paircount(tD, tD, L.RS_EDGES) / 2.0
DR = paircount(tD, tR, L.RS_EDGES)
RR = paircount(tR, tR, L.RS_EDGES) / 2.0
NR = len(rxyz)
ls = (DD / (n * (n - 1) / 2) - 2 * DR / (n * NR) + RR / (NR * (NR - 1) / 2)) / (RR / (NR * (NR - 1) / 2))
# grid estimator on the same n=6000 subsample with matching alpha
Dsub_grid = c.cgrid(D)
Rsub_grid_full = c.Rg  # this is the full 5x randoms grid, alpha mismatched to the n-subsample -- rebuild small
c.Rg_sub = c.cgrid(rxyz)
alpha_grid = n / NR
F = Dsub_grid - alpha_grid * c.Rg_sub
RRg = c.corr_bins(c.Rg_sub)
xi_grid = c.corr_bins(F) / (alpha_grid ** 2 * RRg)
res = dict(n_sub=n, n_rand_sub=NR, r_bin_centers=((L.RS_EDGES[:-1] + L.RS_EDGES[1:]) / 2).tolist(),
           xi_direct_LS=ls.tolist(), xi_grid=xi_grid.tolist(),
           ratio_grid_over_direct=(xi_grid / np.where(ls != 0, ls, np.nan)).tolist(),
           note="Independent 6000-point subsample; both estimators use the SAME random points for their RR/DR terms so "
                "differences are purely NGP-grid vs exact-pair-distance binning + FFT edge effects, not sample variance.")
json.dump(res, open(L.HERE / "paircount_crosscheck.json", "w"), indent=1)
print(json.dumps({k: res[k] for k in ("r_bin_centers", "xi_direct_LS", "xi_grid", "ratio_grid_over_direct")}, indent=1))
