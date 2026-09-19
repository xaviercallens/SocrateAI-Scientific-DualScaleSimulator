"""Test 2 analysis: TDA of rotating-BEC ground states from gpe_rotating.py.

For each Omega:
  R_eff^2 = 2 sqrt(g/pi)/omega', omega' = sqrt(1-Omega^2); n0 = omega' sqrt(g/pi)/g.
  E2a  lower-star path: cmb_tda.betti_curves_from_topology on the open grid
       restricted to r < 0.7 R_eff (qf_common.masked_grid_topology), field
       |psi|^2, sublevel, threshold 0.1 n0 (converted to the function's sigma
       units). b0 there vs N_winding = number of grid plaquettes (all 4
       corners in the disk) with nonzero phase winding of psi.
  E2b  Feynman: winding-vortex density inside r < 0.5 R_eff and r < 0.7 R_eff
       vs Omega/pi.
  E2c  alpha path on winding-vortex positions (plaquette centres) inside
       r < 0.7 R_eff: a_TDA = sqrt(3) median(H1 death) vs a_F =
       sqrt(2 pi/(sqrt3 Omega)); IQR/median of finite H0 deaths; top_bars.
       Cross-check: mean local |psi6| (scipy Delaunay) over vortices not on
       the convex hull.
  E2d  controls: Omega = 0 (no vortex); uniform random points, same N, same
       disk, numpy default_rng(101..110).
Command: prlimit --as=8589934592 -- .venv-tda/bin/python gpe_tda.py
Output: ../results/gpe_tda.json
"""
import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qf_common as qc  # noqa: E402


def interior_psi6(pts):
    from scipy.spatial import ConvexHull
    p6, _ = qc.psi6_local(pts)
    hull = set(ConvexHull(pts).vertices.tolist())
    inner = np.array([j for j in range(len(pts)) if j not in hull])
    return float(np.mean(np.abs(p6[inner]))) if inner.size else None, int(inner.size)


def alpha_stats(pts, aF, web):
    info, by_dim = web.alpha_persistence(pts, (3 * aF) ** 2, "gpe", None)
    d0 = qc.h0_finite_deaths(by_dim)
    h1 = np.array(by_dim[1]) if by_dim[1] else np.empty((0, 2))
    return {"h0": qc.spread_stats(d0),
            "h1_n": int(len(h1)),
            "h1_death_median": float(np.median(h1[:, 1])) if len(h1) else None,
            "a_TDA": float(np.sqrt(3) * np.median(h1[:, 1])) if len(h1) else None,
            "top3_h0_bars": web.top_bars(by_dim, 0, k=3, r_trunc=3 * aF).tolist(),
            "betti_at_truncation": info["betti_numbers_at_truncation"]}


def main():
    cmb, web = qc.cmb(), qc.web()
    out = {}
    for fn in sorted(glob.glob(os.path.join(qc.DATA_ROOT, "gpe", "psi_Omega*.npz"))):
        d = np.load(fn)
        psi = d["psi"].astype(np.complex128); x = d["x"]; g = float(d["g"]); Om = float(d["Omega"])
        dx = x[1] - x[0]
        X, Y = np.meshgrid(x, x, indexing="xy")
        wp = np.sqrt(1 - Om ** 2)
        R = np.sqrt(2 * np.sqrt(g / np.pi) / wp)
        n0 = wp * np.sqrt(g / np.pi) / g
        rr = np.hypot(X, Y)
        dens = np.abs(psi) ** 2
        r = {"Omega": Om, "g": g, "dx": float(dx), "R_eff": float(R), "n0_TF": float(n0),
             "converged": bool(d["converged"]), "steps": int(d["steps"]), "mu": float(d["mu"]),
             "norm_check": float(dens.sum() * dx * dx),
             "max_density": float(dens.max())}
        # winding on all plaquettes
        ph = np.angle(psi)
        q = qc.plaquette_winding(ph)  # periodic roll; plaquettes touching the box edge excluded below
        inside = rr < 0.7 * R
        in4 = inside & np.roll(inside, -1, 1) & np.roll(inside, -1, 0) & np.roll(np.roll(inside, -1, 0), -1, 1)
        in4[-1, :] = False; in4[:, -1] = False
        qv = np.where(in4, q, 0)
        jj, ii = np.nonzero(qv)
        pos = np.stack([x[ii] + dx / 2, x[jj] + dx / 2], 1)
        r["N_winding_r<0.7R"] = int(len(pos)); r["winding_signs"] = {"plus": int((qv > 0).sum()), "minus": int((qv < 0).sum())}
        # E2a lower-star path
        u, e, t = qc.masked_grid_topology(inside)
        f = dens.ravel()
        sig = f[u].std()
        b0, b1, _ = cmb.betti_curves_from_topology(f, u, e, t, np.array([0.1 * n0 / sig, 0.05 * n0 / sig, 0.2 * n0 / sig]), sublevel=True)
        r["lower_star_b0_at_0.1n0"] = int(b0[0]); r["lower_star_b0_at_0.05n0_and_0.2n0"] = [int(b0[1]), int(b0[2])]
        r["lower_star_b1_at_0.1n0"] = int(b1[0])
        r["min_density_in_disk_over_n0"] = float(dens[inside].min() / n0)
        r["E2a_diff_b0_minus_Nwinding"] = int(b0[0]) - int(len(pos))
        if len(pos):
            rv = np.hypot(pos[:, 0], pos[:, 1])
            r["vortex_radii_nearest_disk_edge"] = sorted((0.7 * R - rv).tolist())[:3]
        if Om > 0:
            aF = np.sqrt(2 * np.pi / (np.sqrt(3) * Om))
            nF = Om / np.pi
            rv = np.hypot(pos[:, 0], pos[:, 1])
            n05 = int((rv < 0.5 * R).sum())
            r["feynman"] = {"n_F": nF, "a_F": float(aF),
                            "N_r<0.5R": n05, "density_r<0.5R_over_nF": float(n05 / (np.pi * (0.5 * R) ** 2) / nF),
                            "density_r<0.7R_over_nF": float(len(pos) / (np.pi * (0.7 * R) ** 2) / nF)}
            st = alpha_stats(pos, aF, web)
            st["a_TDA_over_a_F"] = st["a_TDA"] / aF if st["a_TDA"] else None
            st["2median_h0_over_a_F"] = 2 * st["h0"]["median"] / aF
            st["psi6_interior_mean_abs"], st["n_interior"] = interior_psi6(pos)
            r["alpha"] = st
            ctrl = []
            for s in range(101, 111):
                rng = np.random.default_rng(s)
                rad = 0.7 * R * np.sqrt(rng.random(len(pos))); ang = 2 * np.pi * rng.random(len(pos))
                rp = np.stack([rad * np.cos(ang), rad * np.sin(ang)], 1)
                cs = alpha_stats(rp, aF, web)
                p6, ni = interior_psi6(rp)
                ctrl.append({"seed": s, "h0_iqr_over_median": cs["h0"]["iqr_over_median"], "psi6_interior": p6,
                             "a_TDA_over_a_F": cs["a_TDA"] / aF})
            r["random_control"] = {"per_seed": ctrl,
                                   "h0_iqr_over_median_mean": float(np.mean([c["h0_iqr_over_median"] for c in ctrl])),
                                   "psi6_interior_mean": float(np.mean([c["psi6_interior"] for c in ctrl]))}
        np.save(os.path.join(qc.DATA_ROOT, "gpe", f"vortex_positions_Omega{Om:.2f}.npy"), pos)
        out[f"Omega{Om:.2f}"] = r
        print(json.dumps({k: v for k, v in r.items() if k not in ("random_control",)}, default=float)[:1500], flush=True)
    with open(os.path.join(qc.RESULTS, "gpe_tda.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=float)


if __name__ == "__main__":
    main()
