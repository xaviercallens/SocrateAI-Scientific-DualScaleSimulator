"""POST-HOC (not pre-registered) descriptive check for test 3, E3d: the
pre-registered psi6 estimator (global |<exp(6 i theta_b)>| over all Delaunay
bonds, psi6_global_bonds) gave 0.229 at 3 kOe vs the published ~0.035. The
paper's exact normalisation is not recoverable from the text extraction, so
this script reports alternative readings on the SAME minima (saved by
stm_vortex_tda.py), without changing the pre-registered verdict:
  (a) global bond order, bonds longer than 1.6 a_tri removed (hull artefacts)
  (b) per-site mean |psi6_j| over non-hull sites (qf_common.psi6_local)
  (c) |mean_j psi6_j| over non-hull sites (global order of site averages)
Command: prlimit --as=8589934592 -- .venv-tda/bin/python stm_psi6_variants.py
Output: ../results/stm_psi6_variants.json
"""
import glob
import json
import os
import sys

import numpy as np
from scipy.spatial import ConvexHull, Delaunay

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qf_common as qc  # noqa: E402

PHI0 = 2.067833848e-15


def variants(pts, a):
    tri = Delaunay(pts)
    indptr, nb = tri.vertex_neighbor_vertices
    angs = []
    for j in range(len(pts)):
        for k in nb[indptr[j]:indptr[j + 1]]:
            if k > j:
                dv = pts[k] - pts[j]
                if np.hypot(*dv) <= 1.6 * a:
                    angs.append(np.arctan2(dv[1], dv[0]))
    ga = float(np.abs(np.mean(np.exp(6j * np.array(angs)))))
    p6, _ = qc.psi6_local(pts)
    hull = set(ConvexHull(pts).vertices.tolist())
    inner = np.array([j for j in range(len(pts)) if j not in hull])
    return ga, float(np.mean(np.abs(p6[inner]))), float(np.abs(np.mean(p6[inner])))


def main():
    out = {}
    for fn in sorted(glob.glob(os.path.join(qc.DATA_ROOT, "stm", "minima_*kOe_460mK.npz"))):
        H = float(os.path.basename(fn).split("_")[1][:-3])
        a = 1.075 * np.sqrt(PHI0 / (H / 10)) * 1e9
        d = np.load(fn)
        v = np.array([variants(d[k], a) for k in sorted(d.files)])
        out[f"{H:g}kOe"] = {"global_bonds_le_1.6a": float(v[:, 0].mean()), "per_site_mean_abs": float(v[:, 1].mean()),
                            "abs_mean_of_site_psi6": float(v[:, 2].mean()), "n_images": int(len(v))}
    out = dict(sorted(out.items(), key=lambda kv: float(kv[0][:-3])))
    with open(os.path.join(qc.RESULTS, "stm_psi6_variants.json"), "w") as fh:
        json.dump({"note": "POST-HOC descriptive; does not change the pre-registered E3d verdict", "by_field": out}, fh, indent=1)
    for k, x in out.items():
        print(k, {kk: round(vv, 3) for kk, vv in x.items()})


if __name__ == "__main__":
    main()
