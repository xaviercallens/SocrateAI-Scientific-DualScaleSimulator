"""DATASET 4 (breadth, REAL data): quasi-2-D colloidal dispersions of 2.8 um particles.

Source: Zenodo 10.5281/zenodo.14518422, CC-BY-4.0, "Dataset for 'Local area
distribution of 2D colloidal dispersions and its relation to particle diffusion:
a Voronoi tessellation analysis'". Description: "2.8 micron particles confined in
a quasi 2D geometry at different concentrations. File name (in .zip format)
indicates concentration."

Two concentrations, one TopoDB dataset each:
  area fraction 0.89 -- a dense 2-D colloidal crystal (the ordered end)
  area fraction 0.11 -- a dilute 2-D colloidal gas (the disordered end)

KNOWN ANSWER, scale-free so it needs no pixel calibration: for a triangular
lattice of hard discs of diameter d at area fraction phi,
    a / d = sqrt(pi / (2 sqrt(3) phi))
which is 1.0096 at phi = 0.89 and 2.871 at phi = 0.11. d is measured here as
the contact distance (1st percentile of nearest-neighbour distances), which for
hard discs is the particle diameter. The implied pixel size d_um / d_px is then
reported and compared with the 2.8 um particle diameter.

Run: PYTHONPATH=<.pylibs> prlimit --as=8589934592 -- .venv-tda/bin/python compute_colloid.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, "/mnt/disks/disk-socrateai-local-1/wt-topo-qfluid/.pylibs")
import h5py  # noqa: E402

import qf_lib as q  # noqa: E402

SEED = 20260920
N_NULL = 200
N_FRAMES = 20
D_UM = 2.8
ROOT = os.path.join(q.DATA_ROOT, "colloid_zenodo_14518422")
SETS = {
    "0.89": {"file": "0_89_Full_skipped10/01_90k/stack.mat", "key": "stack",
             "zip": "0_89_Full_skipped10.zip", "phi": 0.89,
             "zip_sha256": "714b073f9d9f6c5d82b957fb4ff44a6b2b4d41f1fae38aca5be8ad66ea182e6c",
             "url": "https://zenodo.org/api/records/14518422/files/0_89_Full_skipped10.zip/content"},
    "0.11": {"file": "0_11/B_16_1_2_8um_20X_0_5NA_25FPS_90k_Dilute_490nm_New/track.mat", "key": "tr",
             "zip": "0_11.zip", "phi": 0.11,
             "zip_sha256": "517c8f410f3626d4b5bb2d7cfd6a5e87ca845d58823c9fc45257e1130c41893f",
             "url": "https://zenodo.org/api/records/14518422/files/0_11.zip/content"},
}


def frames(path, key, n_frames):
    with h5py.File(path, "r") as h:
        d = h[key]
        fr = np.asarray(d[4, :])
        uf = np.unique(fr)
        pick = uf[np.linspace(0, len(uf) - 1, n_frames).astype(int)]
        xy = np.asarray(d[0:2, :]).T
        for f in pick:
            yield float(f), xy[fr == f]


def main():
    assert q.assert_alpha_convention()
    rng = np.random.default_rng(SEED)
    out = {"meta": {"seed": SEED, "n_null": N_NULL, "n_frames": N_FRAMES,
                    "particle_diameter_um": D_UM, "command": q.command(),
                    "script": os.path.abspath(__file__),
                    "source": "Zenodo 10.5281/zenodo.14518422 (CC-BY-4.0)",
                    "units": "pixels for the point clouds; lengths also reported in units of the "
                             "measured contact distance, which is scale-free"},
           "by_phi": {}}
    for phi_key, S in SETS.items():
        path = os.path.join(ROOT, S["file"])
        phi = S["phi"]
        per, clouds = [], {}
        for k, (f, pts) in enumerate(frames(path, S["key"], N_FRAMES)):
            if len(pts) < 10:
                continue
            xmin, xmax = pts[:, 0].min(), pts[:, 0].max()
            ymin, ymax = pts[:, 1].min(), pts[:, 1].max()
            area = (xmax - xmin) * (ymax - ymin)
            a_ref = q.a_density(len(pts), area)
            nn_all = q.nn_distance(pts)
            d_contact = float(np.percentile(nn_all, 1))
            bars = q.alpha_bars(pts, (3 * a_ref) ** 2)
            d0, d1 = q.h0_deaths(bars), q.h1_deaths(bars)
            p6, _ = q.psi6_global(pts)
            nn_in = q.nn_distance(pts, (xmin, xmax, ymin, ymax), margin=a_ref)
            per.append({"frame": f, "n": int(len(pts)), "area_px2": float(area),
                        "a_density_px": a_ref,
                        "a_nn_px": float(np.median(nn_in)) if nn_in.size else None,
                        "a_h1_px": float(np.sqrt(3) * np.median(d1)) if d1.size else None,
                        "a_h0_px": float(2 * np.median(d0)) if d0.size else None,
                        "d_contact_px": d_contact,
                        "psi6": p6,
                        "h0_iqr_over_median": q.spread_stats(d0)["iqr_over_median"],
                        "betti0": q.betti_from_bars(bars, 0.55 * a_ref)[0],
                        "betti1": q.betti_from_bars(bars, 0.55 * a_ref)[1]})
            if k == 0:
                first_bars, first_aref = bars, a_ref
            clouds[f"frame{k:02d}"] = pts.astype(np.float32)

        def med(kk):
            v = np.array([p[kk] for p in per if p[kk] is not None], float)
            return float(np.median(v[np.isfinite(v)]))

        n_med = int(med("n")); area = med("area_px2"); a_ref = med("a_density_px")
        L = float(np.sqrt(area))
        nul = {"psi6": [], "h0_iqr_over_median": []}
        for _ in range(N_NULL):
            rp = rng.random((n_med, 2)) * L
            bn = q.alpha_bars(rp, (3 * a_ref) ** 2)
            nul["psi6"].append(q.psi6_global(rp)[0])
            nul["h0_iqr_over_median"].append(q.spread_stats(q.h0_deaths(bn))["iqr_over_median"])
        obs = {"psi6": med("psi6"), "h0_iqr_over_median": med("h0_iqr_over_median")}
        pv = {"psi6": q.rank_p(obs["psi6"], nul["psi6"], "greater"),
              "h0_iqr_over_median": q.rank_p(obs["h0_iqr_over_median"], nul["h0_iqr_over_median"], "less")}

        d_c = med("d_contact_px")
        a_nn = med("a_nn_px")
        a_over_d_pred = float(np.sqrt(np.pi / (2 * np.sqrt(3) * phi)))
        npz = os.path.join(q.RESULTS, f"colloid_phi{phi_key.replace('.','p')}.npz")
        np.savez_compressed(npz, **clouds)
        rec = {
            "phi_published": phi, "n_frames": len(per), "n_particles_median": n_med,
            "area_px2": area, "d_contact_px": d_c,
            "implied_pixel_size_um": D_UM / d_c,
            "a_nn_px": a_nn, "a_density_px": a_ref,
            "a_h1_px": med("a_h1_px"), "a_h0_px": med("a_h0_px"),
            "a_nn_over_d_contact": a_nn / d_c,
            "a_density_over_d_contact": a_ref / d_c,
            "a_over_d_predicted_from_phi": a_over_d_pred,
            "known_answer_rel_error": (a_ref / d_c) / a_over_d_pred - 1,
            "phi_implied_from_geometry": float(n_med * np.pi * (d_c / 2) ** 2 / area),
            "observed": obs, "p_values_rank": pv, "n_null": N_NULL,
            "null_median": {k: float(np.median(v)) for k, v in nul.items()},
            "betti_median": {0: int(med("betti0")), 1: int(med("betti1"))},
            "top_h1_bars": q.top_bars(first_bars, 1, 10),
            "top_h0_bars": q.top_bars(first_bars, 0, 10, r_trunc=3 * first_aref),
            "clouds_npz": npz, "clouds_sha256": q.sha256_file(npz),
            "zip_sha256": S["zip_sha256"], "url": S["url"],
            "per_frame": per,
        }
        out["by_phi"][phi_key] = rec
        print(f"phi={phi} n={n_med} d_contact={d_c:.2f}px -> px={D_UM/d_c*1000:.1f} nm  "
              f"a_density/d={a_ref/d_c:.4f} predicted {a_over_d_pred:.4f} "
              f"(rel {100*rec['known_answer_rel_error']:+.2f}%)  phi_implied={rec['phi_implied_from_geometry']:.3f}  "
              f"psi6={obs['psi6']:.3f}(null {np.median(nul['psi6']):.3f}, p={pv['psi6']:.4f})  "
              f"iqr/med={obs['h0_iqr_over_median']:.3f}(null {np.median(nul['h0_iqr_over_median']):.3f})",
              flush=True)
    print("wrote", q.write_json("colloid.json", out))


if __name__ == "__main__":
    main()
