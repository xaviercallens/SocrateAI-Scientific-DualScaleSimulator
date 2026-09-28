"""STEP 0 diagnostics, as a committed script (CLAUDE.md rule 6: no number without a script).

Four checks whose numbers are quoted in results/step0_verdict.json:
  1. a_fft estimator bias on SYNTHETIC conductance maps with the same geometry
     as the real 20 kOe scans (known answer: the input a).
  2. The six resolved Bragg spots of a real 20 kOe map, each giving a.
  3. The pooled g(r) first peak and the Delaunay bond-length mode over all
     20 images at 20 kOe.
  4. Border over-detection: nearest-neighbour distance distribution, points on
     a border line, and the interior-versus-full areal density.

Run: prlimit --as=8589934592 -- .venv-tda/bin/python step0_diagnostics.py
"""
import glob
import os

import numpy as np
from scipy import ndimage, spatial

import qf_lib as q

SEED = 1
H_KOE = 20.0


def synth_img(a_true, sf, core_nm, noise, n, px, L, rng):
    pts, dy = [], a_true * np.sqrt(3) / 2
    for j in range(int(L / dy) + 3):
        for i in range(int(L / a_true) + 3):
            pts.append([i * a_true + (j % 2) * a_true / 2 - a_true, j * dy - a_true])
    pts = np.array(pts) + rng.normal(0, sf * a_true, (len(pts), 2))
    yy, xx = np.mgrid[0:n, 0:n]
    X, Y = xx * px, yy * px
    im = np.zeros((n, n))
    for p in pts:
        im -= np.exp(-((X - p[0]) ** 2 + (Y - p[1]) ** 2) / (2 * core_nm ** 2))
    return im + rng.normal(0, noise, im.shape)


def main():
    assert q.assert_alpha_convention()
    a = q.a_tri(H_KOE / 10.0)
    fns = sorted(glob.glob(os.path.join(q.STM_ROOT, f"{int(H_KOE)} kOe", "*.sxm")))
    imgs, meta = q.read_sxm(fns[1])
    img = imgs[("Input_7", "fwd")]
    px = meta["range_x_m"] * 1e9 / meta["nx"]
    L = meta["range_x_m"] * 1e9
    n = meta["nx"]
    out = {"a_target_nm": a, "scan_nm": L, "px_nm": px, "nx": n,
           "representative_file": os.path.basename(fns[1]), "seed": SEED,
           "command": q.command()}

    # --- 1. a_fft bias on synthetic maps (known answer)
    rng = np.random.default_rng(SEED)
    bias = []
    for sf in [0.0, 0.05, 0.10, 0.15]:
        for cd in [6.0, 4.0]:
            si = synth_img(a, sf, a / cd, 0.05, n, px, L, rng)
            af, _, _, _ = q.a_fft(si, px, a)
            bias.append({"sigma_over_a": sf, "core_a_over": cd, "a_fft_nm": af,
                         "rel_error": af / a - 1})
    si36 = synth_img(36.0, 0.05, a / 5, 0.05, n, px, L, rng)
    af36, _, _, _ = q.a_fft(si36, px, a)
    out["synthetic_a_fft_bias"] = bias
    out["synthetic_a_true_36nm_recovered"] = af36

    # --- 2. six Bragg spots
    im = q.plane_subtract(img)
    w = np.outer(np.hanning(n), np.hanning(n))
    im = (im - im.mean()) * w
    pad = 8
    F = np.fft.fftshift(np.abs(np.fft.fft2(im, s=(n * pad, n * pad))) ** 2)
    qx = 2 * np.pi * np.fft.fftshift(np.fft.fftfreq(n * pad, d=px))
    QX, QY = np.meshgrid(qx, qx)
    Q = np.sqrt(QX ** 2 + QY ** 2)
    q1e = 4 * np.pi / (np.sqrt(3) * a)
    m = (Q > 0.6 * q1e) & (Q < 1.5 * q1e)
    Fm = np.where(m, F, 0)
    mx = ndimage.maximum_filter(Fm, size=60)
    py, pxx = np.nonzero((Fm == mx) & (Fm > 0.15 * Fm.max()))
    vals = Fm[py, pxx]
    order = np.argsort(vals)[::-1][:6]
    spots = [{"q_inv_nm": float(Q[py[i], pxx[i]]),
              "a_nm": float(4 * np.pi / (np.sqrt(3) * Q[py[i], pxx[i]])),
              "angle_deg": float(np.degrees(np.arctan2(QY[py[i], pxx[i]], QX[py[i], pxx[i]])))}
             for i in order]
    out["bragg_spots"] = spots
    out["bragg_spots_a_median_nm"] = float(np.median([s["a_nm"] for s in spots]))

    # --- 3. g(r) and Delaunay bond mode over all 20 images
    allpts, bonds, dists = [], [], []
    for f in fns:
        ii, _ = q.read_sxm(f)
        p = q.detect_minima(ii[("Input_7", "fwd")], px, a)
        allpts.append(p)
        tri = spatial.Delaunay(p)
        E = set()
        for s in tri.simplices:
            for k in range(3):
                E.add(tuple(sorted((int(s[k]), int(s[(k + 1) % 3])))))
        E = np.array(sorted(E))
        bonds.append(np.linalg.norm(p[E[:, 1]] - p[E[:, 0]], axis=1))
        dd = spatial.distance.pdist(p)
        dists.append(dd[dd < 2.5 * a])
    b = np.concatenate(bonds); b = b[b < 2 * a]
    h, e = np.histogram(b, bins=60, range=(0, 2 * a))
    c = (e[:-1] + e[1:]) / 2
    dall = np.concatenate(dists)
    hg, eg = np.histogram(dall, bins=100, range=(0, 2.5 * a))
    cg = (eg[:-1] + eg[1:]) / 2
    gr = hg / (2 * np.pi * cg * (eg[1] - eg[0]))
    sel = cg > 0.5 * a
    out["delaunay_bond"] = {"median_nm": float(np.median(b)), "mode_nm": float(c[np.argmax(h)]),
                            "mean_nm": float(b.mean()), "n_bonds": int(b.size)}
    out["g_of_r_first_peak_nm"] = float(cg[sel][np.argmax(gr[sel])])

    # --- 4. border over-detection
    nn = np.concatenate([q.nn_distance(p) for p in allpts])
    onb, nin, ntot = [], [], []
    mgn = a / 2
    for p in allpts:
        onb.append(int(np.sum((p[:, 0] == 0) | (p[:, 1] == 0) |
                              (p[:, 0] > L - 2) | (p[:, 1] > L - 2))))
        ii = ((p[:, 0] >= mgn) & (p[:, 0] <= L - mgn) & (p[:, 1] >= mgn) & (p[:, 1] <= L - mgn))
        nin.append(int(ii.sum())); ntot.append(len(p))
    A_in = (L - 2 * mgn) ** 2
    out["border"] = {
        "nn_min_nm": float(nn.min()),
        "nn_percentiles_nm": {str(p): float(v) for p, v in zip([1, 5, 10, 50], np.percentile(nn, [1, 5, 10, 50]))},
        "frac_nn_below_half_a": float((nn < 0.5 * a).mean()),
        "min_pairwise_distance_img00_nm": float(spatial.distance.pdist(allpts[0]).min()),
        "points_exactly_on_a_border_line_median": float(np.median(onb)),
        "n_full_median": float(np.median(ntot)), "n_interior_median": float(np.median(nin)),
        "a_density_full_nm": q.a_density(float(np.median(ntot)), L * L),
        "a_density_interior_nm": q.a_density(float(np.median(nin)), A_in),
    }
    for k in ["bragg_spots_a_median_nm", "g_of_r_first_peak_nm", "synthetic_a_true_36nm_recovered"]:
        print(k, out[k])
    print("delaunay", out["delaunay_bond"])
    print("border", {k: v for k, v in out["border"].items() if not isinstance(v, dict)})
    print("wrote", q.write_json("step0_diagnostics.json", out))


if __name__ == "__main__":
    main()
