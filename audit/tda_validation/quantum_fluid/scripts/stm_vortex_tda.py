"""Test 3: real STS conductance maps of the vortex state in a 20 nm a-Re6Zr film
(Duhan et al., Nat. Commun. 16, 2100 (2025); data Zenodo 10.5281/zenodo.14780459,
'STM data.zip', md5 5b97a8b60660f1f7b521f0f7a526502f, CC-BY-4.0).

.sxm (Nanonis) parser: ASCII header up to ':SCANIT_END:', then the 0x1A 0x04
marker, then big-endian float32 images, channels listed in DATA_INFO, each
'both' channel stored forward then backward, SCAN_PIXELS (nx, ny).
Channel used: 'Input_7' (lock-in conductance, V), forward.
Vortex detector (non-TDA, fixed before looking at any image, see
expectations.json): plane subtraction, Gaussian low-pass in Fourier space with
real-space sigma = a_tri(B)/6, local minima = pixels equal to the minimum of a
disk of radius a_tri(B)/3 (scipy.ndimage.minimum_filter, mode='nearest').
a_tri(B) = 1.075 sqrt(Phi0/B), B = mu0 H (H in kOe -> B = H/10 T).
TDA (pipeline under test): cosmic_web_tda_scaled.alpha_persistence on the
minima (nm), max_alpha_sq = (3 a_tri)^2; top_bars with r_trunc = 3 a_tri.
Cross-check (non-TDA): global |<exp(6 i theta_b)>| over scipy Delaunay bonds.
Negative control: uniform random points, same count, same area,
numpy default_rng(201 + image index).

Command: prlimit --as=8589934592 -- .venv-tda/bin/python stm_vortex_tda.py
Output: ../results/stm_vortex_tda.json and DATA_ROOT/stm/minima_*.npz
"""
import glob
import json
import os
import re
import sys

import numpy as np
from scipy import ndimage

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qf_common as qc  # noqa: E402

PHI0 = 2.067833848e-15  # Wb (paper quotes 2.068e-15)
ROOT = os.path.join(qc.DATA_ROOT, "raw/zenodo_14780459/extracted/STM data/Fig 1 and 2")


def read_sxm(path):
    b = open(path, "rb").read()
    i = b.find(b":SCANIT_END:")
    hdr = b[:i].decode("latin1")
    j = b.find(b"\x1a\x04", i)
    data = b[j + 2:]
    tags = {}
    for m in re.finditer(r":([^:\n]+):\n(.*?)(?=\n:[^:\n]+:\n|\Z)", hdr, re.S):
        tags[m.group(1)] = m.group(2)
    nx, ny = [int(v) for v in tags["SCAN_PIXELS"].split()]
    rx, ry = [float(v) for v in tags["SCAN_RANGE"].split()]
    lines = [ln.strip().split("\t") for ln in tags["DATA_INFO"].strip().split("\n")]
    head, rows = lines[0], lines[1:]
    chans = []
    for r in rows:
        d = dict(zip(head, r))
        dirs = ["fwd", "bwd"] if d["Direction"] == "both" else ["fwd"]
        for dd in dirs:
            chans.append((d["Name"], dd))
    arr = np.frombuffer(data[: len(chans) * nx * ny * 4], dtype=">f4").reshape(len(chans), ny, nx)
    imgs = {c: arr[k].astype(np.float64) for k, c in enumerate(chans)}
    return imgs, dict(nx=nx, ny=ny, range_x_m=rx, range_y_m=ry, scan_dir=tags.get("SCAN_DIR", "").strip(),
                      bias=tags.get("BIAS", "").strip(), file=tags.get("SCAN_FILE", "").strip())


def detect_minima(img, px_nm, a_nm):
    ny, nx = img.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    ok = np.isfinite(img)
    A = np.c_[xx[ok], yy[ok], np.ones(ok.sum())]
    coef, *_ = np.linalg.lstsq(A, img[ok], rcond=None)
    im = img - (coef[0] * xx + coef[1] * yy + coef[2])
    im[~ok] = np.nanmean(im[ok])
    s_px = (a_nm / 6.0) / px_nm
    kx = 2 * np.pi * np.fft.fftfreq(nx); ky = 2 * np.pi * np.fft.fftfreq(ny)
    KX, KY = np.meshgrid(kx, ky)
    f = np.real(np.fft.ifft2(np.fft.fft2(im) * np.exp(-0.5 * s_px ** 2 * (KX ** 2 + KY ** 2))))
    rad = (a_nm / 3.0) / px_nm
    R = int(np.ceil(rad))
    fy, fx = np.mgrid[-R:R + 1, -R:R + 1]
    fp = (fx ** 2 + fy ** 2) <= rad ** 2
    mn = ndimage.minimum_filter(f, footprint=fp, mode="nearest")
    my, mx = np.nonzero(f == mn)
    pts = np.stack([mx * px_nm, my * px_nm], 1)
    return pts


def main():
    web = qc.web()
    res = {"fields": {}}
    folders = sorted(glob.glob(os.path.join(ROOT, "* kOe")), key=lambda p: float(os.path.basename(p).split()[0]))
    os.makedirs(os.path.join(qc.DATA_ROOT, "stm"), exist_ok=True)
    for fo in folders:
        H = float(os.path.basename(fo).split()[0])
        B = H / 10.0
        a_nm = 1.075 * np.sqrt(PHI0 / B) * 1e9
        files = sorted(glob.glob(os.path.join(fo, "*.sxm")))
        d0_all, h1d_all, ps6, counts, exp_counts = [], [], [], [], []
        d0_rand, ps6_rand = [], []
        per_img = []
        allpts = {}
        for k, fn in enumerate(files):
            imgs, meta = read_sxm(fn)
            img = imgs[("Input_7", "fwd")]
            px_nm = meta["range_x_m"] * 1e9 / meta["nx"]
            assert abs(meta["range_x_m"] / meta["nx"] - meta["range_y_m"] / meta["ny"]) < 1e-12
            pts = detect_minima(img, px_nm, a_nm)
            area_m2 = meta["range_x_m"] * meta["range_y_m"]
            n_exp = B * area_m2 / PHI0
            info, by_dim = web.alpha_persistence(pts, (3 * a_nm) ** 2, f"H{H}_{k}", None)
            d0 = qc.h0_finite_deaths(by_dim)
            h1 = np.array(by_dim[1]) if by_dim[1] else np.empty((0, 2))
            p6, nb = qc.psi6_global_bonds(pts) if len(pts) >= 3 else (np.nan, 0)
            rng = np.random.default_rng(201 + k)
            rp = rng.random(pts.shape) * np.array([meta["range_x_m"], meta["range_y_m"]]) * 1e9
            _, bdr = web.alpha_persistence(rp, (3 * a_nm) ** 2, "rand", None)
            d0_rand.append(qc.h0_finite_deaths(bdr)); ps6_rand.append(qc.psi6_global_bonds(rp)[0])
            d0_all.append(d0); h1d_all.append(h1[:, 1] if len(h1) else np.empty(0))
            ps6.append(p6); counts.append(len(pts)); exp_counts.append(n_exp)
            allpts[f"img{k:02d}"] = pts
            per_img.append(dict(file=os.path.basename(fn), n_minima=len(pts), n_expected=n_exp, psi6_global=p6,
                                h0_median=float(np.median(d0)) if d0.size else None))
        np.savez_compressed(os.path.join(qc.DATA_ROOT, "stm", f"minima_{int(H)}kOe_460mK.npz"), **allpts)
        d0c = np.concatenate(d0_all); h1c = np.concatenate(h1d_all); d0r = np.concatenate(d0_rand)
        a_tda = float(np.sqrt(3) * np.median(h1c)) if h1c.size else None
        res["fields"][f"{H:g}kOe"] = {
            "H_kOe": H, "B_T": B, "a_tri_nm": float(a_nm), "n_images": len(files),
            "pixel_nm": px_nm, "scan_nm": meta["range_x_m"] * 1e9, "channel": "Input_7 fwd",
            "count_median": float(np.median(counts)), "expected_count_B_A_over_Phi0": float(exp_counts[0]),
            "count_ratio_median": float(np.median(np.array(counts) / np.array(exp_counts))),
            "h0_pooled": qc.spread_stats(d0c),
            "a_TDA_sqrt3_median_h1_death_nm": a_tda,
            "a_TDA_over_a_tri": a_tda / a_nm if a_tda else None,
            "a_from_2_median_h0_over_a_tri": float(2 * np.median(d0c) / a_nm),
            "psi6_global_mean": float(np.nanmean(ps6)), "psi6_global_std": float(np.nanstd(ps6)),
            "random_control": {"h0_pooled": qc.spread_stats(d0r), "psi6_global_mean": float(np.mean(ps6_rand))},
            "per_image": per_img,
        }
        r = res["fields"][f"{H:g}kOe"]
        print(f"H={H:g} kOe a={a_nm:.1f}nm n={r['count_median']:.0f}/{r['expected_count_B_A_over_Phi0']:.0f} "
              f"aTDA/a={r['a_TDA_over_a_tri']:.3f} H0 IQR/med={r['h0_pooled']['iqr_over_median']:.3f} "
              f"psi6={r['psi6_global_mean']:.3f} | rand IQR/med={r['random_control']['h0_pooled']['iqr_over_median']:.3f} "
              f"psi6={r['random_control']['psi6_global_mean']:.3f}", flush=True)
    res["meta"] = {"Phi0_Wb": PHI0, "detector": "plane-subtract, Gaussian LP sigma=a/6, disk-min radius a/3",
                   "data_root": ROOT}
    with open(os.path.join(qc.RESULTS, "stm_vortex_tda.json"), "w") as fh:
        json.dump(res, fh, indent=1)


if __name__ == "__main__":
    main()
