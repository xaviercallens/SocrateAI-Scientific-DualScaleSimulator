"""STEP 0 (gate): reproduce the published Re6Zr vortex-lattice constant at 20 kOe.

Target (pre-stated in ../expectations.json, committed alone at c003c3b):
  a = 1.075 sqrt(Phi0/B) = 34.566 nm at B = 2.0 T, tolerance 3 per cent.
Primary estimator: a_fft, the first Bragg ring of the azimuthally averaged
|FFT|^2 of the raw conductance map -- independent of both the vortex detector
and the alpha complex.

Also calibrates all five estimators on a synthetic perfect triangular lattice
with matched N and A, and under increasing positional disorder, to show which
way each one is biased (the earlier run's 36.04 nm is predicted to be the
a_h1 disorder bias).

Data: Duhan et al., Nat. Commun. 16, 2100 (2025), Zenodo 10.5281/zenodo.14780459.
Run: prlimit --as=8589934592 -- .venv-tda/bin/python step0_re6zr.py
"""
import glob
import os

import numpy as np

import qf_lib as q

SEED = 20260920
H_KOE = 20.0
TOL = 0.03


def estimators(pts, area, box, a_ref):
    bars = q.alpha_bars(pts, (3 * a_ref) ** 2)
    d0, d1 = q.h0_deaths(bars), q.h1_deaths(bars)
    nn = q.nn_distance(pts, box, margin=a_ref)
    p6, nb = q.psi6_global(pts)
    return {
        "n": int(len(pts)),
        "a_density_nm": q.a_density(len(pts), area),
        "a_nn_nm": float(np.median(nn)) if nn.size else None,
        "a_h1_nm": float(np.sqrt(3) * np.median(d1)) if d1.size else None,
        "a_h0_nm": float(2 * np.median(d0)) if d0.size else None,
        "psi6": p6, "n_bonds": nb,
        "h0_iqr_over_median": q.spread_stats(d0)["iqr_over_median"],
    }


def synth_triangular(n_target, Lx, Ly, a):
    pts = []
    dy = a * np.sqrt(3) / 2
    ny = int(np.ceil(Ly / dy)) + 1
    nx = int(np.ceil(Lx / a)) + 2
    for j in range(ny):
        for i in range(nx):
            x = i * a + (j % 2) * a / 2
            y = j * dy
            if 0 <= x <= Lx and 0 <= y <= Ly:
                pts.append([x, y])
    return np.array(pts)


def main():
    assert q.assert_alpha_convention()
    B = H_KOE / 10.0
    a_target = q.a_tri(B)
    folder = os.path.join(q.STM_ROOT, f"{int(H_KOE)} kOe")
    files = sorted(glob.glob(os.path.join(folder, "*.sxm")))
    assert files, folder

    per_img, a_ffts = [], []
    for fn in files:
        imgs, meta = q.read_sxm(fn)
        img = imgs[("Input_7", "fwd")]
        px_nm = meta["range_x_m"] * 1e9 / meta["nx"]
        af, q1, _, _ = q.a_fft(img, px_nm, a_target)
        pts = q.detect_minima(img, px_nm, a_target)
        Lx, Ly = meta["range_x_m"] * 1e9, meta["range_y_m"] * 1e9
        est = estimators(pts, Lx * Ly, (0, Lx, 0, Ly), a_target)
        est.update(file=os.path.basename(fn), a_fft_nm=af, q1_inv_nm=q1,
                   px_nm=px_nm, scan_nm=Lx, nx=meta["nx"])
        per_img.append(est)
        a_ffts.append(af)
        print(f"{os.path.basename(fn)} a_fft={af:6.2f} a_nn={est['a_nn_nm']:6.2f} "
              f"a_dens={est['a_density_nm']:6.2f} a_h1={est['a_h1_nm']:6.2f} "
              f"a_h0={est['a_h0_nm']:6.2f} n={est['n']} psi6={est['psi6']:.3f}", flush=True)

    def agg(key):
        v = np.array([p[key] for p in per_img if p[key] is not None], dtype=float)
        v = v[np.isfinite(v)]
        return {"median": float(np.median(v)), "q1": float(np.percentile(v, 25)),
                "q3": float(np.percentile(v, 75)), "n": int(v.size)}

    a_fft_med = agg("a_fft_nm")["median"]
    rel = (a_fft_med - a_target) / a_target
    verdict = "PASS" if abs(rel) <= TOL else "FAIL"

    # ---- calibration: perfect lattice + disorder sweep, matched N and A
    rng = np.random.default_rng(SEED)
    Lx = per_img[0]["scan_nm"]
    base = synth_triangular(None, Lx, Lx, a_target)
    calib = []
    for sf in [0.0, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30]:
        pts = base + rng.normal(0, sf * a_target, base.shape)
        e = estimators(pts, Lx * Lx, (0, Lx, 0, Lx), a_target)
        e["sigma_over_a"] = sf
        calib.append(e)
        print(f"  synth sigma/a={sf:.2f} a_nn={e['a_nn_nm']:.2f} a_h1={e['a_h1_nm']:.2f} "
              f"a_h0={e['a_h0_nm']:.2f} psi6={e['psi6']:.3f} iqr/med={e['h0_iqr_over_median']:.3f}",
              flush=True)

    out = {
        "step": "STEP 0 known-answer reproduction",
        "field_kOe": H_KOE, "B_T": B,
        "target_a_nm": a_target,
        "target_source": "a = 1.075 sqrt(Phi0/B), Duhan et al., Nat. Commun. 16, 2100 (2025)",
        "tolerance_prestated": TOL,
        "primary_estimator": "a_fft (first Bragg ring of azimuthally averaged |FFT|^2; detector-free, TDA-free)",
        "a_fft": agg("a_fft_nm"),
        "a_fft_rel_error": float(rel),
        "VERDICT": verdict,
        "secondary": {k: agg(k) for k in ["a_density_nm", "a_nn_nm", "a_h1_nm", "a_h0_nm", "psi6", "h0_iqr_over_median", "n"]},
        "prior_run_value_nm": 36.04,
        "prior_run_estimator": "sqrt(3) * median(H1 alpha death) == a_h1 here",
        "n_images": len(files),
        "seed": SEED,
        "calibration_synthetic": calib,
        "per_image": per_img,
        "command": q.command(),
    }
    p = q.write_json("step0_re6zr_20kOe.json", out)
    print(f"\nTARGET {a_target:.3f} nm   a_fft median {a_fft_med:.3f} nm   "
          f"rel {rel*100:+.2f}%   tol {TOL*100:.0f}%   -> {verdict}")
    print("wrote", p)


if __name__ == "__main__":
    main()
