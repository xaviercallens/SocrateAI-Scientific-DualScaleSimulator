"""DATASET 1: Re6Zr STM vortex maps at every available field (3 - 70 kOe).

One TopoDB dataset per field value; ~20 conductance maps per field.
Alpha persistence on the vortex-core point clouds, plus the detector-free
FFT lattice constant and the psi6 orientational statistic (the inherited
binding limit: the H0 spread measures regular spacing, so it is always
reported paired with psi6).

Null (the published comparison): uniform random points in the same rectangle
at the same count. n_null sets per field, rank p-values.

Data: Duhan et al., Nat. Commun. 16, 2100 (2025), Zenodo 10.5281/zenodo.14780459.
Run: prlimit --as=8589934592 -- .venv-tda/bin/python compute_stm_fields.py
"""
import glob
import os
import sys

import numpy as np

import qf_lib as q

SEED = 20260920
N_NULL = 200
R_BETTI_FRAC = 0.55   # pre-stated Betti threshold, in units of a_tri(B)


def stats_for(pts, area, box, a_ref):
    bars = q.alpha_bars(pts, (3 * a_ref) ** 2)
    d0, d1 = q.h0_deaths(bars), q.h1_deaths(bars)
    p6, _ = q.psi6_global(pts)
    iv1 = bars[1]
    fin = np.isfinite(iv1[:, 1]) if len(iv1) else np.zeros(0, bool)
    tp1 = float(np.sum(iv1[fin, 1] - iv1[fin, 0])) / (a_ref ** 2) if len(iv1) else 0.0
    nn = q.nn_distance(pts, box, margin=a_ref)
    return {
        "bars": bars, "d0": d0, "d1": d1,
        "psi6": p6,
        "h0_iqr_over_median": q.spread_stats(d0)["iqr_over_median"],
        "h1_total_persistence_over_a2": tp1,
        "a_nn_nm": float(np.median(nn)) if nn.size else None,
        "a_h1_nm": float(np.sqrt(3) * np.median(d1)) if d1.size else None,
        "a_h0_nm": float(2 * np.median(d0)) if d0.size else None,
        "a_density_nm": q.a_density(len(pts), area),
        "n": int(len(pts)),
        "betti": q.betti_from_bars(bars, R_BETTI_FRAC * a_ref),
    }


def synth_triangular(Lx, Ly, a):
    pts, dy = [], a * np.sqrt(3) / 2
    for j in range(int(np.ceil(Ly / dy)) + 1):
        for i in range(int(np.ceil(Lx / a)) + 2):
            x, y = i * a + (j % 2) * a / 2, j * dy
            if 0 <= x <= Lx and 0 <= y <= Ly:
                pts.append([x, y])
    return np.array(pts)


def main():
    assert q.assert_alpha_convention()
    rng = np.random.default_rng(SEED)
    folders = sorted(glob.glob(os.path.join(q.STM_ROOT, "* kOe")),
                     key=lambda p: float(os.path.basename(p).split()[0]))
    out = {"meta": {"seed": SEED, "n_null": N_NULL, "r_betti_frac_of_a": R_BETTI_FRAC,
                    "command": q.command(), "script": os.path.abspath(__file__),
                    "detector": "plane-subtract, Gaussian LP sigma=a/6, disk-minimum radius a/3",
                    "source": "Duhan et al., Nat. Commun. 16, 2100 (2025); Zenodo 10.5281/zenodo.14780459",
                    "channel": "Input_7 forward (lock-in conductance)"},
           "fields": {}}
    for fo in folders:
        H = float(os.path.basename(fo).split()[0])
        B = H / 10.0
        a_ref = q.a_tri(B)
        files = sorted(glob.glob(os.path.join(fo, "*.sxm")))
        per, a_ffts = [], []
        pooled_d0, pooled_d1 = [], []
        clouds = {}
        for k, fn in enumerate(files):
            imgs, meta = q.read_sxm(fn)
            img = imgs[("Input_7", "fwd")]
            px = meta["range_x_m"] * 1e9 / meta["nx"]
            Lx, Ly = meta["range_x_m"] * 1e9, meta["range_y_m"] * 1e9
            af, q1, _, _ = q.a_fft(img, px, a_ref)
            pts = q.detect_minima(img, px, a_ref)
            s = stats_for(pts, Lx * Ly, (0, Lx, 0, Ly), a_ref)
            pooled_d0.append(s["d0"]); pooled_d1.append(s["d1"])
            clouds[f"img{k:02d}"] = pts
            per.append({"file": os.path.basename(fn), "a_fft_nm": af, "q1_inv_nm": q1,
                        "n": s["n"], "psi6": s["psi6"],
                        "h0_iqr_over_median": s["h0_iqr_over_median"],
                        "h1_total_persistence_over_a2": s["h1_total_persistence_over_a2"],
                        "a_nn_nm": s["a_nn_nm"], "a_h1_nm": s["a_h1_nm"],
                        "a_h0_nm": s["a_h0_nm"], "a_density_nm": s["a_density_nm"],
                        "betti0": s["betti"][0], "betti1": s["betti"][1],
                        "n_expected_flux": B * (Lx * Ly * 1e-18) / q.PHI0,
                        "scan_nm": Lx, "px_nm": px})
            a_ffts.append(af)
        npz = os.path.join(q.RESULTS, f"stm_cores_{int(H)}kOe.npz")
        np.savez_compressed(npz, **clouds)

        def med(key):
            v = np.array([p[key] for p in per if p[key] is not None], float)
            return float(np.median(v[np.isfinite(v)]))

        Lx = per[0]["scan_nm"]; area = Lx * Lx
        n_med = int(np.median([p["n"] for p in per]))
        # --- null: uniform random points, same count, same rectangle
        nul = {"h0_iqr_over_median": [], "psi6": [], "h1_total_persistence_over_a2": []}
        for _ in range(N_NULL):
            rp = rng.random((n_med, 2)) * Lx
            s = stats_for(rp, area, (0, Lx, 0, Lx), a_ref)
            for kk in nul:
                nul[kk].append(s[kk])
        # --- known-answer control: synthetic perfect triangular lattice, matched A
        sp = synth_triangular(Lx, Lx, a_ref)
        sctl = stats_for(sp, area, (0, Lx, 0, Lx), a_ref)

        obs = {k: med(k) for k in ["h0_iqr_over_median", "psi6", "h1_total_persistence_over_a2"]}
        pv = {"h0_iqr_over_median": q.rank_p(obs["h0_iqr_over_median"], nul["h0_iqr_over_median"], "less"),
              "psi6": q.rank_p(obs["psi6"], nul["psi6"], "greater"),
              "h1_total_persistence_over_a2": q.rank_p(obs["h1_total_persistence_over_a2"],
                                                       nul["h1_total_persistence_over_a2"], "two")}
        rec = {
            "H_kOe": H, "B_T": B, "a_tri_formula_nm": a_ref, "n_images": len(files),
            "scan_nm": Lx, "px_nm": per[0]["px_nm"], "area_nm2": area,
            "n_detected_median": n_med,
            "n_expected_flux": per[0]["n_expected_flux"],
            "count_over_flux": n_med / per[0]["n_expected_flux"],
            "a_fft_nm": {"median": float(np.median(a_ffts)),
                         "q1": float(np.percentile(a_ffts, 25)), "q3": float(np.percentile(a_ffts, 75))},
            "a_nn_nm": med("a_nn_nm"), "a_h1_nm": med("a_h1_nm"),
            "a_h0_nm": med("a_h0_nm"), "a_density_nm": med("a_density_nm"),
            "observed": obs, "p_values_rank": pv, "n_null": N_NULL,
            "null_median": {k: float(np.median(v)) for k, v in nul.items()},
            "betti_median": {0: int(np.median([p["betti0"] for p in per])),
                             1: int(np.median([p["betti1"] for p in per]))},
            "betti_synthetic_perfect": {0: sctl["betti"][0], 1: sctl["betti"][1]},
            "synthetic_control": {k: sctl[k] for k in
                                  ["psi6", "h0_iqr_over_median", "a_nn_nm", "a_h1_nm", "a_h0_nm", "n"]},
            "top_h1_bars": q.top_bars(stats_for(clouds["img00"], area, (0, Lx, 0, Lx), a_ref)["bars"], 1, 10),
            "top_h0_bars": q.top_bars(stats_for(clouds["img00"], area, (0, Lx, 0, Lx), a_ref)["bars"], 0, 10,
                                      r_trunc=3 * a_ref),
            "cores_npz": npz, "cores_sha256": q.sha256_file(npz),
            "per_image": per,
        }
        out["fields"][f"{H:g}kOe"] = rec
        print(f"H={H:5g} kOe a_tri={a_ref:6.2f} a_fft={rec['a_fft_nm']['median']:6.2f} "
              f"({100*(rec['a_fft_nm']['median']/a_ref-1):+5.1f}%) n={n_med}/{rec['n_expected_flux']:.0f} "
              f"psi6={obs['psi6']:.3f}(p={pv['psi6']:.4f}) iqr/med={obs['h0_iqr_over_median']:.3f}"
              f"(p={pv['h0_iqr_over_median']:.4f}) b={rec['betti_median']}", flush=True)

    # --- Abrikosov scaling fit a = C / sqrt(B)
    Bs = np.array([v["B_T"] for v in out["fields"].values()])
    As = np.array([v["a_fft_nm"]["median"] for v in out["fields"].values()])
    C = float(np.exp(np.mean(np.log(As) + 0.5 * np.log(Bs))))
    slope = float(np.polyfit(np.log(Bs), np.log(As), 1)[0])
    out["abrikosov_fit"] = {
        "model": "a_fft = C * B^slope, fitted by least squares in log-log",
        "slope": slope, "slope_expected": -0.5,
        "C_nm_T_half": C,
        "C_expected_nm_T_half": float(1.075 * np.sqrt(q.PHI0) * 1e9),
        "C_ratio": C / (1.075 * np.sqrt(q.PHI0) * 1e9),
        "per_field_ratio": {k: v["a_fft_nm"]["median"] / v["a_tri_formula_nm"] for k, v in out["fields"].items()},
    }
    print("\nAbrikosov fit: slope %.4f (expect -0.5), C ratio %.4f" % (slope, out["abrikosov_fit"]["C_ratio"]))
    print("wrote", q.write_json("stm_fields.json", out))


if __name__ == "__main__":
    sys.exit(main())
