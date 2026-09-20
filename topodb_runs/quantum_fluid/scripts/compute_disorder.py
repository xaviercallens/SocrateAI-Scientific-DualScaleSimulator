"""DATASET 5: a calibrated disorder axis, from 'ordered lattice' to 'random'.

Base: the real 20 kOe Re6Zr vortex-core clouds (the most ordered field measured,
psi6 = 0.641), degraded in controlled steps. This gives TopoDB an axis every
other domain can be compared against.

TWO axes, because they are physically different degradations:
  A) positional noise: every core displaced by an isotropic Gaussian of width
     sigma, sigma/a = 0, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30. Count unchanged.
  B) spurious extra points: a fraction f of uniformly-placed extra points added,
     f = 0, 0.02, 0.05, 0.10, 0.20, 0.40. Positions of the real cores unchanged.
     This axis exists because the Step-0 analysis found the real maps' apparent
     density and their measured spacing disagree, and axis A cannot produce that
     split (positional noise leaves the areal density untouched).

PRE-STATED PREDICTION (expectations.json, committed at c003c3b, before this ran):
  "psi6 degrades FIRST: it should fall to near the uniform-random level by
   sigma/a ~ 0.15, while h0_iqr_over_median stays close to its ordered value
   until sigma/a ~ 0.10 and only then rises toward the Poisson value 0.70."

Run: prlimit --as=8589934592 -- .venv-tda/bin/python compute_disorder.py
"""
import os

import numpy as np

import qf_lib as q

SEED = 20260920
N_NULL = 200
SIGMAS = [0.0, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30]
SPURIOUS = [0.0, 0.02, 0.05, 0.10, 0.20, 0.40]
H_KOE = 20.0


def measure(pts, area, box, a_ref):
    bars = q.alpha_bars(pts, (3 * a_ref) ** 2)
    d0, d1 = q.h0_deaths(bars), q.h1_deaths(bars)
    nn = q.nn_distance(pts, box, margin=a_ref)
    p6, _ = q.psi6_global(pts)
    return {"n": int(len(pts)), "psi6": p6,
            "h0_iqr_over_median": q.spread_stats(d0)["iqr_over_median"],
            "a_nn_nm": float(np.median(nn)) if nn.size else None,
            "a_h1_nm": float(np.sqrt(3) * np.median(d1)) if d1.size else None,
            "a_h0_nm": float(2 * np.median(d0)) if d0.size else None,
            "a_density_nm": q.a_density(len(pts), area),
            "betti0": q.betti_from_bars(bars, 0.55 * a_ref)[0],
            "betti1": q.betti_from_bars(bars, 0.55 * a_ref)[1],
            "bars": bars}


def main():
    assert q.assert_alpha_convention()
    rng = np.random.default_rng(SEED)
    a_ref = q.a_tri(H_KOE / 10.0)
    npz = np.load(os.path.join(q.RESULTS, f"stm_cores_{int(H_KOE)}kOe.npz"))
    L = 353.0
    area = L * L
    box = (0, L, 0, L)

    # --- uniform-random reference level (the 'random' end of the axis)
    nul = {"psi6": [], "h0_iqr_over_median": []}
    n0 = int(np.median([len(npz[k]) for k in npz.files]))
    for _ in range(N_NULL):
        rp = rng.random((n0, 2)) * L
        m = measure(rp, area, box, a_ref)
        nul["psi6"].append(m["psi6"]); nul["h0_iqr_over_median"].append(m["h0_iqr_over_median"])
    ref = {k: float(np.median(v)) for k, v in nul.items()}

    out = {"meta": {"seed": SEED, "n_null": N_NULL, "base": f"Re6Zr 20 kOe cores, {len(npz.files)} images",
                    "a_ref_nm": a_ref, "scan_nm": L, "command": q.command(),
                    "script": os.path.abspath(__file__),
                    "prediction_prestated": "psi6 degrades first; h0_iqr_over_median stays flat until "
                                            "sigma/a ~ 0.10 (expectations.json, commit c003c3b)"},
           "uniform_random_reference": ref, "n_null": N_NULL,
           "axis_A_positional_noise": {}, "axis_B_spurious_points": {}}

    for sf in SIGMAS:
        rows, first = [], None
        for k in npz.files:
            p = npz[k] + rng.normal(0, sf * a_ref, npz[k].shape)
            m = measure(p, area, box, a_ref)
            if first is None:
                first = m
            rows.append(m)
        def med(kk):
            v = np.array([r[kk] for r in rows if r[kk] is not None], float)
            return float(np.median(v[np.isfinite(v)]))
        rec = {kk: med(kk) for kk in ["n", "psi6", "h0_iqr_over_median", "a_nn_nm",
                                      "a_h1_nm", "a_h0_nm", "a_density_nm", "betti0", "betti1"]}
        rec["sigma_over_a"] = sf
        rec["psi6_p_rank"] = q.rank_p(rec["psi6"], nul["psi6"], "greater")
        rec["h0_iqr_p_rank"] = q.rank_p(rec["h0_iqr_over_median"], nul["h0_iqr_over_median"], "less")
        rec["psi6_frac_of_ordered"] = rec["psi6"] / out["axis_A_positional_noise"]["0.00"]["psi6"] \
            if "0.00" in out["axis_A_positional_noise"] else 1.0
        rec["top_h1_bars"] = q.top_bars(first["bars"], 1, 10)
        out["axis_A_positional_noise"][f"{sf:.2f}"] = rec
        print(f"A sigma/a={sf:.2f} n={rec['n']:.0f} psi6={rec['psi6']:.3f}(p={rec['psi6_p_rank']:.4f}) "
              f"iqr/med={rec['h0_iqr_over_median']:.3f}(p={rec['h0_iqr_p_rank']:.4f}) "
              f"a_nn={rec['a_nn_nm']:.2f} a_h1={rec['a_h1_nm']:.2f} a_dens={rec['a_density_nm']:.2f}",
              flush=True)

    for f in SPURIOUS:
        rows, first = [], None
        for k in npz.files:
            base = npz[k]
            nadd = int(round(f * len(base)))
            extra = rng.random((nadd, 2)) * L
            p = np.vstack([base, extra]) if nadd else base
            m = measure(p, area, box, a_ref)
            if first is None:
                first = m
            rows.append(m)
        def med(kk):
            v = np.array([r[kk] for r in rows if r[kk] is not None], float)
            return float(np.median(v[np.isfinite(v)]))
        rec = {kk: med(kk) for kk in ["n", "psi6", "h0_iqr_over_median", "a_nn_nm",
                                      "a_h1_nm", "a_h0_nm", "a_density_nm", "betti0", "betti1"]}
        rec["spurious_fraction"] = f
        rec["psi6_p_rank"] = q.rank_p(rec["psi6"], nul["psi6"], "greater")
        rec["h0_iqr_p_rank"] = q.rank_p(rec["h0_iqr_over_median"], nul["h0_iqr_over_median"], "less")
        rec["top_h1_bars"] = q.top_bars(first["bars"], 1, 10)
        out["axis_B_spurious_points"][f"{f:.2f}"] = rec
        print(f"B f={f:.2f} n={rec['n']:.0f} psi6={rec['psi6']:.3f} "
              f"iqr/med={rec['h0_iqr_over_median']:.3f} a_nn={rec['a_nn_nm']:.2f} "
              f"a_h1={rec['a_h1_nm']:.2f} a_dens={rec['a_density_nm']:.2f} "
              f"(a_h1/a_dens={rec['a_h1_nm']/rec['a_density_nm']:.3f})", flush=True)

    # --- which statistic degrades first?
    A = out["axis_A_positional_noise"]
    p0 = A["0.00"]["psi6"]; i0 = A["0.00"]["h0_iqr_over_median"]
    def first_cross(key, frac_to_random):
        for sf in SIGMAS[1:]:
            r = A[f"{sf:.2f}"]
            if key == "psi6":
                t = p0 + frac_to_random * (ref["psi6"] - p0)
                if r["psi6"] <= t:
                    return sf
            else:
                t = i0 + frac_to_random * (ref["h0_iqr_over_median"] - i0)
                if r["h0_iqr_over_median"] >= t:
                    return sf
        return None
    out["degradation_order"] = {
        "sigma_at_which_psi6_reaches_half_way_to_random": first_cross("psi6", 0.5),
        "sigma_at_which_h0_iqr_reaches_half_way_to_random": first_cross("h0_iqr", 0.5),
        "sigma_at_which_psi6_reaches_90pc_of_the_way_to_random": first_cross("psi6", 0.9),
        "sigma_at_which_h0_iqr_reaches_90pc_of_the_way_to_random": first_cross("h0_iqr", 0.9),
        "ordered_values": {"psi6": p0, "h0_iqr_over_median": i0},
        "random_values": ref,
    }
    print("\ndegradation order:", out["degradation_order"])
    print("wrote", q.write_json("disorder_series.json", out))


if __name__ == "__main__":
    main()
