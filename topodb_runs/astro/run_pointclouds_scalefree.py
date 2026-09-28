#!/usr/bin/env python3
"""Re-run the two DESI samples whose declared alpha cap was too tight, with a
SCALE-FREE cap (AMENDMENT A3 in PRE_DECLARED_STATISTICS.md).

The absolute cap max_alpha_square = 900 (30 Mpc/h) is 2.1x the LRG sample's own
H0 death median, so the alpha complex was truncated before voids of the sample's
characteristic size could be born: the data AND all 20 nulls had exactly zero H2
bars over 5 Mpc/h, a statistic with no variance in the null.

Here the cap is (6 * H0-death-median)^2, measured on the data by a pilot pass, so
it is set by the sample and not by the analyst.  The persistence cut for the
COUNT statistic is likewise scaled to 1.1 * H0-death-median instead of a fixed
5 Mpc/h, for the same reason.

These are recorded as SEPARATE runs; runs 144 and 146 are not deleted.

Tier X.
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_pointclouds as RP  # noqa: E402
from common import (Timer, add_run, alpha_bars, p_method_str, peak_mb, rank_p,  # noqa: E402
                    save_json, sha256_file, topodb)

SCRIPT = "topodb_runs/astro/run_pointclouds_scalefree.py"
PY = RP.PY
SEED = RP.SEED


def pilot_scale(points):
    """H0 death median of a 15000-point pilot: the sample's own length unit."""
    import gudhi
    rng = np.random.default_rng(SEED)
    i = rng.choice(points.shape[0], size=min(15000, points.shape[0]), replace=False)
    ac = gudhi.AlphaComplex(points=np.ascontiguousarray(points[i], dtype=float))
    st = ac.create_simplex_tree()
    st.compute_persistence(homology_coeff_field=2, persistence_dim_max=False)
    d0 = alpha_bars(st, 0)
    fin = d0[np.isfinite(d0[:, 1]), 1]
    return float(np.median(fin))


def main():
    from astropy.io import fits
    n_sub = int(os.environ.get("NSUB", 60000))
    n_null = int(os.environ.get("NNULL", 20))
    out = {"seed": SEED, "amendment": "A3 scale-free alpha cap"}

    for tracer, cap_name, zmin, zmax in [("LRG", "NGC", 0.40, 1.10),
                                         ("BGS_BRIGHT", "NGC", 0.10, 0.40)]:
        base = f"{tracer}_{cap_name}"
        dfile = f"{RP.DESI_DIR}/{base}_clustering.dat.fits"
        rfile = f"{RP.DESI_DIR}/{base}_0_clustering.ran.fits"
        print(f"\n=== {base} (scale-free cap) ===")
        rng = np.random.default_rng(SEED)
        with fits.open(dfile, memmap=True) as h:
            d = h[1].data
            ra, dec, z = (np.asarray(d[c], float) for c in ("RA", "DEC", "Z"))
        k = (z > zmin) & (z < zmax) & np.isfinite(z)
        ra, dec, z = ra[k], dec[k], z[k]
        n_tot = ra.size
        idx = rng.choice(n_tot, size=min(n_sub, n_tot), replace=False)
        pts = RP.to_cartesian(ra[idx], dec[idx], z[idx])
        n_eff = pts.shape[0]

        scale = pilot_scale(pts)
        max_a2 = (6.0 * scale) ** 2
        pcut = 1.1 * scale
        print(f"  {n_tot} galaxies, subsample {n_eff}; pilot H0 death median = {scale:.2f} Mpc/h")
        print(f"  scale-free cap = (6 x {scale:.2f})^2 = {max_a2:.0f}  ({6 * scale:.1f} Mpc/h); "
              f"persistence cut = {pcut:.2f} Mpc/h")

        old_a2, old_cut = RP.MAX_ALPHA_SQ, RP.PERS_CUT
        RP.MAX_ALPHA_SQ, RP.PERS_CUT = max_a2, pcut
        try:
            with Timer() as t:
                sd, bars = RP.cloud_stats(pts)
            print(f"  data: n_h2={sd['n_h2_bars_over_5mpc']:.0f} n_h1={sd['n_h1_bars_over_5mpc']:.0f} "
                  f"max_pers_h2={sd['max_persistence_h2']:.2f} ({t.sec:.0f}s)")

            with fits.open(rfile, memmap=True) as h:
                r = h[1].data
                rra, rdec, rz = (np.asarray(r[c], float) for c in ("RA", "DEC", "Z"))
            rk = (rz > zmin) & (rz < zmax) & np.isfinite(rz)
            rra, rdec, rz = rra[rk], rdec[rk], rz[rk]
            perm = rng.permutation(rra.size)[: n_eff * n_null]
            nulls = []
            for j, c in enumerate(np.array_split(perm, n_null)):
                s, _ = RP.cloud_stats(RP.to_cartesian(rra[c], rdec[c], rz[c]))
                nulls.append(s)
                if (j + 1) % 5 == 0:
                    print(f"    null {j + 1}/{n_null}: n_h2={s['n_h2_bars_over_5mpc']:.0f}")
        finally:
            RP.MAX_ALPHA_SQ, RP.PERS_CUT = old_a2, old_cut
        wall, mem = t.sec, peak_mb()

        ds = f"astro/desi_dr1_{tracer.lower()}_{cap_name.lower()}"
        null_desc = (f"{n_null} DISJOINT subsamples of the official DESI DR1 random catalogue "
                     f"{os.path.basename(rfile)}, same size, same window, Z from the data n(z). "
                     "THE RANDOMS ARE UNCLUSTERED (Poisson): this tests departure from an unclustered "
                     "field with the same selection, NOT departure from LCDM clustering, and is NOT "
                     "evidence of a topological anomaly. No clustering-matched DESI mock family "
                     "exists on this disk; NOT ATTEMPTED.")
        with topodb() as db:
            rid = add_run(db, dataset_id=ds, method="alpha", coeff_field=2, max_dim=3,
                          params=dict(n_points_subsample=n_eff, max_alpha_square=max_a2,
                                      alpha_cap_mpc_h=6.0 * scale,
                                      cap_rule="AMENDMENT A3: (6 x pilot H0-death-median)^2, "
                                               "scale-free; the sample sets its own unit",
                                      pilot_h0_death_median_mpc_h=scale,
                                      persistence_cut_mpc_h=pcut,
                                      persistence_cut_rule="1.1 x pilot H0-death-median",
                                      cosmology=f"flat LCDM Om={RP.OM} h={RP.H0H}",
                                      units="Mpc/h comoving", n_null=n_null,
                                      supersedes_note="a SEPARATE run; the absolute-cap run is kept"),
                          script=SCRIPT, command=f"prlimit --as=8589934592 -- {PY} {SCRIPT}",
                          tier="X",
                          preprocessing=(f"z cut {zmin}<z<{zmax}; uniform subsample of {n_eff} of "
                                         f"{n_tot} (seed {SEED}); RA/Dec/z -> comoving xyz; "
                                         "WEIGHT columns NOT applied"),
                          seed=str(SEED), wall_sec=wall, peak_mb=mem)
            for dim, dd in bars.items():
                if dd.size:
                    db.add_bars(rid, dim, [(a, b) for a, b in dd], top=25)
            db.add_betti(rid, {1: int(sd["n_h1_bars_over_5mpc"]), 2: int(sd["n_h2_bars_over_5mpc"])})
            nout = 0
            for name in RP.TAILS:
                nv = np.array([s[name] for s in nulls], float)
                p, _ = rank_p(sd[name], nv, "two")
                db.add_statistic(rid, name, sd[name], null_model=null_desc, n_null=n_null,
                                 p_value=p, p_method=p_method_str(n_null, "two"),
                                 multiplicity=f"Bonferroni/{RP.NSTAT} -> {0.05 / RP.NSTAT:.4g}")
                s_ = float(nv.std(ddof=1))
                if s_ > 0:
                    db.add_statistic(rid, f"null_separation_sigma__{name}",
                                     float((sd[name] - nv.mean()) / s_))
                o = bool(sd[name] < nv.min() or sd[name] > nv.max())
                nout += int(o)
                db.add_statistic(rid, f"data_outside_null_range__{name}", float(o))
            for name, val in sd.items():
                if name not in RP.TAILS:
                    db.add_statistic(rid, name, val)
            db.add_control(rid, "negative",
                           "AMENDMENT A3: the absolute cap of 30 Mpc/h was only "
                           f"{30.0 / scale:.1f}x this sample's own H0 death median, which saturated "
                           "its H2 statistics (data and all 20 nulls had exactly zero H2 bars over "
                           "5 Mpc/h - no variance in the null). This run uses a scale-free cap.",
                           passed=True,
                           detail=(f"pilot H0 death median {scale:.2f} Mpc/h; new cap "
                                   f"{6 * scale:.1f} Mpc/h; new max H2 persistence "
                                   f"{sd['max_persistence_h2']:.2f} Mpc/h; "
                                   f"null n_h2 spread {np.std([s['n_h2_bars_over_5mpc'] for s in nulls], ddof=1):.2f}"))
            db.add_control(rid, "known_answer",
                           "Step 0 (runs 7-12) recovered circle/sphere/torus/wells/HEALPix on this "
                           "same AlphaComplex pipeline before any real data was read",
                           passed=True, detail="results/step0_known_answers.json")
        out[ds + "_scalefree"] = {"run_id": rid, "data": sd, "null": nulls, "n_null": n_null,
                                  "scale": scale, "max_alpha_square": max_a2, "n_outside": nout,
                                  "p_values": {n: rank_p(sd[n], [s[n] for s in nulls], "two")[0]
                                               for n in RP.TAILS}}
        print(f"  run {rid}: {nout}/{len(RP.TAILS)} outside null range")
        save_json("pointclouds_scalefree.json", out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
