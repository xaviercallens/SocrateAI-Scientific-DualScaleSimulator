#!/usr/bin/env python3
"""3-D comoving point clouds: DESI DR1 LSS catalogues and the SDSS bulk sample.

Method: GUDHI AlphaComplex over Z/2, max_dim 3, on a random subsample that fits
memory.  The filtration value is the squared circumradius; every bar is reported
after sqrt, i.e. in Mpc/h.

Primary statistic is the COUNT of bars above a fixed persistence, per
docs/FUTURE_OBSERVATIONAL_TARGETS.md §1: in three dimensions the spacing spread
that carried the 2-D Re6Zr result never reaches 95 % power, while the count
saturates at the lowest amplitude tested.  The spacing spread is still recorded,
as a declared SECONDARY.

THE NULL IS THE OFFICIAL RANDOM CATALOGUE AND IT IS UNCLUSTERED.  The randoms
share the survey window and n(z) but contain no galaxy clustering, so a p-value
against them answers "does this point set differ topologically from an
unclustered field with the same selection?" -- whose answer is known in advance
to be yes.  It is NOT a test for any anomaly with respect to LCDM.  The memo
states that substituting the randoms for a clustering-matched null "would
manufacture a detection"; the string is therefore carried in `null_model` itself.

Cosmology for comoving distances: flat LCDM, Om = 0.3137, h = 0.6736 (DESI DR1
fiducial / AbacusSummit base).  It sets the unit of length and nothing else.

Tier X (numerics).
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (Timer, add_run, alpha_bars, p_method_str, peak_mb, rank_p,  # noqa: E402
                    save_json, sha256_file, topodb)

SCRIPT = "topodb_runs/astro/run_pointclouds.py"
PY = "/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python"
SEED = 20260920
DESI_DIR = "/mnt/disks/disk-socrateai-local-1/dualscale-data-r3/desi_dr1_lss"
SDSS = "/mnt/disks/disk-socrateai-local-1/dualscale-data-r3/sdss_spec_bulk/sdss_spec_ra100_250_dec-10_50.csv"

OM, H0H = 0.3137, 0.6736          # flat LCDM; distances come out in Mpc/h
PERS_CUT = 5.0                    # Mpc/h, declared in PRE_DECLARED_STATISTICS.md
MAX_ALPHA_SQ = 900.0              # (30 Mpc/h)^2 -- declared hyperparameter, not the full Delaunay


def comoving_mpc_h(z):
    """Comoving distance in Mpc/h for flat LCDM, by Simpson quadrature.

    D_C = (c/H0) * int_0^z dz'/E(z'),  and c/H0 = 2997.92458/h Mpc = 2997.92458 Mpc/h.
    """
    z = np.asarray(z, dtype=float)
    grid = np.linspace(0.0, float(z.max()) * 1.0001 + 1e-9, 4096)
    E = np.sqrt(OM * (1 + grid) ** 3 + (1 - OM))
    integ = np.concatenate([[0.0], np.cumsum(np.diff(grid) * 0.5 * (1 / E[1:] + 1 / E[:-1]))])
    return 2997.92458 * np.interp(z, grid, integ)


def to_cartesian(ra, dec, z):
    d = comoving_mpc_h(z)
    r, dc = np.radians(ra), np.radians(dec)
    return np.c_[d * np.cos(dc) * np.cos(r), d * np.cos(dc) * np.sin(r), d * np.sin(dc)]


def cloud_stats(points):
    """All declared statistics for one 3-D cloud.  Returns (stats, bars)."""
    import gudhi
    ac = gudhi.AlphaComplex(points=np.ascontiguousarray(points, dtype=float))
    st = ac.create_simplex_tree(max_alpha_square=MAX_ALPHA_SQ)
    st.compute_persistence(homology_coeff_field=2, persistence_dim_max=True)
    d0, d1, d2 = (alpha_bars(st, k) for k in (0, 1, 2))
    fin0 = d0[np.isfinite(d0[:, 1]), 1]
    p1 = np.diff(d1, axis=1).ravel() if d1.size else np.zeros(0)
    p2 = np.diff(d2, axis=1).ravel() if d2.size else np.zeros(0)
    p1 = p1[np.isfinite(p1)]
    p2 = p2[np.isfinite(p2)]
    q1, q3 = np.percentile(fin0, [25, 75]) if fin0.size else (0.0, 0.0)
    med = float(np.median(fin0)) if fin0.size else 0.0
    s = {
        "n_h1_bars_over_5mpc": float(np.sum(p1 > PERS_CUT)),
        "n_h2_bars_over_5mpc": float(np.sum(p2 > PERS_CUT)),
        "total_persistence_h1": float(p1.sum()),
        "total_persistence_h2": float(p2.sum()),
        "max_persistence_h2": float(p2.max()) if p2.size else 0.0,
        "h0_death_iqr_over_median": float((q3 - q1) / med) if med else 0.0,
        "n_h1_bars_all": float(p1.size),
        "n_h2_bars_all": float(p2.size),
        "h0_death_median_mpc_h": med,
    }
    return s, {1: d1, 2: d2}


# statistics that are compared against the null, and the tail that means
# "more topological structure than an unclustered field"
TAILS = {"n_h1_bars_over_5mpc": "two", "n_h2_bars_over_5mpc": "two",
         "total_persistence_h1": "two", "total_persistence_h2": "two",
         "max_persistence_h2": "two", "h0_death_iqr_over_median": "two"}
NSTAT = len(TAILS)


def run_sample(db_id, title, data_pts, null_clouds, null_desc, n_sub, dataset_kw,
               notes_extra, out, n_null_label=None):
    """One dataset: one run on the data, statistics against the null clouds."""
    with topodb() as db:
        db.add_dataset(id=db_id, title=title, domain="astro", **dataset_kw)

    with Timer() as t:
        stats, bars = cloud_stats(data_pts)
    wall, mem = t.sec, peak_mb()
    print(f"  data:  {stats}")

    null_stats = []
    for i, c in enumerate(null_clouds):
        s, _ = cloud_stats(c)
        null_stats.append(s)
        print(f"  null {i + 1}/{len(null_clouds)}: n_h2={s['n_h2_bars_over_5mpc']:.0f} "
              f"n_h1={s['n_h1_bars_over_5mpc']:.0f} iqr/med={s['h0_death_iqr_over_median']:.4f}")
    n_null = len(null_stats)

    params = dict(n_points_subsample=int(n_sub), max_alpha_square=MAX_ALPHA_SQ,
                  persistence_cut_mpc_h=PERS_CUT, cosmology=f"flat LCDM Om={OM} h={H0H}",
                  units="Mpc/h comoving", coeff_field=2, n_null=n_null,
                  weights="NOT applied: unweighted galaxy positions (stated, not assumed)")
    cmd = f"prlimit --as=8589934592 -- {PY} {SCRIPT}"
    with topodb() as db:
        rid = add_run(db, dataset_id=db_id, method="alpha", coeff_field=2, max_dim=3,
                      params=params, script=SCRIPT, command=cmd, tier="X",
                      preprocessing=notes_extra, seed=str(SEED), wall_sec=wall, peak_mb=mem)
        for dim, d in bars.items():
            if d.size:
                db.add_bars(rid, dim, [(b, dd) for b, dd in d], top=25)
        # Betti WITHOUT `expected`: no textbook answer exists for a galaxy field,
        # and the gap rule is not a detector (PRE_DECLARED_STATISTICS.md A1).
        db.add_betti(rid, {1: int(stats["n_h1_bars_over_5mpc"]),
                           2: int(stats["n_h2_bars_over_5mpc"])})
        for name, val in stats.items():
            if name in TAILS and n_null > 0:
                nv = [s[name] for s in null_stats]
                p, floor = rank_p(val, nv, TAILS[name])
                db.add_statistic(rid, name, val, null_model=null_desc, n_null=n_null,
                                 p_value=p, p_method=p_method_str(n_null, TAILS[name]),
                                 multiplicity=f"Bonferroni/{NSTAT} -> threshold {0.05 / NSTAT:.4g}")
            else:
                db.add_statistic(rid, name, val)     # no null -> NO p-value
        if n_null >= 3:
            # null calibration: hold one null out, rank it against the rest
            held = [s["n_h2_bars_over_5mpc"] for s in null_stats]
            p_held, _ = rank_p(held[0], held[1:], "two")
            db.add_control(rid, "null_calibration",
                           "a held-out null realisation ranked against the remaining n-1 is not extreme",
                           passed=bool(p_held > 0.05),
                           detail=f"held-out n_h2_bars_over_5mpc rank p = {p_held:.3f} on {n_null - 1} others")
        db.add_control(rid, "known_answer",
                       "Step 0 (runs 7-12) recovered circle/sphere/torus/wells/HEALPix on this same "
                       "AlphaComplex pipeline before any real data was read",
                       passed=True, detail="see topodb_runs/astro/results/step0_known_answers.json")
    out[db_id] = {"run_id": rid, "data": stats, "null": null_stats, "n_null": n_null,
                  "p_values": {n: rank_p(stats[n], [s[n] for s in null_stats], TAILS[n])[0]
                               for n in TAILS} if n_null else {}}
    return rid, stats, null_stats


def desi(tracer, cap, n_sub, n_null, zmin, zmax, out):
    from astropy.io import fits
    base = f"{tracer}_{cap}"
    dfile = f"{DESI_DIR}/{base}_clustering.dat.fits"
    rfile = f"{DESI_DIR}/{base}_0_clustering.ran.fits"
    print(f"\n=== DESI {base} ===")
    rng = np.random.default_rng(SEED)
    with fits.open(dfile, memmap=True) as h:
        d = h[1].data
        ra, dec, z = np.asarray(d["RA"], float), np.asarray(d["DEC"], float), np.asarray(d["Z"], float)
    keep = (z > zmin) & (z < zmax) & np.isfinite(z)
    ra, dec, z = ra[keep], dec[keep], z[keep]
    n_tot = ra.size
    idx = rng.choice(n_tot, size=min(n_sub, n_tot), replace=False)
    pts = to_cartesian(ra[idx], dec[idx], z[idx])
    n_sub_eff = pts.shape[0]
    print(f"  {n_tot} galaxies in {zmin}<z<{zmax}; subsample {n_sub_eff}")

    with fits.open(rfile, memmap=True) as h:
        r = h[1].data
        rra, rdec, rz = np.asarray(r["RA"], float), np.asarray(r["DEC"], float), np.asarray(r["Z"], float)
    rk = (rz > zmin) & (rz < zmax) & np.isfinite(rz)
    rra, rdec, rz = rra[rk], rdec[rk], rz[rk]
    perm = rng.permutation(rra.size)[: n_sub_eff * n_null]       # DISJOINT subsamples
    nulls = [to_cartesian(rra[c], rdec[c], rz[c])
             for c in np.array_split(perm, n_null)]
    print(f"  {rra.size} randoms; {n_null} disjoint subsamples of {nulls[0].shape[0]}")

    null_desc = (f"{n_null} DISJOINT subsamples of the official DESI DR1 random catalogue "
                 f"{os.path.basename(rfile)}, same size, same window, Z drawn from the data n(z). "
                 "THE RANDOMS ARE UNCLUSTERED (Poisson). This p-value tests departure from an "
                 "unclustered field with the same selection function, NOT departure from LCDM "
                 "clustering, and is NOT evidence of a topological anomaly. A clustering-matched "
                 "DESI-footprint mock family does not exist on this disk (see "
                 "docs/FUTURE_OBSERVATIONAL_TARGETS.md §2a) and was NOT ATTEMPTED.")
    kw = dict(source=("DESI DR1 LSS catalogues v1.5, https://data.desi.lbl.gov/public/dr1/survey/catalogs/"
                      "dr1/LSS/iron/LSScats/v1.5/ ; staged copy, sha256 in audit/data_r3_manifest.json"),
              provenance="observation", n_objects=int(n_tot), ambient_dim=3, units="Mpc/h comoving",
              sha256=sha256_file(dfile), local_path=dfile,
              notes=(f"{tracer} {cap}. z cut {zmin}<z<{zmax}. Official randoms "
                     f"{os.path.basename(rfile)} used as the (unclustered) null."))
    return run_sample(f"astro/desi_dr1_{tracer.lower().replace('-', '_').replace('.', '')}_{cap.lower()}",
                      f"DESI DR1 {tracer} {cap} (comoving point cloud)", pts, nulls, null_desc,
                      n_sub_eff, kw,
                      (f"z cut {zmin}<z<{zmax}; uniform random subsample of {n_sub_eff} of {n_tot} "
                       f"galaxies (seed {SEED}); RA/Dec/z -> comoving xyz with flat LCDM Om={OM}; "
                       "WEIGHT columns NOT applied"), out)


def sdss(n_sub, n_null, out):
    print("\n=== SDSS spectroscopic bulk ===")
    rng = np.random.default_rng(SEED)
    d = np.loadtxt(SDSS, delimiter=",", skiprows=1)
    ra, dec, z = d[:, 0], d[:, 1], d[:, 2]
    keep = (z > 0.05) & (z < 0.20) & np.isfinite(z)     # the MGS-dominated slab
    ra, dec, z = ra[keep], dec[keep], z[keep]
    n_tot = ra.size
    idx = rng.choice(n_tot, size=min(n_sub, n_tot), replace=False)
    pts = to_cartesian(ra[idx], dec[idx], z[idx])
    print(f"  {n_tot} galaxies in 0.05<z<0.20 of {d.shape[0]} total; subsample {pts.shape[0]}")

    # NO random catalogue exists for this sample.  The declared null is a
    # SHUFFLED-REDSHIFT null: angular positions kept, z resampled from the
    # sample's own n(z).  It destroys 3-D clustering while keeping the angular
    # mask and the radial selection -- weaker than a mock, and labelled as such.
    nulls = []
    for _ in range(n_null):
        j = rng.permutation(n_tot)[: pts.shape[0]]
        k = rng.permutation(n_tot)[: pts.shape[0]]
        nulls.append(to_cartesian(ra[j], dec[j], z[k]))
    null_desc = (f"{n_null} SHUFFLED-REDSHIFT realisations: angular positions and redshifts drawn "
                 "independently from the sample itself, so the angular footprint and the radial n(z) "
                 "are preserved and all 3-D clustering is destroyed. No random catalogue and no "
                 "completeness map exists for this sample, so this is the only null available; it is "
                 "WEAKER than a mock and is not a clustering-matched null.")
    kw = dict(source=("SDSS SkyServer spectroscopic query, 100 tiles over RA 100-250, Dec -10..50; "
                      "staged as sdss_spec_ra100_250_dec-10_50.csv. PROVENANCE PARTIAL: the original "
                      "query script was not found on disk (audit/data_r3_manifest.json); the staged "
                      "file was cross-matched against an independent DR17 catalogue, 4367/5042 "
                      "matching within 2 arcsec at identical z."),
              provenance="observation", n_objects=int(d.shape[0]), ambient_dim=3, units="Mpc/h comoving",
              sha256=sha256_file(SDSS), local_path=SDSS,
              notes=("HETEROGENEOUS SELECTION: MGS at low z, LRG-like above z~0.15, no zErr column, "
                     "no completeness map, no randoms. The z<0.20 slab used here is MGS-dominated but "
                     "is still not a single well-defined selection."))
    return run_sample("astro/sdss_spec_bulk", "SDSS spectroscopic bulk sample (comoving point cloud)",
                      pts, nulls, null_desc, pts.shape[0], kw,
                      ("0.05<z<0.20 slab; uniform random subsample (seed 20260920); RA/Dec/z -> "
                       f"comoving xyz with flat LCDM Om={OM}; heterogeneous selection, unweighted"), out)


def main():
    n_sub = int(os.environ.get("NSUB", 25000))
    n_null = int(os.environ.get("NNULL", 5))
    out = {"seed": SEED, "n_sub": n_sub, "n_null": n_null}
    which = sys.argv[1:] or ["bgs215ngc", "bgs215sgc", "bgsngc", "lrgngc", "sdss"]
    jobs = {
        "bgs215ngc": lambda: desi("BGS_BRIGHT-21.5", "NGC", n_sub, n_null, 0.10, 0.30, out),
        "bgs215sgc": lambda: desi("BGS_BRIGHT-21.5", "SGC", n_sub, n_null, 0.10, 0.30, out),
        "bgsngc": lambda: desi("BGS_BRIGHT", "NGC", n_sub, n_null, 0.10, 0.40, out),
        "lrgngc": lambda: desi("LRG", "NGC", n_sub, n_null, 0.40, 1.10, out),
        "sdss": lambda: sdss(n_sub, n_null, out),
    }
    for w in which:
        jobs[w]()
        save_json(f"pointclouds_{w}.json", out)
    print("\n" + "=" * 60)
    for k, v in out.items():
        if isinstance(v, dict) and "p_values" in v:
            print(k, "run", v["run_id"], {a: round(b, 4) for a, b in v["p_values"].items()})
    return 0


if __name__ == "__main__":
    sys.exit(main())
