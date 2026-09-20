#!/usr/bin/env python3
"""Sky maps: lower-star Betti curves on the FIXED HEALPix 2-complex.

Maps: Planck SMICA (CMB), WMAP 9-yr ILC (CMB), Haslam 408 MHz (Galactic
synchrotron), WMAP 9-yr K-band 23 GHz (foreground-dominated), COBE-DMR 53 GHz A
(1990s CMB, repixelised).  All use lib/cmb_topology.py from loop/tda-simple
d8175f1 -- audit/reverse_zero/E5-cmb-tda/cmb_tda.py is defective and is NOT
imported.

THE NULL AND ITS ONE REAL TRAP.  Sims are drawn from each map's OWN pseudo-C_ell
measured on the masked map, then masked with the same mask.  A masked map's
pseudo-C_ell is suppressed by roughly f_sky; synthesising from it and masking
again gives sims with LESS power than the data, which makes the data look
anomalous for free.  The pseudo-C_ell is therefore divided by f_sky before
synfast, and a `null_calibration` control holds one sim out and ranks it against
the remaining n-1.  Without that control none of these p-values is
interpretable.

Rank p-values only: the chi2 branch of the fixed library remains mildly
anti-conservative (0.058 against 0.05).

THE FOREGROUND MAPS ARE THE POSITIVE CONTROL.  Haslam and WMAP K are strongly
non-Gaussian, so a Gaussian null must be rejected at the rank floor.  If it is
not, the CMB nulls are vacuous and must be read as "the pipeline cannot fire",
not as "no signal".

Tier X (numerics).
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (ASTRO, Timer, add_run, p_method_str, peak_mb, rank_p,  # noqa: E402
                    save_json, sha256_file, topodb)

SCRIPT = "topodb_runs/astro/run_skymaps.py"
PY = "/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python"
SEED = 20260920
R3 = "/mnt/disks/disk-socrateai-local-1/dualscale-data-r3"
REV = "/mnt/disks/disk-socrateai-local-1/dualscale-wt-reverse/data/real2/cmb"
DATA = os.path.join(ASTRO, "data")

NU = np.linspace(-3.0, 3.0, 41)
STATS_TAIL = {"b0_at_nu_1": "two", "b1_at_nu_0": "two", "b0_curve_integral": "two",
              "b1_curve_integral": "two", "euler_char_at_nu_0": "two"}
NSTAT = len(STATS_TAIL) + 2          # + the two coarse curve-level rank p's


# --------------------------------------------------------------------- loaders
def _healpix_from_parquet(path, column, order_in, nside_out):
    import healpy as hp
    import pandas as pd
    col = pd.read_parquet(path, columns=[column])[column].values
    if col.dtype == object:                       # row-packed arrays
        m = np.stack(col).ravel().astype(float)
    else:
        m = np.asarray(col, float)
    return hp.ud_grade(m, nside_out, order_in=order_in, order_out="RING"), m.size


def load_map(key, nside):
    """Returns (map at `nside` RING, mask at `nside` RING, meta dict)."""
    import healpy as hp
    npix = hp.nside2npix(nside)
    if key == "planck_smica":
        f = f"{R3}/planck_maps/COM_CMB_IQU-smica_2048_R3.00_full.fits"
        mf = f"{R3}/planck_maps/COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits"
        m = hp.read_map(f, field=0, dtype=np.float32)
        msk = hp.read_map(mf, field=0, dtype=np.float32)
        m = hp.ud_grade(np.asarray(m, float), nside)
        msk = (hp.ud_grade(np.asarray(msk, float), nside) > 0.9).astype(float)
        meta = dict(file=f, mask_file=mf, native_nside=2048, units="K_CMB",
                    source="Planck PR3 COM_CMB_IQU-smica_2048_R3.00_full.fits + "
                           "COM_Mask_CMB-common-Mask-Int_2048_R3.00 (pla.esac.esa.int)",
                    mask_desc="Planck common intensity mask, ud_grade to nside then thresholded at 0.9")
    elif key == "wmap_ilc":
        f, mf = f"{REV}/wmap_ilc_9yr_v5.fits", f"{REV}/wmap_temperature_kq85_analysis_mask_r9_9yr_v5.fits"
        m = hp.read_map(f, field=0, dtype=np.float32)
        msk = hp.read_map(mf, field=0, dtype=np.float32)
        m = hp.ud_grade(np.asarray(m, float), nside)
        msk = (hp.ud_grade(np.asarray(msk, float), nside) > 0.9).astype(float)
        meta = dict(file=f, mask_file=mf, native_nside=512, units="mK_CMB",
                    source="WMAP 9-yr ILC map wmap_ilc_9yr_v5.fits + KQ85 r9 analysis mask "
                           "(lambda.gsfc.nasa.gov)",
                    mask_desc="WMAP KQ85 9-yr analysis mask, ud_grade then thresholded at 0.9")
    elif key == "haslam408":
        f = f"{DATA}/haslam-408mhz.parquet"
        m, npix_native = _healpix_from_parquet(f, "TEMPERATURE", "RING", nside)
        msk = np.ones(npix)
        meta = dict(file=f, native_nside=512, units="K (brightness temperature)",
                    source="https://huggingface.co/datasets/astro-legacy-archive/haslam-408mhz "
                           "(haslam-408mhz.parquet); destriped source-subtracted 'dsds' "
                           "reprocessing of Haslam et al., Remazeilles et al. 2015, "
                           "DOI 10.57967/hf/10163",
                    mask_desc="FULL SKY, no mask (this map IS the foreground)",
                    ordering_determined="RING, by the plane/pole contrast test recorded in the control")
    elif key == "wmap_kband":
        f = f"{DATA}/k-band-hdu-1.parquet"
        m, npix_native = _healpix_from_parquet(f, "TEMPERATURE", "NESTED", nside)
        msk = np.ones(npix)
        meta = dict(file=f, native_nside=512, units="mK",
                    source="https://huggingface.co/datasets/astro-legacy-archive/wmap-band-maps-9yr "
                           "(k-band-hdu-1.parquet); WMAP 9-yr K band, 23 GHz",
                    mask_desc="FULL SKY, no mask (foreground-dominated by construction)",
                    ordering_determined="NESTED, by the plane/pole contrast test recorded in the control")
    elif key == "cobe_dmr":
        import pandas as pd
        f = f"{DATA}/DMR_SKYMAP_GALACTIC_53A_4YR.parquet"
        d = pd.read_parquet(f)
        pix = hp.ang2pix(nside, np.asarray(d["RA"], float), np.asarray(d["DEC"], float), lonlat=True)
        acc = np.zeros(npix); cnt = np.zeros(npix)
        np.add.at(acc, pix, np.asarray(d["SIGNAL"], float)); np.add.at(cnt, pix, 1.0)
        m = np.where(cnt > 0, acc / np.maximum(cnt, 1), 0.0)
        msk = (cnt > 0).astype(float)
        meta = dict(file=f, native_pixels=int(d.shape[0]), units="mK (thermodynamic)",
                    source="https://huggingface.co/datasets/astro-legacy-archive/"
                           "cobe-dmr-four-year-sky-maps (DMR_SKYMAP_GALACTIC_53A_4YR.parquet); "
                           "COBE-DMR 4-year 53 GHz channel A, galactic-frame sky map",
                    mask_desc="pixels with no DMR sample are masked out",
                    repixelisation=("the DMR quad-cube 6144 pixels were binned into HEALPix by their "
                                    "tabulated RA/Dec pixel centres and averaged; this is a "
                                    "repixelisation, not a reprojection, and it is lossy"))
    else:
        raise KeyError(key)
    m = np.asarray(m, float)
    m[~np.isfinite(m)] = 0.0
    return m, np.asarray(msk, float), meta


# ------------------------------------------------------------------- pipeline
def curves(temp, unmasked, edges, tris):
    sys.path.insert(0, ASTRO)
    from lib.cmb_topology import betti_curves_from_topology
    return betti_curves_from_topology(temp, unmasked, edges, tris, NU, sublevel=True)


def stats_from(c):
    i1 = int(np.argmin(np.abs(NU - 1.0)))
    i0 = int(np.argmin(np.abs(NU - 0.0)))
    dnu = float(NU[1] - NU[0])
    return {"b0_at_nu_1": float(c["b0"][i1]), "b1_at_nu_0": float(c["b1"][i0]),
            "b0_curve_integral": float(c["b0"].sum() * dnu),
            "b1_curve_integral": float(c["b1"].sum() * dnu),
            "euler_char_at_nu_0": float(c["euler_char_true"][i0])}


def run_one(key, nside, n_sims, out):
    import healpy as hp
    sys.path.insert(0, ASTRO)
    from lib.cmb_topology import build_topology_fixed
    from lib.stats import coarse_stats_fixed

    print(f"\n=== {key} @ nside {nside}, {n_sims} sims ===")
    m, msk, meta = load_map(key, nside)
    f_sky = float(msk.mean())
    print(f"  f_sky {f_sky:.4f}")

    unmasked, edges, tris = build_topology_fixed(msk, nside)
    lmax = 3 * nside - 1
    cl = hp.anafast(m * msk, lmax=lmax) / max(f_sky, 1e-6)     # f_sky-corrected

    with Timer() as t:
        cd = curves(m, unmasked, edges, tris)
    sd = stats_from(cd)
    print(f"  data: {sd}  ({t.sec:.1f}s)")

    rng = np.random.default_rng(SEED)
    sim_b0, sim_b1, sim_stats = [], [], []
    with Timer() as tt:
        for i in range(n_sims):
            np.random.seed(int(rng.integers(1, 2 ** 31 - 1)))
            sm = hp.synfast(cl, nside, lmax=lmax, verbose=False) if "verbose" in \
                hp.synfast.__code__.co_varnames else hp.synfast(cl, nside, lmax=lmax)
            cs = curves(np.asarray(sm, float) * msk, unmasked, edges, tris)
            sim_b0.append(cs["b0"]); sim_b1.append(cs["b1"]); sim_stats.append(stats_from(cs))
            if (i + 1) % 10 == 0:
                print(f"    sim {i + 1}/{n_sims}")
    wall, mem = t.sec + tt.sec, peak_mb()
    sim_b0, sim_b1 = np.array(sim_b0), np.array(sim_b1)

    cs0 = coarse_stats_fixed(sim_b0, cd["b0"])
    cs1 = coarse_stats_fixed(sim_b1, cd["b1"])
    p0, p1 = cs0.get("empirical_rank_p"), cs1.get("empirical_rank_p")
    print(f"  coarse rank p: b0 {p0}  b1 {p1}")

    # null calibration: hold out sim 0, rank against the other n-1
    cal = coarse_stats_fixed(sim_b0[1:], sim_b0[0])
    p_cal = cal.get("empirical_rank_p")
    cal_ok = bool(p_cal is not None and p_cal > 0.05)
    print(f"  null calibration (held-out sim vs {n_sims - 1} others): rank p = {p_cal} "
          f"{'PASS' if cal_ok else 'FAIL'}")

    null_desc = (f"{n_sims} isotropic Gaussian realisations from this map's OWN pseudo-C_ell, "
                 f"measured with anafast on the masked map to lmax={lmax} and DIVIDED BY f_sky="
                 f"{f_sky:.4f} to undo the mask's power suppression, then synfast at nside={nside} "
                 "and multiplied by the same mask. Spectrum-matched; it does NOT model non-Gaussian "
                 "foregrounds, anisotropic noise or the instrument beam.")
    ds_id = f"astro/{key}"
    kw = dict(source=meta["source"], provenance="observation",
              n_objects=int(unmasked.size), ambient_dim=2, units=meta["units"],
              sha256=sha256_file(meta["file"]), local_path=meta["file"],
              notes=" | ".join(f"{k}: {v}" for k, v in meta.items() if k not in ("file", "source")))
    with topodb() as db:
        db.add_dataset(id=ds_id, domain="astro", title=f"{key} sky map", **kw)
        rid = add_run(db, dataset_id=ds_id, method="lower_star_graph", coeff_field=2, max_dim=2,
                      params=dict(nside=nside, n_sims=n_sims, lmax=lmax, f_sky=f_sky,
                                  nu_grid="41 points on [-3,3]", filtration="lower-star sublevel",
                                  normalisation="field divided by sigma of the UNMASKED pixels",
                                  complex="lib/cmb_topology.build_topology_fixed @ d8175f1",
                                  n_vertices=int(unmasked.size), n_edges=int(len(edges)),
                                  n_triangles=int(len(tris)),
                                  euler_of_complex=int(unmasked.size - len(edges) + len(tris)),
                                  f_sky_correction="pseudo-C_ell divided by f_sky before synfast"),
                      script=SCRIPT, command=f"prlimit --as=8589934592 -- {PY} {SCRIPT}",
                      tier="X", preprocessing=f"ud_grade to nside {nside} (RING); {meta['mask_desc']}",
                      seed=str(SEED), wall_sec=wall, peak_mb=mem)
        db.add_betti(rid, {0: int(sd["b0_at_nu_1"]), 1: int(sd["b1_at_nu_0"])})
        pvals = {}
        for name, tail in STATS_TAIL.items():
            nv = [s[name] for s in sim_stats]
            p, _ = rank_p(sd[name], nv, tail)
            pvals[name] = p
            db.add_statistic(rid, name, sd[name], null_model=null_desc, n_null=n_sims,
                             p_value=p, p_method=p_method_str(n_sims, tail),
                             multiplicity=f"Bonferroni/{NSTAT} -> threshold {0.05 / NSTAT:.4g}")
        for nm, cv, pv in (("coarse_rank_p_b0", cs0, p0), ("coarse_rank_p_b1", cs1, p1)):
            if pv is not None:
                db.add_statistic(rid, nm, float(cv["data_chi2_hartlap"]), null_model=null_desc,
                                 n_null=n_sims, p_value=float(pv),
                                 p_method=("pooled-exchangeable rank p over 8 coarse bins "
                                           f"(lib/stats.coarse_stats_fixed); kept {cv['n_bins_kept']} "
                                           f"bins, retained rank {cv['retained_rank']}; floor "
                                           f"{1 / (n_sims + 1):.4g}"),
                                 multiplicity=f"Bonferroni/{NSTAT} -> threshold {0.05 / NSTAT:.4g}")
                pvals[nm] = float(pv)
        db.add_control(rid, "null_calibration",
                       "a held-out null realisation ranked against the remaining n-1 is not extreme "
                       "(guards the f_sky correction: a mis-normalised null makes the data anomalous "
                       "for free)",
                       passed=cal_ok, detail=f"held-out sim b0 coarse rank p = {p_cal}")
        db.add_control(rid, "known_answer",
                       "Step 0 K5: the same fixed complex has Betti (1,0,1) full sky at nside 32 (run 11)",
                       passed=True, detail=f"this run's complex: V-E+F = "
                                           f"{unmasked.size - len(edges) + len(tris)}")
        if key in ("haslam408", "wmap_kband"):
            fired = min(pvals.values()) <= 2.0 / (n_sims + 1)
            db.add_control(rid, "injection",
                           "POSITIVE CONTROL: a strongly non-Gaussian foreground map must reject the "
                           "Gaussian null at (or near) the rank floor; if it does not, every CMB null "
                           "in this campaign is vacuous",
                           passed=bool(fired),
                           detail=f"min rank p over the declared family = {min(pvals.values()):.4g}, "
                                  f"floor {1 / (n_sims + 1):.4g}")
    out[f"{key}_nside{nside}"] = {"run_id": rid, "f_sky": f_sky, "data": sd, "p_values": pvals,
                                  "null_calibration_p": p_cal, "n_sims": n_sims,
                                  "coarse_b0": {k: v for k, v in cs0.items() if k != "kept_bin_indices"},
                                  "wall_sec": wall}
    print(f"  run {rid}: p = { {k: round(v, 4) for k, v in pvals.items()} }")
    return rid, pvals


def main():
    nside = int(os.environ.get("NSIDE", 128))
    n_sims = int(os.environ.get("NSIMS", 50))
    keys = sys.argv[1:] or ["planck_smica", "wmap_ilc", "haslam408", "wmap_kband", "cobe_dmr"]
    out = {"seed": SEED, "nside": nside, "n_sims": n_sims}
    for k in keys:
        ns = min(nside, 32) if k == "cobe_dmr" else nside   # DMR has only 6144 native pixels
        try:
            run_one(k, ns, n_sims, out)
        except Exception as e:
            print(f"  !! {k} FAILED: {type(e).__name__}: {e}")
            out[k] = {"error": f"{type(e).__name__}: {e}"}
        save_json(f"skymaps_nside{nside}.json", out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
