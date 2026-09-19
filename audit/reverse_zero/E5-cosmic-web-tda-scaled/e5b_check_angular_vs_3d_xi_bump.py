#!/usr/bin/env python3
"""
Post-hoc advisor-review check on e5b_lognormal_and_poisson_null.py's full
run (report: e5b_lognormal_and_poisson_null_report.json, generated
2026-09-19). Result already folded into that report's
post_hoc_advisor_review_checks.angular_vs_3d_origin_of_the_xi_bump_at_r_37_45.

WHY: the main run's data xi(r) is non-monotonic (declines from 0.070 at
r=21 Mpc/h to ~0.012 at r=33, rises back to ~0.019 at r=41, then
declines again). Smooth LCDM does not produce this shape (BAO sits at
~105 Mpc/h, not 40). At this shell's comoving distance (~200 Mpc/h), 10
degrees on the sky subtends ~35 Mpc/h transverse -- exactly this range --
so an uncorrected ANGULAR selection/completeness pattern on ~10-degree
scales could produce this bump. The main run's HEALPix mask (nside=64)
is a binary occupied/not-occupied footprint with no completeness
weighting (none is available: the catalogue has only ra, dec, z, zErr).

TEST: recompute the DATA's own Landy-Szalay xi(r) -- same real galaxies,
same r bins -- but with the RANDOM catalogue replaced by
data/real2/cosmic_web/sdss_dr17_random_shuffled_z_ra140_220_dec0_50.csv
(already on disk from an earlier round), restricted to the same shell
window. That random carries the data's EXACT (ra,dec) pairs, so any
PURELY ANGULAR selection/completeness effect divides out of the LS
estimator by construction -- isolating whether the bump is angular or a
genuine 3D feature. This is a SINGLE random catalogue (~2.3x oversample,
vs the main run's 4x/20-seed treatment) -- diagnostic only, tier X.

RESULT (this run, 2026-09-19): the bump PERSISTS (dip near r=33, rise to
~0.010-0.012 at r=39-45) at somewhat lower amplitude than with the main
run's mask-uniform random. CAVEAT on that amplitude difference (added
after a further advisor review): it is EXPECTED and uninformative on its
own -- an angular-position-matched random inherits part of the data's
own real 3D clustering by construction, which generically suppresses the
LS estimator's measured xi(r). Only the PERSISTENCE of the r=37-45 shape
feature, not this amplitude drop, is diagnostic. Persistence under an
angular-position-matched random means the bump is NOT explained by the
nside=64 binary mask -- it is a 3D (or finer-than-nside-64 angular-
completeness-scale) feature the lognormal recipe's single-bias, linear-
theory model does not reproduce. A genuine sub-nside-64 angular effect
(e.g. SDSS fiber-plate boundaries, ~1-2 degree scale) is NOT ruled out by
this test, since it is below this mask's own resolution.

Command:
  timeout 200 prlimit --as=10737418240 -- \
    /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python \
    audit/reverse_zero/E5-cosmic-web-tda-scaled/e5b_check_angular_vs_3d_xi_bump.py
"""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import cosmic_web_tda_scaled as base
import e5b_lognormal_and_poisson_null as m

WT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
SHUF_CSV = os.path.join(WT_ROOT, "data", "real2", "cosmic_web",
                         "sdss_dr17_random_shuffled_z_ra140_220_dec0_50.csv")


def run():
    sel = m.load_real_selection()
    real_xyz = sel["real_xyz"]
    r_lo, r_hi = sel["r_lo"], sel["r_hi"]

    df2 = pd.read_csv(SHUF_CSV, comment="#")
    valid = np.isfinite(df2["z"]) & (df2["z"] > 0)
    df2 = df2[valid].reset_index(drop=True)
    r2_all = base.comoving_r_mpc_over_h(df2["z"].to_numpy(), float(df2["z"].max()))
    mask_win = (r2_all >= r_lo) & (r2_all < r_hi)
    ra2 = df2["ra"].to_numpy()[mask_win]
    dec2 = df2["dec"].to_numpy()[mask_win]
    r2 = r2_all[mask_win]
    print("n available in shuffled-z random (shell window):", len(r2))
    rand_xyz_shuf = base.radec_r_to_xyz(ra2, dec2, r2)

    r_edges_xi = np.linspace(2.0, 70.0, 35)
    r_mid_xi = 0.5 * (r_edges_xi[:-1] + r_edges_xi[1:])
    t0 = time.time()
    rtree_shuf = m.precompute_random_tree(rand_xyz_shuf, r_edges_xi, verbose=True)
    xi_shuf, DD2, RR2, DR2 = m.ls_xi_with_tR(real_xyz, rtree_shuf, verbose=True)
    print("time", time.time() - t0)
    fit_mask = (r_mid_xi >= 20) & (r_mid_xi <= 60)
    print("r    xi_with_angular_matched_random")
    for rr, x2 in zip(r_mid_xi[fit_mask], xi_shuf[fit_mask]):
        print(f"{rr:5.1f}  {x2:9.5f}")
    return {"r_mid": r_mid_xi[fit_mask].tolist(), "xi_angular_matched_random": xi_shuf[fit_mask].tolist()}


if __name__ == "__main__":
    run()
