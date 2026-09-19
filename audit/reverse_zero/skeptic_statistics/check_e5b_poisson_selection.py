#!/usr/bin/env python3
"""Skeptic check: does E5b's Poisson-in-footprint null share the DATA's radial selection and angular
footprint? The data go through base.select_flat_density_subsample (N_TARGET per equal-volume shell);
the Poisson null draws r uniform-in-volume in [r_lo, r_hi] directly at N_FINAL (no per-shell selection).
Checks: N, per-shell counts, two-sample KS of r (data vs Poisson seed 3000, the first committed seed),
angular KS of ra and sin(dec), fraction of Poisson points in HEALPix nside=512 pixels that contain no
data galaxy at all (sub-pixel footprint holes the nside=64 mask cannot see), and the comoving size of
an nside=64 pixel at r_lo/r_hi. Seeds: Poisson seed 3000 (same as committed run). Command (worktree root):
  /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python audit/reverse_zero/skeptic_statistics/check_e5b_poisson_selection.py
"""
import json, os, sys
import numpy as np, healpy as hp
from scipy.stats import ks_2samp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "E5-cosmic-web-tda-scaled"))
import e5b_lognormal_and_poisson_null as E  # noqa: E402
sel = E.load_real_selection()
occ, frac = E.build_mask(sel["df_ra"], sel["df_dec"], E.NSIDE_MASK)
ra_min, ra_max = float(sel["df_ra"].min()), float(sel["df_ra"].max())
dec_min, dec_max = float(sel["df_dec"].min()), float(sel["df_dec"].max())
pra, pdec, pr = E.poisson_footprint_catalogue(occ, E.NSIDE_MASK, ra_min, ra_max, dec_min, dec_max, sel["r_lo"], sel["r_hi"], E.N_FINAL, 3000)
edges = sel["edges"]
cnt = lambda r: [int(((r >= edges[i]) & (r < edges[i + 1])).sum()) for i in E.SELECTION_SHELL_IDX]
# fine-scale footprint holes: nside=512 pixels (~0.11 deg) occupied by ANY galaxy of the full catalogue
nfine = 512
occ_f = np.zeros(hp.nside2npix(nfine), bool)
occ_f[hp.ang2pix(nfine, np.radians(90 - sel["df_dec"]), np.radians(sel["df_ra"]))] = True
pin = occ_f[hp.ang2pix(nfine, np.radians(90 - pdec), np.radians(pra))]
din = occ_f[hp.ang2pix(nfine, np.radians(90 - sel["real_dec"]), np.radians(sel["real_ra"]))]
pix = hp.nside2resol(E.NSIDE_MASK)
out = {"N_data": int(len(sel["real_r"])), "N_poisson": int(len(pr)),
       "shell_counts_data": cnt(sel["real_r"]), "shell_counts_poisson": cnt(pr),
       "ks_r": {"stat": float(ks_2samp(sel["real_r"], pr).statistic), "p": float(ks_2samp(sel["real_r"], pr).pvalue)},
       "ks_ra": {"stat": float(ks_2samp(sel["real_ra"], pra).statistic), "p": float(ks_2samp(sel["real_ra"], pra).pvalue)},
       "ks_sindec": {"stat": float(ks_2samp(np.sin(np.radians(sel["real_dec"])), np.sin(np.radians(pdec))).statistic),
                     "p": float(ks_2samp(np.sin(np.radians(sel["real_dec"])), np.sin(np.radians(pdec))).pvalue)},
       "frac_poisson_in_nside512_pixels_with_no_galaxy": float(1 - pin.mean()),
       "frac_data_in_such_pixels": float(1 - din.mean()),
       "nside64_pixel_deg": float(np.degrees(pix)),
       "nside64_pixel_comoving_mpc_h_at_r_lo_r_hi": [float(pix * sel["r_lo"]), float(pix * sel["r_hi"])],
       "r_lo_r_hi": [sel["r_lo"], sel["r_hi"]],
       "note": ("frac_poisson_in_nside512_pixels_with_no_galaxy is NOT a hole diagnostic: ~193k galaxies over the"
                " footprint give <1 galaxy per nside=512 pixel on average, so emptiness is sparsity-dominated;"
                " frac_data=0 is by construction. ks_ra/ks_sindec differences include genuine angular LSS"
                " (data are clustered), so they are not by themselves a footprint mismatch.")}
json.dump(out, open(os.path.join(HERE, "check_e5b_poisson_selection.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
