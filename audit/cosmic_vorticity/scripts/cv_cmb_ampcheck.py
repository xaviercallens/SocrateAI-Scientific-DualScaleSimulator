#!/usr/bin/env python3
"""Make the CMB injection amplitude interpretable.

The injection grid adds spots of peak amplitude A * sigma_T where sigma_T is
the standard deviation of the PREPPED (unsmoothed) map, but the detector
thresholds on sigma(T_s) of the SMOOTHED map, and a theta_0 = sigma_s Gaussian
spot loses part of its peak to that smoothing.  Without the two conversion
numbers, "smallest detectable A/sigma_T" cannot be read as a distance above the
detection threshold.  This script MEASURES both, on simulations plus the map's
own sigma, and also converts the injected full-sky counts into the effective
in-mask density that the detector actually sees.

Command: timeout 900 prlimit --as=8589934592 -- <venv-tda python> cv_cmb_ampcheck.py
"""
import json
import os
import sys

import numpy as np
import healpy as hp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cv_lib as L  # noqa: E402
import cv_cmb as C  # noqa: E402


def main():
    m, mask = L.load_data("wmap")
    tab, _ = L.disk_table()
    er = L.eroded_mask(mask, tab)
    tp = L.prep_map(m, mask)
    sigma_T = float(tp[mask > 0].std())
    ts = L.smooth_map(tp, mask)
    sigma_Ts = float(ts[mask > 0].std())

    # peak of ONE injected spot after the detector's own smoothing, in sigma(T_s)
    npix = hp.nside2npix(L.NSIDE)
    zero = np.zeros(npix)
    centre = np.array(hp.pix2vec(L.NSIDE, npix // 2))
    rng = np.random.default_rng(0)
    one = C.inject(zero, centre[None, :], 1.0 * sigma_T, rng=rng)   # A = 1 sigma_T
    one_s = L.smooth_map(one, np.ones(npix, np.uint8))
    peak_raw = float(np.abs(one).max())
    peak_smoothed = float(np.abs(one_s).max())

    fsky_mask = float(mask.mean())
    fsky_er = float(er.mean())
    sky_deg2 = 4 * np.pi * (180 / np.pi) ** 2

    out = {
        "purpose": "conversion factors that make the CMB injection amplitudes and densities readable",
        "map": "wmap",
        "sigma_T_prepped_unsmoothed": sigma_T,
        "sigma_Ts_smoothed": sigma_Ts,
        "ratio_sigma_Ts_over_sigma_T": sigma_Ts / sigma_T,
        "one_spot_at_A_equals_1_sigma_T": {
            "theta0_deg": L.SIGMA_S_DEG,
            "peak_before_smoothing_over_sigma_T": peak_raw / sigma_T,
            "peak_after_smoothing_over_sigma_T": peak_smoothed / sigma_T,
            "peak_after_smoothing_over_sigma_Ts": peak_smoothed / sigma_Ts,
            "detector_threshold_in_sigma_Ts": L.NU_PRIMARY,
            "reading": "an A = 1 sigma_T spot reaches %.2f sigma(T_s) after the detector's own "
                       "smoothing, against a threshold of %.1f sigma(T_s); the smoothing loss is the "
                       "factor %.3f." % (peak_smoothed / sigma_Ts, L.NU_PRIMARY, peak_smoothed / peak_raw),
        },
        "footprint": {
            "fsky_mask": fsky_mask, "fsky_after_R_disk_erosion": fsky_er,
            "full_sky_deg2": sky_deg2, "usable_deg2": fsky_er * sky_deg2,
            "note": "registered injection densities are FULL-SKY counts; only the eroded fraction is "
                    "usable, so the effective in-mask count is N_inj * fsky_eroded.",
        },
        "registered_densities_converted": [
            {"N_inj_full_sky": n, "effective_in_mask_count": n * fsky_er,
             "density_per_deg2_full_sky": n / sky_deg2,
             "mean_separation_deg_full_sky": float(np.degrees(np.sqrt(2 / (np.sqrt(3) * (n / (4 * np.pi))))))}
            for n in (100, 300, 1000)
        ],
        "detected_population_for_comparison": {
            "N_data_wmap": int(json.load(open(os.path.join(L.CV, "results", "cmb_data_wmap.json")))
                               ["primary"]["S2_count"]),
            "note": "the injected populations are to be compared with the ~350 candidates the detector "
                    "already finds in the WMAP footprint from the Gaussian field alone; an injected "
                    "population must compete with that background, it is not added to an empty sky.",
        },
        "tier": "X",
    }
    json.dump(out, open(os.path.join(L.CV, "results", "cmb_amplitude_units.json"), "w"), indent=1)
    print(json.dumps({k: out[k] for k in ("sigma_T_prepped_unsmoothed", "sigma_Ts_smoothed",
                                          "ratio_sigma_Ts_over_sigma_T")}, indent=1))
    print(out["one_spot_at_A_equals_1_sigma_T"]["reading"])
    print(json.dumps(out["registered_densities_converted"], indent=1))


if __name__ == "__main__":
    main()
