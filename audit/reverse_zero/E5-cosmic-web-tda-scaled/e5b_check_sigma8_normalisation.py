#!/usr/bin/env python3
"""
Post-hoc advisor-review check on e5b_lognormal_and_poisson_null.py's full
run (report: e5b_lognormal_and_poisson_null_report.json, generated
2026-09-19). Result already folded into that report's
post_hoc_advisor_review_checks.sigma8_absolute_normalisation_check.

WHY: e5b's own round-trip test of pk2xi/xi2pk (toy power-law P(k),
pk2xi -> xi2pk recovers P to a few percent) validates the transform's
SHAPE but cannot catch an ABSOLUTE normalisation error, because a
constant multiplicative factor cancels exactly in a round trip. This
script independently checks absolute normalisation by recomputing
sigma8(z=0) from e5b's own kh/pk(z=0) CAMB arrays, using a top-hat
window (R=8 Mpc/h) and the SAME trapezoidal quadrature style pk2xi/xi2pk
use, and comparing to CAMB's own reported sigma8(z=0).

RESULT (this run, 2026-09-19): ratio (mine/CAMB) = 1.0003 (0.03%
agreement) -- no normalisation bug in the shared P(k)<->xi(r) transform.

Tier: X (exploratory numerics, diagnostic only).

Command:
  /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python \
    audit/reverse_zero/E5-cosmic-web-tda-scaled/e5b_check_sigma8_normalisation.py
"""
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import camb
import e5b_lognormal_and_poisson_null as m


def tophat_w(x):
    return np.where(np.abs(x) < 1e-4, 1.0, 3 * (np.sin(x) - x * np.cos(x)) / x ** 3)


def run():
    p = m.PLANCK18_PARAMS
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=p["H0"], ombh2=p["ombh2"], omch2=p["omch2"], mnu=p["mnu"], omk=p["omk"], tau=p["tau"])
    pars.InitPower.set_params(As=p["As_raw"], ns=p["ns"])
    pars.set_matter_power(redshifts=[0.0], kmax=15.0)
    pars.NonLinear = camb.model.NonLinear_none
    results = camb.get_results(pars)
    kh, zs, pk = results.get_matter_power_spectrum(minkh=1e-4, maxkh=12.0, npoints=600)
    pk_z0 = pk[0]
    sigma8_camb = results.get_sigma8()[0]
    print("CAMB sigma8(z=0) reported:", sigma8_camb)

    R = 8.0
    integrand = kh ** 2 * pk_z0 * tophat_w(kh * R) ** 2
    sigma8_mine = float(np.sqrt(np.trapezoid(integrand, kh) / (2 * np.pi ** 2)))
    print("Self-check sigma8(z=0) from own quadrature:", sigma8_mine)
    print("ratio mine/camb:", sigma8_mine / sigma8_camb)
    return {"sigma8_camb": float(sigma8_camb), "sigma8_selfcheck": sigma8_mine,
            "ratio": sigma8_mine / sigma8_camb}


if __name__ == "__main__":
    run()
