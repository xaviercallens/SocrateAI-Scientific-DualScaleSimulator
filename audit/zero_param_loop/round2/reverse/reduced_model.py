#!/usr/bin/env python3
"""
REVERSE LOOP (round 2), step 1: build the REDUCED model.

Applies BOTH round-2 forward proposals from ../reduce_proposals.json:

  1. a_pot, b_pot DELETED (mechanism: delete_unobservable). The entire
     quintessence sector (compute_dark_energy_observables) is replaced by
     an ANALYTIC flat LCDM block with Omega_m = 1 - omegaLambda = 0.31115
     (omegaLambda := 0.68885, DarkEnergyScale.lean:84) and
     H0 = c / hubbleRadius_m = 67.661122... km/s/Mpc
     (hubbleRadius_m := 1.3672e26 m, SelfDualCutoff.lean:100). This block
     has ZERO free parameters -- it does not call workshopcosmo's
     quintessence ODE at all.

  2. lambda_sym DELETED (mechanism: delete_unobservable). The forward
     round's exact-algebra result (lambda_sym cancels exactly once the
     symmetron ODE is rewritten in psi=phi/phi_0 form; see
     ../lambda_sym_negative_control.py) is APPLIED here by fixing
     lambda_sym at its old nominal value (1.0) and never sweeping it --
     this is NOT re-derived in this script; it is taken from the forward
     round's evidence, exactly as the ground rules for a reverse loop
     require ("apply the proposed reductions ... via harness params").
     The actual psi-form rewrite of run_symmetron_screening_simulation in
     workshopcosmo.py is still NOT implemented (out of scope: this
     worktree may not touch proofs/ or do a full code refactor beyond the
     audit/ and scripts/ areas already used by this loop); holding
     lambda_sym fixed is behaviorally equivalent for every number this
     harness reads out, per the forward round's own negative control.

REMAINING FREE PARAMETERS (2): mu_sym, c4_pta_product.

c4_pta_product is split back into (pta_suppression, c4_c0_ratio) =
(0.005, c4_pta_product/0.005) before calling the untouched
param_loop_sim.compute_pta_observable -- any positive split gives
byte-identical results (round-1 Check C, max diff 1.110e-16); the fixed
denominator 0.005 matches param_loop_sim.DEFAULT_PARAMS["pta_suppression"].

IMPORTANT (load-bearing, per advisor review): mu_sym only feeds the
screening block, and c4_pta_product only feeds the PTA block; NEITHER
appears anywhere in the dark-energy block, which is now a constant. No
real fifth-force/short-range-gravity dataset and no real PTA angular
Gamma(theta) dataset exist in this repo's data/real/ (fifth_force_screening
and pta angular-correlation are both ABSENT, data/real/MANIFEST.json).
So this "2-parameter" reduced model has ZERO parameters that any fetched
dataset in this loop can actually see: chi2_reduced below is a single
computed number, not a fit. Calling this model "a 2-parameter fit to the
data" would be false; it fits the DESI+Pantheon+ data with 0 parameters.
"""
import math
import os
import sys
from typing import Any, Dict

import numpy as np

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
sys.path.insert(0, REPO_ROOT)
sys.path.insert(0, os.path.join(REPO_ROOT, "scripts"))
import param_loop_sim as pls  # noqa: E402

REDUCED_PARAM_NAMES = ["mu_sym", "c4_pta_product"]

# --- Lean constants (quoted verbatim, same pins used by
#     ../decisive_experiment.py: omegaLambda DarkEnergyScale.lean:84,
#     hubbleRadius_m SelfDualCutoff.lean:100) ---
LEAN_OMEGA_LAMBDA = 0.68885
LEAN_HUBBLE_RADIUS_M = 1.3672e26
C_LIGHT_M_S = 299792458.0
MPC_IN_M = 3.0856775814913673e22
LEAN_H0_KM_S_MPC = (C_LIGHT_M_S / LEAN_HUBBLE_RADIUS_M) * MPC_IN_M / 1000.0
OMEGA_M_FROZEN = 1.0 - LEAN_OMEGA_LAMBDA

LAMBDA_SYM_FIXED = pls.DEFAULT_PARAMS["lambda_sym"]  # 1.0, held fixed (deleted, not swept)
PTA_SUPPRESSION_FIXED = pls.DEFAULT_PARAMS["pta_suppression"]  # 0.005, arbitrary split denominator

REDUCED_DEFAULT_PARAMS: Dict[str, float] = {
    "mu_sym": pls.DEFAULT_PARAMS["mu_sym"],
    "c4_pta_product": pls.DEFAULT_PARAMS["c4_c0_ratio"] * pls.DEFAULT_PARAMS["pta_suppression"],  # 0.08035
}

Z_GRID = pls.Z_GRID
Z_POINTS = pls.Z_POINTS


def frozen_lcdm_dark_energy() -> Dict[str, Any]:
    """Analytic flat LCDM(Om=0.31115) block. Zero free parameters. Every
    number below is a deterministic function of Z_GRID/Z_POINTS and the two
    Lean constants; no ODE solve, no a_pot/b_pot."""
    def Ez(z):
        return np.sqrt(OMEGA_M_FROZEN * (1.0 + z) ** 3 + (1.0 - OMEGA_M_FROZEN))

    zg_fine = np.linspace(0.0, 2.45, 4000)
    Ez_fine = Ez(zg_fine)
    inv_E = 1.0 / Ez_fine
    DC_fine = np.concatenate(([0.0], np.cumsum(0.5 * (inv_E[1:] + inv_E[:-1]) * np.diff(zg_fine))))

    z_grid = np.asarray(Z_GRID, dtype=float)
    z_points = np.asarray(Z_POINTS, dtype=float)
    Hz_H0_grid = Ez(z_grid)
    Hz_H0_points = Ez(z_points)
    DM_H0_grid = np.interp(z_grid, zg_fine, DC_fine)
    DM_H0_points = np.interp(z_points, zg_fine, DC_fine)
    DL_H0_grid = (1.0 + z_grid) * DM_H0_grid
    DL_H0_points = (1.0 + z_points) * DM_H0_points
    mu_shape_grid = 5.0 * np.log10(np.clip(DL_H0_grid, 1e-12, None))
    mu_shape_points = 5.0 * np.log10(np.clip(DL_H0_points, 1e-12, None))

    return {
        "success": True,
        "sector": "frozen_flat_LCDM_zero_param",
        "Omega_m_frozen": OMEGA_M_FROZEN,
        "Omega_Lambda_frozen": LEAN_OMEGA_LAMBDA,
        "H0_implied_km_s_Mpc": LEAN_H0_KM_S_MPC,
        "w0_cpl_latetime": -1.0,
        "wa_cpl_latetime": 0.0,
        "z_grid": z_grid.tolist(),
        "w_of_z_grid": [-1.0] * len(z_grid),
        "z_points": z_points.tolist(),
        "w_of_z_points": [-1.0] * len(z_points),
        "H_of_z_over_H0_grid": Hz_H0_grid.tolist(),
        "H_of_z_over_H0_points": Hz_H0_points.tolist(),
        "D_M_times_H0_grid": DM_H0_grid.tolist(),
        "D_M_times_H0_points": DM_H0_points.tolist(),
        "distance_modulus_shape_grid": mu_shape_grid.tolist(),
        "distance_modulus_shape_points": mu_shape_points.tolist(),
    }


# Cached: this block is a pure constant, computed once.
_FROZEN_DE = frozen_lcdm_dark_energy()


def evaluate_reduced_point(params: Dict[str, float]) -> Dict[str, Any]:
    merged = dict(REDUCED_DEFAULT_PARAMS)
    for k in params:
        if k not in REDUCED_DEFAULT_PARAMS:
            raise SystemExit(f"Unknown reduced parameter '{k}'. Valid keys: {sorted(REDUCED_DEFAULT_PARAMS)}")
    merged.update(params)

    mu_sym = float(merged["mu_sym"])
    product = float(merged["c4_pta_product"])
    c4_c0_ratio = product / PTA_SUPPRESSION_FIXED if product > 0 else float("nan")

    screening = pls.compute_screening_observable(mu_sym, LAMBDA_SYM_FIXED)
    pta = pls.compute_pta_observable(PTA_SUPPRESSION_FIXED, c4_c0_ratio)

    return {
        "reduced_params": merged,
        "lambda_sym_fixed_not_free": LAMBDA_SYM_FIXED,
        "dark_energy": _FROZEN_DE,
        "screening": screening,
        "pta": pta,
    }


def _selftest() -> bool:
    """(a) c4_pta_product split-invariance (reuses round-1's proof: any
    positive split gives byte-identical gamma_theta). (b) dark-energy block
    is byte-identical across two different (mu_sym, c4_pta_product)
    points -- i.e. genuinely zero-variance, the load-bearing fact for the
    TDA prediction below."""
    p1 = evaluate_reduced_point({"mu_sym": 1.0, "c4_pta_product": 0.08035})
    p2 = evaluate_reduced_point({"mu_sym": 2.5, "c4_pta_product": 0.5})
    de1, de2 = p1["dark_energy"], p2["dark_energy"]
    de_identical = (de1["H_of_z_over_H0_grid"] == de2["H_of_z_over_H0_grid"]) and \
                   (de1["D_M_times_H0_grid"] == de2["D_M_times_H0_grid"])

    alt = pls.compute_pta_observable(0.001, 0.08035 / 0.001)
    nominal = pls.compute_pta_observable(PTA_SUPPRESSION_FIXED, 0.08035 / PTA_SUPPRESSION_FIXED)
    diff_split = float(np.max(np.abs(np.array(alt["gamma_theta"]) - np.array(nominal["gamma_theta"]))))

    ok = de_identical and (diff_split < 1e-12)
    print(f"[reduced_model round2 selftest] dark_energy block identical across 2 different (mu_sym,c4_pta_product) points: {de_identical}")
    print(f"[reduced_model round2 selftest] diff_split(c4_pta_product invariance)={diff_split:.3e}")
    print(f"[reduced_model round2 selftest] H0_implied_km_s_Mpc={LEAN_H0_KM_S_MPC}, Omega_m_frozen={OMEGA_M_FROZEN}")
    print(f"[reduced_model round2 selftest] ALL_PASS={ok}")
    return ok


if __name__ == "__main__":
    ok = _selftest()
    raise SystemExit(0 if ok else 1)
