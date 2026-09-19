#!/usr/bin/env python3
"""
E3 -- M0 (zero-parameter hypothesis) vs M2 (mu_sym, c4 free) vs fitted LCDM.

M0 is defined by GROUND RULE (task text, quoted): "flat LCDM (Omega_L = 0.68885
frozen, H0 a nuisance as before), GR tensor sector (c4 = 0), NO symmetron sector
(mu_sym removed with the sector)". Per the FRAMING RULE, M0 is a HYPOTHESIS
CHANGE (sectors removed / set to GR values), NOT a derivation from K3 x T2.

Every number in M0 and its source (Lean def / literature / nuisance) is listed
in the "parameters" block of the JSON this script writes; there is no hidden
tuning. M0 has ZERO free theory parameters:
  - Omega_Lambda = 0.68885           -- LeanMaster DarkEnergyScale.lean:84 (Planck 2018, tier L
                                         as physics; tier C to identify with this fit, per
                                         REPORT.md tier note). FROZEN, not fitted.
  - c4_pta_product = 0                -- GR value (Hellings-Downs exactly); mu_sym and
    mu_sym: SECTOR REMOVED              lambda_sym: the symmetron sector is deleted entirely
                                         (not "set to a value"), per the ground-rule wording.
  - H0 (equivalently the BAO amplitude scale u_star, and the SN absolute
    magnitude offset): NUISANCE, marginalized analytically in closed form
    (same convention as round2/round3: chi2_bao_given_base's u_star, and
    chi2_sn_given_mu's additive offset). This is NOT a theory parameter --
    it never multiplies a THEORY prediction that varies with z-shape; it only
    rescales/shifts a shape that is otherwise 100% fixed by Omega_Lambda.

M2 = M0 with mu_sym and c4_pta_product reinstated as free (fitted) parameters.
Because no fifth-force/screening dataset and no PTA angular-correlation dataset
enter this chi2 (both ABSENT, see data manifest), M2's best-fit chi2 on THIS data
is numerically identical to M0's: reinstating the sectors changes nothing that
the DE data can see. This is stated, not glossed over (same conclusion as
round-3 chi2_independent.json "probes_note").

Data:
  DESI DR2 BAO ALL_GCcomb consensus distances + 13x13 covariance
    (data/real2/dark_energy/desi_dr2/desi_gaussian_bao_ALL_GCcomb_{mean,cov}.txt,
     fetched_verified, MANIFEST.json in that directory).
  Pantheon+SH0ES full STAT+SYS covariance (1701x1701), sliced with a boolean
    mask (zHD>0.01, zHD<=2.4, IS_CALIBRATOR==0) to the SAME 1580-SN
    cosmology-only cut convention as rounds 2-3's diagonal-error analysis
    (which used 1590; the 10-SN difference is the added IS_CALIBRATOR==0
    filter, since Cepheid-host calibrator SNe are anchored to CEPH_DIST, not
    a Hubble-flow distance modulus, and must not be fit with the same
    shape+offset cosmological model). The row/column mask is applied
    identically to the data vector and the covariance so correspondence is
    exact. This is the "full Pantheon+ covariance" pipeline action named in
    PRE_REGISTRATION.md P1 (an upgrade from round2/3's diagonal errors). The
    full 1701-row variant INCLUDING calibrators is also computed and reported
    as a disclosed (non-primary) comparison.

Command:
  /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python \
    audit/reverse_zero/E3-M0/m0_model.py
  (run from /mnt/disks/disk-socrateai-local-1/dualscale-wt-reverse)

Seeds: no RNG is used (chi2 is a deterministic closed-form/least-squares
computation); no seed needed, consistent with round2/round3 scripts.
"""
import json
import os

import numpy as np
import pandas as pd

WT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
HERE = os.path.dirname(os.path.abspath(__file__))

# --- frozen constants (Lean defs, quoted verbatim) ---------------------------
OMEGA_LAMBDA = 0.68885          # LeanMaster DarkEnergyScale.lean:84 (tier L physics / tier C ID)
OMEGA_M_FROZEN = 1.0 - OMEGA_LAMBDA
C4_PTA_PRODUCT_M0 = 0.0         # GR value: Gamma(theta) = Hellings-Downs exactly (no free tensor sector)

# --- DESI DR2 ALL_GCcomb (fetched_verified; see data/real2/dark_energy/desi_dr2/MANIFEST.json) ---
DR2_DIR = os.path.join(WT_ROOT, "data", "real2", "dark_energy", "desi_dr2")
_mean_rows = []
with open(os.path.join(DR2_DIR, "desi_gaussian_bao_ALL_GCcomb_mean.txt")) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        z, val, q = line.split()
        _mean_rows.append((float(z), float(val), q))
DESI_Z = np.array([r[0] for r in _mean_rows])
DESI_VAL = np.array([r[1] for r in _mean_rows])
DESI_Q = [r[2] for r in _mean_rows]
DESI_COV = np.loadtxt(os.path.join(DR2_DIR, "desi_gaussian_bao_ALL_GCcomb_cov.txt"))
assert DESI_COV.shape == (len(DESI_VAL), len(DESI_VAL)), DESI_COV.shape
DESI_COV_INV = np.linalg.inv(DESI_COV)
N_BAO = len(DESI_VAL)

# --- Pantheon+SH0ES: same 1701-row file as rounds 1-3, matched to the FULL
#     STAT+SYS covariance (upgrade from round-3's diagonal-only errors) ------
SN_PATH = os.path.join(WT_ROOT, "data", "real", "dark_energy", "pantheon_plus_sh0es.dat")
_sn_full = pd.read_csv(SN_PATH, sep=r"\s+", engine="python")
assert len(_sn_full) == 1701, len(_sn_full)

COV_PATH = os.path.join(WT_ROOT, "data", "real2", "dark_energy", "Pantheon+SH0ES_STAT+SYS.cov")
with open(COV_PATH) as f:
    _n_declared = int(f.readline().strip())
assert _n_declared == 1701, _n_declared
SN_COV_FLAT = np.loadtxt(COV_PATH, skiprows=1)  # 1701*1701 values, one per line (Pantheon+ release format)
assert SN_COV_FLAT.shape == (1701 * 1701,), SN_COV_FLAT.shape
SN_COV_FULL = SN_COV_FLAT.reshape(1701, 1701)
assert np.allclose(SN_COV_FULL, SN_COV_FULL.T, atol=1e-6), "covariance not symmetric"

# Cosmology-only cut (advisor-review fix): the raw 1701-row file includes ~77
# Cepheid-host calibrator SNe (IS_CALIBRATOR==1) whose MU_SH0ES is tied to
# CEPH_DIST, a distance-ladder anchor, NOT a cosmological (Hubble-flow)
# distance modulus -- fitting them with the same shape+offset model as
# Hubble-flow SNe is a likelihood error, not a data upgrade. Rounds 2-3 used
# this same cut (zHD>0.01, zHD<=2.4) and landed on 1590 SNe (no explicit
# IS_CALIBRATOR filter there); adding IS_CALIBRATOR==0 removes 10 more,
# landing on 1580. The covariance is sliced with the SAME boolean mask so
# row/column correspondence with the cut data vector is exact.
SN_CUT_MASK = (
    (_sn_full["zHD"] > 0.01)
    & (_sn_full["zHD"] <= 2.4)
    & (_sn_full["IS_CALIBRATOR"] == 0)
    & np.isfinite(_sn_full["MU_SH0ES"])
)
_sn = _sn_full[SN_CUT_MASK]
SN_Z = _sn["zHD"].to_numpy()
SN_MU = _sn["MU_SH0ES"].to_numpy()
assert np.isfinite(SN_MU).all() and np.all(SN_Z > 0)
SN_COV = SN_COV_FULL[np.ix_(SN_CUT_MASK.to_numpy(), SN_CUT_MASK.to_numpy())]
N_SN = int(SN_CUT_MASK.sum())
SN_COV_INV = np.linalg.inv(SN_COV)

# Disclosed variant kept for comparison: the FULL 1701-row sample including
# calibrators, run with the identical pipeline further below (see
# "variant_full_1701_including_calibrators" in the output JSON).
SN_Z_FULL = _sn_full["zHD"].to_numpy()
SN_MU_FULL = _sn_full["MU_SH0ES"].to_numpy()
SN_COV_INV_FULL = np.linalg.inv(SN_COV_FULL)

N_TOTAL = N_BAO + N_SN


def lcdm_Ez(zg, Om):
    return np.sqrt(Om * (1.0 + zg) ** 3 + (1.0 - Om))


def comoving_distance(zg, Ez):
    inv_E = 1.0 / Ez
    dc = np.concatenate(([0.0], np.cumsum(0.5 * (inv_E[1:] + inv_E[:-1]) * np.diff(zg))))
    return dc, inv_E


ZGRID = np.linspace(0.0, 2.5, 25001)


def bao_base(Om):
    Ez = lcdm_Ez(ZGRID, Om)
    dc, inv_E = comoving_distance(ZGRID, Ez)
    xM = {z: float(np.interp(z, ZGRID, dc)) for z in set(DESI_Z)}
    xH = {z: float(np.interp(z, ZGRID, inv_E)) for z in set(DESI_Z)}
    base = []
    for z, q in zip(DESI_Z, DESI_Q):
        if q == "DM_over_rs":
            base.append(xM[z])
        elif q == "DH_over_rs":
            base.append(xH[z])
        elif q == "DV_over_rs":
            base.append((z * xM[z] ** 2 * xH[z]) ** (1.0 / 3.0))
        else:
            raise ValueError(q)
    return np.array(base)


def chi2_bao(Om):
    base = bao_base(Om)
    denom = float(base @ DESI_COV_INV @ base)
    u_star = float(base @ DESI_COV_INV @ DESI_VAL) / denom  # analytic amplitude nuisance (~ 1/(r_d H0))
    resid = DESI_VAL - u_star * base
    return {"chi2": float(resid @ DESI_COV_INV @ resid), "u_star": u_star}


def sn_mu_shape(Om, z):
    Ez = lcdm_Ez(ZGRID, Om)
    dc, _ = comoving_distance(ZGRID, Ez)
    dc_q = np.interp(z, ZGRID, dc)
    dl_q = (1.0 + z) * dc_q
    return 5.0 * np.log10(np.clip(dl_q, 1e-12, None))


def chi2_sn(Om, z=None, mu=None, cov_inv=None):
    z = SN_Z if z is None else z
    mu = SN_MU if mu is None else mu
    cov_inv = SN_COV_INV if cov_inv is None else cov_inv
    m0 = sn_mu_shape(Om, z)
    resid0 = mu - m0
    ones = np.ones_like(resid0)
    # analytic additive-offset nuisance under the (possibly sliced) covariance:
    # offset = (1^T C^-1 r0) / (1^T C^-1 1)
    denom = float(ones @ cov_inv @ ones)
    offset = float(ones @ cov_inv @ resid0) / denom
    resid = resid0 - offset
    return {"chi2": float(resid @ cov_inv @ resid), "offset_fit": offset}


def total_chi2(Om, sn_variant="cut_1580"):
    b = chi2_bao(Om)
    if sn_variant == "cut_1580":
        s = chi2_sn(Om)
    elif sn_variant == "full_1701":
        s = chi2_sn(Om, z=SN_Z_FULL, mu=SN_MU_FULL, cov_inv=SN_COV_INV_FULL)
    else:
        raise ValueError(sn_variant)
    return {"chi2": b["chi2"] + s["chi2"], "bao": b, "sn": s}


def aic_bic(chi2, k, n):
    return {"AIC": chi2 + 2 * k, "BIC": chi2 + k * np.log(n)}


def main():
    out = {}

    out["data"] = {
        "desi_dr2_all_gccomb": {
            "path": "data/real2/dark_energy/desi_dr2/desi_gaussian_bao_ALL_GCcomb_{mean,cov}.txt",
            "n_rows": N_BAO,
            "z_unique": sorted(set(DESI_Z.tolist())),
        },
        "pantheon_plus_full_cov": {
            "sn_path": "data/real/dark_energy/pantheon_plus_sh0es.dat",
            "cov_path": "data/real2/dark_energy/Pantheon+SH0ES_STAT+SYS.cov",
            "n_sn_primary_cosmology_cut": N_SN,
            "cut": "zHD>0.01 & zHD<=2.4 & IS_CALIBRATOR==0, applied identically to the data vector and the 1701x1701 covariance",
            "n_sn_full_1701_disclosed_variant": 1701,
            "note": "PRIMARY uses the full STAT+SYS covariance (upgrade from round2/3's diagonal-only"
                    " errors) sliced to a cosmology-only 1580-SN sample (round2/3 used 1590 with the"
                    " same z-cut but no IS_CALIBRATOR filter, and diagonal errors only). The full"
                    " 1701-row variant (including 77 Cepheid-host calibrators) is reported separately"
                    " as a disclosed, non-primary comparison -- see variant_full_1701_including_calibrators.",
        },
        "n_total": N_TOTAL,
    }

    # M0: zero free theory parameters. List EVERY number and its source.
    out["parameters"] = {
        "M0": {
            "Omega_Lambda": {"value": OMEGA_LAMBDA, "source": "LeanMaster DarkEnergyScale.lean:84 (Planck 2018; kernel-fixed Lean def, tier L as physics)", "status": "FROZEN"},
            "Omega_m": {"value": OMEGA_M_FROZEN, "source": "= 1 - Omega_Lambda (flat FRW)", "status": "FROZEN"},
            "c4_pta_product": {"value": C4_PTA_PRODUCT_M0, "source": "GR value: PTA angular correlation = Hellings-Downs exactly, no l=4 excess", "status": "FROZEN (sector present, value = GR)"},
            "mu_sym / symmetron sector": {"value": None, "source": "sector deleted entirely (ground-rule wording: 'NO symmetron sector ... removed with the sector')", "status": "REMOVED"},
            "H0 (equiv. u_star, sn_offset)": {"value": "fitted analytically per-dataset", "source": "nuisance, marginalized in closed form (same convention as round2/round3 chi2_bao_given_base/chi2_sn_given_mu)", "status": "NUISANCE, not a theory parameter"},
        },
        "M0_free_theory_parameter_count": 0,
        "M0_free_theory_parameter_count_note": "Omega_Lambda is frozen (not fitted); c4_pta_product is frozen at the GR value 0; mu_sym has no value because its sector is removed. H0/amplitude/offset are nuisances by the same convention already used and disclosed in round2/round3 (they marginalize identically regardless of Omega_m, so 'freezing H0' changes nothing in this fit, as flagged in round2/decisive_experiment.py).",
        "M2": {
            "description": "M0 with mu_sym and c4_pta_product reinstated as FREE (fitted) parameters -- but no fifth-force/screening or PTA-angular dataset enters this chi2 (both ABSENT per the data manifest), so their best-fit values are UNCONSTRAINED and the minimized chi2 on this data is numerically identical to M0's.",
            "k_theory": 2,
        },
        "LCDM_fitted": {
            "description": "Omega_m free (fitted), c4_pta_product=0, no symmetron sector -- same sectors as M0, but Omega_m no longer frozen.",
            "k_theory": 1,
        },
    }

    # --- M0: Omega_m frozen ---
    m0 = total_chi2(OMEGA_M_FROZEN)
    out["M0_result"] = m0
    out["M0_result"]["k"] = 0
    out["M0_result"].update(aic_bic(m0["chi2"], 0, N_TOTAL))

    # --- M2: same sectors, mu_sym & c4 free but unconstrained by this data -> identical chi2 ---
    m2 = dict(m0)
    out["M2_result"] = {
        "chi2": m2["chi2"],
        "k": 2,
        "note": "chi2 identical to M0 because mu_sym, c4_pta_product enter no term in this data's likelihood (no fifth-force/screening or PTA-angular dataset present; both ABSENT).",
        **aic_bic(m2["chi2"], 2, N_TOTAL),
    }

    # --- fitted LCDM: Omega_m free, minimize total chi2 over Om ---
    from scipy.optimize import minimize_scalar
    r = minimize_scalar(lambda Om: total_chi2(Om)["chi2"], bounds=(0.05, 0.95), method="bounded",
                         options={"xatol": 1e-10})
    fitted = total_chi2(float(r.x))
    out["LCDM_fitted_result"] = {
        "Omega_m": float(r.x),
        "chi2": fitted["chi2"],
        "bao_chi2": fitted["bao"]["chi2"],
        "sn_chi2": fitted["sn"]["chi2"],
        "k": 1,
        **aic_bic(fitted["chi2"], 1, N_TOTAL),
    }

    out["delta_chi2_M0_minus_LCDM_fitted"] = m0["chi2"] - fitted["chi2"]
    out["delta_AIC_M0_minus_LCDM_fitted"] = out["M0_result"]["AIC"] - out["LCDM_fitted_result"]["AIC"]
    out["delta_BIC_M0_minus_LCDM_fitted"] = out["M0_result"]["BIC"] - out["LCDM_fitted_result"]["BIC"]
    out["delta_AIC_M2_minus_M0"] = out["M2_result"]["AIC"] - out["M0_result"]["AIC"]
    out["delta_BIC_M2_minus_M0"] = out["M2_result"]["BIC"] - out["M0_result"]["BIC"]

    # pre-registered threshold, quoted verbatim (1 dof case; this is the Omega_m dof)
    out["threshold_preregistered_quote"] = (
        "PRE_REGISTRATION.md Sec 2: 'Thresholds (computed with scipy.stats.chi2.isf): "
        "1 dof: Delta chi2 = 9.00 is 3sigma, and 25.00 is 5sigma; "
        "2 dof: Delta chi2 = 11.83 is 3sigma, and 28.74 is 5sigma.'"
    )
    out["threshold_1dof_3sigma"] = 9.00
    out["threshold_2dof_3sigma"] = 11.83
    out["M0_vs_LCDM_fitted_1dof_verdict"] = (
        "NOT REJECTED" if out["delta_chi2_M0_minus_LCDM_fitted"] <= 9.00 else "REJECTED (>=3sigma, 1 dof)"
    )

    # --- negative controls (Om_L away from 0.68885) ---
    controls = {}
    for om_l in (0.5, 0.72, 0.65):
        om_m = 1.0 - om_l
        c = total_chi2(om_m)
        controls[f"Omega_L={om_l:.2f}"] = {"chi2": c["chi2"], "delta_vs_fitted": c["chi2"] - fitted["chi2"]}
    out["negative_controls"] = controls

    # --- disclosed variant: full 1701-row Pantheon+ sample INCLUDING calibrators ---
    # (kept only as a disclosed comparison per advisor review; NOT the primary result,
    # because IS_CALIBRATOR==1 rows are Cepheid-anchor SNe, not Hubble-flow SNe.)
    m0_full = total_chi2(OMEGA_M_FROZEN, sn_variant="full_1701")
    r_full = minimize_scalar(lambda Om: total_chi2(Om, sn_variant="full_1701")["chi2"],
                              bounds=(0.05, 0.95), method="bounded", options={"xatol": 1e-10})
    fitted_full = total_chi2(float(r_full.x), sn_variant="full_1701")
    n_total_full = N_BAO + 1701
    out["variant_full_1701_including_calibrators"] = {
        "note": "DISCLOSED VARIANT, NOT PRIMARY: includes 77 IS_CALIBRATOR==1 Cepheid-host SNe, "
                "whose MU_SH0ES is tied to CEPH_DIST (a distance-ladder anchor), not a Hubble-flow "
                "distance modulus; fitting them with the same shape+offset model is a likelihood "
                "error. Reported only to show the headline result is stable either way.",
        "n_sn": 1701, "n_total": n_total_full,
        "M0_chi2": m0_full["chi2"], "LCDM_fitted_Omega_m": float(r_full.x),
        "LCDM_fitted_chi2": fitted_full["chi2"],
        "delta_chi2_M0_minus_LCDM_fitted": m0_full["chi2"] - fitted_full["chi2"],
    }

    # --- screening / fifth-force bridge ---
    out["screening_section"] = {
        "fifth_force_data_status": "ABSENT: 2002.11761 e-print tarball fetched and inspected "
            "(data/real2/fifth_force/2002.11761.tar.gz, sha256 in data manifest) -- contains only "
            "FB_ISL_pdf.tex + 9 figure PDFs, NO tabular alpha-lambda exclusion data (grep for "
            "tabular/table environments returned nothing). The Eot-Wash exclusion curve exists only "
            "as a plotted figure; it was not digitised, per ground rules.",
        "unit_bridge_status": "NO BRIDGE. The harness's mu_sym is a dimensionless BVP parameter in "
            "workshopcosmo.run_symmetron_screening_simulation with no documented map to a physical "
            "mass, length, or coupling scale. Confirmed by grep: `grep -n mu_sym workshopcosmo.py` "
            "-> lines 320 (default arg mu_sym: float = 1.0), 361 and 383 (used directly, dimensionless, "
            "inside `((rho/(m_scale**2)) - mu_sym**2)*psi + mu_sym**2*psi**3`); no unit-conversion code "
            "for mu_sym exists anywhere in that file or in scripts/param_loop_sim.py. Therefore even if "
            "a machine-readable Eot-Wash alpha-lambda table existed, "
            "it could not be mapped to a bound on mu_sym today -- 'cannot map' (exactly as pre-flagged "
            "in the ground rules).",
        "does_any_data_prefer_the_symmetron_sector": "No. No fifth-force dataset and no other dataset "
            "in this session's manifest tests mu_sym at all (round3 chi2_independent.json probes_note, "
            "reconfirmed by this script's M2 vs M0 identical chi2). The symmetron sector is UNCONSTRAINED "
            "by every dataset in hand, so its removal in M0 costs nothing and is preferred on parsimony "
            "(AIC/BIC) alone.",
    }

    # --- LeanMaster P2 note (quoted) ---
    out["leanmaster_P2_note"] = {
        "quote_p62": "Within the programme's own T-duality R -> alpha'/R with alpha' = s^2, the "
            "self-dual radius is the unique positive fixed point R = s (p62_selfdual_fixed_point), "
            "where the KK and winding scales coincide (p62_towers_coincide). With the standard KK "
            "normalisation m_n = n/R, this gives kappa = 1: P1 is the programme's prediction, not one "
            "convention among several. (LEANMASTER_CONSTRAINTS, item 'P6.2: T-duality fixes kappa = 1')",
        "quote_p1_excluded": "P1 is excluded. R exceeds each of the three bounds. By the pre-registered "
            "decision rule (T1 fails), P1 is excluded. The failure is not a rounding effect: R > 1.5 x "
            "30 um and R > 1.06 x 44 um. (LEANMASTER_CONSTRAINTS, item 'P1: Extra-dimension prediction')",
        "meaning_for_pre_registered_[21,105]um_bracket": (
            "PRE_REGISTRATION.md P2 registers the bracket [21.0, 105.3] um for s = sqrt(l_P * c/H0), "
            "explicitly because 'the O(1) factor is not derived' (reduced vs non-reduced Planck-mass "
            "convention contributes (8*pi)^(1/4) = 2.239). LeanMaster's p62_selfdual_fixed_point now "
            "shows that WITHIN the programme's OWN T-duality (once alpha' = s^2 is assumed -- itself a "
            "Tier C hypothesis change per the framing rule, not a K3xT2 derivation), kappa is NOT a free "
            "O(1) convention choice: self-duality forces kappa = 1 exactly. This does not shrink the "
            "registered [21,105] um bracket (that bracket already comes from an independent, unrelated "
            "O(1) ambiguity -- the reduced/non-reduced Planck-mass convention, (8*pi)^(1/4) -- which "
            "T-duality's kappa=1 result does not touch or resolve). What it DOES change: it removes the "
            "'is kappa a free tuning knob within T-duality' escape route -- P1's central value (kappa=1, "
            "R=s=47.008 um) is THE unique T-duality prediction, and LeanMaster separately proves that "
            "value is excluded by all three cited bounds (Eot-Wash 30um, Yukawa-range 38.6um, MVV 44um). "
            "The [21,105] um bracket itself is untouched: its lower half (21.0-38.6 um approx, up to the "
            "tightest cited bound) still formally survives per PRE_REGISTRATION.md's own P2 language "
            "('Only the lower part of the bracket survives'), but only as a Planck-mass-convention "
            "artifact, not as anything T-duality favors -- T-duality's OWN preferred point (kappa=1, "
            "47.008 um) sits OUTSIDE the surviving lower half. PRE_REGISTRATION.md is NOT edited by this "
            "script, per ground rules; this is a read-only report of what the two independent facts "
            "(the Lean kappa=1 result and the pre-registered bracket) jointly imply."
        ),
    }

    out["framing_rule_compliance"] = (
        "M0 is a HYPOTHESIS CHANGE: Omega_Lambda is frozen at an imported Planck value (not derived "
        "from K3xT2), c4_pta_product is set to the GR value 0 (not derived), and the symmetron sector "
        "is removed by fiat (not derived). No statement in this script or its output claims K3xT2 "
        "predicts or fixes mu_sym, c4, or Omega_Lambda."
    )

    out_path = os.path.join(HERE, "m0_result.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=1, default=float)
    print(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
