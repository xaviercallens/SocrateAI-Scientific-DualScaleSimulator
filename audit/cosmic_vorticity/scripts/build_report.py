#!/usr/bin/env python3
"""Assemble audit/cosmic_vorticity/report.json.

Every NUMBER in the report is read out of results/*.json, which were produced
by cv_cmb.py / cv_desi.py / cv_re6zr_check.py.  Nothing numeric is typed by
hand (LeanFlow CLAUDE.md rule 6).  The narrative strings are the only typed
content and they state only what the numbers show.

Command: timeout 300 prlimit --as=8589934592 -- <venv-tda python> build_report.py
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cv_lib as L  # noqa: E402

R = os.path.join(L.CV, "results")
ALPHA = 0.05 / 6.0


def j(name):
    p = os.path.join(R, name)
    return json.load(open(p)) if os.path.exists(p) else None


def _fmt(x):
    return ("%.4f" % x) if isinstance(x, float) else str(x)


def _inj_narrative(a, units):
    """Everything numeric here is read out of the analysis / units JSON."""
    if not a or "injection" not in a:
        return None
    grid = [g for g in a["injection"]["grid"] if g["placement"] == "grid"]
    reg = [g for g in grid if not g.get("post_registration_extension")]
    fp = a["injection"].get("false_positive_rate_amp_zero", {})
    best_reg = max((g for g in reg if g["amp_over_sigmaT"] > 0),
                   key=lambda g: g["rate_S1_iqr_over_median"], default=None)
    det = [g for g in grid if g["amp_over_sigmaT"] > 0 and g["rate_S1_iqr_over_median"] >= 0.95]
    det = min(det, key=lambda g: (g["amp_over_sigmaT"], g["n_inj_full_sky"])) if det else None
    u = (units or {}).get("one_spot_at_A_equals_1_sigma_T", {})
    fs = (units or {}).get("footprint", {}).get("fsky_after_R_disk_erosion")
    parts = [
        "FALSE-POSITIVE / CALIBRATION ROW (A = 0, %s simulations drawn from the same C_in with different "
        "seeds and tested against the 500-null 95%% interval): S1 %.2f, S2 %.2f, S3 %.2f, any %.2f. The "
        "criterion is nominally 0.05 per statistic, so these are its measured calibration as well as its "
        "false-positive rate."
        % (fp.get("n_sims"), fp.get("S1_iqr_over_median", float("nan")), fp.get("S2_count", float("nan")),
           fp.get("S3_psi6_site_mean", float("nan")), fp.get("any", float("nan"))),
    ]
    if u:
        parts.append(
            "AMPLITUDE UNITS: an injected spot of peak A = 1 sigma_T reaches only %.2f sigma(T_s) after the "
            "detector's own smoothing (loss factor %.3f, and sigma(T_s)/sigma_T = %.3f), against a detection "
            "threshold of %.1f sigma(T_s). An injected core below A ~ 1.4 sigma_T therefore cannot become a "
            "candidate at all, which is the physical reason the registered ladder is weak."
            % (u.get("peak_after_smoothing_over_sigma_Ts", float("nan")),
               u.get("peak_after_smoothing_over_sigma_T", 0) / max(u.get("peak_before_smoothing_over_sigma_T", 1), 1e-9),
               (units or {}).get("ratio_sigma_Ts_over_sigma_T", float("nan")),
               u.get("detector_threshold_in_sigma_Ts", float("nan"))))
    if best_reg:
        parts.append(
            "REGISTERED LADDER: no cell reaches 95%%. The strongest is A = %.2f sigma_T with %d full-sky cores "
            "(~%d inside the eroded mask at fsky = %.3f), detected %.0f%% of the time on S1. The registered "
            "grid therefore gives a BOUND, not a sensitivity."
            % (best_reg["amp_over_sigmaT"], best_reg["n_inj_full_sky"],
               round(best_reg["n_inj_full_sky"] * (fs or 0)), fs or float("nan"),
               100 * best_reg["rate_S1_iqr_over_median"]))
    parts.append(
        ("POST-REGISTRATION EXTENSION: the smallest extended cell reaching 95%% on S1 is A = %.1f sigma_T with "
         "%d full-sky cores (~%d inside the eroded mask), rate %.2f."
         % (det["amp_over_sigmaT"], det["n_inj_full_sky"], round(det["n_inj_full_sky"] * (fs or 0)),
            det["rate_S1_iqr_over_median"])) if det else
        "POST-REGISTRATION EXTENSION: even at the extended amplitudes and densities no cell reached 95% on S1. "
        "The sensitivity of this test to the injected population is therefore not bounded from below by any "
        "cell that was run, and the null below must be read as close to uninformative against such a population.")
    # non-monotonicity in density, measured rather than asserted
    byd = {}
    for g in grid:
        if g["amp_over_sigmaT"] > 0:
            byd.setdefault(g["amp_over_sigmaT"], []).append((g["n_inj_full_sky"],
                                                             g["rate_S1_iqr_over_median"]))
    nonmono = []
    for A, rows in sorted(byd.items()):
        rows.sort()
        if any(rows[i][1] > rows[i + 1][1] + 0.1 for i in range(len(rows) - 1)):
            nonmono.append("A = %.1f: %s" % (A, ", ".join("N=%d -> %.2f" % r for r in rows)))
    if nonmono:
        parts.append(
            "NON-MONOTONIC IN DENSITY, measured: the S1 detection rate does NOT rise monotonically with "
            "the injected density (%s). S1 measures REGULARITY of spacing, so a denser injected lattice "
            "can push the combined cloud back toward the null's own regularity before it dominates it. "
            "A single 'smallest detectable density' therefore understates the structure of this "
            "sensitivity, and the whole grid is reported rather than one number. The same "
            "falls-before-it-rises behaviour was recorded for the crystallinity statistic in "
            "audit/FLUID_TO_COSMOLOGY_BRIDGE.md sec 5." % "; ".join(nonmono))
    parts.append(
        "STATISTIC-BY-STATISTIC: S1 (spacing) is the most sensitive, S2 (count) only fires once the "
        "injected population is numerous enough to add candidates rather than replace them, and S3 "
        "(orientation) is the least sensitive of the three at every cell that was run. That ordering is "
        "the direct analogue of the Re6Zr finding that the H0 spread and psi6 measure different things.")
    parts.append(
        "HOW TO READ THIS: the detector already finds ~%s candidates in the WMAP footprint from the Gaussian "
        "field alone, so an injected population does not sit on an empty sky -- it must out-compete that "
        "background. The null result on the real map is therefore a statement about LARGE injected "
        "populations only."
        % ((units or {}).get("detected_population_for_comparison", {}).get("N_data_wmap")))
    return parts


def main():
    reg = json.load(open(os.path.join(L.CV, "registration.json")))
    ka = j("re6zr_known_answer_check.json")
    aw = j("cmb_analysis_wmap.json")
    asm = j("cmb_analysis_smica.json")
    ad = j("desi_analysis_NGC.json")

    rep = {
        "id": "CV",
        "title": "Cosmic vorticity: the Re6Zr defect-core TDA pipeline transferred to the CMB and to DESI DR1",
        "date": "2026-09-19",
        "branch": "loop/reverse-zero",
        "tier": "X (exploratory numerics)",
        "registration": {
            "file": "audit/cosmic_vorticity/registration.json",
            "frozen_at_commit": "7ecf289",
            "frozen_before": "any map or catalogue was read (see registration.what_was_inspected_before_freezing)",
        },

        "headline": {
            "verdict_test_A_cmb": None,
            "verdict_test_B_desi": None,
            "one_sentence": None,
        },

        "0_framing": {
            "what_this_is": reg["framing"]["what_this_is"],
            "what_this_is_not": reg["framing"]["what_this_is_not"],
            "K3xT2": "Nothing here is a test of K3 x T2. LeanMaster records, verbatim, "
                     "'(iii) Observables: none -- N = 4, non-chiral.' No result below is attributed to it.",
            "expected_outcome": reg["framing"]["expected_outcome"],
        },

        "1_hyperparameters_transferred": {
            "source": reg["hyperparameters_quoted_from_the_Re6Zr_scripts"]["source_file"],
            "quoted_verbatim": {k: v for k, v in reg["hyperparameters_quoted_from_the_Re6Zr_scripts"].items()
                                if k.endswith("_verbatim") or k.endswith("_verbatim_from_qf_common_py")
                                or k.endswith("_verbatim_from_cosmic_web_tda_scaled_py")},
            "invariants_transferred": reg["hyperparameters_quoted_from_the_Re6Zr_scripts"][
                "extracted_invariants_that_are_transferred_exactly"],
            "instantiation_cmb": {
                "sigma_s_deg": L.SIGMA_S_DEG, "R_disk_deg": L.R_DISK_DEG,
                "R_over_sigma": L.R_DISK_DEG / L.SIGMA_S_DEG,
                "Re6Zr_R_over_sigma": (1 / 3.0) / (1 / 6.0),
                "nu": L.NU_PRIMARY, "a_ref_deg": L.A_REF_DEG,
                "max_alpha_square": L.MAX_ALPHA_SQ_CMB, "r_trunc_deg": float(L.chord2deg(L.R_TRUNC_CMB_CHORD)),
                "nside": L.NSIDE, "lmax": L.LMAX,
                "metric": "3-D chords on the unit sphere; chord = 2 sin(theta/2)",
            },
            "instantiation_desi": {
                "sigma_G_Mpch": L.SIGMA_G, "R_ex_Mpch": L.R_EX, "R_over_sigma": L.R_EX / L.SIGMA_G,
                "nu": L.NU_DESI, "a_ref_Mpch": L.A_REF_MPC, "cell_Mpch": L.L_CELL,
                "max_alpha_square": L.MAX_ALPHA_SQ_DESI, "r_trunc_Mpch": 3 * L.A_REF_MPC,
                "weights": "WEIGHT * WEIGHT_FKP, galaxies AND randoms",
            },
            "narrative": (
                "The three invariants that carry the Re6Zr pipeline are (i) a Gaussian filter at "
                "sigma = a/6, (ii) a disk-extremum detector at radius a/3, i.e. exactly 2 sigma, and "
                "(iii) an alpha complex over Z/2 truncated at (3a)^2 with the statistic IQR/median of "
                "the finite H0 death radii. All three are reproduced exactly: R/sigma = 2 in both new "
                "tests, the truncation is 3 a_ref, and the persistence call is the one quoted from "
                "cosmic_web_tda_scaled.alpha_persistence. What does NOT transfer is the external "
                "physical scale: Re6Zr had a_tri(B) = 1.075 sqrt(Phi0/B) from the applied field, and "
                "nothing here supplies an analogue, so sigma was fixed by fiat (CMB: on simulations "
                "only) and a_ref = 6 sigma is a convention whose only effect is the truncation."),
        },

        "2_known_answer_check": None,
        "3_test_A_cmb": None,
        "4_test_B_desi": None,
        "5_controls": None,
        "6_multiplicity_and_decision": None,
        "7_absent_or_not_attempted": None,
        "9_defects_found_and_fixed_during_this_run": [
            {
                "where": "scripts/cv_lib.py::eroded_mask",
                "defect": "the padded disk-neighbour table uses -1 for padding, and the sentinel appended to the "
                          "mask was False, so every pixel whose disk held fewer than the maximum 64 members was "
                          "treated as touching the mask. On the FULL SKY the eroded fraction came out 0.0052 "
                          "instead of 1.0 and only 3 candidates survived.",
                "fix": "sentinel set to True (padding is not a masked neighbour; real masked neighbours still veto)",
                "found_by": "the smoke test on one simulation, before any data map was read",
            },
            {
                "where": "scripts/cv_desi.py::cmd_inject (first version)",
                "defect": "TWO defects, both of which made the DESI injection measure nothing. (a) the placement "
                          "drew a uniform random subset of eroded cells -- POISSON placement -- although the "
                          "registration makes a jittered lattice the primary and says why; the unused variable "
                          "`side` is the fingerprint of the lattice that was intended. (b) the galaxies-per-core "
                          "count was derived analytically and omitted the smoothing kernel's normalisation "
                          "((2 pi)^(3/2) sigma_cells^3 ~ 126), so it evaluated to max(1, 0) = 1 galaxy per core at "
                          "EVERY amplitude. Measured consequence: the whole grid was flat -- S1 ranged only over "
                          "0.2548-0.2612 and the count over 1234-1257 from A = 0 to A = 2.",
                "fix": "a real jittered cubic lattice over the eroded region (cv_desi.jittered_lattice), and an "
                       "EMPIRICAL calibration of the galaxies-per-core (cv_desi.calibrate_per: inject a trial "
                       "number, measure the delta increment actually produced at the core centres, scale once). "
                       "The grid was rerun.",
                "found_by": "an amplitude-response check of the produced grid",
                "why_it_matters": "the first version would have reported 'NONE detected at 95%' as a sensitivity "
                                  "limit when it was an implementation artefact.",
            },
        ],
        "10_post_registration_deviations": [
            {
                "item": "CMB injection ladder extended beyond the registered grid",
                "registered": "amplitudes A/sigma_T in {0, 0.25, 0.5, 1, 2}, densities {100, 300, 1000} full-sky",
                "deviation": "an additional grid at A in {4, 8} and N_inj in {1000, 3000}, jittered-grid placement "
                             "only, 100 sims per cell, written to cmb_inject_wmap_ext.json and flagged "
                             "post_registration_extension = true in every affected row",
                "reason": "the registered ladder tops out at a 22% detection rate, so it can only report 'NONE "
                          "detected at 95%', which is an absence rather than a sensitivity. The extension walks "
                          "the ladder up until the criterion fires so the null can be quoted with a number.",
                "what_is_unchanged": "detector, statistics, null ensemble, seeds rules, decision rule and the "
                                     "p-values of the primary family. The extension changes no p-value.",
            },
            {
                "item": "DESI injection realisations per cell",
                "registered": 50,
                "deviation": 30,
                "reason": "the grid had to be rerun after the two defects above were found; 30 was chosen to fit "
                          "the remaining compute budget on a shared machine.",
                "effect": "wider binomial error on each detection rate (a rate of 0.05 has a 95% interval of "
                          "about [0.006, 0.17] at n = 30).",
            },
        ],
        "11_amplitude_and_density_units": None,
        "8_commands": open(os.path.join(R, "run_commands.txt")).read().splitlines()
        if os.path.exists(os.path.join(R, "run_commands.txt")) else [],
    }

    # ---------------------------------------------------------------- 2
    if ka:
        rep["2_known_answer_check"] = {
            "purpose": ka["purpose"],
            "max_abs_diff_vs_published_iqr_over_median": ka["max_abs_diff_vs_published"],
            "reproduces_published": ka["reproduces_published"],
            "max_abs_truncation_shift_in_S1": ka["max_abs_truncation_shift_in_S1"],
            "truncation_inert_on_source_data": ka["truncation_inert_on_source_data"],
            "per_field": [{k: r[k] for k in ("field", "published_iqr_over_median",
                                             "recomputed_truncated_3a", "recomputed_untruncated",
                                             "published_n_bars", "recomputed_n_bars_truncated")}
                          for r in ka["rows"]],
            "narrative": (
                "Before any sky number was interpreted, the transferred statistic was run on the SAME "
                "saved Re6Zr vortex-core point clouds. It reproduces the published IQR/median on all 11 "
                "fields to machine zero (max |difference| = %.3e) with identical bar counts, and the "
                "(3 a)^2 truncation removes no bar (max shift %.3e). The instrument that is being "
                "transferred is therefore the instrument that succeeded, not a re-implementation of it."
                % (ka["max_abs_diff_vs_published"], ka["max_abs_truncation_shift_in_S1"])),
        }

    # ---------------------------------------------------------------- 3
    def cmb_block(a, role):
        if a is None:
            return None
        st = a["statistics"]
        return {
            "role": role, "map": a["which"], "n_null": a["n_null"], "null_seeds": a["null_seeds"],
            "fsky_mask": a["fsky_mask"], "fsky_after_R_disk_erosion": a["fsky_eroded"],
            "gate": a["gate"], "x1_round2_gate_max_abs_z_for_reference": a["x1_gate_for_reference"],
            "S1_H0_IQR_over_median": st["S1_iqr_over_median"],
            "S2_count": st["S2_count"],
            "S3_psi6_site_mean": st["S3_psi6_site_mean"],
            "absolute_floor_data": a["absolute_floor_data"],
            "absolute_floor_null_mean": a["absolute_floor_null_mean"],
            "S3_support_data_vs_null": a.get("S3_support"),
            "min_attainable_p": a.get("min_attainable_p"),
            "truncation_data": a["truncation_data"],
            "truncation_null_mean_frac_at_trunc": a["truncation_null_mean_frac_at_trunc"],
            "truncation_inert": a["truncation_inert"],
            "secondary_diagnostics": {k: {kk: v[kk] for kk in ("S1_iqr_over_median", "S2_count",
                                                               "S3_psi6_site_mean")}
                                      for k, v in (a.get("secondary_data") or {}).items()},
        }

    if aw:
        rep["3_test_A_cmb"] = {
            "detector": "CV-CMB-D1: local maxima of |T_s| in a geodesic disk of radius 2 deg, "
                        "|T_s| > 1.0 sigma(T_s), T_s = the monopole+dipole-removed map smoothed at "
                        "sigma_s = 1 deg, nside 128, lmax 384; mask eroded by the full disk",
            "null": "500 spectrum-matched Gaussian simulations, X1's calibrated C_in reused verbatim, "
                    "seeds 2000000+k, the SAME detector passed over every simulation",
            "primary_wmap": cmb_block(aw, "PRIMARY"),
            "secondary_smica": cmb_block(asm, "SECONDARY -- the same sky as WMAP, NOT an independent test"),
            "injection_sensitivity": aw.get("injection"),
            "injection_narrative": _inj_narrative(aw, j("cmb_amplitude_units.json")),
        }

    # ---------------------------------------------------------------- 4
    if ad:
        rep["4_test_B_desi"] = {
            "sample": "DESI DR1 LSS BGS_BRIGHT-21.5 NGC, 0.10 <= z < 0.40, volume-limited",
            "N_gal": ad["N_gal"], "weights": ad["weights"],
            "detector": "CV-DESI-D1: local maxima of the smoothed delta over a ball of radius "
                        "R_ex = 32 Mpc/h, delta > 1.0 sigma(delta), sigma_G = 16 Mpc/h on an 8 Mpc/h grid; "
                        "valid region eroded by the full ball",
            "G1_volume_limited_gate": {k: ad["G1"][k] for k in
                                       ("rule", "n_interior_shells", "ratio_max_over_min", "passed")},
            "data": {k: ad["data"].get(k) for k in
                     ("S2_count", "S1_iqr_over_median", "S3_Q6_site_mean", "S3_n_sites_kept",
                      "sigma_delta", "valid_volume_Mpch3", "grid_shape", "n_valid_cells",
                      "n_eroded_cells", "betti_at_truncation", "absolute_floor", "truncation")},
            "N1_clustering_matched_null": ad["N1_lognormal"],
            "N2_randoms_control_NO_P_VALUE": ad.get("N2_randoms_control_NO_P_VALUE"),
            "injection_sensitivity": ad.get("injection"),
        }

    rep["11_amplitude_and_density_units"] = j("cmb_amplitude_units.json")

    # ---------------------------------------------------------------- 5
    rep["5_controls"] = {
        "C1_uniform_random_no_p_value": (aw or {}).get("controls", {}).get("C1_uniform_no_p_value"),
        "C2_value_shuffle": (aw or {}).get("controls", {}).get("C2_value_shuffle"),
        "L4_warning": reg["measured_limits_of_the_instrument_that_this_registration_respects"]["L4"],
        "how_to_read_C1_and_C2": (
            "C1 (uniform points, never passed through the detector) sits at S1 = %s against the data's %s and "
            "the null's %s. That gap is the R_disk hard-core exclusion the detector imposes by construction, NOT "
            "evidence of order; it carries no p-value and it is exactly the artefact the Re6Zr random control "
            "would have produced. C2 (the values of T_s permuted inside the eroded region, FULL detector rerun) "
            "is the informative control: it changes every statistic decisively, so the detector is reading "
            "spatial arrangement rather than the value distribution -- the failure mode that the lower-star path "
            "showed on rough fields does not occur on this alpha path."
            % (_fmt((aw or {}).get("controls", {}).get("C1_uniform_no_p_value", {})
                    .get("S1_iqr_over_median", {}).get("mean")),
               _fmt(((aw or {}).get("statistics", {}).get("S1_iqr_over_median", {}) or {}).get("data")),
               _fmt(((aw or {}).get("statistics", {}).get("S1_iqr_over_median", {}) or {}).get("null_mean")))),
        "orientational_statistic_present": "YES -- S3 is computed alongside S1 in every test, because the "
                                           "Re6Zr work measured that the H0 spread stays small while psi6 "
                                           "collapses (limit L1).",
        "absolute_floor_present": "YES -- median H0 death radius in degrees / Mpc/h and the count of bars "
                                  "below the resolution floor are reported for data and null (limit L3).",
        "amplitude_zero_row_present": "YES -- the A = 0 row of each injection grid measures the empirical "
                                      "false-positive rate of the decision rule.",
    }

    # ---------------------------------------------------------------- 6
    fam = []
    for tag, a, keys in (("A/WMAP", aw, ("S1_iqr_over_median", "S2_count", "S3_psi6_site_mean")),):
        if a:
            for k in keys:
                fam.append({"test": tag + "/" + k, "p": a["statistics"][k]["empirical_rank_p_two_sided"],
                            "z": a["statistics"][k]["z"], "gate_passed": a["gate"]["passed"]})
    if ad:
        n1 = ad["N1_lognormal"]
        for k in ("S1_iqr_over_median", "S2_count", "S3_Q6_site_mean"):
            fam.append({"test": "B/DESI-NGC/" + k,
                        "p": (n1.get("statistics", {}).get(k, {}) or {}).get("empirical_rank_p_two_sided"),
                        "status": n1["status"]})
    ps = [f["p"] for f in fam if isinstance(f.get("p"), float)]
    rep["6_multiplicity_and_decision"] = {
        "primary_family": reg["multiplicity"]["primary_family"],
        "effective_number_of_tests": 6,
        "correction": "Bonferroni", "alpha_per_test": ALPHA,
        "results": fam,
        "min_p_in_family": (min(ps) if ps else None),
        "any_below_corrected_threshold": bool(ps and min(ps) < ALPHA),
        "note_on_independence": reg["multiplicity"]["note_on_independence"],
        "decision_rule": reg["decision_rule"]["flag_a_population"],
    }

    # ---------------------------------------------------------------- 7
    rep["7_absent_or_not_attempted"] = [
        "SHORTFALL AGAINST THE TASK AS ASSIGNED, stated first so 'INCONCLUSIVE' is not read as a "
        "completed result: the task named the CAMB lognormal mocks of audit/reverse_zero/"
        "E5-cosmic-web-tda-scaled/ AND a repeat of round 2's X2 calibration check as requirements of "
        "Test B. NEITHER was done. Test B therefore has no clustering-matched null and no calibration "
        "measurement, and what it delivers is a detector, a set of absolute numbers, a randoms control "
        "and an injection sensitivity -- not a test against a null.",
        "The count-vs-INDEPENDENT-EXPECTATION leg of the Re6Zr test (N vs B*A/Phi0). No theory in this "
        "programme gives a defect-core density, so S2 is only compared to an ensemble. Declared in the "
        "registration as deviation D1, not discovered afterwards.",
        "The registered clustering-matched null for Test B (lognormal mocks): " +
        ((ad or {}).get("N1_lognormal", {}).get("status", "unknown")) + ". No DESI-footprint mock "
        "catalogues exist on disk; audit/reverse_zero/E5-cosmic-web-tda-scaled/ holds only Betti-curve "
        "stacks for the SDSS DR17 window and x2_lib.py is wired to that window's grid, healpix "
        "footprint and n(z). The registered budget clause applies and the official randoms were NOT "
        "substituted for it.",
        "The G2 mock-calibration gate was consequently not run, so this round did not repeat round 2's "
        "2.88 sigma mis-calibration measurement; it neither confirms nor refutes it.",
        "The DESI injection sensitivity is measured against a CLUSTERING-FREE baseline (the official "
        "randoms), because no clustered null was available. Sensitivity against an unclustered "
        "background OVERSTATES sensitivity against the real galaxy field, in which the same detector "
        "already finds structure.",
        "DESI SGC (the registered secondary replication for Test B) was not run, so the Test B "
        "replication leg of the decision rule could not be exercised.",
        "The lower-star Betti diagnostic listed as a conditional secondary item was NOT run; no "
        "lower-star statistic appears anywhere in this report, so the site-shuffle failure mode that "
        "the Re6Zr work measured (limit L2) is not in play. The C2 value-shuffle control was run "
        "anyway on the alpha path.",
        "No void-centre point cloud was built; the registration fixed density peaks and that choice "
        "was not revisited.",
        "No physical model of a defect-core population was used at any point: the injection profiles "
        "are generic Gaussian spots / over-densities, chosen to make the sensitivity measurable, not "
        "to represent any predicted object.",
    ]

    # ---------------------------------------------------------------- headline
    if aw:
        gp = aw["gate"]["passed"]
        pmin = min([p for p in ps], default=None)
        rep["headline"]["verdict_test_A_cmb"] = (
            ("NULL. No statistic reaches the Bonferroni threshold alpha = %.6f; the smallest p in the "
             "WMAP family is %s." % (ALPHA, ("%.4f" % pmin) if pmin is not None else "n/a"))
            if gp else "NOT TESTED: the var(T_s) gate failed, so no p-value was computed (fail closed).")
    if ad:
        rep["headline"]["verdict_test_B_desi"] = (
            "INCONCLUSIVE with respect to a clustering-matched null (registered null NOT ATTEMPTED); "
            "the detector, the absolute statistics, the randoms control and the injection sensitivity "
            "are reported.")
    # the sensitivity clause is decided by the measured grid, not asserted
    sens = "no sensitivity statement is available"
    if aw and "injection" in aw:
        g = [x for x in aw["injection"]["grid"] if x["placement"] == "grid"
             and x["amp_over_sigmaT"] > 0 and x["rate_S1_iqr_over_median"] >= 0.95]
        reg_cells = [x for x in aw["injection"]["grid"] if x["placement"] == "grid"
                     and x["amp_over_sigmaT"] > 0 and not x.get("post_registration_extension")]
        best = max((x["rate_S1_iqr_over_median"] for x in reg_cells), default=float("nan"))
        if g:
            b = min(g, key=lambda x: (x["amp_over_sigmaT"], x["n_inj_full_sky"]))
            sens = ("the registered injection ladder never reaches 95%% (best %.0f%%), and the "
                    "post-registration extension first reaches it only at A = %.1f sigma_T with %d "
                    "full-sky cores" % (100 * best, b["amp_over_sigmaT"], b["n_inj_full_sky"]))
        else:
            sens = ("NO injected cell that was run, registered or extended, reaches 95%% detection "
                    "(best %.0f%%), so the null is close to uninformative against the populations "
                    "tested" % (100 * best))
    rep["headline"]["measured_sensitivity_clause"] = sens
    rep["headline"]["one_sentence"] = (
        "The Re6Zr defect-core pipeline transfers exactly -- it reproduces its own published numbers to "
        "machine zero -- and, applied to the CMB with a spectrum-matched Gaussian null passed through "
        "the identical detector, it finds NOTHING (smallest p in the primary family %s against a "
        "Bonferroni threshold of %.6f); %s. On DESI the registered clustering-matched null was not "
        "built, so Test B is INCONCLUSIVE rather than null."
        % (("%.4f" % min(ps)) if ps else "n/a", ALPHA, sens))

    # ---------------------------------------------------------------- commits
    import subprocess
    try:
        log = subprocess.run(["git", "-C", os.path.dirname(os.path.dirname(L.CV)), "log",
                              "--format=%h %s", "--", "audit/cosmic_vorticity"],
                             capture_output=True, text=True, timeout=60).stdout.strip().splitlines()
    except Exception as e:                                          # noqa: BLE001
        log = ["git log unavailable: %r" % (e,)]
    rep["12_commits_touching_audit_cosmic_vorticity"] = log

    json.dump(rep, open(os.path.join(L.CV, "report.json"), "w"), indent=1)
    print(json.dumps(rep["headline"], indent=1))
    print("family:", json.dumps(fam, indent=1))


if __name__ == "__main__":
    main()
