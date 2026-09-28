"""Assemble audit/tda_validation/quantum_fluid/report.json from the result JSONs
(every number is read from a results/*.json written by a script; verdict rules
are the pre-stated ones in expectations.json, commit a0fd19a).
Command: prlimit --as=8589934592 -- .venv-tda/bin/python build_report.py
"""
import json
import os
import sys

import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qf_common as qc  # noqa: E402

QF = os.path.abspath(os.path.join(qc.HERE, ".."))


def J(name):
    p = os.path.join(qc.RESULTS, name)
    return json.load(open(p)) if os.path.exists(p) else None


def test3(stm):
    F = stm["fields"]
    f20 = F["20kOe"]
    keys = sorted(F, key=lambda k: F[k]["H_kOe"])
    Hs = [F[k]["H_kOe"] for k in keys]
    iqr = [F[k]["h0_pooled"]["iqr_over_median"] for k in keys]
    p6 = [F[k]["psi6_global_mean"] for k in keys]
    cr = {k: F[k]["count_ratio_median"] for k in keys}
    argmin_iqr = Hs[int(np.argmin(iqr))]
    argmax_p6 = Hs[int(np.argmax(p6))]
    v = {}
    v["E3a"] = {"observed_a_TDA_over_a_tri_20kOe": f20["a_TDA_over_a_tri"], "a_tri_nm": f20["a_tri_nm"],
                "a_TDA_nm": f20["a_TDA_sqrt3_median_h1_death_nm"], "secondary_2median_H0_over_a_tri": f20["a_from_2_median_h0_over_a_tri"],
                "rule": "within +/-10%", "PASS": bool(abs(f20["a_TDA_over_a_tri"] - 1) <= 0.10)}
    sel = [f"{h}kOe" for h in (10, 15, 20, 25, 30)]
    v["E3b"] = {"count_ratio_median_by_field": cr, "rule": "in [0.90,1.25] at 10-30 kOe",
                "PASS": bool(all(0.90 <= cr[k] <= 1.25 for k in sel)),
                "note_3kOe": "published: minima 10-15%% above expected vortex number at 3 kOe; observed ratio %.3f" % cr["3kOe"]}
    v["E3c"] = {"H_kOe": Hs, "h0_iqr_over_median": iqr, "argmin_H_kOe": argmin_iqr,
                "rule": "min at 15/20/25 kOe and value at 3 and 70 kOe larger than at 20 kOe",
                "PASS": bool(argmin_iqr in (15, 20, 25) and F["3kOe"]["h0_pooled"]["iqr_over_median"] > f20["h0_pooled"]["iqr_over_median"]
                             and F["70kOe"]["h0_pooled"]["iqr_over_median"] > f20["h0_pooled"]["iqr_over_median"])}
    v["E3d"] = {"psi6_global_by_field": dict(zip(keys, p6)), "argmax_H_kOe": argmax_p6,
                "published": {"3kOe": 0.035, "20kOe": 0.62, "above_50kOe": "~0"},
                "sub_rules": {"max_at_15_20_25": bool(argmax_p6 in (15, 20, 25)),
                              "value_20kOe_within_0.62+-0.20": bool(abs(F["20kOe"]["psi6_global_mean"] - 0.62) <= 0.20),
                              "value_3kOe_below_0.2": bool(F["3kOe"]["psi6_global_mean"] < 0.2)}}
    v["E3d"]["PASS"] = bool(all(v["E3d"]["sub_rules"].values()))
    rc = f20["random_control"]
    v["E3e"] = {"random_h0_iqr_over_median_by_field": {k: F[k]["random_control"]["h0_pooled"]["iqr_over_median"] for k in keys},
                "random_psi6_by_field": {k: F[k]["random_control"]["psi6_global_mean"] for k in keys},
                "PASS": bool(all(F[k]["random_control"]["h0_pooled"]["iqr_over_median"] >= 0.4 and F[k]["random_control"]["psi6_global_mean"] < 0.1 for k in keys))}
    rho = stats.spearmanr(iqr, p6).correlation
    v["descriptive"] = {"spearman_rho_TDA_spread_vs_psi6_across_11_fields": float(rho),
                        "a_TDA_over_a_tri_all_fields": {k: F[k]["a_TDA_over_a_tri"] for k in keys}}
    return v


def test1(xa):
    v = xa["verdicts"]
    tda = xa["tda"]
    ext = xa["E1a_extrapolation"]
    crossings = {L: {w: tda[L][w]["T_TDA"] for w in ("window_A", "window_B", "window_C")} for L in tda}
    shuf = {L: tda[L]["E1e_classifier_on_shuffled_windowA"]["T_TDA"] for L in tda if "E1e_classifier_on_shuffled_windowA" in tda[L]}
    b0_site = {L: [tda[L]["E1e_shuffle"][k]["real_b0_theta_le_0_per_site"] for k in tda[L]["E1e_shuffle"]] for L in tda if tda[L]["E1e_shuffle"]}
    return {
        "E1a_helicity_crosscheck": {"T_Ups_by_L": dict(zip(map(str, ext["L"]), ext["T_Ups"])),
                                    "T_Ups_err_by_L": {L: xa["E1a_helicity"][L]["T_Ups_crossing_err"] for L in xa["E1a_helicity"]},
                                    "fit": ext["fit"], "T_inf": ext["T_inf"], "reference_T_BKT": 0.89289,
                                    "PASS": xa["E1a_pass"],
                                    "scope": "validates the Monte Carlo and the independent (non-TDA) observable, not the TDA pipeline"},
        "E1b_primary_lower_star_crossover": {"T_TDA_windowA_by_L": {L: crossings[L]["window_A"] for L in crossings},
                                             "T_TDA_boot_std_by_L": {L: tda[L]["window_A"]["T_TDA_boot_std"] for L in tda},
                                             "rule_PASS": v["E1b_pass"],
                                             "VERDICT": "NOT ATTRIBUTABLE (counted as FAIL as a topology validation): the bracket/drift rule is met, but the pre-stated negative control E1e fails and the classifier trained on SITE-SHUFFLED fields reproduces the crossing",
                                             "classifier_on_shuffled_T_TDA_by_L": shuf},
        "E1c_training_window_control": {"crossings_by_L": crossings, "shift_B_minus_A_L64": v["E1c_shift_B_minus_A"],
                                        "shift_C_minus_A_L64": v["E1c_shift_C_minus_A"], "PASS_at_L64": v["E1c_pass"],
                                        "note": "at L=32 window C moves the crossing by %.3f (would fail); rule was stated for L=64" % (crossings["32"]["window_C"] - crossings["32"]["window_A"])},
        "E1d_vortex_betti_identity": {"median_spearman_rho_T_0.70_0.80": v["E1d_median_rho"], "per_T": tda["64"]["E1d_corr"], "PASS": v["E1d_pass"]},
        "E1e_negative_control_site_shuffle": {"min_diff_over_se_T_le_1": v["E1e_min_diff_over_se_T_le_1"],
                                              "b0_theta_le_0_per_site_real_range_L64": [min(b0_site["64"]), max(b0_site["64"])],
                                              "per_T_L64": tda["64"]["E1e_shuffle"], "PASS": v["E1e_pass"]},
        "E1f_alpha_vortex_pairing": {"S_pair_L64": {k: x["S_pair"] for k, x in tda["64"]["E1f_alpha_pairing"].items()},
                                     "S_pair_0.85": v["E1f_S_pair_0.85"], "S_pair_1.50": v["E1f_S_pair_1.50"], "PASS": v["E1f_pass"],
                                     "note": "S_pair is null at T=0.70 because the random null had f_short = 0 in all configurations (division by zero); reported as null, not infinite"},
    }


def test2(gp, edge):
    prim = {k: gp[k] for k in ("Omega0.70", "Omega0.80", "Omega0.90")}
    ext = {k: gp[k] for k in ("Omega0.70_ext", "Omega0.80_ext", "Omega0.90_ext") if k in gp}
    def e2a(sub):
        return {k: {"b0_lower_star": r["lower_star_b0_at_0.1n0"], "N_winding": r["N_winding_r<0.7R"],
                    "diff": r["E2a_diff_b0_minus_Nwinding"], "vortices_just_outside_edge_(0.7R,0.7R+0.25]": edge[k]["n_vortex_centres_in_(0.7R, 0.7R+0.25]"]} for k, r in sub.items()}
    def e2b(sub):
        return {k: {"central_r<0.5R_over_nF": r["feynman"]["density_r<0.5R_over_nF"], "N_central": r["feynman"]["N_r<0.5R"],
                    "r<0.7R_over_nF": r["feynman"]["density_r<0.7R_over_nF"]} for k, r in sub.items()}
    def e2c(sub):
        return {k: {"a_TDA_over_a_F": r["alpha"]["a_TDA_over_a_F"], "a_F": r["feynman"]["a_F"],
                    "2median_H0_over_a_F": r["alpha"]["2median_h0_over_a_F"], "H0_iqr_over_median": r["alpha"]["h0"]["iqr_over_median"],
                    "psi6_interior": r["alpha"]["psi6_interior_mean_abs"], "interior_coordination_counts": r["posthoc_interior_coordination_counts"],
                    "random_H0_iqr_over_median": r["random_control"]["h0_iqr_over_median_mean"], "random_psi6": r["random_control"]["psi6_interior_mean"]} for k, r in sub.items()}
    A, B, C = e2a(prim), e2b(prim), e2c(prim)
    om0 = gp["Omega0.00"]
    return {
        "runs": {k: {"converged": r["converged"], "steps": r["steps"], "mu": r["mu"], "winding_signs": r["winding_signs"]} for k, r in gp.items()},
        "E2a_h0_equals_winding": {"primary": A, "extension": e2a(ext),
                                  "PASS": bool(all(abs(x["diff"]) <= 1 for x in A.values()) and om0["lower_star_b0_at_0.1n0"] == 0 and om0["N_winding_r<0.7R"] == 0),
                                  "note": "every discrepancy is b0 > N_winding and is <= the number of vortex centres lying within 0.25 a_ho OUTSIDE the analysis disk (their low-density cores enter the masked disk; the all-corners-inside winding count excludes them). Omega=0.8 exceeds the pre-stated tolerance of 1."},
        "E2b_feynman": {"primary": B, "extension": e2b(ext),
                        "PASS_extension_states": bool(all(0.85 <= x["central_r<0.5R_over_nF"] <= 1.05 and x["r<0.7R_over_nF"] <= 1.02 for x in e2b(ext).values())),
                        "extension_convergence": {k: r["converged"] for k, r in ext.items()},
                        "PASS": bool(all(0.85 <= x["central_r<0.5R_over_nF"] <= 1.05 and x["r<0.7R_over_nF"] <= 1.02 for x in B.values())),
                        "note": "one vortex in the central disk changes the ratio by ~0.05-0.12 (8-18 central vortices); the Omega=0.8 step-cap value 0.841 becomes 0.925 in the converged extension run. In the extension states all three central densities are inside [0.85, 1.05] (0.915, 0.925, 0.977), but the edge-depletion sub-rule then fails at Omega=0.7 (18 vortices, 1.050 x Omega/pi inside 0.7 R_eff; one vortex above the 1.02 limit). Headline verdict = primary runs (FAIL)."},
        "E2c_lattice_tda": {"primary": C, "extension": e2c(ext),
                            "PASS_TDA_rule": bool(all(abs(x["a_TDA_over_a_F"] - 1) <= 0.10 and x["H0_iqr_over_median"] <= 0.2 for x in C.values())),
                            "PASS_independent_psi6_crosscheck": bool(all(x["psi6_interior"] >= 0.9 for x in C.values())),
                            "reading": "The TDA rule passes (regular spacing at the Feynman scale), but the independent hexatic check fails: the imaginary-time states contain 4-, 5- and 7-coordinated interior vortices (dislocation/disclination defects; runs not converged by the pre-stated criterion except Omega=0.8 extension). The simulated lattice is therefore NOT a clean Abrikosov triangular lattice, and the alpha-path H0 spread did not flag the defects: it measures spacing regularity, not orientational order."},
        "E2d_negative_controls": {"Omega0_N_winding": om0["N_winding_r<0.7R"], "Omega0_b0": om0["lower_star_b0_at_0.1n0"],
                                  "random_H0_iqr_over_median": {k: x["random_H0_iqr_over_median"] for k, x in C.items()},
                                  "random_psi6": {k: x["random_psi6"] for k, x in C.items()},
                                  "PASS": bool(om0["N_winding_r<0.7R"] == 0 and om0["lower_star_b0_at_0.1n0"] == 0 and all(x["random_H0_iqr_over_median"] >= 0.4 and x["random_psi6"] <= 0.6 for x in C.values()))},
    }


def main():
    xa, gp, edge, stm = J("xy_analysis.json"), J("gpe_tda.json"), J("gpe_edge_trace.json"), J("stm_vortex_tda.json")
    p6v, d6, cal = J("stm_psi6_variants.json"), J("domain06_assessment.json"), J("calib_pipeline.json")
    t1, t2, t3 = test1(xa), test2(gp, edge), test3(stm)
    lf = d6["lattice_fit"]
    F = stm["fields"]
    rep = {
        "title": "TDA pipeline validation against established quantum-fluid topology",
        "expectations_commit": "a0fd19a (expectations.json committed alone, before any computation)",
        "pipeline_functions_tested": {
            "lower_star": "audit/reverse_zero/E5-cmb-tda/cmb_tda.py::betti_curves_from_topology (imported via importlib, not reimplemented)",
            "alpha": "audit/reverse_zero/E5-cosmic-web-tda-scaled/cosmic_web_tda_scaled.py::alpha_persistence, top_bars (betti_curve not needed for the pre-stated statistics)"},
        "status_of_tests": "Tests 1 and 2 are simulations of established theory (known physics, not measurements). Test 3 uses real published STM measurements.",
        "calibration": {"perfect_triangular_H0_death_over_a": cal["alpha_path"]["perfect_triangular_a1"]["h0_death_stats"]["median"],
                        "perfect_triangular_H1_death_over_a": cal["alpha_path"]["perfect_triangular_a1"]["h1_death_stats"]["median"],
                        "poisson_H0_iqr_over_median": cal["alpha_path"]["poisson_iqr_over_median_mean"],
                        "lower_star_dips_b0": cal["lower_star_path"]["gaussian_dips_torus64"]["b0"],
                        "torus_euler": cal["lower_star_path"]["torus_euler"]},
        "test1_xy_bkt": {
            "setup": "2-D XY model, J=1, periodic L x L, L = 32, 64, 128; Wolff single-cluster MC (numba), 2000 thermalisation sweep-equivalents, 400 measurements 2 sweep-equivalents apart, 29 temperatures 0.40-1.60 (step 0.05 plus 0.875, 0.925, 0.975, 1.025); seeds 20260919 + 1000 L + round(1000 T). Max integrated autocorrelation time of the energy (measurement units): " + ", ".join(f"L={L}: {xa['E1a_helicity'][L]['tau_int_energy_meas_units_max']:.2f}" for L in xa["E1a_helicity"]) + ". TDA on 200 (L=32, 64) or 100 (L=128) configurations per temperature.",
            **t1,
            "summary": ("Independent cross-check: the helicity-modulus crossings T_Ups(L) = %s give T_inf = %.4f with the (ln L)^-2 form, against T_BKT = 0.89289 (Komura & Okabe 2012): PASS; this validates the simulation, not the TDA. "
                        "Alpha path (vortex point cloud): the H0 short-bar excess over a matched-density random placement falls from %.1f at T=0.85 to %.2f at T=1.50 (L=64), a clean bound-pair to free-vortex signature: PASS. It is not a T_BKT locator and was not claimed as one. "
                        "Lower-star path (angle field, Cole et al.-style classifier): crossings T_TDA = %s drift down with L and sit in the pre-stated bracket, but the pre-stated negative control fails: b0 per site at the median threshold (L=64) is the same for real and site-shuffled fields at every T <= 1.00 (largest |difference| %.1f standard errors; a small excess of up to %.1f SE appears only at T >= 1.45), and the classifier trained on site-shuffled fields gives %s. The crossover is carried by the one-point angle distribution (magnetisation), not by spatial topology, so it is NOT attributable to topology. The vortex-count/Betti identity (E1d) also fails (median Spearman rho = %.3f).") % (
                            ", ".join(f"{v:.4f} (L={L})" for L, v in t1["E1a_helicity_crosscheck"]["T_Ups_by_L"].items()), t1["E1a_helicity_crosscheck"]["T_inf"],
                            t1["E1f_alpha_vortex_pairing"]["S_pair_0.85"], t1["E1f_alpha_vortex_pairing"]["S_pair_1.50"],
                            ", ".join(f"{v:.3f} (L={L})" for L, v in t1["E1b_primary_lower_star_crossover"]["T_TDA_windowA_by_L"].items()),
                            max(abs(x["diff_over_se"]) for k, x in xa["tda"]["64"]["E1e_shuffle"].items() if float(k) <= 1.0),
                            max(x["diff_over_se"] for x in xa["tda"]["64"]["E1e_shuffle"].values()),
                            ", ".join(f"{v:.3f} (L={L})" for L, v in t1["E1b_primary_lower_star_crossover"]["classifier_on_shuffled_T_TDA_by_L"].items()),
                            t1["E1d_vortex_betti_identity"]["median_spearman_rho_T_0.70_0.80"]),
            "caveats": ["BKT finite-size drift is logarithmic; L <= 128 cannot resolve T_BKT to better than a few percent from any crossing statistic.",
                        "The field was gauge-rotated by the magnetisation angle before filtration (pre-stated); Cole et al. used the raw angle on periodic cubical complexes with persistence images, we used a Freudenthal simplicial torus with Betti curves (the pipeline's output). The E1d rule (within-temperature rank correlation) is stricter than Cole et al.'s Fig. 13 (averages at fixed pair number).",
                        "The Betti curves at the median threshold equal those of an i.i.d. (site-shuffled) field of the same histogram: on this graph the raw-angle sublevel topology is dominated by lattice-scale spin-wave noise.",
                        "L=128 shuffle control was an extension run (xy_tda_features.py --shuffle-only 50), not part of the pre-stated L=64 rule.",
                        "The alpha complex of vortex positions is non-periodic (pipeline limitation); pairs across the periodic boundary are split."]},
        "test2_rotating_bec": {
            "setup": "2-D GPE, hbar=m=omega=1, g=1000, 384^2 grid on [-12,12)^2 (dx=0.0625), imaginary-time Strang split-step Fourier with x/y-split rotation operator, dt 0.005 (8000 steps) then 0.002 (up to 12000), seed 7, Omega = 0, 0.7, 0.8, 0.9. Omega=0 converged; the rotating runs hit the step cap (primary, pre-stated); an extension continued them for up to 40000 steps at dt=0.002 (Omega=0.8 converged, 0.7 and 0.9 still drifting in mu at ~1e-5 per 100 steps).",
            **t2,
            "caveats": ["Imaginary-time evolution from a random-phase start gives metastable, defected vortex arrays; the lattice was not annealed to the triangular ground state within budget.",
                        "Vortex positions are plaquette centres (0.0625 a_ho resolution, <= 3.2% of a_F).",
                        "Only 17-35 vortices lie inside the analysis disk; Feynman ratios carry a counting error of one vortex (~0.05-0.12)."]},
        "test3_real_stm_re6zr": {
            "setup": "Duhan et al., Nat. Commun. 16, 2100 (2025); Zenodo 10.5281/zenodo.14780459 (CC-BY-4.0). 20 successive STS conductance maps (channel Input_7, forward, 256x256 px) per field at 460 mK, H = 3-70 kOe; scan size from the .sxm header. Minima detector fixed in expectations.json (Gaussian low-pass sigma = a_tri/6, disk minimum radius a_tri/3).",
            **t3,
            "posthoc_psi6_variants": p6v["by_field"],
            "summary": ("Real data. At 20 kOe the alpha-path lattice constant is %.1f nm against a_tri = 1.075 sqrt(Phi0/B) = %.1f nm (ratio %.3f): PASS. Minima counts match B*A/Phi0 at every field (ratio %.3f-%.3f); at 3 kOe the excess is %.1f%%, matching the paper's stated 10-15%%. The H0 death spread (IQR/median) is smallest at 20 kOe (%.3f) and larger at 3 kOe (%.3f) and 70 kOe (%.3f): the published re-entrant liquid-solid-liquid sequence is recovered by a purely topological statistic: PASS. Independent psi6 peaks at 20 kOe at %.3f (published ~0.62) but is %.3f at 3 kOe (published ~0.035): E3d FAIL on the 3 kOe sub-rule; no post-hoc reading of psi6 brings 3 kOe below 0.2. Random controls: H0 IQR/median %.2f-%.2f, psi6 <= %.3f: PASS. The TDA spread tracks psi6 only moderately across fields (Spearman rho = %.2f): at 50 and 70 kOe psi6 collapses (%.3f, %.3f) while the H0 spread rises only to %.3f and %.3f (20 kOe: %.3f; 3 kOe liquid: %.3f), i.e. the alpha H0 statistic sees positional spacing regularity, not orientational (hexatic) order.") % (
                F["20kOe"]["a_TDA_sqrt3_median_h1_death_nm"], F["20kOe"]["a_tri_nm"], F["20kOe"]["a_TDA_over_a_tri"],
                min(x["count_ratio_median"] for x in F.values()), max(x["count_ratio_median"] for x in F.values()), 100 * (F["3kOe"]["count_ratio_median"] - 1),
                F["20kOe"]["h0_pooled"]["iqr_over_median"], F["3kOe"]["h0_pooled"]["iqr_over_median"], F["70kOe"]["h0_pooled"]["iqr_over_median"],
                F["20kOe"]["psi6_global_mean"], F["3kOe"]["psi6_global_mean"],
                min(x["random_control"]["h0_pooled"]["iqr_over_median"] for x in F.values()), max(x["random_control"]["h0_pooled"]["iqr_over_median"] for x in F.values()),
                max(x["random_control"]["psi6_global_mean"] for x in F.values()),
                t3["descriptive"]["spearman_rho_TDA_spread_vs_psi6_across_11_fields"],
                F["50kOe"]["psi6_global_mean"], F["70kOe"]["psi6_global_mean"], F["50kOe"]["h0_pooled"]["iqr_over_median"], F["70kOe"]["h0_pooled"]["iqr_over_median"],
                F["20kOe"]["h0_pooled"]["iqr_over_median"], F["3kOe"]["h0_pooled"]["iqr_over_median"]),
            "caveats": ["B = mu0 H assumed (thin film).", "Our minima detector is not the authors' (their Fourier filter parameters are not given numerically); the paper's psi6 normalisation is garbled in the text extraction.",
                        "The scan area was chosen by the experimenters to hold ~120 vortices at most fields, so the count check tests detection plus flux quantisation at one density per field."]},
        "local_file_domain06": {
            "file": d6["file"], "sha256": d6["sha256"], "manifest_source_claim": d6["manifest_entry"]["source"],
            "evidence": [
                "all %d points sit on distinct sites of an ideal triangular lattice with spacing exactly 1.0 (i = %d..%d, j = %d..%d; row pitch sqrt(3)/2, alternate-row offset 0.5)" % (lf["n_points"], *lf["i_range"], *lf["j_range"]),
                "points are stored in exact nested-loop raster order (i-major, j-minor): %s" % lf["storage_order_is_i_major_j_minor_raster"],
                "displacements from the ideal sites are isotropic Gaussian: std = %.4f, %.4f; D'Agostino normality p = %.2f, %.2f; max |residual| %.3f" % (*lf["residual_std"], *lf["normaltest_p"], lf["max_abs_residual"]),
                "z column identically %s; lengths dimensionless (no nm/um scale, no field value), no image, no defects, no dislocations" % d6["z_column_unique_values"],
                "the only other array is a bare scalar pinning_strength = %s" % d6["pinning_strength_value"]],
            "note_on_named_source": "The contents of the NIMS SuperCon database and of any 'Scanning SQUID Microscopy Benchmark (YBCO)' were not checked online in this run; the verdict rests only on the file-content evidence above.",
            "verdict": "SYNTHETIC: a generated 20 x 15 triangular lattice (a = 1) with ~0.04a Gaussian jitter, labelled with real-source names. Not a measurement. Not used as validation data.",
            "alpha_signature_for_record": d6["alpha_signature"]},
        "absent_or_not_used": {
            "huggingface_datasets": "searched https://huggingface.co/api/datasets?search= {gross-pitaevskii, superfluid, vortex, bose-einstein, quantum-turbulence, abrikosov, vortex lattice}: no quantum-fluid vortex dataset (only unrelated LLM/text sets).",
            "the_well_polymathic_ai": "https://huggingface.co/api/datasets?author=polymathic-ai lists 21 datasets (active_matter, MHD, shear_flow, rayleigh_benard, ...): none is a GPE / superfluid dataset.",
            "zenodo": "searched zenodo.org/api/records (type=dataset) for vortex-lattice STM, Gross-Pitaevskii vortex, BEC vortex-lattice images, Bitter decoration; used 10.5281/zenodo.14780459 (Duhan et al.). Other candidates seen but not used: 20728008 (classical-field simulation of vortex-lattice melting, no files listed in the API), 5510351 (generalised GP circulation-statistics sample data).",
            "experimental_BEC_vortex_images": "no public raw image set of Abo-Shaeer et al. 2001 found; test 2 is therefore a known-physics simulation."},
        "citations_checked": ["Komura & Okabe, JPSJ 81, 113001 (2012), arXiv:1210.6116: beta_KT = 1.11996(6)",
                              "Hasenbusch, J. Phys. A 38, 5869 (2005), arXiv:cond-mat/0502556: beta_KT = 1.1199",
                              "Cole, Loges & Shiu, PRB 104, 104426 (2021), arXiv:2009.14231: T_XY ~ 0.9 (20x20, logistic regression on persistence images)",
                              "Donato et al., Phys. Rev. E 93, 052138 (2016) (PRE, not PRB): mean-field XY and phi^4 models, not the 2-D XY BKT transition; not used for a number",
                              "Abo-Shaeer, Raman, Vogels & Ketterle, Science 292, 476 (2001): triangular vortex lattices with > 100 vortices",
                              "Duhan et al., Nat. Commun. 16, 2100 (2025), arXiv:2406.07027: a_tri = 1.075 sqrt(Phi0/B), Psi6 ~0.035 (3 kOe) and ~0.62 (20 kOe) at 460 mK, minima 10-15% above vortex number at 3 kOe"],
    }
    verdicts = {
        "E1a_helicity_crosscheck": t1["E1a_helicity_crosscheck"]["PASS"],
        "E1b_lower_star_crossover": "rule met, NOT ATTRIBUTABLE (control reproduces it)",
        "E1c_window_control_L64": t1["E1c_training_window_control"]["PASS_at_L64"],
        "E1d_vortex_betti_identity": t1["E1d_vortex_betti_identity"]["PASS"],
        "E1e_shuffle_control": t1["E1e_negative_control_site_shuffle"]["PASS"],
        "E1f_alpha_pairing": t1["E1f_alpha_vortex_pairing"]["PASS"],
        "E2a_h0_equals_winding": t2["E2a_h0_equals_winding"]["PASS"],
        "E2b_feynman": t2["E2b_feynman"]["PASS"],
        "E2b_feynman_extension_states": t2["E2b_feynman"]["PASS_extension_states"],
        "E2c_lattice_tda_rule": t2["E2c_lattice_tda"]["PASS_TDA_rule"],
        "E2c_psi6_crosscheck": t2["E2c_lattice_tda"]["PASS_independent_psi6_crosscheck"],
        "E2d_controls": t2["E2d_negative_controls"]["PASS"],
        "E3a_lattice_constant": t3["E3a"]["PASS"], "E3b_count": t3["E3b"]["PASS"], "E3c_reentrant_tda": t3["E3c"]["PASS"],
        "E3d_psi6": t3["E3d"]["PASS"], "E3e_random_control": t3["E3e"]["PASS"]}
    rep["verdict_table"] = verdicts
    rep["execution_notes"] = [
        "Every python command ran under prlimit --as=8589934592 and timeout. Jobs longer than the 600 s foreground limit (XY MC for L=128, TDA feature extraction, GPE solves) ran as nohup background processes with timeouts of 3000-10800 s, not in the foreground.",
        "numba (MC only) was installed into the separate venv /mnt/disks/disk-socrateai-local-1/venv-tdaval; all GUDHI/pipeline code ran in .venv-tda unchanged.",
        "Extensions added after the pre-stated runs, labelled as such in scripts and JSON: GPE continuation (--resume, _ext), L=128 shuffle-only control, stm_psi6_variants.py, gpe_edge_trace.py, psi6 perfect-lattice reference in gpe_tda.py. None changes a pre-stated verdict.",
        "The GPE solver was restarted once (after ~7 min) to add checkpointing; the restarted runs are the ones reported."]
    rep["overall_assessment"] = (
        "The alpha-complex path (alpha_persistence/top_bars) recovers established quantum-fluid topology in the cases where the answer is a spacing or pairing scale: "
        "vortex-antivortex binding in the XY model (bound-pair excess vs random), the Feynman-scale vortex spacing in a rotating condensate, and, on real STM data, the Abrikosov lattice constant 1.075 sqrt(Phi0/B) and the published re-entrant order-disorder-order sequence in a-Re6Zr. "
        "Its H0 spread does not see orientational order: it did not flag defected GPE lattices and did not follow the psi6 collapse at 50-70 kOe. "
        "The lower-star path (betti_curves_from_topology) counts components correctly on a smooth field (GPE vortex cores: exact at Omega=0 and 0.9, off by edge-straddling cores at 0.7 and 0.8), but on the raw XY angle field its Betti curves are indistinguishable from a site-shuffled field, so the Cole et al.-style BKT crossover it produces is not attributable to topology. "
        "For the cosmology use: a positive lower-star result on a rough field needs a site-shuffle (one-point-matched) null, which the CMB analysis already has in the form of Gaussian sims of the same spectrum; alpha-path H0 statistics should not be read as evidence of orientational/structural order without an independent orientational check.")
    with open(os.path.join(QF, "report.json"), "w") as fh:
        json.dump(rep, fh, indent=1)
    print(json.dumps(verdicts, indent=1))


if __name__ == "__main__":
    main()
