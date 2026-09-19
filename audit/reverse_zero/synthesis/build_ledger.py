#!/usr/bin/env python3
"""Build audit/reverse_zero/ledger_reverse_zero.json for the reverse-to-zero synthesis.

Exact command (run from anywhere; all paths are absolute below):

  /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python \
    /mnt/disks/disk-socrateai-local-1/dualscale-wt-reverse/audit/reverse_zero/synthesis/build_ledger.py

What it does (no fitting, no randomness, so no seed is needed):
  * reads every number from the committed result JSONs of this branch (loop/reverse-zero);
  * reads the pre-registered thresholds and rules from PRE_REGISTRATION.md (main repo, read-only)
    by line number, asserts the quoted text is present, and recomputes the thresholds with scipy;
  * reads LeanMaster statements VERBATIM from the committed tree at LEANMASTER_COMMIT via
    `git show` (read-only), asserting each quoted substring is on the cited line;
  * reads the external flux-vacua counts from commit ab66bff (k3t2-rigidity-v2 Track E) via
    `git show` in this repository's shared object store (the other worktree is not touched);
  * computes only derived quantities that the synthesis needs (sigma conversions, Bonferroni,
    ratio of the G mu bracket to the imported bound) and writes the ledger.

Framing: M0 is a HYPOTHESIS CHANGE (sectors removed or set to their GR value), not a derivation
from K3 x T2. Nothing here says K3 x T2 predicts or fixes mu_sym, c4_pta_product or Omega_Lambda.
"""
import json
import subprocess
from pathlib import Path

import numpy as np
from scipy import stats

WT = Path("/mnt/disks/disk-socrateai-local-1/dualscale-wt-reverse")
RZ = WT / "audit/reverse_zero"
PREREG = Path("/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/audit/PRE_REGISTRATION.md")
LM = Path("/home/callensxavier_gmail_com/SocrateAI-Scientific-Agora-LeanMaster")
LEANMASTER_COMMIT = "eb791e7"
FLUX_COMMIT = "ab66bff"
OUT = RZ / "ledger_reverse_zero.json"


def load(rel):
    return json.loads((RZ / rel).read_text())


def git(repo, *args):
    return subprocess.run(["git", "-C", str(repo), *args], check=True,
                          capture_output=True, text=True).stdout


# ---------------------------------------------------------------- PRE_REGISTRATION quotes
prereg_lines = PREREG.read_text().splitlines()


def prq(lineno, must):
    line = prereg_lines[lineno - 1]
    assert must in line, (lineno, must, line)
    return {"file_line": f"audit/PRE_REGISTRATION.md:{lineno}", "text": line.strip()}


PRQ = {
    "thresholds_1dof": prq(35, "1 dof: Δχ² = 9.00 is 3σ, and 25.00 is 5σ"),
    "thresholds_2dof": prq(36, "2 dof: Δχ² = 11.83 is 3σ, and 28.74 is 5σ"),
    "R2": prq(28, "Not rejected"),
    "R1": prq(27, "Rejected: Δχ² = +1204.5"),
    "R3": prq(29, "Δχ² = 6.09"),
    "P1_falsified_if": prq(40, "≥ 5σ (Δχ² ≥ 28.74, 2 dof)"),
    "P1_data_1": prq(41, "DESI DR3 or final BAO combined with CMB"),
    "P1_data_2": prq(42, "Euclid"),
    "P1_data_3": prq(43, "Rubin/LSST"),
    "P1_omega_lambda_rule": prq(45, "away from 0.68885 by **≥ 5σ**"),
    "P1_pipeline_action": prq(46, "Pipeline action"),
    "P1_pipeline_action_a": prq(47, "decisive_experiment.py"),
    "P1_pipeline_action_b": prq(48, "report the CPL-vs-Λ Δχ² as it comes out"),
    "P2_bracket": prq(53, "[21.0, 105.3] µm"),
    "P2_torsion": prq(55, "s ≥ 1.5 × 30 µm"),
    "P2_survives": prq(56, "Only the lower part of the bracket survives"),
    "P2_falsified_if": prq(58, "every range ≥ **21.0 µm** at 95% CL"),
    "P3_mu_sym": prq(64, "mu_sym"),
    "P3_c4": prq(65, "c4_pta_product"),
    "P3_gmu": prq(66, "planckGmuBound"),
    "not_derived": prq(80, "None is derived from the K3×T² mathematics"),
}
# threshold re-computation (the rule reads the file's numbers; scipy is the cross-check)
thr = {
    "1dof_3sigma": float(stats.chi2.isf(2 * stats.norm.sf(3), 1)),
    "1dof_5sigma": float(stats.chi2.isf(2 * stats.norm.sf(5), 1)),
    "2dof_3sigma": float(stats.chi2.isf(2 * stats.norm.sf(3), 2)),
    "2dof_5sigma": float(stats.chi2.isf(2 * stats.norm.sf(5), 2)),
}
FILE_THR = {"1dof_3sigma": 9.00, "1dof_5sigma": 25.00, "2dof_3sigma": 11.83, "2dof_5sigma": 28.74}
for k in thr:
    assert abs(thr[k] - FILE_THR[k]) < 0.01, (k, thr[k])


def sigma_equiv(dchi2, dof):
    return float(stats.norm.isf(stats.chi2.sf(dchi2, dof) / 2))


# ---------------------------------------------------------------- LeanMaster quotes (committed tree)
lm_cache = {}


def lmq(path, lineno, must, span=1):
    if path not in lm_cache:
        lm_cache[path] = git(LM, "show", f"{LEANMASTER_COMMIT}:{path}").splitlines()
    lines = lm_cache[path][lineno - 1: lineno - 1 + span]
    joined = " ".join(l.strip() for l in lines)
    assert must in joined, (path, lineno, must, joined)
    return {"file_line": f"{path}:{lineno}" + (f"-{lineno + span - 1}" if span > 1 else ""),
            "verbatim": joined, "commit": LEANMASTER_COMMIT}


S6 = "docs/STREAM6_EXPERIMENT_PLAN.md"
S7 = "docs/STREAM7_HYPOTHESIS_INVENTORY.md"
S8 = "docs/STREAM8_WHICH_K3.md"
LMQ = {
    "P1_excluded": lmq(S6, 74, "By the pre-registered rule, **P1 is excluded**", 4),
    "P62_kappa": lmq(S6, 80, "(the Tier C dual-scale identification of Stream 3)", 4),
    "S6_conclusion": lmq(S6, 87, "does not survive the existing data.", 5),
    "S6_new_hypothesis": lmq(S6, 92, "A new hypothesis would need a new, independently derived observable", 2),
    "CA_excluded": lmq(S7, 83, "**C-A is excluded**", 1),
    "CA_S6_S7_closed": lmq(S7, 85, "both the extra-dimension reading and the Hubble-scale holographic reading", 2),
    "CB_not_testable": lmq(S7, 94, "not falsifiable", 2),
    "S7_where_leaves": lmq(S7, 104, "Every observable the programme can derive without an unconstructed compactification", 4),
    "S7_new_hypothesis": lmq(S7, 108, "The common cause is the dual-scale identification", 4),
    "S7_N4": lmq(S7, 65, "K3 × T² itself gives `N = 4` in four dimensions", 1),
    "S8_forced_not_fitted": lmq(S8, 15, "A choice must be forced, not fitted.", 2),
    "S8_E2_trapping": lmq(S8, 101, "under trapping, T² comes to rest at `(ω, ω)`", 1),
    "S8_E2_observables_none": lmq(S8, 112, "Observables: none", 2),
    "S8_octad_stabilizer_M24": lmq(S8, 134, "`(ℤ₂)⁴ ⋊ A₈`", 2),
    "S8_TW_overarching_A7": lmq(S8, 136, "`(ℤ₂)⁴ ⋊ A₇` (order 40320, a maximal subgroup of `M₂₃`)", 1),
    "S8_E3_observables_none": lmq(S8, 153, "Observables: none (`N = 4`, non-chiral)", 1),
    "S8_synthesis_no_observable": lmq(S8, 272, "No observable follows", 2),
    "S8_nonsymplectic_counterweight": lmq(S8, 318, "Non-symplectic symmetry points the other way", 1),
    "S8_P82_not_vacuum": lmq(S8, 369, "Not for the vacuum.", 1),
    "S8_P82_observables_none": lmq(S8, 374, "Observables: none (`N = 4`)", 1),
    "omegaLambda": lmq("DualScaleCosmology/DarkEnergyScale.lean", 83, "astropy.cosmology.Planck18", 2),
    "planckGmuBound": lmq("DualScaleCosmology/CosmicString.lean", 83, "planckGmuBound : ℝ := 1.5e-7", 1),
    "tadpole_budget_TT": lmq("DualScaleStream2/Flux/Tadpole.lean", 116, "Tripathy–Trivedi eq. (2.3)", 1),
    "tadpole_not_construction": lmq("DualScaleStream2/Flux/Tadpole.lean", 40, "nothing about whether a flux configuration", 3),
}
# theorem presence at the committed tree (names only; no kernel/axiom audit run here)
THEOREMS = {}
for name, path in [
    ("p1_verdict_excluded", "DualScaleCosmology/Stream6Verdict.lean"),
    ("p62_selfdual_fixed_point", "DualScaleCosmology/Stream6Verdict.lean"),
    ("p62_towers_coincide", "DualScaleCosmology/Stream6Verdict.lean"),
    ("ca_verdict_excluded", "DualScaleCosmology/Stream7CA.lean"),
    ("cb_not_testable", "DualScaleCosmology/Stream7CB.lean"),
    ("roots_ww", "DualScaleDyons/SelfDualT2.lean"),
    ("ww_root_lattice_is_A2", "DualScaleDyons/SelfDualT2.lean"),
    ("kummer_roots", "DualScaleDyons/KummerE3.lean"),
    ("octad_is_kummer", "DualScaleDyons/KummerE3.lean"),
    ("order14_not_geometric", "DualScaleDyons/ForgerE4.lean"),
    ("smallest_black_hole", "DualScaleDyons/AttractorCharges.lean"),
    ("tadpole_budget", "DualScaleStream2/Flux/Tadpole.lean"),
]:
    src = git(LM, "show", f"{LEANMASTER_COMMIT}:{path}").splitlines()
    hits = [i + 1 for i, l in enumerate(src) if l.startswith(f"theorem {name}")]
    assert hits, name
    THEOREMS[name] = f"{path}:{hits[0]}"
# P8.4c (K3Enhancement.lean) status at the time of this build: committed or not?
k3e_tracked = git(LM, "ls-files", "DualScaleDyons/K3Enhancement.lean").strip() != ""
lm_status = git(LM, "status", "--porcelain").strip().splitlines()

# ---------------------------------------------------------------- experiment numbers
e1 = load("E1-desi-dr2/e1_desi_dr2_report.json")
e3 = load("E3-M0/m0_result.json")
e2 = load("e2_pta_and_tda_extension/e2_pta_result.json")
chk_chi2 = load("skeptic_statistics/check_e1_e3_chi2.json")
chk_slope = load("skeptic_statistics/check_e2_slope.json")
e5b = load("E5-cosmic-web-tda-scaled/e5b_lognormal_and_poisson_null_report.json")
e5b_r = load("E5-cosmic-web-tda-scaled/e5b_r_range_decomposition.json")
gate_pw = load("skeptic_statistics/check_e5b_gate_power.json")
psel = load("skeptic_statistics/check_e5b_poisson_selection.json")
kac = load("skeptic_statistics/check_e5_known_answer_controls.json")
cmb = load("E5-cmb-tda/e5_cmb_tda_report.json")
cmb_spec = load("skeptic_statistics/check_cmb_null_spectrum.json")
cmb_rule = load("skeptic_statistics/decision_rule_cmb_recal.json")
e4 = load("E4-cosmic-web-tda/e4_cosmic_web_tda_report.json")

dr2 = e1["DR2"]
dchi2_dr2 = dr2["delta_chi2_CPL_vs_Lambda_2dof"]
dchi2_dr1 = e1["DR1"]["delta_chi2_CPL_vs_Lambda_2dof"]
dchi2_m0 = e3["delta_chi2_M0_minus_LCDM_fitted"]

nt = cmb["null_test"]
pvals = {f"{lvl}_{s}": nt[f"statistic_{lvl}"][s]["p_value_chi2_survival"]
         for lvl in ("sublevel", "superlevel") for s in ("b0", "b1", "euler_chi")}
pmin_key = min(pvals, key=pvals.get)
pmin = pvals[pmin_key]
gmu_scan = cmb["cosmic_string_sensitivity"]["gmu_scan_sorted_ascending"]
gmu_frac = {k: v["fraction_exceeding_null_95th_pct_LOO"] for k, v in gmu_scan.items()}
gmu_clear = sorted(float(k) for k, v in gmu_frac.items() if v >= 0.95)
gmu_below = sorted(float(k) for k, v in gmu_frac.items() if v < 0.95)
gmu_lo = max(x for x in gmu_below if x < gmu_clear[0])
planck_gmu = 1.5e-7  # read from LMQ["planckGmuBound"] line; asserted above
assert "1.5e-7" in LMQ["planckGmuBound"]["verbatim"]

flux = json.loads(git(WT, "show", f"{FLUX_COMMIT}:audit/k3t2_rigidity_v2/E-flux/flux_enumeration_results.json"))["results"]
flux_counts = {k: {n: v["with_tadpole_bound_alpha_x_sq_le_24"][n]["total_valid_ordered_pairs"]
                   for n in v["with_tadpole_bound_alpha_x_sq_le_24"]} for k, v in flux.items()}
flux_counts_no_tadpole = {k: {n: v["without_tadpole_bound_negative_control"][n]["total_valid_ordered_pairs"]
                              for n in v["without_tadpole_bound_negative_control"]} for k, v in flux.items()}
# finite list from tadpole (1/2) N_flux + N_D3 = 24, N_flux = 2 alpha^2, 8 | alpha^2 (Track E, tier B)
flux_list = [{"alpha_x_sq": a, "N_flux": 2 * a, "N_D3": 24 - a} for a in range(8, 25, 8)]

ledger = {
    "generated_by": "audit/reverse_zero/synthesis/build_ledger.py",
    "header": "Generated by workflow reverse-to-zero; exploratory (tier X) unless stated; M0 is a hypothesis change, not a derivation.",
    "worktree": str(WT), "branch": "loop/reverse-zero",
    "leanmaster_commit_read": LEANMASTER_COMMIT,
    "leanmaster_working_tree_status_at_build": lm_status,
    "leanmaster_P84c_K3Enhancement_committed": k3e_tracked,
    "prereg_quotes": PRQ,
    "prereg_thresholds_recomputed_scipy": thr,
    "leanmaster_quotes_verbatim": LMQ,
    "leanmaster_theorem_locations_at_commit": THEOREMS,
    "parameter_count": {
        "M0_free_theory_parameters": e3["M0_result"]["k"],
        "M0_profiled_nuisances": ["BAO r_d*h scale (u)", "SN absolute offset (H0/M_B)"],
        "M2_free_theory_parameters": e3["M2_result"]["k"],
        "M2_parameters": ["mu_sym", "c4_pta_product"],
        "status_of_mu_sym_c4": "untested by any dataset in hand: not bounded, not derived, not zero",
        "mechanism": "hypothesis change (symmetron sector deleted; c4_pta_product set to GR value 0; Omega_Lambda frozen at imported Planck18 value)",
        "tier": "X",
    },
    "E1_P1_pipeline_action": {
        "tier": "X",
        "DR2_delta_chi2_CPL_vs_LCDM_2dof": dchi2_dr2,
        "DR2_sigma_equiv": sigma_equiv(dchi2_dr2, 2),
        "DR1_delta_chi2_CPL_vs_LCDM_2dof_fullcov": dchi2_dr1,
        "DR1_sigma_equiv": sigma_equiv(dchi2_dr1, 2),
        "DR1_diagcov_skeptic": chk_chi2["DR1_E1cut_1590_diagcov"]["dchi2_cpl_vs_lcdm_2dof"],
        "DR2_diagcov_skeptic": chk_chi2["DR2_E1cut_1590_diagcov"]["dchi2_cpl_vs_lcdm_2dof"],
        "DR2_1580_cut_skeptic": chk_chi2["DR2_E3cut_1580_fullcov"]["dchi2_cpl_vs_lcdm_2dof"],
        "DR2_w0_wa": [dr2["cpl_om_free"].get("w0"), dr2["cpl_om_free"].get("wa")] if isinstance(dr2.get("cpl_om_free"), dict) else None,
        "DR2_fitted_Omega_L": dr2["fitted_lcdm"]["Om_L"],
        "DR2_Omega_L_distance_sigma": dr2["omega_lambda_distance_from_0.68885_in_sigma"],
        "negative_control_Om_L_0.5_DR2": dr2["negative_control_Om_L_0.5"],
        "falsification_rule_applicable": False,
        "why": "PRE_REGISTRATION.md:40-43 names DESI DR3/final BAO + CMB + SN, Euclid, Rubin. DR2 BAO + SN without CMB is covered only by the pipeline action (46-48).",
        "note_script": "E1 used a new script, not decisive_experiment.py named at PRE_REGISTRATION.md:47",
    },
    "E3_M0_vs_fitted_LCDM": {
        "tier": "X",
        "delta_chi2_1dof": dchi2_m0,
        "sigma_equiv": sigma_equiv(dchi2_m0, 1),
        "chi2_sf": float(stats.chi2.sf(dchi2_m0, 1)),
        "delta_AIC": e3["delta_AIC_M0_minus_LCDM_fitted"],
        "delta_BIC": e3["delta_BIC_M0_minus_LCDM_fitted"],
        "delta_AIC_M2_minus_M0": e3["delta_AIC_M2_minus_M0"],
        "delta_BIC_M2_minus_M0": e3["delta_BIC_M2_minus_M0"],
        "M2_chi2_equals_M0": e3["M2_result"]["chi2"] == e3["M0_result"]["chi2"],
        "variant_1701_delta_chi2": e3["variant_full_1701_including_calibrators"]["delta_chi2_M0_minus_LCDM_fitted"],
        "negative_controls_delta_chi2": {k: v["delta_vs_fitted"] for k, v in e3["negative_controls"].items()},
        "preregistered": False,
        "note": "9.00 comes from the threshold table (PRE_REGISTRATION.md:35); M0 vs fitted LCDM is not itself a registered test (R2 precedent: 'reported, not pre-registered'). AIC/BIC penalty on M2 is for parameters that enter no likelihood term: not evidence.",
    },
    "E2_PTA": {
        "tier": "B (exact linear scan) for slope; data branch NO_DATA",
        "branch": e2["branch"],
        "slope_true": chk_slope["slope_recomputed"],
        "slope_reported_wrong": 6.171,
        "hd_dynamic_range": chk_slope["hd_dynamic_range"],
        "c4_threshold_1pct": chk_slope["threshold_1pct_recomputed"],
        "c4_thresholds_1_5_10pct": [chk_slope["hd_dynamic_range"] * f for f in (0.01, 0.05, 0.10)],
        "parameter_effect": "none",
    },
    "E4_cosmic_web_N400": {"tier": "X", "n_target": e4["n_target"], "note": "real vs z-shuffle Betti tie at N=400; superseded by E5b"},
    "E5b_cosmic_web": {
        "tier": "X",
        "void_control_pass": kac["committed_gate"]["pass"],
        "void_control_n": [kac["committed_gate"]["n_bars_with_death_within_20pct_of_void_radius"], kac["committed_gate"]["n_voids_planted"]],
        "circle_control_ratio": kac["committed_circle_ratio"],
        "b_fit": e5b["xi_r_and_bias"]["b_fit"],
        "gate_rms_z": e5b["lognormal_gate"]["rms_z"],
        "gate_threshold_source": gate_pw["gate_threshold_from_source"],
        "gate_rms_z_for_poisson_xi0": gate_pw["rms_z_poisson_xi0"],
        "gate_passing_scales_of_mock_xi": gate_pw["passing_scale_range"],
        "chi2_diag_mock_vs_poisson_20bins": [gate_pw["chi2_diag_committed_mock"], gate_pw["chi2_diag_poisson"]],
        "n_mocks": e5b["lognormal_mocks_iii"]["n_succeeded"],
        "rank_p_floor": 1.0 / (e5b["lognormal_mocks_iii"]["n_succeeded"] + 1),
        "betti_r_2_20": {h: e5b_r[h]["unvalidated_2_20"] for h in ("H0", "H1", "H2")},
        "betti_r_20_40": {h: e5b_r[h]["gate_validated_20_40"] for h in ("H0", "H1", "H2")},
        "poisson_shell_counts": psel["shell_counts_poisson"],
        "data_shell_counts": psel["shell_counts_data"],
        "poisson_ks_r": psel["ks_r"],
        "nside64_pixel_mpc_h": psel["nside64_pixel_comoving_mpc_h_at_r_lo_r_hi"],
        "verdict": "INCONCLUSIVE; xi(r) gate has no power against an unclustered catalogue, so 'gate passes' is not validation",
    },
    "E5_cmb": {
        "tier": "X",
        "n_sims": nt["n_sims_total"],
        "calibration_ratio_variance": nt["calibration_check"]["ratio_sim_over_data"],
        "p_values_raw": pvals,
        "smallest": {"stat": pmin_key, "p": pmin, "bonferroni_6": min(1.0, 6 * pmin),
                     "bonferroni_2_effective": min(1.0, 2 * pmin)},
        "null_spectrum_bands_data_over_sim": cmb_spec["stepA_pipeline_null"]["bands"],
        "recalibration_gate_passed": cmb_spec["stepB_recalibrated_spectrum_check"]["gate_all_bands_abs_z_lt_3"],
        "recalibration_last_iteration_max_abs_z": float(max(abs(b["z"]) for b in cmb_spec["stepB_recalibrated_spectrum_check"]["iterations"][-1]["bands"])),
        "decision_rule_written_before_results": cmb_rule["written_before_results"],
        "gmu_frac_exceeding": gmu_frac,
        "gmu_95pct_bracket": [gmu_lo, gmu_clear[0]],
        "gmu_bracket_over_planck_bound": [gmu_lo / planck_gmu, gmu_clear[0] / planck_gmu],
        "verdict": "Gaussianity of M0 not rejected by this pipeline, but null sims do not match the data spectrum at ell>150; the non-rejection is not established at the claimed strength (a spectrum-matched null is pending). No topology signal survives controls.",
    },
    "flux_vacua_external": {
        "tier": "B inside stated truncation + L (TT identification); external input",
        "source": f"commit {FLUX_COMMIT} audit/k3t2_rigidity_v2/E-flux/flux_enumeration_results.json (branch loop/k3t2-rigidity; read via git show)",
        "finite_list": flux_list,
        "pair_counts_with_tadpole": flux_counts,
        "pair_counts_without_tadpole": flux_counts_no_tadpole,
        "caveats": ["counts grow with the coefficient box N: finiteness comes from the box, not the tadpole",
                    "index-8 (rank 6) / index-4 (rank 4) subfamily of TT (2.5) quantisation",
                    "TT section 4.1 example family only",
                    "orbit counts are not certified upper bounds",
                    "20 K3 Kahler moduli + T2 Kahler modulus remain free (TT ll. 843-845)"],
    },
    "skeptic_statistics_verdicts": load("skeptic_statistics/reruns_and_verdicts.json")["verdicts"],
    "framing_corrections_applied_in_report": [
        "P1: DR2 BAO+SN without CMB is the pipeline action (PRE_REGISTRATION.md:46-48), not a falsification test; P1 has no 3-sigma bar",
        "M0 vs fitted LCDM: 9.00 is the threshold table, not a registered rule for this comparison",
        "Delta chi2 = 0.746 (1 dof) is 0.864 sigma, not ~0.1 sigma",
        "M24 octad stabilizer is (Z2)^4 x| A8; (Z2)^4 x| A7 (order 40320) is TW's overarching group, maximal in M23",
        "No TDA prediction is attributed to LeanMaster: STREAM8 says 'Observables: none' (lines 112, 153, 374 at eb791e7)",
        "LeanMaster Stream 6 P1 (extra-dimension radius) is not PRE_REGISTRATION P1 (cosmological constant)",
        "P6.2 kappa=1 quoted with its Tier C qualifier (alpha' = s^2)",
        "Omega_Lambda = 0.68885 is an imported Planck18 value (DarkEnergyScale.lean:83-84), tier L, not derived",
        "LeanMaster Streams 6-8 never mention mu_sym or c4_pta_product; their untested status is this project's statement",
        "xi(r) gate pass is non-rejection with no power, not confirmation; GATE_Z_THRESHOLD was fixed in source, not pre-registered",
        "One parameter count: M0 0 theory parameters by hypothesis change (+2 profiled nuisances); M2 2, untested",
    ],
    "absent": [
        {"id": "nanograv_15yr_hd_angular_correlation", "status": "ABSENT",
         "tried": ["nanograv/15yr_stochastic_analysis figure_1", "PTArcade", "Zenodo searches"],
         "recheck_this_session": "/mnt/disks/disk-socrateai-local-1/NANOGrav15yr/NANOGrav15yr_PulsarTiming_v2.1.0/correlations/narrowband/README.nb_correlations: per-pulsar timing-parameter correlations, not an inter-pulsar angular table; raw TOAs present, so HD is reconstructible by a scoped enterprise run (not done)"},
        {"id": "fifth_force_eotwash_2002_11761", "status": "ABSENT", "tried": ["https://arxiv.org/e-print/2002.11761"],
         "note": "source tarball has only .tex + 9 figure PDFs, no table"},
    ],
    "not_attempted": ["CMB likelihood for E1 (so DESI DR2+CMB+SN 3.1 sigma not reproduced)",
                      "literature galaxy-bias table for E5b (b_fit=0.893 not cross-checked)",
                      "kernel/axiom audit of the cited LeanMaster theorems (existence at commit only)"],
}

OUT.write_text(json.dumps(ledger, indent=1, ensure_ascii=False, default=float) + "\n")
print(json.dumps({
    "out": str(OUT),
    "DR2_dchi2": dchi2_dr2, "DR2_sigma": ledger["E1_P1_pipeline_action"]["DR2_sigma_equiv"],
    "M0_dchi2": dchi2_m0, "M0_sigma": ledger["E3_M0_vs_fitted_LCDM"]["sigma_equiv"],
    "cmb_pmin": [pmin_key, pmin, 6 * pmin, 2 * pmin],
    "gmu_bracket": [gmu_lo, gmu_clear[0]], "gmu_ratio": ledger["E5_cmb"]["gmu_bracket_over_planck_bound"],
    "flux_counts": flux_counts, "K3Enhancement_committed": k3e_tracked,
}, indent=1, default=float))
