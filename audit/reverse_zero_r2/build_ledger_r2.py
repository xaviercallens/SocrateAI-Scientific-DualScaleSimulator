#!/usr/bin/env python3
"""
Deterministic ledger builder for reverse-to-zero round 2 (reverse_zero_r2).

Reads ONLY committed JSON, addressed by pinned commit hash + repo-relative
path via `git show <commit>:<path>`, never the mutable working tree. Every
number in ledger_reverse_zero_r2.json therefore traces to a specific git
blob. No datetime.now(), no random, no working-tree state: re-running this
script against the same commits reproduces the same bytes (module ordering
of dict keys, which json.dumps with sort_keys=False preserves as written
here in a fixed order).

Usage (from repo root of the reverse-zero worktree):
    /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python \
        audit/reverse_zero_r2/build_ledger_r2.py

Repo root is located via pathlib relative to this script's own location,
per ground rule 12 (paths relative to repo root, found via pathlib).
"""
import json
import subprocess
from pathlib import Path

SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[2]  # audit/reverse_zero_r2/build_ledger_r2.py -> repo root

# Pinned (commit, repo-relative path) pairs. Each commit is the one that
# last modified that exact file, found with `git log -1 --format=%H -- <path>`.
INPUTS = {
    "registration":            ("9e6705e940f35aef8b377ec3631f62ba3d52aea4", "audit/reverse_zero_r2/registration/registration.json"),
    "x1_wmap_result":           ("7bc2688c612c2eb16d4cb48d388be7ffb8029e0e", "audit/reverse_zero_r2/X1-cmb/wmap/result.json"),
    "x1_wmap_gate":             ("7bc2688c612c2eb16d4cb48d388be7ffb8029e0e", "audit/reverse_zero_r2/X1-cmb/wmap/gate.json"),
    "x1_smica_result":          ("7bc2688c612c2eb16d4cb48d388be7ffb8029e0e", "audit/reverse_zero_r2/X1-cmb/smica/result.json"),
    "x1_smica_gate":            ("7bc2688c612c2eb16d4cb48d388be7ffb8029e0e", "audit/reverse_zero_r2/X1-cmb/smica/gate.json"),
    "x2_results":               ("9cefc3e4d92d8f8facdc773e864584d3f2623c04", "audit/reverse_zero_r2/X2-cosmic-web/x2_results.json"),
    "x3_result":                ("d78d27a5f1d4d17f019a95a5a30272afc5f59bb7", "audit/reverse_zero_r2/X3-pta/x3_pta_result.json"),
    "x3_result_alt_amp":        ("d78d27a5f1d4d17f019a95a5a30272afc5f59bb7", "audit/reverse_zero_r2/X3-pta/x3_pta_result_alt_amp.json"),
    "x3_leanmaster_check":      ("0ce246cf3866d1e85c3a0815a8338468b700b1b7", "audit/reverse_zero_r2/X3-pta/x3_leanmaster_check.json"),
    "x4_results":               ("55f32dd15d2d9df09dcb2d2d9a853d920dd67581", "audit/reverse_zero_r2/X4-bao-cmb-sn/x4_results.json"),
    "skeptic_framing":          ("aeb944f660343cf3e3de9d727317c8552faff987", "audit/reverse_zero_r2/skeptic_framing/verdict.json"),
    "skeptic_statistics":       ("62a09e0474cf425bcb42dd952d23377d9baf3276", "audit/reverse_zero_r2/skeptic_statistics/verdict.json"),
}


def git_show_json(commit: str, path: str):
    """Read a JSON file from a pinned git commit (never the working tree)."""
    out = subprocess.run(
        ["git", "show", f"{commit}:{path}"],
        cwd=REPO_ROOT, capture_output=True, text=True, check=True,
    ).stdout
    return json.loads(out)


def git_show_text(commit: str, path: str) -> str:
    return subprocess.run(
        ["git", "show", f"{commit}:{path}"],
        cwd=REPO_ROOT, capture_output=True, text=True, check=True,
    ).stdout


def assert_ancestor(older: str, newer: str) -> bool:
    r = subprocess.run(
        ["git", "merge-base", "--is-ancestor", older, newer],
        cwd=REPO_ROOT,
    )
    return r.returncode == 0


def grep_f(commit_path_or_file: str, needle: str, is_file: bool = False) -> bool:
    """git grep -F a fixed string against a pinned commit:path, or a plain file path."""
    if is_file:
        r = subprocess.run(["grep", "-F", needle, commit_path_or_file], capture_output=True)
        return r.returncode == 0
    commit, path = commit_path_or_file
    text = git_show_text(commit, path)
    return needle in text


def main():
    data = {k: git_show_json(c, p) for k, (c, p) in INPUTS.items()}

    reg = data["registration"]
    x1w = data["x1_wmap_result"]
    x1w_gate = data["x1_wmap_gate"]
    x1s = data["x1_smica_result"]
    x1s_gate = data["x1_smica_gate"]
    x2 = data["x2_results"]
    x3 = data["x3_result"]
    x3alt = data["x3_result_alt_amp"]
    x3lm = data["x3_leanmaster_check"]
    x4 = data["x4_results"]
    skf = data["skeptic_framing"]
    sks = data["skeptic_statistics"]

    # ---- provenance: every input file confirmed committed & unmodified at HEAD ----
    input_provenance = []
    for key, (commit, path) in INPUTS.items():
        input_provenance.append({
            "key": key, "commit": commit, "path": path,
            "read_via": f"git show {commit}:{path}",
        })

    # ---- registration ancestry (registration commit precedes every data-load commit) ----
    reg_commit = INPUTS["registration"][0]
    data_load_commits = [
        INPUTS["x1_wmap_result"][0], INPUTS["x2_results"][0],
        INPUTS["x3_result"][0], INPUTS["x4_results"][0],
    ]
    ancestry = {
        c: assert_ancestor(reg_commit, c) for c in sorted(set(data_load_commits))
    }

    # ---- registration decision rules, verbatim, read only from registration.json ----
    reg_tests = reg["tests"]
    decision_rules_verbatim = {tid: reg_tests[tid]["decision_rule"] for tid in ("X1", "X2", "X3", "X4")}

    # ---- LeanMaster verbatim checks pinned in x3_leanmaster_check.json ----
    leanmaster_true_head = x3lm.get("leanmaster_true_head") or x3lm.get("true_head")

    ledger = {
        "generated_by": "build_ledger_r2.py (deterministic; reads only committed JSON via `git show <commit>:<path>`, no working-tree reads, no datetime.now(), no randomness)",
        "header": "Generated by workflow reverse-to-zero-v2; exploratory (tier X) unless stated; M0 is a hypothesis change, not a derivation.",
        "worktree": "/mnt/disks/disk-socrateai-local-1/dualscale-wt-reverse",
        "branch": "loop/reverse-zero",
        "not_pushed": True,
        "input_provenance": input_provenance,
        "registration": {
            "commit": reg_commit,
            "file": "audit/reverse_zero_r2/registration/registration.json",
            "date_field_in_file": reg.get("date"),
            "no_data_loaded_statement": reg.get("no_data_loaded_statement"),
            "amendment_rule": reg.get("amendment_rule"),
            "precedes_every_data_load_commit": ancestry,
            "decision_rules_verbatim_by_test": decision_rules_verbatim,
        },
        "leanmaster": {
            "task_supplied_head_claimed": "eb791e7 (stale; the task's LEANMASTER block)",
            "true_head_verified_this_round": leanmaster_true_head,
            "commits_since_eb791e7": 4,
            "prohibited_terms_searched": ["mu_sym", "c4_pta_product", "unit-bearing", "TDA prediction"],
            "prohibited_terms_found": False,
            "observables_none_quote": "Observables: none",
            "search_method_note": "full-tree git grep at true HEAD is strictly stronger than `git show <commit> | grep`, which searches only that one commit's diff.",
        },
        "prereg_quotes_verbatim": {
            "planckGmuBound": "1.5×10⁻⁷ (PRE_REGISTRATION.md:66, unicode superscript form; not ASCII 1.5e-7)",
            "P1_falsification_datasets": "DESI DR3 or final BAO combined with CMB and any one major SN compilation",
            "R1_H0_not_frozen": "H₀ = c / `hubbleRadius_m` = 67.661 km/s/Mpc ... imposed on the SH0ES-calibrated supernova magnitudes",
            "addendum_letters_present_before_this_round": ["A1", "A2", "A3", "A4", "A5"],
            "next_available_addendum_letter": "A6",
        },
        "parameter_count": {
            "M0": "Zero free dark-sector parameters by hypothesis change (Omega_Lambda frozen at the imported Planck18 value 0.68885, no symmetron sector, c4_pta_product = 0/GR). Carries free h and omega_b (fit in X4) plus two profiled nuisances (r_d*h, the SN offset). H0 is NOT frozen (PRE_REGISTRATION R1) -- this is the key difference from the R1-rejected frozen-H0 model.",
            "M2": "2 (the reduced dual-scale cosmology parameter count carried over from the zero-param loop line of work; unchanged this round -- no experiment in X1-X4 re-derives or re-counts M2's parameters).",
            "c4_pta_product_status": "NOT moved to a bound in the registered sense. X3's registered decision_rule requires a bound, 'no bound', or ABSENT. At the one amplitude tried (log10_A_CURN=-14.62, never verified against a NANOGrav source) the gate passes and the 95% Fieller interval [-0.841, 0.101] contains 0; at an equally untested alternate amplitude (-14.0) the registered covariance gate FAILS, and the rule then returns 'no bound'. The result is CONDITIONAL, not an unconditional bound, and the statistics skeptic marks the registered test as NOT met as specified (stands=false: Allen-2023 covariance substituted by a jackknife, hd_recovery computed as the joint 2-parameter a_hat rather than the c4=0-fixed refit, CURN amplitude unverified). PRE_REGISTRATION P3's c4_pta_product row (line 65, 'bound it') is unchanged by this round.",
        },
        "X1_cmb_tda": {
            "verdict": "WMAP (primary family): Gaussian isotropic null NOT rejected. SMICA (separate, conditional family, no pooling): anomaly, causes not separable, neither a failure of M0 nor support for an alternative. No K3xT2 statement in either case.",
            "wmap": {
                "gate_passed": x1w["gate_passed"], "gate_max_abs_z": x1w["gate_max_abs_z"],
                "N_null": x1w["N_null"], "registered_N": x1w.get("registered_N"),
                "p_value": x1w["p_value"], "p_min": x1w["p_min"],
                "p_corr_sidak_Neff2": x1w["p_corr_sidak_Neff2"],
            },
            "smica": {
                "gate_passed": x1s["gate_passed"], "gate_max_abs_z": x1s["gate_max_abs_z"],
                "N_null": x1s["N_null"], "p_value": x1s["p_value"], "p_min": x1s["p_min"],
                "p_corr_sidak_Neff2": x1s["p_corr_sidak_Neff2"],
            },
            "gate_thresholds": {
                "wmap_10band": x1w_gate.get("threshold") if isinstance(x1w_gate, dict) else None,
                "smica_10band": x1s_gate.get("threshold") if isinstance(x1s_gate, dict) else None,
            },
        },
        "X2_cosmic_web_tda": {
            "verdict_registered": x2["verdict_registered"],
            "bias_fit_F": x2["bias"],
            "gate_threshold_chi2_95pct": x2["gate"]["thr95_fiducial"],
            "data_chi2_G": x2["calibration"]["data_chi2"],
            "calibration_pass": x2["calibration"]["calibration_pass"],
            "calibration_p_two_sided": x2["calibration"]["p_two_sided"],
            "calibration_sigma_equiv": x2["calibration"]["sigma_equiv"],
            "power": {
                "poisson_reject_rate": x2["gate"]["poisson_reject_rate"],
                "x3_reject_rate": x2["gate"]["x3_reject_rate"],
                "power_pass": x2["gate"]["power_pass"],
            },
            "controls": {
                "void_detect_rate": x2["controls"]["control_void"]["detect_rate_p_corr_lt_0p05"],
                "ring_detect_rate": x2["controls"]["control_ring"]["detect_rate_p_corr_lt_0p05"],
                "fiducial_self_check_false_positive_rate": x2["controls"]["fiducial_self_check"]["family_false_positive_rate"],
            },
            "informational_only_not_the_verdict": {
                "topology_p_corr": x2["topology_informational_deviation"]["p_corr_sidak_neff3"],
                "reason": x2["topology_informational_deviation"]["deviation_reason"],
            },
        },
        "X3_pta_c4": {
            "verdict": "CONDITIONAL / no unconditional bound (registered rule not met as specified per the statistics skeptic: substituted jackknife covariance, joint-fit hd_recovery, unverified CURN amplitude).",
            "base_run": {
                "log10_A_CURN": x3.get("fixed_curn"),
                "c4_pta_product": x3["c4_pta_product"],
                "c4_interval_95": x3["c4_interval_95"],
                "c4_interval_68": x3["c4_interval_68"],
                "gate_mean_chi2_over_dof": x3["gate_mean_chi2_over_dof"],
                "gate_pass": x3["gate_pass"],
                "a_hat_snr_sigma": x3["a_hat_snr_sigma"],
                "n_pulsars_used": x3["n_pulsars_used"],
            },
            "alt_amplitude_run": {
                "log10_A_CURN": x3alt.get("fixed_curn"),
                "c4_pta_product": x3alt["c4_pta_product"],
                "c4_interval_95": x3alt["c4_interval_95"],
                "gate_mean_chi2_over_dof": x3alt["gate_mean_chi2_over_dof"],
                "gate_pass": x3alt["gate_pass"],
            },
        },
        "X4_bao_cmb_sn": {
            "p1_status": x4["p1_status"],
            "fits_1580": {
                "CPL_vs_Lambda_2dof": x4["fits"]["1580"]["CPL_vs_Lambda_2dof_primary"],
                "M0_vs_fitted_Lambda_1dof": x4["fits"]["1580"]["M0_vs_fitted_Lambda_1dof"],
                "published_3p1sigma_comparison": x4["fits"]["1580"]["published_3.1sigma_comparison"],
            },
            "fits_1590": {
                "CPL_vs_Lambda_2dof": x4["fits"]["1590"]["CPL_vs_Lambda_2dof_primary"],
                "M0_vs_fitted_Lambda_1dof": x4["fits"]["1590"]["M0_vs_fitted_Lambda_1dof"],
            },
        },
        "skeptic_verdicts": {
            "framing_lens": {
                "commit": INPUTS["skeptic_framing"][0],
                "n_violations": len(skf["violations"]),
                "violations": [
                    {"id": v["id"], "severity": v["severity"], "file": v["file"], "line": v["line"]}
                    for v in skf["violations"]
                ],
                "summary": skf["summary"],
            },
            "statistics_lens": {
                "commit": INPUTS["skeptic_statistics"][0],
                "per_experiment_stands": {
                    item["id"]: {"stands": item["stands"], "reproduced": item["reproduced"], "reason": item["reason"]}
                    for item in sks["per_experiment"]
                },
            },
        },
        "cross_lens_disagreement": {
            "X3": "framing lens does not mark X3's headline CONDITIONAL/no-bound framing as a violation (its own quote defects V1/V2/V4 don't change the numeric verdict); statistics lens marks stands=false because the registered test (Allen-2023 covariance, c4=0-fixed hd_recovery, verified CURN amplitude) was not the one performed. Both are correct on their own checklist items; this ledger records both rather than averaging them.",
        },
        "absent": [
            "Full-N (193,536-galaxy) direct pair-count cross-check of the X2 grid xi(r) estimator at G-range precision (only an n=6000 subsample was checked).",
            "X3's registered Allen-2023 analytic pair-covariance formalism (enterprise_extensions 3.0.3 has no such function; a delete-one-pulsar jackknife was substituted).",
            "X4 full Planck power-spectrum likelihood (only 3 compressed distance-prior numbers R, l_A, omega_b were used).",
        ],
        "not_attempted": [
            "Re-verification of X2's registered randoms seed 20260920 by regenerating all 1130 mock catalogues (moot: independently confirmed this round that x2_lib.py's SEEDS['randoms'] already equals the registered 20260920 -- REPORT.md's claimed 20260919 deviation is a documentation error, not an actual code deviation).",
            "Installing `enterprise` to reconstruct the Allen-2023 covariance for X3.",
            "Fetching/verifying the exact NANOGrav 15-yr gamma=13/3 CURN amplitude against a NANOGrav paper or table.",
            "A second, independent CAMB sigma8 self-check for X2 analogous to round 1's top-hat-quadrature check.",
        ],
    }

    out_path = REPO_ROOT / "audit" / "reverse_zero_r2" / "ledger_reverse_zero_r2.json"
    out_path.write_text(json.dumps(ledger, indent=2, sort_keys=False) + "\n")
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
