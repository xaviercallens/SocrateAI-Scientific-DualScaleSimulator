"""Backfill source (a2): the fixed-library validation.

Source worktree : dualscale-wt-tdasimple, branch loop/tda-simple
Source files    : audit/tda_validation/tda_fixed/{fix_expectations.json,
                  validation_results.json}
Pre-registration: fix_expectations.json, "written before any line of
                  tda_fixed/{stats.py,cmb_topology.py,pointcloud.py} was written
                  and before any acceptance number was computed" (commit 1bc7ebc
                  per validation_results.json `pre_registration`).
Tier            : the source states "tier_of_every_observed_number: X
                  (numerics)"; every run here is tier X. The expected full-sky
                  Betti vector (1, 0, 1) and Euler characteristic 2 of a closed
                  2-sphere are literature (tier L) and are carried as `expected`.

This source is the repair of the two defects the simple suite found, so the DB
gets both sides: the ORIGINAL complex, whose full-sky Betti vector at nside 64 is
(1, 0, 49147), and the FIXED complex, whose Betti vector is (1, 0, 1) with Euler
characteristic 2. Both are ingested as separate runs, the first with
`expected = (1, 0, 1)` and therefore `matches = 0`.

The chi2 p-values in D1_summary_table are calibration p-values of a TEST (a KS
test of 100 chi2 p-values against uniform), not p-values of a physical result.
They are carried with the ensemble that produced them (100 F3 realisations) as
the null. Where the source labels something post-hoc - notably
`gate_4_all_200_maps`, which `gate_4_note` calls "the post-hoc 200-map
extension" - no expectation is recorded.
"""
from __future__ import annotations

from . import _common as C

DIR = "audit/tda_validation/tda_fixed"
SCRIPT = f"{DIR}/run_validation.py"


def ingest(db) -> dict:
    src = C.Source(db, "a2_tda_fixed", "tdasimple")
    res = C.load(src.rel(DIR, "validation_results.json"))
    exp = C.load(src.rel(DIR, "fix_expectations.json"))
    parts = res.get("parts") or {}
    base = dict(pre_registration=res.get("pre_registration"),
                pre_registration_written=exp.get("written"),
                pre_registration_rules=exp.get("rules"),
                git_head_at_assemble_recorded_in_source=res.get("git_head_at_assemble"),
                tier_of_every_observed_number=res.get("tier_of_every_observed_number"),
                never_claimed=res.get("never_claimed"),
                versions=res.get("versions"),
                fixed_library_sha256=res.get("fixed_library_sha256"))

    # ============================ D2: the HEALPix surface complex, before/after
    did = "synthetic/tda_fixed/healpix_sphere_complex"
    src.dataset(id=did, domain="synthetic",
                title="the HEALPix sphere as a simplicial complex: the original build_topology and "
                      "the fixed cmb_topology.build_topology_fixed, at nside 8/16/32/64",
                source="a closed 2-sphere has Betti (1, 0, 1) and Euler characteristic 2 "
                       "(tier L, standard topology)",
                provenance="synthetic_control", ambient_dim=2, units="HEALPix pixels",
                local_path=f"{DIR}/validation_results.json")
    topo = parts.get("topology") or {}
    d2 = res.get("D2_acceptance") or {}
    last_fixed_run = None
    for nside, blk in (topo.get("levels") or {}).items():
        params = dict(base, leg="D2/fixed_complex", nside=int(nside), complex="fixed",
                      combinatorics={k: v for k, v in blk.items() if k.startswith(("n_", "predicted_",
                                                                                   "dedup_", "degenerate_"))},
                      gate=topo.get("gate"),
                      command_recorded_in_source=topo.get("_command"))
        rid = src.run(dataset_id=did, method="lower_star_graph", coeff_field=2, max_dim=2,
                      params=params, preprocessing="fixed HEALPix surface triangulation, full sky",
                      script=f"{DIR}/cmb_topology.py",
                      command=str(topo.get("_command", SCRIPT + " --step topology")), tier="X",
                      wall_sec=blk.get("build_sec"), peak_mb=topo.get("_peak_rss_mb"))
        last_fixed_run = rid   # remembered, never looked up: the database is shared
        # Betti of the fixed complex is asserted by the D2 gate 2 table, for nside 32 and 64.
        if str(nside) in (d2.get("gate_2_betti_1_0_1") or {}):
            src.betti(rid, {0: 1, 1: 0, 2: 1}, {0: 1, 1: 0, 2: 1})
        for k in ("n_vertices", "n_edges", "n_triangles", "euler_char_V_minus_E_plus_F", "npix"):
            src.stat(rid, k, blk.get(k))
        src.control(rid, "known_answer",
                    f"D2 gate 1 (nside {nside}): the fixed complex is a surface with the predicted "
                    "corner, edge and face counts and Euler characteristic 2",
                    blk.get("gate_1_combinatorial_pass"),
                    detail=f"chi = {blk.get('euler_char_V_minus_E_plus_F')}, V/E/F = "
                           f"{blk.get('n_vertices')}/{blk.get('n_edges')}/{blk.get('n_triangles')}")
        if str(nside) in (d2.get("gate_2_betti_1_0_1") or {}):
            src.control(rid, "known_answer",
                        f"D2 gate 2 (nside {nside}): the fixed complex has full-sky Betti (1, 0, 1)",
                        (d2.get("gate_2_betti_1_0_1") or {}).get(str(nside)),
                        detail="the pre-registered target for a closed 2-sphere")
    for nside, blk in (d2.get("before") or {}).items():
        params = dict(base, leg="D2/original_complex", nside=int(nside), complex="original",
                      note=blk.get("note"))
        rid = src.run(dataset_id=did, method="lower_star_graph", coeff_field=2, max_dim=2,
                      params=params,
                      preprocessing="ORIGINAL cmb_tda.build_topology, unchanged, for the before column",
                      script="audit/reverse_zero/E5-cmb-tda/cmb_tda.py",
                      command=str(topo.get("_command", SCRIPT + " --step topology"))
                              + f"  [before column, nside {nside}]",
                      tier="X", wall_sec=blk.get("wall_sec"))
        src.betti(rid, C.betti_map(blk.get("gudhi_betti")), {0: 1, 1: 0, 2: 1})
        for k in ("n_vertices", "n_edges", "n_triangles", "euler_char_V_minus_E_plus_F"):
            src.stat(rid, k, blk.get(k))
        src.control(rid, "known_answer",
                    f"the ORIGINAL complex at nside {nside} has full-sky Betti (1, 0, 1)",
                    blk.get("gudhi_betti") == [1, 0, 1],
                    detail=f"observed {blk.get('gudhi_betti')} with Euler characteristic "
                           f"{blk.get('euler_char_V_minus_E_plus_F')}")
    before64 = (d2.get("before") or {}).get("64") or {}
    src.finding(dataset_id=did, tier="X", verdict="failed",
                claim=f"D2: the ORIGINAL sky complex at nside 64 returns Betti "
                      f"{before64.get('gudhi_betti')} and Euler characteristic "
                      f"{before64.get('euler_char_V_minus_E_plus_F')} for a closed sphere, where the "
                      f"fixed complex returns (1, 0, 1) and chi = 2 at every nside tested",
                caveat="this is the defect the suite's F3 topology diagnostic exposed; both complexes "
                       "are ingested so the before/after pair is visible",
                reference=f"{DIR}/validation_results.json")
    src.note_discrepancy(
        "tda_fixed D2",
        f"the original complex reports b2 = {(before64.get('gudhi_betti') or [None, None, None])[2]} "
        "on the full sky where a closed 2-sphere has b2 = 1; the number is ingested as recorded, with "
        "the pre-stated expectation attached, so its betti row carries matches = 0")
    if last_fixed_run is not None:
        # gate 3 is a property of the fixed complex as a whole, not of one nside
        src.control(last_fixed_run, "known_answer",
                    "D2 gate 3: a disk mask on the fixed complex gives a disc",
                    d2.get("gate_3_disk_mask"), detail="pre-registered gate 3")

    # ================================== D1: the statistic's calibration, 5 columns
    did_f3 = "synthetic/tda_fixed/F3_gaussian_ensemble_curves"
    src.dataset(id=did_f3, domain="synthetic",
                title="the simple suite's committed F3 curves: 200 Gaussian isotropic HEALPix nside-64 "
                      "realisations (100 ensemble + 100 test)",
                source="exchangeability: p-values of an independent realisation against a same-C_ell "
                       "ensemble are uniform on [0,1] (tier L, construction)",
                provenance="synthetic_control", n_objects=49152, ambient_dim=2,
                units="HEALPix pixels, nside 64", local_path=f"{DIR}/validation_results.json",
                notes="the 100 test p-values share ONE 100-member ensemble, so they are correlated; "
                      "the source measures a variance inflation of 1.67-1.74 from this")
    d1a = res.get("D1_acceptance") or {}
    cal = res.get("calibration_context") or {}
    COLUMNS = {
        "A_before": "the original statistic on the committed F3 curves (before the fix)",
        "B_minimal": "the minimal enumerated fix only (dead-bin drop, df = kept bins, eigen conditioning)",
        "C_full": "the full fixed statistic on the committed F3 curves (the case the acceptance names)",
        "D0_fixed_topology_old_stats": "the fixed complex with the OLD statistic",
        "D_full_fixed_pipeline": "the fixed complex AND the fixed statistic",
    }
    for column, description in COLUMNS.items():
        params = dict(base, leg=f"D1/{column}", column=column, description=description,
                      gate=d1a.get("gate"), binomial_context=d1a.get("binomial_context"),
                      known_caveat=d1a.get("known_caveat"),
                      curves_from=(parts.get("d1new") or parts.get("d1old") or {}).get("curves_from"))
        rid = src.run(dataset_id=did_f3, method="lower_star_graph", coeff_field=2, max_dim=1,
                      params=params, preprocessing=description, script=f"{DIR}/stats.py",
                      command=str((parts.get("d1new") or {}).get("_command", SCRIPT + " --step d1new")),
                      tier="X", seed="30000..30199",
                      peak_mb=(parts.get("d1new") or {}).get("_peak_rss_mb"))
        null_desc = ("the 100-realisation same-C_ell Gaussian ensemble of the suite's F3 design; the "
                     "100 test realisations are scored against it, and the KS p reported here is the "
                     "calibration test of those 100 p-values against U(0,1)")
        for curve in ("b0", "b1", "chi"):
            blk = ((res.get("D1_summary_table") or {}).get(curve) or {}).get(column) or {}
            ks = blk.get("ks_p")
            src.stat(rid, f"{curve}_chi2_branch_KS_p_vs_uniform", ks, null_model=null_desc,
                     n_null=100, p_value=ks, p_method="chi2",
                     multiplicity="three curves (b0, b1, Euler) are gated; the source applies no "
                                  "multiplicity correction across them")
            rks = blk.get("rank_ks_p")
            src.stat(rid, f"{curve}_rank_branch_KS_p_vs_uniform", rks, null_model=null_desc,
                     n_null=100, p_value=rks, p_method="rank")
            src.stat(rid, f"{curve}_n_below_0.05", blk.get("n_below_0.05"))
            src.stat(rid, f"{curve}_rank_n_below_0.05", blk.get("rank_n_below_0.05"))
            src.stat(rid, f"{curve}_df", blk.get("df"))
            src.control(rid, "null_calibration",
                        f"D1 gate on {curve} [{column}]: {d1a.get('gate')}",
                        blk.get("gate_pass"),
                        detail=f"KS p = {ks}, #{{p<0.05}} = {blk.get('n_below_0.05')}, df = "
                               f"{blk.get('df')}")
    src.finding(dataset_id=did_f3, tier="X", verdict="recovered",
                claim="D1 on the declared case (column C, the fixed statistic on the suite's "
                      "committed F3 curves): the gate passes for b0, b1 and the Euler curve "
                      + str((res.get("headline_findings") or {}).get(
                          "1_D1_acceptance_met_on_the_declared_case"))[:300],
                caveat=str(d1a.get("known_caveat")), reference=f"{DIR}/validation_results.json")
    hf = res.get("headline_findings") or {}
    src.finding(dataset_id=did_f3, tier="X", verdict="failed",
                claim="D1 columns B and D do NOT meet the gate: "
                      + str(hf.get("2_the_minimal_enumerated_fix_is_not_sufficient"))[:300]
                      + " " + str(hf.get("3_full_fixed_pipeline_b1_still_exceeds_the_count_clause"))[:400],
                caveat="the threshold was not moved after the fact; the source instead measures the "
                       "gate's own false-failure rate",
                reference=f"{DIR}/validation_results.json")

    # the calibration experiment that measures the gate's false-failure rate
    did_cal = "synthetic/tda_fixed/calibration_multivariate_normal"
    src.dataset(id=did_cal, domain="synthetic",
                title="calibration experiment: multivariate normal, p=6, 100 independent ensembles of "
                      "100 with 100 test vectors each (10000 test vectors)",
                source=f"exact uniform target P(p < 0.05) = "
                       f"{cal.get('exact_uniform_target_P_rank_p_below_0.05')} (tier L, construction)",
                provenance="synthetic_control", n_objects=10000,
                local_path=f"{DIR}/validation_results.json", notes=str(cal.get("question")))
    rid = src.run(dataset_id=did_cal, method="lower_star_graph", coeff_field=2, max_dim=0,
                  params=dict(base, leg="calibcheck", construction=cal.get("construction"),
                              chi2_p_branch=cal.get("chi2_p_branch")),
                  preprocessing=str(cal.get("construction")), script=f"{DIR}/stats.py",
                  command=str((parts.get("calibcheck") or {}).get("_command", SCRIPT + " --step calibcheck")),
                  tier="X", seed="777", wall_sec=cal.get("wall_sec"),
                  peak_mb=(parts.get("calibcheck") or {}).get("_peak_rss_mb"))
    for k in ("exact_uniform_target_P_rank_p_below_0.05", "pooled_rate_over_all_10000_test_vectors",
              "per_ensemble_count_mean", "per_ensemble_count_sd", "binomial_sd_if_independent",
              "variance_inflation_factor", "empirical_P_count_ge_11", "empirical_P_count_ge_12",
              "empirical_P_count_ge_13"):
        src.stat(rid, k, cal.get(k))
    for k, v in (cal.get("chi2_p_branch") or {}).items():
        if isinstance(v, (int, float)) and not isinstance(v, bool):
            src.stat(rid, f"chi2_branch_{k}", v)
    src.control(rid, "null_calibration",
                "on exactly Gaussian data the pooled rank-p rejection rate matches the exact uniform "
                "target",
                abs((cal.get("pooled_rate_over_all_10000_test_vectors") or 0)
                    - (cal.get("exact_uniform_target_P_rank_p_below_0.05") or 0)) < 0.005,
                detail=f"pooled {cal.get('pooled_rate_over_all_10000_test_vectors')} against exact "
                       f"{cal.get('exact_uniform_target_P_rank_p_below_0.05')}; the chi2 branch "
                       f"rejects at "
                       f"{(cal.get('chi2_p_branch') or {}).get('pooled_rate_P_chi2_p_below_0.05')}")
    src.finding(run_id=rid, dataset_id=did_cal, tier="X", verdict="artefact",
                claim=str(hf.get("5_the_gate_itself_has_a_false_failure_rate"))[:500],
                caveat=str(hf.get("4_residual_limitation_of_the_chi2_branch"))[:500],
                reference=f"{DIR}/validation_results.json")

    # ------------------------------------------------ the negative control
    neg = parts.get("negctrl") or {}
    did_neg = "synthetic/tda_fixed/negative_control_beam6deg"
    src.dataset(id=did_neg, domain="synthetic",
                title="negative control: 20 maps with a 6 deg beam scored against the 3 deg ensemble",
                source="construction: a different beam width must be rejected (tier L)",
                provenance="synthetic_control", n_objects=20,
                local_path=f"{DIR}/validation_results.json", notes=str(neg.get("construction")))
    rid = src.run(dataset_id=did_neg, method="lower_star_graph", coeff_field=2, max_dim=1,
                  params=dict(base, leg="negctrl", construction=neg.get("construction"),
                              declared_gate=neg.get("declared_gate")),
                  preprocessing=str(neg.get("construction")), script=f"{DIR}/stats.py",
                  command=str(neg.get("_command", SCRIPT + " --step negctrl")), tier="X",
                  seed="40000..40019", peak_mb=neg.get("_peak_rss_mb"))
    for curve, blk in (neg.get("keys") or {}).items():
        med = blk.get("median_chi2_p")
        src.stat(rid, f"{curve}_median_chi2_p", med,
                 null_model="the same 100-map 3-degree ensemble, with the FIXED complex and the "
                            "FIXED statistic",
                 n_null=100, p_value=med, p_method="chi2")
        src.stat(rid, f"{curve}_n_chi2_p_below_0.05", blk.get("n_chi2_p_below_0.05"))
        src.stat(rid, f"{curve}_n_rank_p_below_0.05", blk.get("n_rank_p_below_0.05"))
        src.control(rid, "negative", f"{neg.get('declared_gate')} ({curve})",
                    blk.get("gate_pass_15_of_20"),
                    detail=f"{blk.get('n_chi2_p_below_0.05')} of {blk.get('n_control_maps')} maps "
                           f"rejected; median chi2 p = {med}")
    src.finding(dataset_id=did_neg, tier="X", verdict="recovered",
                claim=str(hf.get("6_negative_controls_do_reject"))[:400],
                caveat="rejecting a 6-degree beam against a 3-degree ensemble is a large effect; "
                       "this shows the statistic CAN fail, not that it is sensitive to small ones",
                reference=f"{DIR}/validation_results.json")

    # ------------------------------------- the smooth-field gate and its extension
    smooth = parts.get("smooth") or {}
    rid = src.run(dataset_id=did_f3, method="lower_star_graph", coeff_field=2, max_dim=1,
                  params=dict(base, leg="D2/gate_4_smooth_field", gate=smooth.get("gate"),
                              declared_tolerance=smooth.get("declared_tolerance_0.15_of_peak"),
                              nu_grid=smooth.get("nu_grid"), n_maps=1,
                              gate_4_note=d2.get("gate_4_note")),
                  preprocessing="one smooth Gaussian realisation (seed 30000) through the old and the "
                                "fixed complex",
                  script=f"{DIR}/cmb_topology.py",
                  command=str(smooth.get("_command", SCRIPT + " --step smooth")), tier="X",
                  seed="30000", peak_mb=smooth.get("_peak_rss_mb"))
    for k in ("max_abs_delta_b0", "max_abs_delta_b1", "max_b0_old", "max_b1_old",
              "rms_delta_b0", "rms_delta_b1", "peak_nu_b0_old", "peak_nu_b0_fixed",
              "peak_nu_b1_old", "peak_nu_b1_fixed", "grid_step", "euler_check_at_final_nu"):
        src.stat(rid, k, smooth.get(k))
    src.control(rid, "known_answer",
                f"D2 gate 4 (pre-registered, n = 1): {smooth.get('gate')}",
                smooth.get("gate_4_pass"), detail=str(d2.get("gate_4_note"))[:400])
    g4 = parts.get("gate4full") or {}
    all200 = d2.get("gate_4_all_200_maps") or {}
    rid = src.run(dataset_id=did_f3, method="lower_star_graph", coeff_field=2, max_dim=1,
                  params=dict(base, leg="D2/gate_4_all_200_maps", post_hoc=True,
                              post_hoc_note="gate_4_note calls this 'the post-hoc 200-map extension'; "
                                            "it is arithmetic on committed curves, is stronger than "
                                            "the pre-registered n=1 gate and does not hold for every "
                                            "realisation. No expectation is stored for it.",
                              n_maps=all200.get("n_maps"),
                              declared_tolerance=g4.get("declared_tolerance"),
                              same_seeds_verified=g4.get("same_seeds_verified")),
                  preprocessing="all 200 F3 realisations through both complexes (arithmetic on the "
                                "committed curves, no new simulation)",
                  script=f"{DIR}/run_validation.py",
                  command=str(g4.get("_command", SCRIPT + " --step gate4full")), tier="X",
                  seed="30000..30199", peak_mb=g4.get("_peak_rss_mb"))
    for grp in ("frac_delta_b0", "frac_delta_b1", "peak_shift_b0", "peak_shift_b1"):
        for k, v in (all200.get(grp) or {}).items():
            src.stat(rid, f"{grp}_{k}", v)
    src.control(rid, "known_answer",
                "post-hoc extension of gate 4 to all 200 realisations",
                all200.get("gate_4_all_maps_pass"),
                detail=f"max fractional b0 difference {(all200.get('frac_delta_b0') or {}).get('max')}, "
                       f"{(all200.get('frac_delta_b0') or {}).get('n_over_0.15')} maps over the 15% "
                       f"tolerance; b0 peak shifted by more than one grid step on "
                       f"{(all200.get('peak_shift_b0') or {}).get('n_over_one_grid_step')} maps")
    src.finding(dataset_id=did_f3, tier="X", verdict="inconclusive",
                claim=str(hf.get("8_gate_4_is_n_equals_1_as_pre_registered_and_the_200_map_extension_is_not_uniform"))[:600],
                caveat="the pre-registered gate is n = 1 and it passes; the stronger post-hoc "
                       "extension does not, and the source reports both rather than choosing one",
                reference=f"{DIR}/validation_results.json")

    # ------------------------------------------------ P2 and F1 regressions
    p2 = parts.get("p2") or {}
    did_p2 = "synthetic/tda_fixed/P2_sphere_regression"
    src.dataset(id=did_p2, domain="synthetic",
                title="P2 regression: 5000 points on S2 through the ORIGINAL alpha pipeline and the "
                      "fixed wrapper",
                source="H_*(S2) = (Z, 0, Z); Hatcher 2002 Ch. 2 (tier L)",
                provenance="synthetic_control", n_objects=p2.get("n_points"), ambient_dim=3,
                local_path=f"{DIR}/validation_results.json")
    rid = src.run(dataset_id=did_p2, method="alpha", coeff_field=2, max_dim=2,
                  params=dict(base, leg="P2_regression", gate=p2.get("gate"),
                              numbers_produced_by=p2.get("numbers_produced_by"),
                              top_persistence_H0=p2.get("top_persistence_H0"),
                              top_persistence_H1=p2.get("top_persistence_H1"),
                              top_persistence_H2=p2.get("top_persistence_H2"),
                              cap=p2.get("cap"), P_all=p2.get("P_all"),
                              wrapper_bars_identical_to_original=p2.get("wrapper_bars_identical_to_original")),
                  preprocessing="alpha complex on 5000 points on the unit sphere",
                  script=f"{DIR}/pointcloud.py",
                  command=str(p2.get("_command", SCRIPT + " --step p2")), tier="X",
                  wall_sec=p2.get("wall_sec_original"), peak_mb=p2.get("_peak_rss_mb"))
    src.betti(rid, C.betti_map(p2.get("observed_betti_vector")),
              C.betti_map(p2.get("expected_betti")))
    src.stat(rid, "cap", p2.get("cap"))
    src.stat(rid, "P_all", p2.get("P_all"))
    src.control(rid, "known_answer",
                "P2 regression: the fixed wrapper produces bars identical to the original pipeline "
                "and the observed Betti vector still matches the committed suite result",
                p2.get("pass"),
                detail=f"observed {p2.get('observed_betti_vector')}, committed suite "
                       f"{p2.get('committed_suite_observed_betti_vector')}, "
                       f"wrapper_bars_identical_to_original = "
                       f"{p2.get('wrapper_bars_identical_to_original')}")
    f1 = parts.get("f1") or {}
    did_f1 = "synthetic/tda_fixed/F1_gate5_regression"
    src.dataset(id=did_f1, domain="synthetic",
                title="D2 gate 5: the F1 grid field (256x256, 7 Gaussian wells, seed 21) through the "
                      "original and the fixed Freudenthal function",
                source="a 256x256 disc has Euler characteristic 1 (tier L)",
                provenance="synthetic_control", n_objects=65536, ambient_dim=2,
                local_path=f"{DIR}/validation_results.json")
    for which, blk in (("original", f1.get("original_function_rerun") or {}),
                       ("fixed", f1.get("fixed_function") or {})):
        rid = src.run(dataset_id=did_f1, method="lower_star_graph", coeff_field=2, max_dim=1,
                      params=dict(base, leg=f"D2/gate_5/{which}", gate=f1.get("gate"),
                                  committed_suite_result=f1.get("committed_suite_result"),
                                  curves_identical_b0=f1.get("curves_identical_b0"),
                                  curves_identical_b1=f1.get("curves_identical_b1"),
                                  euler_note=f1.get("euler_note")),
                      preprocessing=f"{which} Freudenthal sublevel function on the F1 grid field",
                      script=f"{DIR}/cmb_topology.py",
                      command=str(f1.get("_command", SCRIPT + " --step f1")), tier="X", seed="21",
                      wall_sec=blk.get("wall_sec"), peak_mb=f1.get("_peak_rss_mb"))
        src.betti(rid, {0: blk.get("b0_at_max_nu"), 1: blk.get("b1_at_max_nu")})
        src.stat(rid, "max_nu_b0", blk.get("max_nu_b0"))
        src.stat(rid, "max_nu_b1", blk.get("max_nu_b1"))
        src.stat(rid, "euler_char_true_at_final_nu", f1.get("euler_char_true_at_final_nu"))
        src.stat(rid, "b0_minus_b1_at_final_nu", f1.get("b0_minus_b1_at_final_nu"))
        src.control(rid, "known_answer", f"D2 gate 5: {f1.get('gate')}", f1.get("gate_5_pass"),
                    detail=f"curves identical b0 {f1.get('curves_identical_b0')}, b1 "
                           f"{f1.get('curves_identical_b1')}; {f1.get('euler_note')}")

    # ---------------------------------------------------- the regression guard
    guard = res.get("regression_guard") or {}
    did_g = "synthetic/tda_fixed/regression_guard"
    src.dataset(id=did_g, domain="synthetic",
                title="the regression guard test suite: 30 tests that must pass on the fixed library "
                      "and fail on the originals",
                source="construction: a guard that cannot fail is not a guard (tier L)",
                provenance="synthetic_control", n_objects=30,
                local_path=f"{DIR}/validation_results.json")
    for mode, blk in (guard.get("modes") or {}).items():
        rid = src.run(dataset_id=did_g, method="lower_star_graph", coeff_field=2, max_dim=2,
                      params=dict(base, leg=f"regression_guard/{mode}", gate=guard.get("gate"),
                                  summary_line=blk.get("summary_line"),
                                  failed_tests=blk.get("failed_tests"),
                                  exit_code=blk.get("exit_code")),
                      preprocessing=mode.replace("_", " "), script=f"{DIR}/tests",
                      command=str(blk.get("command")), tier="X", wall_sec=blk.get("wall_sec"),
                      peak_mb=guard.get("_peak_rss_mb"))
        src.stat(rid, "n_failed_tests", blk.get("n_failed"))
        src.control(rid, "known_answer", str(guard.get("gate")),
                    (blk.get("n_failed") == 0) if mode == "against_fixed_library"
                    else (blk.get("n_failed") > 0),
                    detail=str(blk.get("summary_line")))
    src.finding(dataset_id=did_g, tier="X",
                verdict="recovered" if guard.get("guard_pass") else "failed",
                claim=str(guard.get("verdict")),
                caveat=f"the 19 tests that fail on the originals are: {guard.get('modes', {}).get('against_original_implementations', {}).get('failed_tests')}",
                reference=f"{DIR}/validation_results.json")

    # ---------------------------------------------- what the fix does to E5
    e5p = parts.get("e5") or {}
    did_e5 = "astro/wmap_ilc9_e5_superseded"
    if db.con.execute("SELECT 1 FROM dataset WHERE id=?", (did_e5,)).fetchone() is None:
        src.dataset(id=did_e5, domain="astro",
                    title="WMAP 9-yr ILC map as used by the superseded E5 CMB TDA run",
                    source="NASA LAMBDA WMAP 9-yr ILC", provenance="observation", ambient_dim=2,
                    notes="SUPERSEDED; see the d2_reverse_zero backfill module")
    rid = src.run(dataset_id=did_e5, method="lower_star_graph", coeff_field=2, max_dim=1,
                  params=dict(base, leg="e5_impact", what=e5p.get("what"),
                              ensemble=e5p.get("ensemble"), n_sims=e5p.get("n_sims"),
                              nside_of_that_run=e5p.get("nside_of_that_run"), SCOPE=e5p.get("SCOPE"),
                              recomputes="the E5 p-values under the OLD, the MINIMAL and the NEW "
                                         "statistic on the same committed curves"),
                  preprocessing=str(e5p.get("SCOPE"))[:400], script=f"{DIR}/stats.py",
                  command=str(e5p.get("_command", SCRIPT + " --step e5")), tier="X",
                  peak_mb=e5p.get("_peak_rss_mb"))
    e5i = res.get("e5_impact") or {}
    nulldesc = (f"the E5 run's own {e5p.get('n_sims')}-realisation Gaussian ensemble "
                f"({e5p.get('ensemble')}), rescored")
    for stat_name, blk in e5i.items():
        for branch, key, method in (("OLD", "OLD_p_chi2", "chi2"), ("MINIMAL", "MINIMAL_p_chi2", "chi2"),
                                    ("NEW", "NEW_p_chi2", "chi2"), ("OLD_rank", "OLD_rank_p", "rank"),
                                    ("NEW_rank", "NEW_rank_p", "rank")):
            p = blk.get(key)
            src.stat(rid, f"{stat_name}_{branch}_p", p, null_model=nulldesc, n_null=e5p.get("n_sims"),
                     p_value=p, p_method=method)
        src.stat(rid, f"{stat_name}_delta_p_chi2", blk.get("delta_p_chi2"))
        src.stat(rid, f"{stat_name}_NEW_n_bins_kept", blk.get("NEW_n_bins_kept"))
    crossed = [k for k, v in e5i.items() if v.get("crosses_0.05")]
    src.finding(run_id=rid, dataset_id=did_e5, tier="X", verdict="artefact",
                claim=f"re-scoring the E5 CMB statistics with the fixed statistic moves "
                      f"{len(crossed)} of {len(e5i)} of them across p = 0.05 ({crossed}); the "
                      "largest change is "
                      + str(max(((abs(v.get('delta_p_chi2') or 0), k) for k, v in e5i.items()),
                                default=(0, ''))[1]),
                caveat=str(e5p.get("SCOPE"))[:500],
                reference=f"{DIR}/validation_results.json")

    src.skip("tda_fixed persistence diagrams",
             "this source works on Betti CURVES over a nu grid and on combinatorial counts of the "
             "sky complex; the only birth/death data it records is the P2 regression's top-3 "
             "persistence values per dimension, which are persistences, not birth/death pairs, so "
             "no bars could be stored")
    src.skip("tda_fixed f3_curve_chunks",
             "the 200 committed F3 Betti curves are stored chunked for provenance; they are the "
             "input to the D1 columns already ingested, not separate runs")
    return src.report()
