"""Backfill source (d2): the reverse-loop CMB and cosmic-web TDA runs.

Source worktree : dualscale-wt-reverse
Source files    :
  audit/reverse_zero_r2/X1-cmb/{wmap,smica}/{gate.json, calibration.json, result.json}
  audit/reverse_zero/E5-cmb-tda/e5_cmb_tda_report.json
  audit/reverse_zero/E5-cosmic-web-tda-scaled/e5_cosmic_web_tda_scaled_report.json

SUPERSEDED RUNS. The E5-* runs are ingested and MARKED AS SUPERSEDED, because a
later run in this same repository declares their statistics function defective:
audit/cosmic_vorticity/results/cmb_lowerstar_wmap.json records
  modules_used.NOT_used = "audit/reverse_zero/E5-cmb-tda/cmb_tda.py (the
  defective originals)"
and audit/tda_validation/tda_fixed/ exists precisely to re-do them. The schema's
`run` table has no `notes` column, so the mark is carried in THREE places:
  * params_json["superseded"] on every E5 run,
  * the dataset `notes` field,
  * the `caveat` of every E5 finding.
Their numbers are still ingested verbatim, because the task is a record of what
was computed, not a record of what survived.

X1-cmb is the round-2 CMB test and is NOT superseded. It reports empirical-rank
p-values with N_null = 500 stated in the same file, so the p-values are carried
with their null. It also records `N_shortfall_vs_registered: true`
(registered_N = 2000, run N = 500); that shortfall is carried in params and in
the finding caveat.
"""
from __future__ import annotations

from . import _common as C

X1 = "audit/reverse_zero_r2/X1-cmb"
E5C = "audit/reverse_zero/E5-cmb-tda"
E5W = "audit/reverse_zero/E5-cosmic-web-tda-scaled"

SUPERSEDED = (
    "SUPERSEDED: audit/cosmic_vorticity/results/cmb_lowerstar_wmap.json lists "
    "audit/reverse_zero/E5-cmb-tda/cmb_tda.py under modules_used.NOT_used as 'the defective "
    "originals', and audit/tda_validation/tda_fixed/ was written to re-do this work with a fixed "
    "statistics function. The numbers below are recorded as this run produced them and should not "
    "be read as current results."
)

MAPS = {
    "wmap": dict(id="astro/wmap_ilc9_nside128", title="WMAP 9-yr ILC temperature map, degraded to nside 128",
                 source="NASA LAMBDA WMAP 9-yr ILC"),
    "smica": dict(id="astro/planck_smica_nside128", title="Planck SMICA temperature map, degraded to nside 128",
                  source="Planck Legacy Archive SMICA"),
}


def _x1(src, which: str) -> None:
    res = C.load(src.rel(X1, which, "result.json"))
    gate = C.load(src.rel(X1, which, "gate.json"))
    cal = C.load(src.rel(X1, which, "calibration.json"))
    did = MAPS[which]["id"]
    params = dict(
        run="X1-cmb (reverse-zero round 2)", which=which,
        nu_grid=res.get("nu"), family=res.get("family"),
        N_null=res.get("N_null"), null_seeds=res.get("null_seeds"),
        registered_N=res.get("registered_N"),
        N_shortfall_vs_registered=res.get("N_shortfall_vs_registered"),
        p_min_resolution=res.get("p_min_resolution"),
        gate=gate.get("gate"), gate_threshold=gate.get("threshold"),
        gate_max_abs_z=res.get("gate_max_abs_z"), gate_passed=res.get("gate_passed"),
        gate_var_ratio_sim_over_data=res.get("gate_var_ratio_sim_over_data"),
        calibration_note=cal.get("note"), calibration_seeds=cal.get("calibration_seeds"),
        lmax=cal.get("lmax"), nside=cal.get("nside"), fsky=cal.get("fsky"),
        data_curves=res.get("data_curves"), null_mean=res.get("null_mean"),
        decision=res.get("decision"),
        euler_uncorrected=res.get("euler_uncorrected"),
        note_no_bars="the pipeline stores Betti CURVES b(nu) on a 41-point nu grid, not birth/death "
                     "pairs, so this run has no bars",
    )
    rid = src.run(dataset_id=did, method="lower_star_graph", coeff_field=2, max_dim=1, params=params,
                  preprocessing="masked HEALPix graph at nside 128; sublevel and superlevel "
                                "filtrations of T/sigma on a nu grid from -4 to 4 in 41 steps",
                  script=f"{X1}/x1_cmb.py", command=f"{X1}/x1_cmb.py --which {which}", tier="X",
                  seed=str(res.get("null_seeds")))
    null_desc = (f"{res.get('N_null')} Gaussian isotropic realisations (seeds {res.get('null_seeds')}) "
                 "from the iteratively calibrated input spectrum, through the same mask and the same "
                 "filtration; empirical rank p-values")
    for name, p in (res.get("p_value") or {}).items():
        src.stat(rid, f"p_{name}", p, null_model=null_desc, n_null=res.get("N_null"),
                 p_value=p, p_method="rank",
                 multiplicity=f"Sidak over N_eff = {res.get('diag_neff_from_null_T_correlation')} "
                              f"effective tests; p_corr = {res.get('p_corr_sidak_Neff2')}")
    for name, t in (res.get("T_data") or {}).items():
        src.stat(rid, f"T_data_{name}", t)
    for name, v in (res.get("max_abs_std_residual_per_curve") or {}).items():
        src.stat(rid, f"max_abs_std_residual_{name}", v)
    src.stat(rid, "p_min", res.get("p_min"))
    src.stat(rid, "p_corr_sidak_Neff2", res.get("p_corr_sidak_Neff2"))
    src.stat(rid, "diag_neff_from_null_T_correlation", res.get("diag_neff_from_null_T_correlation"))
    eu = res.get("euler_uncorrected") or {}
    # The Euler statistic is reported UNCORRECTED and excluded from the family by the source,
    # so its p-value is carried with its null and its exclusion is stated in `multiplicity`.
    src.stat(rid, "p_euler_chi_uncorrected_excluded_from_family", eu.get("p"),
             null_model=null_desc, n_null=res.get("N_null"), p_value=eu.get("p"), p_method="rank",
             multiplicity="NONE: the source reports this statistic uncorrected and excludes it from "
                          "the multiplicity family (" + str(eu.get("note")) + ")")
    src.stat(rid, "T_euler_chi_uncorrected", eu.get("T"))
    src.control(rid, "null_calibration",
                f"X1 gate: {gate.get('gate')}", gate.get("passed"),
                detail=f"max|z| = {gate.get('max_abs_z')} against threshold {gate.get('threshold')}; "
                       f"pixel variance ratio sim/data = {gate.get('var_ratio_sim_over_data')}; "
                       f"decision: {gate.get('decision')}")
    hist = cal.get("history") or []
    if hist:
        src.control(rid, "null_calibration",
                    "X1 spectrum calibration: 6 fixed iterations of the per-ell correction rule",
                    None,
                    detail=f"max|z| over ell bands went from {hist[0].get('max_abs_z_bands')} "
                           f"(iteration {hist[0].get('iteration')}) to {hist[-1].get('max_abs_z_bands')} "
                           f"(iteration {hist[-1].get('iteration')}); passed=None because the source "
                           "records no pass/fail flag for the calibration itself")
    ss = res.get("string_sensitivity") or {}
    verdict = "null" if "not rejected" in str(res.get("decision", "")) else "inconclusive"
    src.finding(run_id=rid, dataset_id=did, tier="X", verdict=verdict,
                claim=f"X1-cmb / {which}: over the family {res.get('family')} the smallest rank "
                      f"p-value is {res.get('p_min')} and the Sidak-corrected p is "
                      f"{res.get('p_corr_sidak_Neff2')}. Source decision: {res.get('decision', '')[:400]}",
                caveat=f"N_null = {res.get('N_null')} against a registered N of {res.get('registered_N')} "
                       f"(N_shortfall_vs_registered = {res.get('N_shortfall_vs_registered')}), so the "
                       f"smallest attainable p is {res.get('p_min_resolution')}; the Euler statistic is "
                       "reported uncorrected and outside the family. String-injection sensitivity model: "
                       + str(ss.get("model", ""))[:300],
                reference=f"{X1}/{which}/result.json")
    src.skip(f"X1-cmb {which} per-null curves",
             "null_chunk_*.npz and string_a*_chunk_*.npz hold Betti CURVES with shape "
             "(100, 4, 41) over the nu grid, not birth/death pairs; no persistence diagram exists "
             "to store, so only the data curves, the null mean/std and the statistics are carried")


def ingest(db) -> dict:
    src = C.Source(db, "d2_reverse_zero", "reverse")
    for key, meta in MAPS.items():
        src.dataset(id=meta["id"], domain="astro", title=meta["title"], source=meta["source"],
                    provenance="observation", n_objects=196608, ambient_dim=2,
                    units="HEALPix pixels, nside 128",
                    notes="used by both the round-2 X1-cmb run and the superseded E5 run")
    for which in ("wmap", "smica"):
        _x1(src, which)

    # ------------------------------------------------------- E5 CMB (superseded)
    e5 = C.load(src.rel(E5C, "e5_cmb_tda_report.json"))
    nt = e5.get("null_test") or {}
    did = "astro/wmap_ilc9_e5_superseded"
    src.dataset(id=did, domain="astro",
                title="WMAP 9-yr ILC map as used by the SUPERSEDED E5 CMB TDA run",
                source="NASA LAMBDA WMAP 9-yr ILC", provenance="observation",
                n_objects=nt.get("n_unmasked_pixels_work_res"), ambient_dim=2,
                units="HEALPix pixels (unmasked, working resolution)",
                notes=SUPERSEDED)
    for direction in ("sublevel", "superlevel"):
        params = dict(run="E5 CMB TDA", superseded=SUPERSEDED, direction=direction,
                      experiment=e5.get("experiment"), framing=e5.get("framing"),
                      script=e5.get("script"),
                      n_sims_total=nt.get("n_sims_total"), nu_grid=nt.get("nu_grid_sigma"),
                      fsky=nt.get("fsky"), lmax_used=nt.get("lmax_used"),
                      n_edges=nt.get("n_edges"), n_triangles=nt.get("n_triangles"),
                      n_coarse_bins=nt.get("n_coarse_bins"),
                      calibration_check=nt.get("calibration_check"),
                      data_betti_curve=(nt.get("data_betti_curve") or {}).get(direction),
                      note_no_bars="Betti curves over a nu grid, not birth/death pairs")
        rid = src.run(dataset_id=did, method="lower_star_graph", coeff_field=2, max_dim=1,
                      params=params,
                      preprocessing=f"{direction} filtration on the masked HEALPix graph",
                      script=f"{E5C}/cmb_tda.py",
                      command=(e5.get("commands_exact") or ["see e5_cmb_tda_report.json"])[0],
                      tier="X", seed=str((nt.get("seeds_all") or [None])[0]))
        stats = nt.get(f"statistic_{direction}") or {}
        null_desc = (f"{nt.get('n_sims_total')} Gaussian isotropic hp.synfast realisations through the "
                     "same mask and the same filtration")
        for stat_name, blk in stats.items():
            n_sims = blk.get("n_sims") or nt.get("n_sims_total")
            src.stat(rid, f"{stat_name}_empirical_rank_p", blk.get("empirical_rank_p"),
                     null_model=null_desc, n_null=n_sims, p_value=blk.get("empirical_rank_p"),
                     p_method="rank",
                     multiplicity=f"the source tests {(e5.get('multiple_testing_note') or {}).get('n_statistics_tested')} "
                                  f"statistics; Bonferroni-corrected p = "
                                  f"{(e5.get('multiple_testing_note') or {}).get('bonferroni_corrected_p')}")
            src.stat(rid, f"{stat_name}_chi2_survival_p", blk.get("p_value_chi2_survival"),
                     null_model=null_desc + " (chi2 survival branch with the Hartlap factor "
                                            f"{blk.get('hartlap_factor')})",
                     n_null=n_sims, p_value=blk.get("p_value_chi2_survival"), p_method="chi2")
            src.stat(rid, f"{stat_name}_data_chi2_hartlap", blk.get("data_chi2_hartlap"))
        cc = nt.get("calibration_check") or {}
        src.control(rid, "null_calibration",
                    "E5: the simulated ensemble reproduces the data's unmasked pixel variance",
                    None,
                    detail=f"ratio sim/data = {cc.get('ratio_sim_over_data')} "
                           f"(data {cc.get('data_unmasked_pixel_variance')}, "
                           f"sims {cc.get('sim_ensemble_mean_variance')}); the source records no "
                           f"pass/fail threshold, so no verdict is stored. {cc.get('note')}")
    mt = e5.get("multiple_testing_note") or {}
    src.finding(dataset_id=did, tier="X", verdict="null",
                claim=f"E5 CMB TDA: of {mt.get('n_statistics_tested')} statistics the smallest raw "
                      f"p-value is {(mt.get('smallest_raw_p_value') or {}).get('p_raw')} "
                      f"({(mt.get('smallest_raw_p_value') or {}).get('statistic')}); the "
                      f"Bonferroni-corrected p is {mt.get('bonferroni_corrected_p')}, above 0.05",
                caveat=SUPERSEDED + " Residual shape: " + str(mt.get("residual_curve_shape")),
                reference=f"{E5C}/e5_cmb_tda_report.json")
    css = e5.get("cosmic_string_sensitivity") or {}
    src.finding(dataset_id=did, tier="X", verdict="inconclusive",
                claim=f"E5 cosmic-string injection sensitivity: the smallest Gmu clearing 95% of the "
                      f"scanned values is {css.get('smallest_gmu_clearing_95pct_of_scanned_values')}, "
                      f"using {css.get('statistic_used')}",
                caveat=SUPERSEDED + " " + str(css.get("interpretation"))[:400],
                reference=f"{E5C}/e5_cmb_tda_report.json")
    src.skip("E5 CMB per-simulation Betti curves (chunks/*.npz)",
             "the null ensemble is stored as Betti curves over the nu grid, not as diagrams; "
             "no bars exist to ingest")

    # ------------------------------------- E5 cosmic web (superseded, HAS bars)
    e5w = C.load(src.rel(E5W, "e5_cosmic_web_tda_scaled_report.json"))
    sc, ac = e5w.get("sample_construction") or {}, e5w.get("alpha_complex") or {}
    rr, ps = e5w.get("real_result") or {}, e5w.get("pre_stated_comparison_statistic") or {}
    nc = e5w.get("null_construction") or {}
    didw = "astro/sdss_dr17_cosmic_web_e5_superseded"
    src.dataset(id=didw, domain="astro",
                title="SDSS DR17 spectroscopic galaxies, ra 140-220, dec 0-50, 4 equal-comoving-volume "
                      "shells downsampled to 25000 points",
                source=str(sc.get("catalogue")), provenance="observation",
                n_objects=rr.get("n_points"), ambient_dim=3, units="Mpc/h",
                local_path=f"{E5W}/e5_cosmic_web_tda_scaled_report.json",
                notes=SUPERSEDED)
    params = dict(run="E5 cosmic-web TDA (scaled)", superseded=SUPERSEDED,
                  framing_rule=e5w.get("framing_rule"),
                  background_cosmology=e5w.get("background_cosmology"),
                  max_alpha_square=ac.get("max_alpha_square"),
                  r_max_persistence_mpc_over_h=ac.get("r_max_persistence_mpc_over_h"),
                  alpha_note=ac.get("note"),
                  selected_shells=sc.get("selected_shells"), n_final=sc.get("n_final"),
                  null_construction=nc.get("method"), n_null_realizations=nc.get("n_realizations"),
                  pre_stated_statistic=ps.get("definition"),
                  bars_note="the source committed only the TOP THREE bars per dimension")
    rid = src.run(dataset_id=didw, method="alpha", coeff_field=2, max_dim=2, params=params,
                  preprocessing="equal-comoving-volume shell selection then random downsample to "
                                "6250 per shell", script=f"{E5W}/cosmic_web_tda_scaled.py",
                  command=f"{E5W}/cosmic_web_tda_scaled.py", tier="X",
                  seed=str(sc.get("seed_subsample_real")), wall_sec=rr.get("runtime_sec"))
    src.betti(rid, C.betti_map(rr.get("betti_numbers_at_truncation_r50")))
    for dim in (0, 1, 2):
        src.bars(rid, dim, [(b[0], b[1]) for b in rr.get(f"top_bars_H{dim}") or []])
    for k, v in (ps.get("l2_real_vs_null_mean") or {}).items():
        # An L2 distance to the null mean, with no rank/p machinery in the source: no p-value.
        src.stat(rid, f"l2_real_vs_null_mean_{k}", v)
    for k, blk in (ps.get("r_points_outside_95pct_null_envelope") or {}).items():
        src.stat(rid, f"n_r_points_outside_95pct_envelope_{k}", blk.get("n_r_points_outside_95pct_envelope"))
        src.stat(rid, f"n_r_points_total_in_range_{k}", blk.get("n_r_points_total_in_range"))
    for k, blk in (ps.get("bottleneck_real_vs_null_top_bars") or {}).items():
        for q in ("mean", "std", "min", "max"):
            src.stat(rid, f"bottleneck_real_vs_null_{k}_{q}", blk.get(q))
    kac = e5w.get("known_answer_controls") or {}
    circ = kac.get("circle_H1") or {}
    src.control(rid, "known_answer",
                "planted circle: the pipeline detects a planted loop in H1",
                circ.get("dominant_to_second_H1_ratio") is not None,
                detail=f"n_points {circ.get('n_points')}, dominant/second H1 persistence ratio "
                       f"{circ.get('dominant_to_second_H1_ratio')}; the source records "
                       f"betti_numbers {circ.get('betti_numbers')} at truncation and no pass flag, "
                       "so 'passed' here records only that the ratio was computed")
    pv = kac.get("planted_voids_H2") or {}
    gate = pv.get("gate") or {}
    src.control(rid, "known_answer",
                f"planted voids: {gate.get('n_bars_with_death_within_20pct_of_void_radius')} H2 bars die "
                f"within 20% of the planted radius {pv.get('void_radius_mpc_or_length_unit')}, matching "
                f"the {gate.get('n_voids_planted')} planted voids",
                gate.get("pass"),
                detail=f"persistence gap bar6 -> bar7 = {gate.get('persistence_gap_bar6_to_bar7')}")
    n_out = {k: (v or {}).get("n_r_points_outside_95pct_envelope")
             for k, v in (ps.get("r_points_outside_95pct_null_envelope") or {}).items()}
    src.finding(run_id=rid, dataset_id=didw, tier="X", verdict="inconclusive",
                claim=f"E5 cosmic web: the real Betti curves fall outside the 95% null envelope at "
                      f"{n_out} r-points of 77 in the pre-stated range 2-40 Mpc/h (the excursions sit "
                      f"at small r); bottleneck distances of the top-3 bars to the 20 null realisations "
                      f"are H0 {(ps.get('bottleneck_real_vs_null_top_bars') or {}).get('H0', {}).get('mean')}, "
                      f"H1 {(ps.get('bottleneck_real_vs_null_top_bars') or {}).get('H1', {}).get('mean')}, "
                      f"H2 {(ps.get('bottleneck_real_vs_null_top_bars') or {}).get('H2', {}).get('mean')}",
                caveat=SUPERSEDED + " The source itself states a mismatch at small r is 'unsurprising "
                       "and not diagnostic on its own', because SDSS fiber-collision incompleteness "
                       "suppresses close pairs and no fiber-collision weights were applied. The null "
                       "is a redshift permutation, not a lognormal mock (deferred). No p-value is "
                       "reported by the source for any of these statistics, and none is invented here.",
                reference=f"{E5W}/e5_cosmic_web_tda_scaled_report.json")
    src.note_discrepancy(
        "E5 cosmic web",
        "the real Betti curves leave the 95% null envelope at 9 (H0), 18 (H1) and 23 (H2) of 77 "
        "r-points, yet the run reports no p-value and its own interpretation says such an excursion "
        "at small r is not diagnostic; both the counts and that caveat are ingested")
    src.skip("E5 cosmic web full diagrams and the 20 null realisations' diagrams",
             "only the top three bars per dimension of the real sample were committed; the per-null "
             "diagrams were summarised into bottleneck mean/std/min/max and not saved")
    src.skip("E5-cosmic-web e5b lognormal/Poisson null follow-up "
             "(e5b_lognormal_and_poisson_null_report.json, e5b_r_range_decomposition.json)",
             "these follow-ups report correlation-function and Betti-stack diagnostics rather than a "
             "persistent-homology run with its own Betti numbers or bars, so no run record was created")
    return src.report()
