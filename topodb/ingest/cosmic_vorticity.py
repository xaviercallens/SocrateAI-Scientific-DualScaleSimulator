"""Backfill source (d1): the cosmic-vorticity run - the Re6Zr defect-core TDA
pipeline transferred to the CMB (WMAP, SMICA) and to DESI DR1.

Source worktree : dualscale-wt-reverse, branch loop/reverse-zero
Source files    : audit/cosmic_vorticity/{registration.json, report.json,
                  results/*.json}
Pre-registration: registration.json, frozen at commit 7ecf289 before any data
                  was read (report.json `registration`).
Tier            : report.json states "X (exploratory numerics)" for the whole run.

No bars exist. The H0 bars are consumed in memory; what the source commits is
their summary (n, median, q1, q3, IQR/median, mean, cv), a Betti vector at one
truncation radius, and floor diagnostics. Those are ingested as Betti numbers
plus statistics, and the absence of diagrams is recorded under `skipped`.

Two legs deliberately carry NO p-value, and none is invented:
  * `controls.C1_uniform_no_p_value` - the source's own key name says so;
  * the whole DESI leg - `N1_clustering_matched_null.status = "NOT ATTEMPTED"`,
    and report.json's multiplicity table records `"p": null, "status": "NOT
    ATTEMPTED"` for all three DESI tests.
"""
from __future__ import annotations

from . import _common as C

DIR = "audit/cosmic_vorticity"
R = f"{DIR}/results"

CMB_DATASETS = {
    "wmap": dict(id="astro/cv/wmap_ilc9_defect_cores",
                 title="WMAP 9-yr ILC: smoothed-temperature defect-core peak set on the sphere"),
    "smica": dict(id="astro/cv/planck_smica_defect_cores",
                  title="Planck SMICA: smoothed-temperature defect-core peak set on the sphere"),
}


def _h0_block(src, rid, blk, prefix=""):
    h0 = blk.get("h0") or {}
    for k in ("n", "median", "q1", "q3", "iqr_over_median", "mean", "cv"):
        src.stat(rid, f"{prefix}h0_{k}", h0.get(k))
    for k, v in (blk.get("absolute_floor") or {}).items():
        src.stat(rid, f"{prefix}absolute_floor_{k}", v)
    for k, v in (blk.get("truncation") or {}).items():
        src.stat(rid, f"{prefix}truncation_{k}", v)
    for k in ("S1_iqr_over_median", "S2_count", "S3_psi6_site_mean", "S3_Q6_site_mean",
              "S3_n_sites_kept", "S3_n_bonds", "sigma_Ts", "var_Ts_eroded", "sigma_delta"):
        if k in blk:
            src.stat(rid, f"{prefix}{k}", blk.get(k))


def ingest(db) -> dict:
    src = C.Source(db, "d1_cosmic_vorticity", "reverse")
    rep = C.load(src.rel(DIR, "report.json"))
    reg = C.load(src.rel(DIR, "registration.json"))
    commands = rep.get("8_commands") or []
    base = dict(run="cosmic_vorticity (CV)", tier_note=rep.get("tier"),
                registration=rep.get("registration"),
                headline=rep.get("headline"),
                hyperparameters_transferred=(rep.get("1_hyperparameters_transferred") or {}).get(
                    "invariants_transferred"),
                note_no_bars="the H0 barcode is consumed in memory; the source commits only its "
                             "summary statistics, a Betti vector at one truncation radius and floor "
                             "diagnostics, so this run has no bars")

    # --------------------------------------- 2) the Re6Zr known-answer check
    kat = C.load(src.rel(R, "re6zr_known_answer_check.json"))
    did = "quantum_fluid/cv/re6zr_stm_published_point_clouds"
    src.dataset(id=did, domain="quantum_fluid",
                title="Re6Zr STM vortex point clouds at 11 fields, as published "
                      "(the source of the transferred hyperparameters)",
                source=str(kat.get("published_reference") or kat.get("source_point_clouds")),
                provenance="observation", ambient_dim=2, units="nm",
                local_path=f"{R}/re6zr_known_answer_check.json")
    rid = src.run(dataset_id=did, method="alpha", coeff_field=2, max_dim=0,
                  params=dict(base, leg="known_answer_check", statistic=kat.get("statistic"),
                              purpose=kat.get("purpose")),
                  preprocessing="alpha H0 death radii, truncated at (3 a_tri)^2",
                  script=f"{DIR}/scripts/cv_lib.py", command="cv known-answer check", tier="X")
    for row in kat.get("rows") or []:
        f = row.get("field")
        src.stat(rid, f"published_iqr_over_median_{f}", row.get("published_iqr_over_median"))
        src.stat(rid, f"recomputed_truncated_3a_{f}", row.get("recomputed_truncated_3a"))
        src.stat(rid, f"truncation_shift_in_S1_{f}", row.get("truncation_shift_in_S1"))
    src.stat(rid, "max_abs_diff_vs_published", kat.get("max_abs_diff_vs_published"))
    src.stat(rid, "max_abs_truncation_shift_in_S1", kat.get("max_abs_truncation_shift_in_S1"))
    src.control(rid, "known_answer",
                "the transferred detector reproduces the published IQR/median of the Re6Zr H0 death "
                "radii at all 11 fields, and the truncation is inert on the source data",
                bool(kat.get("reproduces_published")) and bool(kat.get("truncation_inert_on_source_data")),
                detail=f"max |difference vs published| = {kat.get('max_abs_diff_vs_published')}, "
                       f"max |truncation shift in S1| = {kat.get('max_abs_truncation_shift_in_S1')}")
    src.finding(run_id=rid, dataset_id=did, tier="X", verdict="recovered",
                claim="the CV detector reproduces the published Re6Zr IQR/median exactly "
                      f"(max absolute difference {kat.get('max_abs_diff_vs_published')}) at all 11 "
                      "fields, and the 3 a_tri truncation shifts S1 by "
                      f"{kat.get('max_abs_truncation_shift_in_S1')}",
                caveat="this checks the transfer of the statistic, not that the statistic means "
                       "anything on the sky",
                reference=f"{DIR}/report.json section 2_known_answer_check")

    # ------------------------------------------------------- 3) Test A: CMB
    for which in ("wmap", "smica"):
        ana = C.load(src.rel(R, f"cmb_analysis_{which}.json"))
        data = C.load(src.rel(R, f"cmb_data_{which}.json"))
        meta = CMB_DATASETS[which]
        did = meta["id"]
        src.dataset(id=did, domain="astro", title=meta["title"],
                    source=str(data.get("map_file")), provenance="observation",
                    n_objects=((data.get("primary") or {}).get("h0") or {}).get("n"),
                    ambient_dim=3, units="chord length on the unit sphere",
                    local_path=f"{R}/cmb_data_{which}.json",
                    notes=("secondary, the SAME sky as WMAP and NOT an independent test "
                           "(report.json 3_test_A_cmb.secondary_smica.role)") if which == "smica" else None)
        gate = ana.get("gate") or {}
        params = dict(base, leg=f"test_A_cmb/{which}",
                      role="PRIMARY" if which == "wmap" else "SECONDARY (the same sky as WMAP)",
                      n_null=ana.get("n_null"), null_seeds=ana.get("null_seeds"),
                      fsky_mask=ana.get("fsky_mask"), fsky_eroded=ana.get("fsky_eroded"),
                      gate=gate, map_file=data.get("map_file"), mask_file=data.get("mask_file"),
                      truncation=(data.get("primary") or {}).get("truncation"),
                      min_attainable_p=ana.get("min_attainable_p"))
        cmd = next((c for c in commands if f"--which {which}" in c and "null" in c), f"cv_cmb.py --which {which}")
        rid = src.run(dataset_id=did, method="alpha", coeff_field=2, max_dim=1, params=params,
                      preprocessing="smoothed map, peak detection, R_disk erosion of the mask, "
                                    "alpha H0 death radii of the peak set",
                      script=f"{DIR}/scripts/cv_cmb.py", command=cmd, tier="X",
                      seed=str(ana.get("null_seeds")))
        prim = data.get("primary") or {}
        src.betti(rid, C.betti_map(prim.get("betti_at_truncation")))
        _h0_block(src, rid, prim)
        null_desc = (f"{ana.get('n_null')} Gaussian isotropic realisations (seeds "
                     f"{ana.get('null_seeds')}) through the same mask, the same smoothing, the same "
                     "peak detector and the same alpha filtration")
        for sname, blk in (ana.get("statistics") or {}).items():
            p = blk.get("empirical_rank_p_two_sided")
            src.stat(rid, f"{sname}", blk.get("data"))
            src.stat(rid, f"p_{sname}", p, null_model=null_desc + "; convention "
                                                    + str(blk.get("p_value_convention")),
                     n_null=ana.get("n_null"), p_value=p, p_method="rank",
                     multiplicity=f"Bonferroni over {(rep.get('6_multiplicity_and_decision') or {}).get('effective_number_of_tests')} "
                                  f"tests, alpha per test "
                                  f"{(rep.get('6_multiplicity_and_decision') or {}).get('alpha_per_test')}")
            src.stat(rid, f"z_{sname}", blk.get("z"))
            src.stat(rid, f"null_mean_{sname}", blk.get("null_mean"))
            src.stat(rid, f"null_std_{sname}", blk.get("null_std"))
        src.control(rid, "null_calibration",
                    f"CV gate: {gate.get('quantity')} of the data lies inside the simulated "
                    f"distribution (|z| below {gate.get('threshold')})",
                    gate.get("passed"),
                    detail=f"z = {gate.get('z')}; decision: {gate.get('decision')}")
        src.control(rid, "known_answer",
                    "truncation inertness: the fraction of H0 deaths at or above the truncation "
                    "radius is zero in the data",
                    bool(ana.get("truncation_inert")),
                    detail=f"r_trunc = {(prim.get('truncation') or {}).get('r_trunc_deg')} deg; "
                           f"fraction at or above r_trunc = "
                           f"{(prim.get('truncation') or {}).get('frac_h0_deaths_at_or_above_r_trunc')}")
        # secondary variants of the detector, all on the same map
        for var, blk in (data.get("secondary") or {}).items():
            vparams = dict(params, leg=f"test_A_cmb/{which}/secondary/{var}",
                           secondary_variant=var,
                           secondary_note="a registered secondary detector variant, outside the "
                                          "pre-registered multiplicity family; the source reports "
                                          "no p-value for it, so none is stored")
            vrid = src.run(dataset_id=did, method="alpha", coeff_field=2, max_dim=1, params=vparams,
                           preprocessing=f"secondary detector variant '{var}'",
                           script=f"{DIR}/scripts/cv_cmb.py",
                           command=cmd + f"  [secondary variant {var}]", tier="X")
            src.betti(vrid, C.betti_map(blk.get("betti_at_truncation")))
            _h0_block(src, vrid, blk)

        if which == "wmap":
            for cname, ckind in (("C1_uniform_no_p_value", "negative"), ("C2_value_shuffle", "shuffle")):
                cblk = (ana.get("controls") or {}).get(cname) or {}
                if not cblk:
                    continue
                crid = src.run(dataset_id=did, method="alpha", coeff_field=2, max_dim=1,
                               params=dict(params, leg=f"test_A_cmb/wmap/control/{cname}",
                                           control_note=cblk.get("note"),
                                           n_realisations=cblk.get("n"), seeds=cblk.get("seeds"),
                                           no_p_value_by_design=(cname == "C1_uniform_no_p_value")),
                               preprocessing=str(cblk.get("note"))[:300],
                               script=f"{DIR}/scripts/cv_cmb.py",
                               command=f"cv_cmb.py controls --which wmap", tier="X",
                               seed=str(cblk.get("seeds")))
                for sname in ("S1_iqr_over_median", "S2_count", "S3_psi6_site_mean"):
                    sblk = cblk.get(sname) or {}
                    src.stat(crid, f"{sname}_mean", sblk.get("mean"))
                    src.stat(crid, f"{sname}_std", sblk.get("std"))
                    p = sblk.get("empirical_rank_p_vs_null")
                    if p is not None:
                        src.stat(crid, f"p_{sname}_vs_null", p, null_model=null_desc,
                                 n_null=ana.get("n_null"), p_value=p, p_method="rank")
                if cname == "C1_uniform_no_p_value":
                    src.control(crid, "negative",
                                "C1: a uniform random point set of the same size, NOT passed through "
                                "the detector; the source stores no p-value for it by design",
                                None, detail=str(cblk.get("note"))[:400])
                else:
                    fires = all((cblk.get(s) or {}).get("empirical_rank_p_vs_null") is not None
                                for s in ("S1_iqr_over_median", "S2_count", "S3_psi6_site_mean"))
                    src.control(crid, "shuffle",
                                "C2: permuting the smoothed temperature values inside the eroded "
                                "region and rerunning the FULL detector separates it from the null "
                                "on all three statistics",
                                fires,
                                detail="rank p vs null = "
                                       + str({s: (cblk.get(s) or {}).get("empirical_rank_p_vs_null")
                                              for s in ("S1_iqr_over_median", "S2_count",
                                                        "S3_psi6_site_mean")}))
            inj = ana.get("injection") or {}
            fp = inj.get("false_positive_rate_amp_zero") or {}
            irid = src.run(dataset_id=did, method="alpha", coeff_field=2, max_dim=1,
                           params=dict(params, leg="test_A_cmb/wmap/injection",
                                       criterion=inj.get("criterion"),
                                       null_interval=inj.get("null_interval"),
                                       grid=inj.get("grid"),
                                       post_hoc_marker="rows carry post_registration_extension; "
                                                       "report.json 10_post_registration_deviations "
                                                       "records that the ladder was extended beyond "
                                                       "the registered grid"),
                           preprocessing="Gaussian spot injection at a ladder of amplitude and density",
                           script=f"{DIR}/scripts/cv_cmb.py", command="cv_cmb.py inject --which wmap",
                           tier="X")
            for sname, rate in fp.items():
                src.stat(irid, f"false_positive_rate_amp_zero_{sname}", rate)
            for cell in inj.get("grid") or []:
                tag = (f"amp{cell.get('amp_over_sigmaT')}_n{cell.get('n_inj_full_sky')}"
                       f"_{cell.get('placement')}"
                       + ("_postreg" if cell.get("post_registration_extension") else ""))
                for sname in ("S1_iqr_over_median", "S2_count", "S3_psi6_site_mean", "any"):
                    src.stat(irid, f"detection_rate_{sname}_{tag}", cell.get(f"rate_{sname}"))
            sd = inj.get("smallest_detected_at_95pc_S1_grid_placement") or {}
            src.control(irid, "injection",
                        "injection sensitivity: the smallest amplitude/density cell reaching 95% "
                        "detection on S1",
                        bool(sd),
                        detail=f"{sd}; false-positive rate at amplitude zero: {fp}")
            src.finding(run_id=irid, dataset_id=did, tier="X", verdict="inconclusive",
                        claim="CMB injection sensitivity: the detector reaches 95% detection on S1 "
                              f"only at {sd.get('amp_over_sigmaT')} sigma_T with "
                              f"{sd.get('n_inj_full_sky')} full-sky injections "
                              f"(rate {sd.get('rate_S1_iqr_over_median')}); the false-positive rate "
                              f"at zero amplitude is {fp.get('any')} on the union of the three statistics",
                        caveat=f"the cell that reaches 95% is marked post_registration_extension = "
                               f"{sd.get('post_registration_extension')} in the source, i.e. it lies "
                               "outside the registered ladder",
                        reference=f"{DIR}/report.json section 10_post_registration_deviations")

        ss = ana.get("statistics") or {}
        min_p = min((b.get("empirical_rank_p_two_sided") for b in ss.values()
                     if b.get("empirical_rank_p_two_sided") is not None), default=None)
        src.finding(dataset_id=did, tier="X", verdict="null",
                    claim=f"Test A / {which}: none of S1 (H0 IQR/median), S2 (count) or S3 (psi6) "
                          f"separates the map from the {ana.get('n_null')}-realisation Gaussian null; "
                          f"the smallest two-sided rank p is {min_p}",
                    caveat=("SECONDARY: the same sky as WMAP, not an independent confirmation "
                            "(report.json corrections_to_earlier_commit_messages)" if which == "smica"
                            else "the smallest attainable p at N=500 is "
                                 f"{(ana.get('min_attainable_p') or {}).get('value')}, above the "
                                 "Bonferroni threshold, so a detection could not have been resolved "
                                 "by this run even in principle"),
                    reference=f"{DIR}/report.json section 3_test_A_cmb")

    # -------------------------------- the lower-star secondary diagnostic
    ls = C.load(src.rel(R, "cmb_lowerstar_wmap.json"))
    did = CMB_DATASETS["wmap"]["id"]
    lrid = src.run(dataset_id=did, method="lower_star_graph", coeff_field=2, max_dim=1,
                   params=dict(base, leg="test_A_cmb/wmap/lower_star_secondary",
                               n_null=ls.get("n_null"), null_seeds=ls.get("null_seeds"),
                               n_shuffle=ls.get("n_shuffle"), shuffle_seeds=ls.get("shuffle_seeds"),
                               nu_grid=ls.get("nu_grid"), filtration=ls.get("filtration"),
                               modules_used=ls.get("modules_used"),
                               complex_info=ls.get("complex_info"),
                               status="registered conditional SECONDARY item, outside the "
                                      "multiplicity family"),
                   preprocessing="lower-star persistence of the smoothed temperature on the masked "
                                 "HEALPix graph; Betti curves b0(nu), b1(nu)",
                   script=f"{DIR}/scripts/cv_lowerstar.py",
                   command=next((c for c in commands if "cv_lowerstar" in c), "cv_lowerstar.py"),
                   tier="X", seed=str(ls.get("null_seeds")))
    ci = ls.get("complex_info") or {}
    for k in ("n_vertices", "n_edges", "n_triangles", "euler_char_V_minus_E_plus_F",
              "n_full_sky_adjacency_edges", "n_full_sky_diagonals"):
        src.stat(lrid, k, ci.get(k))
    lnull = (f"{ls.get('n_null')} Gaussian isotropic realisations (seeds {ls.get('null_seeds')}), "
             "the SAME null ensemble as the alpha path")
    for curve, blk in (ls.get("data_vs_null") or {}).items():
        p = blk.get("empirical_rank_p")
        src.stat(lrid, f"p_{curve}_rank", p, null_model=lnull, n_null=ls.get("n_null"),
                 p_value=p, p_method="rank",
                 multiplicity="outside the pre-registered multiplicity family (registered "
                              "conditional secondary)")
        pc = blk.get("p_value_chi2_survival")
        src.stat(lrid, f"p_{curve}_chi2_DIAGNOSTIC_ONLY", pc,
                 null_model=lnull + f"; chi2 survival with Hartlap factor {blk.get('hartlap_factor')} "
                                    f"and df {blk.get('df')} - the source labels this branch "
                                    "DIAGNOSTIC ONLY",
                 n_null=ls.get("n_null"), p_value=pc, p_method="chi2")
        src.stat(lrid, f"{curve}_data_chi2_hartlap", blk.get("data_chi2_hartlap"))
        src.stat(lrid, f"{curve}_n_bins_kept", blk.get("n_bins_kept"))
    cv = ls.get("control_verdict") or {}
    for curve, blk in (ls.get("shuffle_control_vs_null") or {}).items():
        src.stat(lrid, f"shuffle_{curve}_median_rank_p", blk.get("median_rank_p"))
        src.stat(lrid, f"shuffle_{curve}_fraction_below_0.05", blk.get("fraction_below_0.05"))
    src.control(lrid, "shuffle",
                "value-shuffled maps are separated from the null by the lower-star Betti curves "
                "(the diagnostic has power)",
                bool(cv.get("diagnostic_stands")),
                detail=f"fraction of shuffles separated: b0 "
                       f"{cv.get('fraction_of_shuffles_separated_from_the_null_b0')}, b1 "
                       f"{cv.get('fraction_b1')}; reading: {str(cv.get('reading'))[:300]}")
    src.finding(run_id=lrid, dataset_id=did, tier="X", verdict="null",
                claim="lower-star secondary diagnostic on WMAP: the data Betti curves are consistent "
                      f"with the {ls.get('n_null')}-realisation Gaussian null (rank p = "
                      f"{(ls.get('data_vs_null') or {}).get('b0', {}).get('empirical_rank_p')} for b0 "
                      f"and {(ls.get('data_vs_null') or {}).get('b1', {}).get('empirical_rank_p')} for "
                      "b1), while 100% of the value-shuffled controls ARE separated",
                caveat="this leg is a registered conditional secondary, outside the multiplicity "
                       "family; it deliberately does NOT use "
                       "audit/reverse_zero/E5-cmb-tda/cmb_tda.py, which it names as 'the defective "
                       "originals'",
                reference=f"{R}/cmb_lowerstar_wmap.json")

    # ------------------------------------------------------ 4) Test B: DESI
    dana = C.load(src.rel(R, "desi_analysis_NGC.json"))
    ddata = C.load(src.rel(R, "desi_data_NGC.json"))
    g1 = C.load(src.rel(R, "desi_gate1_NGC.json"))
    did = "astro/cv/desi_dr1_ngc_peaks"
    src.dataset(id=did, domain="astro",
                title="DESI DR1 NGC galaxies, 0.1 <= z <= 0.4: smoothed-density peak set in 3-D",
                source="DESI DR1 LSS catalogue (NGC), WEIGHT * WEIGHT_FKP for galaxies and randoms",
                provenance="observation", n_objects=ddata.get("N_gal"), ambient_dim=3, units="Mpc/h",
                local_path=f"{R}/desi_data_NGC.json",
                notes="Test B carries NO p-value: its clustering-matched null was NOT ATTEMPTED and "
                      "the official randoms are explicitly not substituted for it")
    params = dict(base, leg="test_B_desi", N_gal=ddata.get("N_gal"), N_ran=ddata.get("N_ran"),
                  weights=ddata.get("weights"), gate1=g1,
                  N1_clustering_matched_null=dana.get("N1_lognormal"),
                  no_p_value_by_design="the clustering-matched null was NOT ATTEMPTED; report.json's "
                                       "multiplicity table records p = null, status = NOT ATTEMPTED "
                                       "for all three DESI tests")
    rid = src.run(dataset_id=did, method="alpha", coeff_field=2, max_dim=2, params=params,
                  preprocessing="volume-limited shell gate, smoothed density peaks, alpha H0 death "
                                "radii in Mpc/h",
                  script=f"{DIR}/scripts/cv_desi.py",
                  command=next((c for c in commands if "cv_desi" in c), "cv_desi.py --cap NGC"),
                  tier="X")
    prim = ddata.get("primary") or {}
    src.betti(rid, C.betti_map(prim.get("betti_at_truncation")))
    _h0_block(src, rid, prim)
    src.control(rid, "known_answer",
                f"G1 volume-limited gate: {g1.get('rule')}", g1.get("passed"),
                detail=f"{g1.get('n_interior_shells')} interior shells, max/min ratio "
                       f"{g1.get('ratio_max_over_min')}")
    n2 = dana.get("N2_randoms_control_NO_P_VALUE") or {}
    crid = src.run(dataset_id=did, method="alpha", coeff_field=2, max_dim=2,
                   params=dict(params, leg="test_B_desi/N2_randoms_control_NO_P_VALUE",
                               control_note=n2.get("note"), seeds=n2.get("seeds"), n=n2.get("n"),
                               no_p_value_by_design="the source's own key name says NO_P_VALUE"),
                   preprocessing="official randoms, shell-matched (clustering-free)",
                   script=f"{DIR}/scripts/cv_desi.py", command="cv_desi.py control --cap NGC",
                   tier="X", seed=str(n2.get("seeds")))
    for sname in ("S1_iqr_over_median", "S2_count", "S3_Q6_site_mean"):
        sblk = n2.get(sname) or {}
        src.stat(crid, f"{sname}_mean", sblk.get("mean"))
        src.stat(crid, f"{sname}_std", sblk.get("std"))
    src.control(crid, "negative",
                "N2: the official randoms, shell-matched and clustering-free, are NOT used as a "
                "null and carry no p-value",
                None, detail=str(n2.get("note"))[:400])
    inj = dana.get("injection") or {}
    src.finding(run_id=rid, dataset_id=did, tier="X", verdict="inconclusive",
                claim="Test B (DESI DR1 NGC): the volume-limited gate passes and the detector runs, "
                      "but no p-value exists because the clustering-matched null was not attempted",
                caveat=str((dana.get("N1_lognormal") or {}).get("verdict"))[:500],
                reference=f"{DIR}/report.json section 4_test_B_desi")
    sens = (rep.get("4_test_B_desi") or {}).get("injection_sensitivity") or {}
    src.finding(dataset_id=did, tier="X", verdict="inconclusive",
                claim="DESI injection sensitivity: "
                      + str(sens.get("smallest_detected_at_95pc_S1_iqr_over_median"))[:300],
                caveat=f"union false-positive rate {sens.get('false_positive_rate_union')}; "
                       "report.json section 12 records a commit 'S1 is the WRONG statistic in 3-D'",
                reference=f"{DIR}/report.json section 4_test_B_desi")

    # ----------------------------------------------- multiplicity and the rest
    mult = rep.get("6_multiplicity_and_decision") or {}
    src.finding(tier="X", verdict="inconclusive",
                claim=f"CV overall: over the pre-registered family of {mult.get('effective_number_of_tests')} "
                      f"tests ({mult.get('correction')}, alpha per test {mult.get('alpha_per_test')}) "
                      f"the minimum p is {mult.get('min_p_in_family')} and "
                      f"any_below_corrected_threshold = {mult.get('any_below_corrected_threshold')}; "
                      "three of the six tests report p = null with status NOT ATTEMPTED",
                caveat=str(mult.get("note_on_independence"))[:400],
                reference=f"{DIR}/report.json section 6_multiplicity_and_decision")
    for item in rep.get("7_absent_or_not_attempted") or []:
        src.skip("cosmic_vorticity: absent or not attempted", str(item))
    for d in rep.get("9_defects_found_and_fixed_during_this_run") or []:
        src.note_discrepancy("cosmic_vorticity defect found during the run",
                             f"{d.get('where')}: {d.get('defect')}")
    for d in rep.get("10_post_registration_deviations") or []:
        src.skip("cosmic_vorticity post-registration deviation",
                 f"registered: {d.get('registered')}; item: {d.get('item')}")
    for c in rep.get("corrections_to_earlier_commit_messages") or []:
        src.note_discrepancy("cosmic_vorticity correction to an earlier commit message",
                             f"{c.get('commit_subject')}: said {str(c.get('what_it_said'))[:200]} - "
                             f"correction: {str(c.get('correction'))[:300]}")
    src.skip("cosmic_vorticity persistence diagrams",
             "no birth/death pairs are committed anywhere in audit/cosmic_vorticity: the H0 bars are "
             "consumed in memory and only their summary (n, median, q1, q3, IQR/median, mean, cv), a "
             "Betti vector at one truncation radius and floor diagnostics were saved")
    src.skip("cosmic_vorticity 500-realisation null rows (cmb_null_*.json)",
             "each of the 500 null realisations records the same summary statistics rather than a "
             "diagram; they are the null behind the stored p-values, not separate runs")
    return src.report()
