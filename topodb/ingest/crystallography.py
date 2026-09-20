"""Backfill source (b3): the A4/C12 discrete-symmetry lens on WMAP and SMICA.

Source worktree : dualscale-wt-tdaval
Source files    : audit/tda_validation/crystallography/{lens_spec.json, report.json,
                  tda_secondary_wmap.json, pointcloud_kat.json,
                  real_map_result_{wmap,planck}.json, groups_check.json}
Pre-registration: lens_spec.json, committed before any lens code was run and
                  before any real map was read (report.json
                  `design_note_and_preregistration`).

WHAT IS AND IS NOT INGESTED HERE
--------------------------------
The GATED statistic of this run is NOT persistent homology. lens_spec.json
section 7 defines the primary as `C_B(G, R)`, the fraction of band power in the
G-invariant spherical-harmonic subspace, and its p-values come from 200 Gaussian
realisations. That is a harmonic statistic; no persistence enters it. It is
therefore NOT ingested as a run, and it is recorded under `skipped` with its
verdict carried as a dataset-level finding so the DB does not lose the fact that
both maps came out NULL.

What IS ingested is the two legs where persistence was actually computed:
  * tda_secondary_wmap.json - lower-star Betti curves b0(nu), b1(nu) of the
    invariant field m_G, the residual r_G and the map T, plus gudhi
    CubicalComplex bar COUNTS on a gnomonic projection. The file itself records
    `gated: false` and `note: "secondary, reported not gated (lens_spec.json
    section 11)"`, so no expectation is stored and no p-value exists.
  * pointcloud_kat.json - a synthetic known-answer test in which the Betti
    numbers of an A4-orbit union and of a matched random control are the counts
    of gudhi AlphaComplex persistence intervals.

Both legs are tier X, as their files state.

NO BARS. `grep` over every JSON in this directory finds no birth, death,
top_bars, persistence_intervals or dgm key: the source keeps only interval
COUNTS and Betti curves. Nothing could be stored in the `bar` table.
"""
from __future__ import annotations

from . import _common as C

DIR = "audit/tda_validation/crystallography"
SCRIPT = f"{DIR}/discrete_symmetry_lens.py"


def ingest(db) -> dict:
    src = C.Source(db, "b3_crystallography", "tdaval")
    rep = C.load(src.rel(DIR, "report.json"))
    spec = C.load(src.rel(DIR, "lens_spec.json"))
    tda = C.load(src.rel(DIR, "tda_secondary_wmap.json"))
    kat = C.load(src.rel(DIR, "pointcloud_kat.json"))
    wmap = C.load(src.rel(DIR, "real_map_result_wmap.json"))
    planck = C.load(src.rel(DIR, "real_map_result_planck.json"))
    groups = C.load(src.rel(DIR, "groups_check.json"))
    base = dict(pre_registration=rep.get("design_note_and_preregistration"),
                tiers=rep.get("tiers"), seeds=rep.get("seeds"),
                hardware=rep.get("hardware"),
                what_the_lens_computes=(rep.get("what_the_lens_computes") or {}).get("tda_path"),
                gated=False,
                gating_note="lens_spec.json section 7 makes the TDA path a SECONDARY, NOT-GATED "
                            "report; the gated primary is the harmonic statistic C_B(G, R), which "
                            "is not persistent homology and is not ingested as a run")

    # --------------------------------- the lower-star secondary on WMAP
    did = "astro/crystallography/wmap_ilc9_lens_fields"
    meta = (rep.get("real_map_results") or {}).get("wmap") or {}
    src.dataset(id=did, domain="astro",
                title="WMAP ILC9 and the two fields the A4/C12 lens derives from it: the invariant "
                      "projection m_G and the residual r_G",
                source=str(meta.get("map_path") or "NASA LAMBDA WMAP ILC9"),
                provenance="observation", ambient_dim=2, units="HEALPix pixels",
                sha256=meta.get("map_sha256"), local_path=f"{DIR}/tda_secondary_wmap.json",
                notes=f"complex under test: {str(tda.get('complex'))[:400]} | "
                      f"{str(tda.get('tda_fixed_status'))[:300]}")
    for group, curves in (tda.get("curves") or {}).items():
        for field in ("mG", "rG", "T"):
            params = dict(base, leg=f"tda_secondary/{group}/{field}", group=group, field=field,
                          nu_grid=curves.get("nu_grid"),
                          b0_curve=curves.get(f"b0_{field}"), b1_curve=curves.get(f"b1_{field}"),
                          complex=tda.get("complex"),
                          tda_fixed_status=tda.get("tda_fixed_status"),
                          note=tda.get("note"),
                          b2_note="the source restricts itself to b0 and b1 and refuses to call "
                                  "b0 - b1 an Euler characteristic, because the simple suite records "
                                  "a hollow-tetrahedron defect in this complex that produces "
                                  "spurious H2 classes")
            rid = src.run(dataset_id=did, method="lower_star_graph", coeff_field=2, max_dim=1,
                          params=params,
                          preprocessing=f"lower-star persistence of {field} on the masked HEALPix "
                                        "graph over a 41-point nu grid",
                          script=SCRIPT, command=f"{SCRIPT} --stage tda --map wmap",
                          tier=tda.get("tier", "X"), peak_mb=tda.get("peak_rss_mb"))
            src.stat(rid, f"max_b1_{field}", curves.get(f"max_b1_{field}"))
            if field == "mG":
                src.stat(rid, "power_fraction_in_mG", curves.get("power_fraction_in_mG"))
                src.stat(rid, "mG_rG_normalised_overlap_on_unmasked_pixels",
                         curves.get("mG_rG_normalised_overlap_on_unmasked_pixels"))
            cub = curves.get(f"cubical_gnomonic_{field}")
            if isinstance(cub, dict):
                for k, v in cub.items():
                    src.stat(rid, f"cubical_gnomonic_{k}", v)
                src.control(rid, "known_answer",
                            "independent cross-check: gudhi CubicalComplex lower-star persistence on "
                            "a gnomonic projection centred on the north galactic pole",
                            None,
                            detail=f"n_H0_bars {cub.get('n_H0_bars')}, n_H1_bars "
                                   f"{cub.get('n_H1_bars')} over {cub.get('n_finite')} finite pixels "
                                   f"at {cub.get('reso_arcmin')} arcmin; the source states no "
                                   "pass/fail criterion for this cross-check, so no verdict is stored")
            src.control(rid, "negative",
                        f"{field}: the source reports this curve NOT GATED (secondary), so it carries "
                        "no decision and no p-value",
                        None, detail=str(tda.get("note")))
    a4, c12 = (tda.get("curves") or {}).get("A4") or {}, (tda.get("curves") or {}).get("C12") or {}
    src.finding(dataset_id=did, tier="X", verdict="inconclusive",
                claim=f"the lower-star Betti curves of the A4- and C12-invariant projections of WMAP "
                      f"reach max b1 of {a4.get('max_b1_mG')} (A4 m_G) and {c12.get('max_b1_mG')} "
                      f"(C12 m_G) against {a4.get('max_b1_T')} for the map itself; the invariant "
                      f"subspaces hold {a4.get('power_fraction_in_mG')} (A4) and "
                      f"{c12.get('power_fraction_in_mG')} (C12) of the power",
                caveat="this leg is explicitly NOT gated and has no null ensemble of its own, so no "
                       "p-value exists and none is invented; the complex under test also carries a "
                       "known hollow-tetrahedron defect recorded in the simple suite, which is why "
                       "the source uses only b0 and b1",
                reference=f"{DIR}/tda_secondary_wmap.json")
    src.note_discrepancy(
        "crystallography TDA secondary",
        "the source declares that audit/tda_validation/tda_fixed/ 'does not exist on branch "
        "loop/tda-simple', so this leg was run on the UNFIXED complex; the tda_fixed backfill "
        "(source a2) shows that complex returns b2 = 49147 on the full sky where a closed sphere "
        "has b2 = 1. This leg uses only b0 and b1, which the defect does not touch, and says so.")

    # ---------------------------------- the synthetic point-set known answer
    did_pc = "synthetic/crystallography/a4_orbit_pointcloud"
    src.dataset(id=did_pc, domain="synthetic",
                title=str(kat.get("construction"))[:300],
                source="construction: an exact A4-orbit union folds to 60 distinct points under A4 "
                       "and a random set does not (tier L, group theory)",
                provenance="synthetic_control", n_objects=720, ambient_dim=3,
                local_path=f"{DIR}/pointcloud_kat.json",
                notes=str(kat.get("not_run")))
    for leg in ("A4_orbit_union", "random_control"):
        blk = kat.get(leg) or {}
        params = dict(base, leg=f"pointcloud_kat/{leg}", construction=kat.get("construction"),
                      n_points=blk.get("n_points"),
                      n_distinct_after_A4_folding=blk.get("n_distinct_after_A4_folding"),
                      expected_n_distinct_if_exact_orbits=blk.get("expected_n_distinct_if_exact_orbits"),
                      betti_note="betti_of_folded_set is the COUNT of gudhi AlphaComplex persistence "
                                 "intervals per dimension; the intervals themselves were not saved")
        rid = src.run(dataset_id=did_pc, method="alpha", coeff_field=2, max_dim=2, params=params,
                      preprocessing=f"{leg}: fold the point set by the A4 action, then an alpha "
                                    "complex on the folded set",
                      script=SCRIPT, command=f"{SCRIPT} --stage pointcloud",
                      tier=kat.get("tier", "X"), seed=str(kat.get("seed")),
                      peak_mb=kat.get("peak_rss_mb"))
        src.betti(rid, {int(k): int(v) for k, v in (blk.get("betti_of_folded_set") or {}).items()})
        for k in ("A4_orbit_residual_mean", "A4_orbit_residual_median",
                  "C12_orbit_residual_mean", "C12_orbit_residual_median",
                  "n_distinct_after_A4_folding"):
            src.stat(rid, k, blk.get(k))
        if leg == "A4_orbit_union":
            src.control(rid, "known_answer",
                        "an exact A4-orbit union folds to the expected number of distinct points "
                        "and has zero residual under A4",
                        blk.get("n_distinct_after_A4_folding") ==
                        blk.get("expected_n_distinct_if_exact_orbits"),
                        detail=f"{blk.get('n_distinct_after_A4_folding')} distinct against an "
                               f"expected {blk.get('expected_n_distinct_if_exact_orbits')}; A4 "
                               f"residual {blk.get('A4_orbit_residual_mean')}, C12 residual "
                               f"{blk.get('C12_orbit_residual_mean')}")
        else:
            src.control(rid, "negative",
                        "a matched random set does NOT fold to the expected number of distinct points",
                        blk.get("n_distinct_after_A4_folding") !=
                        blk.get("expected_n_distinct_if_exact_orbits"),
                        detail=f"{blk.get('n_distinct_after_A4_folding')} distinct against an "
                               f"expected {blk.get('expected_n_distinct_if_exact_orbits')} for an "
                               "exact orbit set")
    gs = kat.get("group_specific") or {}
    src.finding(dataset_id=did_pc, tier="X", verdict="recovered",
                claim=f"the folding statistic is group-specific: the A4-symmetric set has residual "
                      f"{gs.get('target_residual_symmetric_over_random')} under A4 and "
                      f"{gs.get('control_group_residual_on_symmetric_set_is_NOT_small')} under C12, "
                      f"and the Betti numbers of the folded sets differ sharply "
                      f"({(kat.get('A4_orbit_union') or {}).get('betti_of_folded_set')} vs "
                      f"{(kat.get('random_control') or {}).get('betti_of_folded_set')})",
                caveat=str(kat.get("not_run")), reference=f"{DIR}/report.json")

    # ------------------------------------- the gated harmonic lens: NOT a run
    for name, res, sha_key in (("wmap", wmap, "wmap"), ("planck", planck, "planck")):
        tests = res.get("tests") or {}
        src.finding(dataset_id=did if name == "wmap" else None, tier="X", verdict="null",
                    claim=f"the GATED A4/C12 crystallinity lens on {name} is a NULL result: rank p = "
                          f"{(tests.get('A4') or {}).get('rank_p')} (A4) and "
                          f"{(tests.get('C12') or {}).get('rank_p')} (C12) against "
                          f"{res.get('n_null')} Gaussian realisations, both far above the Bonferroni "
                          f"threshold {res.get('alpha_bonferroni')}; "
                          f"detection_at_bonferroni is "
                          f"{(tests.get('A4') or {}).get('detection_at_bonferroni')} / "
                          f"{(tests.get('C12') or {}).get('detection_at_bonferroni')}",
                    caveat="this statistic is the band-power fraction in the G-invariant harmonic "
                           "subspace, NOT persistent homology, so it is recorded as a finding and "
                           "not as a run in this database. The pre-registered sensitivity ladder "
                           "found no amplitude reaching 95% power "
                           f"({(rep.get('sensitivity_preregistered') or {}).get('smallest_amplitude_with_power_0.95')}), "
                           "and WMAP and Planck observe the SAME sky, so the two p-values are "
                           "correlated rather than independent confirmations",
                    reference=f"{DIR}/real_map_result_{name}.json")
    src.skip("crystallography gated A4/C12 harmonic lens (real_map_result_*.json, injection_*.json, "
             "probe_wmap.json, bandfrac_wmap.json)",
             "the gated statistic C_B(G, R) is the fraction of band power in the G-invariant "
             "spherical-harmonic subspace; no persistent homology enters it, so no run record was "
             "created. Its NULL verdict and its p-values (with n_null = 200 and the null-model prose "
             "from lens_spec.json section 10) are carried as findings instead",
             n_null=wmap.get("n_null"),
             null_model=str((spec.get("10_nulls") or {}).get("sky"))[:400])
    src.skip("crystallography tier-B group checks (groups_check.json)",
             "exact group-theoretic checks (|Aut(D_4)| = "
             f"{(groups.get('D4_lattice') or {}).get('Aut_order_exhaustive')} with negative controls, "
             "Hurwitz closure, the 24->12 double cover, invariant dimensions by two methods, Wigner "
             "matrices) with all_checks_pass = "
             f"{groups.get('all_checks_pass')}. These are tier B but are not homology computations, "
             "so no run record was created",
             all_checks_pass=groups.get("all_checks_pass"))
    src.skip("crystallography persistence diagrams",
             "no birth, death, top_bars, persistence_intervals or dgm key exists in any JSON in this "
             "directory; the code computes persistence and keeps only interval COUNTS "
             "(n_H0_bars / n_H1_bars, betti_of_folded_set) and Betti curves")
    src.skip("crystallography cosmic-web application of the point-set lens",
             str(kat.get("not_run")))
    for item in rep.get("limitations") or []:
        src.skip("crystallography stated limitation", str(item)[:400])
    for c in rep.get("corrections_to_earlier_commit_messages") or []:
        src.note_discrepancy("crystallography correction to an earlier commit message",
                             f"{c.get('commit_subject')}: {str(c.get('correction'))[:300]}")
    bb = rep.get("sensitivity_broadband_extra_rungs_RUN_AFTER_THE_DATA") or {}
    if bb:
        src.skip("crystallography broadband sensitivity extra rungs",
                 "the section is named RUN_AFTER_THE_DATA and carries an ORDERING_WARNING: "
                 + str(bb.get("ORDERING_WARNING"))[:400])
    return src.report()
