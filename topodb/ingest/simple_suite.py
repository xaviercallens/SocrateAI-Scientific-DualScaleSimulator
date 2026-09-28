"""Backfill source (a1): the 20-case known-answer suite on simple spaces.

Source worktree : dualscale-wt-tdasimple, branch loop/tda-simple
Source files    : audit/tda_validation/simple_suite/{expectations.json, results.json, report.json}
Pre-registration: expectations.json, written 2026-09-19 before any number in
                  results.json existed; it carries `expected_betti` per case,
                  so `expected=` IS passed for the cases it covers.
Tiers           : the suite states "expected topology: L (Hatcher 2002; Kunneth;
                  construction). Every observed number: X (numerics)." Every run
                  is therefore ingested at tier X; the literature expectation is
                  carried as `expected` on the betti rows, and the L-tier source
                  for it is recorded in the dataset `source` field.

Known limitations carried across:
 * `top_bars_birth_death_capped` in results.json holds only the SIX longest bars
   per dimension, and the source has already CAPPED the death of essential bars
   at the filtration maximum. Those bars are stored as they appear, and the
   separately recorded `n_infinite` is stored as a statistic so a reader can see
   how many of the stored deaths are caps rather than real deaths.
 * p-values: report.json's `p_value_note` states they compare the weakest
   expected bar to the MAX bar of each of 50 matched uniform-box nulls (the P9
   cases; P9a `seeds` is [5000, 5049], i.e. 50). n_null=50 is carried from that
   statement, never inferred from 1/p.
"""
from __future__ import annotations

from . import _common as C

DIR = "audit/tda_validation/simple_suite"
SCRIPT = f"{DIR}/simple_suite.py"

# space -> dataset metadata. `source` carries the tier-L literature reference the
# suite's expectations.json states for the expected topology.
SPACES = {
    "P1": dict(id="synthetic/simple_suite/P1_circle_S1", title="circle S1 (unit radius) in R3, n=2000, noise 0.05",
               n_objects=2000, ambient_dim=3),
    "P2": dict(id="synthetic/simple_suite/P2_sphere_S2", title="sphere S2 (unit radius) in R3, N=5000, noise 0.02",
               n_objects=5000, ambient_dim=3),
    "P3a": dict(id="synthetic/simple_suite/P3a_torus_T2_R3", title="torus T2 in R3 (R=2.5, r=1), N=8000, noise 0.02",
                n_objects=8000, ambient_dim=3),
    "P3b": dict(id="synthetic/simple_suite/P3b_clifford_torus_R4", title="flat (Clifford) torus in R4, N=8000, no noise",
                n_objects=8000, ambient_dim=4),
    "P4": dict(id="synthetic/simple_suite/P4_figure_eight", title="figure-eight (wedge of two circles) in R3, N=2000",
               n_objects=2000, ambient_dim=3),
    "P5": dict(id="synthetic/simple_suite/P5_klein_bottle_R4", title="Klein bottle in R4 (R=2, r=1), N=8000",
               n_objects=8000, ambient_dim=4),
    "P6": dict(id="synthetic/simple_suite/P6_RP2", title="real projective plane RP2 (R4 embedding / R6 Veronese)",
               n_objects=6000, ambient_dim=4),
    "P7": dict(id="synthetic/simple_suite/P7_flat_T3_R6", title="flat 3-torus T3 in R6, N swept 250..16000",
               n_objects=16000, ambient_dim=6),
    "P8": dict(id="synthetic/simple_suite/P8_planted_voids", title="6 planted spherical voids (radius 15) in a box of side 100",
               n_objects=13727, ambient_dim=3),
    "P9a": dict(id="synthetic/simple_suite/P9a_poisson_null_N2000", title="Poisson null, N=2000, bounding box of the P1 cloud, 50 seeds",
                n_objects=2000, ambient_dim=3),
    "P9b": dict(id="synthetic/simple_suite/P9b_poisson_null_N5000", title="Poisson null, N=5000, bounding box of the P2 cloud, 50 seeds",
                n_objects=5000, ambient_dim=3),
    "P9c": dict(id="synthetic/simple_suite/P9c_poisson_null_N8000", title="Poisson null, N=8000, bounding box of the P3a cloud, 50 seeds",
                n_objects=8000, ambient_dim=3),
    "P10a": dict(id="synthetic/simple_suite/P10a_two_concentric_spheres", title="two concentric spheres, radii 1 and 2.5, N=1500+4500",
                 n_objects=6000, ambient_dim=3),
    "P10b": dict(id="synthetic/simple_suite/P10b_sphere_plus_circle", title="unit sphere plus a disjoint unit circle at (4,0,0)",
                 n_objects=4000, ambient_dim=3),
    "N1": dict(id="synthetic/simple_suite/N1_gaussian_blob", title="negative control: isotropic Gaussian blob in R3, N=2000",
               n_objects=2000, ambient_dim=3),
    "N2": dict(id="synthetic/simple_suite/N2_column_shuffled_torus", title="negative control: P3a torus with x,y,z columns independently permuted",
               n_objects=8000, ambient_dim=3),
    "C0": dict(id="synthetic/simple_suite/C0_collapse_edges_check", title="method check: 300-point subsample of the P3a cloud, L=1.2",
               n_objects=300, ambient_dim=3),
    "F1": dict(id="synthetic/simple_suite/F1_seven_gaussian_wells", title="256x256 grid field with 7 equal-depth Gaussian wells",
               n_objects=65536, ambient_dim=2, units="grid pixels"),
    "F1shuf": dict(id="synthetic/simple_suite/F1shuf_shuffled_grid", title="F1 field with its pixel values shuffled (post-hoc negative control TDA-N3)",
                   n_objects=65536, ambient_dim=2, units="grid pixels"),
    "F3": dict(id="synthetic/simple_suite/F3_healpix_gaussian_ensemble", title="HEALPix nside 64 full sky, 200 Gaussian isotropic realisations (100 ensemble + 100 test)",
               n_objects=49152, ambient_dim=2, units="HEALPix pixels, nside 64"),
}


def _dataset_kwargs(case_id, exp_case):
    meta = dict(SPACES[case_id])
    meta.setdefault("domain", "synthetic")
    meta.setdefault("provenance", "synthetic_control")
    src = []
    if exp_case:
        for k in ("space", "source", "construction"):
            if exp_case.get(k):
                src.append(f"{k}: {exp_case[k]}")
    meta["source"] = "generated by simple_suite.py; " + " | ".join(src) if src else "generated by simple_suite.py"
    meta["local_path"] = f"{DIR}/results.json"
    return meta


def _per_dim(src, run_id, per_dim, prefix=""):
    """Store bars, n_bars, n_infinite and the dominance-rule value for each dimension."""
    for dim_s, info in sorted((per_dim or {}).items(), key=lambda kv: int(kv[0])):
        dim = int(dim_s)
        bars = info.get("top_bars_birth_death_capped") or []
        src.bars(run_id, dim, [(b[0], b[1]) for b in bars if isinstance(b, list) and len(b) == 2])
        src.stat(run_id, f"{prefix}n_bars_H{dim}", info.get("n_bars"))
        src.stat(run_id, f"{prefix}n_infinite_H{dim}", info.get("n_infinite"))
        v = C.rule_value(info.get("rule"))
        if v is not None:
            src.stat(run_id, f"{prefix}dominance_rule_value_H{dim}", v)


def _rules(per_dim) -> dict:
    return {f"H{d}": (i or {}).get("rule") for d, i in sorted((per_dim or {}).items())}


def ingest(db) -> dict:
    src = C.Source(db, "a1_simple_suite", "tdasimple")
    res = C.load(src.rel(DIR, "results.json"))
    rep = C.load(src.rel(DIR, "report.json"))
    exp = C.load(src.rel(DIR, "expectations.json"))
    exp_by_id = {c["id"]: c for c in exp["cases"]}
    cases = res["cases"]
    head_at_aggregate = res.get("git_head_at_aggregate")
    base_params = {
        "git_head_at_aggregate_recorded_in_source": head_at_aggregate,
        "expectations_pre_registered": "audit/tda_validation/simple_suite/expectations.json, written "
                                       + str(exp.get("written")) + " on branch " + str(exp.get("branch")),
        "tier_note": rep.get("tiers"),
        "bars_note": "top_bars_birth_death_capped holds only the 6 longest bars per dimension and the "
                     "source caps the death of essential bars at the filtration maximum",
    }

    # ---- the seven two-engine cases -------------------------------------
    for cid in ["P1", "P2", "P3a", "P3b", "P4", "P10a", "P10b"]:
        case = cases[cid]
        e = exp_by_id.get(cid, {})
        did = src.dataset(**_dataset_kwargs(cid, e))
        expected = C.betti_map(case.get("expected_betti"))
        for engine, key_obs, key_dim, method, fnkey in [
            ("pipeline", "observed_betti_pipeline", "per_dim_pipeline", "alpha", "function"),
            ("rips_reference", "observed_betti_rips_reference", "per_dim_reference", "rips", "reference_function"),
        ]:
            params = dict(base_params, case=cid, engine=engine, function=case.get(fnkey),
                          dominance_rules=_rules(case.get(key_dim)),
                          pass_rule_pre_registered=e.get("pass"),
                          expected_betti_source=e.get("source"), expected_betti_tier=e.get("tier_expected"))
            rid = src.run(dataset_id=did, method=method, coeff_field=2, max_dim=len(case[key_obs]) - 1,
                          params=params, preprocessing=e.get("construction"), script=SCRIPT,
                          command=case.get("command", "see results.json"), tier="X",
                          seed=str(e.get("construction", ""))[:200] or None,
                          wall_sec=case.get("wall_sec_total"), peak_mb=case.get("peak_rss_mb"))
            obs = C.betti_map(case[key_obs])
            # The pre-registration states one expected Betti vector for the case; both engines
            # were gated against it, so both carry `expected`.
            src.betti(rid, obs, expected if expected else None)
            _per_dim(src, rid, case.get(key_dim))
            passed = case.get("pass") if engine == "pipeline" else case.get("reference_pass")
            src.control(rid, "known_answer",
                        f"{cid}: observed Betti equals the pre-stated {case.get('expected_betti')} "
                        f"under the pre-registered dominance rules",
                        passed, detail=str(_rules(case.get(key_dim))))
            if engine == "pipeline" and case.get("p_values"):
                for name, p in case["p_values"].items():
                    src.stat(rid, f"p_{name}", p, null_model=(
                        f"50 matched uniform-box Poisson nulls (case {case.get('null')}); the weakest "
                        f"expected bar is compared to the MAX bar of each null (report.json p_value_note)"),
                        n_null=50, p_value=p, p_method="rank")
                src.stat(rid, "p_value_gate_pass", 1.0 if case.get("p_value_gate_pass") else 0.0)
            verdict = "recovered" if passed else "failed"
            src.finding(run_id=rid, dataset_id=did, tier="X",
                        claim=(f"{cid} ({e.get('space', cid)}): {engine} observed Betti {case[key_obs]} "
                               f"against the pre-stated {case.get('expected_betti')}"),
                        verdict=verdict,
                        caveat=("single construction and seed; observed numbers are tier X (numerics), the "
                                "expectation is tier L (" + str(e.get("source")) + ")"),
                        reference=f"{DIR}/report.json")

    # P2/P3a/P3b/P10a/P10b: the Rips reference disagrees with the pipeline. Record it.
    for cid in ["P2", "P3a", "P3b", "P10a", "P10b"]:
        c = cases[cid]
        if not c.get("reference_pass"):
            src.note_discrepancy(
                f"simple_suite {cid}",
                f"the Rips reference engine observed {c['observed_betti_rips_reference']} where the "
                f"alpha pipeline observed {c['observed_betti_pipeline']} and the pre-registration "
                f"expected {c['expected_betti']}; both engines are ingested as separate runs")

    # ---- P5 and P6: variants x parts x field ----------------------------
    for cid in ["P5", "P6"]:
        case = cases[cid]
        e = exp_by_id.get(cid, {})
        did = src.dataset(**_dataset_kwargs(cid, e))
        exp_by_field = case.get("expected_betti") or {}
        for variant, vinfo in (case.get("variants") or {}).items():
            for part, pinfo in (vinfo.get("parts") or {}).items():
                if pinfo.get("status") != "ok":
                    src.skip(f"{cid} variant {variant} part {part}",
                             f"the source records status={pinfo.get('status')!r} error={pinfo.get('error')!r}: "
                             f"no Betti numbers were produced",
                             wall_sec=pinfo.get("wall_sec_total"), peak_rss_mb=pinfo.get("peak_rss_mb"))
                    continue
                for key, block in pinfo.items():
                    if not isinstance(block, dict) or "observed" not in block:
                        continue
                    field = 3 if "Z/3" in key else 2
                    method = "rips" if key.startswith("rips") else "alpha"
                    exp_vec = exp_by_field.get("Z/3" if field == 3 else "Z/2")
                    params = dict(base_params, case=cid, variant=variant, part=part, engine_key=key,
                                  function=pinfo.get("function"),
                                  dominance_rules=_rules(block.get("per_dim")),
                                  n_simplices=pinfo.get("direct_alpha_n_simplices"),
                                  expected_betti_by_field=exp_by_field,
                                  expected_betti_source=e.get("source"))
                    rid = src.run(dataset_id=did, method=method, coeff_field=field,
                                  max_dim=len(block["observed"]) - 1, params=params,
                                  preprocessing=e.get("construction"), script=SCRIPT,
                                  command=pinfo.get("command", "see results.json"), tier="X",
                                  wall_sec=pinfo.get("wall_sec_total"), peak_mb=pinfo.get("peak_rss_mb"))
                    obs = C.betti_map(block["observed"])
                    expected = C.betti_map(exp_vec) if exp_vec else None
                    # the observed vector can be longer than the expectation (ambient H3); only
                    # the pre-stated dimensions carry `expected`.
                    src.betti(rid, obs, expected)
                    _per_dim(src, rid, block.get("per_dim"))
                    src.control(rid, "known_answer",
                                f"{cid} [{variant}/{part}/{key}] against the pre-stated {exp_vec}",
                                block.get("pass"), detail=str(_rules(block.get("per_dim"))))
                    src.finding(run_id=rid, dataset_id=did, tier="X",
                                claim=f"{cid} [{variant}/{part}/{key}] observed Betti {block['observed']} "
                                      f"against the pre-stated {exp_vec}",
                                verdict="recovered" if block.get("pass") else "failed",
                                caveat=str(case.get("pass_note")) if case.get("pass_note") else None,
                                reference=f"{DIR}/report.json")
    for kc in res.get("killed_or_crashed_runs") or []:
        src.skip("simple_suite killed/crashed run",
                 f"exit_code={kc.get('exit_code')}: no result was produced",
                 args=kc.get("args"), command=kc.get("command"))

    # ---- P7: the N sweep plus the post-hoc landmark runs -----------------
    case = cases["P7"]
    e = exp_by_id["P7"]
    did = src.dataset(**_dataset_kwargs("P7", e))
    expected = C.betti_map(case.get("expected_betti"))
    for row in case.get("sweep") or []:
        params = dict(base_params, case="P7", N=row.get("N"), function=case.get("function"),
                      dominance_rules=row.get("rules"),
                      simplices_after_expansion=row.get("simplices_after_expansion"),
                      max_edge_length=1.6, expansion_dim=4)
        rid = src.run(dataset_id=did, method="rips", coeff_field=2, max_dim=3, params=params,
                      preprocessing=f"N={row.get('N')} uniform angles, seed 7", script=SCRIPT,
                      command=f"{SCRIPT} --case P7 (N={row.get('N')})", tier="X", seed="7",
                      wall_sec=row.get("wall_sec"), peak_mb=row.get("peak_rss_mb"))
        src.betti(rid, C.betti_map(row.get("observed")), expected)
        for d, rule in (row.get("rules") or {}).items():
            v = C.rule_value(rule)
            if v is not None:
                src.stat(rid, f"dominance_rule_value_H{d}", v)
        src.control(rid, "known_answer", f"P7 T3 at N={row.get('N')} against the pre-stated [1,3,3,1]",
                    row.get("pass"), detail=str(row.get("rules")))
    for row in case.get("post_hoc_landmarks") or []:
        params = dict(base_params, case="P7", post_hoc=True,
                      post_hoc_note="landmark subsampling was NOT pre-registered; no `expected` is "
                                    "recorded for these runs, so `matches` stays NULL",
                      n_landmarks=row.get("n_landmarks"), pool=row.get("pool"),
                      dominance_rules=row.get("rules"))
        rid = src.run(dataset_id=did, method="rips", coeff_field=2, max_dim=3, params=params,
                      preprocessing=f"post hoc: {row.get('n_landmarks')} landmarks from a pool of {row.get('pool')}",
                      script=SCRIPT, command=f"{SCRIPT} --case P7 post-hoc landmarks", tier="X",
                      wall_sec=row.get("wall_sec"), peak_mb=row.get("peak_rss_mb"))
        src.betti(rid, C.betti_map(row.get("observed")))       # post hoc: no `expected`
        for d, rule in (row.get("rules") or {}).items():
            v = C.rule_value(rule)
            if v is not None:
                src.stat(rid, f"dominance_rule_value_H{d}", v)
    src.finding(dataset_id=did, tier="X",
                claim="P7: the flat 3-torus T3 in R6 was not recovered within the run budget; the "
                      f"source records smallest_N_recovered = {case.get('smallest_N_recovered')!r} over "
                      f"N = 250..16000 and two post-hoc landmark runs",
                verdict="failed",
                caveat="a budget limit of Rips in R6, not a disagreement with the expected Betti vector "
                       "[1,3,3,1]; H3 was recovered at N>=2000 while H1/H2 were not",
                reference=f"{DIR}/report.json")

    # ---- P8: planted voids ----------------------------------------------
    case = cases["P8"]
    e = exp_by_id["P8"]
    did = src.dataset(**_dataset_kwargs("P8", e))
    s42 = case["seed42"]
    params = dict(base_params, case="P8", function=case.get("function"),
                  expected=case.get("expected"), max_alpha_sq=37.5 ** 2, r_trunc=37.5,
                  n_points=s42.get("n_points"), n_voids_planted=s42.get("n_voids_planted"))
    rid = src.run(dataset_id=did, method="alpha", coeff_field=2, max_dim=2, params=params,
                  preprocessing=e.get("construction"), script=SCRIPT,
                  command=case.get("command", "see results.json"), tier="X", seed="42",
                  wall_sec=case.get("wall_sec_total"), peak_mb=case.get("peak_rss_mb"))
    src.bars(rid, 2, [(b[0], b[1]) for b in s42.get("top8_H2_bars_birth_death") or []])
    for i, v in enumerate(s42.get("death_rel_err_top6") or []):
        src.stat(rid, f"H2_death_relative_error_bar{i + 1}", v)
    src.stat(rid, "ratio_p6_over_p7", s42.get("ratio_p6_over_p7"))
    src.stat(rid, "pipeline_gate_20pct_count", s42.get("pipeline_gate_20pct_count"))
    src.stat(rid, "pipeline_gate_gap", s42.get("pipeline_gate_gap"))
    src.stat(rid, "robustness_seeds_passing_of_5", case.get("robustness_pass_count_of_5"))
    src.control(rid, "known_answer",
                "P8: 6 H2 bars whose death radius is within 5% of the planted void radius 15, and the 7th is not",
                s42.get("tda5_rule_pass"),
                detail=f"top-6 relative death errors {s42.get('death_rel_err_top6')}; "
                       f"p6/p7 = {s42.get('ratio_p6_over_p7')}")
    src.control(rid, "known_answer", "P8 robustness across 5 seeds",
                case.get("robustness_pass_count_of_5") == 5,
                detail=f"{case.get('robustness_pass_count_of_5')} of 5 seeds passed")
    src.finding(run_id=rid, dataset_id=did, tier="X",
                claim="P8: the alpha pipeline recovers the 6 planted voids with H2 deaths within 5% of the "
                      f"planted radius 15 (deaths {s42.get('top8_H2_bars_birth_death', [[None, None]])[0][1]:.4f} "
                      "and five more), and the 7th bar is far from 15",
                verdict="recovered" if case.get("pass") else "failed",
                caveat="report.json records the ungated margin ratio_p6_over_p7 = "
                       f"{case.get('seed42', {}).get('ratio_p6_over_p7')}, which is below the 5x threshold "
                       "used elsewhere in the suite but was not a gate for this case",
                reference=f"{DIR}/report.json")

    # ---- P9a/b/c: the matched Poisson nulls ------------------------------
    e9 = exp_by_id["P9"]
    for cid in ["P9a", "P9b", "P9c"]:
        case = cases[cid]
        did = src.dataset(**_dataset_kwargs(cid, e9))
        params = dict(base_params, case=cid, matched_to=case.get("matched_to"), N=case.get("N"),
                      box_lo=case.get("box_lo"), box_hi=case.get("box_hi"), seeds=case.get("seeds"),
                      function=case.get("function"), pass_rule_pre_registered=e9.get("pass"))
        rid = src.run(dataset_id=did, method="alpha", coeff_field=2, max_dim=2, params=params,
                      preprocessing=f"50 seeds {case.get('seeds')}, homogeneous Poisson in the "
                                    f"{case.get('matched_to')} bounding box", script=SCRIPT,
                      command=case.get("command", "see results.json"), tier="X",
                      seed=str(case.get("seeds")), wall_sec=case.get("wall_sec_total"),
                      peak_mb=case.get("peak_rss_mb"))
        for dim in (1, 2):
            for q, v in (case.get(f"max_persistence_H{dim}") or {}).items():
                src.stat(rid, f"max_persistence_H{dim}_{q}", v)
            for q, v in (case.get(f"ratio_p1_p2_H{dim}") or {}).items():
                src.stat(rid, f"ratio_p1_p2_H{dim}_{q}", v)
        src.stat(rid, "n_seeds_no_dominant_bar", case.get("n_seeds_no_dominant_bar"))
        src.control(rid, "null_calibration",
                    f"{cid}: no dominant bar in H1 or H2 in at least 48 of 50 Poisson seeds",
                    case.get("pass"),
                    detail=f"{case.get('n_seeds_no_dominant_bar')} of 50 seeds had no dominant bar")
        src.finding(run_id=rid, dataset_id=did, tier="X",
                    claim=f"{cid}: the matched Poisson null produces no dominant H1 or H2 bar in "
                          f"{case.get('n_seeds_no_dominant_bar')} of 50 seeds",
                    verdict="null" if case.get("pass") else "failed",
                    caveat="this is the null the P1/P2/P3a p-values are taken against; "
                           "a uniform box does not match non-uniform marginals (see case N2)",
                    reference=f"{DIR}/report.json")

    # ---- N1 and N2: negative controls ------------------------------------
    case, e = cases["N1"], exp_by_id["N1"]
    did = src.dataset(**_dataset_kwargs("N1", e))
    rid = src.run(dataset_id=did, method="alpha", coeff_field=2, max_dim=2,
                  params=dict(base_params, case="N1", function=case.get("function"),
                              pass_rule_pre_registered=e.get("pass")),
                  preprocessing=e.get("construction"), script=SCRIPT,
                  command=case.get("command", ""), tier="X", seed="31",
                  wall_sec=case.get("wall_sec_total"), peak_mb=case.get("peak_rss_mb"))
    for dim in (1, 2):
        blk = case.get(f"H{dim}") or {}
        src.stat(rid, f"p1_H{dim}", blk.get("p1"))
        src.stat(rid, f"p2_H{dim}", blk.get("p2"))
        src.stat(rid, f"ratio_p1_p2_H{dim}", blk.get("ratio"))
    src.control(rid, "negative", "N1: a Gaussian blob shows no dominant H1 or H2 bar (p1 < 5 p2)",
                case.get("pass"),
                detail=f"H1 ratio {case.get('H1', {}).get('ratio')}, H2 ratio {case.get('H2', {}).get('ratio')}")
    src.finding(run_id=rid, dataset_id=did, tier="X",
                claim="N1: the isotropic Gaussian blob produces no dominant H1 or H2 bar",
                verdict="null", caveat="one seed", reference=f"{DIR}/report.json")

    case, e = cases["N2"], exp_by_id["N2"]
    did = src.dataset(**_dataset_kwargs("N2", e))
    rid = src.run(dataset_id=did, method="alpha", coeff_field=2, max_dim=2,
                  params=dict(base_params, case="N2", function=case.get("function"),
                              pass_rule_pre_registered=e.get("pass"),
                              H1_top5_persistence=case.get("H1_top5"), H2_top3_persistence=case.get("H2_top3")),
                  preprocessing=e.get("construction"), script=SCRIPT, command=case.get("command", ""),
                  tier="X", seed="32,33,34", wall_sec=case.get("wall_sec_total"),
                  peak_mb=case.get("peak_rss_mb"))
    src.stat(rid, "ratio_p2_p3_H1", case.get("H1_p2_over_p3"))
    p = case.get("p_value_second_H1_bar_vs_P9c")
    src.stat(rid, "p_second_H1_bar_vs_P9c", p, null_model=(
        "the 50-seed matched uniform-box Poisson null P9c (N=8000 in the P3a bounding box); "
        "rank of the 2nd H1 bar against the max bar of each null"), n_null=50, p_value=p, p_method="rank")
    src.control(rid, "shuffle", "N2: independently permuting the x,y,z columns of the P3a torus destroys "
                                "its two dominant H1 bars AND makes the 2nd H1 bar insignificant vs P9c",
                case.get("pass"),
                detail=f"H1 p2/p3 = {case.get('H1_p2_over_p3')} (the ratio half passes); "
                       f"p = {p} against P9c (the p-value half fails)")
    src.finding(run_id=rid, dataset_id=did, tier="X",
                claim="N2: the column shuffle destroyed the two dominant H1 bars of the P3a torus "
                      f"(p2/p3 = {case.get('H1_p2_over_p3')}), but the combined pre-registered rule is recorded as FAIL",
                verdict="failed",
                caveat=str(next((r.get("note") for r in rep.get("table", []) if r.get("id") == "N2"), None)),
                reference=f"{DIR}/report.json")

    # ---- C0: collapse_edges preserves the diagram -------------------------
    case, e = cases["C0"], exp_by_id["C0"]
    did = src.dataset(**_dataset_kwargs("C0", e))
    rid = src.run(dataset_id=did, method="rips", coeff_field=2, max_dim=2,
                  params=dict(base_params, case="C0", function=case.get("function"),
                              no_collapse=case.get("no_collapse"), collapse=case.get("collapse"),
                              versions=case.get("versions"),
                              git_head_recorded_in_source=case.get("git_head"),
                              script_sha256=case.get("script_sha256")),
                  preprocessing=e.get("construction"), script=SCRIPT, command=case.get("command", ""),
                  tier="X", seed="40", wall_sec=case.get("wall_sec_total"), peak_mb=case.get("peak_rss_mb"))
    for d, v in (case.get("bottleneck_by_dim") or {}).items():
        src.stat(rid, f"bottleneck_collapse_vs_no_collapse_H{d}", v)
    src.control(rid, "known_answer",
                "C0: collapse_edges leaves the Rips diagram unchanged (bottleneck <= 1e-12 in dims 0,1,2)",
                case.get("pass"), detail=str(case.get("bottleneck_by_dim")))
    src.finding(run_id=rid, dataset_id=did, tier="X",
                claim="C0: edge collapse changed the Rips diagram by a bottleneck distance at or below "
                      "1.6e-308 in dims 0,1,2 while cutting the expanded complex from 18441 to 1270 simplices",
                verdict="recovered", caveat="one point cloud, one L",
                reference=f"{DIR}/report.json")

    # ---- F1 and F1shuf: cubical / Freudenthal on a grid field ------------
    for cid, kind, post_hoc in [("F1", "known_answer", False), ("F1shuf", "negative", True)]:
        case = cases[cid]
        e = exp_by_id.get("F1", {})
        did = src.dataset(**_dataset_kwargs(cid, e if cid == "F1" else {}))
        for engine, method, block in [("cubical", "cubical", case.get("cubical") or {}),
                                      ("freudenthal_lower_star", "lower_star_graph", case.get("pipeline") or {})]:
            params = dict(base_params, case=cid, engine=engine, function=block.get("function"),
                          grid=case.get("grid"), centres=case.get("centres"), sigma_px=case.get("sigma_px"),
                          post_hoc=post_hoc or None,
                          post_hoc_note=("F1shuf is labelled post_hoc=true in results.json (TDA-N3); "
                                         "no pre-stated Betti expectation is recorded for it")
                          if post_hoc else None,
                          direct_local_minima_8nbr=case.get("direct_local_minima_8nbr"),
                          git_head_recorded_in_source=case.get("git_head"),
                          versions=case.get("versions"))
            rid = src.run(dataset_id=did, method=method, coeff_field=2, max_dim=1, params=params,
                          preprocessing="sublevel filtration of the grid field", script=SCRIPT,
                          command=case.get("command", ""), tier="X", seed="21",
                          wall_sec=block.get("wall_sec"), peak_mb=case.get("peak_rss_mb"))
            if engine == "cubical":
                src.stat(rid, "n_H0_bars_positive_persistence", block.get("n_H0_bars_positive_persistence"))
                src.stat(rid, "n_H0_bars_total", block.get("n_H0_bars_total"))
                obs = {0: block.get("n_H0_bars_total")}
                src.betti(rid, obs, {0: 7} if cid == "F1" else None)
            else:
                for k in ("n_edges", "n_triangles", "nu_points", "max_nu_b0", "max_nu_b1",
                          "b0_at_max_nu", "b1_at_max_nu"):
                    src.stat(rid, k, block.get(k))
                src.betti(rid, {0: block.get("b0_at_max_nu"), 1: block.get("b1_at_max_nu")})
            src.stat(rid, "direct_local_minima_8nbr", case.get("direct_local_minima_8nbr"))
            src.control(rid, kind,
                        (f"{cid}: sublevel H0 bar count equals the 7 local minima of the field"
                         if cid == "F1" else
                         f"{cid}: on the value-shuffled field the H0 count differs from 7 "
                         f"({case.get('pass_rule')})"),
                        case.get("pass"),
                        detail=f"cubical H0 bars {(case.get('cubical') or {}).get('n_H0_bars_total')}, "
                               f"direct minima {case.get('direct_local_minima_8nbr')}, "
                               f"pipeline max_nu b0 {(case.get('pipeline') or {}).get('max_nu_b0')}")
        verdict = "recovered" if cid == "F1" else "null"
        src.finding(dataset_id=did, tier="X",
                    claim=(f"{cid}: CubicalComplex counts {(case.get('cubical') or {}).get('n_H0_bars_total')} "
                           f"sublevel H0 bars and the direct 8-neighbour detector "
                           f"{case.get('direct_local_minima_8nbr')} local minima; the Freudenthal pipeline "
                           f"reports max_nu b0 = {(case.get('pipeline') or {}).get('max_nu_b0')}"),
                    verdict=verdict if case.get("pass") else "failed",
                    caveat=("post-hoc negative control (results.json marks post_hoc=true)" if post_hoc else
                            "the pipeline and the cubical engine agree here; on the shuffled field (F1shuf) "
                            "they do not (6340 vs 7341)"),
                    reference=f"{DIR}/report.json")
    if cases["F1shuf"]["pipeline"]["max_nu_b0"] != cases["F1shuf"]["cubical"]["n_H0_bars_total"]:
        src.note_discrepancy(
            "simple_suite F1shuf",
            f"on the shuffled field the Freudenthal pipeline reports max_nu b0 = "
            f"{cases['F1shuf']['pipeline']['max_nu_b0']} while CubicalComplex and the direct detector both "
            f"report {cases['F1shuf']['cubical']['n_H0_bars_total']}; both numbers are ingested, on "
            "separate runs, exactly as the source records them")

    # ---- F2: periodic cubical tori T1..T4 --------------------------------
    case, e = cases["F2"], exp_by_id["F2"]
    for name, t in (case.get("tori") or {}).items():
        n = len(t.get("grid") or [])
        did = src.dataset(id=f"synthetic/simple_suite/F2_periodic_cubical_{name}", domain="synthetic",
                          title=f"periodic cubical torus {name} on a {t.get('grid')} grid, field 1 + 1e-3 N(0,1)",
                          source=f"expected beta_k = binomial({n}, k); Kunneth (tier L per expectations.json)",
                          provenance="synthetic_control", n_objects=None, ambient_dim=n,
                          local_path=f"{DIR}/results.json")
        for periodic, method_note, betti_key, expected in [
            (True, "gudhi.PeriodicCubicalComplex", "betti_numbers", t.get("expected")),
            (False, "gudhi.CubicalComplex (non-periodic negative control)", "nonperiodic_control_betti", None),
        ]:
            params = dict(base_params, case="F2", torus=name, grid=t.get("grid"), periodic=periodic,
                          function=case.get("function"), engine=method_note,
                          expected_betti_pre_registered=t.get("expected"),
                          infinite_bars_per_dim=t.get("infinite_bars_per_dim") if periodic else None,
                          git_head_recorded_in_source=case.get("git_head"), versions=case.get("versions"))
            rid = src.run(dataset_id=did, method="cubical", coeff_field=2, max_dim=n, params=params,
                          preprocessing=e.get("construction"), script=SCRIPT,
                          command=case.get("command", ""), tier="X", seed=str(22 + n),
                          wall_sec=t.get("wall_sec"), peak_mb=case.get("peak_rss_mb"))
            src.betti(rid, C.betti_map(t.get(betti_key)), C.betti_map(expected) if expected else None)
            if periodic:
                for d, v in enumerate(t.get("infinite_bars_per_dim") or []):
                    src.stat(rid, f"n_infinite_bars_H{d}", v)
                src.control(rid, "known_answer",
                            f"F2 {name}: periodic cubical Betti numbers equal binomial({n}, k) = {t.get('expected')}",
                            t.get("pass"), detail=f"observed {t.get(betti_key)}")
            else:
                src.control(rid, "negative",
                            f"F2 {name}: the NON-periodic complex on the same field is contractible",
                            t.get("nonperiodic_control_betti") == [1] + [0] * n,
                            detail=f"observed {t.get('nonperiodic_control_betti')}")
        src.finding(dataset_id=did, tier="X",
                    claim=f"F2 {name}: PeriodicCubicalComplex returns Betti {t.get(betti_key)} against the "
                          f"pre-stated binomial({n}, k) = {t.get('expected')}, and the non-periodic control "
                          f"on the same field returns {t.get('nonperiodic_control_betti')}",
                    verdict="recovered" if t.get("pass") else "failed",
                    caveat="one field realisation per torus", reference=f"{DIR}/report.json")

    # ---- F3: the CMB p-value calibration ---------------------------------
    case, e = cases["F3"], exp_by_id["F3"]
    did = src.dataset(**_dataset_kwargs("F3", e))
    f3rep = next((r for r in rep.get("table", []) if r.get("id") == "F3"), {})
    params = dict(base_params, case="F3", function=case.get("function"), n_maps=case.get("n_maps"),
                  expected=e.get("expected"),
                  topology_diagnostic_betti_full_sky=case.get("topology_diagnostic"),
                  branches="two p-value branches: 'hartlap' (chi2 with the Hartlap-corrected covariance) "
                           "and 'empirical_rank'")
    rid = src.run(dataset_id=did, method="lower_star_graph", coeff_field=2, max_dim=1, params=params,
                  preprocessing="cmb_tda.build_topology(full-sky mask, nside 64) + "
                                "betti_curves_from_topology(sublevel) + coarse_stats",
                  script=SCRIPT, command=f"{SCRIPT} --case F3", tier="X", seed="30000..30199")
    for branch in ("hartlap", "empirical_rank"):
        blk = f3rep.get(branch) or {}
        for stat_name, vals in blk.items():
            src.stat(rid, f"{branch}_{stat_name}_ks_p", vals.get("ks_p"))
            src.stat(rid, f"{branch}_{stat_name}_n_below_0.05", vals.get("n_below_0.05"))
            if "ks_stat" in vals:
                src.stat(rid, f"{branch}_{stat_name}_ks_stat", vals.get("ks_stat"))
    src.control(rid, "null_calibration",
                "F3 (TDA-7): p-values of an independent realisation against a same-C_ell ensemble are "
                "uniform (KS p >= 0.01 and at most 9 of 100 below 0.05) for b0, b1 and chi",
                f3rep.get("pass"),
                detail="hartlap branch: " + str({k: {kk: vv for kk, vv in v.items() if kk in
                                                     ("ks_p", "n_below_0.05", "gate_pass")}
                                                 for k, v in (f3rep.get("hartlap") or {}).items()}))
    src.finding(run_id=rid, dataset_id=did, tier="X",
                claim="F3: the chi2/Hartlap p-value branch is not uniform for b1 "
                      f"(KS p = {(f3rep.get('hartlap') or {}).get('b1', {}).get('ks_p')}, "
                      f"{(f3rep.get('hartlap') or {}).get('b1', {}).get('n_below_0.05')} of 100 below 0.05), "
                      "while the empirical-rank branch on the same maps is consistent with uniform "
                      f"(b1 KS p = {(f3rep.get('empirical_rank') or {}).get('b1', {}).get('ks_p')})",
                verdict="failed",
                caveat="the full-sky topology diagnostic in the same case reports Betti "
                       f"{case.get('topology_diagnostic')} for a closed sphere, which is not (1,0,1); "
                       "this is recorded as the suite found it",
                reference=f"{DIR}/report.json")
    src.note_discrepancy(
        "simple_suite F3",
        f"the full-sky topology diagnostic records Betti {case.get('topology_diagnostic')} where a closed "
        "2-sphere has (1,0,1); the source reports this as a defect of the sublevel construction under test, "
        "and the number is ingested as recorded")

    for dev in rep.get("deviations_from_expectations") or []:
        src.skip("simple_suite deviation from the pre-registration",
                 "recorded verbatim from report.json deviations_from_expectations", detail=dev)
    for f in rep.get("other_failures") or []:
        src.skip(f"simple_suite {f.get('case')}", f"{f.get('result')}: {f.get('cause')}")
    src.skip("simple_suite per-null persistence diagrams",
             "the 50 Poisson seeds per P9 case are summarised by quantiles in results.json; no per-seed "
             "diagram was saved, so only the quantile statistics could be carried")
    src.skip("simple_suite full persistence diagrams",
             "results.json stores only the 6 longest bars per dimension "
             "(top_bars_birth_death_capped); the full diagrams were not committed")
    return src.report()
