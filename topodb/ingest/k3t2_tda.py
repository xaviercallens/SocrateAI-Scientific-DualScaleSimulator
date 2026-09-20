"""Backfill source (c): chain-level homology of the Kummer/K3 construction, and
the Part B point-sample study.

Source worktree : dualscale-wt-tdak3t2, branch loop/tda-k3t2
Source files    : audit/tda_validation/k3t2/{expectations.json, expectations_A5b.json,
                  results.json, report.json}
Pre-registration: expectations.json states `written_before_computation: true`,
                  date 2026-09-19 (report.json cites commit 529bf2b); the A5b
                  predictions are in expectations_A5b.json, written after A5 and
                  before A5b was run. Every `expected` vector carried into the DB
                  comes from a partA/partB table row whose source marks it as the
                  pre-stated value; where a row records
                  "not pre-specified (chi=8 only)" NO expectation is stored, so
                  `matches` stays NULL for it.
Tiers           : report.json states "Part A tier B (exact arithmetic mod p with
                  negative controls) with one tier-L input in A4; Part B tier X."
                  Part A rows carry the row's own `tier` field (B); Part B rows
                  carry X, as the task requires for sampled results.

The honest headline this module carries into the DB: Part B's K3 runs settle on a
beta_2 PLATEAU OF 27, not the expected 22, in all three seeds. That is ingested
as an observed Betti number of 27 against a pre-stated expectation of 22, so the
`matches` flag is 0 and the finding's verdict is 'failed', not 'recovered'.
"""
from __future__ import annotations

import re

from . import _common as C

DIR = "audit/tda_validation/k3t2"
FIELD = {"F2": 2, "F3": 3, "F5": 5}

# Which script produced each Part A family, from report.json `commands`.
PARTA_SCRIPT = {
    "A1": f"{DIR}/A1_tori_and_engine_controls.py",
    "A2": f"{DIR}/A2_orbifold.py",
    "A3": f"{DIR}/A3_orbifold_x_T2.py",
    "A4": f"{DIR}/A4_resolved_K3xT2.py",
    "A5b": f"{DIR}/A5b_kummer_code_F2.py",
    "A5": f"{DIR}/A5_negative_controls.py",
}
PARTB_METHOD = {"alpha": "alpha", "rips": "rips", "witness": "witness"}


def _slug(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", s.lower()).strip("_")


def _family(complex_name: str) -> str:
    for k in ("A5b", "A1", "A2", "A3", "A4", "A5"):
        if complex_name.startswith(k + " "):
            return k
    return "A?"


def ingest(db) -> dict:
    src = C.Source(db, "c_k3t2_chain_and_sampling", "tdak3t2")
    res = C.load(src.rel(DIR, "results.json"))
    rep = C.load(src.rel(DIR, "report.json"))
    exp = C.load(src.rel(DIR, "expectations.json"))
    runtime = res.get("partA_runtime_memory") or {}
    base = {
        "branch_recorded_in_source": res.get("branch"),
        "hardware_recorded_in_source": res.get("hardware"),
        "expectations_pre_registered": "audit/tda_validation/k3t2/expectations.json "
                                       f"(written_before_computation={exp['meta'].get('written_before_computation')}, "
                                       f"{exp['meta'].get('date')}); report.json cites commit 529bf2b",
        "kunneth_note": exp["meta"].get("note"),
        "tier_note": rep.get("tiers"),
    }

    # ================================================================ Part A
    seen_datasets: set[str] = set()
    for row in res.get("partA_table") or []:
        name = row["complex"]
        fam = _family(name)
        did = f"mathematics/k3t2/partA/{_slug(name)}"
        if did not in seen_datasets:
            src.dataset(id=did, domain="mathematics",
                        title=f"{name} (cell complex, {row.get('n_cells')} cells)",
                        source="cell complex constructed by the Part A scripts; expected homology "
                               "pre-registered in expectations.json (tier L input, theory)",
                        provenance="synthetic_control", n_objects=row.get("n_cells"),
                        ambient_dim=max(0, len(row.get("cells_per_dim") or [1]) - 1),
                        units="cells", local_path=f"{DIR}/results.json",
                        notes="exact cellular chain complex; no filtration, so this run carries Betti "
                              "numbers and no bars")
            seen_datasets.add(did)
        field = FIELD[row["field"]]
        expected_raw = row.get("expected")
        pre_stated = isinstance(expected_raw, list)
        params = dict(base, part="A", complex=name, family=fam, field=row["field"],
                      cells_per_dim=row.get("cells_per_dim"), n_cells=row.get("n_cells"),
                      ranks_of_boundary_maps=row.get("ranks_d_k"), chi=row.get("chi"),
                      n_fixed_cells=row.get("n_fixed_cells"),
                      lefschetz_number_chain_level=row.get("lefschetz_number_chain_level"),
                      resolved_bits=row.get("resolved_bits"),
                      engine1_agrees=row.get("engine1_agrees"), engine2_agrees=row.get("engine2_agrees"),
                      gudhi_agrees=row.get("gudhi_agrees"),
                      expected_verbatim=expected_raw,
                      tier_verbatim_in_source=row.get("tier"),
                      expectation_status=("pre-stated in expectations.json" if pre_stated else
                                          f"the source records expected = {expected_raw!r}, so no "
                                          "expectation is stored and `matches` stays NULL"))
        # Some A4 rows record tier as "B (+1 tier-L input: class of phi_2)". The schema takes a
        # single letter, so the leading letter is used and the full string is kept in params.
        tier = str(row.get("tier", "B")).strip()[:1] or "B"
        rid = src.run(dataset_id=did, method="chain_complex", coeff_field=field,
                      max_dim=max(0, len(row.get("computed") or [1]) - 1), params=params,
                      preprocessing="exact sparse elimination over F_p (two independent rank engines)",
                      script=PARTA_SCRIPT.get(fam, f"{DIR}/make_results.py"),
                      command=f"{PARTA_SCRIPT.get(fam, 'make_results.py')} ({name}, {row['field']})",
                      tier=tier, seed=None,
                      wall_sec=row.get("sec") or (runtime.get(fam) or [None])[0],
                      peak_mb=(runtime.get(fam) or [None, None])[1])
        src.betti(rid, C.betti_map(row.get("computed")),
                  C.betti_map(expected_raw) if pre_stated else None)
        src.stat(rid, "euler_characteristic", row.get("chi"))
        src.stat(rid, "n_cells", row.get("n_cells"))
        src.stat(rid, "n_fixed_cells", row.get("n_fixed_cells"))
        src.stat(rid, "lefschetz_number_chain_level", row.get("lefschetz_number_chain_level"))
        if pre_stated:
            src.control(rid, "known_answer",
                        f"{name} over {row['field']}: computed homology equals the pre-stated {expected_raw}",
                        row.get("PASS"), detail=f"computed {row.get('computed')}")
        else:
            src.control(rid, "known_answer",
                        f"{name} over {row['field']}: only the pre-stated Euler characteristic "
                        f"chi = {row.get('chi')} was gated (the Betti vector was not pre-specified)",
                        row.get("PASS"), detail=f"computed {row.get('computed')}, chi {row.get('chi')}")
        for engine, agrees in (("second independent rank engine", row.get("engine2_agrees")),
                               ("first independent rank engine", row.get("engine1_agrees")),
                               ("gudhi periodic cubical complex", row.get("gudhi_agrees"))):
            if agrees is not None:
                src.control(rid, "known_answer", f"cross-check: {engine} agrees", bool(agrees),
                            detail=f"wall_sec engine2 {row.get('sec_engine2')}")
    src.finding(tier="B", verdict="recovered",
                claim="Part A: all 119 table rows agree with their pre-registered expectation "
                      "(or, for the three F2 orbifold rows, with the pre-stated chi = 8 alone); the "
                      "resolved Kummer complex gives K3 = (1,0,22,0,1) over F2, F3 and F5 and "
                      "(K3 x T^2) = (1,2,23,44,23,2,1) over the same three fields, at chain level, "
                      "with Kunneth used only to state the expectation and never inside a computation",
                caveat="report.json states one tier-L input in A4: that the Kummer disc bundle "
                       "D(O(-2)) is the mapping cylinder of RP^3 -> S^2 whose pull-back of [S^2] "
                       "generates H^2(RP^3;Z) = Z/2. Over odd primes that input is not tested (the "
                       "wrong-class control also gives (1,0,22,0,1)); only F2 discriminates it, where "
                       "the wrong class gives (1,0,27,5,1)",
                reference=f"{DIR}/report.json")

    # ================================================================ Part B
    spaces = {
        "T2": dict(title="flat T^2, product of unit circles in R^4", ambient_dim=4),
        "T3": dict(title="flat T^3, product of unit circles in R^6", ambient_dim=6),
        "T4": dict(title="flat T^4, product of unit circles in R^8", ambient_dim=8),
        "orbifold": dict(title="T^4/Z2 in the Z2-invariant R^14 embedding", ambient_dim=14),
        "K3": dict(title="K3 as the Fermat quartic, Hermitian projector embedding in R^16", ambient_dim=16),
        "K3xT2": dict(title="K3 x T^2: K3 projector (R^16) x flat T^2 of radius 0.5 (R^4)", ambient_dim=20),
        "null4m": dict(title="null: uniform 4-cube in R^16, density-matched to the K3 sample", ambient_dim=16),
        "quadric": dict(title="positive control: Fermat quadric x^2+y^2+z^2+w^2=0 (S^2 x S^2), "
                              "same projector embedding as K3", ambient_dim=16),
    }
    for key, meta in spaces.items():
        src.dataset(id=f"mathematics/k3t2/partB/{key}", domain="mathematics", title=meta["title"],
                    source="point sample generated by B_driver.py / B_explore.py; expected Betti "
                           "vector pre-registered in expectations.json partB",
                    provenance="synthetic_control", ambient_dim=meta["ambient_dim"],
                    local_path=f"{DIR}/B_summary.json",
                    notes="Part B is tier X (numerics on a finite sample), unlike the exact Part A "
                          "complexes of the same spaces")

    for row in res.get("partB_table") or []:
        space = row["space"]
        did = f"mathematics/k3t2/partB/{space}"
        at = row.get("at_largest_N") or {}
        params = dict(base, part="B", study=row.get("study"), space=space,
                      embedding=row.get("embedding"), tau=row.get("tau"),
                      sampling=row.get("sampling"),
                      expected_betti_pre_registered=row.get("expected"),
                      criterion=exp["partB"]["criterion_recovered"]["definition"],
                      window_ratio_R=exp["partB"]["criterion_recovered"]["R"],
                      N_min_all_3_seeds_F3=row.get("N_min_all_3_seeds_F3"),
                      largest_N_completed=row.get("largest_N_completed"),
                      stopped_because=row.get("stopped_because"),
                      n_simplices_at_largest_N=at.get("n_simplices"))
        rid = src.run(dataset_id=did, method=PARTB_METHOD[row["method"]], coeff_field=3,
                      max_dim=max(0, len(row.get("expected") or [1]) - 1), params=params,
                      preprocessing=row.get("sampling"), script=f"{DIR}/B_driver.py",
                      command=f"{DIR}/B_driver.py --study {row.get('study')}", tier=row.get("tier", "X"),
                      seed="0,1,2", wall_sec=at.get("seconds"), peak_mb=at.get("maxrss_MB"))
        # The source records a window criterion, not an observed Betti vector, for these rows,
        # so no Betti numbers are stored here; only the statistics it reports.
        src.stat(rid, "N_min_all_3_seeds_F3", row.get("N_min_all_3_seeds_F3"))
        src.stat(rid, "largest_N_completed", row.get("largest_N_completed"))
        src.stat(rid, "best_window_ratio_F3_at_largest_N", at.get("best_window_ratio_F3"))
        src.stat(rid, "n_simplices_at_largest_N", at.get("n_simplices"))
        src.control(rid, "known_answer",
                    f"{row.get('study')}: the pre-registered window criterion (the full expected Betti "
                    f"vector {row.get('expected')} holds on a window with e2 >= 1.5 e1, all 3 seeds)",
                    row.get("recovered"),
                    detail=f"N_min = {row.get('N_min_all_3_seeds_F3')}, largest N completed "
                           f"{row.get('largest_N_completed')}, stopped because: {row.get('stopped_because')}")
        if not row.get("recovered"):
            src.skip(f"k3t2 Part B {row.get('study')} observed Betti vector",
                     "the run did not satisfy the pre-registered window criterion and the source "
                     "records no single observed Betti vector for it, only the window statistics",
                     stopped_because=row.get("stopped_because"))

    # -- the K3 beta_2 = 27 plateau (the headline Part B result) ------------
    for seed_key, blk in (res.get("partB_K3") or {}).items():
        if seed_key == "N3000_tau0.9":
            continue
        plateau = blk.get("beta2_plateau(value,e1,e2,ratio)") or [None] * 4
        tail = blk.get("F2_beta_curve_tail")
        params = dict(base, part="B", study="partB_K3_plateau", seed_key=seed_key,
                      N=blk.get("N"), tau=blk.get("tau"),
                      expected_betti_pre_registered=exp["partB"]["B3_K3_Fermat_quartic"]["expected_H0_H2"],
                      beta2_plateau_value=plateau[0], plateau_window=[plateau[1], plateau[2]],
                      plateau_window_ratio=plateau[3],
                      beta2_at=blk.get("beta2_at"),
                      beta2_ever_equals_22_on_grid=blk.get("beta2_ever_equals_22_on_grid"),
                      F2_beta_curve_tail=tail,
                      sampling="farthest-point subsample" if seed_key.endswith("fps") else "iid uniform")
        rid = src.run(dataset_id="mathematics/k3t2/partB/K3", method="rips", coeff_field=2,
                      max_dim=2, params=params,
                      preprocessing=f"N={blk.get('N')} points, tau={blk.get('tau')}, "
                                    f"{'farthest-point subsample' if seed_key.endswith('fps') else 'iid uniform'}",
                      script=f"{DIR}/B_summarize.py",
                      command=f"{DIR}/B_driver.py K3 ({seed_key})", tier="X", seed=seed_key)
        # The plateau IS an observed Betti reading, so it is stored against the pre-stated (1,0,22).
        if plateau[0] is not None:
            src.betti(rid, {0: 1, 1: 0, 2: int(plateau[0])}, {0: 1, 1: 0, 2: 22})
        src.stat(rid, "beta2_plateau_value", plateau[0])
        src.stat(rid, "beta2_plateau_window_ratio", plateau[3])
        src.stat(rid, "n_H2_bars", blk.get("n_H2_bars"))
        src.stat(rid, "n_H2_bars_alive_at_tau", blk.get("n_H2_bars_alive_at_tau"))
        for eps, v in (blk.get("beta2_at") or {}).items():
            src.stat(rid, f"beta2_at_eps_{eps}", v)
        src.control(rid, "known_answer",
                    "K3 point sample: beta_2 equals the pre-stated 22 somewhere on the scanned "
                    "epsilon grid", not blk.get("beta2_ever_equals_22_on_grid") is False and
                    blk.get("beta2_ever_equals_22_on_grid"),
                    detail=f"beta2_ever_equals_22_on_grid = {blk.get('beta2_ever_equals_22_on_grid')}; "
                           f"the plateau value is {plateau[0]} on [{plateau[1]}, {plateau[2]}]")
    lw = res.get("partB_K3_longest_beta2_ge1_window") or {}
    src.finding(dataset_id="mathematics/k3t2/partB/K3", tier="X", verdict="failed",
                claim="K3 point samples (Fermat quartic, projector embedding, N=4000, tau=0.8): beta_0 = 1 "
                      "and beta_1 = 0 are recovered, but beta_2 settles on a plateau of 27 in all three "
                      "iid seeds (28 in the farthest-point-subsampled seed) and never equals the "
                      "pre-registered 22 anywhere on the scanned epsilon grid",
                caveat="the same-dimension positive control (the Fermat quadric S^2 x S^2, b2 = 2) never "
                       "left the noise regime at the sizes that fit in 8 GiB, so it neither reproduces "
                       "nor rules out the over-count; report.json states the origin of the 27 is "
                       "undetermined by these runs. The density-matched null does NOT produce the "
                       f"plateau (beta_2 = 0 on [0.64, 0.8]), but its longest beta_2>=1 windows "
                       f"(ratio ~1.7-2.2) are comparable to the K3 runs' ({lw})",
                reference=f"{DIR}/report.json")
    src.note_discrepancy(
        "k3t2 Part B, K3",
        "the exact Part A chain complex gives b2(K3) = 22 over F2/F3/F5, while the Part B point sample "
        "of the same space reads a stable b2 plateau of 27 (28 under farthest-point subsampling). Both "
        "numbers are ingested: the 22 as a tier-B chain_complex run, the 27 as a tier-X rips run whose "
        "betti row carries expected=22 and matches=0")

    # -- the density-matched null false-positive check ---------------------
    for key, blk in (res.get("partB_null_false_positive_check") or {}).items():
        n_str, seed = key.split("_seed")
        params = dict(base, part="B", study="partB_null_false_positive_check", N=int(n_str[1:]),
                      seed=seed, tau=0.8,
                      expected_betti_pre_registered=exp["partB"]["B4_null_random_4d"]["expected_H0_H2"],
                      false_positive_rule=exp["partB"]["B4_null_random_4d"]["false_positive_rule"])
        rid = src.run(dataset_id="mathematics/k3t2/partB/null4m", method="rips", coeff_field=2, max_dim=2,
                      params=params,
                      preprocessing="uniform 4-cube in R^16, median 10-NN distance matched to the K3 sample",
                      script=f"{DIR}/B_summarize.py", command=f"{DIR}/B_driver.py null4m ({key})",
                      tier="X", seed=seed)
        src.stat(rid, "max_window_ratio_beta2_ge1", blk.get("max_window_ratio_beta2_ge1"))
        src.stat(rid, "beta2_at_eps_0.7", blk.get("beta2_at_0.7"))
        triggered = blk.get("false_positive_rule_(>=1.5)")
        src.control(rid, "null_calibration",
                    "the pre-registered rule 'any beta_2>=1 window of ratio >= 1.5 is a false positive' "
                    "is NOT triggered on the density-matched null",
                    (not triggered) if triggered is not None else None,
                    detail=f"max window ratio {blk.get('max_window_ratio_beta2_ge1')}, "
                           f"beta_2 at eps=0.7 is {blk.get('beta2_at_0.7')}")
    n_trig = sum(1 for b in (res.get("partB_null_false_positive_check") or {}).values()
                 if b.get("false_positive_rule_(>=1.5)"))
    n_tot = len(res.get("partB_null_false_positive_check") or {})
    src.finding(dataset_id="mathematics/k3t2/partB/null4m", tier="X", verdict="artefact",
                claim=f"the density-matched uniform null triggers the pre-registered "
                      f"'beta_2 >= 1 window of ratio >= 1.5' rule in {n_trig} of {n_tot} runs "
                      "(every run at N >= 1000), through overlapping short noise bars",
                caveat="so 'some persistent H2' in this setting is not evidence of anything; what the "
                       "null does not produce is a plateau (its beta_2 is 0 on [0.64, 0.8])",
                reference=f"{DIR}/report.json")

    for idx, blk in (res.get("partB_failures_recorded") or {}).items():
        src.skip(f"k3t2 Part B run index {idx} ({blk.get('space')}, N={blk.get('N')}, tau={blk.get('tau')})",
                 f"the source records reason={blk.get('reason')!r}: the run produced no barcode")
    fit = res.get("partB_scaling_fit") or {}
    src.finding(dataset_id="mathematics/k3t2/partB/T4", tier="X", verdict="inconclusive",
                claim=f"a two-point fit N_min = {fit.get('A')} * {fit.get('B_per_dimension')}^d over the "
                      f"recovered flat tori ({fit.get('points')}) extrapolates to about "
                      f"{fit.get('extrapolated_N_min_T4')} points for T^4",
                caveat=str(fit.get("caveat")) + f"; T^4 status: {fit.get('T4_status')}",
                reference=f"{DIR}/report.json")
    src.finding(dataset_id="mathematics/k3t2/partB/K3xT2", tier="X", verdict="null",
                claim="K3 x T^2 point samples at N = 4000, 8000, 16000 and 32000 (99M simplices, 6.7 GB) "
                      "show beta_1 in the thousands at mid scales and nothing resembling the pre-stated "
                      "(1, 2, 23); N = 64000 timed out",
                caveat="a compute bound, not a statement that no method could recover it; report.json: "
                       "'no claim that K3 x T2 topology is visible in cosmological point data is "
                       "supported by point-sample TDA'",
                reference=f"{DIR}/report.json")

    src.skip("k3t2 full barcodes (B_*.json, B_explore.json)",
             "the per-run JSONs record Betti CURVES on an epsilon grid rather than birth/death pairs, "
             "so no persistence diagram could be stored for the Part B runs; the curves are carried as "
             "beta2_at_eps_* statistics instead")
    src.skip("k3t2 Part A bars",
             "Part A computes exact cellular homology of an unfiltered chain complex; there is no "
             "filtration and therefore no barcode to store")
    for lim in rep.get("limitations") or []:
        src.skip("k3t2 stated limitation", lim)
    return src.report()
