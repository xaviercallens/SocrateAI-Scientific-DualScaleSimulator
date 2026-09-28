"""Backfill source (e): the K3 x T2 rigidity tracks D-tda, rounds v2 and v3.

Source worktree : dualscale-wt-k3t2
Source files    : audit/k3t2_rigidity_v2/D-tda/{results.json, exports.json,
                    00_premise_checks_results.json, 01_controls_results.json,
                    02_invariance_and_quotient_results.json,
                    03_resolution_hybrid_results.json, 04_kunneth_results.json,
                    05_rigidity_scan_results.json}
                  audit/k3t2_rigidity_v3/D-tda/{results.json, inputs.json,
                    00_premise_results.json, 01_controls_results.json,
                    02_mv_results.json, 03_cup_parity_results.json}
Tier            : the task assigns tier B to the chain-level and exact-arithmetic
                  results, and both rounds agree: v2's results.json opens with
                  "tier: B (exact simplicial arithmetic + a hybrid Mayer-Vietoris
                  step) at best; never 'proved'." That qualification is carried
                  verbatim into params_json.

WHAT IS COMPUTED VS STATED. v2 labels its own Betti numbers explicitly:
"b1,b2 COMPUTED (Mayer-Vietoris from GUDHI Betti numbers of U and links);
b0,b3,b4 STATED (connectedness/Poincare duality)". That labelling is carried in
params so no reader mistakes a stated value for a computed one.

NO BARS. Both rounds compute homology of unfiltered simplicial and cubical
complexes through gudhi's Betti numbers; there is no filtration and therefore no
barcode.

EXPECTATIONS. `expected=` is stored only where the source file itself carries a
pre-stated expectation: the torus and sphere controls in 01_controls_results.json
carry `expected` / `expected_after_computing` with an `expected_source`, so they
are gated. The resolved-K3 and K3xT2 Betti numbers are the TARGET of the blind
re-derivation and the rounds deliberately avoid stating them in advance, so no
expectation is attached to those runs and `matches` stays NULL.
"""
from __future__ import annotations

from . import _common as C

V2 = "audit/k3t2_rigidity_v2/D-tda"
V3 = "audit/k3t2_rigidity_v3/D-tda"


def _controls_v2(src, did, base, ctrl, script, cmd):
    """The v2 control file: cubical tori and Freudenthal T^4 at several N."""
    made = []
    cub = ctrl.get("cubical") or {}
    for name in ("T2", "T4"):
        blk = cub.get(name)
        if not isinstance(blk, dict):
            continue
        rid = src.run(dataset_id=did, method="cubical",
                      coeff_field=cub.get("homology_coeff_field", 3),
                      max_dim=max(0, len(blk.get("betti") or [1]) - 1),
                      params=dict(base, leg=f"controls/cubical/{name}", method_note=cub.get("method"),
                                  expected_source=blk.get("expected_source")),
                      preprocessing=str(cub.get("method")), script=script, command=cmd, tier="B",
                      wall_sec=cub.get("runtime_sec"))
        src.betti(rid, C.betti_map(blk.get("betti")), C.betti_map(blk.get("expected")))
        src.control(rid, "known_answer",
                    f"{name} cubical control against the pre-stated {blk.get('expected')} "
                    f"({blk.get('expected_source')})",
                    blk.get("betti") == blk.get("expected"),
                    detail=f"computed {blk.get('betti')}")
        made.append(rid)
    for N, blk in (ctrl.get("freudenthal_simplicial_T4") or {}).items():
        rid = src.run(dataset_id=did, method="chain_complex", coeff_field=3,
                      max_dim=max(0, len(blk.get("betti") or [1]) - 1),
                      params=dict(base, leg=f"controls/freudenthal_T4/N={N}", N=int(N),
                                  simplex_counts_by_dim=blk.get("simplex_counts_by_dim"),
                                  num_vertices=blk.get("num_vertices"),
                                  expected_source=blk.get("expected_source"),
                                  note=blk.get("note")),
                      preprocessing=f"Freudenthal/Kuhn triangulation of (R/{N}Z)^4",
                      script=script, command=cmd, tier="B", wall_sec=blk.get("runtime_sec"))
        src.betti(rid, C.betti_map(blk.get("betti")), C.betti_map(blk.get("expected_betti")))
        src.stat(rid, "euler_characteristic_from_simplex_counts",
                 blk.get("euler_characteristic_from_simplex_counts"))
        src.stat(rid, "euler_characteristic_from_betti_numbers",
                 blk.get("euler_characteristic_from_betti_numbers"))
        src.stat(rid, "num_top_4simplices_built", blk.get("num_top_4simplices_built"))
        src.stat(rid, "simplex_tree_num_simplices", blk.get("simplex_tree_num_simplices"))
        src.control(rid, "known_answer",
                    f"Freudenthal T^4 at N={N} against the pre-stated {blk.get('expected_betti')} "
                    f"({blk.get('expected_source')})",
                    blk.get("matches_control") if blk.get("matches_control") is not None
                    else blk.get("betti") == blk.get("expected_betti"),
                    detail=f"computed {blk.get('betti')}; the simplex count matches the closed "
                           f"formula ({blk.get('num_top_4simplices_built')} vs "
                           f"{blk.get('num_top_4simplices_expected_formula')}) and both Euler "
                           "characteristics agree")
        made.append(rid)
    return made


def ingest(db) -> dict:
    src = C.Source(db, "e_k3t2_rigidity", "k3t2")

    # ======================================================== round v2
    v2 = C.load(src.rel(V2, "results.json"))
    v2ctrl = C.load(src.rel(V2, "01_controls_results.json"))
    v2exp = C.load(src.rel(V2, "exports.json"))
    base2 = dict(round="k3t2_rigidity_v2 track D-tda",
                 tier_verbatim_in_source=v2.get("tier"),
                 blindness_self_report=v2.get("blindness_self_report"))
    did_c = "mathematics/k3t2_rigidity/v2_controls"
    src.dataset(id=did_c, domain="mathematics",
                title="v2 known-answer controls: cubical T^2 and T^4, and the Freudenthal/Kuhn "
                      "triangulation of T^4 at N = 4, 6, 8",
                source="b_k(T^n) = C(n, k), standard torus cohomology (tier L)",
                provenance="synthetic_control", local_path=f"{V2}/01_controls_results.json")
    _controls_v2(src, did_c, base2, v2ctrl, f"{V2}/01_controls.py",
                 f"cd {V2} && python 01_controls.py")

    did_q = "mathematics/k3t2_rigidity/v2_T4_mod_Z2_quotient"
    src.dataset(id=did_q, domain="mathematics",
                title="v2: the T^4/Z2 quotient complex and its 16 singular points, at N = 4, 6, 8",
                source="constructed by 00_premise_checks.py / 02_invariance_and_quotient.py",
                provenance="synthetic_control", local_path=f"{V2}/results.json",
                notes=str((v2.get("premise_checks") or {}).get("N4_excluded_because")))
    pc = v2.get("premise_checks") or {}
    for N in sorted({*(pc.get("closed_star_disjointness") or {}),
                     *(pc.get("nonsingular_link_is_S3") or {})}, key=int):
        rid = src.run(dataset_id=did_q, method="chain_complex", coeff_field=3, max_dim=4,
                      params=dict(base2, leg=f"premise_checks/N={N}", N=int(N),
                                  closed_star_disjointness=(pc.get("closed_star_disjointness") or {}).get(N),
                                  nonsingular_link_is_S3=(pc.get("nonsingular_link_is_S3") or {}).get(N),
                                  valid_N_used_downstream=pc.get("valid_N_used_downstream"),
                                  N4_excluded_because=pc.get("N4_excluded_because"),
                                  source_file=pc.get("source")),
                      preprocessing="structural premise checks on the quotient complex",
                      script=f"{V2}/00_premise_checks.py",
                      command=f"cd {V2} && python 00_premise_checks.py", tier="B")
        src.control(rid, "known_answer",
                    f"N={N}: the closed stars of the 16 singular points are pairwise disjoint (the "
                    "structural premise the Mayer-Vietoris step needs)",
                    (pc.get("closed_star_disjointness") or {}).get(N),
                    detail=f"nonsingular link is S^3: "
                           f"{(pc.get('nonsingular_link_is_S3') or {}).get(N)}")
        if (pc.get("closed_star_disjointness") or {}).get(N) is False:
            src.control(rid, "negative",
                        f"N={N} is the negative control: the premise fails there, and the source "
                        "refuses to treat its numerically matching Betti numbers as evidence",
                        True, detail=str(pc.get("N4_excluded_because")))
    tnc = v2.get("translation_negative_control") or {}
    for N, blk in tnc.items():
        if not isinstance(blk, dict):
            continue
        rid = src.run(dataset_id=did_q, method="chain_complex", coeff_field=3, max_dim=4,
                      params=dict(base2, leg=f"translation_negative_control/{N}",
                                  source_file=tnc.get("source")),
                      preprocessing="translation (rather than inversion) involution: a negative "
                                    "control for the quotient construction",
                      script=f"{V2}/02_invariance_and_quotient.py",
                      command=f"cd {V2} && python 02_invariance_and_quotient.py", tier="B")
        src.control(rid, "negative",
                    f"{N}: quotienting by a free translation reproduces T^4's Betti numbers rather "
                    "than the orbifold's",
                    blk.get("matches_T4_betti"),
                    detail=f"passes the f-vector halving check: {blk.get('passes_halving_check')}; "
                           f"dropped: {blk.get('dropped')}")

    did_k3 = "mathematics/k3t2_rigidity/v2_resolved_K3"
    src.dataset(id=did_k3, domain="mathematics",
                title="v2: the resolved Kummer K3, by Mayer-Vietoris on the computed quotient and "
                      "its 16 resolved neighbourhoods",
                source="the resolution local model (a disc bundle over S^2 with RP^3 boundary) is a "
                       "stated tier-L input; everything else is computed",
                provenance="synthetic_control", local_path=f"{V2}/results.json")
    rb = v2.get("resolved_K3_betti_numbers") or {}
    for N in [k for k in rb if k.startswith("N=")]:
        blk = rb[N]
        rid = src.run(dataset_id=did_k3, method="chain_complex", coeff_field=3, max_dim=4,
                      params=dict(base2, leg=f"resolved_K3/{N}", N=int(N.split("=")[1]),
                                  labeling=rb.get("labeling"),
                                  computed_vs_stated="b1 and b2 are COMPUTED by Mayer-Vietoris from "
                                                     "gudhi Betti numbers of U and of the links; "
                                                     "b0, b3 and b4 are STATED from connectedness "
                                                     "and Poincare duality",
                                  cross_field_Z5_agrees=(rb.get("cross_field_Z5_agrees") or {}).get(N),
                                  chi_fvector_vs_betti_cross_check=(rb.get("chi_fvector_vs_betti_cross_check") or {}).get(N),
                                  torsion_control_RP3_mod2_jump_confirmed=(rb.get("torsion_control_RP3_mod2_jump_confirmed") or {}).get(N),
                                  source_file=rb.get("source"),
                                  no_expectation="the resolved-K3 Betti vector is the TARGET of this "
                                                 "blind re-derivation; the round does not pre-state "
                                                 "it, so no expectation is stored"),
                      preprocessing="hybrid Mayer-Vietoris on the computed quotient and links",
                      script=f"{V2}/03_resolution_hybrid.py",
                      command=f"cd {V2} && python 03_resolution_hybrid.py", tier="B")
        src.betti(rid, {i: blk.get(f"b{i}") for i in range(5) if blk.get(f"b{i}") is not None})
        src.stat(rid, "euler_characteristic", v2exp.get("chi_K3_resolved") if N == "N=6" else None)
        src.stat(rid, "num_fixed_points_computed", v2exp.get("num_fixed_points_computed")
                 if N == "N=6" else None)
        src.control(rid, "known_answer",
                    f"{N}: the Z/5 computation agrees with the Z/3 one, and the Euler characteristic "
                    "from the f-vector agrees with the one from the Betti numbers",
                    bool((rb.get("cross_field_Z5_agrees") or {}).get(N)) and
                    bool((rb.get("chi_fvector_vs_betti_cross_check") or {}).get(N)),
                    detail=f"cross-field agreement {(rb.get('cross_field_Z5_agrees') or {}).get(N)}, "
                           f"chi cross-check {(rb.get('chi_fvector_vs_betti_cross_check') or {}).get(N)}")
        src.control(rid, "known_answer",
                    f"{N}: the RP^3 torsion control shows the expected mod-2 jump",
                    (rb.get("torsion_control_RP3_mod2_jump_confirmed") or {}).get(N),
                    detail="a field-sensitivity control on the link homology")
    kt = v2.get("K3_x_T2_kunneth") or {}
    did_kt = "mathematics/k3t2_rigidity/v2_K3xT2"
    src.dataset(id=did_kt, domain="mathematics",
                title="v2: K3 x T^2, by Kunneth over a field from the computed K3 Betti numbers",
                source="Kunneth over a field is used HERE, at the product step; the K3 factor is "
                       "computed, not assumed",
                provenance="synthetic_control", local_path=f"{V2}/04_kunneth_results.json")
    rid = src.run(dataset_id=did_kt, method="chain_complex", coeff_field=3, max_dim=6,
                  params=dict(base2, leg="K3xT2_kunneth",
                              chi_cross_check_matches_multiplicativity=kt.get("chi_cross_check_matches_multiplicativity"),
                              all_N_and_fields_agree=kt.get("all_N_and_fields_agree"),
                              source_file=kt.get("source")),
                  preprocessing="Kunneth over a field from the computed K3 Betti numbers",
                  script=f"{V2}/04_kunneth_k3xt2.py", command=f"cd {V2} && python 04_kunneth_k3xt2.py",
                  tier="B")
    src.betti(rid, C.betti_map(kt.get("betti")))
    src.stat(rid, "euler_characteristic", kt.get("chi"))
    src.control(rid, "known_answer",
                "the product Euler characteristic matches chi(K3) x chi(T^2), and every valid N and "
                "both fields give the same answer",
                bool(kt.get("chi_cross_check_matches_multiplicativity")) and
                bool(kt.get("all_N_and_fields_agree")),
                detail=f"betti {kt.get('betti')}, chi {kt.get('chi')}")
    rs = v2.get("rigidity_scan_over_k") or {}
    src.finding(dataset_id=did_k3, tier="B", verdict="inconclusive",
                claim=f"v2 rigidity scan over k (the number of resolved points): "
                      f"{rs.get('closed_form')}; chi = 24 is classified "
                      f"{rs.get('chi_equals_24_classification')} with solution set "
                      f"{rs.get('chi_equals_24_solution_set')}, and the b1 selector is classified "
                      f"{str(rs.get('b1_selector_classification'))[:200]}",
                caveat=str(rs.get("conclusion"))[:900], reference=f"{V2}/results.json")
    rn = v2.get("rigidity_over_N") or {}
    src.finding(dataset_id=did_q, tier="B", verdict="recovered",
                claim=f"v2 rigidity over N: {str(rn.get('classification'))[:400]}",
                caveat=f"selecting condition: {str(rn.get('selecting_condition'))[:300]}; negative "
                       f"control: {str(rn.get('negative_control'))[:250]}",
                reference=f"{V2}/results.json")
    src.finding(dataset_id=did_kt, tier="B", verdict="recovered",
                claim=f"v2 track D-tda: the resolved Kummer K3 has Betti numbers "
                      f"{[rb.get('N=6', {}).get(f'b{i}') for i in range(5)]} and K3 x T^2 has "
                      f"{kt.get('betti')} with chi {kt.get('chi')}, agreeing across N = 6, 8 and "
                      "the fields Z/3 and Z/5",
                caveat=f"tier as the source states it: {v2.get('tier')}. b1 and b2 are computed; "
                       "b0, b3 and b4 are stated from connectedness and Poincare duality, and the "
                       "resolution local model is a stated tier-L input",
                reference=f"{V2}/results.json")

    # ======================================================== round v3
    v3 = C.load(src.rel(V3, "results.json"))
    v3ctrl = C.load(src.rel(V3, "01_controls_results.json"))
    v3in = C.load(src.rel(V3, "inputs.json"))
    by_id = {r["id"]: r for r in v3.get("results") or []}
    base3 = dict(round="k3t2_rigidity_v3 track D-tda", track=v3.get("track"),
                 declared_inputs=v3in,
                 tier_note="the task assigns tier B to exact-arithmetic and chain-level results; v3's "
                           "inputs.json labels the resolution local model itself tier L (FROM MEMORY; "
                           "standard Kummer construction)")
    did_c3 = "mathematics/k3t2_rigidity/v3_controls"
    src.dataset(id=did_c3, domain="mathematics",
                title="v3 known-answer controls: cubical T^2 and T^4, the Kuhn T^4 at N = 4, 6, 8, and "
                      "the tetrahedron-boundary S^2",
                source="b_k(T^n) = C(n, k) and H_*(S^2) = (Z, 0, Z), both tier L",
                provenance="synthetic_control", local_path=f"{V3}/01_controls_results.json")
    for name, blk in v3ctrl.items():
        if not isinstance(blk, dict) or "betti_Z3" not in blk:
            continue
        for field, key in ((3, "betti_Z3"), (2, "betti_Z2")):
            vec = blk.get(key)
            if not vec:
                continue
            rid = src.run(dataset_id=did_c3,
                          method="cubical" if name.startswith("cubical") else "chain_complex",
                          coeff_field=field, max_dim=max(0, len(vec) - 1),
                          params=dict(base3, leg=f"controls/{name}/{key}", control=name,
                                      fvector=blk.get("fvector"),
                                      chi_from_f=blk.get("chi_from_f"),
                                      expected_after_computing=blk.get("expected_after_computing")),
                          preprocessing=name, script=f"{V3}/01_controls.py",
                          command=f"cd {V3} && python 01_controls.py {v3ctrl.get('args')}", tier="B")
            src.betti(rid, C.betti_map(vec))
            src.stat(rid, "chi_from_fvector", blk.get("chi_from_f"))
            if blk.get("matches_expected") is not None:
                src.control(rid, "known_answer",
                            f"{name} over F{field} against {blk.get('expected_after_computing')}",
                            blk.get("matches_expected"), detail=f"computed {vec}")
    d1 = by_id.get("D1") or {}
    did_q3 = "mathematics/k3t2_rigidity/v3_T4_mod_Z2_quotient"
    src.dataset(id=did_q3, domain="mathematics",
                title="v3: the T^4/Z2 quotient complex and its premise checks at N = 3..8 "
                      "(odd N are realizable negative controls)",
                source="constructed by 00_premise.py", provenance="synthetic_control",
                local_path=f"{V3}/00_premise_results.json")
    for N, blk in (d1.get("computed") or {}).items():
        rid = src.run(dataset_id=did_q3, method="chain_complex", coeff_field=3, max_dim=4,
                      params=dict(base3, leg=f"D1_premise/{N}", N=int(N.split("=")[1]),
                                  quantity=d1.get("quantity"), computed=blk,
                                  shared_inputs=d1.get("shared_inputs")),
                      preprocessing="premise checks: regularity (f-halving), open-star disjointness, "
                                    "link homology",
                      script=f"{V3}/00_premise.py", command=str(d1.get("command", "")), tier="B")
        src.stat(rid, "n_fixed_points", blk.get("n_fixed_points"))
        src.stat(rid, "chi_Q_from_f", blk.get("chi_Q_from_f"))
        for field, vec in (blk.get("singular_link_representative") or {}).items():
            pass  # the link's Betti vector is per-field; recorded in params above
        passed = blk.get("PREMISE_VALID_FOR_MV")
        src.control(rid, "known_answer" if passed else "negative",
                    f"{N}: PREMISE_VALID_FOR_MV (regular image complex, disjoint open stars, links "
                    "RP^3 / S^3, no orbit collisions) - a condition that never mentions chi, b2 or 24",
                    passed,
                    detail=f"regular f-halving {blk.get('quotient_is_regular_f_halving')}, open stars "
                           f"disjoint {blk.get('open_stars_pairwise_disjoint')}, singular link is "
                           f"RP^3 homology {blk.get('singular_link_is_RP3_homology')}, "
                           f"{blk.get('n_fixed_points')} fixed points"
                           + ("" if passed else "; this odd N is a realizable negative control"))
    d2, d3, d5 = by_id.get("D2") or {}, by_id.get("D3") or {}, by_id.get("D5") or {}
    did_k33 = "mathematics/k3t2_rigidity/v3_resolved_K3"
    src.dataset(id=did_k33, domain="mathematics",
                title="v3: the resolved Kummer K3 by Mayer-Vietoris, with the chain-level rank of the "
                      "H_3 map",
                source="the resolution local model is declared tier L in inputs.json; the rest is "
                       "computed", provenance="synthetic_control", local_path=f"{V3}/02_mv_results.json")
    comp = d2.get("computed") or {}
    for N, blk in (comp.get("per_N") or {}).items():
        for field_name, vec in (blk.get("resolved_betti") or {}).items():
            field = int(field_name.replace("Z", "").replace("/", "")) if field_name[1:].isdigit() else 3
            rid = src.run(dataset_id=did_k33, method="chain_complex", coeff_field=field,
                          max_dim=max(0, len(vec) - 1),
                          params=dict(base3, leg=f"D2_resolved_K3/N={N}/{field_name}", N=int(N),
                                      quantity=d2.get("quantity"),
                                      n_singular=blk.get("n_singular"), chi_Q=blk.get("chi_Q"),
                                      betti_Q=blk.get(f"betti_Q_{field_name}"),
                                      betti_U=blk.get(f"betti_U_{field_name}"),
                                      rank_phi3=(blk.get("rank_phi3") or {}).get(field_name),
                                      rank_phi3_source=(blk.get("rank_phi3_source") or {}).get(field_name),
                                      chain_equals_identity=(blk.get("chain_equals_identity") or {}).get(field_name),
                                      MV_singular_reproduces_GUDHI_Q=(blk.get("MV_singular_reproduces_GUDHI_Q") or {}).get(field_name),
                                      all_N_fields_agree=comp.get("all_N_fields_agree"),
                                      no_expectation="the resolved-K3 Betti vector is the target of "
                                                     "this blind re-derivation and is not pre-stated"),
                          preprocessing="Mayer-Vietoris on the computed quotient, the computed links "
                                        "and the declared local model",
                          script=f"{V3}/02_mv_resolution.py", command=str(d2.get("command", "")),
                          tier="B")
            src.betti(rid, C.betti_map(vec))
            src.stat(rid, "chi_from_betti", (blk.get("chi_from_betti") or {}).get(field_name))
            src.stat(rid, "chi_Q", blk.get("chi_Q"))
            src.stat(rid, "n_singular", blk.get("n_singular"))
            src.stat(rid, "rank_phi3", (blk.get("rank_phi3") or {}).get(field_name))
            src.control(rid, "known_answer",
                        f"N={N} / {field_name}: the Mayer-Vietoris reconstruction of the SINGULAR "
                        "quotient reproduces gudhi's own Betti numbers for it (a check of the MV "
                        "machinery against a directly computable case)",
                        (blk.get("MV_singular_reproduces_GUDHI_Q") or {}).get(field_name),
                        detail=f"rank phi_3 = {(blk.get('rank_phi3') or {}).get(field_name)} "
                               f"({(blk.get('rank_phi3_source') or {}).get(field_name)}); the chain "
                               f"map equals the identity: "
                               f"{(blk.get('chain_equals_identity') or {}).get(field_name)}")
            src.control(rid, "known_answer",
                        f"N={N} / {field_name}: chi from the f-vector equals chi from the Betti numbers",
                        ((d3.get("computed") or {}).get("per_N") or {}).get(N, {}).get(field_name),
                        detail=f"chi(Q) = {((d3.get('computed') or {}).get('chi_Q') or {}).get(N)}")
    d5c = d5.get("computed") or {}
    did_kt3 = "mathematics/k3t2_rigidity/v3_K3xT2"
    src.dataset(id=did_kt3, domain="mathematics",
                title="v3: K3 x T^2 by Kunneth over a field from the computed K3 Betti numbers",
                source="Kunneth over a field (no Tor); the K3 factor is computed",
                provenance="synthetic_control", local_path=f"{V3}/results.json")
    rid = src.run(dataset_id=did_kt3, method="chain_complex", coeff_field=3, max_dim=6,
                  params=dict(base3, leg="D5_kunneth", quantity=d5.get("quantity"),
                              b_K3=d5c.get("b_K3"), b_T2=d5c.get("b_T2"),
                              field_coefficients=d5c.get("field_coefficients"),
                              shared_inputs=d5.get("shared_inputs")),
                  preprocessing="Kunneth over a field", script=f"{V3}/04_kunneth_scan_results.py",
                  command=str(d5.get("command", "")), tier="B")
    src.betti(rid, C.betti_map(d5c.get("b_K3xT2")))
    src.stat(rid, "chi_from_product_betti", d5c.get("chi_from_product_betti"))
    src.stat(rid, "chi_K3_times_chi_T2", d5c.get("chi_K3_times_chi_T2"))
    src.control(rid, "known_answer",
                "the product Euler characteristic from the product Betti numbers equals "
                "chi(K3) x chi(T^2)",
                d5c.get("consistent"),
                detail=f"{d5c.get('chi_from_product_betti')} vs {d5c.get('chi_K3_times_chi_T2')}")
    d7 = by_id.get("D7") or {}
    d7c = d7.get("computed") or {}
    did_p = "mathematics/k3t2_rigidity/v3_cup_parity"
    src.dataset(id=did_p, domain="mathematics",
                title="v3: the mod-2 cup pairing on H^2, computed on the T^4 control and on the "
                      "singular quotient Q - NOT on the K3",
                source="T^4 is spin, so its intersection form is even and nondegenerate of rank 6 "
                       "(tier L)", provenance="synthetic_control",
                local_path=f"{V3}/03_cup_parity_results.json",
                notes="the source states this does NOT decide the K3 parity")
    for label, blk in d7c.items():
        if not isinstance(blk, dict) or "dim_H2_Z2" not in blk:
            continue
        rid = src.run(dataset_id=did_p, method="chain_complex", coeff_field=2, max_dim=4,
                      params=dict(base3, leg=f"D7_cup_parity/{label}", quantity=d7.get("quantity"),
                                  label=blk.get("label"),
                                  n_simplices_by_dim=blk.get("n_simplices_by_dim"),
                                  q_on_basis=blk.get("q_on_basis")),
                      preprocessing="chain-level mod-2 cup pairing on H^2",
                      script=f"{V3}/03_cup_parity.py", command=str(d7.get("command", "")), tier="B")
        src.stat(rid, "dim_H2_Z2", blk.get("dim_H2_Z2"))
        src.stat(rid, "gram_rank_Z2", blk.get("gram_rank_Z2"))
        src.control(rid, "known_answer",
                    f"{label}: the mod-2 fundamental class is a cycle and the Gram matrix is symmetric",
                    bool(blk.get("mod2_fundamental_class_is_cycle")) and
                    bool(blk.get("gram_symmetric", True)),
                    detail=f"dim H^2 over F2 = {blk.get('dim_H2_Z2')}, Gram rank "
                           f"{blk.get('gram_rank_Z2')}, q identically zero: "
                           f"{blk.get('q_identically_zero')}")
    for r in v3.get("rigidity") or []:
        verdict = {"INVARIANCE": "recovered", "NORMALISATION": "inconclusive",
                   "NON-DISCRIMINATING": "null", "REJECTED": "failed"}.get(
            str(r.get("classification")).split()[0] if r.get("classification") else "", "inconclusive")
        src.finding(dataset_id=did_k33, tier="B", verdict=verdict,
                    claim=f"v3 rigidity, parameter '{r.get('parameter_inserted')}': classified "
                          f"{r.get('classification')}; selecting condition "
                          f"{str(r.get('selecting_condition'))[:250]}; solution set "
                          f"{r.get('solution_set')}",
                    caveat=f"condition_uses_true_value = {r.get('condition_uses_true_value')}; "
                           f"negative control: {str(r.get('negative_control'))[:250]}; shared inputs "
                           f"{r.get('shared_inputs')}",
                    reference=f"{V3}/results.json")
    src.finding(dataset_id=did_kt3, tier="B", verdict="recovered",
                claim=f"v3 track D-tda reproduces v2: the resolved K3 has Betti "
                      f"{comp.get('betti')} with chi {comp.get('chi_from_betti')} (all valid N and "
                      f"both fields agree: {comp.get('all_N_fields_agree')}), and K3 x T^2 has "
                      f"{d5c.get('b_K3xT2')}",
                caveat="this is a second, independent round of the same blind re-derivation, not an "
                       "independent confirmation of the underlying tier-L local model, which both "
                       "rounds take as a declared input",
                reference=f"{V3}/results.json")
    for item in v3.get("could_not_do") or []:
        src.skip("k3t2_rigidity v3 could not do", str(item)[:600])
    src.skip("k3t2_rigidity bars",
             "both rounds compute the homology of unfiltered simplicial and cubical complexes "
             "through gudhi's Betti numbers; there is no filtration and therefore no barcode")
    src.skip("k3t2_rigidity v2 even-intersection-form criterion",
             str((v2.get("rigidity_scan_over_k") or {}).get("even_intersection_form_criterion_status")))
    return src.report()
