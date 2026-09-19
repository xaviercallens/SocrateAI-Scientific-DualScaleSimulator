"""Guard-verification switch.

A regression guard that has never been seen to fail is not a guard
(CLAUDE.md rule 2: a check that cannot fail is not a check).  With

    TDA_GUARD_AGAINST_ORIGINAL=1 <venv>/bin/python -m pytest \
        audit/tda_validation/tda_fixed/tests -q

the three fixed entry points are rebound to the ORIGINAL implementations in
``audit/reverse_zero/E5-cmb-tda/cmb_tda.py`` (with their return values mapped
onto the fixed API, keeping the original SEMANTICS: df is always the nominal
bin count, no bin is ever dropped, and 'euler' means b0 - b1).  The tests that
assert the fixed behaviour must then FAIL.  The measured outcome of both modes
is recorded in ``validation_results.json`` under ``regression_guard``.

Without the variable the suite runs normally against the fixed library.
"""
import importlib.util
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
WT = os.path.abspath(os.path.join(PKG, "..", "..", ".."))
sys.path.insert(0, os.path.dirname(PKG))


def _original():
    spec = importlib.util.spec_from_file_location(
        "cmb_tda_original_for_conftest",
        os.path.join(WT, "audit/reverse_zero/E5-cmb-tda/cmb_tda.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def pytest_configure(config):
    if os.environ.get("TDA_GUARD_AGAINST_ORIGINAL") != "1":
        return
    from tda_fixed import cmb_topology as ct
    from tda_fixed import stats as fs
    orig = _original()

    def coarse_stats_shim(sim_curves, data_curve, n_bins=8, **kw):
        r = dict(orig.coarse_stats(sim_curves, data_curve, n_bins=n_bins))
        r.update({"df": r["n_bins"], "retained_rank": r["n_bins"],
                  "n_bins_kept": r["n_bins"], "n_bins_nominal": n_bins,
                  "kept_bin_indices": list(range(n_bins)),
                  "kept_bin_indices_pooled_MISSING": True,
                  "dropped_bin_indices": [], "dropped_bin_reasons": {},
                  "dropped_bin_ensemble_values": {}, "test_differs_in_dropped_bin": {},
                  "kept_block_condition_number": float("nan")})
        return r

    def betti_curves_shim(temp, unmasked, edges_arr, tris_arr, nu_grid, sublevel=True):
        b0, b1, chi = orig.betti_curves_from_topology(temp, unmasked, edges_arr, tris_arr,
                                                      nu_grid, sublevel=sublevel)
        n = np.asarray(nu_grid).size
        return {"b0": b0, "b1": b1,
                "euler_char_true": chi,          # the ORIGINAL mislabelling: b0 - b1
                "b0_minus_b1": b0 - b1, "b2_implied": np.zeros(n, dtype=int),
                "n_vertices_at_nu": np.zeros(n, dtype=int),
                "n_edges_at_nu": np.zeros(n, dtype=int),
                "n_triangles_at_nu": np.zeros(n, dtype=int)}

    def build_topology_shim(mask, nside, return_info=False):
        unm, e, tr = orig.build_topology(mask, nside)
        if return_info:
            V, E, F = int(unm.size), int(e.shape[0]), int(tr.shape[0])
            return unm, e, tr, {"n_vertices": V, "n_edges": E, "n_triangles": F,
                                "n_quad_corners": -1, "n_tri_corners": -1,
                                "euler_char_V_minus_E_plus_F": V - E + F}
        return unm, e, tr

    fs.coarse_stats_fixed = coarse_stats_shim
    ct.build_topology_fixed = build_topology_shim
    ct.betti_curves_from_topology = betti_curves_shim
    print("\n[conftest] TDA_GUARD_AGAINST_ORIGINAL=1: fixed entry points rebound to the "
          "ORIGINAL cmb_tda implementations; the fixed-behaviour tests must now FAIL.")
