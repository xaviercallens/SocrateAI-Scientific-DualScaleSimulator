"""tda_fixed -- corrected implementations of the TDA code that the
known-answer suite on branch loop/tda-simple found to be defective.

The originals in ``audit/reverse_zero/E5-cmb-tda/cmb_tda.py`` and
``audit/reverse_zero/E5-cosmic-web-tda-scaled/`` are the committed record of
past runs and are NOT changed here: other branches' results were produced with
them.  Each affected function in the originals carries a marked comment block
naming the defect, the date and the replacement below.

  stats.coarse_stats_fixed          replaces cmb_tda.coarse_stats
                                    (defect D1: miscalibrated chi2 p-values)
  cmb_topology.build_topology_fixed replaces cmb_tda.build_topology
                                    (defect D2: every 4-clique filled, so the
                                     full sky had b2 = 49147 instead of 1)
  cmb_topology.betti_curves_from_topology
                                    replaces cmb_tda.betti_curves_from_topology
                                    (same filtration; returns the TRUE Euler
                                     characteristic under 'euler_char_true'
                                     and the old b0-b1 under 'b0_minus_b1')
  pointcloud.alpha_persistence      NO defect found; carried here only to give
                                    run_validation.py one import surface

The fixes and every acceptance gate were declared in ``fix_expectations.json``
(commit 1bc7ebc) before this package was written.  Measured before/after
numbers are in ``validation_results.json``, produced by ``run_validation.py``.
Unit tests that FAIL on the originals and PASS here are in ``tests/``.

Tier of every number produced by this package: X (numerics).
"""

__all__ = ["stats", "cmb_topology", "pointcloud"]
