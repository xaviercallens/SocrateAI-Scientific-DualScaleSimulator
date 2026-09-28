"""Vendored copy of the FIXED TDA library.

SOURCE (not modified here, byte-identical copies):
  branch loop/tda-simple, worktree /mnt/disks/disk-socrateai-local-1/dualscale-wt-tdasimple
  audit/tda_validation/tda_fixed/{stats.py,cmb_topology.py}
  source commit d8175f1db4d395edf7151125a631c44009d500ca (range 1bc7ebc..d8175f1)
  sha256 stats.py        9ae09e5489a3d1d474c4d2be44d341a968b26d049e76a9b1e791021545af8a61
  sha256 cmb_topology.py 1f21825c75e560b98f7cc5546178aeb23db9cffb4b2ee29db44555f09dbc5a39

These REPLACE the defective originals in audit/reverse_zero/E5-cmb-tda/cmb_tda.py:
  stats.coarse_stats_fixed          <- cmb_tda.coarse_stats   (D1: miscalibrated chi2)
  cmb_topology.build_topology_fixed <- cmb_tda.build_topology (D2: b2 = 49147 on the full sky)
  cmb_topology.betti_curves_from_topology (true Euler characteristic)

cmb_tda.py is NOT imported anywhere under audit/cosmic_vorticity/.
The RANK p-value is preferred; the chi2 branch remains mildly anti-conservative
(0.058 vs 0.05) and is reported as a diagnostic only.

Tier of every number produced by this package: X (numerics).
"""
__all__ = ["stats", "cmb_topology"]
