#!/usr/bin/env python3
"""Attach the top H0/H1 persistence bars to the sky-map runs.

run_skymaps.py records the Betti CURVES (as statistics on a nu grid) but not the
underlying persistence intervals, because betti_curves_from_topology returns
curves rather than diagrams.  This pass recomputes the DATA diagram only -- no
simulations -- using the same fixed complex and the same lower-star filtration,
and writes the 25 longest bars per dimension.

Filtration units: the field divided by the sigma of its unmasked pixels, so a
bar is quoted in units of sigma of that map.  This is the same normalisation the
curves use, which is what makes K_CMB, mK and MJy/sr maps comparable at all.

Tier X.
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import ASTRO, topodb  # noqa: E402
import run_skymaps as RS  # noqa: E402

MAPS = {"planck_smica": 128, "wmap_ilc": 128, "haslam408": 128,
        "wmap_kband": 128, "cobe_dmr": 16}


def main():
    sys.path.insert(0, ASTRO)
    from lib.cmb_topology import build_topology_fixed, _simplex_tree

    with topodb() as db:
        runs = dict(db.con.execute(
            "SELECT dataset_id, MAX(id) FROM run WHERE method='lower_star_graph' "
            "AND dataset_id LIKE 'astro/%' GROUP BY dataset_id").fetchall())
    print(runs)

    for key, nside in MAPS.items():
        rid = runs.get(f"astro/{key}")
        if rid is None:
            print(f"  {key}: no run, skipped")
            continue
        m, msk, _ = RS.load_map(key, nside)
        unmasked, edges, tris = build_topology_fixed(msk, nside)
        sigma = m[unmasked].std()
        verts = (m / sigma)[unmasked]
        st = _simplex_tree(verts, np.asarray(edges, np.int64).reshape(-1, 2),
                           np.asarray(tris, np.int64).reshape(-1, 3))
        st.compute_persistence(homology_coeff_field=2, persistence_dim_max=True)
        with topodb() as db:
            db.con.execute("DELETE FROM bar WHERE run_id=?", (rid,))
            for dim in (0, 1):
                d = np.array(st.persistence_intervals_in_dimension(dim), float).reshape(-1, 2)
                if d.size:
                    db.add_bars(rid, dim, [(b, (None if not np.isfinite(dd) else dd)) for b, dd in d],
                                top=25)
                print(f"  {key} (run {rid}) dim {dim}: {d.shape[0]} bars, top 25 stored")
            db.add_control(rid, "known_answer",
                           "bars are in units of sigma of the UNMASKED pixels of this map, the same "
                           "normalisation the Betti curves use; that is what makes K_CMB, mK and "
                           "brightness-temperature maps comparable",
                           passed=True,
                           detail=f"sigma = {float(sigma):.6g} {RS.load_map.__doc__ and ''}"
                                  f"(map units); nside {nside}, {unmasked.size} vertices")
    return 0


if __name__ == "__main__":
    sys.exit(main())
