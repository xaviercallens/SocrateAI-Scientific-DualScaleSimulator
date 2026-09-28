"""Labelled synthetic control: /mnt/disks/.../bio_datasets/01_hic_*.npy.

Those local arrays carry real-source names but are SYNTHETIC.  The earlier genetics
validation established this; the diagnostics are RECOMPUTED here rather than cited,
so the TopoDB record rests on numbers this script produced:

  * all 23 files in that directory were written within a fraction of a second;
  * 01_hic_contact_map is an exactly monotone function of the pairwise distances of
    01_hic_chromatin_coords (Spearman -1.0), which no measured Hi-C matrix is;
  * it has no zero off-diagonal entries, and a constant diagonal.

It is ingested with provenance 'synthetic_control' so that the database records why
it must never be used as biological data, and so a future query for "Hi-C" cannot
pick it up as an observation.  Its persistence is computed too, which shows what a
noiseless synthetic ring looks like next to the three real bacterial Hi-C runs.

Usage: python run_synthetic_control.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.stats import spearmanr

sys.path.insert(0, "/mnt/disks/disk-socrateai-local-1/wt-topo-bio/topodb_runs/biology")
from bio_common import (RESULTS, Block, S_stat, dominance, finite, max_pers,  # noqa: E402
                        rips_collapsed_from_distance, sha256, top_bars)

SEED = 20260920
SRC = Path("/mnt/disks/disk-socrateai-local-1/bio_datasets")
SCRIPT = "topodb_runs/biology/run_synthetic_control.py"
COMMAND = ("timeout 300 prlimit --as=8589934592 -- "
           "/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python "
           "topodb_runs/biology/run_synthetic_control.py")


def main() -> int:
    C = np.load(SRC / "01_hic_contact_map.npy")
    X = np.load(SRC / "01_hic_chromatin_coords.npy")
    n = C.shape[0]

    # --- recomputed evidence that this is synthetic ---
    D = squareform(pdist(X))
    iu = np.triu_indices(n, 1)
    rho = float(spearmanr(C[iu], D[iu]).statistic)
    mtimes = sorted(p.stat().st_mtime for p in SRC.glob("*"))
    ev = {"n_bins": int(n),
          "spearman_contact_vs_coordinate_distance": round(rho, 6),
          "frac_offdiag_exactly_zero": float(np.mean(C[iu] == 0)),
          "diagonal_unique_values": sorted(set(np.round(np.diag(C), 9).tolist()))[:5],
          "n_files_in_directory": len(mtimes),
          "mtime_span_seconds": round(mtimes[-1] - mtimes[0], 4),
          "coords_shape": list(X.shape)}
    print(json.dumps(ev, indent=1))

    # --- persistence of the same array, for the record ---
    Cs = C.copy()
    np.fill_diagonal(Cs, 0.0)
    off = ~np.eye(n, dtype=bool)
    pos = Cs[off & (Cs > 0)]
    Cs[off & (Cs <= 0)] = 0.5 * float(pos.min()) if pos.size else 1.0
    Dm = np.zeros_like(Cs)
    Dm[off] = Cs[off] ** (-1.0 / 3.0)
    t = time.time()
    dg = rips_collapsed_from_distance(Dm, max_hom_dim=1)
    wall = round(time.time() - t, 2)
    dom = dominance(dg[1])

    blk = Block("synthetic_control", SCRIPT, COMMAND)
    blk.dataset(id="synthetic/bio_datasets_01_hic_contact_map", domain="synthetic",
                title="LABELLED SYNTHETIC: bio_datasets/01_hic_contact_map.npy "
                      "(named as Hi-C, generated from coordinates)",
                source=f"local file {SRC / '01_hic_contact_map.npy'}; no upstream URL exists",
                provenance="synthetic_control", n_objects=int(n), ambient_dim=int(n),
                units="arbitrary 'contact' units; distance = contact^(-1/3)",
                sha256=sha256(SRC / "01_hic_contact_map.npy"),
                local_path=str(SRC / "01_hic_contact_map.npy"),
                notes="NOT DATA. Recomputed evidence: " + json.dumps(ev) + ". A measured Hi-C "
                      "matrix is never an exactly monotone function of a coordinate set's "
                      "pairwise distances, never has zero off-diagonal zeros, and does not have a "
                      "constant diagonal. Ingested only so that a query for Hi-C in this database "
                      "meets it labelled as a synthetic control rather than as an observation.")
    blk.run(dataset_id="synthetic/bio_datasets_01_hic_contact_map", method="rips", coeff_field=2,
            max_dim=1,
            params={"metric": "contact^(-1/3) fed directly to Rips (NO embedding)",
                    "exponent": 1 / 3, "max_hom_dim": 1, "edge_collapse": True,
                    "expansion_dim": 2, "n_bins": int(n)},
            preprocessing="diagonal zeroed; no other cleaning was needed (there are no "
                          "off-diagonal zeros, which is itself the tell)",
            seed=SEED, tier="X", wall_sec=wall, diagrams=dg,
            betti={1: 1 if dom >= 2 else 0},
            stats=[{"name": "h1_dominance_P1_over_P2", "value": dom},
                   {"name": "h1_S_longest_over_total", "value": S_stat(dg[1])},
                   {"name": "h1_max_persistence", "value": max_pers(dg[1])},
                   {"name": "h1_finite_bar_count", "value": float(len(finite(dg[1])))},
                   {"name": "spearman_contact_vs_coordinate_distance", "value": rho},
                   {"name": "mtime_span_of_all_23_files_seconds",
                    "value": ev["mtime_span_seconds"]}],
            controls=[{"kind": "negative",
                       "description": "provenance check: a real Hi-C contact matrix cannot be an "
                                      "exactly monotone function of a coordinate set's distances",
                       "passed": False,
                       "detail": f"Spearman(contact, coordinate distance) = {rho:.6f}; "
                                 f"{ev['frac_offdiag_exactly_zero']:.3f} of off-diagonal entries "
                                 f"are zero; the whole directory was written within "
                                 f"{ev['mtime_span_seconds']} s"}],
            findings=[{"claim": f"bio_datasets/01_hic_contact_map.npy is synthetic, not Hi-C: it is "
                                f"an exactly monotone function of 01_hic_chromatin_coords "
                                f"(Spearman {rho:.4f}), it has no off-diagonal zeros, and all "
                                f"{ev['n_files_in_directory']} files in that directory were written "
                                f"within {ev['mtime_span_seconds']} s. Its Rips H1 dominance is "
                                f"{dom:.2f}. It is recorded as a synthetic control and is used as "
                                f"data nowhere in this block.",
                       "verdict": "artefact", "tier": "X",
                       "caveat": "the mtime span is filesystem metadata, which is suggestive rather "
                                 "than conclusive on its own; the Spearman -1 relation is the "
                                 "decisive evidence",
                       "reference": "/mnt/disks/disk-socrateai-local-1/dualscale-wt-tdaval/audit/"
                                    "tda_validation/genetics/local_datasets_assessment.json"}])
    (RESULTS / "synthetic_control.json").write_text(json.dumps(
        {"evidence": ev, "h1_dominance": dom, "seed": SEED}, indent=1))
    print(f"dominance {dom:.3f}, top bars {top_bars(dg[1], 3)}")
    blk.write()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
