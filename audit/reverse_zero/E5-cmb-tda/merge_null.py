#!/usr/bin/env python
"""
Pool several --dump-npz null chunks (from cmb_tda.py run_null) produced with
disjoint --seed-offset values into ONE null ensemble, verify the chunks are
consistent (same data curves, same topology counts, disjoint seeds), and
write:
  - a merged .npz (same schema as a single-chunk dump, for --null-cache)
  - a JSON report with the full coarse_stats (chi2, p-value, LOO empirical
    rank p, ensemble mean/std curves, standardized residuals) for b0/b1/chi
    in BOTH sublevel and superlevel directions.

Usage:
  .venv-tda/bin/python merge_null.py chunk0.npz chunk1.npz ... \
      --out-npz null_merged.npz --out-json null_report.json
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cmb_tda as C


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("chunks", nargs="+")
    ap.add_argument("--out-npz", required=True)
    ap.add_argument("--out-json", required=True)
    args = ap.parse_args()

    loaded = [np.load(c) for c in args.chunks]

    # consistency checks
    ref = loaded[0]
    for i, c in enumerate(loaded[1:], start=1):
        assert np.allclose(c["b0_d"], ref["b0_d"]) and np.allclose(c["chi_d"], ref["chi_d"]), \
            f"chunk {args.chunks[i]} has different DATA curves than chunk 0 -- not comparable"
        assert np.allclose(c["nu_grid"], ref["nu_grid"]), f"chunk {args.chunks[i]} has different nu_grid"
        assert int(c["n_unmasked"]) == int(ref["n_unmasked"]), f"chunk {args.chunks[i]} has different topology"
    all_seeds = np.concatenate([c["seeds"] for c in loaded])
    assert len(set(all_seeds.tolist())) == len(all_seeds), f"seed collision across chunks: {sorted(all_seeds.tolist())}"
    assert len(all_seeds) >= 100, (
        f"ground rules require >=100 null sims; only {len(all_seeds)} pooled from "
        f"{len(args.chunks)} chunk(s) -- add more chunks before merging"
    )

    def cat(key):
        return np.concatenate([c[key] for c in loaded], axis=0)

    sim_b0 = cat("sim_b0"); sim_b1 = cat("sim_b1"); sim_chi = cat("sim_chi")
    sim_b0_s = cat("sim_b0_super"); sim_b1_s = cat("sim_b1_super"); sim_chi_s = cat("sim_chi_super")
    sim_var = cat("sim_var")
    n_total = sim_chi.shape[0]

    b0_d = ref["b0_d"]; b1_d = ref["b1_d"]; chi_d = ref["chi_d"]
    b0_d_s = ref["b0_d_super"]; b1_d_s = ref["b1_d_super"]; chi_d_s = ref["chi_d_super"]
    nu_grid = ref["nu_grid"]
    data_var = float(ref["data_var"])
    fsky = float(ref["fsky"])

    sim_var_mean = float(sim_var.mean())
    var_ratio = sim_var_mean / data_var if data_var else float("nan")

    np.savez(args.out_npz,
             sim_b0=sim_b0, sim_b1=sim_b1, sim_chi=sim_chi,
             sim_b0_super=sim_b0_s, sim_b1_super=sim_b1_s, sim_chi_super=sim_chi_s,
             b0_d=b0_d, b1_d=b1_d, chi_d=chi_d,
             b0_d_super=b0_d_s, b1_d_super=b1_d_s, chi_d_super=chi_d_s,
             nu_grid=nu_grid, seeds=all_seeds, data_var=data_var,
             sim_var=sim_var, fsky=fsky,
             n_unmasked=int(ref["n_unmasked"]), n_edges=int(ref["n_edges"]), n_tri=int(ref["n_tri"]),
             lmax=int(ref["lmax"]))

    report = {
        "tier": "X (exploratory numerics)",
        "framing": "M0 predicts Gaussian statistics (no topological defects). A null "
                    "result here is consistent with M0 and plain LCDM equally; it is "
                    "NOT evidence for K3xT2 and changes the free-parameter count by 0.",
        "n_sims_total": n_total,
        "chunks_merged": args.chunks,
        "seeds_all": sorted(int(x) for x in all_seeds.tolist()),
        "nu_grid_sigma": nu_grid.tolist(),
        "fsky": fsky, "n_unmasked_pixels_work_res": int(ref["n_unmasked"]),
        "n_edges": int(ref["n_edges"]), "n_triangles": int(ref["n_tri"]),
        "lmax_used": int(ref["lmax"]),
        "calibration_check": {
            "data_unmasked_pixel_variance": data_var,
            "sim_ensemble_mean_variance": sim_var_mean,
            "ratio_sim_over_data": var_ratio,
            "note": "ratio should be close to 1; computed over the pooled N_sims_total ensemble",
        },
        "data_betti_curve": {
            "sublevel": {"b0": b0_d.tolist(), "b1": b1_d.tolist(), "euler_chi": chi_d.tolist()},
            "superlevel": {"b0": b0_d_s.tolist(), "b1": b1_d_s.tolist(), "euler_chi": chi_d_s.tolist()},
        },
        "statistic_sublevel": {
            "b0": C.coarse_stats(sim_b0, b0_d),
            "b1": C.coarse_stats(sim_b1, b1_d),
            "euler_chi": C.coarse_stats(sim_chi, chi_d),
        },
        "statistic_superlevel": {
            "b0": C.coarse_stats(sim_b0_s, b0_d_s),
            "b1": C.coarse_stats(sim_b1_s, b1_d_s),
            "euler_chi": C.coarse_stats(sim_chi_s, chi_d_s),
        },
        "n_coarse_bins": C.N_BINS,
    }
    with open(args.out_json, "w") as f:
        json.dump(report, f, indent=2)
    print(f"merged {len(args.chunks)} chunks, n_sims_total={n_total}, wrote {args.out_npz} and {args.out_json}", file=sys.stderr)


if __name__ == "__main__":
    main()
