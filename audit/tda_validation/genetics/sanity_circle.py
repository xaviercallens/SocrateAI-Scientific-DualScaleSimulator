"""Known-answer check of the pipeline code on its own synthetic helper.
Command: python sanity_circle.py  (seed fixed inside build_circle_cloud: seed=0)"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from tda_common import PIPE, alpha_h1_summary, rips_h1_summary, pipeline_sha256, dump
xyz = PIPE.build_circle_cloud(200, noise_sigma=0.05, seed=0)
a = alpha_h1_summary(xyz, "circle")
r = rips_h1_summary(points=xyz)
rng = np.random.RandomState(1)
blob = rng.normal(size=(200, 3))
b = alpha_h1_summary(blob, "gauss_blob")
out = {"pipeline_file_sha256": pipeline_sha256(),
       "circle_n200_sigma0.05_seed0": {"alpha": a, "rips": {k: v for k, v in r.items() if k != "all_h1_persistence_sorted"}},
       "gaussian_blob_n200_seed1_alpha": b,
       "criterion": "one H1 bar with death ~ radius 1 and P1/P2 >= 3",
       "passed": bool(a["P1_over_P2"] >= 3 and 0.7 <= a["top_h1_bars_birth_death"][0][1] <= 1.1)}
dump(out, "sanity_circle.json")
print(out["passed"], a["top_h1_bars_birth_death"][:3], a["P1_over_P2"], a["S"], "| blob", b["P1_over_P2"], b["S"], "| rips", r["top_h1_bars_birth_death"][:2])
