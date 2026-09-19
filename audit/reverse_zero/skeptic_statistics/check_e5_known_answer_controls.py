#!/usr/bin/env python3
"""Skeptic check: re-execute E5's two known-answer controls (circle H1, planted voids H2) through the
committed cosmic_web_tda_scaled.py functions, twice in-process, and compare with the committed report
(e5_cosmic_web_tda_scaled_report.json known_answer_controls). Seeds as in the committed script
(SEED_SUBSAMPLE_REAL=42). Side files from alpha_persistence go to --scratch. Command (worktree root):
  P audit/reverse_zero/skeptic_statistics/check_e5_known_answer_controls.py --scratch <dir>
"""
import argparse, json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
E5 = os.path.join(HERE, "..", "E5-cosmic-web-tda-scaled")
sys.path.insert(0, E5)
import cosmic_web_tda_scaled as B  # noqa: E402
ap = argparse.ArgumentParser(); ap.add_argument("--scratch", required=True); a = ap.parse_args()
committed = json.load(open(os.path.join(E5, "e5_cosmic_web_tda_scaled_report.json")))["known_answer_controls"]
runs = []
for rep in range(2):
    cx = B.build_circle_cloud(n=2000, seed=B.SEED_SUBSAMPLE_REAL)
    _, bdc = B.alpha_persistence(cx, max_alpha_sq=100.0, label="ctrl_circle_%d" % rep, out_dir=a.scratch)
    h1 = B.top_bars(bdc, 1, k=3, r_trunc=10.0)
    ratio = float((h1[0, 1] - h1[0, 0]) / max(h1[1, 1] - h1[1, 0], 1e-9))
    vx, _ = B.build_planted_void_cloud(n_fill=15000, box=100.0, void_radius=15.0, n_voids=6, seed=B.SEED_SUBSAMPLE_REAL)
    vt = 2.5 * 15.0
    _, bdv = B.alpha_persistence(vx, max_alpha_sq=vt ** 2, label="ctrl_voids_%d" % rep, out_dir=a.scratch)
    h2 = B.top_bars(bdv, 2, k=8, r_trunc=vt)
    deaths = [float(x[1]) for x in h2[:6]]
    runs.append({"circle_ratio": ratio, "void_n_points": int(len(vx)), "top8_H2": [[float(x[0]), float(x[1])] for x in h2],
                 "n_top6_death_within_20pct": int(sum(abs(d - 15.0) <= 3.0 for d in deaths)),
                 "max_abs_death_err_pct_top6": float(max(abs(d - 15.0) / 15.0 * 100 for d in deaths))})
out = {"runs": runs, "committed_circle_ratio": committed["circle_H1"]["dominant_to_second_H1_ratio"],
       "committed_top8_H2": committed["planted_voids_H2"]["top_H2_bars_birth_death"],
       "committed_gate": committed["planted_voids_H2"]["gate"],
       "identical_to_committed": bool(runs[0]["top8_H2"] == committed["planted_voids_H2"]["top_H2_bars_birth_death"]
                                      and runs[0]["circle_ratio"] == committed["circle_H1"]["dominant_to_second_H1_ratio"]),
       "rep0_equals_rep1": runs[0] == runs[1],
       "density_note": "void control density 15000/100^3 (Mpc/h)^-3 in a clean periodic-free cube with perfectly empty spheres; not matched to the data's density or footprint"}
json.dump(out, open(os.path.join(HERE, "check_e5_known_answer_controls.json"), "w"), indent=1)
print(json.dumps({k: v for k, v in out.items() if k != "committed_top8_H2"}, indent=1))
