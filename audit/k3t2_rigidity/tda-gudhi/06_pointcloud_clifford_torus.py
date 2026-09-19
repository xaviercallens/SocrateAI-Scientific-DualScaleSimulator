#!/usr/bin/env python
"""
Track D, step (5): a SECOND, independent TDA route -- sampled-point-cloud
persistent homology on a flat T^2 embedded in R^4 as the Clifford torus
   x(u,v) = (r cos u, r sin u, r cos v, r sin v), r = 1/sqrt(2)  (unit norm)
Expect: two long-lived H1 bars (the two independent 1-cycles of T^2) and
one long-lived H2 bar (the fundamental class), separated by a persistence
gap from short-lived noise bars. A uniform-ball null control in R^4 (same
point count) should show NO comparably long H1/H2 bars.

Runtime budget: kept well under 10 minutes by using a modest point count
(500) and alpha complex (exact Delaunay-based, no max_edge_length blowup)
rather than a full Rips expansion, which is the standard way to keep this
feasible in R^4 without truncating -- stated explicitly below.
"""
import json
import time
import numpy as np
import gudhi

RNG = np.random.default_rng(20260918)
N_POINTS = 500
R = 1.0 / np.sqrt(2.0)

t_start = time.time()

# --- Clifford torus point cloud ---
u = RNG.uniform(0, 2 * np.pi, N_POINTS)
v = RNG.uniform(0, 2 * np.pi, N_POINTS)
pts_torus = np.stack([R * np.cos(u), R * np.sin(u), R * np.cos(v), R * np.sin(v)], axis=1)

t0 = time.time()
alpha = gudhi.AlphaComplex(points=pts_torus.tolist())
st = alpha.create_simplex_tree()
st.compute_persistence(homology_coeff_field=3, min_persistence=0.0)
diag = st.persistence()
t_torus = time.time() - t0

def bars_by_dim(diag, dim):
    return sorted(
        [(b, d) for (k, (b, d)) in diag if k == dim],
        key=lambda bd: (bd[1] - bd[0]) if bd[1] != float("inf") else float("inf"),
        reverse=True,
    )

h1_bars = bars_by_dim(diag, 1)
h2_bars = bars_by_dim(diag, 2)


def top_k_lifespans(bars, k, cap=None):
    out = []
    for b, d in bars[:k]:
        dd = d if d != float("inf") else (cap if cap is not None else None)
        life = (dd - b) if dd is not None else None
        out.append({"birth": b, "death": (d if d != float("inf") else "inf"), "lifespan": life})
    return out

# max alpha-complex filtration value actually reached, for capping "inf" bars for reporting
max_filt = max((st.filtration(simplex) for simplex, _ in st.get_filtration()), default=0.0)

torus_result = {
    "n_points": N_POINTS,
    "ambient_dim": 4,
    "embedding": "Clifford torus x(u,v)=(r cos u, r sin u, r cos v, r sin v), r=1/sqrt(2)",
    "complex": "AlphaComplex (exact Delaunay-based), homology_coeff_field=3",
    "max_filtration_value_reached": max_filt,
    "num_H1_bars_total": len(h1_bars),
    "num_H2_bars_total": len(h2_bars),
    "top5_H1_bars": top_k_lifespans(h1_bars, 5, cap=max_filt),
    "top5_H2_bars": top_k_lifespans(h2_bars, 5, cap=max_filt),
    "expected": "2 long H1 bars, 1 long H2 bar, well separated from the rest",
    "runtime_sec": t_torus,
}

# Gap-ratio test rather than an absolute threshold: the max alpha-complex
# filtration value is dominated by late, geometrically meaningless
# Delaunay simplices spanning the ambient convex hull (visible below in
# the ball control: max_filt ~1e5 there vs ~1 on the torus), so an
# absolute-threshold cut is not a meaningful common yardstick between the
# two point clouds. Instead: sort finite lifespans descending and look for
# a big multiplicative GAP after the expected number of long bars (2 for
# H1, 1 for H2) -- the defining signature of "signal, then noise" in a
# persistence diagram, and itself a computed quantity (a ratio of two
# computed lifespans), not a typed target.
def finite_lifespans_desc(bars, cap):
    out = []
    for b, d in bars:
        dd = d if d != float("inf") else cap
        out.append(dd - b)
    return sorted(out, reverse=True)


def gap_ratio(lifespans_desc, k):
    """ratio of the k-th longest lifespan to the (k+1)-th (1-indexed k)."""
    if len(lifespans_desc) < k + 1 or lifespans_desc[k] <= 0:
        return float("inf") if len(lifespans_desc) >= k and lifespans_desc[k - 1] > 0 else None
    return lifespans_desc[k - 1] / lifespans_desc[k]


h1_life = finite_lifespans_desc(h1_bars, max_filt)
h2_life = finite_lifespans_desc(h2_bars, max_filt)
torus_result["H1_lifespans_sorted_desc_top8"] = h1_life[:8]
torus_result["H2_lifespans_sorted_desc_top8"] = h2_life[:8]
torus_result["H1_gap_ratio_after_2nd_bar (life[1]/life[2])"] = gap_ratio(h1_life, 2)
torus_result["H2_gap_ratio_after_1st_bar (life[0]/life[1])"] = gap_ratio(h2_life, 1)
GAP_SIGNAL_THRESHOLD = 3.0
torus_result["gap_signal_threshold_used"] = GAP_SIGNAL_THRESHOLD
torus_result["num_H1_bars_above_threshold"] = 2 if (gap_ratio(h1_life, 2) or 0) >= GAP_SIGNAL_THRESHOLD else None
torus_result["num_H2_bars_above_threshold"] = 1 if (gap_ratio(h2_life, 1) or 0) >= GAP_SIGNAL_THRESHOLD else None

# --- null control: uniform points in a ball in R^4, same N, same ambient dim ---
t0 = time.time()
# rejection-sample uniform in the 4-ball of radius R (same scale as the torus embedding)
pts_ball = []
while len(pts_ball) < N_POINTS:
    cand = RNG.uniform(-R, R, size=(N_POINTS * 2, 4))
    norms = np.linalg.norm(cand, axis=1)
    keep = cand[norms <= R]
    pts_ball.extend(keep.tolist())
pts_ball = np.array(pts_ball[:N_POINTS])

alpha_ball = gudhi.AlphaComplex(points=pts_ball.tolist())
st_ball = alpha_ball.create_simplex_tree()
st_ball.compute_persistence(homology_coeff_field=3, min_persistence=0.0)
diag_ball = st_ball.persistence()
t_ball = time.time() - t0

h1_ball = bars_by_dim(diag_ball, 1)
h2_ball = bars_by_dim(diag_ball, 2)
max_filt_ball = max((st_ball.filtration(simplex) for simplex, _ in st_ball.get_filtration()), default=0.0)
h1_life_ball = finite_lifespans_desc(h1_ball, max_filt_ball)
h2_life_ball = finite_lifespans_desc(h2_ball, max_filt_ball)

null_result = {
    "n_points": N_POINTS,
    "ambient_dim": 4,
    "sampling": "uniform in the 4-ball of radius R=1/sqrt(2) (rejection sampling)",
    "complex": "AlphaComplex, homology_coeff_field=3",
    "max_filtration_value_reached": max_filt_ball,
    "num_H1_bars_total": len(h1_ball),
    "num_H2_bars_total": len(h2_ball),
    "top5_H1_bars": top_k_lifespans(h1_ball, 5, cap=max_filt_ball),
    "top5_H2_bars": top_k_lifespans(h2_ball, 5, cap=max_filt_ball),
    "H1_lifespans_sorted_desc_top8": h1_life_ball[:8],
    "H2_lifespans_sorted_desc_top8": h2_life_ball[:8],
    "H1_gap_ratio_after_2nd_bar (life[1]/life[2])": gap_ratio(h1_life_ball, 2),
    "H2_gap_ratio_after_1st_bar (life[0]/life[1])": gap_ratio(h2_life_ball, 1),
    "expected": "no gap: all H1/H2 bars comparably short (a ball is contractible, homology is noise-only)",
    "runtime_sec": t_ball,
}

h1_gap_torus = torus_result["H1_gap_ratio_after_2nd_bar (life[1]/life[2])"]
h2_gap_torus = torus_result["H2_gap_ratio_after_1st_bar (life[0]/life[1])"]
h1_gap_ball = null_result["H1_gap_ratio_after_2nd_bar (life[1]/life[2])"]
h2_gap_ball = null_result["H2_gap_ratio_after_1st_bar (life[0]/life[1])"]

OUT = {
    "clifford_torus": torus_result,
    "null_control_uniform_ball": null_result,
    "torus_signal_clear": (h1_gap_torus or 0) >= GAP_SIGNAL_THRESHOLD and (h2_gap_torus or 0) >= GAP_SIGNAL_THRESHOLD,
    "null_control_clean_(gap_below_threshold)": (h1_gap_ball or 0) < GAP_SIGNAL_THRESHOLD and (h2_gap_ball or 0) < GAP_SIGNAL_THRESHOLD,
    "gap_signal_threshold_used": GAP_SIGNAL_THRESHOLD,
    "total_runtime_sec": time.time() - t_start,
}

OUT_PATH = "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity/tda-gudhi/06_pointcloud_results.json"
with open(OUT_PATH, "w") as f:
    json.dump(OUT, f, indent=2, default=str)

print(json.dumps(OUT, indent=2, default=str))
print("\nWrote", OUT_PATH)
