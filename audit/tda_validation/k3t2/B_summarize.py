#!/usr/bin/env python
"""
Part B summary (tier X): reads B_<study>.json, B_explore.json, B2_injectivity.json; writes
B_summary.json and figure B3_K3_H2.png (H2 persistence diagram and beta_2(eps) curve, K3 vs density-matched null).
Run: prlimit --as=8589934592 -- <venv-python> B_summarize.py
"""
import json, math, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
L = lambda f: json.load(open(os.path.join(HERE, f)))
studies = ["T2_alpha", "T2_rips", "T3_rips", "T3_rips_fps", "T4_rips", "T4_rips_fps",
           "orbifold_rips", "orbifold_rips_fps", "K3_rips", "K3_rips_fps", "null4m_rips"]
S = {}
for s in studies:
    d = L(f"B_{s}.json")
    ok_runs = [r for r in d["runs"] if not r.get("failed")]
    last = max(ok_runs, key=lambda r: r["spec"]["N"])
    S[s] = {"space": d["config"]["space"], "method": d["config"]["method"], "tau": d["config"]["tau"],
            "fps": bool(d["config"].get("fps")), "expected": d["config"]["expected"],
            "N_min_all3seeds_F3": d["N_min"], "largest_N_completed": d["largest_N_completed"],
            "stopped_because": d["stopped_because"],
            "at_largest_N": {"n_simplices": last["info"]["n_simplices"], "seconds": last["wall"], "maxrss_MB": last["maxrss_MB"],
                             "best_window_ratio_F3": last["by_field"]["F3"]["best_window_ratio"],
                             "recovered_F2": last["by_field"]["F2"]["recovered"]},
            "per_N_seed": [{"N": r["spec"]["N"], "seed": r["spec"]["seed"],
                            **({"failed": r["reason"]} if r.get("failed") else
                               {"F3_ratio": round(r["by_field"]["F3"]["best_window_ratio"], 3),
                                "F3_recovered": r["by_field"]["F3"]["recovered"],
                                "F2_recovered": r["by_field"]["F2"]["recovered"],
                                "n_simplices": r["info"]["n_simplices"], "sec": r["wall"]})} for r in d["runs"]]}

# scaling fit on the same method (Rips, tau=1.2, iid uniform): N_min = A * B^d
fit = {}
pts = [(2, S["T2_rips"]["N_min_all3seeds_F3"]), (3, S["T3_rips"]["N_min_all3seeds_F3"])]
if all(n for _, n in pts):
    (d1, n1), (d2, n2) = pts
    B = (n2 / n1) ** (1 / (d2 - d1)); A = n1 / B ** d1
    fit = {"method": "Rips tau=1.2, iid uniform, grid of doublings (so N_min resolution is a factor 2)",
           "points": pts, "A": A, "B_per_dimension": B, "extrapolated_N_min_T4": A * B ** 4,
           "T4_status": f"not recovered up to N={S['T4_rips']['largest_N_completed']}; next grid N failed: {S['T4_rips']['stopped_because']}",
           "caveat": "two points only; extrapolation is indicative (tier X)"}
    pf = [(3, S["T3_rips_fps"]["N_min_all3seeds_F3"])]
    fit["fps_T3_N_min"] = pf[0][1]

# B_explore.json (8.9 MB, uncommitted) -> B_explore_trimmed.json (committed): identical except that the full H2
# diagrams of the K3xT2 runs (tens of thousands of noise bars, unused below) are replaced by their bar counts.
if os.path.exists(os.path.join(HERE, "B_explore.json")):
    full = L("B_explore.json")
    for k, r in full["runs"].items():
        if r["spec"]["space"] == "K3xT2" and not r.get("failed"):
            for f, v in r["by_field"].items():
                if "diagram_H2" in v:
                    v["diagram_H2_trimmed_n_bars"] = len(v.pop("diagram_H2"))
    full["trimmed_note"] = "K3xT2 diagram_H2 lists removed (counts kept); everything else identical to B_explore.json"
    json.dump(full, open(os.path.join(HERE, "B_explore_trimmed.json"), "w"), indent=1)
E = L("B_explore_trimmed.json")["runs"]


def curve(diag, tau, grid):
    b = np.array([x[0] for x in diag]); d = np.array([x[1] if x[1] is not None else np.inf for x in diag])
    return np.array([int(((b <= e) & (d > e)).sum()) for e in grid])


def windows_ge1(diag, tau):
    """longest window (ratio) on which beta >= 1 (union of bar intervals)."""
    iv = sorted((b, min(d if d is not None else tau, tau)) for b, d in diag if b < tau)
    best, cur = 0.0, None
    for b, d in iv:
        if cur and b <= cur[1]:
            cur[1] = max(cur[1], d)
        else:
            if cur and cur[0] > 0:
                best = max(best, cur[1] / cur[0])
            cur = [b, d]
    if cur and cur[0] > 0:
        best = max(best, cur[1] / cur[0])
    return best


def plateau(c, grid):
    """longest run of constant nonzero beta on the grid: (value, e_start, e_end, ratio)."""
    best = (0, None, None, 0.0)
    i = 0
    while i < len(c):
        j = i
        while j + 1 < len(c) and c[j + 1] == c[i]:
            j += 1
        if c[i] > 0 and grid[i] > 0:
            r = grid[j] / grid[i]
            if r > best[3]:
                best = (int(c[i]), float(grid[i]), float(grid[j]), float(r))
        i = j + 1
    return best


tau = 0.8
grid = np.linspace(0.3, 0.799, 500)
k3 = {}
for key in ("0", "1", "2", "3"):
    r = E[key]
    diag = r["by_field"]["F3"]["diagram_H2"]
    c = curve(diag, tau, grid)
    k3[f"seed{r['spec']['seed']}{'_fps' if r['spec'].get('fps') else ''}"] = {
        "N": r["spec"]["N"], "tau": tau, "n_H2_bars": len(diag), "n_H2_bars_alive_at_tau": sum(1 for x in diag if x[1] is None),
        "beta2_plateau(value,e1,e2,ratio)": plateau(c, grid),
        "beta2_ever_equals_22_on_grid": bool((c == 22).any()),
        "beta2_at": {str(e): int(curve(diag, tau, [e])[0]) for e in (0.6, 0.65, 0.7, 0.75, 0.79)},
        "F2_beta_curve_tail": (list(r["by_field"]["F2"]["beta_curve"].items())[-4:] if "F2" in r["by_field"] else None)}
r4 = E["4"]
k3["N3000_tau0.9"] = {"beta_curve": r4["by_field"]["F3"]["beta_curve"]}
nulls = {}
for r in L("B_null4m_rips.json")["runs"]:
    diag = r["by_field"]["F3"]["diagram_H2"]
    nulls[f"N{r['spec']['N']}_seed{r['spec']['seed']}"] = {
        "max_window_ratio_beta2_ge1": windows_ge1(diag, tau),
        "false_positive_rule_(>=1.5)": windows_ge1(diag, tau) >= 1.5,
        "beta2_at_0.7": int(curve(diag, tau, [0.7])[0])}
k3_fp = {k: windows_ge1(E[k]["by_field"]["F3"]["diagram_H2"], tau) for k in ("0", "1", "2")}
fails = {k: {"space": v["spec"]["space"], "N": v["spec"]["N"], "tau": v["spec"]["tau"], "reason": v["reason"]}
         for k, v in E.items() if v.get("failed")}
k3xt2 = {k: {"N": v["spec"]["N"], "n_simplices": v["info"]["n_simplices"], "sec": v["wall"], "maxrss_MB": v["maxrss_MB"],
             "best_window_ratio": v["by_field"]["F3"]["best_window_ratio"],
             "beta_curve_tail": list(v["by_field"]["F3"]["beta_curve"].items())[-5:]}
         for k, v in E.items() if v["spec"]["space"] == "K3xT2" and not v.get("failed")}

out = {"tier": "X", "criterion": "window ratio >= 1.5 with exact expected Betti vector (pre-registered)",
       "studies": S, "scaling_fit": fit, "B2_injectivity": L("B2_injectivity.json"),
       "B3_K3": k3, "B3_K3_longest_beta2>=1_window_ratio": k3_fp, "B3_K3xT2": k3xt2,
       "B4_null_density_matched": nulls, "explore_failures": fails,
       "explore_orbifold_T4_curves": {"orbifold_N12800": E["15"]["by_field"]["F3"]["beta_curve"],
                                      "T4_N12800": E["16"]["by_field"]["F3"]["beta_curve"]}}
json.dump(out, open(os.path.join(HERE, "B_summary.json"), "w"), indent=1)

# ---------------- figure
BLUE, ORANGE, GRAY = "#2a78d6", "#eb6834", "#8a8a85"
fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
dk = E["0"]["by_field"]["F3"]["diagram_H2"]
dn = [r for r in L("B_null4m_rips.json")["runs"] if r["spec"]["N"] == 4000 and r["spec"]["seed"] == 0][0]["by_field"]["F3"]["diagram_H2"]
for diag, col, mk, lab in ((dn, ORANGE, "x", "density-matched 4-cube null"), (dk, BLUE, "o", "Fermat K3 sample")):
    b = np.array([x[0] for x in diag]); d = np.array([x[1] if x[1] is not None else tau * 1.03 for x in diag])
    ax[0].scatter(b, d, s=10, c=col, marker=mk, label=lab, alpha=0.7, linewidths=0.8)
ax[0].plot([0.2, 0.85], [0.2, 0.85], color=GRAY, lw=1)
ax[0].axhline(tau * 1.03, color=GRAY, lw=0.8, ls=":")
ax[0].text(0.21, tau * 1.03 + 0.005, "alive at truncation tau=0.8", color="#555", fontsize=8)
ax[0].set_xlabel("birth (Rips edge length)"); ax[0].set_ylabel("death")
ax[0].set_title("H2 persistence diagram, F3, N=4000, seed 0", fontsize=10)
ax[0].legend(frameon=False, fontsize=8, loc="lower right")
for key, ls in (("0", "-"), ("1", "--"), ("2", "-.")):
    c = curve(E[key]["by_field"]["F3"]["diagram_H2"], tau, grid)
    ax[1].plot(grid, c, color=BLUE, lw=2, ls=ls, label=f"K3 seed {E[key]['spec']['seed']}")
cn = curve(dn, tau, grid)
ax[1].plot(grid, cn, color=ORANGE, lw=2, label="null seed 0")
ax[1].axhline(22, color=GRAY, lw=1, ls="--"); ax[1].text(0.44, 15, "b2(K3) = 22", color="#555", fontsize=8)
ax[1].set_ylim(0, 120); ax[1].set_xlabel("scale eps (Rips edge length)"); ax[1].set_ylabel("beta_2(eps)")
ax[1].set_title("beta_2 vs scale, F3, N=4000 (tau=0.8)", fontsize=10)
ax[1].legend(frameon=False, fontsize=8)
for a in ax:
    a.spines[["top", "right"]].set_visible(False); a.grid(alpha=0.2)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "B3_K3_H2.png"), dpi=130)
print(json.dumps({k: (v["N_min_all3seeds_F3"], v["largest_N_completed"], v["stopped_because"]) for k, v in S.items()}, indent=0))
print("fit", fit)
print("K3", {k: v.get("beta2_plateau(value,e1,e2,ratio)") for k, v in k3.items()})
print("null FP", {k: round(v["max_window_ratio_beta2_ge1"], 3) for k, v in nulls.items()})
print("K3 beta2>=1 windows", k3_fp)
