#!/usr/bin/env python
"""
SKEPTIC audit, Track B part 3 (dyons/part3_polar.py).

Two questions:
 (Q1) part3 scans N with M fixed at the literal 648 (line 205) and then M with N
      fixed at the literal 324 (line 250). Is (324,648) unique on the JOINT 2-D
      grid N in 320..328 x M in 640..656 when neither is pinned?
 (Q2) The unexplained remainder -324 at (n,l)=(4,-5), D=-9, was hidden by a
      D>=-4 window. Hypothesis: build_A21's s<=-1 branch drops the
      q^{s^2+s} prefactor of the term q^{s^2+s} y^{2s+1}/(1-q^s y)^2:
        correct:  n = s^2+s - 2s - s*k = s^2 - s - s*k
        part3:    n = -2s - s*k            (identical only for s=-1)
      For s=-2,k=0 part3 puts a +1 at (4,-5); the correct position is (6,-5).
      Test: rebuild A21 with the prefactor and rescan on the FULL window
      (n<=QCHK-2, |l|<=20, no D filter).

Reuses part3's own inputs (theta cache A, B, E4; Hurwitz numbers; its 1/A
recursion and G2/A) by executing part3's source only up to its "Step 5"
marker, so no part3 output file is written. Must be run with cwd=dyons/.
Output: skeptic/s2_polar_joint_scan_results.json
"""
import json
import os
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
DYONS = os.path.join(os.path.dirname(HERE), "dyons")
os.chdir(DYONS)
import sys
sys.path.insert(0, DYONS)

src = open(os.path.join(DYONS, "part3_polar.py")).read()
cut = src.index("# ---- Step 5")
ns = {"__name__": "part3_prefix"}
exec(compile(src[:cut], "part3_polar.py[:Step5]", "exec"), ns)

G2_over_A = ns["G2_over_A"]
A = ns["A"]
E4_2d = ns["E4_2d"]
Hhat = ns["Hhat"]
mul, add, scal = ns["mul"], ns["add"], ns["scal"]
QCHK, YCAP = ns["QCHK"], ns["YCAP"]
A21_part3 = ns["A21"]


def build_A21_corrected(qchk, ycap):
    """sum_s q^{s^2+s} y^{2s+1} (1-q^s y)^{-2}, expanded for |q|<|y|<1."""
    out = {}
    smax = int(qchk ** 0.5) + 3
    for s in range(1, smax + 1):
        for k in range(0, qchk + 1):
            n = s * s + s + s * k
            l = 2 * s + 1 + k
            if n <= qchk and abs(l) <= ycap:
                out[(n, l)] = out.get((n, l), 0) + (k + 1)
    for m in range(1, ycap + 1):  # s = 0
        out[(0, m)] = out.get((0, m), 0) + m
    for s in range(-smax, 0):
        # (1-x)^{-2} = x^{-2} (1-1/x)^{-2}, x = q^s y, |x|>1
        for k in range(0, 4 * qchk + 4):
            n = s * s + s - 2 * s - s * k      # prefactor q^{s^2+s} KEPT
            l = 2 * s + 1 - 2 - k
            if 0 <= n <= qchk and abs(l) <= ycap:
                out[(n, l)] = out.get((n, l), 0) + (k + 1)
    return {k: Fr(v) for k, v in out.items() if v != 0}


A21_fixed = build_A21_corrected(QCHK, YCAP)
diff_A21 = {str(k): [str(A21_part3.get(k, 0)), str(A21_fixed.get(k, 0))]
            for k in sorted(set(A21_part3) | set(A21_fixed))
            if A21_part3.get(k, 0) != A21_fixed.get(k, 0)
            and k[0] <= QCHK - 2 and abs(k[1]) <= 20}

threeE4A = scal(3, mul(E4_2d, A, QCHK, None))


def mismatches(A21, N, M, dfilter):
    lhs = add(G2_over_A, scal(-N, A21))
    rhs = add(threeE4A, scal(-M, Hhat))
    keys = {k for k in set(lhs) | set(rhs)
            if k[0] <= QCHK - 2 and abs(k[1]) <= 20
            and (not dfilter or 4 * k[0] - k[1] ** 2 >= -4)}
    return {k: (lhs.get(k, 0), rhs.get(k, 0)) for k in keys
            if lhs.get(k, 0) != rhs.get(k, 0)}


res = {"QCHK": QCHK, "YCAP": YCAP,
       "A21_part3_vs_corrected_differences_in_window": diff_A21}
for label, A21 in (("part3_A21", A21_part3), ("corrected_A21", A21_fixed)):
    for dfilter in (True, False):
        sols, counts = [], {}
        for N in range(320, 329):
            for M in range(640, 657):
                mm = mismatches(A21, N, M, dfilter)
                counts[f"{N},{M}"] = len(mm)
                if not mm:
                    sols.append([N, M])
        key = f"{label}__{'D>=-4_window' if dfilter else 'full_window_no_D_filter'}"
        best = min(counts, key=counts.get)
        entry = {"joint_solutions": sols, "grid_points": len(counts),
                 "min_mismatch_point": best, "min_mismatch_count": counts[best]}
        if not sols:
            mm = mismatches(A21, 324, 648, dfilter)
            entry["mismatches_at_324_648"] = {str(k): [str(a), str(b)] for k, (a, b) in mm.items()}
        # negative controls: nearest neighbours of the solution
        entry["mismatch_counts_neighbours"] = {p: counts[p] for p in
                                               ("323,648", "325,648", "324,647", "324,649", "323,647", "325,649")}
        res[key] = entry

out = os.path.join(HERE, "s2_polar_joint_scan_results.json")
with open(out, "w") as f:
    json.dump(res, f, indent=1)
print(json.dumps(res, indent=1))
