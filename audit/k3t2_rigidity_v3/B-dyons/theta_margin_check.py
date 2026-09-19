"""
Truncation-margin check: rerun the Jacobi-form construction at QMAX=40 and QMAX=60 (caches from
`theta_forms.py 40` and `theta_forms.py 60`) and compare (a) A, B, E4, E6 coefficients at n <= QCHK,
(b) cB(D) for D <= 4*4*QCHK, (c) the solved G_1..G_4 coefficient vectors (part2) and the (N,c,M) solution (part3).

Run: cd audit/k3t2_rigidity_v3/B-dyons && <venv-python> theta_margin_check.py 40 60 8
  args: QMAX_small QMAX_large QCHK   (committed run: 40 60 8; part2 uses QCHK=7 there, part3 QCHK=8; both are compared at QCHK passed here / QCHK-1)
Writes theta_margin_check_results.json.
"""
import json, sys
from fractions import Fraction as Fr
from common import HERE, load_cache, load2d, load1d, load_cB, pf, load_inputs
import part2_psi_m as P2
import part3_immortal as P3

if __name__ == "__main__":
    q1, q2, QCHK = map(int, sys.argv[1:4])
    c1, c2 = load_cache(q1), load_cache(q2)
    same = {}
    for key in ("A_series", "B_series"):
        a, b = load2d(c1, key), load2d(c2, key)
        same[key] = all(a.get((n, l), Fr(0)) == b.get((n, l), Fr(0)) for n in range(QCHK + 1) for l in range(-40, 41))
    for key in ("E4_series", "E6_series"):
        a, b = load1d(c1, key), load1d(c2, key)
        same[key] = all(a.get(n, Fr(0)) == b.get(n, Fr(0)) for n in range(QCHK + 1))
    ca, cb = load_cB(c1), load_cB(c2)
    Dmax = 16 * QCHK
    same["cB_D_le_16QCHK"] = all(ca.get(D) == cb.get(D) for D in range(-3, Dmax + 1))
    knorm = Fr(load_inputs()["K3_elliptic_genus_factor_k"]["value"])
    r1 = P2.compute(q1, QCHK - 1, knorm, with_control=False)
    r2 = P2.compute(q2, QCHK - 1, knorm, with_control=False)
    same["part2_G1_to_G4_solutions_equal"] = all(r1[k]["solution"] == r2[k]["solution"] for k in r1)
    with open(HERE / "hurwitz_results.json") as f:
        Hh = {int(D): pf(v) for D, v in json.load(f)["H_table"].items()}
    s1 = P3.run(q1, QCHK, 20, 60, knorm, Hh, control=False)
    s2 = P3.run(q2, QCHK, 20, 60, knorm, Hh, control=False)
    same["part3_NcM_equal"] = s1["solution_N_c_M"] == s2["solution_N_c_M"]
    out = {"args": {"QMAX_small": q1, "QMAX_large": q2, "QCHK": QCHK}, "agree": same, "all_agree": all(same.values()),
           "part2_solutions": {k: r1[k]["solution"] for k in r1}, "part3_solution": s1["solution_N_c_M"]}
    with open(HERE / "theta_margin_check_results.json", "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))
