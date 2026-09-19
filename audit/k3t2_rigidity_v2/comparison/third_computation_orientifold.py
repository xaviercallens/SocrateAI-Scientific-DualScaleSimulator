"""Third, minimal computation for the DISAGREE rows of the sealed-v2 comparison
(orientifold / D3-tadpole block).

Run (from a clean checkout, no arguments):
  /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python \
      audit/k3t2_rigidity_v2/comparison/third_computation_orientifold.py
(the script resolves every path relative to its own location, so the working
directory does not matter).  Writes third_computation_orientifold_results.json
next to itself.

What it computes (exact integer arithmetic only):
  (1) Fixed points of the involution x -> -x on the real torus R^n / Z^n, for
      n = 2 (the T^2 of K3 x T^2/Z2) and n = 4 (the T^4 of the Kummer T^4/Z2),
      by brute-force enumeration of the half-period grid (1/2)Z^n / Z^n and a
      direct test 2x == 0 mod Z^n.  No closed form 2^n is used.
  (2) The 7-brane (D7) charge budget of the K3 x T^2/Z2 orientifold: the O7
      planes sit at the fixed points of the involution on T^2 (they wrap K3), so
      their number is the n=2 count from (1).  With the tier-L per-plane charge
      Q(O7^-) = -4 in D7 units (input, stated as such), the number of D7 branes
      needed for cancellation is solved for, not assumed.
  (3) The D3 tadpole target chi(X)/24 for two candidate X: X = K3 (what the
      sealed LeanMaster statement d3_tadpole_target_is_one uses) and
      X = K3 x K3 (the M-theory / F-theory 4-fold dual of the K3 x T^2/Z2
      orientifold, tier L).  chi(K3 x K3) is computed by the Kunneth formula
      from Track D's BLIND GUDHI Betti numbers of K3 (read from
      D-tda/results.json), not typed.
  (4) Negative control: using the T^4 fixed-point count (16) for the O7 count
      reproduces the sealed LeanMaster TadpoleCancellation totals -64 / +64
      only if the D7 count is doubled to 32 AND the D7 unit charge is 2 --
      i.e. the sealed numbers are the T^4 (wrong torus) bookkeeping.

Tier: B for the arithmetic; L for the physics identifications (O7 at T^2 fixed
points, Q(O7^-) = -4, the F-theory dual on K3 x K3).  Nothing here is 'proved'.
"""
import itertools
import json
import os
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
V2 = os.path.dirname(HERE)


def fixed_points_of_negation(n):
    """Enumerate x in (1/2)Z^n / Z^n (a superset of all candidates, since
    -x == x mod Z^n  <=>  2x in Z^n) and keep those with 2x == 0 mod Z^n.
    Also scan a finer grid (1/4)Z^n/Z^n to show no other point is fixed."""
    fixed = []
    grid = [Fraction(k, 4) for k in range(4)]
    for x in itertools.product(grid, repeat=n):
        if all(((-xi) - xi) % 1 == 0 for xi in x):
            fixed.append(tuple(str(xi) for xi in x))
    return fixed


def kunneth_betti(b_x, b_y):
    out = [0] * (len(b_x) + len(b_y) - 1)
    for i, bx in enumerate(b_x):
        for j, by in enumerate(b_y):
            out[i + j] += bx * by
    return out


def euler(b):
    return sum((-1) ** k * v for k, v in enumerate(b))


def main():
    res = {"script": os.path.relpath(__file__, V2),
           "command": "/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python "
                      "audit/k3t2_rigidity_v2/comparison/third_computation_orientifold.py"}

    fp2 = fixed_points_of_negation(2)
    fp4 = fixed_points_of_negation(4)
    res["fixed_points_T2"] = {"count": len(fp2), "points": fp2,
                              "grid_scanned": "(1/4)Z^2/Z^2 (16 points)"}
    res["fixed_points_T4"] = {"count": len(fp4), "grid_scanned": "(1/4)Z^4/Z^4 (256 points)"}

    # (2) 7-brane charge
    q_o7 = -4  # tier-L INPUT: O7^- charge in D7 units
    n_o7 = len(fp2)
    total_o7 = n_o7 * q_o7
    q_d7 = 1   # unit
    n_d7_needed = Fraction(-total_o7, q_d7)
    res["seven_brane_charge"] = {
        "inputs_tier_L": {"Q_O7_minus_in_D7_units": q_o7, "Q_D7": q_d7,
                          "O7_location": "fixed points of z->-z on T^2 (O7 wraps K3)"},
        "n_O7_computed": n_o7,
        "total_O7_charge_D7_units": total_o7,
        "n_D7_needed_for_cancellation_computed": str(n_d7_needed),
        "total_D7_charge_D7_units": int(n_d7_needed) * q_d7,
        "net": total_o7 + int(n_d7_needed) * q_d7,
    }
    # (4) negative control: wrong torus
    n_o7_wrong = len(fp4)
    res["negative_control_wrong_torus_T4"] = {
        "n_O7_if_T4_fixed_points_used": n_o7_wrong,
        "total_O7_charge_D7_units": n_o7_wrong * q_o7,
        "D7_needed_in_unit_charge": -n_o7_wrong * q_o7,
        "reading": "the sealed LeanMaster TadpoleCancellation value -64 equals the T^4 (Kummer) fixed-point count x (-4); "
                   "the orientifold of K3 x T^2 has its O7 planes at the T^2 fixed points",
    }

    # (3) D3 tadpole target from Track D's blind Betti numbers
    with open(os.path.join(V2, "D-tda", "results.json")) as f:
        dres = json.load(f)
    bk3 = [dres["resolved_K3_betti_numbers"]["N=6"][f"b{k}"] for k in range(5)]
    bk3_8 = [dres["resolved_K3_betti_numbers"]["N=8"][f"b{k}"] for k in range(5)]
    assert bk3 == bk3_8
    chi_k3 = euler(bk3)
    b_k3k3 = kunneth_betti(bk3, bk3)
    chi_k3k3 = euler(b_k3k3)
    res["d3_tadpole"] = {
        "input_betti_K3_from": "D-tda/results.json resolved_K3_betti_numbers (N=6 == N=8)",
        "betti_K3": bk3,
        "chi_K3": chi_k3,
        "betti_K3xK3_kunneth": b_k3k3,
        "chi_K3xK3": chi_k3k3,
        "chi_K3xK3_mod_24": chi_k3k3 % 24,
        "target_chi_K3xK3_over_24": Fraction(chi_k3k3, 24).__str__(),
        "target_chi_K3_over_24_as_in_sealed_lean": Fraction(chi_k3, 24).__str__(),
        "orientifold_D3_bookkeeping_tier_L": {
            "induced_D3_per_O7": 2, "induced_D3_per_D7": 1,
            "total": n_o7 * 2 + int(n_d7_needed) * 1,
            "note": "per-object induced charges are tier-L inputs (TT lines 164-171); the counts n_O7, n_D7 are the ones computed above"},
    }
    res["verdicts"] = {
        "d3_target_from_K3xK3": int(Fraction(chi_k3k3, 24)),
        "d3_target_from_bookkeeping": n_o7 * 2 + int(n_d7_needed),
        "two_routes_agree": int(Fraction(chi_k3k3, 24)) == n_o7 * 2 + int(n_d7_needed),
        "sealed_lean_d3_target_is_one_reproduced_only_by_chi_K3_over_24": Fraction(chi_k3, 24) == 1,
    }
    out = os.path.join(HERE, "third_computation_orientifold_results.json")
    with open(out, "w") as f:
        json.dump(res, f, indent=1)
    print(json.dumps(res["verdicts"]))
    print(json.dumps(res["seven_brane_charge"]))


if __name__ == "__main__":
    main()
