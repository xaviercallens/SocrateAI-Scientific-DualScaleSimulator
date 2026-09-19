"""Skeptic (math lens): Hurwitz class numbers by a route independent of Track B's reduced-form count.

Route (tier L inputs: Dirichlet's class number formula for a fundamental discriminant D0<0,
 h(D0) = -(w/(2|D0|)) * sum_{a=1}^{|D0|} chi_{D0}(a) a, with chi_{D0} the Kronecker symbol;
 and the conductor formula H(|D0| f^2) = (2 h(D0)/w(D0)) * sum_{d|f} mu(d) chi_{D0}(d) sigma_1(f/d);
 convention H(0) = -1/12, H(n) = 0 for n < 0 or n = 1,2 mod 4):
Kronecker symbol computed from scratch (Jacobi reciprocity + (D/2) rule). Exact Fractions.
Then the Lean-side combinations of rows B-dyons#13 (h12 = 12H), #18, #19 are evaluated, and only AFTER
computing, Track B's table and the sealed Lean lists are loaded into 'expected' fields for comparison.
Negative controls: (i) dropping the 2/w weight (plain class number h instead of Hurwitz H) changes H(3), H(4);
(ii) dropping the conductor correction mu(d)chi(d) (sum sigma_1(f) only) changes H(12), H(28).
Command: cd audit/k3t2_rigidity_v2/skeptic_math && \
  /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python skB_hurwitz_dirichlet.py
"""
import json, os
from fractions import Fraction as Fr


def jacobi(a, n):  # n odd positive
    a %= n; r = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5): r = -r
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3: r = -r
        a %= n
    return r if n == 1 else 0


def kronecker(D, m):  # (D/m) for m >= 1, D = 0,1 mod 4
    r = 1
    while m % 2 == 0:
        m //= 2
        if D % 2 == 0: return 0
        r *= 1 if D % 8 in (1, 7) else -1
    return r * jacobi(D, m) if m > 1 else r


def is_fundamental(D):  # D < 0
    if D % 4 == 1:
        m = -D
        return all(m % (p * p) for p in range(2, int(m ** 0.5) + 1))
    if D % 4 == 0:
        m = -D // 4
        if m % 4 not in (1, 2): return False
        return all(m % (p * p) for p in range(2, int(m ** 0.5) + 1))
    return False


def w_of(D, uniform=False):
    if uniform: return 2
    return {-3: 6, -4: 4}.get(D, 2)


def h_fund(D0, uniform=False):
    s = sum(kronecker(D0, a) * a for a in range(1, -D0 + 1))
    return Fr(-w_of(D0, uniform) * s, 2 * (-D0))


def mu(n):
    r, p = 1, 2
    while p * p <= n:
        if n % p == 0:
            n //= p
            if n % p == 0: return 0
            r = -r
        p += 1
    return -r if n > 1 else r


def sigma1(n):
    return sum(d for d in range(1, n + 1) if n % d == 0)


def H(n, uniform=False, nocond=False):
    if n == 0: return Fr(-1, 12)
    if n < 0 or n % 4 in (1, 2): return Fr(0)
    D = -n
    for f in range(int(n ** 0.5), 0, -1):  # largest f with D/f^2 fundamental
        if n % (f * f) == 0 and (D // (f * f)) % 4 in (0, 1) and is_fundamental(D // (f * f)):
            D0 = D // (f * f); break
    pref = h_fund(D0) if uniform else 2 * h_fund(D0) / w_of(D0)
    if nocond: return pref * sigma1(f)
    return pref * sum(mu(d) * kronecker(D0, d) * sigma1(f // d) for d in range(1, f + 1) if f % d == 0)


out = {"H_0_to_40": {n: str(H(n)) for n in range(0, 41)}}
h12 = lambda D: 12 * H(D)
row13 = [h12(D) for D in [0, 3, 4, 7, 8, 11, 12, 15]]
row18 = [-(h12(D) + (2 * h12(D // 4) if D >= 0 and D % 4 == 0 else 0)) for D in [-4, -1, 0, 4, 7, 8, 12, 15, 16]]
row19 = [-(h12(D) + (3 * h12(D // 9) if D >= 0 and D % 9 == 0 else 0)) for D in [-9, -4, -1, 0, 3, 8, 11, 12, 15]]
out["computed"] = {"B-dyons#13_h12": [str(x) for x in row13],
                   "B-dyons#18_dmz911": [str(x) for x in row18],
                   "B-dyons#19_dmz912": [str(x) for x in row19]}
out["negative_control_no_2_over_w"] = {"H(3)": str(H(3, True)), "H(4)": str(H(4, True)),
                                       "differs_from_correct": H(3, True) != H(3) and H(4, True) != H(4)}
out["negative_control_no_conductor_correction"] = {"H(12)": str(H(12, nocond=True)), "H(28)": str(H(28, nocond=True)),
                                       "differs_from_correct": H(12, nocond=True) != H(12) and H(28, nocond=True) != H(28)}
# ---- expected values loaded only after computing ----
here = os.path.dirname(os.path.abspath(__file__))
bt = json.load(open(os.path.join(here, "..", "B-dyons", "hurwitz_class_numbers_results.json")))
def find(d, key):
    if isinstance(d, dict):
        for k, v in d.items():
            if k == key: return v
            r = find(v, key)
            if r is not None: return r
    return None
tabB = find(bt, "H_table")
agreeB = None
if isinstance(tabB, dict):
    agreeB = all(Fr(str(tabB[k]["H"])) == H(int(k)) for k in tabB if int(k) <= 40)
    out["trackB_table_size_compared"] = sum(1 for k in tabB if int(k) <= 40)
out["expected"] = {
    "trackB_H_table_agrees_0_40": agreeB,
    "trackB_source": "B-dyons/hurwitz_class_numbers_results.json:H_table",
    "lean_sealed_B-dyons#13": [-1, 4, 6, 12, 12, 12, 16, 24],
    "lean_sealed_B-dyons#18": [0, 0, 3, -6, -12, -12, -24, -24, -30],
    "lean_sealed_B-dyons#19": [0, 0, 0, 4, -4, -12, -12, -16, -24],
    "lean_source": "comparison.json rows B-dyons#13/#18/#19 lean_value (sealed v3.20.0)",
}
out["match"] = {"#13": row13 == out["expected"]["lean_sealed_B-dyons#13"],
                "#18": row18 == out["expected"]["lean_sealed_B-dyons#18"],
                "#19": row19 == out["expected"]["lean_sealed_B-dyons#19"]}
json.dump(out, open(os.path.join(here, "skB_hurwitz_dirichlet_results.json"), "w"), indent=1)
print(json.dumps({k: out[k] for k in ("computed", "negative_control_no_2_over_w", "negative_control_no_conductor_correction", "match")}, indent=1),
      "trackB_agree", agreeB)
