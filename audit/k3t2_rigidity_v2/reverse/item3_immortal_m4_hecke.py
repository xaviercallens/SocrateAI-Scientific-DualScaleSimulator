"""
Reverse-pass item 3: LeanMaster's `ImmortalHigher.lean` proves the DMZ (9.11)/(9.13)
identity `(psi_{0,m+1}^opt - 12^{m+1}*A*A_{2,m})/A = -12^{m+1}*(H|V_m)` for m = 2, 3 (both
PRIME), where `H|V_m` is the naive Hecke-like operator `12H(D)+m*12H(D/m^2)*[m^2|D]`
(Lean's `hurwitzVSer`). Lean's own docstring: "Not proved: m>=4; the general statement
(9.13), DMZ Sec.10" -- explicitly flagged as open because m=4 is NOT prime.

This is a DIFFERENT statement from the one v1 already tested (v1's P3/P4 extended the
CRUDER pole-coefficient identity `polarDefect`/`immortal_exact_m23`, i.e.
`G_{m+1}(tau,0) = p24(m+1)`, to m=4,5 -- excluded here per the task).

WHAT THE HECKE-LIKE OPERATOR SHOULD BE FOR m=4 (non-prime, m=2^2): DMZ Sec.4.4, eq.
(4.37), define V_{k,t}: c(phi|V_{k,t};n,r) = sum_{d | gcd(n,r,t)} d^{k-1} c(phi; nt/d^2, r/d).
For prime m the divisors of gcd(n,r,m) are only {1,m} (when m | gcd(n,r,m)) or {1}, which
is exactly Lean's two-term formula. For m=4 the divisors of gcd(n,r,4) can also be {1,2,4}
or {1,2} -- an intermediate d=2 term that Lean's/the naive two-term "prime pattern" MISSES.
DMZ's paper (`papers/foundations/1208_4074.txt`, l.3851 area, quoted below) states this
explicitly and gives the m=4 fix:

    Phi_{2,4}^opt = -H|V_4 + 2*H|U_2                                    (quoted, Tier L)

where U_2 (eq. 4.36) sends c(n,r) -> written at (n, 2r) (an index-raising re-embedding, NOT
a divisor sum). This script tests THIS statement, independently, against the K3-elliptic-
genus/DMVV-computed side, using Table 2 of the same paper (l.3899-3907) for the companion
weak Jacobi form:

    psi_{0,5}^opt = B^5 - 10*E4*A^2*B^3 + 20*E6*A^3*B^2 - 15*E4^2*A^4*B + 4*E4*E6*A^5
                                                                          (quoted, Tier L)

PREDICTION (fixed before `computed` is touched; sign fixed against the m=2,3 Lean pattern,
where `dmz_911_verified` gives `(psi0,3opt-12^3 A A2,2)/A = -144*hurwitzVSer(2,2,true) =
-12^3*(H|V2) = +12^3*Phi_2,2opt` since `Phi_2,2opt=-H|V2`, i.e. the LHS is `+12^{m+1}*Phi`,
not `-12^{m+1}*Phi`):
  (psi_{0,5}^opt - 12^5*A*A_{2,4})/A  ==  +12^5 * Phi_{2,4}^opt
                                       ==  -12^5 * (H|V4) + 2*12^5 * (H|U2)
computed here from G_5 = [p^5] of the DMVV product (independent of the psi_{0,5}^opt/A_2,4
ansatz -- G_5 comes from the K3 elliptic genus alone), through q^QCHK.

NEGATIVE CONTROLS (both structural, not built from any target value):
  (a) the NAIVE two-term "prime pattern" H|V4_naive(n,r) := H(D) + 4*H(D/16)*[16|D] (which
      is what you'd write if you (wrongly) assumed the m=2,3 two-term formula generalises
      literally to m=4, omitting the d=2 divisor) must FAIL to match the computed side at
      some (n,l) where the d=2 term is nonzero but the d=4 term is not (i.e. 2 | gcd(n,r,4)
      but 4 does not).
  (b) the paper's OWN alternative formula (9.14), Phi_{2,m}^opt = -2*H|V^{(1)}_{2,m} for
      prime powers, expands (via (4.38), m1=1 in V^{(1)}) to -2*H|V4 + 4*H|U2 for m=4 --
      DIFFERENT by a factor of 2 in each term from the m=4-specific formula above. Both are
      tested against the computed side; this is reported honestly as a genuine textual
      ambiguity/typo risk in the source, resolved by which one the independent computation
      agrees with (agreement of two independent computations, not literal fiat).

Data: reads B-dyons/theta_forms_cache.json (QMAX=40; regenerate with
`cd ../B-dyons && python theta_forms.py 40`) and
B-dyons/hurwitz_class_numbers_results.json (DMAX=40; regenerate with
`cd ../B-dyons && python hurwitz_class_numbers.py`). Both files are already committed.

Run: /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python item3_immortal_m4_hecke.py
Writes item3_immortal_m4_hecke_results.json.
"""
import json
import sys
from fractions import Fraction as Fr
from math import gcd

sys.path.insert(0, "../B-dyons")
from series import mul, add, scal  # noqa: E402  (n,l)-keyed exact Laurent-series engine

with open("../B-dyons/theta_forms_cache.json") as f:
    cache = json.load(f)
with open("../B-dyons/hurwitz_class_numbers_results.json") as f:
    Hdata_full = json.load(f)
H_TABLE = Hdata_full["H_table"]
DMAX_H = Hdata_full["DMAX"]


def pf(s):
    if s is None:
        return None
    if isinstance(s, (int, float)):
        return Fr(s)
    if "/" in s:
        a, b = s.split("/")
        return Fr(int(a), int(b))
    return Fr(int(s))


def load2d(key):
    out = {}
    for k, v in cache[key].items():
        n, l = k.split(",")
        out[(int(n), int(l))] = pf(v)
    return out


def load1d(key):
    return {int(k): pf(v) for k, v in cache[key].items()}


QMAX_theta = cache["QMAX"]
A = load2d("A_series")
B = load2d("B_series")
E4 = load1d("E4_series")
E6 = load1d("E6_series")
cD = {int(k): pf(v) for k, v in cache["cD_table"].items()}

QCHK = 6     # truncation order for the m=4 identity check
YCAP = 40    # y-window cap
assert QMAX_theta >= 6 * QCHK, (
    f"theta cache QMAX={QMAX_theta} too small: need >=6x QCHK={QCHK} margin for a "
    f"degree-5 identity (part2_psi_m.py found 4x adequate for its degree-3 case; we use a "
    f"larger margin here since G_5 is degree-5 in the DMVV product r-grading)"
)

E4_2d = {(n, 0): v for n, v in E4.items() if n <= QCHK}
E6_2d = {(n, 0): v for n, v in E6.items() if n <= QCHK}


def A_r(r, qchk):
    out = {}
    for s in range(0, qchk + 1):
        tmax = int((4 * r * s + 1) ** 0.5) + 2
        for t in range(-tmax, tmax + 1):
            D = 4 * r * s - t * t
            if D in cD:
                out[(s, t)] = out.get((s, t), Fr(0)) + cD[D]
    return {k: v for k, v in out.items() if v != 0}


def subst_qk_yk(series, k, qchk):
    out = {}
    for (n, l), v in series.items():
        if k * n > qchk:
            continue
        out[(k * n, k * l)] = out.get((k * n, k * l), Fr(0)) + v
    return out


KMAX_P = 5  # need G_5 (m=4 -> G_{m+1}=G_5); extends part2_psi_m.py's KMAX_P=3
L = {n: {} for n in range(1, KMAX_P + 1)}
for r in range(1, KMAX_P + 1):
    Ar = A_r(r, QCHK)
    for k in range(1, KMAX_P // r + 1):
        pk = r * k
        term = scal(Fr(1, k), subst_qk_yk(Ar, k, QCHK))
        L[pk] = add(L[pk], term)

G = {0: {(0, 0): Fr(1)}}
for n in range(1, KMAX_P + 1):
    acc = {}
    for j in range(1, n + 1):
        term = scal(j, mul(L.get(j, {}), G[n - j], QCHK, None))
        acc = add(acc, term)
    G[n] = scal(Fr(1, n), acc)
G5 = G[4 + 1]

# ---- validation: G_5(q^0, y=1) should equal p24(5) = 176256 (Lean's `p24_values`, and
# v1's already-Tier-A range; this is a GATE, not the deliverable -- we do NOT extrapolate
# p24 here, we only check our G_5 engine reproduces the already-established value at q^0).
def p24_via_partitions(n, kmax):
    s = {0: Fr(1)}
    for r in range(1, kmax + 1):
        factor = {}
        jmax = kmax // r
        for j in range(jmax + 1):
            c = Fr(1)
            for i in range(j):
                c *= (24 + i)
            if j > 0:
                fact = 1
                for i in range(1, j + 1):
                    fact *= i
                c /= fact
            factor[r * j] = c
        news = {}
        for e1, v1 in s.items():
            for e2, v2 in factor.items():
                e = e1 + e2
                if e > kmax:
                    continue
                news[e] = news.get(e, Fr(0)) + v1 * v2
        s = news
    return s.get(n, Fr(0))


p24_5 = p24_via_partitions(5, 5)
G5_at_q0_y1 = sum(v for (n, l), v in G5.items() if n == 0)
gate_p24_5 = (G5_at_q0_y1 == p24_5 == 176256)

# ---- psi_{0,5}^opt (Table 2, column "5", quoted Tier L -- see docstring)
B5 = mul(mul(B, B, QCHK), mul(B, mul(B, B, QCHK), QCHK), QCHK)
A2 = mul(A, A, QCHK)
A3 = mul(A2, A, QCHK)
A4 = mul(A2, A2, QCHK)
A5 = mul(A4, A, QCHK)
B2 = mul(B, B, QCHK)
B3 = mul(B2, B, QCHK)
psi05opt = add(
    B5,
    scal(-10, mul(E4_2d, mul(A2, B3, QCHK), QCHK)),
    scal(20, mul(E6_2d, mul(A3, B2, QCHK), QCHK)),
    scal(-15, mul(mul(E4_2d, E4_2d, QCHK), mul(A4, B, QCHK), QCHK)),
    scal(4, mul(mul(E4_2d, E6_2d, QCHK), A5, QCHK)),
)

# ---- A_{2,4} (build_A2m, general m, copied verbatim from part3_polar.py -- generic, not
# specific to m=1 or to any target value)
def build_A2m(m, qchk, ycap):
    out = {}
    s = 1
    while m * s * s + s <= qchk:
        base_n = m * s * s + s
        kmax = (qchk - base_n) // s
        for k in range(kmax + 1):
            n = base_n + s * k
            l = 2 * m * s + 1 + k
            if abs(l) <= ycap:
                out[(n, l)] = out.get((n, l), Fr(0)) + (k + 1)
        s += 1
    for l in range(1, ycap + 1):
        out[(0, l)] = out.get((0, l), Fr(0)) + l
    t = 1
    while m * t * t + t <= qchk:
        base_n = m * t * t + t
        kmax = (qchk - base_n) // t
        for k in range(kmax + 1):
            n = base_n + t * k
            l = -2 * m * t - 1 - k
            if abs(l) <= ycap:
                out[(n, l)] = out.get((n, l), Fr(0)) + (k + 1)
        t += 1
    return {k: v for k, v in out.items() if v != 0}


A24 = build_A2m(4, QCHK, YCAP)
A_times_A24 = mul(A, A24, QCHK, YCAP)

# ---- invA (generic; copied construction pattern from part3_polar.py)
def slice_n(d2d, n):
    return {l: v for (nn, l), v in d2d.items() if nn == n}


def conv_y(s1, s2, ycap):
    out = {}
    for l1, v1 in s1.items():
        for l2, v2 in s2.items():
            l = l1 + l2
            if abs(l) > ycap:
                continue
            out[l] = out.get(l, Fr(0)) + v1 * v2
    return {k: v for k, v in out.items() if v != 0}


A_by_n = {n: slice_n(A, n) for n in range(QCHK + 1)}
invA_by_n = {0: {m: Fr(m) for m in range(1, YCAP + 1)}}
R0 = invA_by_n[0]
for n in range(1, QCHK + 1):
    acc = {}
    for i in range(1, n + 1):
        term = conv_y(A_by_n[i], invA_by_n[n - i], YCAP)
        for l, v in term.items():
            acc[l] = acc.get(l, Fr(0)) + v
    invA_by_n[n] = {l: -v for l, v in conv_y(R0, acc, YCAP).items() if v != 0}


def to_2d(seq_by_n):
    out = {}
    for n, s in seq_by_n.items():
        for l, v in s.items():
            out[(n, l)] = v
    return out


invA = to_2d(invA_by_n)


def mul_by_invA(s2d, qchk):
    """s2d * invA, n-graded convolution (both already truncated at qchk)."""
    s_by_n = {n: slice_n(s2d, n) for n in range(qchk + 1)}
    out = {}
    for n in range(qchk + 1):
        acc = {}
        for i in range(n + 1):
            term = conv_y(s_by_n.get(i, {}), invA_by_n.get(n - i, {}), YCAP)
            for l, v in term.items():
                acc[l] = acc.get(l, Fr(0)) + v
        for l, v in acc.items():
            if v != 0:
                out[(n, l)] = v
    return out


numerator5 = add(G5, scal(-(12 ** 5), A_times_A24))
finite_part5 = mul_by_invA(numerator5, QCHK)  # (G5 - 12^5*A*A24)/A candidate

lhs2 = mul_by_invA(add(psi05opt, scal(-(12 ** 5), A_times_A24)), QCHK)  # (psi05opt-12^5*A*A24)/A

# ---- H(D) generating function and Hecke-like operators
def Hval(D):
    if D < 0:
        return Fr(0)
    if D == 0:
        return Fr(-1, 12)
    if D % 4 not in (0, 3):
        return Fr(0)
    key = str(D)
    if key in H_TABLE:
        return pf(H_TABLE[key]["H"])
    return None  # out of tabulated range


def c_H(n, r):
    D = 4 * n - r * r
    return Hval(D)


def divisors(m):
    return [d for d in range(1, m + 1) if m % d == 0]


def c_H_V4_full(n, r):
    """The CORRECT structural Hecke-like V_4 (eq. 4.37), full divisor sum d|gcd(n,r,4)."""
    g = gcd(gcd(n, abs(r)), 4)
    total = Fr(0)
    used_d2 = False
    for d in divisors(g):
        if (4 * n) % (d * d) != 0 or r % d != 0:
            return None
        n2 = (4 * n) // (d * d)
        r2 = r // d
        v = c_H(n2, r2)
        if v is None:
            return None
        total += d * v
        if d == 2:
            used_d2 = True
    return total, used_d2


def c_H_V4_naive_two_term(n, r):
    """WRONG generalisation: only d in {1,4} (the literal m=2,3-prime pattern), skipping d=2."""
    D = 4 * n - r * r
    v0 = Hval(D)
    if v0 is None:
        return None
    total = v0
    if n % 4 == 0 and r % 4 == 0 and D % 16 == 0:
        v16 = Hval(D // 16)
        if v16 is None:
            return None
        total += 4 * v16
    return total


def c_H_U2(n, r):
    if r % 2 != 0:
        return Fr(0)
    return c_H(n, r // 2)


# ---- build target series over the (n,l) grid where finite_part5 is defined
CONST = 12 ** 5
grid = sorted(set(finite_part5) | set(lhs2))

target_correct = {}   # 12^5*(H|V4_full) - 2*12^5*(H|U2)   [Phi_2,4^opt = -H|V4+2H|U2 quoted formula]
target_naive = {}     # 12^5*(H|V4_naive_two_term) - 2*12^5*(H|U2)   [WRONG generalisation]
target_914 = {}       # 2*12^5*(H|V4_full) - 4*12^5*(H|U2)   [(9.14)-style -2H|V^(1)_{2,4}, factor-2 variant]
skipped_out_of_range = []
for (n, r) in grid:
    vfull = c_H_V4_full(n, r)
    vu2 = c_H_U2(n, r)
    if vfull is None or vu2 is None:
        skipped_out_of_range.append((n, r))
        continue
    vfull_val, used_d2 = vfull
    # sign check against the m=2,3 Lean pattern: dmz_911_verified gives
    # (psi0,3opt-12^3 A A2,2)/A = -144*hurwitzVSer(2,2,true) = -12^3*(H|V2) = +12^3*Phi_2,2opt
    # (Phi_2,2opt=-H|V2 per DMZ (9.11)), i.e. the LHS equals +12^{m+1}*Phi_2,m^opt, not -.
    target_correct[(n, r)] = -CONST * vfull_val + 2 * CONST * vu2       # = +12^5 * Phi_2,4opt
    target_914[(n, r)] = -2 * CONST * vfull_val + 4 * CONST * vu2       # = +12^5 * (9.14)-variant
    vnaive = c_H_V4_naive_two_term(n, r)
    if vnaive is not None:
        target_naive[(n, r)] = -CONST * vnaive + 2 * CONST * vu2

keep = [k for k in grid if k not in skipped_out_of_range and k[0] <= QCHK - 1 and abs(k[1]) <= 15]


def compare(lhs, target, keys):
    mism = {}
    for k in keys:
        a = lhs.get(k, Fr(0))
        b = target.get(k, Fr(0))
        if a != b:
            mism[str(k)] = [str(a), str(b)]
    return mism


mism_G5_correct = compare(finite_part5, target_correct, keep)
mism_G5_naive = compare(finite_part5, target_naive, keep)
mism_G5_914 = compare(finite_part5, target_914, keep)
mism_psi_correct = compare(lhs2, target_correct, keep)

holds_correct_via_G5 = (len(mism_G5_correct) == 0)
holds_correct_via_psi = (len(mism_psi_correct) == 0)
holds_naive = (len(mism_G5_naive) == 0)
holds_914 = (len(mism_G5_914) == 0)

# ---- negative control (a): find a (n,r) where the naive two-term pattern differs from the
# full divisor sum (i.e. the d=2 term is genuinely nonzero there) -- demonstrates the two
# formulas are not accidentally identical over the tested window.
d2_matters_examples = []
for (n, r) in keep:
    vfull = c_H_V4_full(n, r)
    vnaive = c_H_V4_naive_two_term(n, r)
    if vfull is not None and vnaive is not None:
        vfull_val, used_d2 = vfull
        if used_d2 and vfull_val != vnaive:
            d2_matters_examples.append({"n": n, "r": r, "full": str(vfull_val), "naive_two_term": str(vnaive)})

result = {
    "prediction": {
        "claim": "(psi_{0,5}^opt - 12^5*A*A_{2,4})/A == -12^5*Phi_{2,4}^opt == "
                 "12^5*(H|V4_full) - 2*12^5*(H|U2), with H|V4_full the FULL divisor-sum "
                 "Hecke-like operator (d in {1,2,4}, not just {1,4}), computed independently "
                 "from G_5 = [p^5] of the DMVV product (K3 elliptic genus alone).",
        "source": "Dabholkar-Murthy-Zagier arXiv:1208.4074, papers/foundations/1208_4074.txt, "
                   "quoted line 'Phi_2,4^opt = -H|V4+2H|U2' (~l.3851) and Table 2 col '5' "
                   "(~l.3899-3907) for psi_{0,5}^opt.",
        "lean_gap_closed": "DualScaleDyons/ImmortalHigher.lean states m=2,3 only "
                            "('Not proved: m>=4'); this tests the paper's own explicit m=4 fix.",
    },
    "validation_gate": {
        "G5_at_q0_y1": str(G5_at_q0_y1),
        "p24_5_independently_computed": str(p24_5),
        "matches_known_p24_5_176256 (Lean p24_values, v1's already-established range)": gate_p24_5,
    },
    "computed": {
        "QCHK": QCHK, "QMAX_theta": QMAX_theta, "DMAX_hurwitz": DMAX_H,
        "num_keys_compared": len(keep),
        "num_keys_skipped_out_of_hurwitz_range": len(skipped_out_of_range),
        "holds_correct_formula_via_G5_DMVV": holds_correct_via_G5,
        "holds_correct_formula_via_psi05opt_ansatz": holds_correct_via_psi,
        "mismatches_correct_via_G5_sample": dict(list(mism_G5_correct.items())[:10]),
        "holds_naive_two_term_wrong_generalisation": holds_naive,
        "mismatches_naive_sample": dict(list(mism_G5_naive.items())[:10]),
        "holds_paper_eq_9_14_style_factor2_variant": holds_914,
        "honesty_note_on_9_14_ambiguity": (
            "DMZ's own eq. (9.14), Phi_2,m^opt = -2*H|V^{(1)}_{2,m} for prime powers, "
            "expands at m=4 (V^{(1)}_{2,4} = V_{2,4} - 2*U_2 by (4.38) with mu(2)=-1, "
            "s^{k-1}=2^1) to -2*H|V4+4*H|U2 -- twice the m=4-specific formula quoted just "
            "above it in the same paper. We test both; the computation, not either quote, "
            "decides which normalisation is the one that actually matches the K3-elliptic-"
            "genus side."
        ),
    },
    "negative_controls": {
        "naive_two_term_generalisation_fails": (not holds_naive) if len(keep) > 0 else None,
        "d2_term_matters_examples (full formula's d=2 contribution is nonzero and the two "
        "formulas differ there)": d2_matters_examples[:10],
        "num_examples_where_d2_matters": len(d2_matters_examples),
    },
    "tier": "B (arithmetic) / L (the two quoted ansatze themselves)",
    "command": "cd audit/k3t2_rigidity_v2/reverse && "
               "/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python item3_immortal_m4_hecke.py",
    "prerequisite_data_regeneration_commands": [
        "cd audit/k3t2_rigidity_v2/B-dyons && "
        "/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python theta_forms.py 40",
        "cd audit/k3t2_rigidity_v2/B-dyons && "
        "/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python hurwitz_class_numbers.py",
    ],
}

with open("item3_immortal_m4_hecke_results.json", "w") as f:
    json.dump(result, f, indent=1)

print(json.dumps({
    "gate_p24_5_ok": gate_p24_5,
    "holds_correct_via_G5": holds_correct_via_G5,
    "holds_correct_via_psi": holds_correct_via_psi,
    "holds_naive (should be False)": holds_naive,
    "holds_9_14_variant (should differ from correct unless normalisations coincide)": holds_914,
    "num_d2_matters_examples": len(d2_matters_examples),
    "num_keys_compared": len(keep),
}, indent=1))
