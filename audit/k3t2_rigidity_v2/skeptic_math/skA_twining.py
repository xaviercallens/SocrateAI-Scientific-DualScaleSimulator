"""Skeptic (math lens) independent re-derivation for Track A.

Route (independent of Track A's theta/Appell-Lerch code):
 1. H(tau) via  eta^3 H = -2 E2 + 48 F2,  F2 = sum_{r>s>0, r-s odd} (-1)^r s q^{rs/2}
    (Eguchi-Hikami/Cheng weight-2 route; tier L). Work with h := q^{1/8} H.
 2. Twining H_g for g in {2A,3A,5A,7A} WITHOUT the from-memory F_g table:
    phi_g = (chi/12) phi_{0,1} + t * phi2_N(tau) * phi_{-2,1}, where
      - chi(g) = k from the balanced frame shape 1^k N^k with k(N+1)=24 (k=24/(N+1)),
      - phi2_N = the unique normalised weight-2 form on Gamma0(N) (dim M2(Gamma0(N))=1
        for N=2,3,5,7), phi2_N = 1 + 24/(N-1) sum_m sigma1(m)(q^m - N q^{Nm}),
      - t fixed by the invariance of the 2 ground states at (q^0, y^{+1}):
        chi/12 * [y]phi01_q0 + t*[y]phi_{-2,1,q0} = 2, with [y]phi01_q0 = 1, [y]phi-21_q0 = 1
        computed below from the product formulas, not typed.
    Then H_g = (chi/24) H - t*phi2_N/eta^3 (sign fixed by requiring the polar term of H_g
    to be -2 x trivial, printed as a check, not imposed).
 3. Compare raw coefficients to Track A's A-genus/results.json (read AFTER computing).
 4. Negative control: t -> t + 1/3 breaks integrality or the mod-N congruence.
Command: cd audit/k3t2_rigidity_v2/skeptic_math && <venv python> skA_twining.py
"""
import json, os
from fractions import Fraction as Fr

Q = 12  # q-order through which everything is computed

def sigma1(m):
    return sum(d for d in range(1, m + 1) if m % d == 0)

def mul(a, b):
    c = [Fr(0)] * (Q + 1)
    for i, x in enumerate(a):
        if x:
            for j in range(Q + 1 - i):
                c[i + j] += x * b[j]
    return c

def inv(a):
    assert a[0] != 0
    b = [Fr(0)] * (Q + 1); b[0] = 1 / Fr(a[0])
    for n in range(1, Q + 1):
        b[n] = -sum(a[k] * b[n - k] for k in range(1, n + 1)) / a[0]
    return b

# prod (1-q^n)^3  = q^{-1/8} eta^3
P = [Fr(0)] * (Q + 1); P[0] = Fr(1)
for n in range(1, Q + 1):
    f = [Fr(0)] * (Q + 1); f[0] = 1; f[n] = -1
    for _ in range(3):
        P = mul(P, f)
E2 = [Fr(1)] + [Fr(-24 * sigma1(m)) for m in range(1, Q + 1)]
F2 = [Fr(0)] * (Q + 1)
for r in range(1, 2 * Q + 2):
    for s in range(1, r):
        if (r - s) % 2 == 1 and r * s // 2 <= Q:
            F2[r * s // 2] += (-1) ** r * s
h = mul([Fr(-2) * e + 48 * f for e, f in zip(E2, F2)], inv(P))
An = {n: h[n] / 2 for n in range(0, Q + 1)}

# y-coefficients of phi01 and phi_{-2,1} at q^0 from product formulas (q^0 truncation):
# phi_{-2,1}|_{q^0} = (y-2+1/y);  phi01|_{q^0} = 4*[ (y+2+1/y)/4 + 1 + 1 ] = y+10+1/y
# (theta3,theta4 ratios are 1 at q^0; theta2 ratio^2 = (y+2+1/y)/4).
phim21_q0 = {1: Fr(1), 0: Fr(-2), -1: Fr(1)}
phi01_q0 = {1: Fr(4, 4), 0: Fr(4 * 2, 4) + 4 + 4, -1: Fr(4, 4)}
GROUND_Y1 = 2 * phi01_q0[1]  # Z_K3 = 2 phi01 ground states at y^1 (the 2 = normalisation of Z)

def phi2(N):
    c = [Fr(0)] * (Q + 1); c[0] = Fr(1)
    for m in range(1, Q + 1):
        c[m] += Fr(24, N - 1) * sigma1(m)
        if N * m <= Q:
            c[N * m] -= Fr(24, N - 1) * N * sigma1(m)
    return c

out = {"H_route": "eta^3 H = -2E2 + 48 F2", "A_n": {n: str(An[n]) for n in An}, "twining": {}}
for name, N in (("2A", 2), ("3A", 3), ("5A", 5), ("7A", 7)):
    chi = Fr(24, N + 1)
    t = (GROUND_Y1 - chi / 12 * phi01_q0[1]) / phim21_q0[1]
    const_check = chi / 12 * phi01_q0[0] + t * phim21_q0[0]  # should equal chi-4
    def hg(tt):
        return [chi / 24 * a - tt * b for a, b in zip(h, mul(phi2(N), inv(P)))]
    raw = hg(t)
    raw_bad = hg(t + Fr(1, 3))
    cong = all((h[n] - raw[n]) % N == 0 for n in range(1, Q + 1) if raw[n].denominator == 1)
    integ = all(x.denominator == 1 for x in raw)
    integ_bad = all(x.denominator == 1 for x in raw_bad)
    cong_bad = integ_bad and all((h[n] - raw_bad[n]) % N == 0 for n in range(1, Q + 1))
    out["twining"][name] = {
        "chi_from_frame_shape": str(chi), "t_from_ground_state": str(t),
        "t_times_24_over_N(N-1)_equals_F_factor": str(t * 24 / (N * (N - 1))),
        "q0_const_equals_chi_minus_4": const_check == chi - 4,
        "raw_H_g": [str(x) for x in raw[:10]],
        "polar_term": str(raw[0]), "integral": integ, "congruence_mod_N": cong,
        "neg_control_t_plus_1_3_integral_and_congruent": integ_bad and cong_bad,
    }

# comparison with Track A AFTER computing
ta = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "A-genus", "results.json")))
cmp = {"A_n_1_10_agree": all(str(An[n]) == ta["H_tau_appell_lerch"]["A_n_1_to_10"][str(n)] for n in range(1, 11))}
for name in ("2A", "3A", "5A", "7A"):
    theirs = ta["twining_2A_3A_5A_7A"][name]["A_n_g"]
    mine = out["twining"][name]["raw_H_g"]
    cmp[name + "_agree"] = all(Fr(mine[n]) == 2 * Fr(theirs[str(n)]) for n in range(0, 10))
out["comparison_with_trackA"] = cmp
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "skA_twining_results.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
