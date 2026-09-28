"""
Part 2: Delta*psi_m = G_{m+1}/A for m = -1..3, with G_k the p^k coefficient of the DMVV product, computed by
DIRECT expansion of prod_{r>=1,s>=0,t} (1 - p^r q^s y^t)^(-c(4rs-t^2)) (c = k_norm * cB; v3: no log/exp recurrence).
The coefficients of G_k in the weak-Jacobi monomial basis are SOLVED exactly (never typed):
    G_k = sum_i x_i * E4^i E6^j A^a B^b   (a+b=k, -2a+4i+6j=0)      [ansatz = declared input weak_Jacobi_ring_monomials]
as an over-determined linear system over all (n,l), n<=QCHK. Reported: number of equations, unknowns, rank(M),
rank([M|b]) and the residual (b - M x, exact). m=3 (G_4) is the "DMZ-style line" extension.

Run: cd audit/k3t2_rigidity_v3/B-dyons && <venv-python> part2_psi_m.py 40 7 2
  args: QMAX (theta cache), QCHK (q-order of the identities), K_NORM (must equal the declared input k, asserted)
  (committed run: 40 7 2)
Writes part2_psi_m_results.json.

Negative control (same parameter as the selection: the coefficient vector x): every neighbour x + e_i/d, x - e_i/d
(d = common denominator of the solved x; i.e. the neighbouring integer in the integer-normalised coordinates
d*x) is substituted and the number of failing (n,l) equations is reported.
k-scan: the same solve for k_norm in 1..6 (every integer) -> tells whether the SHAPE selects k (it does not:
weight-0 index-m weak Jacobi forms are closed under the DMVV construction for any multiple k).
"""
import json, sys, itertools
from fractions import Fraction as Fr
from math import lcm
import sympy as sp
from series import mul, add, scal
from common import HERE, load_cache, load2d, load1d, load_cB, load_inputs

KP = 4  # highest p order


def expand_product(cD, QCHK, KP):
    """Direct expansion; keys (kp, n, l) -> int."""
    S = {(0, 0, 0): 1}
    for r in range(1, KP + 1):
        for s in range(0, QCHK + 1):
            top = 4 * r * s
            for t in range(-int((top + 1) ** 0.5) - 1, int((top + 1) ** 0.5) + 2):
                e = cD.get(top - t * t, Fr(0))
                if e == 0:
                    continue
                e = int(e)
                fac = {(0, 0, 0): 1}
                c = Fr(1)
                j = 1
                while r * j <= KP and s * j <= QCHK:
                    c = c * (e + j - 1) / j
                    assert c.denominator == 1
                    if c != 0:
                        fac[(r * j, s * j, t * j)] = int(c)
                    j += 1
                new = {}
                for (a1, b1, l1), v1 in S.items():
                    for (a2, b2, l2), v2 in fac.items():
                        a, b = a1 + a2, b1 + b2
                        if a > KP or b > QCHK:
                            continue
                        key = (a, b, l1 + l2)
                        new[key] = new.get(key, 0) + v1 * v2
                S = {k: v for k, v in new.items() if v != 0}
    return S


def G_of_k(S, k):
    return {(n, l): Fr(v) for (a, n, l), v in S.items() if a == k}


def monomials(m):
    """list of (name, (a,b,i,j)) with a+b=m, -2a+4i+6j=0, i,j>=0."""
    out = []
    for a in range(0, m + 1):
        b = m - a
        w = 2 * a
        for i in range(0, w // 4 + 1):
            rem = w - 4 * i
            if rem >= 0 and rem % 6 == 0:
                out.append((f"E4^{i} E6^{rem // 6} A^{a} B^{b}", (a, b, i, rem // 6)))
    return out


def power(s, k, QCHK):
    out = {(0, 0): Fr(1)}
    for _ in range(k):
        out = mul(out, s, QCHK, None)
    return out


def solve_G(G, m, A, B, E4, E6, QCHK):
    mons = monomials(m)
    E4_2d = {(n, 0): v for n, v in E4.items() if n <= QCHK}
    E6_2d = {(n, 0): v for n, v in E6.items() if n <= QCHK}
    cols = []
    for name, (a, b, i, j) in mons:
        s = mul(power(A, a, QCHK), power(B, b, QCHK), QCHK, None)
        s = mul(s, mul(power(E4_2d, i, QCHK), power(E6_2d, j, QCHK), QCHK, None), QCHK, None)
        cols.append(s)
    keys = sorted(set(G) | set().union(*[set(c) for c in cols]))
    keys = [k for k in keys if k[0] <= QCHK]
    M = sp.Matrix([[sp.Rational(c.get(k, Fr(0)).numerator, c.get(k, Fr(0)).denominator) for c in cols] for k in keys])
    bvec = sp.Matrix([sp.Rational(G.get(k, Fr(0)).numerator, G.get(k, Fr(0)).denominator) for k in keys])
    rM = M.rank()
    rMb = M.row_join(bvec).rank()
    sol = None
    resid_nonzero = None
    if rM == len(mons):
        # least-squares-free exact: solve normal equations (exact), then check residual
        x = (M.T * M).LUsolve(M.T * bvec)
        res = bvec - M * x
        resid_nonzero = sum(1 for v in res if v != 0)
        sol = [Fr(int(v.p), int(v.q)) for v in x]
    return {"monomials": [n for n, _ in mons], "num_equations": len(keys), "num_unknowns": len(mons),
            "rank_M": rM, "rank_augmented": rMb, "consistent": rM == rMb,
            "num_nonzero_residual_equations": resid_nonzero,
            "solution": None if sol is None else [str(v) for v in sol]}, sol, M, bvec, keys


def neighbour_control(sol, M, bvec):
    d = 1
    for v in sol:
        d = lcm(d, v.denominator)
    base = [sp.Rational(v.numerator, v.denominator) for v in sol]
    res = {}
    for i in range(len(sol)):
        for sgn in (+1, -1):
            x = list(base)
            x[i] += sp.Rational(sgn, d)
            r = bvec - M * sp.Matrix(x)
            res[f"x[{i}]{'+' if sgn > 0 else '-'}1/{d}"] = sum(1 for v in r if v != 0)
    return d, res


def compute(QMAX, QCHK, knorm, ms=(1, 2, 3, 4), with_control=True):
    cache = load_cache(QMAX)
    A = load2d(cache, "A_series")
    B = load2d(cache, "B_series")
    E4 = load1d(cache, "E4_series")
    E6 = load1d(cache, "E6_series")
    cB = load_cB(cache)
    assert 4 * KP * QCHK <= 4 * QMAX - 3
    cD = {D: knorm * v for D, v in cB.items()}
    S = expand_product(cD, QCHK, KP)
    out = {}
    for m in ms:
        G = G_of_k(S, m)
        info, sol, M, bvec, keys = solve_G(G, m, A, B, E4, E6, QCHK)
        if sol is not None and with_control:
            d, ctrl = neighbour_control(sol, M, bvec)
            info["neighbour_denominator_d"] = d
            info["neighbour_perturbations_num_failing_equations"] = ctrl
            info["all_neighbours_fail"] = all(v > 0 for v in ctrl.values())
        out[f"G_{m}"] = info
    return out


if __name__ == "__main__":
    QMAX, QCHK, knorm_arg = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    inp = load_inputs()
    knorm = Fr(inp["K3_elliptic_genus_factor_k"]["value"])
    assert knorm == knorm_arg, "third argument must equal the declared input k"
    main = compute(QMAX, QCHK, knorm)
    # Delta psi_m = G_{m+1}/A : Delta psi_{-1} = 1/A (G_0 = 1), Delta psi_0 = G_1/A ...
    scan = {}
    for k in range(1, 7):
        r = compute(QMAX, QCHK, Fr(k), ms=(1, 2), with_control=False)
        scan[str(k)] = {"G_1_solution": r["G_1"]["solution"], "G_1_consistent": r["G_1"]["consistent"],
                        "G_2_solution": r["G_2"]["solution"], "G_2_consistent": r["G_2"]["consistent"],
                        "monomials_G_2": r["G_2"]["monomials"]}
    # expected fields filled AFTER computing (FROM MEMORY / v2 note, tier L)
    expected = {
        "G_2 = (9/4) B^2 + (3/4) E4 A^2  [4 Dpsi_1 = 9B^2/A + 3E4A]": ["9/4", "3/4"],
        "G_3 = (50/27) B^3 + (48/27) E4 A^2 B + (10/27) E6 A^3": ["50/27", "48/27", "10/27"],
        "source": "DMZ arXiv:1208.4074 eq 5.16, FROM MEMORY / v2 note; never used in a computation path",
    }
    out = {
        "args": {"QMAX": QMAX, "QCHK": QCHK, "k_norm": str(knorm)},
        "G_solutions_exact": main,
        "delta_psi_lines": {
            "m=-1": "Delta psi_{-1} = 1/A (G_0=1, empty product)",
            "m=0": "Delta psi_0 = G_1/A, G_1 solved coefficient(s) in main['G_1']",
            "m=1": "Delta psi_1 = G_2/A", "m=2": "Delta psi_2 = G_3/A", "m=3": "Delta psi_3 = G_4/A (DMZ-style extension)",
        },
        "expected_after_computing": expected,
        "k_scan_1_to_6": scan,
        "shared_inputs": ["K3_elliptic_genus_factor_k", "DMVV_product_formula", "Jacobi_theta_definitions",
                          "weak_Jacobi_ring_monomials"],
        "rigidity": {
            "parameter_inserted": "coefficient vector x of G_k in the weak-Jacobi monomial basis",
            "selecting_condition": "G_k(DMVV, own expansion) = sum x_i mono_i for all (n,l), n<=QCHK (linear system; target-free)",
            "condition_uses_true_value": False,
            "classification": "RIGID_GIVEN_DEFINITION",
            "solution_set": "unique solution when rank_M = num_unknowns (see G_solutions_exact); given declared k and the monomial ansatz",
            "control_perturbs_same_parameter": True,
            "negative_control": "neighbour_perturbations_num_failing_equations (each coordinate +-1/d)",
            "k_dependence": "solvable for every k in 1..6: the SHAPE does not select k (see k_scan)",
        },
    }
    with open(HERE / "part2_psi_m_results.json", "w") as f:
        json.dump(out, f, indent=1)
    for m in (1, 2, 3, 4):
        g = main[f"G_{m}"]
        print(m, g["monomials"], g["solution"], "rank", g["rank_M"], g["rank_augmented"], "resid", g["num_nonzero_residual_equations"], "eqs", g["num_equations"], "allfail", g.get("all_neighbours_fail"))
    for k, v in scan.items():
        print("k", k, v["G_1_solution"], v["G_2_solution"], v["G_2_consistent"])
