"""
Track C, item (3): O(d,d;Z) structure checks, d = 2, 3, exact (sympy Rational),
200 seeded random integer samples per check per d.

eta = [[0, I_d], [I_d, 0]]  (the O(d,d) bilinear form, block form).

Checked to preserve eta (g^T eta g == eta), exactly:
  (i)  Theta-shifts   g = [[I, Theta], [0, I]],  Theta antisymmetric integer.
  (ii) Basis changes  g = [[A, 0], [0, A^{-T}]], A in GL(d,Z) (det A = +-1).
  (iii) Factorized duality K_i (swap coordinate i between the two D-dim blocks,
        identity elsewhere) -- also checked to square to the identity.

Negative control: (i) repeated with Theta SYMMETRIC and explicitly forced
nonzero -- must FAIL to preserve eta.

Generalized metric identity: eta * H(G,0) * eta == H(G^{-1},0), where
H(G,B) = [[G - B G^{-1} B, B G^{-1}], [-G^{-1} B, G^{-1}]]; checked with B=0
so H(G,0) = diag(G, G^{-1}), for random SPD-ish rational G (both integer and
non-integer entries, always exact).
"""
import json
import random
import sys

import sympy as sp

sys.path.insert(0, "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity/lattices-duality")

N_SAMPLES = 200


def eta_matrix(d):
    Z = sp.zeros(d, d)
    I = sp.eye(d)
    return sp.Matrix(sp.BlockMatrix([[Z, I], [I, Z]]))


def random_antisymmetric(d, rng, force_nonzero=False):
    T = sp.zeros(d, d)
    for i in range(d):
        for j in range(i + 1, d):
            v = rng.randint(-5, 5)
            T[i, j] = v
            T[j, i] = -v
    if force_nonzero and T == sp.zeros(d, d):
        T[0, 1] = 1
        T[1, 0] = -1
    return T


def random_symmetric_nonzero(d, rng):
    S = sp.zeros(d, d)
    for i in range(d):
        for j in range(i, d):
            v = rng.randint(-5, 5)
            S[i, j] = v
            S[j, i] = v
    if S == sp.zeros(d, d):
        S[0, 0] = 1  # force nonzero, still symmetric
    return S


def random_GLdZ(d, rng, n_moves=8):
    """Build A in GL(d,Z) with det = +-1 exactly, by composing random
    elementary unimodular row operations onto the identity."""
    A = sp.eye(d)
    for _ in range(n_moves):
        kind = rng.choice(["shear", "swap", "negate"])
        if kind == "shear" and d > 1:
            i, j = rng.sample(range(d), 2)
            k = rng.choice([-2, -1, 1, 2])
            A[i, :] = A[i, :] + k * A[j, :]
        elif kind == "swap" and d > 1:
            i, j = rng.sample(range(d), 2)
            A[i, :], A[j, :] = A[j, :].copy(), A[i, :].copy()
        else:
            i = rng.randrange(d)
            A[i, :] = -A[i, :]
    return A


def factorized_duality(d, i):
    """K_i: identity 2d x 2d matrix, but swap coordinate i between the two
    D-dim blocks (row/col i <-> row/col d+i). This is the 'factorized'
    (single-direction) T-duality; eta itself is the product of all K_i."""
    K = sp.eye(2 * d)
    K[i, i] = 0
    K[d + i, d + i] = 0
    K[i, d + i] = 1
    K[d + i, i] = 1
    return K


def random_spd(d, rng):
    M = sp.Matrix(d, d, lambda i, j: rng.randint(-3, 3))
    G = M.T * M + (d + 1) * sp.eye(d)  # SPD, integer entries, exact
    return G


def generalized_metric(G, B):
    Ginv = G.inv()
    top_left = G - B * Ginv * B
    top_right = B * Ginv
    bot_left = -Ginv * B
    bot_right = Ginv
    return sp.Matrix(sp.BlockMatrix([[top_left, top_right], [bot_left, bot_right]]))


def run_for_d(d):
    rng = random.Random(20260918 + d)  # seeded, reproducible
    eta = eta_matrix(d)
    out = {"d": d}

    # (i) theta-shift preserves eta
    theta_pass = 0
    for _ in range(N_SAMPLES):
        T = random_antisymmetric(d, rng)
        g = sp.Matrix(sp.BlockMatrix([[sp.eye(d), T], [sp.zeros(d, d), sp.eye(d)]]))
        if sp.simplify(g.T * eta * g - eta) == sp.zeros(2 * d, 2 * d):
            theta_pass += 1
    out["theta_shift_preserves_eta"] = {"passed": theta_pass, "total": N_SAMPLES,
                                         "all_passed": theta_pass == N_SAMPLES}

    # negative control: symmetric, forced-nonzero Theta must FAIL
    neg_fail_count = 0  # count of samples that correctly FAIL to preserve eta
    for _ in range(N_SAMPLES):
        S = random_symmetric_nonzero(d, rng)
        g = sp.Matrix(sp.BlockMatrix([[sp.eye(d), S], [sp.zeros(d, d), sp.eye(d)]]))
        preserved = sp.simplify(g.T * eta * g - eta) == sp.zeros(2 * d, 2 * d)
        if not preserved:
            neg_fail_count += 1
    out["negative_control_symmetric_theta"] = {
        "num_samples": N_SAMPLES,
        "num_correctly_failing": neg_fail_count,
        "control_passes": neg_fail_count == N_SAMPLES,  # every nonzero-symmetric Theta must break eta
    }

    # (ii) basis change [[A,0],[0,A^{-T}]] preserves eta, A in GL(d,Z)
    basis_pass = 0
    for _ in range(N_SAMPLES):
        A = random_GLdZ(d, rng)
        detA = A.det()
        assert detA in (1, -1), f"sampler produced non-unimodular A, det={detA}"
        AinvT = A.inv().T
        g = sp.Matrix(sp.BlockMatrix([[A, sp.zeros(d, d)], [sp.zeros(d, d), AinvT]]))
        if sp.simplify(g.T * eta * g - eta) == sp.zeros(2 * d, 2 * d):
            basis_pass += 1
    out["basis_change_preserves_eta"] = {"passed": basis_pass, "total": N_SAMPLES,
                                          "all_passed": basis_pass == N_SAMPLES}

    # (iii) factorized duality: preserves eta, and squares to identity
    fact_preserve_pass = 0
    fact_square_pass = 0
    total_fact_checks = 0
    for _ in range(N_SAMPLES):
        i = rng.randrange(d)
        K = factorized_duality(d, i)
        total_fact_checks += 1
        if sp.simplify(K.T * eta * K - eta) == sp.zeros(2 * d, 2 * d):
            fact_preserve_pass += 1
        if sp.simplify(K * K - sp.eye(2 * d)) == sp.zeros(2 * d, 2 * d):
            fact_square_pass += 1
    out["factorized_duality_preserves_eta"] = {
        "passed": fact_preserve_pass, "total": total_fact_checks,
        "all_passed": fact_preserve_pass == total_fact_checks,
    }
    out["factorized_duality_squares_to_identity"] = {
        "passed": fact_square_pass, "total": total_fact_checks,
        "all_passed": fact_square_pass == total_fact_checks,
    }
    # also explicitly check the FULL eta itself (product of all K_i) squares to 1
    eta_sq_is_id = sp.simplify(eta * eta - sp.eye(2 * d)) == sp.zeros(2 * d, 2 * d)
    out["full_eta_squares_to_identity"] = bool(eta_sq_is_id)

    # generalized metric identity: eta H(G,0) eta == H(G^{-1},0)
    gm_pass = 0
    for _ in range(N_SAMPLES):
        G = random_spd(d, rng)
        Ginv = G.inv()
        H_G = generalized_metric(G, sp.zeros(d, d))
        H_Ginv = generalized_metric(Ginv, sp.zeros(d, d))
        lhs = eta * H_G * eta
        if sp.simplify(lhs - H_Ginv) == sp.zeros(2 * d, 2 * d):
            gm_pass += 1
    out["generalized_metric_identity_eta_H_eta_eq_H_Ginv"] = {
        "passed": gm_pass, "total": N_SAMPLES, "all_passed": gm_pass == N_SAMPLES,
    }

    return out


results = {
    "track": "C",
    "item": "3_oddz_checks",
    "method": "Exact sympy Rational/Integer arithmetic; eta=[[0,I],[I,0]] block form; "
               "200 seeded random samples per check per d (seed=20260918+d, "
               "python random.Random, reproducible); GL(d,Z) elements built by "
               "composing random unimodular row operations onto the identity "
               "(guarantees exact det=+-1, not sampled-and-filtered).",
    "n_samples_per_check": N_SAMPLES,
    "results_by_d": {str(d): run_for_d(d) for d in (2, 3)},
}

out_path = "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity/lattices-duality/03_oddz_checks_result.json"
with open(out_path, "w") as f:
    json.dump(results, f, indent=2)

print(json.dumps(results, indent=2))
