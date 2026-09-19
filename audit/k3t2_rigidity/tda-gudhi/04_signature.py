#!/usr/bin/env python
"""
Track D, step (3) continued: signature ingredients b2 = b2^+ + b2^-.

Two independent pieces, both computed:

(A) The Z_2-INVARIANT part of H^2(T^4;Q). Since the action is x -> -x on
    H^1(T^4;Q) = Q^4 (acts as -1), it acts on H^2(T^4;Q) = Lambda^2(Q^4)
    by (-1)*(-1) = +1, i.e. trivially -- so ALL of H^2(T^4;Q) survives to
    the quotient (this is an independent check against the computed
    quotient Betti number b2=6 from 02_invariance_and_quotient.py: dim
    Lambda^2(Q^4) = C(4,2) = 6, computed below by literally building the
    basis, matches).

    The cup product / wedge pairing on this 6-dim space,
        (e_i^e_j) . (e_k^e_l) = sign(perm taking (i,j,k,l)->(1,2,3,4)) if
                                 {i,j,k,l}={1,2,3,4}, else 0
    is built EXACTLY with sympy Rationals/integers (Levi-Civita symbol),
    not typed from memory, and diagonalized. Eigenvalue SIGNS are read
    off as floats (permitted by the ground rules for sign/trace checks);
    the matrix itself is exact.

(B) The 16 exceptional (-2)-curves from resolving the 16 isolated
    A1 (RP^3-linked) singular points. STATED, not computed here: each
    exceptional P^1 has self-intersection -2 and the 16 are mutually
    orthogonal and orthogonal to the pulled-back H^2(T^4) classes -- this
    is the standard local structure of an isolated ODP/A1 resolution
    (one (-2)-curve per point, no interaction between distinct points'
    exceptional divisors since the points are disjoint). We did NOT
    verify orthogonality by direct intersection-form computation (that
    needs the actual resolved 4-manifold's cup product, out of reach of
    the link/complement homology computed in 03); it is asserted as a
    standard fact and flagged as such.

Honesty limits (explicit, per ground rules):
  - This gives the signature of the INTERSECTION FORM OVER Q (rank +
    signs), not the integral unimodular lattice structure. Establishing
    the even unimodular lattice II_{3,19} needs more than eigenvalue
    signs (it needs the actual integral Gram matrix / overlattice
    argument), which is OUT OF SCOPE here and NOT claimed.
  - c1=0 (trivial canonical class) is NOT computed anywhere in this
    track; Noether's formula chi=24 <=> 12 chi(O) with c1=0 is used only
    as a stated consistency identity in 07_rigidity_scan.py, never
    derived.
"""
import json
import itertools
import sys
import time
import sympy as sp

t0 = time.time()

# --- build Lambda^2(Q^4) basis and the wedge/cup pairing matrix, exactly ---
basis = list(itertools.combinations(range(1, 5), 2))  # (i,j), i<j ; 6 elements
n = len(basis)
assert n == 6

def levi_civita_sign(seq):
    """Sign of the permutation `seq` (a rearrangement of 1..4), 0 if repeated."""
    if len(set(seq)) != len(seq):
        return 0
    # count inversions
    s = 0
    seq = list(seq)
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                s += 1
    return (-1) ** s

M = sp.zeros(n, n)
for a, (i, j) in enumerate(basis):
    for b, (k, l) in enumerate(basis):
        quad = (i, j, k, l)
        if len(set(quad)) == 4:
            M[a, b] = sp.Integer(levi_civita_sign(quad))
        else:
            M[a, b] = sp.Integer(0)

# sanity: M should be symmetric (cup product on H^2 of a 4-manifold is symmetric)
is_symmetric = (M == M.T)

# eigen-decomposition, exact
eigen_exact = M.eigenvals()  # dict eigenvalue -> multiplicity, exact (sympy) values
eig_list = []
for val, mult in eigen_exact.items():
    for _ in range(mult):
        eig_list.append(val)

# signs as floats (allowed for sign determination)
eig_signs = []
for val in eig_list:
    fv = float(val)
    eig_signs.append(1 if fv > 1e-9 else (-1 if fv < -1e-9 else 0))

b2_plus_T4_part = sum(1 for s in eig_signs if s > 0)
b2_minus_T4_part = sum(1 for s in eig_signs if s < 0)
b2_zero_T4_part = sum(1 for s in eig_signs if s == 0)

# --- combine with the 16 exceptional (-2)-curves ---
num_exceptional = 16
b2_plus_total = b2_plus_T4_part  # exceptional curves are all negative-definite, add nothing to +
b2_minus_total = b2_minus_T4_part + num_exceptional
b2_total = b2_plus_total + b2_minus_total
signature = b2_plus_total - b2_minus_total

OUT = {
    "invariant_H2_T4_basis_size": n,
    "invariant_H2_T4_basis_size_expected": "C(4,2)=6, computed as len(combinations(range(1,5),2))",
    "cup_product_matrix_is_symmetric": is_symmetric,
    "cup_product_matrix_exact": [[str(M[a, b]) for b in range(n)] for a in range(n)],
    "eigenvalues_exact": [str(v) for v in eig_list],
    "eigenvalue_signs": eig_signs,
    "b2_plus_from_T4_invariant_part_COMPUTED": b2_plus_T4_part,
    "b2_minus_from_T4_invariant_part_COMPUTED": b2_minus_T4_part,
    "b2_zero_eigenvalues_T4_invariant_part": b2_zero_T4_part,
    "cross_check_vs_quotient_b2": {
        "quotient_b2_computed_in_step2": 6,
        "b2_plus_plus_minus_plus_zero": b2_plus_T4_part + b2_minus_T4_part + b2_zero_T4_part,
        "matches": (b2_plus_T4_part + b2_minus_T4_part + b2_zero_T4_part) == 6,
    },
    "num_exceptional_minus2_curves": num_exceptional,
    "exceptional_curves_orthogonality_and_self_intersection_status": "STATED (standard fact for isolated A1/ODP resolutions), NOT computed from an explicit intersection-form computation in this track",
    "b2_plus_total": b2_plus_total,
    "b2_minus_total": b2_minus_total,
    "b2_total_(should_match_resolved_K3_b2_from_step3)": b2_total,
    "b2_total_matches_step3_computation": b2_total == 22,
    "signature_b2plus_minus_b2minus": signature,
    "signature_expected_source": "K3 signature -16, standard fact (from memory, unverified) -- compared here only as a cross-check, not used as an input to the computation above",
    "signature_matches_expected": signature == -16,
    "honesty_limits": [
        "Only the intersection form's SIGNS/rank over Q are computed here, not the integral unimodular lattice II_{3,19}.",
        "Orthogonality and self-intersection -2 of the 16 exceptional curves is a stated standard fact, not independently computed by intersection-form arithmetic in this track.",
        "c1=0 / Noether's formula is not derived anywhere in this track; see 07_rigidity_scan.py for how chi=24 is used only as a stated consistency identity.",
    ],
    "runtime_sec": time.time() - t0,
}

OUT_PATH = "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity/tda-gudhi/04_signature_results.json"
with open(OUT_PATH, "w") as f:
    json.dump(OUT, f, indent=2)

print(json.dumps(OUT, indent=2))
print("\nWrote", OUT_PATH)
