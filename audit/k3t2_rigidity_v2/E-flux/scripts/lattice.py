"""
Shared lattice utilities for Track E (flux vacua on K3 x T2/Z2, Tripathy-Trivedi
hep-th/0301139).

Basis and metric come from Tripathy-Trivedi Appendix A, eqs. (A.6)-(A.7)
(source: audit/k3t2_rigidity_v2/sources/hep-th_0301139_TripathyTrivedi.txt,
lines 1849-1868). This is the EXPLICIT T4/Z2 basis {e1,...,e6} of the rank-6
sublattice Gamma_{3,3} subset H^2(K3,Z), diagonal metric

    (e1,e1) = (e2,e2) = (e3,e3) = +2
    (e4,e4) = (e5,e5) = (e6,e6) = -2
    all other inner products zero.

IMPORTANT (flagged, not smoothed over): TT also give a DIFFERENT presentation
of the same abstract lattice Gamma_{3,3} in eqs. (A.2)-(A.3) (lines 1778-1797),
as three hyperbolic planes U (Gram matrix [[0,1],[1,0]] per block, det -1,
unimodular). The diagonal basis (A.6)-(A.7) has Gram matrix diag(2,2,2,-2,-2,-2)
per basis vector, determinant -2^6 = -64 (NOT unimodular). These two Gram
matrices are related by TT's remark "we can make a change of basis such that
H_{3,3} = 2*eta_{3,3}" (lines 1815-1817) -- that statement is only a REAL
linear-algebra fact (both are indefinite quadratic forms of signature (3,3)
over R), not an INTEGRAL equivalence: U and "2*eta" are not the same Z-lattice
(discriminant -1 vs -64). We verified this is consequential: TT's own worked
example (4.15) alpha_x = 2 e1 requires (e1,e1) = 2 (from A.6) to reproduce
their number alpha_x^2 = 8 (line 878-885). Under the (A.2)/(A.3) hyperbolic
U-basis, e1 would be a null vector (norm 0), and TT's own arithmetic (4.15)-
(4.16) would NOT reproduce. We therefore use the (A.6)/(A.7) diagonal basis
throughout this track, and we do NOT use the "U+U" slice suggested as an
example truncation in our instructions, because it fails to reproduce TT's
own checked numbers. Every reported count in this track states which basis
it is computed in.

Flux quantisation: TT eq. (2.5) (lines 199-208) plus the remark after
(2.6)-(2.7) (lines 217-224) requires that in the ei basis, flux vector
coefficients are EVEN integers: alpha_x = sum_i a_i e_i with a_i in 2*Z.
We therefore represent each flux vector by an integer tuple n = (n_1,...,n_r)
with a_i = 2*n_i, n_i in Z.
"""
from fractions import Fraction
import itertools

# sign convention: +1 for e1,e2,e3 (positive-norm block), -1 for e4,e5,e6
# (negative-norm block), per (A.6)-(A.7).
FULL_SIGNS = (1, 1, 1, -1, -1, -1)
FULL_LABELS = ("e1", "e2", "e3", "e4", "e5", "e6")


def truncation(indices):
    """Return (signs, labels) for a sub-basis given by 0-based indices into
    the full 6-vector basis {e1,...,e6}."""
    signs = tuple(FULL_SIGNS[i] for i in indices)
    labels = tuple(FULL_LABELS[i] for i in indices)
    return signs, labels


def norm(n, signs):
    """Exact integer norm v^2 for v = sum_i (2*n_i) e_i, i.e.
    v^2 = sum_i signs_i * (2 n_i)^2 = 8 * sum_i signs_i * n_i^2."""
    return 8 * sum(s * ni * ni for s, ni in zip(signs, n))


def dot(n, m, signs):
    """Exact integer inner product of v = sum (2 n_i) e_i and
    w = sum (2 m_i) e_i: v.w = 8 * sum_i signs_i * n_i * m_i."""
    return 8 * sum(s * ni * mi for s, ni, mi in zip(signs, n, m))


def enumerate_vectors(rank, bound):
    """All integer tuples n of length `rank` with |n_i| <= bound (inclusive),
    EXCLUDING the all-zero vector."""
    for n in itertools.product(range(-bound, bound + 1), repeat=rank):
        if any(n):
            yield n


def coordinate_automorphisms(signs):
    """Lattice automorphisms of the diagonal form that are (a) independent
    permutations of coordinates WITHIN the positive-norm block and WITHIN the
    negative-norm block, and (b) independent sign flips of each coordinate.
    Both preserve the diagonal quadratic form sum signs_i x_i^2 exactly,
    because the metric has no cross terms and equal-sign entries share the
    same norm value (+2 or -2). This is a FINITE, explicitly justified
    SUBGROUP of the full (infinite, arithmetic) automorphism group
    Aut(Gamma_{3,3}); it is NOT the full duality group, and orbit counts
    under it are upper bounds on the count of physically distinct
    configurations, not exact orbit counts under the full duality group.

    Returns a list of functions g: tuple -> tuple.
    """
    pos_idx = [i for i, s in enumerate(signs) if s == 1]
    neg_idx = [i for i, s in enumerate(signs) if s == -1]
    n = len(signs)

    group_elems = []
    for pos_perm in itertools.permutations(pos_idx):
        for neg_perm in itertools.permutations(neg_idx):
            # combined permutation array: where does output slot i's value
            # come from in the input
            perm = [None] * n
            for src, dst in zip(pos_idx, pos_perm):
                perm[dst] = src
            for src, dst in zip(neg_idx, neg_perm):
                perm[dst] = src
            for flips in itertools.product((1, -1), repeat=n):
                group_elems.append((tuple(perm), flips))

    def make_g(perm, flips):
        def g(v):
            return tuple(flips[i] * v[perm[i]] for i in range(n))
        return g

    return [make_g(perm, flips) for perm, flips in group_elems]


def group_order(rank_pos, rank_neg):
    return (__import__("math").factorial(rank_pos) * (2 ** rank_pos)) * (
        __import__("math").factorial(rank_neg) * (2 ** rank_neg)
    )
