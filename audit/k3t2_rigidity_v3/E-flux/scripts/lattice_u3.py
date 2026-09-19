"""Lattice model in TT's UNIMODULAR presentation (A.2)-(A.3): Gamma_{3,3} = U + U + U.

Coordinates: lambda = sum_i x_i f_i + y_i g_i, ordered (f1,g1,f2,g2,f3,g3), Gram = H33 (parsed from
the pinned TT text, not retyped).  Flux vectors are alpha = 2*lambda with lambda in Gamma
(TT eq. 2.5 and the App. A remark 'lattice vectors with even coefficients').

The 'group G' is the finite group of integer matrices  (Z2 x Z2) wr S3  acting on the three U blocks:
  block permutations, x_i<->y_i in a block, (x_i,y_i)->-(x_i,y_i) in a block.
"""
import itertools
import numpy as np
import sympy as sp
from common import parse_H33

H33 = parse_H33()
H33_np = np.array(H33, dtype=np.int64)

def group_G():
    """All 384 matrices M (acting on column coordinate vectors) of the wreath-product group."""
    mats = []
    for perm in itertools.permutations(range(3)):
        for swaps in itertools.product((0, 1), repeat=3):
            for signs in itertools.product((1, -1), repeat=3):
                M = np.zeros((6, 6), dtype=np.int64)
                for blk in range(3):
                    dst = perm[blk]
                    s = signs[blk]
                    if swaps[blk]:
                        M[2 * dst + 0, 2 * blk + 1] = s   # y_blk -> x_dst
                        M[2 * dst + 1, 2 * blk + 0] = s   # x_blk -> y_dst
                    else:
                        M[2 * dst + 0, 2 * blk + 0] = s
                        M[2 * dst + 1, 2 * blk + 1] = s
                mats.append(M)
    return mats

def det_on_positive_plane(M):
    """Exact determinant of M restricted to the positive 3-plane P = span(f_i+g_i)
    (M maps P to itself for every element of G; asserted)."""
    P = sp.Matrix([[1 if (j == 2 * i or j == 2 * i + 1) else 0 for j in range(6)] for i in range(3)]).T  # 6x3
    MP = sp.Matrix(M.tolist()) * P
    # solve P * R = MP
    R = (P.T * P).inv() * P.T * MP
    assert P * R == MP, "M does not preserve the positive plane"
    return R.det()

def orient_plus(M):
    return det_on_positive_plane(M) == 1

def is_isometry(M, H=None):
    H = H33_np if H is None else H
    return np.array_equal(M.T @ H @ M, H)

def in_D(v):
    """v in D (the sublattice spanned by TT's (A.6) e_i): x_i == y_i mod 2 in each block."""
    return all((v[2 * i] - v[2 * i + 1]) % 2 == 0 for i in range(3))
