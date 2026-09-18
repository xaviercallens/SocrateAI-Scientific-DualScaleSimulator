"""
Shared exact-arithmetic helpers for the K3 x T^2 rigidity track.

Builds a Freudenthal (Kuhn) triangulation of the periodic grid Z_N^d,
i.e. a simplicial triangulation of the d-torus (R/N Z)^d, and provides
group-action utilities (negation mod N, translation mod N) used to
build quotient complexes as abstract simplicial complexes for GUDHI's
SimplexTree.

All vertex/simplex bookkeeping is done with Python ints and tuples
(exact, no floats). Floats only ever appear in GUDHI's persistence
computation (filtration values, which are all 0 for these simplicial
complexes) and in eigenvalue-sign checks elsewhere.
"""
from __future__ import annotations
import itertools
from fractions import Fraction


def freudenthal_top_simplices(N: int, d: int):
    """
    Top-dimensional (d-simplex) list for the Freudenthal/Kuhn triangulation
    of the periodic grid Z_N^d (i.e. the d-torus (R/NZ)^d).

    For each base vertex v in {0,...,N-1}^d and each permutation pi of
    {0,...,d-1}, the simplex has vertices
        v,
        v + e_{pi(0)},
        v + e_{pi(0)} + e_{pi(1)},
        ...,
        v + e_{pi(0)} + ... + e_{pi(d-1)}
    with all coordinate arithmetic taken mod N (periodic identification).

    Returns a list of frozensets... no: returns a list of tuples of length d+1
    of vertex-tuples (each vertex is a tuple of d ints in [0,N)).
    Degenerate simplices (repeated vertex, which happens if N is too small)
    are dropped and counted; caller should check none were dropped.
    """
    assert N >= 2
    simplices = []
    degenerate = 0
    base_vertices = list(itertools.product(range(N), repeat=d))
    perms = list(itertools.permutations(range(d)))
    for v in base_vertices:
        for pi in perms:
            verts = [v]
            cur = list(v)
            for k in pi:
                cur[k] = (cur[k] + 1) % N
                verts.append(tuple(cur))
            if len(set(verts)) != d + 1:
                degenerate += 1
                continue
            simplices.append(tuple(verts))
    return simplices, degenerate


def all_faces(top_simplices):
    """All faces (of every dimension >=0) of a list of top simplices, as a set
    of frozensets of vertices."""
    faces = set()
    for simp in top_simplices:
        n = len(simp)
        for r in range(1, n + 1):
            for combo in itertools.combinations(simp, r):
                faces.add(frozenset(combo))
    return faces


def euler_characteristic_from_faces(faces):
    """chi = sum_{k>=0} (-1)^k * (#k-simplices), k = dim = |face|-1."""
    from collections import Counter
    counts = Counter(len(f) - 1 for f in faces)
    chi = 0
    for dim, cnt in counts.items():
        chi += ((-1) ** dim) * cnt
    return chi, dict(sorted(counts.items()))


def negate_mod(v, N):
    return tuple((-x) % N for x in v)


def translate_mod(v, N, shift):
    return tuple((x + s) % N for x, s in zip(v, shift))


def apply_map_to_simplex(simp, vmap):
    return tuple(vmap(v) for v in simp)


def build_simplex_tree(top_simplices, SimplexTree):
    """Insert all top simplices (and hence all faces) into a fresh SimplexTree
    with filtration value 0 everywhere."""
    st = SimplexTree()
    for simp in top_simplices:
        # gudhi wants vertices as ints; encode each d-tuple vertex as a single
        # int via a bijection so the SimplexTree (which expects hashable ints)
        # can use it as a vertex label. We build the encoding once by caller.
        st.insert(list(simp), filtration=0.0)
    return st


def make_vertex_encoder(N, d):
    """Bijection from {0,...,N-1}^d to {0,...,N^d - 1} (mixed-radix), and back."""
    def encode(v):
        idx = 0
        for x in v:
            idx = idx * N + x
        return idx

    def decode(idx):
        coords = []
        for _ in range(d):
            coords.append(idx % N)
            idx //= N
        return tuple(reversed(coords))

    return encode, decode
