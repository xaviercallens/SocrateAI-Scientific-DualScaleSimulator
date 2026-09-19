"""Shared exact-arithmetic helpers for Track D (v3). Relative paths only.

Freudenthal/Kuhn triangulation of the periodic grid Z_N^4 (N even), the
canonical-representative quotient by negation, sparse GF(p) linear algebra
and a Mayer-Vietoris bookkeeping function. All counts are Python ints.
"""
from __future__ import annotations
import itertools
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]          # repo (worktree) root: .../audit/k3t2_rigidity_v3/D-tda -> root
D = 4


def kuhn_top(N, d=D):
    tops = []
    for v in itertools.product(range(N), repeat=d):
        for pi in itertools.permutations(range(d)):
            verts = [v]
            cur = list(v)
            for k in pi:
                cur[k] = (cur[k] + 1) % N
                verts.append(tuple(cur))
            assert len(set(verts)) == d + 1, "degenerate simplex: N too small"
            tops.append(tuple(verts))
    return tops


def neg(v, N):
    return tuple((-x) % N for x in v)


def build_quotient(N):
    """Return (tops_T, tops_Q(set of frozenset), fixed_points, n_lost_to_orbit_collision)."""
    tops = kuhn_top(N)
    canon = lambda v: min(v, neg(v, N))
    q = set()
    collide = 0
    for s in tops:
        c = tuple(canon(v) for v in s)
        if len(set(c)) != len(c):
            collide += 1
            continue
        q.add(frozenset(c))
    fixed = sorted({v for s in tops for v in s if neg(v, N) == v})
    return tops, q, fixed, collide


def all_faces_by_dim(tops):
    """dict dim -> set of frozenset faces (all faces of all tops)."""
    out = {}
    for s in tops:
        s = tuple(s)
        for r in range(1, len(s) + 1):
            for c in itertools.combinations(s, r):
                out.setdefault(r - 1, set()).add(frozenset(c))
    return out


def fvector(faces):
    return [len(faces[k]) for k in sorted(faces)]


def chi_from_f(f):
    return sum((-1) ** k * c for k, c in enumerate(f))


# ---------------- GUDHI helpers ----------------
def gudhi_betti(top_faces, field):
    """top_faces: iterable of iterables of hashable vertices (facets)."""
    import gudhi
    verts = sorted({v for s in top_faces for v in s})
    idx = {v: i for i, v in enumerate(verts)}
    st = gudhi.SimplexTree()
    for s in top_faces:
        st.insert(sorted(idx[v] for v in s), filtration=0.0)
    st.compute_persistence(homology_coeff_field=field, min_persistence=0, persistence_dim_max=True)
    return list(st.betti_numbers())


def pad(b, n=5):
    b = list(b)
    return b + [0] * (n - len(b))


# ---------------- sparse GF(p) elimination ----------------
class SparseRank:
    """Incremental rank over GF(p). Columns are dict row->val (p prime)."""

    def __init__(self, p):
        self.p = p
        self.piv = {}
        self.rank = 0

    def add(self, col):
        p = self.p
        col = {r: v % p for r, v in col.items() if v % p}
        while col:
            r = max(col)
            if r in self.piv:
                pc = self.piv[r]
                f = col[r]  # pc[r] normalised to 1
                for rr, vv in pc.items():
                    nv = (col.get(rr, 0) - f * vv) % p
                    if nv:
                        col[rr] = nv
                    else:
                        col.pop(rr, None)
            else:
                inv = pow(col[r], -1, p)
                self.piv[r] = {rr: (vv * inv) % p for rr, vv in col.items()}
                self.rank += 1
                return True
        return False


def nullspace_gfp(cols, p):
    """cols: list of dict row->val. Return basis of {c: sum c_i col_i = 0} as list of dict i->val."""
    piv = {}
    null = []
    for i, col in enumerate(cols):
        col = {r: v % p for r, v in col.items() if v % p}
        comb = {i: 1}
        while col:
            r = max(col)
            if r in piv:
                pc, pcomb = piv[r]
                f = col[r]
                for rr, vv in pc.items():
                    nv = (col.get(rr, 0) - f * vv) % p
                    if nv: col[rr] = nv
                    else: col.pop(rr, None)
                for ii, vv in pcomb.items():
                    nv = (comb.get(ii, 0) - f * vv) % p
                    if nv: comb[ii] = nv
                    else: comb.pop(ii, None)
            else:
                inv = pow(col[r], -1, p)
                piv[r] = ({rr: vv * inv % p for rr, vv in col.items()},
                          {ii: vv * inv % p for ii, vv in comb.items()})
                break
        if not col:
            null.append(comb)
    return null


def boundary_col(simplex_sorted, index_of_face):
    col = {}
    for i in range(len(simplex_sorted)):
        f = simplex_sorted[:i] + simplex_sorted[i + 1:]
        col[index_of_face[f]] = (-1) ** i
    return col


# ---------------- Mayer-Vietoris ----------------
def mv_betti(bA, bB, bAB, ranks, maxdeg=4):
    """Field-coefficient Mayer-Vietoris  ... -> H_n(AB) -phi_n-> H_n(A)+H_n(B) -> H_n(X) -> H_{n-1}(AB) -> ...
    bA,bB,bAB: lists of dims (padded). ranks: dict n -> rank phi_n (must be supplied for every n with
    bAB[n] > 0). If bAB[n]==0 the rank is 0 automatically. Returns b_n(X) and the pieces."""
    out = []
    for n in range(0, maxdeg + 1):
        def rk(m):
            if m < 0:
                return 0, 0
            dom = bAB[m] if m < len(bAB) else 0
            if dom == 0:
                return dom, 0
            assert m in ranks, f"rank of phi_{m} not supplied but H_{m}(A cap B) has dim {dom}"
            return dom, ranks[m]
        codim = (bA[n] if n < len(bA) else 0) + (bB[n] if n < len(bB) else 0)
        _, rn = rk(n)
        dom1, r1 = rk(n - 1)
        out.append({"n": n, "coker_phi_n": codim - rn, "ker_phi_n-1": dom1 - r1, "b_n": codim - rn + dom1 - r1})
    return out
