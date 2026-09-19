"""
lib_cells.py -- exact cellular chain complexes and homology over Z/p.

Written for audit/tda_validation/k3t2 (branch loop/tda-k3t2). Independent of
audit/k3t2_rigidity_v2/D-tda (no code copied from lib_freudenthal.py: Part A here
uses CUBICAL CW complexes, not Freudenthal simplicial ones).

Contents
  * CC: a finite chain complex of a CW complex, cells with dims and INTEGER
    boundary coefficients (exact Python ints).
  * cubical(...): cubical CW complexes on Z_N (periodic) or integer boxes.
      cell (v, S): prod_{i in S} [v_i, v_i+1] x prod_{i notin S} {v_i}
      d(v,S) = sum_{j in S} (-1)^{pos(j,S)} [ (v+e_j, S\\j) - (v, S\\j) ]
  * quotient_by_involution: cellular chains of X/G for G = Z/2 acting
    cellularly, freely on every cell except cells fixed POINTWISE (checked).
    C(X/G) = coinvariants C(X)_G, one generator per orbit.
  * product: cellular chains of X x Y,
      d(a x b) = da x b + (-1)^{|a|} a x db.
  * glue_cylinders: algebraic mapping cylinder of a chain map phi: C(L) -> C(B)
    glued along a subcomplex L of X:  d(c x I) = (dc) x I + (-1)^{|c|} (phi(c) - c).
  * homology_mod_p (engine 1): pair elimination (algebraic reduction of the
    chain complex, Markowitz-type pivot order), all dimensions at once.
  * homology_mod_p_colred (engine 2): classical per-dimension column reduction
    with lowest-row pivots and clearing. Independent code path, used to
    cross-check engine 1.
  * dense_rank_mod_p: dense Gaussian elimination mod p (small cross-checks).

Nothing here knows any expected Betti number.
"""
from __future__ import annotations

import heapq
import itertools
import resource
import time


# ----------------------------------------------------------------------------
# chain complex container
# ----------------------------------------------------------------------------
class CC:
    """Finite cellular chain complex with integer boundary coefficients."""

    def __init__(self):
        self.dims = []      # dim of cell i
        self.bd = []        # dict face_id -> int coefficient (nonzero)
        self.label = []     # optional geometric label

    def add(self, dim, bd, label=None):
        self.dims.append(dim)
        self.bd.append({k: v for k, v in bd.items() if v != 0})
        self.label.append(label)
        return len(self.dims) - 1

    def __len__(self):
        return len(self.dims)

    def fvector(self):
        top = max(self.dims) if self.dims else -1
        f = [0] * (top + 1)
        for d in self.dims:
            f[d] += 1
        return f

    def chi(self):
        return sum((-1) ** k * x for k, x in enumerate(self.fvector()))

    def check_d2(self):
        """Exact integer check that d o d = 0. Returns list of violations (cell, face, coeff)."""
        bad = []
        for c, b in enumerate(self.bd):
            acc = {}
            for f, x in b.items():
                if self.dims[f] != self.dims[c] - 1:
                    bad.append((c, f, "dim"))
                for g, y in self.bd[f].items():
                    acc[g] = acc.get(g, 0) + x * y
            for g, v in acc.items():
                if v != 0:
                    bad.append((c, g, v))
        return bad


# ----------------------------------------------------------------------------
# cubical complexes
# ----------------------------------------------------------------------------
def cubical(sizes, periodic, keep=None):
    """
    Cubical CW complex.
    sizes[i]: if periodic[i], coordinate i lives in Z_{sizes[i]} (circle with
              sizes[i] vertices and sizes[i] edges); else vertices 0..sizes[i]
              (an interval subdivided into sizes[i] unit edges).
    keep(v, S) -> bool: optional predicate selecting a SUBCOMPLEX (caller must
              make it closed under faces; checked).
    Returns (cc, index) where index maps (v, S) -> cell id, cc.label[i] = (v, S).
    """
    n = len(sizes)
    cells = []
    for S_bits in range(1 << n):
        S = tuple(i for i in range(n) if S_bits >> i & 1)
        ranges = []
        for i in range(n):
            if periodic[i]:
                ranges.append(range(sizes[i]))
            else:
                ranges.append(range(sizes[i]) if i in S else range(sizes[i] + 1))
        for v in itertools.product(*ranges):
            if keep is None or keep(v, S):
                cells.append((len(S), v, S))
    cells.sort(key=lambda t: (t[0], t[1], t[2]))
    index = {(v, S): i for i, (_, v, S) in enumerate(cells)}
    cc = CC()
    for dim, v, S in cells:
        bd = {}
        for pos, j in enumerate(S):
            Sm = tuple(i for i in S if i != j)
            up = list(v)
            up[j] = up[j] + 1
            if periodic[j]:
                up[j] %= sizes[j]
            sgn = -1 if pos % 2 else 1
            for w, s in ((tuple(up), sgn), (v, -sgn)):
                key = (w, Sm)
                if key not in index:
                    raise ValueError(f"keep() not closed under faces: {key} missing (face of {(v, S)})")
                fid = index[key]
                bd[fid] = bd.get(fid, 0) + s
        cc.add(dim, bd, (v, S))
    return cc, index


def cubical_torus(N, n):
    return cubical([N] * n, [True] * n)


def cubical_sphere(m):
    """S^{m-1} = boundary of [-1,1]^m, unit-subdivided (coords shifted by +1: 0..2)."""
    def keep(v, S):
        return any((i not in S) and (v[i] in (0, 2)) for i in range(m))
    return cubical([2] * m, [False] * m, keep)


# ----------------------------------------------------------------------------
# involution quotients
# ----------------------------------------------------------------------------
def quotient_by_involution(cc, index, act):
    """
    act(v, S) -> ((v', S), sign): the cellular involution on labelled cells and the
    orientation sign with which it maps cell (v,S) onto cell (v',S).
    Checks: act is an involution; a cell mapped to itself must have sign +1
    (fixed POINTWISE-compatible orientation) -- here only vertices can be fixed.
    Returns (qcc, info). qcc.label[o] = representative (v,S) of orbit o.
    """
    n = len(cc)
    img = [None] * n
    sgn = [0] * n
    fixed = []
    for i, lab in enumerate(cc.label):
        (w, S2), s = act(*lab)
        j = index[(w, S2)]
        img[i], sgn[i] = j, s
    for i in range(n):
        if img[img[i]] != i or sgn[i] * sgn[img[i]] != 1:
            raise ValueError(f"not an involution at cell {i}")
        if img[i] == i:
            fixed.append(i)
            if sgn[i] != 1 or cc.dims[i] != 0:
                raise ValueError(f"cell {cc.label[i]} fixed setwise with dim {cc.dims[i]} sign {sgn[i]}")
    # orbits: representative = smaller id
    orbit_of = [None] * n
    osign = [0] * n
    reps = []
    for i in range(n):
        j = img[i]
        if j < i:
            continue
        o = len(reps)
        reps.append(i)
        orbit_of[i], osign[i] = o, 1
        if j != i:
            # act(rep) = sgn * (cell j)  =>  cell j = sgn * act(rep)  ->  sgn * [orbit]
            orbit_of[j], osign[j] = o, sgn[i]
    q = CC()
    for o, r in enumerate(reps):
        bd = {}
        for f, x in cc.bd[r].items():
            fo = orbit_of[f]
            bd[fo] = bd.get(fo, 0) + x * osign[f]
        q.add(cc.dims[r], bd, cc.label[r])
    info = {"n_cells_cover": n, "n_orbits": len(reps), "fixed_cells": [cc.label[i] for i in fixed],
            "n_fixed_cells": len(fixed), "fixed_cell_dims": sorted({cc.dims[i] for i in fixed})}
    return q, info, orbit_of, osign


def neg_action(N):
    """v -> -v on Z_N^n: (v,S) -> (v',S) with v'_i = -v_i-1 (i in S), -v_i (else); sign (-1)^|S|."""
    def act(v, S):
        w = tuple(((-x - 1) if i in S else (-x)) % N for i, x in enumerate(v))
        return (w, S), (-1) ** len(S)
    return act


def antipodal_action_box(m):
    """x -> -x on the boundary of [-1,1]^m in shifted coords 0..2 (center 1)."""
    def act(v, S):
        w = tuple((2 - x - 1) if i in S else (2 - x) for i, x in enumerate(v))
        return (w, S), (-1) ** len(S)
    return act


def klein_action(Mx, Ny):
    """Klein bottle = T^2 (Z_{2Mx} x Z_Ny) / tau, tau(x,y) = (x+Mx, -y)."""
    def act(v, S):
        x, y = v
        w = ((x + Mx) % (2 * Mx), ((-y - 1) if 1 in S else (-y)) % Ny)
        return (w, S), (-1) ** (1 if 1 in S else 0)
    return act


# ----------------------------------------------------------------------------
# products
# ----------------------------------------------------------------------------
def product(A, B):
    """Cellular chain complex of A x B: cells (a,b), d(a x b) = da x b + (-1)^|a| a x db."""
    nB = len(B)
    P = CC()
    # id(a,b) = a*nB + b
    for a in range(len(A)):
        da = A.dims[a]
        s = -1 if da % 2 else 1
        for b in range(nB):
            bd = {}
            for f, x in A.bd[a].items():
                k = f * nB + b
                bd[k] = bd.get(k, 0) + x
            for g, y in B.bd[b].items():
                k = a * nB + g
                bd[k] = bd.get(k, 0) + s * y
            P.dims.append(da + B.dims[b])
            P.bd.append({k: v for k, v in bd.items() if v})
            P.label.append((a, b))
    return P


# ----------------------------------------------------------------------------
# subcomplexes and algebraic mapping cylinders
# ----------------------------------------------------------------------------
def subcomplex(cc, cell_ids):
    """Restrict to a set of cells (must be closed under faces; checked). Returns (sub, old->new map)."""
    ids = sorted(cell_ids)
    new = {c: i for i, c in enumerate(ids)}
    S = CC()
    for c in ids:
        bd = {}
        for f, x in cc.bd[c].items():
            if f not in new:
                raise ValueError(f"cell set not closed under faces: {f} (face of {c})")
            bd[new[f]] = x
        S.add(cc.dims[c], bd, cc.label[c])
    return S, new


def glue_cylinders(X, keep_cells, pieces):
    """
    Build the chain complex of  Y := X[keep_cells]  union_{L_i}  Cyl(phi_i).
      keep_cells : cells of X forming the subcomplex Y (closed under faces).
      pieces     : list of dicts {"L": list of cell ids of X inside keep_cells,
                                  "base": list of (dim, bd_dict_in_base_local_ids),
                                  "phi": function(cell_id_in_X, dim) -> dict base_local_id -> int}
    Cylinder cells: c x I for c in L (dim+1), plus the base cells.
      d(c x I) = sum_f [dc:f] (f x I) + (-1)^{|c|} (phi(c) - c)
    Returns (Z, info) with Z.label entries ('Y', old_id) / ('cyl', i, c) / ('base', i, j).
    """
    Y, ymap = subcomplex(X, keep_cells)
    Z = CC()
    Z.dims = list(Y.dims)
    Z.bd = [dict(b) for b in Y.bd]
    Z.label = [("Y", c) for c in sorted(keep_cells)]
    info = []
    for i, pc in enumerate(pieces):
        L = pc["L"]
        Lset = set(L)
        base_ids = []
        for j, (d, bdl) in enumerate(pc["base"]):
            bd = {base_ids[k]: v for k, v in bdl.items()}
            base_ids.append(Z.add(d, bd, ("base", i, j)))
        cyl_id = {}
        for c in sorted(L, key=lambda c: X.dims[c]):
            k = X.dims[c]
            bd = {}
            for f, x in X.bd[c].items():
                if f not in Lset:
                    raise ValueError("L not a subcomplex")
                bd[cyl_id[f]] = bd.get(cyl_id[f], 0) + x
            s = -1 if k % 2 else 1
            for bj, y in pc["phi"](c, k).items():
                bd[base_ids[bj]] = bd.get(base_ids[bj], 0) + s * y
            yc = ymap[c]
            bd[yc] = bd.get(yc, 0) - s
            cyl_id[c] = Z.add(k + 1, bd, ("cyl", i, c))
        info.append({"n_L_cells": len(L), "n_base_cells": len(base_ids)})
    return Z, info


def relative_complex(cc, drop_cells):
    """Quotient chain complex C(cc)/C(sub) where drop_cells is a subcomplex."""
    drop = set(drop_cells)
    keep = [c for c in range(len(cc)) if c not in drop]
    new = {c: i for i, c in enumerate(keep)}
    R = CC()
    for c in keep:
        R.add(cc.dims[c], {new[f]: x for f, x in cc.bd[c].items() if f in new}, cc.label[c])
    return R


# ----------------------------------------------------------------------------
# homology engine 1: pair elimination mod p
# ----------------------------------------------------------------------------
def homology_mod_p(cc, p, return_stats=False):
    """
    Betti numbers of cc over F_p by algebraic reduction of the chain complex:
    repeatedly pick a face a with a nonzero incidence [db:a] (a unit mod p),
    replace dc -> dc - ([dc:a]/[db:a]) db for the other cofaces c of a, delete a, b.
    The reduced complex is chain-homotopy equivalent; at the end all boundaries
    vanish (asserted) and b_k = number of surviving k-cells.
    Pivot order: face with fewest cofaces first (free faces => no fill-in).
    """
    t0 = time.time()
    n = len(cc)
    bd = [dict() for _ in range(n)]
    cob = [dict() for _ in range(n)]
    for c in range(n):
        for f, x in cc.bd[c].items():
            x %= p
            if x:
                bd[c][f] = x
                cob[f][c] = x
    alive = [True] * n
    heap = [(len(cob[a]), a) for a in range(n) if cob[a]]
    heapq.heapify(heap)
    n_pairs = 0
    fill = 0
    while heap:
        key, a = heapq.heappop(heap)
        if not alive[a] or not cob[a]:
            continue
        if len(cob[a]) != key:
            heapq.heappush(heap, (len(cob[a]), a))
            continue
        # choose coface b with the smallest boundary
        b = min(cob[a], key=lambda c: len(bd[c]))
        u = bd[b][a]
        inv = pow(u, p - 2, p) if p > 2 else 1
        fb = bd[b]
        for c in list(cob[a]):
            if c == b:
                continue
            f = (bd[c][a] * inv) % p
            bc = bd[c]
            for x, cx in fb.items():
                nv = (bc.get(x, 0) - f * cx) % p
                if nv:
                    if x not in bc:
                        fill += 1
                    bc[x] = nv
                    cob[x][c] = nv
                else:
                    if x in bc:
                        del bc[x]
                        del cob[x][c]
        # remove b
        for x in fb:
            del cob[x][b]
        for d in cob[b]:
            del bd[d][b]
        cob[b] = {}
        bd[b] = {}
        # remove a
        for y in bd[a]:
            del cob[y][a]
        bd[a] = {}
        cob[a] = {}
        alive[a] = alive[b] = False
        n_pairs += 1
        for x in fb:
            if alive[x] and cob[x]:
                heapq.heappush(heap, (len(cob[x]), x))
    top = max(cc.dims) if n else -1
    betti = [0] * (top + 1)
    for c in range(n):
        if alive[c]:
            assert not bd[c] and not cob[c], "reduction incomplete"
            betti[cc.dims[c]] += 1
    stats = {"engine": "pair-elimination", "p": p, "n_cells": n, "n_pairs": n_pairs,
             "fill_in_entries": fill, "seconds": round(time.time() - t0, 3),
             "maxrss_MB": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1)}
    return (betti, stats) if return_stats else betti


# ----------------------------------------------------------------------------
# homology engine 2: column reduction per dimension with clearing
# ----------------------------------------------------------------------------
def homology_mod_p_colred(cc, p, return_stats=False):
    """
    Independent engine: rank of each boundary matrix d_k over F_p by the
    standard 'lowest nonzero row' column reduction (as in persistence
    algorithms), processed from the top dimension down with clearing.
    b_k = n_k - rank d_k - rank d_{k+1}.
    """
    t0 = time.time()
    n = len(cc)
    top = max(cc.dims) if n else -1
    by_dim = [[] for _ in range(top + 1)]
    for c in range(n):
        by_dim[cc.dims[c]].append(c)
    rank = [0] * (top + 2)
    cleared = set()
    for k in range(top, 0, -1):
        pivot_col = {}  # low row -> reduced column dict
        r = 0
        for c in by_dim[k]:
            if c in cleared:
                continue
            col = {f: x % p for f, x in cc.bd[c].items() if x % p}
            while col:
                low = max(col)
                if low in pivot_col:
                    other = pivot_col[low]
                    f = (col[low] * pow(other[low], p - 2, p)) % p if p > 2 else 1
                    for x, v in other.items():
                        nv = (col.get(x, 0) - f * v) % p
                        if nv:
                            col[x] = nv
                        else:
                            col.pop(x, None)
                else:
                    pivot_col[low] = col
                    r += 1
                    break
        rank[k] = r
        cleared = set(pivot_col.keys())
    nk = [len(x) for x in by_dim]
    betti = [nk[k] - rank[k] - rank[k + 1] for k in range(top + 1)]
    stats = {"engine": "column-reduction", "p": p, "n_cells": n, "ranks_d_k": rank[1:top + 1],
             "seconds": round(time.time() - t0, 3),
             "maxrss_MB": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1)}
    return (betti, stats) if return_stats else betti


# ----------------------------------------------------------------------------
# dense helpers (small complexes only)
# ----------------------------------------------------------------------------
def dense_rank_mod_p(rows, p):
    """rows: list of lists of ints. Gaussian elimination mod p, pure Python."""
    M = [[x % p for x in r] for r in rows]
    if not M:
        return 0
    m, ncol = len(M), len(M[0])
    rk = 0
    for col in range(ncol):
        piv = next((i for i in range(rk, m) if M[i][col]), None)
        if piv is None:
            continue
        M[rk], M[piv] = M[piv], M[rk]
        inv = pow(M[rk][col], p - 2, p) if p > 2 else 1
        M[rk] = [(x * inv) % p for x in M[rk]]
        for i in range(m):
            if i != rk and M[i][col]:
                f = M[i][col]
                M[i] = [(a - f * b) % p for a, b in zip(M[i], M[rk])]
        rk += 1
        if rk == m:
            break
    return rk


def boundary_matrix_dense(cc, k):
    """Dense matrix of d_k : C_k -> C_{k-1} (rows = (k-1)-cells)."""
    rows_ids = [c for c in range(len(cc)) if cc.dims[c] == k - 1]
    col_ids = [c for c in range(len(cc)) if cc.dims[c] == k]
    ri = {c: i for i, c in enumerate(rows_ids)}
    M = [[0] * len(col_ids) for _ in rows_ids]
    for j, c in enumerate(col_ids):
        for f, x in cc.bd[c].items():
            M[ri[f]][j] = x
    return M, rows_ids, col_ids


def betti_dense(cc, p):
    top = max(cc.dims)
    nk = cc.fvector()
    rank = [0] * (top + 2)
    for k in range(1, top + 1):
        M, _, _ = boundary_matrix_dense(cc, k)
        rank[k] = dense_rank_mod_p(M, p) if M and M[0] else 0
    return [nk[k] - rank[k] - rank[k + 1] for k in range(top + 1)]


def nullspace_mod_p(rows, ncols, p):
    """Basis of {x : M x = 0 mod p} for M given as list of rows (ints). Returns list of vectors."""
    M = [[x % p for x in r] for r in rows]
    pivcols = []
    rk = 0
    m = len(M)
    for col in range(ncols):
        piv = next((i for i in range(rk, m) if M[i][col]), None)
        if piv is None:
            continue
        M[rk], M[piv] = M[piv], M[rk]
        inv = pow(M[rk][col], p - 2, p) if p > 2 else 1
        M[rk] = [(x * inv) % p for x in M[rk]]
        for i in range(m):
            if i != rk and M[i][col]:
                f = M[i][col]
                M[i] = [(a - f * b) % p for a, b in zip(M[i], M[rk])]
        pivcols.append(col)
        rk += 1
    free = [c for c in range(ncols) if c not in set(pivcols)]
    basis = []
    for fc in free:
        v = [0] * ncols
        v[fc] = 1
        for i, pc in enumerate(pivcols):
            v[pc] = (-M[i][fc]) % p
        basis.append(v)
    return basis


def rss_mb():
    return round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1)
