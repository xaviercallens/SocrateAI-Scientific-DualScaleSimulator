"""A genuine 2-complex of the HEALPix sphere, and lower-star Betti/Euler
curves on it (fixed replacement for cmb_tda.build_topology and
cmb_tda.betti_curves_from_topology).

DEFECT FIXED (D2, found 2026-09-19 by the known-answer suite on branch
loop/tda-simple, case F3topo in
audit/tda_validation/simple_suite/results.json):
``audit/reverse_zero/E5-cmb-tda/cmb_tda.py::build_topology`` took the
8-neighbour pixel adjacency graph and inserted EVERY 3-clique.  Each group of
four pixels meeting at a grid corner is a 4-clique; filling all four of its
triangles builds a hollow tetrahedron, i.e. a 2-sphere.  Measured at nside 64,
full sky: V = 49152, E = 196596, F = 196592, V - E + F = 49148,
Betti = (1, 0, 49147).  A triangulated S2 has Betti (1, 0, 1) and chi = 2.

Third clause of the same defect: ``betti_curves_from_topology`` returned
``chi = b0 - b1`` under the name ``euler_chi``.  For a 2-complex the Euler
characteristic is ``b0 - b1 + b2``, so that key was mislabelled (and on the
OLD complex b2 was 49147, so it was wrong by a lot).  This module returns the
true Euler characteristic under ``euler_char_true`` and keeps the old quantity
under the explicit key ``b0_minus_b1``.

CONSTRUCTION (declared in fix_expectations.json, commit 1bc7ebc, before this
file was written).  Measured properties of healpy 1.20.0 that it rests on,
verified at nside 8/16/32/64 before the rule was written:

  * ``hp.get_all_neighbours`` returns ``[SW, W, NW, N, NE, E, SE, S]``.
  * The EVEN indices (SW, NW, NE, SE) are the EDGE-sharing neighbours: a
    sampled pair shares 2 pixel corners (160/160 sampled pairs at nside 8).
    The ODD indices (W, N, E, S) are corner-only: they share 1 corner
    (159/159).  This is the opposite of the naive reading of the direction
    names -- HEALPix pixels are diamonds, so their edges face the diagonals.
  * Exactly 24 neighbour entries are ``-1`` at every nside tested, all at odd
    (corner-only) positions: the 8 points of the grid where only three pixels
    meet, each seen from its 3 member pixels.

So, for each pixel ``p`` and each ``k`` in 0..3, the four pixels meeting at
one corner of ``p`` are ``G = {p, nb[2k], nb[2k+1], nb[(2k+2) % 8]}`` with the
cycle ``p - nb[2k] - nb[2k+1] - nb[2k+2] - p``.  Groups are canonicalised and
de-duplicated (each corner is enumerated once per member pixel).  A 4-pixel
corner becomes TWO triangles split along ONE canonical diagonal (through the
smallest pixel index of the group and its opposite in the cycle); a corner
whose ``nb[2k+1]`` is ``-1`` has only three pixels and becomes ONE triangle.
No 4-clique is ever filled.

Consequence for the 1-skeleton, stated because it is a real change: the fixed
complex has the 4-adjacency edges plus ONE diagonal per corner, where the old
complex had both diagonals.  Both diagonals cannot be kept: with both present
and only two of the four triangles filled, the quad carries an unfilled
1-cycle, so no surface with the full 8-adjacency 1-skeleton exists without
adding a vertex at the pixel centre.  The effect on the b0/b1 curves is
measured in ``run_validation.py`` against a declared 15% tolerance.

Combinatorics (predicted in fix_expectations.json before being computed):
with V = npix, T3 = 8 three-pixel corners and Q = V - 6 four-pixel corners,
E = 2V + Q and F = 2Q + T3, hence V - E + F = -V + Q + T3 = 2.

Tier of every number this module produces: X (numerics).
"""
from __future__ import annotations

import numpy as np

__all__ = [
    "healpix_faces",
    "build_topology_fixed",
    "betti_curves_from_topology",
    "complex_betti",
]

NEIGHBOUR_ORDER = ["SW", "W", "NW", "N", "NE", "E", "SE", "S"]
EDGE_SHARING_INDICES = (0, 2, 4, 6)   # SW, NW, NE, SE  -- share 2 corners
CORNER_ONLY_INDICES = (1, 3, 5, 7)    # W, N, E, S      -- share 1 corner


# ------------------------------------------------------------------ faces
def healpix_faces(nside):
    """Enumerate the corners of the HEALPix grid as faces of a 2-complex.

    Returns ``(quads, tris3, info)``:
      quads  (Q, 4) int64, each row a 4-pixel corner in CYCLIC order
      tris3  (T3, 3) int64, each row a 3-pixel corner
      info   counts and the combinatorial checks

    Pure combinatorics on ``hp.get_all_neighbours``; no geometry, no gudhi.
    """
    import healpy as hp

    npix = hp.nside2npix(nside)
    nb = np.asarray(hp.get_all_neighbours(nside, np.arange(npix)), dtype=np.int64)
    assert nb.shape == (8, npix), nb.shape

    # Report where the -1s are (a property of the grid, not an assumption).
    neg_by_index = {NEIGHBOUR_ORDER[k]: int((nb[k] < 0).sum()) for k in range(8)}
    neg_at_edge_positions = int(sum((nb[k] < 0).sum() for k in EDGE_SHARING_INDICES))

    p = np.arange(npix, dtype=np.int64)
    quad_rows, tri_rows = [], []
    for k in range(4):
        a = nb[2 * k]                 # edge-sharing neighbour
        c = nb[2 * k + 1]             # corner-only neighbour (may be -1)
        b = nb[(2 * k + 2) % 8]       # the other edge-sharing neighbour
        ok4 = (a >= 0) & (b >= 0) & (c >= 0)
        ok3 = (a >= 0) & (b >= 0) & (c < 0)
        if ok4.any():
            quad_rows.append(np.c_[p[ok4], a[ok4], c[ok4], b[ok4]])
        if ok3.any():
            tri_rows.append(np.c_[p[ok3], a[ok3], b[ok3]])

    Q = np.vstack(quad_rows) if quad_rows else np.zeros((0, 4), dtype=np.int64)
    T = np.vstack(tri_rows) if tri_rows else np.zeros((0, 3), dtype=np.int64)

    # Each corner is enumerated once per member pixel: de-duplicate on the
    # sorted vertex set, keeping one representative row (which carries the
    # cyclic order needed to pick a diagonal).
    n_quad_raw, n_tri_raw = Q.shape[0], T.shape[0]
    degenerate_quads = int(np.sum([len(set(r)) != 4 for r in Q])) if Q.size else 0
    degenerate_tris = int(np.sum([len(set(r)) != 3 for r in T])) if T.size else 0

    _, iq = np.unique(np.sort(Q, axis=1), axis=0, return_index=True) if Q.size else (None, np.array([], int))
    _, it = np.unique(np.sort(T, axis=1), axis=0, return_index=True) if T.size else (None, np.array([], int))
    quads = Q[np.sort(iq)] if Q.size else Q
    tris3 = T[np.sort(it)] if T.size else T

    info = {
        "nside": int(nside),
        "npix": int(npix),
        "neighbour_order": NEIGHBOUR_ORDER,
        "n_missing_neighbour_entries": int((nb < 0).sum()),
        "missing_neighbour_entries_by_direction": neg_by_index,
        "n_missing_at_edge_sharing_positions": neg_at_edge_positions,
        "n_quad_corner_rows_before_dedup": int(n_quad_raw),
        "n_tri_corner_rows_before_dedup": int(n_tri_raw),
        "n_quad_corners": int(quads.shape[0]),
        "n_tri_corners": int(tris3.shape[0]),
        "degenerate_quad_rows": degenerate_quads,
        "degenerate_tri_rows": degenerate_tris,
        "dedup_factor_quads": (n_quad_raw / quads.shape[0]) if quads.shape[0] else None,
        "dedup_factor_tris": (n_tri_raw / tris3.shape[0]) if tris3.shape[0] else None,
        "predicted_n_quad_corners_V_minus_6": int(npix - 6),
        "predicted_n_tri_corners": 8,
    }
    return quads, tris3, info


def build_topology_fixed(mask, nside, return_info=False):
    """Fixed replacement for ``cmb_tda.build_topology``.

    Same call shape: ``(mask, nside) -> (unmasked, edges_arr, tris_arr)``, so
    it is a drop-in for ``betti_curves_from_topology``.  ``mask`` is a HEALPix
    map of 0/1 with ``mask > 0`` meaning "keep".  Indices in ``edges_arr`` and
    ``tris_arr`` are positions in ``unmasked`` (as in the original).

    With a mask, the result is the INDUCED subcomplex: every simplex all of
    whose vertices are unmasked is kept.
    """
    import healpy as hp

    npix = hp.nside2npix(nside)
    mask = np.asarray(mask)
    assert mask.size == npix, (mask.size, npix)
    unmasked = np.where(mask > 0)[0]
    idx = -np.ones(npix, dtype=np.int64)
    idx[unmasked] = np.arange(unmasked.size)

    quads, tris3, info = healpix_faces(nside)

    # --- triangles: split each quad along ONE canonical diagonal.
    # Cyclic order of a quad row is (v0, v1, v2, v3); opposite pairs are
    # (v0, v2) and (v1, v3).  The diagonal is the one through the smallest
    # pixel index in the group -- canonical, independent of which member
    # pixel enumerated the corner.
    if quads.size:
        amin = np.argmin(quads, axis=1)
        r = np.arange(quads.shape[0])
        d0 = quads[r, amin]
        d1 = quads[r, (amin + 2) % 4]
        o0 = quads[r, (amin + 1) % 4]
        o1 = quads[r, (amin + 3) % 4]
        # two triangles: (d0, d1, o0) and (d0, d1, o1)
        tri_from_quads = np.vstack([np.c_[d0, d1, o0], np.c_[d0, d1, o1]])
        diagonals = np.c_[d0, d1]
    else:
        tri_from_quads = np.zeros((0, 3), dtype=np.int64)
        diagonals = np.zeros((0, 2), dtype=np.int64)

    all_tris = np.vstack([tri_from_quads, tris3]) if tris3.size else tri_from_quads

    # --- edges: the 4-adjacency (edge-sharing) pairs + one diagonal per quad.
    nb = np.asarray(hp.get_all_neighbours(nside, np.arange(npix)), dtype=np.int64)
    p = np.arange(npix, dtype=np.int64)
    e_pairs = []
    for k in EDGE_SHARING_INDICES:
        q = nb[k]
        ok = q >= 0
        e_pairs.append(np.c_[p[ok], q[ok]])
    adj_edges = np.vstack(e_pairs)
    all_edges = np.vstack([adj_edges, diagonals]) if diagonals.size else adj_edges
    all_edges = np.unique(np.sort(all_edges, axis=1), axis=0)

    n_adj_edges_full = int(np.unique(np.sort(adj_edges, axis=1), axis=0).shape[0])

    # --- restrict to the unmasked induced subcomplex, remap to local indices
    ei = idx[all_edges]
    keep_e = (ei >= 0).all(axis=1)
    edges_arr = np.unique(ei[keep_e], axis=0).astype(np.int64) if keep_e.any() else np.zeros((0, 2), np.int64)

    ti = idx[all_tris]
    keep_t = (ti >= 0).all(axis=1)
    tris_arr = np.unique(np.sort(ti[keep_t], axis=1), axis=0).astype(np.int64) if keep_t.any() else np.zeros((0, 3), np.int64)

    V, E, F = int(unmasked.size), int(edges_arr.shape[0]), int(tris_arr.shape[0])
    info.update({
        "construction": "HEALPix quad-corner 2-complex: 4-adjacency edges + one canonical "
                        "diagonal per 4-pixel corner; 2 triangles per 4-pixel corner, "
                        "1 per 3-pixel corner; no clique filling",
        "diagonal_rule": "through the smallest pixel index of the corner group and its opposite in the cycle",
        "n_full_sky_adjacency_edges": n_adj_edges_full,
        "n_full_sky_diagonals": int(diagonals.shape[0]),
        "predicted_full_sky_edges_2V_plus_Q": int(2 * info["npix"] + info["n_quad_corners"]),
        "predicted_full_sky_faces_2Q_plus_T3": int(2 * info["n_quad_corners"] + info["n_tri_corners"]),
        "n_vertices": V, "n_edges": E, "n_triangles": F,
        "euler_char_V_minus_E_plus_F": V - E + F,
        "masked": bool(V != info["npix"]),
        "tier": "X",
    })
    if return_info:
        return unmasked, edges_arr, tris_arr, info
    return unmasked, edges_arr, tris_arr


# ------------------------------------------------------- filtration / curves
def _simplex_tree(verts_f, edges_arr, tris_arr):
    import gudhi
    st = gudhi.SimplexTree()
    n = verts_f.size
    vs = np.arange(n, dtype=np.int32)[None, :]
    try:
        st.insert_batch(vs, verts_f.astype(np.float64))
        if edges_arr.size:
            ef = np.maximum(verts_f[edges_arr[:, 0]], verts_f[edges_arr[:, 1]])
            st.insert_batch(edges_arr.T.astype(np.int32), ef.astype(np.float64))
        if tris_arr.size:
            tf = np.max(verts_f[tris_arr], axis=1)
            st.insert_batch(tris_arr.T.astype(np.int32), tf.astype(np.float64))
    except (AttributeError, TypeError):  # older gudhi: fall back to the original's loop
        for i in range(n):
            st.insert([int(i)], filtration=float(verts_f[i]))
        if edges_arr.size:
            ef = np.maximum(verts_f[edges_arr[:, 0]], verts_f[edges_arr[:, 1]])
            for i in range(edges_arr.shape[0]):
                st.insert([int(edges_arr[i, 0]), int(edges_arr[i, 1])], filtration=float(ef[i]))
        if tris_arr.size:
            tf = np.max(verts_f[tris_arr], axis=1)
            for i in range(tris_arr.shape[0]):
                st.insert([int(tris_arr[i, 0]), int(tris_arr[i, 1]), int(tris_arr[i, 2])],
                          filtration=float(tf[i]))
    return st


def betti_curves_from_topology(temp, unmasked, edges_arr, tris_arr, nu_grid, sublevel=True):
    """Lower-star Betti and TRUE Euler curves.

    Same signature as ``cmb_tda.betti_curves_from_topology``, so it also runs
    on an externally supplied complex (the F1 grid Freudenthal triangulation
    uses this path, not the HEALPix one).

    The filtration is unchanged from the original: vertices carry T/sigma
    (negated when ``sublevel=False``), an edge or triangle carries the max of
    its vertices.

    Returns a DICT, not a 3-tuple, so the Euler key cannot be misread:
      b0, b1              persistent Betti curves on nu_grid
      euler_char_true     the TRUE Euler characteristic of the filtered
                          subcomplex, #V(nu) - #E(nu) + #F(nu)
      b0_minus_b1         the quantity the ORIGINAL returned under the name
                          'euler_chi'.  It equals the Euler characteristic
                          only where b2(nu) = 0.
      b2_implied          euler_char_true - b0_minus_b1, i.e. b2(nu)
    """
    temp = np.asarray(temp, dtype=float)
    unmasked = np.asarray(unmasked)
    edges_arr = np.asarray(edges_arr, dtype=np.int64).reshape(-1, 2)
    tris_arr = np.asarray(tris_arr, dtype=np.int64).reshape(-1, 3)
    nu_grid = np.asarray(nu_grid, dtype=float)

    sigma = temp[unmasked].std()
    t = temp / sigma
    if not sublevel:
        t = -t
    verts_f = t[unmasked]

    st = _simplex_tree(verts_f, edges_arr, tris_arr)
    st.make_filtration_non_decreasing()  # a no-op for a max-filtration; kept for parity
    st.compute_persistence(persistence_dim_max=True)
    diag0 = np.array(st.persistence_intervals_in_dimension(0), dtype=float).reshape(-1, 2)
    diag1 = np.array(st.persistence_intervals_in_dimension(1), dtype=float).reshape(-1, 2)

    def curve(diag):
        if diag.size == 0:
            return np.zeros(nu_grid.size, dtype=int)
        return np.array([int(np.sum((diag[:, 0] <= nu) &
                                    ((diag[:, 1] > nu) | np.isinf(diag[:, 1]))))
                         for nu in nu_grid])

    b0 = curve(diag0)
    b1 = curve(diag1)

    # TRUE Euler characteristic of the filtered subcomplex, by counting
    # simplices with filtration <= nu (signed by dimension).
    fv = np.sort(verts_f)
    fe = np.sort(np.maximum(verts_f[edges_arr[:, 0]], verts_f[edges_arr[:, 1]])) if edges_arr.size else np.zeros(0)
    ft = np.sort(np.max(verts_f[tris_arr], axis=1)) if tris_arr.size else np.zeros(0)
    nV = np.searchsorted(fv, nu_grid, side="right")
    nE = np.searchsorted(fe, nu_grid, side="right") if fe.size else np.zeros(nu_grid.size, int)
    nF = np.searchsorted(ft, nu_grid, side="right") if ft.size else np.zeros(nu_grid.size, int)
    euler_true = nV - nE + nF

    return {
        "b0": b0,
        "b1": b1,
        "euler_char_true": euler_true,
        "b0_minus_b1": b0 - b1,
        "b2_implied": euler_true - (b0 - b1),
        "n_vertices_at_nu": nV,
        "n_edges_at_nu": nE,
        "n_triangles_at_nu": nF,
        "euler_key_note": "euler_char_true is #V-#E+#F of the filtered subcomplex. "
                          "The original cmb_tda returned b0-b1 under the key 'euler_chi'; "
                          "that quantity is here under 'b0_minus_b1'.",
    }


def complex_betti(n_vertices, edges_arr, tris_arr):
    """Betti numbers of the complex with a constant filtration (no field)."""
    verts_f = np.zeros(n_vertices, dtype=float)
    st = _simplex_tree(verts_f, np.asarray(edges_arr, np.int64).reshape(-1, 2),
                       np.asarray(tris_arr, np.int64).reshape(-1, 3))
    st.compute_persistence(persistence_dim_max=True)
    b = list(st.betti_numbers())
    while len(b) < 3:
        b.append(0)
    return b[:3]
