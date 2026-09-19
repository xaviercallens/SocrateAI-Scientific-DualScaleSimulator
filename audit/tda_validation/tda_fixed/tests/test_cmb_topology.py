"""Regression guard for defect D2 (malformed HEALPix 2-complex).

Every test here FAILS against
``audit/reverse_zero/E5-cmb-tda/cmb_tda.py::build_topology`` and PASSES
against ``tda_fixed.cmb_topology.build_topology_fixed``.  The
``test_OLD_*`` tests pin the defect down, so the suite is known to be able to
fail (CLAUDE.md rule 2).

Run:
  OMP_NUM_THREADS=1 timeout 590 prlimit --as=8589934592 -- \
    <venv>/bin/python -m pytest audit/tda_validation/tda_fixed/tests -q
"""
import importlib.util
import math
import os
import sys
from collections import Counter

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
WT = os.path.abspath(os.path.join(PKG, "..", "..", ".."))
sys.path.insert(0, os.path.dirname(PKG))

from tda_fixed import cmb_topology as ct  # noqa: E402

hp = pytest.importorskip("healpy")


@pytest.fixture(scope="module")
def original():
    spec = importlib.util.spec_from_file_location(
        "cmb_tda_original_for_topology_tests",
        os.path.join(WT, "audit/reverse_zero/E5-cmb-tda/cmb_tda.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def full_sky(nside):
    return np.ones(hp.nside2npix(nside), np.uint8)


# ------------------------------------------------------------------ D2.1
@pytest.mark.parametrize("nside", [8, 16, 32])
def test_fixed_full_sky_has_euler_characteristic_2(nside):
    """Combinatorial, computed WITHOUT gudhi so a dedup or ordering error
    cannot hide behind a homology computation."""
    unm, e, tr, info = ct.build_topology_fixed(full_sky(nside), nside, return_info=True)
    V, E, F = int(unm.size), int(e.shape[0]), int(tr.shape[0])
    assert V == hp.nside2npix(nside)
    assert info["n_tri_corners"] == 8
    assert info["n_quad_corners"] == V - 6
    assert E == 2 * V + info["n_quad_corners"]
    assert F == 2 * info["n_quad_corners"] + 8
    assert V - E + F == 2


@pytest.mark.parametrize("nside", [8, 16, 32])
def test_OLD_full_sky_euler_characteristic_is_not_2(original, nside):
    unm, e, tr = original.build_topology(full_sky(nside), nside)
    V, E, F = int(unm.size), int(e.shape[0]), int(tr.shape[0])
    assert V - E + F != 2, "the original is expected to be wrong here"
    assert V - E + F == V - 4      # measured: 49148 at nside 64, 12284 at nside 32


@pytest.mark.parametrize("nside", [8, 16, 32])
def test_fixed_full_sky_betti_is_1_0_1(nside):
    unm, e, tr = ct.build_topology_fixed(full_sky(nside), nside)
    assert ct.complex_betti(int(unm.size), e, tr) == [1, 0, 1]


@pytest.mark.parametrize("nside", [8, 16])
def test_OLD_full_sky_betti_2_is_enormous(original, nside):
    unm, e, tr = original.build_topology(full_sky(nside), nside)
    b = ct.complex_betti(int(unm.size), e, tr)
    assert b[0] == 1 and b[1] == 0
    assert b[2] == int(unm.size) - 5, "the original fills every 4-clique (hollow tetrahedra)"


# ------------------------------------------------------------------ D2.2
def test_fixed_complex_is_a_surface_every_edge_in_at_most_two_triangles():
    """The structural reason the old complex was wrong: filling all four
    triangles of a 4-clique puts three triangles on some edges."""
    nside = 8
    unm, e, tr = ct.build_topology_fixed(full_sky(nside), nside)
    c = Counter()
    for a, b, d in tr:
        c[(a, b)] += 1
        c[(a, d)] += 1
        c[(b, d)] += 1
    assert max(c.values()) == 2
    assert min(c.values()) == 2, "a closed surface: every edge in exactly two triangles"
    assert len(c) == int(e.shape[0]), "every edge is a face of the complex"


def test_OLD_complex_is_not_a_surface(original):
    nside = 8
    unm, e, tr = original.build_topology(full_sky(nside), nside)
    c = Counter()
    for a, b, d in tr:
        c[(min(a, b), max(a, b))] += 1
        c[(min(a, d), max(a, d))] += 1
        c[(min(b, d), max(b, d))] += 1
    assert max(c.values()) > 2, "the original is expected to put >2 triangles on some edge"


def test_no_four_clique_is_completely_filled():
    """The exact defect: a 4-clique with all four of its triangles present is
    a hollow tetrahedron and adds 1 to b2."""
    nside = 8
    unm, e, tr = ct.build_topology_fixed(full_sky(nside), nside)
    tris = set(map(tuple, np.sort(tr, axis=1)))
    filled = 0
    for a, b, d in list(tris)[:4000]:
        # any fourth vertex forming the other three triangles?
        cand = {x for (p, q, x) in tris if (p, q) == (a, b)} | \
               {x for (p, x, q) in tris if (p, q) == (a, d)}
        for x in cand:
            quad = tuple(sorted({a, b, d, x}))
            if len(quad) != 4:
                continue
            faces = [tuple(sorted(quad[:3])), tuple(sorted(quad[1:])),
                     tuple(sorted((quad[0], quad[1], quad[3]))),
                     tuple(sorted((quad[0], quad[2], quad[3])))]
            if all(f in tris for f in faces):
                filled += 1
    assert filled == 0


# ------------------------------------------------------------------ D2.3
def test_disk_mask_gives_a_disc(original):
    nside = 16
    vec = np.array(hp.pix2vec(nside, np.arange(hp.nside2npix(nside))))
    mask = (vec[2] >= math.cos(math.radians(30.0))).astype(np.uint8)
    unm, e, tr = ct.build_topology_fixed(mask, nside)
    assert ct.complex_betti(int(unm.size), e, tr) == [1, 0, 0]
    assert int(unm.size) - int(e.shape[0]) + int(tr.shape[0]) == 1
    unm_o, e_o, tr_o = original.build_topology(mask, nside)
    b_o = ct.complex_betti(int(unm_o.size), e_o, tr_o)
    assert b_o[2] > 0, "the original is expected to find spurious 2-cycles in a disc"


# ------------------------------------------------------------------ D2.4
def test_euler_key_is_the_true_euler_characteristic_not_b0_minus_b1(original):
    """The third clause of D2: the original returned b0 - b1 under the name
    'euler_chi'.  On the full sphere that is 1, not 2."""
    nside = 8
    rng = np.random.default_rng(99)
    temp = rng.normal(size=hp.nside2npix(nside))
    unm, e, tr = ct.build_topology_fixed(full_sky(nside), nside)
    nu = np.linspace(-4, 4, 41)
    r = ct.betti_curves_from_topology(temp, unm, e, tr, nu, sublevel=True)
    assert r["b0_minus_b1"][-1] == 1
    assert r["euler_char_true"][-1] == 2, "a 2-sphere has chi = 2"
    assert r["b2_implied"][-1] == 1
    # and the true Euler characteristic equals #V - #E + #F at every level
    assert np.array_equal(r["euler_char_true"],
                          r["n_vertices_at_nu"] - r["n_edges_at_nu"] + r["n_triangles_at_nu"])
    # the original returns b0 - b1 as its third value, mislabelled euler_chi
    _, _, chi_old = original.betti_curves_from_topology(temp, unm, e, tr, nu, sublevel=True)
    assert chi_old[-1] == 1 != r["euler_char_true"][-1]


# ------------------------------------------------------------------ D2.5
def test_filtration_code_is_unchanged_on_an_external_complex(original):
    """The F1 path: a grid Freudenthal triangulation supplied from outside.
    The topology fix must not move it."""
    n = 64
    rng = np.random.default_rng(7)
    f = rng.normal(size=(n, n))
    idx = np.arange(n * n).reshape(n, n)
    e = np.r_[np.c_[idx[:, :-1].ravel(), idx[:, 1:].ravel()],
              np.c_[idx[:-1, :].ravel(), idx[1:, :].ravel()],
              np.c_[idx[:-1, :-1].ravel(), idx[1:, 1:].ravel()]].astype(np.int64)
    tr = np.r_[np.c_[idx[:-1, :-1].ravel(), idx[:-1, 1:].ravel(), idx[1:, 1:].ravel()],
               np.c_[idx[:-1, :-1].ravel(), idx[1:, :-1].ravel(), idx[1:, 1:].ravel()]].astype(np.int64)
    temp = f.ravel()
    v = np.arange(n * n)
    nu = np.linspace(temp.min() / temp.std(), temp.max() / temp.std(), 201)
    b0_o, b1_o, _ = original.betti_curves_from_topology(temp, v, e, tr, nu, sublevel=True)
    r = ct.betti_curves_from_topology(temp, v, e, tr, nu, sublevel=True)
    assert np.array_equal(b0_o, r["b0"])
    assert np.array_equal(b1_o, r["b1"])


def test_neighbour_parity_assumption_holds(original):
    """The construction rests on 'even indices of get_all_neighbours are the
    edge-sharing neighbours'.  Check it against hp.boundaries rather than
    trusting the direction names."""
    nside = 8
    nb = np.asarray(hp.get_all_neighbours(nside, np.arange(hp.nside2npix(nside))), dtype=np.int64)
    corners = lambda q: {tuple(np.round(hp.boundaries(nside, int(q), step=1)[:, k], 9))  # noqa: E731
                         for k in range(4)}
    rng = np.random.default_rng(0)
    for p in rng.choice(hp.nside2npix(nside), size=12, replace=False):
        cp = corners(p)
        for k in range(8):
            q = nb[k, p]
            if q < 0:
                assert k in ct.CORNER_ONLY_INDICES, "a missing neighbour must be a corner one"
                continue
            shared = len(cp & corners(q))
            assert shared == (2 if k in ct.EDGE_SHARING_INDICES else 1), (p, k, shared)
