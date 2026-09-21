"""
Known-answer tests against Stream 1's Lean 4 formalization.

    repo   SocrateAI-DualScaleTopologicalUniverseModel-LeanProposal  ("Stream 1")
    commit bb74acb56f386a97e433f94eb0b2632ed03bc4ca
    files  Agora/Geometry/{ModularAction,SelfDual,MnLattice,SymSquareForms}.lean
    status 0 sorry; axioms: Lean's three only.

Each test names the Lean theorem it checks.  These are NOT re-proofs -- the Lean
kernel is the judge of the theorems.  They are regression guards that this
repository's constants and code paths stay consistent with what is proved there,
and negative controls for the two claims this repository got wrong
(``audit/STREAM1_BRIDGE.md``, findings S1-F1 and S1-F2).

Scope, from Stream 1 ``ModularAction.lean:47``: "NO physics."  Nothing here is a
physical claim.
"""

import math
import warnings

import numpy as np
import pytest

warnings.filterwarnings("ignore")

from leanflow.core.gamma0n_plus import (  # noqa: E402
    fricke_involution,
    fricke_matrix,
    gram_G0N,
    gram_U_plus_2N,
    height,
    is_invariant_under,
    period,
    rho,
    rho_AL,
    self_dual_tau,
    sym2,
)

# The level this repository actually works at.  `specs/spec phase 1.md` carries an
# "EXTREMAL LEVEL-12 ETA-QUOTIENT" sub-article, and `workshopcosmo.FRICKE_Y` is
# 1/sqrt(12).  See S1-F1.
N_REPO = 12


# ---------------------------------------------------------------------------
# ModularAction.lean -- the exact integer representation of Gamma_0(N)+
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("N", [1, 2, 7, 12, 49])
def test_rhoAL_fricke_is_the_N_independent_integer_matrix(N):
    """``rhoAL_fricke (N : K) : rhoAL N 0 (-1) 1 0 = !![0,-1,0; -1,0,0; 0,0,-1]``.

    The Fricke involution is an EXACT INTEGER matrix whose entries do not
    involve N or sqrt(N): the sqrt(N) of the Atkin-Lehner element cancels in the
    symmetric square.  That is why it is a lattice isometry at all.
    """
    F = np.array([[0, -1, 0], [-1, 0, 0], [0, 0, -1]], dtype=np.int64)
    assert np.array_equal(rho_AL(N, 0, -1, 1, 0), F)
    assert np.array_equal(fricke_matrix(N), F)
    assert F.dtype == np.int64


@pytest.mark.parametrize("N", [1, 2, 7, 12, 49])
def test_fricke_is_an_isometry_and_an_involution(N):
    """``rhoAL_isometry`` at ``(Nad - bc)^2 = 1``, and ``swap_involution``."""
    F, T = fricke_matrix(N), gram_U_plus_2N(N)
    assert np.array_equal(F.T @ T @ F, T)          # rhoAL_isometry
    assert np.array_equal(F @ F, np.eye(3, dtype=np.int64))


@pytest.mark.parametrize("N", [1, 2, 7, 12])
def test_rhoAL_det_is_the_determinant_cube(N):
    """``rhoAL_det : (rhoAL N a b c d).det = (N*a*d - b*c)^3``.

    At the Fricke element this is ``+1``: Gamma_0(N)+ acts by ORIENTATION-PRESERVING
    isometries, i.e. lands in SO(2,1), not merely O(2,1).
    """
    for (a, b, c, d) in [(0, -1, 1, 0), (1, 4, -2, -1)]:
        assert round(np.linalg.det(rho_AL(N, a, b, c, d))) == (N * a * d - b * c) ** 3


@pytest.mark.parametrize("N", [2, 7, 12])
def test_rho_isometry_homomorphism_det_and_trace(N):
    """``rho_isometry``, ``rho_mul``, ``rho_det``, ``rho_trace``."""
    T = gram_U_plus_2N(N)
    gs = [(1, 0, 0, 1), (1, 1, 0, 1), (1, 0, 1, 1), (1, -1, 1, 1 - N)]
    for (a, b, c, d) in gs:
        R = rho(N, a, b, c, d)
        det_g = a * d - N * b * c
        if det_g ** 2 == 1:
            assert np.array_equal(R.T @ T @ R, T)                     # rho_isometry
        assert round(np.linalg.det(R)) == det_g ** 3                  # rho_det
        assert np.trace(R) == (a + d) ** 2 - det_g                    # rho_trace (Sym^2 character)
    # rho_mul (verbatim parametrisation):
    #   rho N a b c d * rho N a' b' c' d'
    #     = rho N (a*a' + N*b*c') (a*b' + b*d') (c*a' + d*c') (N*c*b' + d*d')
    for (a, b, c, d), (a2, b2, c2, d2) in [((1, 1, 0, 1), (1, 0, 1, 1)),
                                           ((2, 1, 1, 3), (1, -1, 2, 1))]:
        lhs = rho(N, a, b, c, d) @ rho(N, a2, b2, c2, d2)
        rhs = rho(N, a * a2 + N * b * c2, a * b2 + b * d2,
                  c * a2 + d * c2, N * c * b2 + d * d2)
        assert np.array_equal(lhs, rhs)


def test_sym2_lift_is_contravariant_not_covariant():
    """``sym2_contravariant`` with the negative control ``sym2_not_covariant``.

    Sym^2 acts by substitution into the form, Q -> Q o M, and substitution
    reverses composition.  The covariant form is FALSE.  LeanMaster v3.44.0
    (ede49f0) adopted this after re-deriving it; it is binding on any future
    Sym^2 code in this repository (S1-F5).
    """
    M = np.array([[1, 1], [0, 1]], dtype=np.int64)
    Mp = np.array([[1, 0], [1, 1]], dtype=np.int64)
    assert np.array_equal(sym2(M @ Mp), sym2(Mp) @ sym2(M))       # sym2_contravariant
    assert not np.array_equal(sym2(M @ Mp), sym2(M) @ sym2(Mp))   # sym2_not_covariant

    # sym2_det, sym2_trace, sym2_isometry_general on the form lattice G0 = diag-ish
    for X in (M, Mp, M @ Mp, np.array([[2, 1], [1, 1]], dtype=np.int64)):
        detX = int(round(np.linalg.det(X)))
        assert round(np.linalg.det(sym2(X))) == detX ** 3          # sym2_det
        assert np.trace(sym2(X)) == np.trace(X) ** 2 - detX        # sym2_trace


def test_the_two_rank3_lattices_are_not_isometric():
    """``G0N_det_ne_TN_det`` / ``no_isometry_G0N_TN`` (S1-F6).

    ``<1> + U(2N)`` has det ``-4N^2``; ``U + <2N>`` has det ``-2N``.  An isometry has
    determinant +-1 and so preserves the Gram determinant, hence none exists.
    """
    for N in (1, 2, 7, 12):
        T, G = gram_U_plus_2N(N), gram_G0N(N)
        assert round(np.linalg.det(T)) == -2 * N          # TN_det
        assert round(np.linalg.det(G)) == -4 * N * N      # G0N_det
        if N != 1:  # at N = 1 the determinants coincide and the argument says nothing
            assert round(np.linalg.det(T)) != round(np.linalg.det(G))


# ---------------------------------------------------------------------------
# SelfDual.lean -- the self-dual locus and the dual-scale height
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("N", [1, 2, 7, 12, 49])
def test_self_dual_locus_is_N_tau_squared_eq_minus_one(N):
    """``root_orthogonal_iff_selfdual : bN N rootK (period N tau) = 0 <-> N tau^2 = -1``."""
    tau = self_dual_tau(N)
    assert abs(N * tau ** 2 + 1) < 1e-12
    # the (-2)-root e - f, paired with the period through the form of U + <2N>
    w = period(N, tau)
    pairing = w[0] * (-1) + w[1] * 1          # bN N ![1,-1,0] omega  (= -(N tau^2 + 1))
    assert abs(pairing) < 1e-12
    # and it is a genuine (-2)-root
    r = np.array([1, -1, 0], dtype=np.int64)
    assert r @ gram_U_plus_2N(N) @ r == -2    # root_norm


@pytest.mark.parametrize("N", [1, 2, 7, 12, 49])
def test_fricke_involution_fixes_the_self_dual_point_and_S_does_not(N):
    """``rhoAL_mulVec_period_field`` at ``(0,-1,1,0)``: ``tau -> -1/(N tau)``.

    Negative control: ``S : tau -> -1/tau`` fixes it only at ``N = 1``.  This is
    exactly the error in S1-F1/S1-F2 -- the code's level-12 point is not fixed by
    the SL(2,Z) generator the fold uses.
    """
    tau = self_dual_tau(N)
    assert abs(fricke_involution(tau, N) - tau) < 1e-14
    S = -1.0 / tau
    if N == 1:
        assert abs(S - tau) < 1e-14
    else:
        assert abs(S - tau) > 1e-6


def test_height_bound_fricke_invariance_and_equality_locus():
    """``height_ge_two``, ``height_fricke``, ``height_eq_two_iff``.

    ``height_ge_two`` is LeanMaster's circle dual-scale bound invoked directly, so
    this also guards agreement with LeanMaster ede49f0 (v3.44.0).
    """
    import random
    rng = random.Random(2)
    N = N_REPO
    worst = 0.0
    for _ in range(20000):
        t = 10.0 ** rng.uniform(-3, 3)
        assert height(N, t) >= 2.0 - 1e-12                       # height_ge_two
        worst = max(worst, abs(height(N, 1.0 / (N * t)) - height(N, t)))  # height_fricke
    assert worst < 1e-8
    t_sd = 1.0 / math.sqrt(N)
    assert abs(height(N, t_sd) - 2.0) < 1e-14                    # height_eq_two_iff (=>)
    assert height(N, t_sd * 1.01) > 2.0                          # height_eq_two_iff (<=)


# ---------------------------------------------------------------------------
# Against THIS repository's constants and code paths
# ---------------------------------------------------------------------------

def test_repo_fricke_constant_is_the_level_12_self_dual_point():
    """S1-F1.  ``FRICKE_Y = 1/sqrt(12)`` is ``tau = i/sqrt(12)``, i.e. level N = 12.

    If this fails, either the constant or the level has been changed and
    ``audit/STREAM1_BRIDGE.md`` plus the manuscript must be updated with it.
    """
    from workshopcosmo import FRICKE_X, FRICKE_Y
    assert FRICKE_X == 0.0
    assert abs(1.0 / FRICKE_Y ** 2 - N_REPO) < 1e-12
    assert abs(complex(FRICKE_X, FRICKE_Y) - self_dual_tau(N_REPO)) < 1e-15


def test_potential_is_stationary_at_the_level_12_point_and_not_at_tau_i():
    """S1-F1, the evidence that the CODE is right and the manuscript was wrong."""
    import workshopcosmo as wc
    _, dx, dy = wc.compute_potential(0.0, 1.0 / math.sqrt(N_REPO))
    assert abs(dx) < 1e-12 and abs(dy) < 1e-12          # stationary at i/sqrt(12)
    _, dx_i, dy_i = wc.compute_potential(0.0, 1.0)
    assert abs(dy_i) > 1.0                               # NOT stationary at tau = i
    _, dxo, dyo = wc.compute_potential(0.5, math.sqrt(3) / 2)
    assert abs(dxo) < 1e-12 and abs(dyo) < 1e-12        # stationary at the orbifold point


def test_potential_is_T_invariant_but_not_S_or_Fricke_invariant():
    """S1-F2, the negative control.

    Folding is only valid under a genuine symmetry.  ``compute_potential`` is
    invariant under ``T : tau -> tau + 1`` and under NEITHER ``S : tau -> -1/tau``
    (which ``modular_domain_fold`` applies) nor Fricke ``tau -> -1/(N tau)``.
    """
    import workshopcosmo as wc
    V = lambda x, y: wc.compute_potential(x, y)[0]

    ok_T, dev_T = is_invariant_under(V, lambda z: z + 1, samples=500, atol=1e-9)
    assert ok_T, f"T should be a symmetry, deviation {dev_T}"

    ok_S, dev_S = is_invariant_under(V, lambda z: -1.0 / z, samples=500)
    assert not ok_S and dev_S > 1.0

    ok_F, dev_F = is_invariant_under(V, lambda z: -1.0 / (N_REPO * z), samples=500)
    assert not ok_F and dev_F > 1.0


def test_sl2z_fold_moves_the_repo_own_fricke_saddle():
    """S1-F2.  The SL(2,Z) fold teleports the model's own critical point.

    ``|i/sqrt(12)| = 0.2887 < 1`` so the point is outside the SL(2,Z) fundamental
    domain and ``S`` maps it to ``i*sqrt(12)``.
    """
    from leanflow.core.projections import modular_domain_fold
    from workshopcosmo import FRICKE_Y

    # With the S-step enabled -- the behaviour before the S1-F2 fix, kept as the
    # regression guard for why it is now opt-in.
    x_f, y_f, folds = modular_domain_fold(0.0, FRICKE_Y, apply_S=True)
    assert folds >= 1
    assert abs(y_f - math.sqrt(N_REPO)) < 1e-9          # -> i*sqrt(12), not fixed
    assert abs(y_f - FRICKE_Y) > 1.0

    # Default is now T-only, which leaves the saddle where it is.
    x_d, y_d, folds_d = modular_domain_fold(0.0, FRICKE_Y)
    assert folds_d == 0
    assert abs(y_d - FRICKE_Y) < 1e-15


# ---------------------------------------------------------------------------
# S1-F8 -- the solver's state-layout contract
# ---------------------------------------------------------------------------

def test_solver_skips_projections_on_an_undeclared_state_layout():
    """S1-F8.  The projections used to assume (Re tau, Im tau) = indices (0, 1).

    For ``workshopcosmo.cosmology_rhs`` the state is ``[a, x, y, u, v]``: index 0 is
    the SCALE FACTOR and Im tau is at index 2.  Guessing (0, 1) there folded the
    scale factor as if it were Re tau and clamped Re tau -- which is legitimately
    0 at the Fricke point and 0.5 at the orbifold point -- to 1/sqrt(12).
    An undeclared layout must now skip, not guess.
    """
    import numpy as np
    import workshopcosmo as wc
    from leanflow.core.solver import Solver

    y0 = np.array([1e-10, 0.001, wc.FRICKE_Y + 0.001, 0.0, 0.0])
    rhs = lambda t, y: wc.cosmology_rhs(t, y)

    with pytest.warns(RuntimeWarning, match="modulus_indices"):
        r = Solver(method="Radau").solve(rhs, y0, (0.0, 50.0))
    assert r.telemetry.modular_folds_count == 0
    assert r.telemetry.projections_applied == 0


def test_solver_projects_the_right_components_when_told():
    """S1-F8.  With ``modulus_indices=(1, 2)`` the cosmology state projects correctly."""
    import numpy as np
    import workshopcosmo as wc
    from leanflow.core.solver import Solver

    y0 = np.array([1e-10, 0.001, wc.FRICKE_Y + 0.001, 0.0, 0.0])
    rhs = lambda t, y: wc.cosmology_rhs(t, y)
    r = Solver(method="Radau", modulus_indices=(1, 2)).solve(rhs, y0, (0.0, 50.0))

    # Re(tau) is T-folded into [-1/2, 1/2]; Im(tau) respects the floor; a is untouched.
    assert r.y[1].min() >= -0.5 - 1e-9 and r.y[1].max() <= 0.5 + 1e-9
    assert r.y[2].min() >= wc.FRICKE_Y - 1e-9
    assert r.y[0].max() > 1.0                       # the scale factor grew, unprojected
    assert r.telemetry.modular_folds_count > 0


def test_two_and_four_component_layouts_keep_the_old_default():
    """S1-F8.  The layouts the projection docstrings name are unchanged: (0, 1)."""
    import numpy as np
    from leanflow.core.solver import Solver

    rhs4 = lambda t, y: np.array([y[2], y[3], 0.0, 0.0])
    r = Solver(method="RK45").solve(rhs4, np.array([3.2, 1.5, 0.0, 0.0]), (0.0, 1.0))
    assert abs(r.y[0, -1] - 0.2) < 1e-9             # T-folded 3.2 -> 0.2
    assert r.telemetry.modular_folds_count > 0
