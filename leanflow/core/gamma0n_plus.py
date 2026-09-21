"""
=============================================================================
Gamma_0(N)+ : the exact integer representation, transcribed from Stream 1's
              Lean 4 formalization
=============================================================================

Source of every identity below (quoted in ``audit/STREAM1_BRIDGE.md``):

    repo   SocrateAI-DualScaleTopologicalUniverseModel-LeanProposal  ("Stream 1")
    commit bb74acb56f386a97e433f94eb0b2632ed03bc4ca
    files  Agora/Geometry/ModularAction.lean, Agora/Geometry/SelfDual.lean,
           Agora/Geometry/MnLattice.lean
    status 0 sorry; axioms: Lean's three only.

This module is a TRANSCRIPTION, not a re-derivation. It exists so that the
kernel-proved identities are executable here and guarded by regression tests
(``tests/test_stream1_bridge.py``). Every function names the theorem it mirrors.

SCOPE, stated because it is easy to overstate (Stream 1, ModularAction.lean:47):

    "NO physics. That `Gamma_0(N)+ = O+(U+<2N>)/+-1` is Dolgachev's theorem
     (literature, not proved here): we prove the inclusion `>=` constructively,
     which is the direction the applications need."

Nothing here is a physical claim, and nothing here licenses folding a
trajectory: a fold is only valid for a function that is actually invariant
under the group. Use :func:`is_invariant_under` to check before folding.
Stream 1's README states program-wide that the Sym^2 relation supplies no
physical coupling and that the Fricke/Narain matrix coincidence is a fact about
a lattice isometry, not a physical identification.
=============================================================================
"""

from __future__ import annotations

import math
from typing import Callable

import numpy as np

__all__ = [
    "gram_U_plus_2N",
    "gram_G0N",
    "sym2",
    "rho",
    "rho_AL",
    "fricke_matrix",
    "period",
    "fricke_involution",
    "self_dual_tau",
    "height",
    "is_invariant_under",
]


def gram_U_plus_2N(N: int) -> np.ndarray:
    """Gram matrix of ``U + <2N>`` in the basis ``(e, f, w)``.

    Mirrors ``MnLattice.TN`` / ``ModularAction.TNR``::

        TNR (N : K) : Matrix (Fin 3) (Fin 3) K := !![0, 1, 0; 1, 0, 0; 0, 0, 2 * N]

    Determinant ``-2N``.  NOTE (Stream 1 ``no_isometry_G0N_TN``): this is NOT the
    lattice ``<1> + U(2N)`` of determinant ``-4N^2``.  The two are both rank 3 and
    both of signature (2,1), and they are NOT isometric -- an isometry has
    determinant +-1 and so preserves the Gram determinant.  Do not conflate them.
    """
    return np.array([[0, 1, 0], [1, 0, 0], [0, 0, 2 * N]], dtype=np.int64)


def gram_G0N(N: int) -> np.ndarray:
    """Gram matrix of the discriminant lattice of ``Gamma_0(N)``-forms, ``b^2 - 4Nac``.

    Mirrors ``SymSquareForms.G0N``::

        def G0N (N : Z) : Matrix (Fin 3) (Fin 3) Z := !![0, 0, -(2 * N); 0, 1, 0; -(2 * N), 0, 0]

    Determinant ``-4N^2`` (``G0N_det``).  This is the OTHER rank-3 lattice: see
    :func:`gram_U_plus_2N` and ``no_isometry_G0N_TN``.
    """
    return np.array([[0, 0, -(2 * N)], [0, 1, 0], [-(2 * N), 0, 0]], dtype=np.int64)


def sym2(M: np.ndarray) -> np.ndarray:
    """The ``Sym^2`` action on the coefficient triple ``(a, b, c)`` of ``a x^2 + b xy + c y^2``.

    Mirrors ``SymSquareForms.sym2`` exactly -- induced by the substitution
    ``(x, y) -> (alpha x + beta y, gamma x + delta y)`` for ``M = !![alpha, beta; gamma, delta]``::

        !![M 0 0 ^ 2,           M 0 0 * M 1 0,                 M 1 0 ^ 2;
           2 * (M 0 0 * M 0 1), M 0 0 * M 1 1 + M 0 1 * M 1 0, 2 * (M 1 0 * M 1 1);
           M 0 1 ^ 2,           M 0 1 * M 1 1,                 M 1 1 ^ 2]

    CONTRAVARIANT (``sym2_contravariant``): ``sym2(M M') = sym2(M') sym2(M)``.
    The covariant form is false, with the negative control ``sym2_not_covariant``
    on ``M = !![1,1;0,1]``, ``M' = !![1,0;1,1]``.
    """
    a, b = M[0, 0], M[0, 1]
    c, d = M[1, 0], M[1, 1]
    return np.array(
        [
            [a * a, a * c, c * c],
            [2 * (a * b), a * d + b * c, 2 * (c * d)],
            [b * b, b * d, d * d],
        ],
        dtype=M.dtype,
    )


def rho(N: int, a: int, b: int, c: int, d: int) -> np.ndarray:
    """``rho(g)`` for ``g = [[a, b], [Nc, d]]`` in ``Gamma_0(N)``.

    Mirrors ``ModularAction.rho``.  It is an isometry of ``U + <2N>`` as soon as
    ``(ad - Nbc)^2 = 1`` (``rho_isometry``), a homomorphism (``rho_mul``), has
    ``det rho = (ad - Nbc)^3`` (``rho_det``) and character
    ``tr rho = (a+d)^2 - (ad - Nbc)`` (``rho_trace``) -- the character of ``Sym^2``
    of the standard representation.

    WARNING (``sym2_contravariant``): the ``Sym^2`` LIFT of matrices is an
    ANTI-homomorphism, ``sym2(M M') = sym2(M') sym2(M)``.  ``rho`` as defined here
    is the group homomorphism on ``Gamma_0(N)``; do not compose raw ``Sym^2`` lifts
    as if they were covariant.
    """
    return np.array(
        [
            [d * d, -(N * c * c), 2 * N * c * d],
            [-(N * b * b), a * a, -(2 * N * a * b)],
            [b * d, -(a * c), a * d + N * b * c],
        ],
        dtype=np.int64,
    )


def rho_AL(N: int, a: int, b: int, c: int, d: int) -> np.ndarray:
    """``rho(W)`` for the Atkin-Lehner element ``W = [[Na, b], [Nc, Nd]]/sqrt(N)``.

    Mirrors ``ModularAction.rhoAL``::

        !![N * d ^ 2, -(c ^ 2), 2 * N * c * d;
           -(b ^ 2), N * a ^ 2, -(2 * N * a * b);
           b * d, -(a * c), N * a * d + b * c]

    Stream 1's docstring: "The `sqrt N` cancels in the symmetric square, so this
    is again an INTEGER matrix -- which is why Atkin-Lehner involutions are
    lattice isometries at all."  Isometry when ``(Nad - bc)^2 = 1``
    (``rhoAL_isometry``); ``det = (Nad - bc)^3`` (``rhoAL_det``), so the whole of
    ``Gamma_0(N)+``, Fricke included, acts by ORIENTATION-PRESERVING isometries.
    """
    return np.array(
        [
            [N * d * d, -(c * c), 2 * N * c * d],
            [-(b * b), N * a * a, -(2 * N * a * b)],
            [b * d, -(a * c), N * a * d + b * c],
        ],
        dtype=np.int64,
    )


def fricke_matrix(N: int | None = None) -> np.ndarray:
    """The Fricke involution as an exact integer matrix -- independent of ``N``.

    Mirrors ``ModularAction.rhoAL_fricke``::

        theorem rhoAL_fricke (N : K) : rhoAL N 0 (-1) 1 0 = !![0, -1, 0; -1, 0, 0; 0, 0, -1]

    and ``rhoAL_fricke_eq_neg_swap (N : Z) : rhoAL N 0 (-1) 1 0 = -swap``.

    This is the content of the observation that the Fricke involution "remains an
    exact integer matrix isometry": the entries carry no ``N`` and no ``sqrt(N)``.
    ``N`` is accepted and ignored so that call sites can pass it for clarity; when
    given it is checked against :func:`rho_AL`.
    """
    F = np.array([[0, -1, 0], [-1, 0, 0], [0, 0, -1]], dtype=np.int64)
    if N is not None:
        assert np.array_equal(rho_AL(N, 0, -1, 1, 0), F), "rhoAL_fricke violated"
    return F


def period(N: int, tau: complex) -> np.ndarray:
    """``omega(tau) = e - N tau^2 f + tau w``, in coordinates ``(e, f, w)``.

    Mirrors ``MnLattice.period`` / ``ModularAction.periodR``.  Isotropic for every
    ``tau`` (``period_isotropic``).  Source recorded by Stream 1 as Dolgachev
    (1996) section 7, the tube-domain realization of the period domain of
    ``U + <2N>`` as the upper half-plane.
    """
    return np.array([1.0, -N * tau * tau, tau], dtype=np.complex128)


def fricke_involution(tau: complex, N: int) -> complex:
    """``tau -> -1/(N tau)``.

    Mirrors ``ModularAction.rhoAL_mulVec_period_field`` specialized at
    ``(a,b,c,d) = (0,-1,1,0)``: "Over a field: `tau -> (Na tau + b)/(N(c tau + d))`.
    At `(a,b,c,d) = (0,-1,1,0)` this is `tau -> -1/(N tau)`, the Fricke involution".

    This is the involution that fixes the level-``N`` self-dual point.  It is NOT
    ``S : tau -> -1/tau``, which is the ``N = 1`` case.
    """
    return -1.0 / (N * tau)


def self_dual_tau(N: int) -> complex:
    """The self-dual point ``tau = i/sqrt(N)``, the fixed point of Fricke at level ``N``.

    Mirrors ``SelfDual.root_orthogonal_iff_selfdual``::

        theorem root_orthogonal_iff_selfdual (N tau : K) :
            bN N rootK (period N tau) = 0 <-> N * tau ^ 2 = -1

    i.e. the self-dual locus is exactly where the ``(-2)``-root ``e - f`` becomes
    orthogonal to the period.  ``N tau^2 = -1`` gives ``tau = i/sqrt(N)``.
    """
    return complex(0.0, 1.0 / math.sqrt(N))


def height(N: float, t: float) -> float:
    """``H_N(t) = N t^2 + 1/(N t^2)`` on the imaginary axis ``tau = it``.

    Mirrors ``SelfDual.height``, with three proved properties:

    * ``height_ge_two``   -- ``H_N >= 2``; Stream 1 obtains this by invoking
      LeanMaster's ``circle_effective_scale_ge_two`` at ``R = N t^2``, so this IS
      the programme's dual-scale bound, not an analogue of it.
    * ``height_fricke``   -- ``H_N(1/(Nt)) = H_N(t)``.
    * ``height_eq_two_iff`` -- equality holds exactly on ``N t^2 = 1``.
    """
    x = N * t * t
    return x + 1.0 / x


def is_invariant_under(
    f: Callable[[float, float], float],
    transform: Callable[[complex], complex],
    samples: int = 2000,
    atol: float = 1e-9,
    seed: int = 0,
    x_range: tuple[float, float] = (-0.5, 0.5),
    log10_y_range: tuple[float, float] = (-0.6, 0.8),
) -> tuple[bool, float]:
    """Check numerically whether ``f(x, y)`` is invariant under ``transform`` on ``tau``.

    Returns ``(invariant, max_abs_deviation)``.

    Folding a trajectory back into a fundamental domain is only meaningful for a
    function that is genuinely invariant under the group being used.  This repo's
    ``compute_potential`` is invariant under ``T : tau -> tau + 1`` but NOT under
    ``S : tau -> -1/tau`` and NOT under Fricke ``tau -> -1/(N tau)``; see
    ``audit/STREAM1_BRIDGE.md`` finding S1-F2.  Call this before folding.
    """
    import random

    rng = random.Random(seed)
    worst = 0.0
    for _ in range(samples):
        x = rng.uniform(*x_range)
        y = 10.0 ** rng.uniform(*log10_y_range)
        tau_t = transform(complex(x, y))
        worst = max(worst, abs(f(tau_t.real, tau_t.imag) - f(x, y)))
    return worst <= atol, worst
