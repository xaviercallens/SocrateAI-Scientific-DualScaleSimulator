"""
=============================================================================
A genuinely Gamma_0(N)+ -invariant potential on the upper half-plane
=============================================================================

WHY THIS EXISTS -- spec item L1 (`specs/LEANFLOW_ARCHITECTURE.md`).

The potential this repository actually integrates, `workshopcosmo.compute_potential`,
is a hand-built double well whose two stationary points are placed by construction:

    V = a*cos^2(pi x) + b*sin^2(pi x) + cos^2(pi x)*(y-y_F)^2 + sin^2(pi x)*(y-y_O)^2

It has NO modular symmetry. Measured (audit/STREAM1_BRIDGE.md S1-F2):
T : tau -> tau+1 holds to 6.7e-16, but S : tau -> -1/tau deviates by 8.77 and the
level-12 Fricke involution tau -> -1/(12 tau) by 35.81. Every fold, projection and
"native duality" claim in this repository presupposed an invariance that was absent,
which is why the S-fold had to be disabled by default.

THE CONSTRUCTION. j is SL(2,Z)-invariant, so under the Fricke involution
W_N : tau -> -1/(N tau),

    j(tau)  ->  j(-1/(N tau)) = j(N tau)        [ -1/(N tau) = S(N tau) ]
    j(N tau) ->  j(-1/tau)     = j(tau)

so W_N SWAPS the two terms and their sum is W_N-invariant. For
g = [[a,b],[Nc,d]] in Gamma_0(N): j(g tau) = j(tau) since g is in SL(2,Z), and
N*(g tau) = g'(N tau) with g' = [[a, bN],[c, d]], det g' = ad - Nbc = 1, so
j(N g tau) = j(N tau). Hence

    F_N(tau) = j(tau) + j(N tau)   is Gamma_0(N)+ -invariant.

Verified numerically to ~1e-23 under W_N and under Gamma_0(N) generators
(`tests/test_modular_potential.py`).

STATUS AND SCOPE. This module is **opt-in and changes no default**. It does not
touch `workshopcosmo.compute_potential`, any committed telemetry, or the
manuscript. Switching the simulation over is a separate decision with a cost:
every reported trajectory would change. What it does is discharge L1's acceptance
criterion, which unblocks L2 (fold by the right group), L3 (Fricke reflection
instead of clamping) and L5 (give the level an observable consequence).

NO PHYSICS is claimed. Stream 1 records program-wide that no exact physical
observable exists anywhere in this programme. This is a modular function, not a
scalar field of our universe, and the level N is an input -- `audit/K3_SELECTION.md`
records that this repository's N = 12 is fitted, not forced.
=============================================================================
"""

from __future__ import annotations

import math
from typing import Tuple

import mpmath

__all__ = [
    "j_invariant",
    "gamma0_plus_invariant",
    "self_dual_tau",
    "modular_potential",
    "modular_potential_xy",
]

# mpmath's kleinj is the normalized J = j/1728.
_J_NORM = 1728


def j_invariant(tau: complex, dps: int = 25) -> mpmath.mpc:
    """Klein's j, normalized so that j(i) = 1728 and j(exp(2 pi i/3)) = 0."""
    with mpmath.workdps(dps):
        return _J_NORM * mpmath.kleinj(mpmath.mpc(tau))


def gamma0_plus_invariant(tau: complex, N: int, dps: int = 25) -> mpmath.mpc:
    """``F_N(tau) = j(tau) + j(N tau)`` -- invariant under all of ``Gamma_0(N)+``.

    Invariant under `Gamma_0(N)` because j is `SL(2,Z)`-invariant, and under the
    Fricke involution because `W_N` swaps the two summands. See the module
    docstring for the two-line proof.
    """
    with mpmath.workdps(dps):
        t = mpmath.mpc(tau)
        return _J_NORM * (mpmath.kleinj(t) + mpmath.kleinj(N * t))


def self_dual_tau(N: int) -> complex:
    """The level-N self-dual point ``tau = i/sqrt(N)``, fixed by ``W_N``.

    Stream 1 `root_orthogonal_iff_selfdual`: the self-dual locus is `N tau^2 = -1`.
    """
    return complex(0.0, 1.0 / math.sqrt(N))


def modular_potential(tau: complex, N: int = 12, mode: str = "ratio",
                      scale: float = 10.0, dps: int = 25) -> float:
    """A bounded, ``Gamma_0(N)+``-invariant potential, zero on the orbit of ``i/sqrt(N)``.

    Invariant because `F_N` is, and any function of an invariant is invariant.
    Two forms, both in `[0, 1]`, chosen by `mode`. Neither is a physical potential;
    both are modular functions.

    ``mode="ratio"`` (default) -- **scale-free**::

        V = |F_N(tau) - F_0| / (|F_N(tau)| + |F_0|),     F_0 = F_N(i/sqrt N)

    Measured: `V ~ r^2` around the self-dual point with fitted slope **2.000 over
    four decades** (`r = 1e-3 .. 1e-7`), and gradient `(0, -5e-9)` there, i.e. a
    genuine nondegenerate minimum. **Limitation, stated:** `|F_N|` outgrows `|F_0|`
    so fast that `V` saturates to `1.0` away from the well -- there is a clean basin
    on the `Gamma_0(N)+` orbit and an almost forceless plateau elsewhere.

    ``mode="log"`` -- **spread across the domain**::

        u = |F_N(tau) - F_0|,   w = log(1+u),   V = w / (scale + w)

    Measured: `V` runs `0 -> 0.97` across the domain, so it is not flat. **Limitation,
    stated:** the minimum is needle-like. `F_N` is so steep that the quadratic core
    extends only to `r ~ 1e-7`; outside it the well is logarithmic, and a finite
    difference at `h = 1e-5` already reports a nonzero gradient at the minimum.

    Why the log is needed at all in that mode: `j(N tau) ~ exp(2 pi N y)`, so the
    raw `u/(1+u)` saturates to `1.0` to machine precision almost everywhere -- it
    reads `1.0` at `tau = i`, at `0.3+0.9i` and at `0.1+5i` alike.
    """
    if mode not in ("ratio", "log"):
        raise ValueError(f"mode must be 'ratio' or 'log', got {mode!r}")
    with mpmath.workdps(dps):
        ref = gamma0_plus_invariant(self_dual_tau(N), N, dps=dps)
        val = gamma0_plus_invariant(tau, N, dps=dps)
        u = abs(val - ref)
        if mode == "ratio":
            return float(u / (abs(val) + abs(ref)))
        w = mpmath.log(1 + u)
        return float(w / (scale + w))


def modular_potential_xy(x: float, y: float, N: int = 12, mode: str = "ratio",
                         scale: float = 10.0, dps: int = 25) -> float:
    """`modular_potential` in the ``(Re tau, Im tau)`` signature the projections use."""
    return modular_potential(complex(x, y), N=N, mode=mode, scale=scale, dps=dps)
