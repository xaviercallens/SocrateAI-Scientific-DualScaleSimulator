"""
Spec item L1 (`specs/LEANFLOW_ARCHITECTURE.md`): a potential that is actually
Gamma_0(N)+ -invariant.

These tests ARE the acceptance criterion written in that spec. They also carry the
negative controls the criterion demands, so that "invariant" cannot be satisfied
by a test that could not fail.

No physics is claimed. Stream 1 records program-wide that no exact physical
observable exists anywhere in this programme.
"""

import math
import warnings

import mpmath
import numpy as np
import pytest

warnings.filterwarnings("ignore")

from leanflow.core.gamma0n_plus import is_invariant_under  # noqa: E402
from leanflow.core.modular_potential import (  # noqa: E402
    gamma0_plus_invariant,
    j_invariant,
    modular_potential,
    modular_potential_xy,
    self_dual_tau,
)

N = 12  # the level this repository is implicitly at; see audit/K3_SELECTION.md


# ---------------------------------------------------------------------------
# The j-invariant, against known values
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("tau,expected", [
    (complex(0, 1), 1728),                                   # j(i) = 1728
    (complex(0, 2), 287496),                                 # j(2i) = 66^3
    (complex(0, math.sqrt(2)), 8000),                        # j(i sqrt 2) = 20^3
])
def test_j_invariant_known_values(tau, expected):
    assert abs(complex(j_invariant(tau)) - expected) < 1e-6 * max(1, abs(expected))


def test_j_vanishes_at_the_corner():
    """j(exp(2 pi i / 3)) = 0 -- the order-3 elliptic point."""
    rho = complex(math.cos(2 * math.pi / 3), math.sin(2 * math.pi / 3))
    assert abs(complex(j_invariant(rho))) < 1e-15


# ---------------------------------------------------------------------------
# F_N(tau) = j(tau) + j(N tau) is Gamma_0(N)+ -invariant
# ---------------------------------------------------------------------------

_DPS = 30


@pytest.mark.parametrize("z", [(0.13, 0.9), (-0.31, 0.45), (0.42, 1.7)])
def test_F_is_fricke_invariant(z):
    """W_N swaps j(tau) and j(N tau), so their sum is fixed.

    The transformation is done IN mpmath: `F_N` is so steep (|F| reaches 1e55 in
    this domain) that building `-1/(N tau)` in float64 first costs ~13 digits of
    the comparison. See test_float64_input_costs_precision, which records that.
    """
    with mpmath.workdps(_DPS):
        tau = mpmath.mpc(*z)
        img = -1 / (N * tau)          # transformation computed at high precision
    a = gamma0_plus_invariant(tau, N, dps=_DPS)
    b = gamma0_plus_invariant(img, N, dps=_DPS)
    assert abs(b - a) / abs(a) < 1e-25


@pytest.mark.parametrize("g", [(1, 1, 0, 1), (1, 0, 1, 1), (5, 1, 2, 5)])
def test_F_is_gamma0N_invariant(g):
    """For g = [[a,b],[Nc,d]] with ad - Nbc = 1."""
    a, b, c, d = g
    assert a * d - N * b * c == 1, "test matrix must lie in Gamma_0(N)"
    for z in ((0.13, 0.9), (-0.31, 0.45)):
        with mpmath.workdps(_DPS):
            tau = mpmath.mpc(*z)
            gt = (a * tau + b) / (N * c * tau + d)
        x = gamma0_plus_invariant(tau, N, dps=_DPS)
        y = gamma0_plus_invariant(gt, N, dps=_DPS)
        assert abs(y - x) / abs(x) < 1e-25


def test_float64_input_costs_precision():
    """RECORDED, not a defect of the invariance: the identity is exact to working
    precision (~1e-30 at dps=30), but `F_N` is steep enough that a float64 `tau`
    degrades the check to ~1e-13. Anyone wiring this into a float64 solver must
    know that the invariance is not usable at double precision for large Im tau.
    """
    z = complex(0.13, 0.9)
    with mpmath.workdps(_DPS):
        exact = mpmath.mpc(z.real, z.imag)
        img_hi = -1 / (N * exact)                 # high-precision transformation
    img_lo = -1 / (N * z)                         # the same map done in float64
    base = gamma0_plus_invariant(exact, N, dps=_DPS)
    hi = abs(gamma0_plus_invariant(img_hi, N, dps=_DPS) - base) / abs(base)
    lo = abs(gamma0_plus_invariant(img_lo, N, dps=_DPS) - base) / abs(base)
    assert float(hi) < 1e-25, "the identity itself is exact"
    assert float(lo) > 1e-16, "float64 input measurably degrades it"


def test_F_is_NOT_invariant_under_S_negative_control():
    """NEGATIVE CONTROL. S : tau -> -1/tau is the N = 1 Fricke map and must NOT fix
    F_N for N > 1 -- otherwise the invariance test above is checking nothing."""
    tau = complex(0.13, 0.9)
    a = gamma0_plus_invariant(tau, N)
    b = gamma0_plus_invariant(-1 / tau, N)
    assert abs(b - a) / abs(a) > 1e-3


# ---------------------------------------------------------------------------
# L1 ACCEPTANCE CRITERION, verbatim from the spec
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("mode", ["ratio", "log"])
def test_L1_acceptance_criterion(mode):
    """spec L1: both checks must return True, using gamma0n_plus.is_invariant_under."""
    V = lambda x, y: modular_potential_xy(x, y, N=N, mode=mode)
    ok_T, dev_T = is_invariant_under(V, lambda z: z + 1, samples=80, atol=1e-9, seed=3)
    ok_F, dev_F = is_invariant_under(V, lambda z: -1 / (N * z), samples=80, atol=1e-9, seed=3)
    assert ok_T, f"T invariance failed, deviation {dev_T}"
    assert ok_F, f"Fricke invariance failed, deviation {dev_F}"


def test_L1_negative_control_the_old_potential_still_fails():
    """The criterion must be able to fail. The hand-built double well this
    repository actually integrates is T-invariant and NOT Fricke-invariant
    (audit/STREAM1_BRIDGE.md S1-F2); if this ever passes, the check is broken."""
    import workshopcosmo as wc
    V0 = lambda x, y: wc.compute_potential(x, y)[0]
    ok_T, _ = is_invariant_under(V0, lambda z: z + 1, samples=80, atol=1e-9, seed=3)
    ok_F, dev_F = is_invariant_under(V0, lambda z: -1 / (N * z), samples=80, atol=1e-9, seed=3)
    assert ok_T, "the old potential IS T-invariant"
    assert not ok_F and dev_F > 1.0, "the old potential must still fail Fricke"


def test_L1_negative_control_a_perturbation_breaks_invariance():
    """Spec L1: 'a deliberately non-invariant perturbation must fail the same check'."""
    V_bad = lambda x, y: modular_potential_xy(x, y, N=N) + 0.05 * y
    ok_F, dev_F = is_invariant_under(V_bad, lambda z: -1 / (N * z), samples=60, atol=1e-9, seed=3)
    assert not ok_F and dev_F > 1e-3


# ---------------------------------------------------------------------------
# Shape of the well
# ---------------------------------------------------------------------------

def test_potential_vanishes_on_the_whole_orbit_not_just_the_point():
    """V = 0 at i/sqrt(N) AND at its images under T and W_N -- the signature of a
    genuine invariant rather than a well dug at one place by hand."""
    t0 = self_dual_tau(N)
    for image in (t0, t0 + 1, t0 - 1, -1 / (N * t0)):
        assert modular_potential(image, N) < 1e-12


def test_ratio_mode_minimum_is_quadratic_and_stationary():
    """Measured slope 2.000 over four decades, and a vanishing gradient."""
    t0 = self_dual_tau(N)
    def radial(r):
        return sum(modular_potential(t0 + complex(r * math.cos(a), r * math.sin(a)), N)
                   for a in (0.0, 1.1, 2.3, 3.7)) / 4.0
    v3, v5 = radial(1e-3), radial(1e-5)
    slope = math.log(v3 / v5) / math.log(1e-3 / 1e-5)
    assert abs(slope - 2.0) < 0.05, f"expected V ~ r^2, got r^{slope:.3f}"
    h = 1e-6
    dx = (modular_potential(t0 + h, N) - modular_potential(t0 - h, N)) / (2 * h)
    dy = (modular_potential(t0 + 1j * h, N) - modular_potential(t0 - 1j * h, N)) / (2 * h)
    assert abs(dx) < 1e-6 and abs(dy) < 1e-6


def test_bounded_in_unit_interval():
    for tau in (complex(0, 1), complex(0.3, 0.9), complex(0.1, 5.0), complex(0.45, 0.3)):
        for mode in ("ratio", "log"):
            v = modular_potential(tau, N, mode=mode)
            assert 0.0 <= v <= 1.0


def test_rejects_unknown_mode():
    with pytest.raises(ValueError):
        modular_potential(complex(0, 1), N, mode="nope")


# ---------------------------------------------------------------------------
# L5 precondition: does the LEVEL have any consequence at all?
# ---------------------------------------------------------------------------

def test_level_changes_the_landscape_unlike_the_old_potential():
    """`audit/K3_SELECTION.md` §6.3 recorded the obstruction: changing `N` moved
    `FRICKE_Y` and nothing else, so no observable depended on the level and no
    level could be falsified. This test is the regression guard for that.

    The old potential does not read `N` at runtime at all -- `N` enters only through
    a module-level constant -- so its values are bit-identical across levels. The
    modular potential depends on `N` everywhere.
    """
    import random
    import workshopcosmo as wc

    rng = random.Random(11)
    pts = [complex(rng.uniform(-0.5, 0.5), 10 ** rng.uniform(-0.5, 0.5)) for _ in range(8)]

    # Old: it cannot depend on the level, because it has no level parameter.
    # (Asserting that two calls agree would only test determinism.)
    import inspect
    params = inspect.signature(wc.compute_potential).parameters
    assert "N" not in params and "level" not in params, \
        f"compute_potential grew a level argument: {list(params)}"

    # New: the level moves the landscape at every probe.
    devs = [abs(modular_potential(t, 7, mode="log") - modular_potential(t, 12, mode="log"))
            for t in pts]
    assert min(devs) > 1e-3, f"level dependence vanished somewhere: min={min(devs)}"
    assert sum(devs) / len(devs) > 0.01

    # MOST OF THAT IS A REPARAMETRIZATION, not structure. In log mode V is
    # essentially a function of N*Im(tau); rescaling so that 12*Im(tau') = 7*Im(tau)
    # removes ~76% of the apparent difference. Recorded so the number in
    # audit/K3_SELECTION.md 6.3 is not read as stronger than it is.
    resc = [abs(modular_potential(t, 7, mode="log")
                - modular_potential(complex(t.real, t.imag * 7 / 12), 12, mode="log"))
            for t in pts]
    assert sum(resc) / len(resc) < 0.5 * sum(devs) / len(devs)

    # Ratio mode is different: the rescaling does NOT remove its level dependence,
    # so that one is genuine structure rather than a stretched axis.
    d_r = [abs(modular_potential(t, 7) - modular_potential(t, 12)) for t in pts]
    r_r = [abs(modular_potential(t, 7)
               - modular_potential(complex(t.real, t.imag * 7 / 12), 12)) for t in pts]
    assert sum(r_r) / len(r_r) > 0.8 * sum(d_r) / len(d_r)

    # ... and the wells sit at different places, as they must.
    assert abs(self_dual_tau(7) - self_dual_tau(12)) > 0.08
    assert modular_potential(self_dual_tau(7), 7) < 1e-12
    assert modular_potential(self_dual_tau(7), 12) > 1e-3   # not a well at the wrong level
