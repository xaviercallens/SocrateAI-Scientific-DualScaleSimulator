"""lambda_sym deletion: equivalence test and negative control.

Claim: the psi-form solver (workshopcosmo.run_symmetron_screening_simulation, no
lambda_sym) returns the same retained screening observables as the legacy phi-form
solver (legacy_phi_form.py, with lambda_sym) for every lambda_sym.

Negative control: a psi-form solver with a deliberately lambda-dependent cubic term
(mu^2 * lambda * psi^3, i.e. a wrong substitution) must DISAGREE with the legacy
solver for lambda != 1. If it agreed, the comparison would have no power.

Run: .venv-tda/bin/python -m pytest -q audit/zero_param_loop/lambda_deletion/test_lambda_deletion.py
"""
import inspect
import os
import sys

import numpy as np
import pytest
from scipy.integrate import solve_bvp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)

import workshopcosmo as wc  # noqa: E402
from legacy_phi_form import legacy_phi_form  # noqa: E402

LAMBDAS = [0.01, 0.1, 1.0, 10.0, 100.0]
MUS = [0.3, 1.0, 2.0]
KEYS = ["screening_suppression_factor", "phi_center_ratio", "phi_surface_ratio"]
# solve_bvp runs with tol=1e-4; the legacy phi-form scales its residual by phi_0,
# so the two solvers stop at slightly different iterates. Tolerance: 1e-3 relative
# (plus an absolute floor for ratios that are ~1e-14).
RTOL, ATOL = 1e-3, 1e-12


def _close(a, b):
    return abs(a - b) <= ATOL + RTOL * max(abs(a), abs(b))


def test_lambda_is_gone_from_the_model():
    sig = inspect.signature(wc.run_symmetron_screening_simulation)
    assert "lambda_sym" not in sig.parameters
    src = inspect.getsource(wc.run_symmetron_screening_simulation)
    code = "\n".join(l for l in src.splitlines() if not l.strip().startswith("#"))
    body = code.split('"""')[-1]  # after the docstring
    assert "lambda" not in body


@pytest.mark.parametrize("mu", MUS)
@pytest.mark.parametrize("lam", LAMBDAS)
def test_psi_form_equals_phi_form(mu, lam):
    new = wc.run_symmetron_screening_simulation(mu_sym=mu)
    old = legacy_phi_form(mu_sym=mu, lambda_sym=lam)
    for k in KEYS:
        assert np.isfinite(new[k]) and np.isfinite(old[k])
        assert _close(new[k], old[k]), (k, mu, lam, new[k], old[k])


def test_mu_is_live():
    # positive control: the remaining parameter really moves the observable
    a = wc.run_symmetron_screening_simulation(mu_sym=1.0)["screening_suppression_factor"]
    b = wc.run_symmetron_screening_simulation(mu_sym=0.3)["screening_suppression_factor"]
    assert not _close(a, b)


def _wrong_psi_form(mu, lam, r_max=10.0, n_points=400, rho_in=1000.0, rho_out=0.01):
    """Deliberately wrong substitution: keeps a lambda in the cubic term."""
    r_core, eps = 1.0, 1e-4
    r = np.concatenate([
        np.linspace(eps, r_core - 0.2, n_points // 3, endpoint=False),
        np.linspace(r_core - 0.2, r_core + 0.2, n_points // 3, endpoint=False),
        np.linspace(r_core + 0.2, r_max, n_points - 2 * (n_points // 3)),
    ])
    rho = lambda x: rho_out + (rho_in - rho_out) / (1.0 + np.exp((x - r_core) / 0.05))

    def f(x, y):
        return np.vstack((y[1], (rho(x) - mu ** 2) * y[0] + mu ** 2 * lam * y[0] ** 3 - 2.0 / x * y[1]))

    g = np.zeros((2, n_points))
    g[0] = np.where(r < r_core, 0.001, 1.0 - np.exp(-(r - r_core)))
    g[1] = np.gradient(g[0], r)
    s = solve_bvp(f, lambda ya, yb: np.array([ya[1], yb[0] - 1.0]), r, g, max_nodes=5000, tol=1e-4)
    i = int(np.argmin(np.abs(s.x - r_core)))
    return {"screening_suppression_factor": float(s.y[0][i] ** 2), "success": bool(s.success)}


def test_negative_control_wrong_substitution_is_detected():
    disagreements = 0
    for lam in (0.01, 100.0):
        wrong = _wrong_psi_form(1.0, lam)
        old = legacy_phi_form(mu_sym=1.0, lambda_sym=lam)
        if wrong["success"] and not _close(wrong["screening_suppression_factor"], old["screening_suppression_factor"]):
            disagreements += 1
    assert disagreements >= 1, "negative control failed: a lambda-dependent psi-form was not detected"
