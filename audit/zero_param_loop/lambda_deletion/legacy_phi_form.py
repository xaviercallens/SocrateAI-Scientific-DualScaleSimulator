"""Verbatim copy of workshopcosmo.run_symmetron_screening_simulation as of commit
10d3ed2 (phi-form, with lambda_sym), renamed legacy_phi_form. Kept only as the
reference for the psi-form equivalence test. Do not use in the model."""
import math
from typing import Any, Dict
import numpy as np
from scipy.integrate import solve_bvp


def legacy_phi_form(
    r_max: float = 10.0,
    n_points: int = 400,
    rho_in: float = 1000.0,
    rho_out: float = 0.01,
    mu_sym: float = 1.0,
    lambda_sym: float = 1.0,
    m_scale: float = 1.0,
) -> Dict[str, Any]:
    """
    Solves the non-linear radial Symmetron equation:
    d²phi/dr² + (2/r) dphi/dr = dV_eff/dphi
    V_eff(phi) = 0.5 * (rho(r)/M² - mu²) phi² + 0.25 * lambda * phi⁴

    Inside body (r <= 1.0): rho = rho_in >> mu² M² => phi -> 0 (symmetry restored).
    Outside body (r > 1.0): rho = rho_out << mu² M² => phi -> phi_0 = mu/sqrt(lambda).
    Boundary conditions: dphi/dr(0) = 0, phi(r_max) = phi_0.
    """
    r_core = 1.0
    phi_0 = mu_sym / math.sqrt(lambda_sym)  # vacuum expectation value (VEV)

    # Regularized radial coordinate to avoid 1/r singularity at r=0
    eps = 1e-4
    # Refined grid: higher density around r_core to capture sharp screening transitions
    r_grid = np.concatenate([
        np.linspace(eps, r_core - 0.2, n_points // 3, endpoint=False),
        np.linspace(r_core - 0.2, r_core + 0.2, n_points // 3, endpoint=False),
        np.linspace(r_core + 0.2, r_max, n_points - 2 * (n_points // 3))
    ])

    def rho_profile(r: np.ndarray) -> np.ndarray:
        # Smooth step across surface
        width = 0.05
        return rho_out + (rho_in - rho_out) / (1.0 + np.exp((r - r_core) / width))

    def ode_system(r: np.ndarray, y: np.ndarray) -> np.ndarray:
        phi = y[0]
        dphi_dr = y[1]
        rho = rho_profile(r)
        # dV_eff/dphi = (rho/M² - mu²) phi + lambda * phi³
        dv_dphi = ((rho / (m_scale ** 2)) - (mu_sym ** 2)) * phi + lambda_sym * (phi ** 3)
        d2phi_dr2 = dv_dphi - (2.0 / r) * dphi_dr
        return np.vstack((dphi_dr, d2phi_dr2))

    def bc(ya: np.ndarray, yb: np.ndarray) -> np.ndarray:
        # At r = eps: dphi/dr = 0
        # At r = r_max: phi = phi_0
        return np.array([ya[1], yb[0] - phi_0])

    # Initial guess
    y_guess = np.zeros((2, n_points))
    # Step function transition from 0 to phi_0 outside core
    y_guess[0] = np.where(r_grid < r_core, 0.001 * phi_0, phi_0 * (1.0 - np.exp(-(r_grid - r_core))))
    y_guess[1] = np.gradient(y_guess[0], r_grid)

    res = solve_bvp(ode_system, bc, r_grid, y_guess, max_nodes=5000, tol=1e-4)

    if not res.success:
        # Refined relaxation fallback
        phi_sol = y_guess[0]
        for _ in range(50):
            dphi = np.gradient(phi_sol, r_grid)
            d2phi = np.gradient(dphi, r_grid)
            rho = rho_profile(r_grid)
            res_val = d2phi + (2.0 / r_grid) * dphi - (((rho / (m_scale ** 2)) - (mu_sym ** 2)) * phi_sol + lambda_sym * (phi_sol ** 3))
            phi_sol -= 0.01 * res_val
            phi_sol[0] = phi_sol[1]
            phi_sol[-1] = phi_0
        r_eval = r_grid
        phi_eval = np.clip(phi_sol, 0.0, phi_0)
        dphi_eval = np.gradient(phi_eval, r_eval)
    else:
        r_eval = res.x
        phi_eval = res.y[0]
        dphi_eval = res.y[1]

    # Compute physical screening parameters
    idx_center = 0
    idx_surface = int(np.argmin(np.abs(r_eval - r_core)))
    phi_center_ratio = float(phi_eval[idx_center] / phi_0)
    phi_surface_ratio = float(phi_eval[idx_surface] / phi_0)

    # Fifth force ratio F_phi / F_N = (dphi/dr * phi / M) / (G M_enc / r²)
    # For thin-shell screening: Delta R / R = (phi_out - phi_in) / (6 beta M_Pl Phi_N)
    # Fifth-force suppression factor ~ (phi_surface / phi_0)²
    screening_factor = float(phi_surface_ratio ** 2)
    is_screened = (phi_center_ratio < 0.05) and (screening_factor < 5e-4)

    return {
        "success": True,
        "r": r_eval.tolist(),
        "phi_ratio": (phi_eval / phi_0).tolist(),
        "dphi_dr": dphi_eval.tolist(),
        "phi_0": phi_0,
        "phi_center_ratio": phi_center_ratio,
        "phi_surface_ratio": phi_surface_ratio,
        "screening_suppression_factor": screening_factor,
        "is_screened": is_screened,
        "cassini_bound_satisfied": screening_factor < 5e-4,
    }
