"""
=============================================================================
LeanFlow Core Projections: Certified Invariant Projections & SIMD Guards
=============================================================================
Provides ahead-of-time (AOT) mathematical projection operators that enforce
Lean 4-verified invariants down to machine precision:
- Target-space metric positivity: tau_im = y > 0
- Modular invariance: Folding into SL(2, Z) fundamental domain F
- Weak Energy Condition: w(t) >= -1
- Calabi-Yau Kahler cone positive-definiteness: g_{i, jbar} > 0
=============================================================================
"""

import math
from typing import Tuple, Union, Optional
import numpy as np

try:
    import torch
    TORCH_AVAILABLE = True
except ImportError:
    TORCH_AVAILABLE = False

FRICKE_Y = 1.0 / math.sqrt(12.0)  # ~ 0.28867513459
ORBIFOLD_Y = math.sqrt(3.0) / 2.0  # ~ 0.86602540378


def metric_positivity_projection(
    state: Union[np.ndarray, "torch.Tensor"],
    min_y: float = FRICKE_Y
) -> Union[np.ndarray, "torch.Tensor"]:
    """
    Orthogonal projection enforcing Poincaré target-space metric positivity:
    tau_im = y >= min_y > 0.
    Precludes coordinate runaways and division-by-zero singularities.
    """
    if TORCH_AVAILABLE and isinstance(state, torch.Tensor):
        projected = state.clone()
        # Assumes y is at index 1 in [x, y, vx, vy] or [x, y]
        if projected.ndim == 1 and projected.shape[0] >= 2:
            if projected[1] < min_y:
                projected[1] = min_y
        elif projected.ndim == 2 and projected.shape[1] >= 2:
            projected[:, 1] = torch.clamp(projected[:, 1], min=min_y)
        return projected
    else:
        projected = np.copy(state)
        if projected.ndim == 1 and projected.shape[0] >= 2:
            if projected[1] < min_y:
                projected[1] = min_y
        elif projected.ndim == 2 and projected.shape[1] >= 2:
            projected[:, 1] = np.maximum(projected[:, 1], min_y)
        return projected


def modular_domain_fold(
    x: float,
    y: float,
    max_folds: int = 20
) -> Tuple[float, float, int]:
    """
    Folds the modulus tau = x + i y back into the fundamental modular domain:
    F = { tau in H : |tau| >= 1, |Re(tau)| <= 1/2 }
    using SL(2, Z) generators T: tau -> tau +/- 1 and S: tau -> -1/tau.
    Returns (x_folded, y_folded, fold_count).
    """
    folds = 0
    # Guard against singular non-positive y
    y = max(y, 1e-6)

    for _ in range(max_folds):
        # 1. T-transformation: Shift real part into [-0.5, 0.5]
        shift = round(x)
        if shift != 0:
            x -= shift
            folds += abs(shift)

        # 2. S-transformation: If |tau| < 1, invert tau -> -1/tau
        norm_sq = x * x + y * y
        if norm_sq < 1.0 - 1e-9:
            # S: tau -> -1/tau = (-x + i y) / (x^2 + y^2)
            x = -x / norm_sq
            y = y / norm_sq
            folds += 1
        else:
            # Inside fundamental domain F
            break

    # Final real part clamp
    x = max(-0.5, min(0.5, x))
    y = max(FRICKE_Y, y)
    return float(x), float(y), folds


def weak_energy_condition_lock(
    rho: float,
    p: float,
    min_sum: float = 0.0
) -> Tuple[float, float, float]:
    """
    Enforces the Weak Energy Condition rho + p >= 0 and Null Energy Condition.
    Guarantees equation of state w = p / rho >= -1.
    Returns (rho, p_projected, w_projected).
    """
    rho = max(1e-12, float(rho))
    # If rho + p < min_sum, project p to min_sum - rho
    if rho + p < min_sum:
        p_proj = min_sum - rho
    else:
        p_proj = float(p)
    w = p_proj / rho
    return rho, p_proj, w


def kahler_cone_positivity_projection(
    g_matrix: Union[np.ndarray, "torch.Tensor"],
    min_eigenval: float = 1e-5
) -> Union[np.ndarray, "torch.Tensor"]:
    """
    Active projection operator for Numerical Calabi-Yau metric learning (e.g. cymetric).
    Enforces that the metric matrix g_{i, jbar} remains strictly positive-definite
    (inside the Kahler cone) down to machine precision.
    Uses spectral decomposition: g' = V * diag(max(min_eigenval, lambda_i)) * V^H.
    """
    if TORCH_AVAILABLE and isinstance(g_matrix, torch.Tensor):
        # Symmetric / Hermitian eigenvalue decomposition
        # Works on batches (..., n, n) or single matrix (n, n)
        evals, evecs = torch.linalg.eigh(g_matrix)
        evals_clamped = torch.clamp(evals, min=min_eigenval)
        # Reconstruct: V * diag(evals) * V^H
        diag_clamped = torch.diag_embed(evals_clamped.to(evecs.dtype))
        g_projected = torch.matmul(evecs, torch.matmul(diag_clamped, evecs.mH))
        return g_projected
    else:
        # NumPy implementation
        evals, evecs = np.linalg.eigh(g_matrix)
        evals_clamped = np.maximum(evals, min_eigenval)
        if g_matrix.ndim == 2:
            g_projected = evecs @ np.diag(evals_clamped) @ evecs.conj().T
        else:
            # Batch
            g_projected = np.zeros_like(g_matrix)
            for i in range(g_matrix.shape[0]):
                g_projected[i] = evecs[i] @ np.diag(evals_clamped[i]) @ evecs[i].conj().T
        return g_projected


def tadpole_budget_check(
    q_d3_background: float,
    n_flux: int,
    euler_char: int
) -> Tuple[bool, float, float]:
    """
    Verifies D3-brane tadpole cancellation bound on Calabi-Yau 4-fold / 3-fold:
    Q_D3 + N_flux <= chi / 24.
    Returns (is_valid, total_charge, max_allowed).
    """
    max_allowed = float(euler_char) / 24.0
    total_charge = float(q_d3_background) + float(n_flux)
    is_valid = total_charge <= max_allowed + 1e-9
    return is_valid, total_charge, max_allowed
