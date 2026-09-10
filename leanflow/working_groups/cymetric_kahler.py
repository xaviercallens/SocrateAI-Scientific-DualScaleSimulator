"""
=============================================================================
Working Group 1: Numerical Calabi-Yau Metric Stabilization (cymetric / PINNs)
=============================================================================
Provides active projection operators and loss regularizers that enforce:
- Kahler cone positive-definiteness: g_{i, jbar} > 0
- Non-vanishing Monge-Ampere volume forms: det(g) > 0
- Machine-precision Ricci-flatness stabilization for neural metric surrogates.
=============================================================================
"""

from typing import Union, Tuple, Optional
import numpy as np

try:
    import torch
    import torch.nn as nn
    TORCH_AVAILABLE = True
except ImportError:
    TORCH_AVAILABLE = False


def project_kahler_metric(
    g: Union[np.ndarray, "torch.Tensor"],
    min_eigenval: float = 1e-4
) -> Union[np.ndarray, "torch.Tensor"]:
    """
    Stabilizes candidate Calabi-Yau metric tensors by projecting them
    onto the interior of the Kahler cone (strictly positive-definite).
    Precludes metric degeneracies in neural networks and PINN solvers.
    """
    if TORCH_AVAILABLE and isinstance(g, torch.Tensor):
        # Hermitian / symmetric eigendecomposition
        evals, evecs = torch.linalg.eigh(g)
        evals_proj = torch.clamp(evals, min=min_eigenval)
        return torch.matmul(evecs, torch.matmul(torch.diag_embed(evals_proj), evecs.mH))
    else:
        evals, evecs = np.linalg.eigh(g)
        evals_proj = np.maximum(evals, min_eigenval)
        if g.ndim == 2:
            return evecs @ np.diag(evals_proj) @ evecs.conj().T
        else:
            g_out = np.zeros_like(g)
            for i in range(g.shape[0]):
                g_out[i] = evecs[i] @ np.diag(evals_proj[i]) @ evecs[i].conj().T
            return g_out


def compute_monge_ampere_loss(
    g: Union[np.ndarray, "torch.Tensor"],
    omega_wedge_omega_bar: Union[np.ndarray, "torch.Tensor"]
) -> Union[float, "torch.Tensor"]:
    """
    Computes the Monge-Ampere Ricci-flatness residual:
    L_MA = || det(g) / (Omega ^ Omega_bar) - 1 ||^2
    with protected positive determinants.
    """
    if TORCH_AVAILABLE and isinstance(g, torch.Tensor):
        det_g = torch.linalg.det(g).real
        det_g_safe = torch.clamp(det_g, min=1e-12)
        ratio = det_g_safe / torch.clamp(omega_wedge_omega_bar, min=1e-12)
        return torch.mean((ratio - 1.0) ** 2)
    else:
        det_g = np.real(np.linalg.det(g))
        det_g_safe = np.maximum(det_g, 1e-12)
        ratio = det_g_safe / np.maximum(omega_wedge_omega_bar, 1e-12)
        return float(np.mean((ratio - 1.0) ** 2))


if TORCH_AVAILABLE:
    class CymetricKahlerRegularizer(nn.Module):
        """
        PyTorch loss regularizer for Calabi-Yau metric learning (cymetric).
        Penalizes metric eigenvalues leaving the Kahler cone.
        """
        def __init__(self, min_eigenval: float = 1e-4, weight: float = 10.0):
            super().__init__()
            self.min_eigenval = min_eigenval
            self.weight = weight

        def forward(self, g_pred: torch.Tensor) -> torch.Tensor:
            evals = torch.linalg.eigvalsh(g_pred)
            violations = torch.clamp(self.min_eigenval - evals, min=0.0)
            loss = torch.sum(violations ** 2) * self.weight
            return loss
