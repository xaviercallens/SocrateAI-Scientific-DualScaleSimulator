"""
=============================================================================
Neural Calabi-Yau Metric Learning PINN for the Quintic Threefold P4[5]
=============================================================================
Defines a Physics-Informed Neural Network approximating the Ricci-flat metric
on the Fermat Quintic threefold:
    z0^5 + z1^5 + z2^5 + z3^5 + z4^5 = 0
Equipped with LeanFlow's active spectral projection onto the Kahler cone:
    g'_{i, jbar} = V * diag(max(min_eval, lambda_i)) * V^H
=============================================================================
"""

import math
from typing import Tuple, Optional
import numpy as np
import torch
import torch.nn as nn

from leanflow.core.projections import kahler_cone_positivity_projection


class QuinticMetricPINN(nn.Module):
    """
    Physics-Informed Neural Network predicting the 3x3 Hermitian metric tensor
    g_{i, jbar} on the tangent bundle of the Quintic threefold in local charts.
    """
    def __init__(self, input_dim: int = 6, hidden_dim: int = 64, num_layers: int = 3):
        super().__init__()
        self.input_dim = input_dim
        layers = []
        in_d = input_dim
        for _ in range(num_layers):
            layers.append(nn.Linear(in_d, hidden_dim))
            layers.append(nn.GELU())
            in_d = hidden_dim
        self.backbone = nn.Sequential(*layers)
        
        # Metric output: 3 real diagonal components + 3 complex off-diagonal components
        # Total parameters for 3x3 Hermitian matrix = 3 + 2 * 3 = 9 real values
        self.head_diag = nn.Linear(hidden_dim, 3)
        self.head_offdiag_re = nn.Linear(hidden_dim, 3)
        self.head_offdiag_im = nn.Linear(hidden_dim, 3)
        self.to(torch.float64)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Raw unconstrained forward pass.
        May produce indefinite matrices (drift outside Kahler cone) during early training.
        Returns batch of 3x3 complex Hermitian matrices.
        """
        x = x.to(torch.float64)
        batch_size = x.shape[0]
        h = self.backbone(x)
        
        # Diagonal elements (real)
        d = self.head_diag(h)  # (B, 3)
        
        # Off-diagonal elements (real and imaginary parts for indices (0,1), (0,2), (1,2))
        re = self.head_offdiag_re(h)  # (B, 3)
        im = self.head_offdiag_im(h)  # (B, 3)
        
        # Construct complex Hermitian batch: g[b, i, j]
        g = torch.zeros((batch_size, 3, 3), dtype=torch.complex128, device=x.device)
        
        # Diagonals
        g[:, 0, 0] = torch.complex(d[:, 0], torch.zeros_like(d[:, 0]))
        g[:, 1, 1] = torch.complex(d[:, 1], torch.zeros_like(d[:, 1]))
        g[:, 2, 2] = torch.complex(d[:, 2], torch.zeros_like(d[:, 2]))
        
        # (0, 1) and (1, 0)
        c01 = torch.complex(re[:, 0], im[:, 0])
        g[:, 0, 1] = c01
        g[:, 1, 0] = c01.conj()
        
        # (0, 2) and (2, 0)
        c02 = torch.complex(re[:, 1], im[:, 1])
        g[:, 0, 2] = c02
        g[:, 2, 0] = c02.conj()
        
        # (1, 2) and (2, 1)
        c12 = torch.complex(re[:, 2], im[:, 2])
        g[:, 1, 2] = c12
        g[:, 2, 1] = c12.conj()
        
        return g

    def forward_guarded(
        self,
        x: torch.Tensor,
        min_eigenval: float = 1e-3
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Forward pass protected by LeanFlow's active Kahler cone spectral projection:
        Returns: (stabilized_metric, minimum_eigenvalues)
        Guarantees strict positive-definiteness g_{i, jbar} > 0.
        """
        g_raw = self.forward(x)
        # Apply LeanFlow orthogonal spectral projection
        g_stabilized = kahler_cone_positivity_projection(g_raw, min_eigenval=min_eigenval)
        
        # Compute eigenvalues of stabilized metric
        evals, _ = torch.linalg.eigh(g_stabilized)
        min_evals = evals[:, 0]
        return g_stabilized, min_evals

    def compute_monge_ampere_loss(
        self,
        x: torch.Tensor,
        omega_sq: torch.Tensor,
        use_guardrail: bool = True,
        min_eigenval: float = 1e-3
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Evaluates the normalized Monge-Ampère loss:
            L_MA = (1 / N) * sum | det(g) / omega_sq - 1 |^2
        """
        if use_guardrail:
            g, min_evals = self.forward_guarded(x, min_eigenval=min_eigenval)
        else:
            g = self.forward(x)
            evals, _ = torch.linalg.eigh(g)
            min_evals = evals[:, 0]
            
        # det(g) for 3x3 Hermitian matrix is real
        det_g = torch.linalg.det(g).real
        
        # Monge-Ampère residual
        ratio = det_g / (omega_sq + 1e-12)
        loss = torch.mean((ratio - 1.0) ** 2)
        return loss, min_evals


def sample_quintic_chart_points(num_points: int = 100, seed: int = 42) -> Tuple[torch.Tensor, torch.Tensor]:
    """
    Generates synthetic evaluation points in an affine chart of the Quintic threefold.
    Returns: (x_coords, omega_sq)
    """
    torch.manual_seed(seed)
    # Affine coordinates: 3 complex variables = 6 real dimensions
    x = torch.randn((num_points, 6), dtype=torch.float64) * 0.5
    
    # Holomorphic volume form approximation |Omega|^2
    r_sq = torch.sum(x**2, dim=1)
    omega_sq = 1.0 / (1.0 + r_sq)**4
    return x, omega_sq
