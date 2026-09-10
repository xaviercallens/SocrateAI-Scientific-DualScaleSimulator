"""
=============================================================================
LeanFlow: Generalized Geometry & Double Field Theory (DFT) Framework
=============================================================================
Provides full 2D x 2D generalized metric construction, O(D, D; Z) invariant
verification, Courant / Dorfman brackets, and non-geometric flux chain tracking:
    H_{ijk} ---> f^i_{jk} ---> Q^{ij}_k ---> R^{ijk}
=============================================================================
"""

from typing import Tuple, Dict, Any, Optional, List
import numpy as np


class GeneralizedMetric:
    """
    Represents the 2D x 2D Generalized Metric H_{MN} on the generalized
    tangent bundle TM + T*M in Double Field Theory (DFT).
    """
    def __init__(self, g: np.ndarray, b: Optional[np.ndarray] = None):
        """
        Initializes the generalized metric from target space metric g (D x D)
        and Kalb-Ramond 2-form B (D x D, antisymmetric).
        """
        g = np.asarray(g, dtype=np.float64)
        if g.ndim != 2 or g.shape[0] != g.shape[1]:
            raise ValueError(f"Metric g must be a square matrix, got shape {g.shape}")
        
        self.d = g.shape[0]
        
        # Verify metric symmetry
        if not np.allclose(g, g.T, atol=1e-8):
            raise ValueError("Target space metric g must be symmetric")
            
        # Verify positive-definiteness
        evals = np.linalg.eigvalsh(g)
        if np.min(evals) <= 0:
            raise ValueError(f"Target space metric g must be positive-definite, min eval = {np.min(evals)}")
            
        self.g = g
        self.g_inv = np.linalg.inv(g)

        if b is None:
            self.b = np.zeros((self.d, self.d), dtype=np.float64)
        else:
            b = np.asarray(b, dtype=np.float64)
            if b.shape != (self.d, self.d):
                raise ValueError(f"B-field must have shape ({self.d}, {self.d}), got {b.shape}")
            if not np.allclose(b, -b.T, atol=1e-8):
                raise ValueError("B-field must be antisymmetric (B = -B^T)")
            self.b = b

        # Construct O(D, D) invariant metric: eta = [[0, I], [I, 0]]
        i_d = np.eye(self.d, dtype=np.float64)
        z_d = np.zeros((self.d, self.d), dtype=np.float64)
        self.eta = np.block([[z_d, i_d], [i_d, z_d]])

        # Construct 2D x 2D Generalized Metric:
        # H = [[g - B g^{-1} B,  B g^{-1}],
        #      [-g^{-1} B,       g^{-1}  ]]
        b_ginv = self.b @ self.g_inv
        ginv_b = self.g_inv @ self.b
        top_left = self.g - (self.b @ self.g_inv @ self.b)
        top_right = b_ginv
        bottom_left = -ginv_b
        bottom_right = self.g_inv

        self.h = np.block([
            [top_left, top_right],
            [bottom_left, bottom_right]
        ])

    def verify_odd_invariance(self, tol: float = 1e-8) -> bool:
        """
        Verifies the fundamental DFT group constraint:
        H^T * eta * H = eta
        """
        prod = self.h.T @ self.eta @ self.h
        return bool(np.allclose(prod, self.eta, atol=tol))

    def verify_inverse_relation(self, tol: float = 1e-8) -> bool:
        """
        Verifies that the generalized inverse is given by:
        H^{-1} = eta * H * eta
        """
        h_inv = np.linalg.inv(self.h)
        eta_h_eta = self.eta @ self.h @ self.eta
        return bool(np.allclose(h_inv, eta_h_eta, atol=tol))

    def t_duality_transform(self, direction: int) -> "GeneralizedMetric":
        """
        Performs continuous Buscher / O(D, D; Z) inversion along a specified
        isometry direction (0 <= direction < D).
        """
        if direction < 0 or direction >= self.d:
            raise ValueError(f"Direction {direction} out of bounds for D={self.d}")

        # Construct O(D, D; Z) reflection matrix for this direction
        t_mat = np.eye(2 * self.d, dtype=np.float64)
        # Swap components (direction) and (direction + D)
        i = direction
        j = direction + self.d
        t_mat[i, i] = 0.0
        t_mat[j, j] = 0.0
        t_mat[i, j] = 1.0
        t_mat[j, i] = 1.0

        # Transform generalized metric: H' = T^T * H * T
        h_prime = t_mat.T @ self.h @ t_mat
        
        # Extract new metric g' and B-field B'
        bottom_right_prime = h_prime[self.d:, self.d:]
        g_prime = np.linalg.inv(bottom_right_prime)
        
        bottom_left_prime = h_prime[self.d:, :self.d]
        b_prime = -g_prime @ bottom_left_prime

        # Symmetrize g' and antisymmetrize B' to remove float precision noise
        g_prime = 0.5 * (g_prime + g_prime.T)
        b_prime = 0.5 * (b_prime - b_prime.T)

        return GeneralizedMetric(g=g_prime, b=b_prime)


class FluxChain:
    """
    Tracks the complete string T-duality non-geometric flux chain:
        H_{ijk} ---> f^i_{jk} ---> Q^{ij}_k ---> R^{ijk}
    under successive T-duality operations along compact cycles.
    """
    def __init__(
        self,
        dim: int = 3,
        h_flux: Optional[np.ndarray] = None,
        f_geom: Optional[np.ndarray] = None,
        q_flux: Optional[np.ndarray] = None,
        r_flux: Optional[np.ndarray] = None
    ):
        self.dim = dim
        self.h = h_flux if h_flux is not None else np.zeros((dim, dim, dim), dtype=np.float64)
        self.f = f_geom if f_geom is not None else np.zeros((dim, dim, dim), dtype=np.float64)
        self.q = q_flux if q_flux is not None else np.zeros((dim, dim, dim), dtype=np.float64)
        self.r = r_flux if r_flux is not None else np.zeros((dim, dim, dim), dtype=np.float64)

    @classmethod
    def from_standard_three_torus_h_flux(cls, h_val: float = 1.0) -> "FluxChain":
        """
        Initializes a standard 3-torus background with uniform NS-NS 3-form flux:
        H_{012} = h_val.
        """
        chain = cls(dim=3)
        # Fully antisymmetric Levi-Civita component
        for i, j, k in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
            chain.h[i, j, k] = h_val
        for i, j, k in [(0, 2, 1), (2, 1, 0), (1, 0, 2)]:
            chain.h[i, j, k] = -h_val
        return chain

    def apply_t_duality(self, cycle: int) -> "FluxChain":
        """
        Applies T-duality along direction `cycle` (0, 1, or 2), advancing the flux state:
        - If H_{ijk} has index `cycle`, it transforms into geometric flux f^i_{jk}.
        - If f^i_{jk} has lower index `cycle`, it transforms into non-geometric Q^{ik}_j.
        - If Q^{ij}_k has lower index `cycle`, it transforms into non-local R^{ijk}.
        """
        new_chain = FluxChain(dim=self.dim)
        new_chain.h = self.h.copy()
        new_chain.f = self.f.copy()
        new_chain.q = self.q.copy()
        new_chain.r = self.r.copy()

        c = cycle
        # 1. H_{c, j, k} -> f^c_{j, k}
        for j in range(self.dim):
            for k in range(self.dim):
                if abs(self.h[c, j, k]) > 1e-9:
                    val = self.h[c, j, k]
                    new_chain.f[c, j, k] += val
                    new_chain.h[c, j, k] = 0.0

        # 2. f^i_{c, k} -> Q^{i, c}_k
        for i in range(self.dim):
            for k in range(self.dim):
                if abs(self.f[i, c, k]) > 1e-9:
                    val = self.f[i, c, k]
                    new_chain.q[i, c, k] += val
                    new_chain.f[i, c, k] = 0.0

        # 3. Q^{i, j}_c -> R^{i, j, c}
        for i in range(self.dim):
            for j in range(self.dim):
                if abs(self.q[i, j, c]) > 1e-9:
                    val = self.q[i, j, c]
                    new_chain.r[i, j, c] += val
                    new_chain.q[i, j, c] = 0.0

        return new_chain

    def active_flux_type(self) -> str:
        """Identifies the dominant flux nature of the current background."""
        h_norm = float(np.linalg.norm(self.h))
        f_norm = float(np.linalg.norm(self.f))
        q_norm = float(np.linalg.norm(self.q))
        r_norm = float(np.linalg.norm(self.r))

        if r_norm > 1e-6:
            return "Non-Geometric Non-Local (R-flux)"
        elif q_norm > 1e-6:
            return "Non-Geometric Locally-Geometric (Q-flux T-Fold)"
        elif f_norm > 1e-6:
            return "Geometric Twisted Torus (f-flux)"
        elif h_norm > 1e-6:
            return "Standard Geometric (H-flux)"
        else:
            return "Unfluxed Geometry"


def compute_t_fold_monodromy(q_val: float, radius: float = 1.0) -> np.ndarray:
    """
    Computes the O(2, 2; Z) monodromy matrix M experienced when encircling
    the base S^1 cycle of a non-geometric T-fold.
    """
    # Monodromy element in O(2, 2; Z):
    # M = [[I, beta], [0, I]] where beta^{12} = -beta^{21} = q_val
    m = np.eye(4, dtype=np.float64)
    m[0, 3] = q_val
    m[1, 2] = -q_val
    return m


def dorfman_bracket(
    u: np.ndarray,
    v: np.ndarray,
    structure_constants: Optional[np.ndarray] = None,
    h_flux: Optional[np.ndarray] = None
) -> np.ndarray:
    """
    Computes the Dorfman bracket (Courant algebroid anchor / derived bracket)
    u circ_D v on sections of TM + T*M.
    u, v are 2D-vectors: u = [X, xi], v = [Y, eta].
    """
    u = np.asarray(u, dtype=np.float64)
    v = np.asarray(v, dtype=np.float64)
    tot_dim = len(u)
    d = tot_dim // 2
    x = u[:d]
    xi = u[d:]
    y = v[:d]
    eta = v[d:]

    # Vector part: [X, Y]^k = f^k_{ij} X^i Y^j
    res_vec = np.zeros(d, dtype=np.float64)
    if structure_constants is not None:
        # f[k, i, j]
        res_vec = np.einsum('kij,i,j->k', structure_constants, x, y)

    # 1-form part: L_X eta - i_Y d xi + i_Y i_X H
    res_form = np.zeros(d, dtype=np.float64)
    if structure_constants is not None:
        # L_X eta: (L_X eta)_k = -f^i_{jk} X^j eta_i
        res_form -= np.einsum('ijk,j,i->k', structure_constants, x, eta)
        # i_Y d xi: (i_Y d xi)_k = f^i_{jk} Y^j xi_i
        res_form -= np.einsum('ijk,j,i->k', structure_constants, y, xi)

    if h_flux is not None:
        # (i_Y i_X H)_k = H_{ijk} X^i Y^j
        res_form += np.einsum('ijk,i,j->k', h_flux, x, y)

    return np.concatenate([res_vec, res_form])


def courant_bracket(
    u: np.ndarray,
    v: np.ndarray,
    structure_constants: Optional[np.ndarray] = None,
    h_flux: Optional[np.ndarray] = None
) -> np.ndarray:
    """
    Computes the Courant bracket [u, v]_C on sections of TM + T*M:
        [u, v]_C = 1/2 (u circ_D v - v circ_D u)
    The Courant bracket is manifestly skew-symmetric: [u, v]_C = -[v, u]_C.
    """
    d_uv = dorfman_bracket(u, v, structure_constants, h_flux)
    d_vu = dorfman_bracket(v, u, structure_constants, h_flux)
    return 0.5 * (d_uv - d_vu)

