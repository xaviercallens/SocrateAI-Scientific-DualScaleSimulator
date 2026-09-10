"""
=============================================================================
LeanFlow Datasets: Kreuzer-Skarke (KS) 4D Reflexive Polyhedra Ingestion
=============================================================================
Parses Kreuzer-Skarke reflexive polyhedra definitions for Calabi-Yau threefolds,
extracts Hodge numbers (h11, h21), computes Euler characteristic chi, and derives
the D-brane tadpole bound.
=============================================================================
"""

from dataclasses import dataclass
from typing import List, Dict, Any, Optional
import numpy as np


@dataclass
class KreuzerSkarkePolytope:
    ks_id: str
    vertices: np.ndarray             # Shape (V, 4) in Z^4
    h11: int
    h21: int
    euler_characteristic: int       # chi = 2 * (h11 - h21)
    picard_number: int              # rho = h11
    is_reflexive: bool
    tadpole_bound: float            # Q_D3 <= |chi| / 24

    def summary(self) -> Dict[str, Any]:
        return {
            "ks_id": self.ks_id,
            "num_vertices": self.vertices.shape[0],
            "h11": self.h11,
            "h21": self.h21,
            "chi": self.euler_characteristic,
            "picard_number": self.picard_number,
            "is_reflexive": self.is_reflexive,
            "tadpole_bound": self.tadpole_bound,
        }


def parse_ks_polytope(
    vertices: List[List[int]],
    h11: int,
    h21: int,
    ks_id: str = "KS_Custom",
    is_reflexive: bool = True
) -> KreuzerSkarkePolytope:
    """
    Parses a 4D lattice polyhedron defining a toric Calabi-Yau threefold hypersurface.
    """
    V = np.array(vertices, dtype=int)
    assert V.shape[1] == 4, f"Expected 4D vertices, got shape {V.shape}"

    chi = 2 * (h11 - h21)
    tadpole_bound = abs(chi) / 24.0

    return KreuzerSkarkePolytope(
        ks_id=ks_id,
        vertices=V,
        h11=h11,
        h21=h21,
        euler_characteristic=chi,
        picard_number=h11,
        is_reflexive=is_reflexive,
        tadpole_bound=tadpole_bound
    )
