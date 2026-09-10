"""
=============================================================================
LeanFlow Datasets: Standard String Theory Manifolds Repository
=============================================================================
Provides canonical, pre-indexed Calabi-Yau threefolds and fourfolds used by
computational string theory groups (CICY, Kreuzer-Skarke, K3 x T2).
=============================================================================
"""

from typing import Dict, List, Optional
from leanflow.datasets.cicy import CICYConfiguration, parse_cicy_matrix
from leanflow.datasets.kreuzer_skarke import KreuzerSkarkePolytope, parse_ks_polytope


# 1. Canonical CICY Manifolds
_PRESET_CICY: Dict[str, CICYConfiguration] = {
    "quintic": parse_cicy_matrix(
        matrix_data=[[5]],
        ambient_spaces=[4],
        name="Quintic_P4[5]",
        h11=1,
        h21=101
    ),
    "cicy_7887": parse_cicy_matrix(
        matrix_data=[
            [1, 1, 0],
            [1, 0, 1],
            [0, 1, 1],
            [0, 0, 2]
        ],
        ambient_spaces=[1, 1, 1, 1],
        name="CICY_7887",
        h11=4,
        h21=68
    ),
    "tian_yau": parse_cicy_matrix(
        matrix_data=[
            [3, 0, 1],
            [0, 3, 1]
        ],
        ambient_spaces=[3, 3],
        name="Tian_Yau_P3xP3",
        h11=14,
        h21=23
    ),
    "bicubic": parse_cicy_matrix(
        matrix_data=[
            [3],
            [3]
        ],
        ambient_spaces=[2, 2],
        name="Bicubic_P2xP2",
        h11=2,
        h21=83
    ),
    "k3_x_t2": parse_cicy_matrix(
        matrix_data=[
            [4, 0],
            [0, 2]
        ],
        ambient_spaces=[3, 1],
        name="K3_x_T2_Fibration",
        h11=3,
        h21=3
    )
}

# 2. Canonical Kreuzer-Skarke Polyhedra
_PRESET_KS: Dict[str, KreuzerSkarkePolytope] = {
    "ks_quintic": parse_ks_polytope(
        vertices=[
            [1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
            [-1, -1, -1, -1]
        ],
        h11=1,
        h21=101,
        ks_id="KS_Quintic_4D",
        is_reflexive=True
    ),
    "ks_bicubic": parse_ks_polytope(
        vertices=[
            [1, 0, 0, 0],
            [-1, -1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, -1, -1],
            [0, 1, 0, 1]
        ],
        h11=2,
        h21=83,
        ks_id="KS_Bicubic_4D",
        is_reflexive=True
    ),
    "ks_k3_surface": parse_ks_polytope(
        vertices=[
            [1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [-1, -1, -1, 0]
        ],
        h11=20,
        h21=0,
        ks_id="KS_K3_Orbifold",
        is_reflexive=True
    )
}


def load_cicy(identifier: str) -> CICYConfiguration:
    """
    Loads a canonical CICY Calabi-Yau geometry by name or ID.
    Available: 'quintic', 'cicy_7887', 'tian_yau', 'bicubic', 'k3_x_t2'.
    """
    key = identifier.lower().strip()
    if key in _PRESET_CICY:
        return _PRESET_CICY[key]
    raise ValueError(f"Unknown CICY identifier '{identifier}'. Available: {list(_PRESET_CICY.keys())}")


def load_ks(identifier: str) -> KreuzerSkarkePolytope:
    """
    Loads a canonical Kreuzer-Skarke 4D reflexive polytope.
    Available: 'ks_quintic', 'ks_bicubic', 'ks_k3_surface'.
    """
    key = identifier.lower().strip()
    if key in _PRESET_KS:
        return _PRESET_KS[key]
    raise ValueError(f"Unknown KS identifier '{identifier}'. Available: {list(_PRESET_KS.keys())}")


def list_available_manifolds() -> Dict[str, List[str]]:
    """Returns a dictionary listing all pre-indexed Calabi-Yau geometries."""
    return {
        "cicy": list(_PRESET_CICY.keys()),
        "kreuzer_skarke": list(_PRESET_KS.keys())
    }
