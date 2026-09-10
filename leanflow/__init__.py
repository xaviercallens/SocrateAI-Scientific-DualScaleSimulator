"""
=============================================================================
LeanFlow: Guaranteed-Invariant Scientific Machine Learning & String Toolkit
=============================================================================
Phase 2 Endorsement Toolkit for Computational String Theory & SciML:
- Python & PyTorch Interoperable Stiff ODE/PDE Solver
- String Theory Databases (Kreuzer-Skarke, CICY)
- Declarative Invariant Specification DSL (@leanflow.guardrail)
- Working Group Toolkits (cymetric Kahler stabilization, Swampland MCMC filter)
- Asynchronous Lean 4 Proof Kernel Verification Bridge
=============================================================================
"""

__version__ = "2.0.0"

from leanflow.core.solver import Solver, SolverResult, SolverTelemetry
from leanflow.core.dsl import guardrail, InvariantRegistry, GuardrailTelemetry
from leanflow.core.projections import (
    metric_positivity_projection,
    modular_domain_fold,
    weak_energy_condition_lock,
    kahler_cone_positivity_projection,
    tadpole_budget_check,
)
from leanflow.datasets.standard_manifolds import (
    load_cicy,
    load_ks,
    list_available_manifolds,
)
from leanflow.datasets.cicy import CICYConfiguration, parse_cicy_matrix
from leanflow.datasets.kreuzer_skarke import KreuzerSkarkePolytope, parse_ks_polytope
from leanflow.bridge.lean_ipc import LeanVerificationClient
from leanflow.working_groups.cymetric_kahler import (
    project_kahler_metric,
    compute_monge_ampere_loss,
)
from leanflow.working_groups.swampland_mcmc import SwamplandMCMCFilter
from leanflow.working_groups.topological_defects import (
    TopologicalDefectExtractor,
    DefectSummary,
)
from leanflow.core.supergravity import (
    SupergravityModuli,
    SupergravityParameters,
    evaluate_sugra_landscape_grid,
)
from leanflow.core.generalized_geometry import (
    GeneralizedMetric,
    FluxChain,
    compute_t_fold_monodromy,
)
from leanflow.core.dbrane_inflow import (
    DpBraneBoundaryState,
    CallanHarveyInflow,
)


def extract_topological_defects(
    phi_grid,
    energy_threshold: float = 0.65
) -> DefectSummary:
    """
    Convenience function to extract topological defect 1-skeletons from a scalar field grid.
    """
    extractor = TopologicalDefectExtractor(energy_threshold=energy_threshold)
    return extractor.extract_from_grid(phi_grid)


__all__ = [
    "__version__",
    "Solver",
    "SolverResult",
    "SolverTelemetry",
    "guardrail",
    "InvariantRegistry",
    "GuardrailTelemetry",
    "metric_positivity_projection",
    "modular_domain_fold",
    "weak_energy_condition_lock",
    "kahler_cone_positivity_projection",
    "tadpole_budget_check",
    "load_cicy",
    "load_ks",
    "list_available_manifolds",
    "CICYConfiguration",
    "parse_cicy_matrix",
    "KreuzerSkarkePolytope",
    "parse_ks_polytope",
    "LeanVerificationClient",
    "project_kahler_metric",
    "compute_monge_ampere_loss",
    "SwamplandMCMCFilter",
    "TopologicalDefectExtractor",
    "DefectSummary",
    "extract_topological_defects",
    "SupergravityModuli",
    "SupergravityParameters",
    "evaluate_sugra_landscape_grid",
    "GeneralizedMetric",
    "FluxChain",
    "compute_t_fold_monodromy",
    "DpBraneBoundaryState",
    "CallanHarveyInflow",
]
