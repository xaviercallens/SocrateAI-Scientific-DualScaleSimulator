"""
LeanFlow Core Module
"""
from leanflow.core.projections import (
    metric_positivity_projection,
    modular_domain_fold,
    weak_energy_condition_lock,
    kahler_cone_positivity_projection,
    tadpole_budget_check,
)
from leanflow.core.dsl import guardrail, InvariantRegistry, GuardrailTelemetry
from leanflow.core.solver import Solver, SolverResult, SolverTelemetry

__all__ = [
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
]
