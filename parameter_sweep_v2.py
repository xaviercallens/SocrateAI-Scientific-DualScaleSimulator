#!/usr/bin/env python3
"""
Parameter Sensitivity Sweep for DualScale Simulator - Version 2
==============================================================

Systematically varies each free parameter (±10%) and measures impact on:
- Observables: w0, wa, screening_factor, boom_action
- Constraints: Cassini bound, DESI consistency, axiom audit passing rate

Reports sensitivity metrics: output change per 1% parameter variation.
"""

import os
import sys
import json
import math
import numpy as np
from typing import Dict, Any, List, Tuple, Callable
from dataclasses import dataclass
import copy

# Add project root to path
sys.path.insert(0, os.path.dirname(__file__))


@dataclass
class SweepPoint:
    """Results at one parameter value."""
    param_name: str
    param_value: float
    w0: float
    wa: float
    screening_factor: float
    boom_action: float
    cassini_satisfied: bool
    desi_consistent: bool
    axiom_audit_passing: bool
    chi2_dual_scale: float
    delta_chi2: float


@dataclass
class ParameterSensitivity:
    """Sensitivity analysis for one parameter."""
    param_name: str
    nominal_value: float
    tested_range: Tuple[float, float]

    # Observable values at three test points
    w0_sweep: List[float]
    wa_sweep: List[float]
    screening_sweep: List[float]
    boom_action_sweep: List[float]

    # Sensitivity: change in output per 1% parameter change
    w0_sensitivity: float
    wa_sensitivity: float
    screening_sensitivity: float
    boom_action_sensitivity: float

    # Constraint status at each point
    cassini_status: List[bool]
    desi_status: List[bool]
    axiom_status: List[bool]

    # Which constraints depend on this parameter
    constrains_cassini: bool
    constrains_desi: bool
    constrains_axiom: bool


def make_potential_function(a_pot: float, b_pot: float) -> Callable:
    """
    Factory for potential function with specified a_pot, b_pot.
    Returns function with signature: compute_potential(x, y) -> (V, dV/dx, dV/dy)
    """
    FRICKE_X = 0.0
    FRICKE_Y = 1.0 / math.sqrt(12.0)
    ORBIFOLD_X = 0.5
    ORBIFOLD_Y = math.sqrt(3.0) / 2.0

    def compute_potential(x: float, y: float) -> Tuple[float, float, float]:
        """Potential with custom a_pot, b_pot."""
        pi = math.pi
        cos_pi_x = math.cos(pi * x)
        sin_pi_x = math.sin(pi * x)
        cos2 = cos_pi_x * cos_pi_x
        sin2 = sin_pi_x * sin_pi_x
        sin2x = math.sin(2.0 * pi * x)

        dyf = y - FRICKE_Y
        dyo = y - ORBIFOLD_Y

        v = a_pot * cos2 + b_pot * sin2 + cos2 * (dyf ** 2) + sin2 * (dyo ** 2)
        dv_dx = -pi * sin2x * (a_pot - b_pot + (dyf ** 2) - (dyo ** 2))
        dv_dy = 2.0 * cos2 * dyf + 2.0 * sin2 * dyo

        return v, dv_dx, dv_dy

    return compute_potential


def run_parametrized_simulation(
    a_pot: float = 1.0,
    b_pot: float = 0.01,
    mu_sym: float = 1.0,
    lambda_sym: float = 1.0,
    pta_suppression: float = 0.005,
    c4_c0_ratio: float = 16.07,
) -> SweepPoint:
    """
    Run simulation with given parameters.
    Returns observables and constraint satisfaction.
    """
    # Import workshopcosmo dynamically to allow patching
    import workshopcosmo as ws

    # Patch the compute_potential function in the module
    original_compute_potential = ws.compute_potential
    ws.compute_potential = make_potential_function(a_pot, b_pot)

    try:
        # Run quintessence simulation
        quint_res = ws.run_quintessence_simulation(t_max=70.0, num_points=500)

        # Run symmetron screening
        sym_res = ws.run_symmetron_screening_simulation(
            r_max=10.0,
            n_points=400,
            mu_sym=mu_sym,
            lambda_sym=lambda_sym,
        )

        # Run observables analysis
        obs_res = ws.run_observables_analysis(quint_res)

        # Run nanograv simulation
        nano_res = ws.run_nanograv_hexadecapole_simulation(c4_c0_ratio=c4_c0_ratio)

        # Run axiom audit
        lean_res = ws.run_lean_verification()

        # Extract observables
        w0 = obs_res.get("w0_fit", -0.5)
        wa = obs_res.get("wa_fit", 0.0)
        screening_factor = sym_res.get("screening_suppression_factor", 1.0)

        # Boom action: use chi2 from observables analysis
        chi2_ds = obs_res.get("chi2_dual_scale", 1.0)
        delta_chi2 = obs_res.get("delta_chi2", 0.0)
        boom_action = abs(chi2_ds) + abs(delta_chi2)

        # Constraint satisfaction
        cassini = sym_res.get("cassini_bound_satisfied", False)
        desi_pref = obs_res.get("prefers_dual_scale", False)
        axiom_ok = lean_res.get("verified_via_axiom_audit", False)

        # Build parameter summary string for identification
        param_id = f"a{a_pot:.2f}_b{b_pot:.3f}_mu{mu_sym:.2f}_l{lambda_sym:.2f}"

        return SweepPoint(
            param_name=param_id,
            param_value=0.0,  # Will be overridden by caller
            w0=float(w0),
            wa=float(wa),
            screening_factor=float(screening_factor),
            boom_action=float(boom_action),
            cassini_satisfied=bool(cassini),
            desi_consistent=bool(desi_pref),
            axiom_audit_passing=bool(axiom_ok),
            chi2_dual_scale=float(chi2_ds),
            delta_chi2=float(delta_chi2),
        )

    finally:
        # Restore original function
        ws.compute_potential = original_compute_potential


def compute_observable_sensitivity(
    values: List[float],
    param_variation_percent: float = 10.0
) -> float:
    """
    Compute sensitivity: change in observable per 1% parameter change.

    Args:
        values: [nominal-10%, nominal, nominal+10%]
        param_variation_percent: 10.0 for ±10% variation

    Returns:
        sensitivity metric (change per 1% parameter variation)
    """
    if len(values) < 3 or any(np.isnan(v) for v in values):
        return 0.0

    v_minus = values[0]
    v_nominal = values[1]
    v_plus = values[2]

    # Output change over full 20% parameter range
    output_change = abs(v_plus - v_minus)

    # Normalize to per 1% change
    # 20% parameter change -> divide by 20
    sensitivity = output_change / (2.0 * param_variation_percent)

    return float(sensitivity)


def sweep_parameter(
    param_name: str,
    nominal_value: float,
    param_variation: float = 10.0,
) -> ParameterSensitivity:
    """
    Sweep one free parameter over ±param_variation% range.
    Measure sensitivity of each observable to this parameter.
    """
    print(f"\n{'='*70}")
    print(f"Sweeping: {param_name}")
    print(f"Nominal: {nominal_value:.6f}")
    print(f"Range: [{nominal_value*(1-param_variation/100):.6f}, {nominal_value*(1+param_variation/100):.6f}]")
    print(f"{'='*70}")

    # Test points
    frac = param_variation / 100.0
    v_minus = nominal_value * (1.0 - frac)
    v_nominal = nominal_value
    v_plus = nominal_value * (1.0 + frac)
    test_points = [v_minus, v_nominal, v_plus]

    # Default parameters
    base_params = {
        "a_pot": 1.0,
        "b_pot": 0.01,
        "mu_sym": 1.0,
        "lambda_sym": 1.0,
        "pta_suppression": 0.005,
        "c4_c0_ratio": 16.07,
    }

    # Run at each test point
    w0_values = []
    wa_values = []
    screening_values = []
    boom_values = []
    cassini_results = []
    desi_results = []
    axiom_results = []

    for i, test_val in enumerate(test_points):
        print(f"\n  Point {i+1}/3: {param_name} = {test_val:.6f}")

        params = copy.copy(base_params)
        params[param_name] = test_val

        try:
            result = run_parametrized_simulation(**params)

            w0_values.append(result.w0)
            wa_values.append(result.wa)
            screening_values.append(result.screening_factor)
            boom_values.append(result.boom_action)
            cassini_results.append(result.cassini_satisfied)
            desi_results.append(result.desi_consistent)
            axiom_results.append(result.axiom_audit_passing)

            print(f"    w0={result.w0:.4f}, wa={result.wa:.4f}")
            print(f"    screening={result.screening_factor:.6f}, boom={result.boom_action:.4f}")
            print(f"    Cassini={result.cassini_satisfied}, DESI={result.desi_consistent}, Axiom={result.axiom_audit_passing}")

        except Exception as e:
            print(f"    ERROR: {str(e)[:100]}")
            # Use NaN placeholders
            w0_values.append(np.nan)
            wa_values.append(np.nan)
            screening_values.append(np.nan)
            boom_values.append(np.nan)
            cassini_results.append(False)
            desi_results.append(False)
            axiom_results.append(False)

    # Compute sensitivities
    w0_sens = compute_observable_sensitivity(w0_values, param_variation)
    wa_sens = compute_observable_sensitivity(wa_values, param_variation)
    screen_sens = compute_observable_sensitivity(screening_values, param_variation)
    boom_sens = compute_observable_sensitivity(boom_values, param_variation)

    # Detect constraint dependence
    cassini_dependent = any(c1 != c2 for c1, c2 in zip(cassini_results[:-1], cassini_results[1:]))
    desi_dependent = any(d1 != d2 for d1, d2 in zip(desi_results[:-1], desi_results[1:]))
    axiom_dependent = any(a1 != a2 for a1, a2 in zip(axiom_results[:-1], axiom_results[1:]))

    print(f"\n  Sensitivities:")
    print(f"    w0:        {w0_sens:.6f} per 1% change")
    print(f"    wa:        {wa_sens:.6f} per 1% change")
    print(f"    screening: {screen_sens:.6f} per 1% change")
    print(f"    boom:      {boom_sens:.6f} per 1% change")
    print(f"  Constraint dependence:")
    print(f"    Cassini: {cassini_dependent}, DESI: {desi_dependent}, Axiom: {axiom_dependent}")

    return ParameterSensitivity(
        param_name=param_name,
        nominal_value=nominal_value,
        tested_range=(v_minus, v_plus),
        w0_sweep=w0_values,
        wa_sweep=wa_values,
        screening_sweep=screening_values,
        boom_action_sweep=boom_values,
        w0_sensitivity=w0_sens,
        wa_sensitivity=wa_sens,
        screening_sensitivity=screen_sens,
        boom_action_sensitivity=boom_sens,
        cassini_status=cassini_results,
        desi_status=desi_results,
        axiom_status=axiom_results,
        constrains_cassini=cassini_dependent,
        constrains_desi=desi_dependent,
        constrains_axiom=axiom_dependent,
    )


def main():
    """Run full parameter sweep."""
    print("\n" + "="*70)
    print("PARAMETER SENSITIVITY SWEEP - DualScale Cosmological Simulator")
    print("="*70)

    # Free parameters and nominal values
    free_params = [
        ("a_pot", 1.0),
        ("b_pot", 0.01),
        ("mu_sym", 1.0),
        ("lambda_sym", 1.0),
        ("pta_suppression", 0.005),
        ("c4_c0_ratio", 16.07),
    ]

    results = {}
    for param_name, nominal_val in free_params:
        result = sweep_parameter(param_name, nominal_val, param_variation=10.0)
        results[param_name] = result

    # Identify most constraining parameters
    most_sensitive_cassini = max(
        (p for p in results.values() if p.constrains_cassini),
        key=lambda r: r.screening_sensitivity,
        default=None
    )
    most_sensitive_desi = max(
        (p for p in results.values() if p.constrains_desi),
        key=lambda r: (r.w0_sensitivity + r.wa_sensitivity) / 2.0,
        default=None
    )
    most_sensitive_axiom = max(
        (p for p in results.values() if p.constrains_axiom),
        key=lambda r: r.screening_sensitivity,
        default=None
    )

    # Save results
    output_file = "parameter_sweep_results.json"
    json_results = {}
    for param_name, sens in results.items():
        json_results[param_name] = {
            "nominal_value": float(sens.nominal_value),
            "tested_range": [float(sens.tested_range[0]), float(sens.tested_range[1])],
            "w0_sensitivity": float(sens.w0_sensitivity),
            "wa_sensitivity": float(sens.wa_sensitivity),
            "screening_sensitivity": float(sens.screening_sensitivity),
            "boom_action_sensitivity": float(sens.boom_action_sensitivity),
            "constrains_cassini": bool(sens.constrains_cassini),
            "constrains_desi": bool(sens.constrains_desi),
            "constrains_axiom": bool(sens.constrains_axiom),
            "cassini_status": [bool(c) for c in sens.cassini_status],
            "desi_status": [bool(d) for d in sens.desi_status],
            "axiom_status": [bool(a) for a in sens.axiom_status],
        }

    json_results["summary"] = {
        "most_constraining_cassini": most_sensitive_cassini.param_name if most_sensitive_cassini else "none",
        "most_constraining_desi": most_sensitive_desi.param_name if most_sensitive_desi else "none",
        "most_constraining_axiom": most_sensitive_axiom.param_name if most_sensitive_axiom else "none",
    }

    with open(output_file, "w") as f:
        json.dump(json_results, f, indent=2)

    print(f"\nResults saved to {output_file}")
    return results, json_results


if __name__ == "__main__":
    results_objs, json_results = main()
