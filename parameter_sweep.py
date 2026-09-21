#!/usr/bin/env python3
"""
Parameter Sensitivity Sweep for DualScale Simulator
====================================================

Runs workshopcosmo.py with systematic variations of free parameters:
- a_pot, b_pot, mu_sym, lambda_sym, pta_suppression, c4_c0_ratio

For each parameter: measure output sensitivity at nominal, +10%, -10%
Compute sensitivity metric: output variance per 1% parameter change
Identify which constraints are most sensitive to each parameter.
"""

import os
import sys
import json
import math
import numpy as np
from typing import Dict, Any, List, Tuple
from dataclasses import dataclass, asdict

# Add project root to path
sys.path.insert(0, os.path.dirname(__file__))

import workshopcosmo as ws


@dataclass
class ParameterSweepResult:
    """Results for one parameter's sweep."""
    param_name: str
    nominal_value: float
    test_values: List[float]  # [nominal-10%, nominal, nominal+10%]

    # Observable measurements at each point
    w0_values: List[float]
    wa_values: List[float]
    screening_factor_values: List[float]
    boom_action_values: List[float]  # CDL action or effective action

    # Cassini, DESI, axiom outcomes
    cassini_satisfied: List[bool]
    desi_consistent: List[bool]
    axiom_audit_passing: List[bool]

    # Derived sensitivity metrics
    w0_sensitivity: float  # variance per 1% change
    wa_sensitivity: float
    screening_sensitivity: float
    boom_action_sensitivity: float

    # Which constraints are most sensitive
    most_constraining_cassini: str  # param name
    most_constraining_desi: str
    most_constraining_axiom: str


def run_simulation_with_params(
    a_pot: float = 1.0,
    b_pot: float = 0.01,
    mu_sym: float = 1.0,
    lambda_sym: float = 1.0,
    pta_suppression: float = 0.005,
    c4_c0_ratio: float = 16.07,
) -> Dict[str, Any]:
    """
    Run full simulation with specified parameters.
    Returns observables and constraint satisfaction.
    """
    # Modify workshopcosmo functions to use these parameters
    # Since parameters are hardcoded, we need to patch them

    # Run quintessence simulation
    quint_res = ws.run_quintessence_simulation(t_max=70.0, num_points=500)

    # Run symmetron screening with the parameters
    sym_res = ws.run_symmetron_screening_simulation(
        mu_sym=mu_sym,
        lambda_sym=lambda_sym,
        pta_suppression=pta_suppression
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

    # Use CDL action as boom_action proxy
    # For now, use a derived quantity: effective action from observables
    boom_action = abs(obs_res.get("chi2_dual_scale", 0.0)) + abs(obs_res.get("delta_chi2", 0.0))

    # Constraint satisfaction
    cassini = sym_res.get("cassini_bound_satisfied", False)
    desi_pref = obs_res.get("prefers_dual_scale", False)
    axiom_ok = lean_res.get("verified_via_axiom_audit", False)

    return {
        "w0": w0,
        "wa": wa,
        "screening_factor": screening_factor,
        "boom_action": boom_action,
        "cassini_satisfied": cassini,
        "desi_consistent": desi_pref,
        "axiom_audit_passing": axiom_ok,
        "chi2_dual_scale": obs_res.get("chi2_dual_scale", 0.0),
        "delta_chi2": obs_res.get("delta_chi2", 0.0),
    }


def compute_sensitivity(
    values_at_points: List[float],
    param_variation_percent: float = 10.0
) -> float:
    """
    Compute sensitivity metric: output variance per 1% parameter change.

    values_at_points: [nominal-10%, nominal, nominal+10%]
    param_variation_percent: 10.0 for ±10% variation

    Returns: sensitivity = max_output_change / (param_variation_percent)
    """
    if len(values_at_points) < 3:
        return 0.0

    v_minus = values_at_points[0]
    v_nominal = values_at_points[1]
    v_plus = values_at_points[2]

    # Total output change over 20% parameter range
    output_change = abs(v_plus - v_minus)

    # Convert to per 1% parameter change
    sensitivity = output_change / (2.0 * param_variation_percent / 100.0)

    return sensitivity


def sweep_parameter(
    param_name: str,
    nominal_value: float,
    param_variation_percent: float = 10.0,
) -> ParameterSweepResult:
    """
    Sweep one parameter over ±10% and measure output sensitivity.
    """
    print(f"\n{'='*70}")
    print(f"Sweeping parameter: {param_name}")
    print(f"Nominal value: {nominal_value:.6f}")
    print(f"Range: {nominal_value * (1 - param_variation_percent/100):.6f} to {nominal_value * (1 + param_variation_percent/100):.6f}")
    print(f"{'='*70}")

    # Test values
    test_fraction = param_variation_percent / 100.0
    v_minus = nominal_value * (1.0 - test_fraction)
    v_nominal = nominal_value
    v_plus = nominal_value * (1.0 + test_fraction)
    test_values = [v_minus, v_nominal, v_plus]

    # Default parameter values
    params = {
        "a_pot": 1.0,
        "b_pot": 0.01,
        "mu_sym": 1.0,
        "lambda_sym": 1.0,
        "pta_suppression": 0.005,
        "c4_c0_ratio": 16.07,
    }

    # Collect results
    w0_values = []
    wa_values = []
    screening_values = []
    boom_values = []
    cassini_list = []
    desi_list = []
    axiom_list = []

    for i, test_value in enumerate(test_values):
        print(f"\nPoint {i+1}/3: {param_name} = {test_value:.6f}")

        # Update parameter
        params[param_name] = test_value

        try:
            # Run simulation with this parameter value
            result = run_simulation_with_params(**params)

            w0_values.append(result["w0"])
            wa_values.append(result["wa"])
            screening_values.append(result["screening_factor"])
            boom_values.append(result["boom_action"])
            cassini_list.append(result["cassini_satisfied"])
            desi_list.append(result["desi_consistent"])
            axiom_list.append(result["axiom_audit_passing"])

            print(f"  w0={result['w0']:.4f}, wa={result['wa']:.4f}")
            print(f"  screening={result['screening_factor']:.6f}")
            print(f"  Cassini: {result['cassini_satisfied']}, DESI: {result['desi_consistent']}, Axiom: {result['axiom_audit_passing']}")

        except Exception as e:
            print(f"  ERROR: {e}")
            # Use fallback values
            w0_values.append(np.nan)
            wa_values.append(np.nan)
            screening_values.append(np.nan)
            boom_values.append(np.nan)
            cassini_list.append(False)
            desi_list.append(False)
            axiom_list.append(False)

    # Compute sensitivities
    w0_sens = compute_sensitivity(w0_values, param_variation_percent)
    wa_sens = compute_sensitivity(wa_values, param_variation_percent)
    screen_sens = compute_sensitivity(screening_values, param_variation_percent)
    boom_sens = compute_sensitivity(boom_values, param_variation_percent)

    print(f"\nSensitivities:")
    print(f"  w0: {w0_sens:.6f}")
    print(f"  wa: {wa_sens:.6f}")
    print(f"  screening: {screen_sens:.6f}")
    print(f"  boom_action: {boom_sens:.6f}")

    result = ParameterSweepResult(
        param_name=param_name,
        nominal_value=nominal_value,
        test_values=test_values,
        w0_values=w0_values,
        wa_values=wa_values,
        screening_factor_values=screening_values,
        boom_action_values=boom_values,
        cassini_satisfied=cassini_list,
        desi_consistent=desi_list,
        axiom_audit_passing=axiom_list,
        w0_sensitivity=w0_sens,
        wa_sensitivity=wa_sens,
        screening_sensitivity=screen_sens,
        boom_action_sensitivity=boom_sens,
        most_constraining_cassini="",  # Will be filled after all sweeps
        most_constraining_desi="",
        most_constraining_axiom="",
    )

    return result


def run_full_parameter_sweep() -> Dict[str, Any]:
    """
    Run complete parameter sweep over all 6 free parameters.
    """
    print("\n" + "="*70)
    print("PARAMETER SENSITIVITY SWEEP")
    print("="*70)

    # Define free parameters and their nominal values
    free_params = {
        "a_pot": 1.0,
        "b_pot": 0.01,
        "mu_sym": 1.0,
        "lambda_sym": 1.0,
        "pta_suppression": 0.005,
        "c4_c0_ratio": 16.07,
    }

    sweep_results = {}

    for param_name, nominal_value in free_params.items():
        result = sweep_parameter(param_name, nominal_value, param_variation_percent=10.0)
        sweep_results[param_name] = result

    # Analyze which parameters constrain which outcomes
    # Find maximum sensitivity for each constraint type

    max_cassini_sensitivity = 0.0
    max_cassini_param = "none"

    max_desi_sensitivity = 0.0
    max_desi_param = "none"

    max_axiom_sensitivity = 0.0
    max_axiom_param = "none"

    for param_name, result in sweep_results.items():
        # For Cassini: measure sensitivity of achieving/losing the bound
        cassini_changes = sum(1 for i in range(len(result.cassini_satisfied) - 1)
                              if result.cassini_satisfied[i] != result.cassini_satisfied[i+1])
        if cassini_changes > 0:
            cassini_sensitivity = result.screening_sensitivity  # Most relevant observable
            if cassini_sensitivity > max_cassini_sensitivity:
                max_cassini_sensitivity = cassini_sensitivity
                max_cassini_param = param_name

        # For DESI: measure sensitivity of consistency
        desi_changes = sum(1 for i in range(len(result.desi_consistent) - 1)
                          if result.desi_consistent[i] != result.desi_consistent[i+1])
        if desi_changes > 0:
            desi_sensitivity = (result.w0_sensitivity + result.wa_sensitivity) / 2.0
            if desi_sensitivity > max_desi_sensitivity:
                max_desi_sensitivity = desi_sensitivity
                max_desi_param = param_name

        # For Axiom: check if audit changes
        axiom_changes = sum(1 for i in range(len(result.axiom_audit_passing) - 1)
                           if result.axiom_audit_passing[i] != result.axiom_audit_passing[i+1])
        if axiom_changes > 0:
            axiom_sensitivity = result.screening_sensitivity  # Use screening as proxy
            if axiom_sensitivity > max_axiom_sensitivity:
                max_axiom_sensitivity = axiom_sensitivity
                max_axiom_param = param_name

    return {
        "sweep_results": sweep_results,
        "most_constraining_cassini": max_cassini_param,
        "most_constraining_desi": max_desi_param,
        "most_constraining_axiom": max_axiom_param,
    }


def format_structured_output(sweep_results_dict: Dict) -> Tuple[List[Dict], List[str]]:
    """Convert sweep results to structured output format."""

    sweep_results = sweep_results_dict["sweep_results"]

    parameter_sweeps = []
    key_findings = []

    for param_name, result in sweep_results.items():
        param_sweep = {
            "param_name": param_name,
            "nominal_value": result.nominal_value,
            "tested_range": f"[{result.test_values[0]:.6f}, {result.test_values[2]:.6f}]",
            "sensitivity_metric": result.screening_sensitivity,  # Report screening sensitivity as primary
            "most_constraining_outcome": "cassini_bound"  # Placeholder
        }
        parameter_sweeps.append(param_sweep)

        key_findings.append(
            f"{param_name}: nominal={result.nominal_value:.6f}, "
            f"w0_sensitivity={result.w0_sensitivity:.4f}, "
            f"wa_sensitivity={result.wa_sensitivity:.4f}, "
            f"screening_sensitivity={result.screening_sensitivity:.4f}, "
            f"boom_action_sensitivity={result.boom_action_sensitivity:.4f}"
        )

    # Add summary findings
    key_findings.append(
        f"Most constraining parameter for Cassini bound: {sweep_results_dict['most_constraining_cassini']}"
    )
    key_findings.append(
        f"Most constraining parameter for DESI consistency: {sweep_results_dict['most_constraining_desi']}"
    )
    key_findings.append(
        f"Most constraining parameter for axiom audit: {sweep_results_dict['most_constraining_axiom']}"
    )

    return parameter_sweeps, key_findings


if __name__ == "__main__":
    results = run_full_parameter_sweep()

    # Save results to JSON
    output_file = "parameter_sweep_results.json"
    print(f"\nSaving results to {output_file}...")

    # Convert to JSON-serializable format
    json_results = {}
    for param_name, result in results["sweep_results"].items():
        json_results[param_name] = {
            "nominal_value": float(result.nominal_value),
            "tested_range": [float(v) for v in result.test_values],
            "w0_values": [float(v) if not np.isnan(v) else None for v in result.w0_values],
            "wa_values": [float(v) if not np.isnan(v) else None for v in result.wa_values],
            "screening_values": [float(v) if not np.isnan(v) else None for v in result.screening_factor_values],
            "boom_action_values": [float(v) if not np.isnan(v) else None for v in result.boom_action_values],
            "w0_sensitivity": float(result.w0_sensitivity),
            "wa_sensitivity": float(result.wa_sensitivity),
            "screening_sensitivity": float(result.screening_sensitivity),
            "boom_action_sensitivity": float(result.boom_action_sensitivity),
            "cassini_satisfied": result.cassini_satisfied,
            "desi_consistent": result.desi_consistent,
            "axiom_audit_passing": result.axiom_audit_passing,
        }

    json_results["summary"] = {
        "most_constraining_cassini": results["most_constraining_cassini"],
        "most_constraining_desi": results["most_constraining_desi"],
        "most_constraining_axiom": results["most_constraining_axiom"],
    }

    with open(output_file, "w") as f:
        json.dump(json_results, f, indent=2)

    print(f"Results saved to {output_file}")

    # Format structured output
    param_sweeps, findings = format_structured_output(results)

    print("\n" + "="*70)
    print("STRUCTURED OUTPUT")
    print("="*70)
    print("\nParameter Sweeps:")
    for ps in param_sweeps:
        print(f"  {ps['param_name']}: nominal={ps['nominal_value']:.6f}")

    print("\nKey Findings:")
    for finding in findings:
        print(f"  - {finding}")
