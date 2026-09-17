#!/usr/bin/env python3
"""
Tests for cross-consistency checker.

Positive control: consistent fixtures → exit 0
Negative control: one mutated number → exit 1 naming the quantity
"""

import json
import tempfile
import shutil
import subprocess
import sys
import os
from pathlib import Path


def create_consistent_fixture(tmpdir: str) -> None:
    """Create a fixture with consistent values across all files."""

    # Create necessary subdirectories
    os.makedirs(f"{tmpdir}/audit", exist_ok=True)
    os.makedirs(f"{tmpdir}/papers/T-dulaity alone", exist_ok=True)
    os.makedirs(f"{tmpdir}/scripts", exist_ok=True)

    # Create simulation_results.json
    sim_results = {
        "axe3_formal_proof": {
            "theorems_verified_count": 9
        },
        "axe4_moonshine_voa": {
            "theorems_verified_count": 16
        },
        "axe6_kummer_tda_langevin": {
            "mapper_nodes": 187,
            "mapper_edges": 557,
            "final_string_count": 101,
            "theorems_verified_count": 11
        },
        "observables": {
            "w0": -0.9999963005031299,
            "wa": -5.46990569108541e-06,
            "delta_chi2_vs_lcdm": -0.0004271558779258555,
            "ln_bayes_factor": 0.00021357793896292776
        }
    }
    with open(f"{tmpdir}/simulation_results.json", "w") as f:
        json.dump(sim_results, f, indent=2)

    # Create tachyon_condensation_summary.json
    tachyon = {
        "soliton_width": 1.5398970458808714
    }
    with open(f"{tmpdir}/tachyon_condensation_summary.json", "w") as f:
        json.dump(tachyon, f, indent=2)

    # Create swampland_geodesic_summary.json
    swampland = {
        "final_mass_gap": 0.007157187468270595
    }
    with open(f"{tmpdir}/swampland_geodesic_summary.json", "w") as f:
        json.dump(swampland, f, indent=2)

    # Create LaTeX file with matching values
    latex_content = r"""
\documentclass{article}
\begin{document}

% Mapper nodes and edges: 187 nodes, 557 edges, $\beta_1 = 376$
(187 nodes, 557 edges, $\beta_1 = 376$ persistent loops)

% Betti number
$\beta_1 = 376$

% Screening factor: 1.24 × 10^{-4}
$\Delta R/R \approx 1.24 \times 10^{-4}$

% Gamma - 1: 2.48 × 10^{-6}
$|\gamma - 1| \approx 2.48 \times 10^{-6}$

% Bubble radius
$R_c \approx 2.82$

% Bounce action
$S_E = 137.15$

\end{document}
"""
    with open(f"{tmpdir}/papers/T-dulaity alone/T_duality_Alone.tex", "w") as f:
        f.write(latex_content)

    # Create parameters.json for the fixture
    params = {
        "quantities": {
            "mapper_nodes": {
                "description": "Number of nodes in TDA Mapper 1-skeleton graph",
                "canonical_source": {
                    "type": "json_file",
                    "path": "simulation_results.json",
                    "key_path": "axe6_kummer_tda_langevin.mapper_nodes"
                },
                "canonical_value": 187,
                "comparison": "exact_integer",
                "occurrences": [
                    {
                        "file": "papers/T-dulaity alone/T_duality_Alone.tex",
                        "regex": r"\(([0-9]+)\s+nodes.*edges.*\$\\beta_1"
                    }
                ]
            },
            "mapper_edges": {
                "description": "Number of edges in TDA Mapper 1-skeleton graph",
                "canonical_source": {
                    "type": "json_file",
                    "path": "simulation_results.json",
                    "key_path": "axe6_kummer_tda_langevin.mapper_edges"
                },
                "canonical_value": 557,
                "comparison": "exact_integer",
                "occurrences": [
                    {
                        "file": "papers/T-dulaity alone/T_duality_Alone.tex",
                        "regex": r"nodes,\s+([0-9]+)\s+edges.*\$\\beta_1"
                    }
                ]
            },
            "beta_1": {
                "description": "First Betti number",
                "canonical_source": {
                    "type": "computed",
                    "formula": "|E| - |V| + |C|",
                    "expected_value": 376
                },
                "canonical_value": 376,
                "comparison": "exact_integer",
                "occurrences": [
                    {
                        "file": "papers/T-dulaity alone/T_duality_Alone.tex",
                        "regex": r"\$\\beta_1\s*=\s*([0-9]+)"
                    }
                ]
            },
            "final_string_count": {
                "description": "Number of topological cosmic strings",
                "canonical_source": {
                    "type": "json_file",
                    "path": "simulation_results.json",
                    "key_path": "axe6_kummer_tda_langevin.final_string_count"
                },
                "canonical_value": 101,
                "comparison": "exact_integer",
                "occurrences": []
            },
            "soliton_width": {
                "description": "Width of Sen soliton",
                "canonical_source": {
                    "type": "json_file",
                    "path": "tachyon_condensation_summary.json",
                    "key_path": "soliton_width"
                },
                "canonical_value": 1.5398970458808714,
                "comparison": "relative_tolerance",
                "tolerance": 1e-9,
                "occurrences": []
            },
            "sdl_final_mass_gap": {
                "description": "Final SDC mass gap",
                "canonical_source": {
                    "type": "json_file",
                    "path": "swampland_geodesic_summary.json",
                    "key_path": "final_mass_gap"
                },
                "canonical_value": 0.007157187468270595,
                "comparison": "relative_tolerance",
                "tolerance": 1e-9,
                "occurrences": []
            },
            "screening_factor": {
                "description": "Symmetron screening factor",
                "canonical_source": {
                    "type": "latex_constant",
                    "source": "papers/T-dulaity alone/T_duality_Alone.tex"
                },
                "canonical_value": 1.24e-4,
                "comparison": "relative_tolerance",
                "tolerance": 1e-9,
                "occurrences": [
                    {
                        "file": "papers/T-dulaity alone/T_duality_Alone.tex",
                        "regex": r"≈\s*([0-9.]+)\s*\\times\s*10\^\\{-4\\}"
                    }
                ]
            },
            "gamma_minus_1": {
                "description": "PPN parameter deviation",
                "canonical_source": {
                    "type": "latex_constant",
                    "source": "papers/T-dulaity alone/T_duality_Alone.tex"
                },
                "canonical_value": 2.48e-6,
                "comparison": "relative_tolerance",
                "tolerance": 1e-9,
                "occurrences": [
                    {
                        "file": "papers/T-dulaity alone/T_duality_Alone.tex",
                        "regex": r"≈\s*([0-9.]+)\s*\\times\s*10\^\\{-6\\}"
                    }
                ]
            },
            "w0": {
                "description": "w0 parameter",
                "canonical_source": {
                    "type": "json_file",
                    "path": "simulation_results.json",
                    "key_path": "observables.w0"
                },
                "canonical_value": -0.9999963005031299,
                "comparison": "relative_tolerance",
                "tolerance": 1e-9,
                "occurrences": []
            },
            "wa": {
                "description": "wa parameter",
                "canonical_source": {
                    "type": "json_file",
                    "path": "simulation_results.json",
                    "key_path": "observables.wa"
                },
                "canonical_value": -5.46990569108541e-06,
                "comparison": "relative_tolerance",
                "tolerance": 1e-9,
                "occurrences": []
            },
            "delta_chi2": {
                "description": "Delta chi2",
                "canonical_source": {
                    "type": "json_file",
                    "path": "simulation_results.json",
                    "key_path": "observables.delta_chi2_vs_lcdm"
                },
                "canonical_value": -0.0004271558779258555,
                "comparison": "relative_tolerance",
                "tolerance": 1e-9,
                "occurrences": []
            },
            "ln_bayes_factor": {
                "description": "ln Bayes factor",
                "canonical_source": {
                    "type": "json_file",
                    "path": "simulation_results.json",
                    "key_path": "observables.ln_bayes_factor"
                },
                "canonical_value": 0.00021357793896292776,
                "comparison": "relative_tolerance",
                "tolerance": 1e-9,
                "occurrences": []
            },
            "lean_theorem_count_tadpole": {
                "description": "Lean theorem count tadpole",
                "canonical_source": {
                    "type": "json_file",
                    "path": "simulation_results.json",
                    "key_path": "axe3_formal_proof.theorems_verified_count"
                },
                "canonical_value": 9,
                "comparison": "exact_integer",
                "occurrences": []
            },
            "lean_theorem_count_moonshine": {
                "description": "Lean theorem count moonshine",
                "canonical_source": {
                    "type": "json_file",
                    "path": "simulation_results.json",
                    "key_path": "axe4_moonshine_voa.theorems_verified_count"
                },
                "canonical_value": 16,
                "comparison": "exact_integer",
                "occurrences": []
            },
            "lean_theorem_count_kummer_tda": {
                "description": "Lean theorem count kummer_tda",
                "canonical_source": {
                    "type": "json_file",
                    "path": "simulation_results.json",
                    "key_path": "axe6_kummer_tda_langevin.theorems_verified_count"
                },
                "canonical_value": 11,
                "comparison": "exact_integer",
                "occurrences": []
            },
            "bubble_radius": {
                "description": "Bubble radius",
                "canonical_source": {
                    "type": "latex_constant",
                    "source": "papers/T-dulaity alone/T_duality_Alone.tex"
                },
                "canonical_value": 2.82,
                "comparison": "relative_tolerance",
                "tolerance": 1e-9,
                "occurrences": [
                    {
                        "file": "papers/T-dulaity alone/T_duality_Alone.tex",
                        "regex": r"R_c\s*≈\s*([0-9.]+)"
                    }
                ]
            },
            "bounce_action": {
                "description": "Bounce action",
                "canonical_source": {
                    "type": "latex_constant",
                    "source": "papers/T-dulaity alone/T_duality_Alone.tex"
                },
                "canonical_value": 137.15,
                "comparison": "relative_tolerance",
                "tolerance": 1e-9,
                "occurrences": [
                    {
                        "file": "papers/T-dulaity alone/T_duality_Alone.tex",
                        "regex": r"S_E\s*=\s*([0-9.]+)"
                    }
                ]
            }
        },
        "grep_commands_used": []
    }
    with open(f"{tmpdir}/audit/parameters.json", "w") as f:
        json.dump(params, f, indent=2)


def test_consistent_fixture():
    """Positive control: consistent fixture should exit 0."""
    tmpdir = tempfile.mkdtemp()
    try:
        create_consistent_fixture(tmpdir)

        # Copy the checker script to the fixture
        checker_src = Path(__file__).parent.parent / "scripts" / "cross_consistency_check.py"
        checker_dst = Path(tmpdir) / "scripts" / "cross_consistency_check.py"
        shutil.copy(checker_src, checker_dst)

        # Run the checker
        result = subprocess.run(
            [sys.executable, str(checker_dst)],
            cwd=tmpdir,
            capture_output=True,
            text=True
        )

        print("=== POSITIVE CONTROL OUTPUT ===")
        print(result.stdout)
        print(result.stderr)

        # Should exit 0 (consistent)
        assert result.returncode == 0, f"Expected exit code 0, got {result.returncode}"

    finally:
        shutil.rmtree(tmpdir)


def test_mutated_number():
    """Negative control: one mutated number should exit 1 and name the quantity."""
    tmpdir = tempfile.mkdtemp()
    try:
        create_consistent_fixture(tmpdir)

        # Mutate mapper_nodes in the JSON
        with open(f"{tmpdir}/simulation_results.json", "r") as f:
            sim_results = json.load(f)

        # Change mapper_nodes from 187 to 999
        sim_results["axe6_kummer_tda_langevin"]["mapper_nodes"] = 999

        with open(f"{tmpdir}/simulation_results.json", "w") as f:
            json.dump(sim_results, f, indent=2)

        # Copy the checker script
        checker_src = Path(__file__).parent.parent / "scripts" / "cross_consistency_check.py"
        checker_dst = Path(tmpdir) / "scripts" / "cross_consistency_check.py"
        shutil.copy(checker_src, checker_dst)

        # Run the checker
        result = subprocess.run(
            [sys.executable, str(checker_dst)],
            cwd=tmpdir,
            capture_output=True,
            text=True
        )

        print("=== NEGATIVE CONTROL OUTPUT ===")
        print(result.stdout)
        print(result.stderr)

        # Should exit 1 (inconsistent)
        assert result.returncode == 1, f"Expected exit code 1, got {result.returncode}"

        # Should mention mapper_nodes in the output
        output = result.stdout + result.stderr
        assert "mapper_nodes" in output, "Expected 'mapper_nodes' in output when it has a mismatch"

    finally:
        shutil.rmtree(tmpdir)


if __name__ == "__main__":
    print("Running positive control test...")
    test_consistent_fixture()
    print("PASS: Positive control\n")

    print("Running negative control test...")
    test_mutated_number()
    print("PASS: Negative control\n")

    print("All tests passed!")
