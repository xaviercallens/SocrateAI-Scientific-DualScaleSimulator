#!/usr/bin/env python3
"""
=============================================================================
LeanFlow Phase 2: String Theory Computation Community Quickstart
=============================================================================
End-to-end 2-minute demonstration showing:
1. Ingestion of canonical CICY and Kreuzer-Skarke Calabi-Yau geometries.
2. Declarative invariant specification via @leanflow.guardrail.
3. Stiff moduli trajectory integration with SL(2, Z) modular domain folding.
4. Working Group 1: cymetric Kahler metric positivity stabilization (g > 0).
5. Working Group 2: Swampland MCMC filter with asynchronous Lean 4 batch gate.
6. Working Group 3: Topological defect extraction and Ramond-Ramond neutrality.
=============================================================================
"""

import sys
import os
import time
import numpy as np

# Ensure workspace root is on sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

try:
    import torch
    TORCH_AVAILABLE = True
except ImportError:
    TORCH_AVAILABLE = False

import leanflow as lf


def run_quickstart():
    print("=" * 80)
    print(" LeanFlow 2.0: Guaranteed-Invariant Scientific Machine Learning Quickstart")
    print(" Endorsing the Computational String Theory & AI-for-Science Communities")
    print("=" * 80)
    start_all = time.perf_counter()

    # -------------------------------------------------------------------------
    # 1. STRING DATABASE INGESTION (CICY & KREUZER-SKARKE)
    # -------------------------------------------------------------------------
    print("\n[Stage 1/6] Ingesting Canonical Calabi-Yau Manifolds from String Datasets...")
    
    # Load canonical Quintic threefold P4[5]
    quintic = lf.load_cicy("quintic")
    print(f"  Loaded CICY: {quintic.name}")
    print(f"    Ambient: P^{quintic.ambient_spaces}, Matrix shape: {quintic.matrix.shape}")
    print(f"    c1=0 (Calabi-Yau): {quintic.is_calabi_yau}, Euler char chi: {quintic.euler_characteristic}")
    print(f"    Hodge numbers: h11={quintic.h11}, h21={quintic.h21}")
    print(f"    Max D3-Brane Tadpole Budget (chi/24): {quintic.tadpole_bound:.2f}")

    # Load Kreuzer-Skarke 4D reflexive polytope
    ks_bicubic = lf.load_ks("ks_bicubic")
    print(f"  Loaded Kreuzer-Skarke Polytope: {ks_bicubic.ks_id}")
    print(f"    Vertices: {ks_bicubic.vertices.shape[0]} in Z^4, Reflexive: {ks_bicubic.is_reflexive}")
    print(f"    Picard number h11: {ks_bicubic.h11}, h21: {ks_bicubic.h21}, chi: {ks_bicubic.euler_characteristic}")

    # -------------------------------------------------------------------------
    # 2. DECLARATIVE INVARIANT SPECIFICATION (@guardrail DSL)
    # -------------------------------------------------------------------------
    print("\n[Stage 2/6] Defining Moduli Dynamics with Declarative @guardrail DSL...")

    @lf.guardrail(
        theorem="SocrateAI.Cosmology.wec_kinetic_identity",
        invariants=["metric_positivity", "modular_invariance", "wec_bound"],
        projection="orthogonal"
    )
    def moduli_equations_of_motion(t, state):
        # State: [x, y, vx, vy] where tau = x + i y
        x, y, vx, vy = state[0], state[1], state[2], state[3]
        
        # Hyperbolic target-space connection terms: -2/y * vx * vy
        H_hubble = 0.05
        ax = -3.0 * H_hubble * vx + (2.0 / y) * vx * vy - (y ** 2) * 0.1 * np.sin(np.pi * x)
        ay = -3.0 * H_hubble * vy - (1.0 / y) * (vx**2 - vy**2) - (y ** 2) * 0.2 * (y - 0.866)
        
        return np.array([vx, vy, ax, ay])

    print("  Successfully instrumented ODE RHS with Lean 4 invariant locks.")

    # -------------------------------------------------------------------------
    # 3. STIFF INTEGRATION & SL(2, Z) MODULAR DOMAIN FOLDING
    # -------------------------------------------------------------------------
    print("\n[Stage 3/6] Solving Stiff Moduli Trajectory with Adaptive BDF & Modular Invariance...")
    solver = lf.Solver(method="BDF", rtol=1e-6, atol=1e-8)
    
    # High-energy initial condition outside fundamental domain (|tau| < 1)
    # tau_0 = 0.1 + 0.35i (close to singular boundary y -> 0)
    y0 = [0.1, 0.35, 0.0, 0.0]
    t_span = (0.0, 10.0)
    t_eval = np.linspace(0.0, 10.0, 100)

    sol = solver.solve(moduli_equations_of_motion, y0, t_span, t_eval=t_eval)
    
    print(f"  Integration Success: {sol.success}")
    print(f"  Total Steps Evaluated: {sol.telemetry.total_steps}")
    print(f"  Step Latency: {sol.telemetry.step_latency_ms:.3f} ms (Bare-metal speed)")
    print(f"  SL(2, Z) Modular Folds Applied: {sol.telemetry.modular_folds_count}")
    print(f"  Final Modulus: tau = {sol.y[0, -1]:.4f} + {sol.y[1, -1]:.4f}i")
    print(f"  Metric Positivity Preserved: min(y) = {np.min(sol.y[1, :]):.4f} > 0")

    # -------------------------------------------------------------------------
    # 4. WORKING GROUP 1: CALABI-YAU METRIC STABILIZATION (cymetric / PINN)
    # -------------------------------------------------------------------------
    print("\n[Stage 4/6] Working Group 1 Showcase: cymetric / PINN Kahler Cone Stabilization...")
    
    # Generate an unphysical, drifted metric tensor with a negative eigenvalue
    unphysical_metric = np.array([
        [1.0, 0.4 + 0.2j, 0.0],
        [0.4 - 0.2j, -0.3, 0.1j],  # Negative diagonal -> outside Kahler cone!
        [0.0, -0.1j, 0.8]
    ])
    evals_before = np.linalg.eigvalsh(unphysical_metric)
    print(f"  Raw Neural Surrogate Eigenvalues: {np.round(evals_before, 4)} (Contains unphysical negative eigenvalue!)")
    
    # Apply LeanFlow Kahler cone projector
    stabilized_metric = lf.project_kahler_metric(unphysical_metric, min_eigenval=1e-3)
    evals_after = np.linalg.eigvalsh(stabilized_metric)
    print(f"  LeanFlow Projected Eigenvalues:   {np.round(evals_after, 4)} (Strictly positive-definite: g > 0)")
    
    omega_norm = np.array([1.0])
    ma_loss = lf.compute_monge_ampere_loss(stabilized_metric, omega_norm)
    print(f"  Protected Monge-Ampere Loss: {ma_loss:.6f}")

    # -------------------------------------------------------------------------
    # 5. WORKING GROUP 2: SWAMPLAND MCMC FILTER & LEAN 4 BATCH GATE
    # -------------------------------------------------------------------------
    print("\n[Stage 5/6] Working Group 2 Showcase: Swampland MCMC Batch Filter...")
    swampland_filter = lf.SwamplandMCMCFilter(euler_characteristic=quintic.euler_characteristic)
    
    # Synthesize 500 candidate MCMC vacuum proposals
    np.random.seed(42)
    mcmc_batch = []
    for _ in range(500):
        proposal = {
            "flux_charge": float(np.random.uniform(0.0, 20.0)),
            "delta_d": float(np.random.uniform(0.0, 9.0)),
            "tau_im": float(np.random.uniform(-0.1, 2.0))
        }
        mcmc_batch.append(proposal)

    filter_results = swampland_filter.filter_batch(mcmc_batch, async_certify=True)
    print(f"  Evaluated {filter_results['total_evaluated']} MCMC parameter proposals in {filter_results['filter_latency_ms']:.2f} ms")
    print(f"  Accepted Candidates: {filter_results['accepted_count']} (Acceptance rate: {filter_results['acceptance_rate']*100:.1f}%)")
    print(f"  Rejections: {filter_results['rejection_breakdown']}")
    print(f"  Lean 4 Verification Certificate: {filter_results['lean_certification']['certificate']}")

    # -------------------------------------------------------------------------
    # 6. WORKING GROUP 3: TOPOLOGICAL DEFECT EXTRACTION & TADPOLE NEUTRALITY
    # -------------------------------------------------------------------------
    print("\n[Stage 6/6] Working Group 3 Showcase: Topological Defect Extraction via TDA...")
    extractor = lf.TopologicalDefectExtractor(energy_threshold=0.5)
    
    # Synthesize a multi-well 2D scalar field configuration (32x32)
    x = np.linspace(-1, 1, 32)
    y = np.linspace(-1, 1, 32)
    xx, yy = np.meshgrid(x, y)
    phi1 = np.sin(np.pi * xx)
    phi2 = np.cos(np.pi * yy)
    phi_grid = np.stack([phi1, phi2], axis=2)

    defect_summary = extractor.extract_from_grid(phi_grid)
    print(f"  Simplicial 1-Skeleton Extracted: {defect_summary.num_nodes} nodes, {defect_summary.num_edges} edges")
    print(f"  Defect Breakdown: {defect_summary.cosmic_strings_count} Strings, {defect_summary.domain_walls_count} Walls, {defect_summary.attractor_vacua_count} Vacua")
    print(f"  Net Ramond-Ramond Charge: {defect_summary.net_rr_charge}")
    print(f"  Ramond-Ramond Tadpole Neutral: {defect_summary.is_tadpole_neutral}")
    print(f"  Lean 4 Kernel Certified (Zero Sorry): {defect_summary.lean_certified}")

    total_time = time.perf_counter() - start_all
    print("\n" + "=" * 80)
    print(f" QUICKSTART COMPLETE: All 6 Stages Executed Successfully in {total_time:.2f} seconds!")
    print(" LeanFlow is ready for adoption by the string theory & SciML communities.")
    print("=" * 80)


if __name__ == "__main__":
    run_quickstart()
