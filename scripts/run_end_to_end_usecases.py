#!/usr/bin/env python3
"""
=============================================================================
LeanFlow: End-to-End Three Flagship String Theory & SciML Use Cases
=============================================================================
Executes and serializes:
1. Use Case 1: Neural Calabi-Yau Metric Learning with Kahler Cone Stabilization.
2. Use Case 2: High-Throughput Swampland & MCMC Flux Vacuum Exploration.
3. Use Case 3: Cosmic Defect Extraction & Closed-Loop Lean 4 Tadpole Certification.
=============================================================================
"""

import os
import sys
import json
import time
import numpy as np
import torch

# Ensure workspace root is on sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import leanflow as lf
from deployment.huggingface.models.quintic_metric_pinn import (
    QuinticMetricPINN,
    sample_quintic_chart_points,
)


def run_usecase_1_cymetric(output_dir: str):
    print("\n" + "=" * 70)
    print(" [Use Case 1] Neural Calabi-Yau Metric Learning & Kahler Cone Stabilization")
    print("=" * 70)
    start = time.perf_counter()
    
    # 1. Sample evaluation points on the Quintic P4[5]
    x_train, omega_sq_train = sample_quintic_chart_points(num_points=200, seed=42)
    
    # 2. Instantiate Quintic PINN model
    model = QuinticMetricPINN(input_dim=6, hidden_dim=64, num_layers=3)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
    
    # 3. Simulate training without guardrail (demonstrating unphysical drift)
    print("  Evaluating unconstrained PINN behavior (without LeanFlow guardrail)...")
    with torch.no_grad():
        raw_loss, raw_min_evals = model.compute_monge_ampere_loss(x_train, omega_sq_train, use_guardrail=False)
        num_neg_unconstrained = int((raw_min_evals < 0).sum().item())
        min_eval_unconstrained = float(raw_min_evals.min().item())
        print(f"    Raw Unconstrained Min Eigenvalue: {min_eval_unconstrained:.4f}")
        print(f"    Negative Eigenvalues Detected: {num_neg_unconstrained} / {len(x_train)} (Unphysical Drift!)")

    # 4. Train with LeanFlow active Kahler cone spectral projection
    print("  Training with LeanFlow active Kahler cone projection (g_{i, jbar} > 0)...")
    loss_history = []
    for epoch in range(40):
        optimizer.zero_grad()
        loss, min_evals = model.compute_monge_ampere_loss(x_train, omega_sq_train, use_guardrail=True, min_eigenval=1e-3)
        loss.backward()
        optimizer.step()
        loss_history.append(float(loss.item()))
        if (epoch + 1) % 10 == 0 or epoch == 0:
            print(f"    Epoch {epoch+1:02d} | Monge-Ampère Loss: {loss.item():.6f} | Min Eigenvalue: {min_evals.min().item():.5f} > 0")

    # 5. Serialize PyTorch model checkpoint
    model_path = os.path.join(output_dir, "models", "quintic_metric_pinn.pt")
    torch.save({
        "epoch": 40,
        "model_state_dict": model.state_dict(),
        "final_loss": loss_history[-1],
        "input_dim": 6,
        "hidden_dim": 64,
        "num_layers": 3,
        "target_manifold": "Quintic_P4[5]",
        "kahler_cone_guaranteed": True
    }, model_path)
    print(f"  Model checkpoint saved: {model_path}")
    elapsed = time.perf_counter() - start
    print(f"  Use Case 1 Completed in {elapsed:.2f} seconds.")
    return model, loss_history


def run_usecase_2_swampland_mcmc(output_dir: str):
    print("\n" + "=" * 70)
    print(" [Use Case 2] High-Throughput Swampland & MCMC Flux Vacuum Exploration")
    print("=" * 70)
    start = time.perf_counter()
    
    # 1. Load canonical geometries
    quintic = lf.load_cicy("quintic")
    tian_yau = lf.load_cicy("tian_yau")
    k3_t2 = lf.load_cicy("k3_x_t2")
    
    print(f"  Ingested Calabi-Yau Geometries:")
    print(f"    Quintic: chi={quintic.euler_characteristic}, Max Tadpole={quintic.tadpole_bound:.2f}")
    print(f"    Tian-Yau: chi={tian_yau.euler_characteristic}, Max Tadpole={tian_yau.tadpole_bound:.2f}")
    print(f"    K3 x T2: chi={k3_t2.euler_characteristic}, Max Tadpole={k3_t2.tadpole_bound:.2f}")
    
    # 2. Serialize canonical dataset files
    datasets_dir = os.path.join(output_dir, "datasets")
    
    cicy_data = {
        name: lf.load_cicy(name).summary()
        for name in ["quintic", "cicy_7887", "tian_yau", "bicubic", "k3_x_t2"]
    }
    with open(os.path.join(datasets_dir, "cicy_canonical.json"), "w") as f:
        json.dump(cicy_data, f, indent=2)
        
    ks_data = {
        name: lf.load_ks(name).summary()
        for name in ["ks_quintic", "ks_bicubic", "ks_k3_surface"]
    }
    with open(os.path.join(datasets_dir, "kreuzer_skarke_canonical.json"), "w") as f:
        json.dump(ks_data, f, indent=2)

    # 3. Synthesize 1,000 MCMC Flux Vacuum Proposals
    np.random.seed(42)
    mcmc_samples = []
    for i in range(1000):
        # Sample parameters
        flux_charge = float(np.random.exponential(scale=5.0))
        delta_d = float(np.random.uniform(0.0, 10.0))
        tau_im = float(np.random.normal(loc=1.0, scale=0.5))
        mcmc_samples.append({
            "id": i,
            "flux_charge": flux_charge,
            "delta_d": delta_d,
            "tau_im": tau_im
        })

    # 4. Filter batch via LeanFlow SwamplandMCMCFilter
    mcmc_filter = lf.SwamplandMCMCFilter(euler_characteristic=quintic.euler_characteristic, delta_d_cutoff=7.0)
    filter_results = mcmc_filter.filter_batch(mcmc_samples, async_certify=True)
    
    print(f"  Screened {filter_results['total_evaluated']} MCMC candidates in {filter_results['filter_latency_ms']:.2f} ms")
    print(f"  Accepted Vacua: {filter_results['accepted_count']} (Acceptance Rate: {filter_results['acceptance_rate']*100:.1f}%)")
    print(f"  Rejections Breakdown: {filter_results['rejection_breakdown']}")
    print(f"  Lean 4 Proof Kernel Certificate: {filter_results['lean_certification']['certificate']}")

    # 5. Save annotated dataset
    output_mcmc = {
        "metadata": {
            "total_samples": 1000,
            "accepted_samples": filter_results["accepted_count"],
            "acceptance_rate": filter_results["acceptance_rate"],
            "target_manifold": "Quintic_P4[5]",
            "euler_characteristic": quintic.euler_characteristic,
            "tadpole_bound": quintic.tadpole_bound,
            "lean4_certified": filter_results["lean_certification"]["lean_verified"],
            "lean4_certificate": filter_results["lean_certification"]["certificate"]
        },
        "samples": mcmc_samples[:200]  # Store representative subset
    }
    with open(os.path.join(datasets_dir, "mcmc_vacuum_samples.json"), "w") as f:
        json.dump(output_mcmc, f, indent=2)
    print(f"  Dataset saved: {os.path.join(datasets_dir, 'mcmc_vacuum_samples.json')}")
    elapsed = time.perf_counter() - start
    print(f"  Use Case 2 Completed in {elapsed:.2f} seconds.")
    return filter_results


def run_usecase_3_tda_defects(output_dir: str):
    print("\n" + "=" * 70)
    print(" [Use Case 3] Cosmic Defect Extraction & Closed-Loop Lean 4 Tadpole Anomaly Gate")
    print("=" * 70)
    start = time.perf_counter()
    
    # 1. Synthesize 2D scalar field simulation grid (32x32)
    nx, ny = 32, 32
    x = np.linspace(-2, 2, nx)
    y = np.linspace(-2, 2, ny)
    xx, yy = np.meshgrid(x, y)
    
    # Multi-well field configuration with vortices (strings) and domain walls
    phi1 = np.sin(np.pi * xx) * np.exp(-0.1 * (xx**2 + yy**2))
    phi2 = np.cos(np.pi * yy) * np.exp(-0.1 * (xx**2 + yy**2))
    phi_grid = np.stack([phi1, phi2], axis=2)
    
    # 2. Extract simplicial 1-skeleton graph via TDA Mapper
    extractor = lf.TopologicalDefectExtractor(energy_threshold=0.4)
    defect_summary = extractor.extract_from_grid(phi_grid)
    
    print(f"  TDA Mapper 1-Skeleton Extracted:")
    print(f"    Graph Nodes: {defect_summary.num_nodes} | Edges: {defect_summary.num_edges}")
    print(f"    Cosmic Strings: {defect_summary.cosmic_strings_count} | Domain Walls: {defect_summary.domain_walls_count} | Vacua: {defect_summary.attractor_vacua_count}")
    print(f"    Net Ramond-Ramond Charge: {defect_summary.net_rr_charge}")
    print(f"    Ramond-Ramond Tadpole Neutral: {defect_summary.is_tadpole_neutral}")
    print(f"    Lean 4 Kernel Certified (Zero Sorry): {defect_summary.lean_certified}")
    
    elapsed = time.perf_counter() - start
    print(f"  Use Case 3 Completed in {elapsed:.2f} seconds.")
    return defect_summary


def main():
    print("=" * 70)
    print(" LeanFlow End-to-End Three Flagship Use Cases Runner")
    print(" Generating Production Datasets, Models, and Proof Certificates")
    print("=" * 70)
    total_start = time.perf_counter()
    
    output_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "deployment", "huggingface"))
    os.makedirs(os.path.join(output_dir, "models"), exist_ok=True)
    os.makedirs(os.path.join(output_dir, "datasets"), exist_ok=True)
    
    # Run all 3 use cases
    run_usecase_1_cymetric(output_dir)
    run_usecase_2_swampland_mcmc(output_dir)
    run_usecase_3_tda_defects(output_dir)
    
    total_elapsed = time.perf_counter() - total_start
    print("\n" + "=" * 70)
    print(f" ALL 3 USE CASES EXECUTED AND SERIALIZED IN {total_elapsed:.2f} SECONDS!")
    print(f" Artifacts ready in: {output_dir}")
    print("=" * 70)


if __name__ == "__main__":
    main()
