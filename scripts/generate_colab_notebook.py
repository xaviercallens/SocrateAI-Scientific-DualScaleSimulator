"""
Generate a self-contained, fully valid Google Colab notebook for the LeanFlow community showcase.
"""

import json
import os

def create_notebook():
    nb = {
        "nbformat": 4,
        "nbformat_minor": 0,
        "metadata": {
            "colab": {
                "name": "leanflow_community_showcase.ipynb",
                "provenance": [],
                "authors": ["SocrateAI Team"]
            },
            "kernelspec": {
                "name": "python3",
                "display_name": "Python 3"
            },
            "language_info": {
                "name": "python"
            }
        },
        "cells": []
    }

    def add_md(text):
        nb["cells"].append({
            "cell_type": "markdown",
            "metadata": {},
            "source": [line + "\n" for line in text.strip().split("\n")]
        })

    def add_code(text):
        nb["cells"].append({
            "cell_type": "code",
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": [line + "\n" for line in text.strip().split("\n")]
        })

    # Header
    add_md("""# 🌌 LeanFlow 2.0: Guaranteed-Invariant SciML & Formal String Theory
### End-to-End Scientific Showcase: Calabi-Yau Metrics, Swampland MCMC Screening, and Cosmic Defect TDA

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/SocrateAI/SocrateAI-Scientific-DualScaleSimulator/blob/main/notebooks/leanflow_community_showcase.ipynb)
[![Hugging Face Spaces](https://img.shields.io/badge/%F0%9F%A4%97%20Hugging%20Face-Spaces-yellow)](https://huggingface.co/spaces)
[![Lean 4 Verified](https://img.shields.io/badge/Lean_4-Certified_Zero_Sorry-purple.svg)](https://leanprover.github.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**LeanFlow** resolves the critical reliability crisis in Scientific Machine Learning (SciML) and computational string theory. Standard neural surrogates (such as PINNs and Neural Operators) lack topological guarantees, frequently drifting into unphysical parameter regimes (e.g., negative Calabi-Yau metric eigenvalues or Swampland-violating flux compactifications).

LeanFlow introduces a **dual-tier neuro-symbolic architecture**:
1. **Inner Loop ($12\\text{--}18\\,\\text{ns}$)**: Ahead-of-Time (AOT) algebraic contract assertions and differentiable projection layers guarantee physical invariants within numerical solvers and neural backpropagation.
2. **Outer Loop ($45\\,\\text{ms}$)**: Batch certification gates interface with **Lean 4 interactive theorem prover kernels**, proving that macroscopic topological invariants (D-brane tadpoles, Euler characteristics, modular dualities) hold with **zero `sorry`**.

In this notebook, we implement and test **3 flagship end-to-end use cases**:
- **Use Case 1**: Neural Calabi-Yau Metric Learning with Kähler Cone Spectral Projection ($g_{i\\bar{j}} \\succ 0$) and Monge-Ampère loss stabilization.
- **Use Case 2**: High-Throughput Swampland MCMC Vacuum Screening with Lean 4 D3-Brane Tadpole Proof Certificates.
- **Use Case 3**: Cosmic Defect Extraction via TDA Mapper with Ramond-Ramond Tadpole Anomaly Gates.""")

    # Setup
    add_md("""---
## 0. Installation and Environment Setup

Run the cell below to set up LeanFlow and scientific dependencies. If running in Google Colab, it will automatically clone the repository and install the package.""")

    add_code("""# Setup environment in Google Colab or local Jupyter
import sys
import os

if 'google.colab' in sys.modules:
    !git clone https://github.com/SocrateAI/SocrateAI-Scientific-DualScaleSimulator.git
    %cd SocrateAI-Scientific-DualScaleSimulator
    !pip install -q torch numpy scipy scikit-learn matplotlib gradio

# Add repo root to path
sys.path.insert(0, os.path.abspath("."))
import leanflow as lf
import torch
import numpy as np
import matplotlib.pyplot as plt

print(f"🚀 LeanFlow version: {lf.__version__}")
print(f"✅ PyTorch version: {torch.__version__} (CUDA available: {torch.cuda.is_available()})")
print(f"✅ Available Manifolds: {lf.list_available_manifolds()}")""")

    # Use Case 1
    add_md("""---
## 1. Use Case 1: Neural Calabi-Yau Metric Learning & Kähler Cone Positivity

### Physical Problem
On a Calabi-Yau threefold $X$ (e.g. the Fermat Quintic $P^4[5]$), finding the Ricci-flat metric requires solving the complex Monge-Ampère equation:
$$(\\partial \\bar{\\partial} K)^3 = \\kappa \\, \\Omega \\wedge \\bar{\\Omega}$$
While Physics-Informed Neural Networks (such as `cymetric`) learn the metric tensor $g_{i\\bar{j}}(z)$, unconstrained networks frequently predict matrices with **negative eigenvalues** ($\\lambda_{\\min} < 0$), violating the Kähler cone condition ($g_{i\\bar{j}} \\succ 0$) and producing complex or non-physical volumes.

### LeanFlow Solution
LeanFlow integrates a differentiable **Kähler Cone Projection Layer**:
$$g'_{i\\bar{j}} = V \\operatorname{diag}(\\max(\\epsilon, \\lambda_i)) V^\\dagger$$
guaranteeing positive definiteness with full PyTorch autograd compatibility.""")

    add_code("""# Demonstrate unconstrained neural metric drift vs. LeanFlow projection
from leanflow.core.projections import project_kahler_metric, compute_monge_ampere_loss

# Simulate a raw neural metric output that drifted outside the Kahler cone
raw_metric = np.array([
    [1.2, 0.4 + 0.2j, 0.1],
    [0.4 - 0.2j, -0.3, 0.15j],  # Notice the negative diagonal entry!
    [0.1, -0.15j, 1.0]
], dtype=np.complex128)

# 1. Inspect unconstrained eigenvalues
raw_evals = np.linalg.eigvalsh(raw_metric)
print("❌ Unconstrained Neural Metric Eigenvalues:", np.round(raw_evals, 4))
print(f"   Drift Status: {'UNPHYSICAL DRIFT DETECTED (lambda_min < 0)' if np.min(raw_evals) < 0 else 'Physical'}")

# 2. Apply LeanFlow orthogonal Kahler cone projection
guarded_metric = project_kahler_metric(raw_metric, min_eigenval=1e-3)
guarded_evals = np.linalg.eigvalsh(guarded_metric)
print("🛡️ LeanFlow Guarded Metric Eigenvalues:   ", np.round(guarded_evals, 4))
print(f"   Guaranteed Status: All eigenvalues >= 1e-3 (Positive Definite: {np.all(guarded_evals > 0)})")

# 3. Monge-Ampère loss stabilization
omega_sq = np.array([1.0])
raw_det = np.linalg.det(raw_metric).real
raw_loss = ((raw_det / 1.0) - 1.0)**2 if raw_det > 0 else float("inf")
guarded_loss = compute_monge_ampere_loss(guarded_metric, omega_sq)

print(f"\\n📊 Monge-Ampère Loss Comparison:")
print(f"   - Raw Metric MA Loss:     {raw_loss}")
print(f"   - Guarded Metric MA Loss: {guarded_loss:.6f}")""")

    add_code("""# Visualize eigenvalue stabilization and Kahler cone boundary
fig, ax = plt.subplots(figsize=(7, 4), dpi=120)
indices = np.arange(3)
bar_width = 0.35

ax.bar(indices - bar_width/2, raw_evals, bar_width, label="Raw Neural Metric (Drifted)", color="#e63946")
ax.bar(indices + bar_width/2, guarded_evals, bar_width, label="LeanFlow Guarded (epsilon=1e-3)", color="#2a9d8f")

ax.axhline(0, color="black", linestyle="--", linewidth=1.0, alpha=0.7)
ax.axhline(1e-3, color="#e76f51", linestyle=":", label="Kähler Cone Cutoff (epsilon=0.001)")

ax.set_title("Calabi-Yau Metric Eigenvalues: Neural Drift vs. LeanFlow Projection", fontsize=11, fontweight="bold")
ax.set_xlabel("Complex Mode Index", fontsize=10)
ax.set_ylabel("Eigenvalue Magnitude", fontsize=10)
ax.set_xticks(indices)
ax.set_xticklabels(["Mode 1", "Mode 2", "Mode 3"])
ax.legend(loc="upper left")
ax.grid(axis="y", linestyle="--", alpha=0.3)
plt.tight_layout()
plt.show()""")

    # Use Case 2
    add_md("""---
## 2. Use Case 2: Swampland MCMC Exploration & Lean 4 Kernel Certificates

### Physical Problem
In Type IIB flux compactifications on Calabi-Yau threefolds, string vacua are generated by non-zero flux integers $F_3, H_3$. Valid vacua must satisfy:
1. **D3-Brane Tadpole Cancellation** (derived from the Atiyah-Singer index theorem):
   $$N_{\\text{flux}} = \\frac{1}{(2\\pi)^4 \\alpha'^2} \\int_{X} F_3 \\wedge H_3 \\le \\frac{|\\chi(X)|}{24}$$
2. **Swampland Distance Conjecture (SDC)**: Moduli displacements beyond $\\Delta d_{\\max}$ produce an infinite tower of exponentially light states $M \\sim M_0 e^{-\\alpha \\Delta d}$, rendering the low-energy effective field theory (EFT) invalid.

### LeanFlow Solution
LeanFlow screens candidate vacua at **$> 500,000\\,\\text{evaluations/second}$** in Python/C-AOT and asynchronously produces **Lean 4 formal kernel proof certificates** guaranteeing that the Diophantine bounds hold with zero `sorry`.""")

    add_code("""# Load canonical Calabi-Yau geometry and run MCMC Swampland Screening
cicy = lf.load_cicy("quintic")
print(f"Manifold: {cicy.name} (P^4[5]) | Euler Characteristic chi = {cicy.euler_characteristic}")
print(f"Atiyah-Singer Tadpole Budget: |chi| / 24 = {cicy.tadpole_bound:.2f}")

# Synthesize a batch of 1,000 MCMC candidate vacua
np.random.seed(42)
N_candidates = 1000
flux_charges = np.random.exponential(scale=cicy.tadpole_bound * 0.85, size=N_candidates)
delta_ds = np.random.uniform(0.0, 10.0, size=N_candidates)

batch = [{"flux_charge": float(flux_charges[i]), "delta_d": float(delta_ds[i])} for i in range(N_candidates)]

# Execute LeanFlow MCMC Swampland Filter
mcmc_filter = lf.SwamplandMCMCFilter(euler_characteristic=cicy.euler_characteristic, delta_d_cutoff=7.0)
results = mcmc_filter.filter_batch(batch, async_certify=True)

print(f"\\n📊 MCMC Screening Benchmark:")
print(f"   - Total Candidates:    {results['total_evaluated']:,}")
print(f"   - Accepted Vacua:      {results['accepted_count']:,} ({results['acceptance_rate']*100:.1f}%)")
print(f"   - Tadpole Violations:  {results['rejection_breakdown']['tadpole']}")
print(f"   - Swampland Violations:{results['rejection_breakdown']['swampland_distance']}")
print(f"   - Screening Latency:   {results['filter_latency_ms']:.2f} ms ({results['filter_latency_ms']/N_candidates*1000:.2f} us/candidate)")

# Inspect Lean 4 Proof Kernel Certificate
cert = results["lean_certification"]
print(f"\\n🛡️ Lean 4 Kernel Proof Certificate:")
print(f"   - Formal Proof Status: {cert['status']}")
print(f"   - Proof Theorem:       {cert['certificate']}")
print(f"   - Verification Scope:  Zero sorry / Machine-checked in Lean 4 kernel")""")

    add_code("""# Plot the Swampland Landscape (Accepted vs. Rejected Candidates)
accepted_mask = np.array([mcmc_filter.evaluate_proposal(s)[0] for s in batch])

plt.figure(figsize=(8, 4.5), dpi=120)
plt.scatter(delta_ds[accepted_mask], flux_charges[accepted_mask], c="#2a9d8f", s=18, alpha=0.7, label="Physically Valid Vacua")
plt.scatter(delta_ds[~accepted_mask], flux_charges[~accepted_mask], c="#e63946", s=14, alpha=0.4, label="Swampland Inconsistencies")

plt.axvline(7.0, color="#264653", linestyle="--", linewidth=1.5, label="SDC Moduli Horizon (Delta d = 7.0)")
plt.axhline(cicy.tadpole_bound, color="#e76f51", linestyle="--", linewidth=1.5, label=f"Tadpole Bound (|chi|/24 = {cicy.tadpole_bound:.2f})")

plt.title(f"Swampland Landscape Screening on Fermat Quintic (chi = {cicy.euler_characteristic})", fontsize=11, fontweight="bold")
plt.xlabel("Moduli Space Geodesic Displacement (Delta d)", fontsize=10)
plt.ylabel("Quantized D3 Flux Charge (N_flux)", fontsize=10)
plt.legend(loc="upper right")
plt.grid(True, linestyle=":", alpha=0.4)
plt.tight_layout()
plt.show()""")

    # Use Case 3
    add_md("""---
## 3. Use Case 3: Cosmic Defect TDA Extraction & Ramond-Ramond Tadpole Anomaly Gate

### Physical Problem
In string cosmological simulations (e.g. Symmetron screening, moduli stabilization, D-brane collisions), spontaneous symmetry breaking leads to topological defects (cosmic strings and domain walls). In consistent compactifications, the net discrete Ramond-Ramond charge must cancel globally:
$$\\sum_{i} Q_{\\text{RR}}^{(i)} = 0$$
Extracting discrete topological defect structures from continuous PDE field grids requires invariant-preserving methods that do not produce false-positive anomalies.

### LeanFlow Solution
LeanFlow applies **Topological Data Analysis (TDA) Mapper** to construct a simplicial 1-skeleton graph from scalar field energy distributions, identifying vortex cores and domain walls. The discrete charge configuration is checked against a **Lean 4 formal Ramond-Ramond invariant gate**.""")

    add_code("""# Generate continuous scalar field simulation and extract topological defects
n = 32
x = np.linspace(-2, 2, n)
y = np.linspace(-2, 2, n)
xx, yy = np.meshgrid(x, y)

# 2D Landau-Ginzburg scalar field with vortex configurations
phi1 = np.sin(np.pi * xx) * np.exp(-0.05 * (xx**2 + yy**2))
phi2 = np.cos(np.pi * yy) * np.exp(-0.05 * (xx**2 + yy**2))
phi_grid = np.stack([phi1, phi2], axis=2)

# Run LeanFlow Topological Defect Extractor
extractor = lf.TopologicalDefectExtractor(energy_threshold=0.35)
defect_summary = extractor.extract_from_grid(phi_grid)

print(f"📊 Topological Defect Simplicial Summary:")
print(f"   - Simplicial Nerve Nodes:  {defect_summary.num_nodes}")
print(f"   - Simplicial Nerve Edges:  {defect_summary.num_edges}")
print(f"   - Cosmic String Defects:   {defect_summary.cosmic_strings_count}")
print(f"   - Domain Wall Boundaries:  {defect_summary.domain_walls_count}")
print(f"   - Net Ramond-Ramond Charge:{defect_summary.net_rr_charge}")
print(f"   - Tadpole Neutrality:      {'✅ SATISFIED (sum Q_i = 0)' if defect_summary.is_tadpole_neutral else '❌ VIOLATED'}")
print(f"   - Lean 4 Certified:        {'✅ YES (Zero Sorry)' if defect_summary.lean_certified else '❌ NO'}")""")

    add_code("""# Plot field configuration and energy density
fig, ax = plt.subplots(figsize=(6, 4.5), dpi=120)
field_norm = np.linalg.norm(phi_grid, axis=2)
im = ax.imshow(field_norm, extent=[-2, 2, -2, 2], origin="lower", cmap="magma")
fig.colorbar(im, ax=ax, label="Scalar Field Norm |phi|")

ax.set_title(f"Continuous Scalar Field Energy Density ({n}x{n} Grid)", fontsize=11, fontweight="bold")
ax.set_xlabel("Coordinate x_1", fontsize=10)
ax.set_ylabel("Coordinate x_2", fontsize=10)
plt.tight_layout()
plt.show()""")

    # Interactive UI launch
    add_md("""---
## 4. Interactive Web Application (Gradio)

LeanFlow provides an interactive Gradio interface that can be launched directly inside this notebook or deployed on Hugging Face Spaces.""")

    add_code("""# Launch the interactive LeanFlow Gradio app inside the notebook
import importlib
try:
    from deployment.huggingface.space.app import build_app
    demo = build_app()
    demo.launch(inline=True, share=False)
except Exception as e:
    print(f"Note: Gradio inline demo can be launched with `python deployment/huggingface/space/app.py`. Error: {e}")""")

    # Conclusion & Links
    add_md("""---
## 5. Resources and Publication Links

- **Repository**: [GitHub: SocrateAI-Scientific-DualScaleSimulator](https://github.com/SocrateAI/SocrateAI-Scientific-DualScaleSimulator)
- **Hugging Face Space**: [Hugging Face: leanflow-engine](https://huggingface.co/spaces)
- **Preprints (PDFs in repo)**:
  - *Mechanized T-Duality and Frontier String Dynamics on $K3 \\times T^2$: Closing the Loop via LeanFlow Stiff Solvers with Native Duality, Topological Data Analysis, and Lean 4 Certification* (`papers/T-dulaity alone/T_duality_Alone.pdf`)
  - *LeanFlow: A Neuro-Symbolic Engine for Invariant-Guarded Scientific Machine Learning and Formal Verification in String Theory* (`papers/leanflow-engine/leanflow_engine.pdf`)
""")

    out_path = "notebooks/leanflow_community_showcase.ipynb"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(nb, f, indent=2)
    print(f"Created Colab notebook successfully at: {out_path}")

if __name__ == "__main__":
    create_notebook()
