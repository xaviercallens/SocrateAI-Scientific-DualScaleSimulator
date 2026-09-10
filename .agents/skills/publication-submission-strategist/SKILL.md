---
name: publication-submission-strategist
description: >-
  Strategizes publication venue selection (JHEP, CPC, JCAP, SciPost, NeurIPS/ICLR),
  drafts tailored cover letters to journal editors, optimizes title/abstract framing
  for crossover theoretical physics and AI-for-Science audiences, and conducts
  pre-submission technical audits. Use this skill when preparing submission packages,
  deciding target journals, or writing cover letters.
---

# Publication Submission Strategist: String Theory & AI-for-Science

This skill guides the strategic preparation, venue selection, cover letter drafting, and submission packaging for the manuscript on LeanFlow, topological T-duality, and cosmological simulations.

---

## 1. Publication Venue Decision Matrix

Select the primary target venue based on the research objective and emphasis:

| Venue | Primary Focus | Best Framing Angle | Key Strengths to Emphasize | Potential Concerns to Preempt |
| :--- | :--- | :--- | :--- | :--- |
| **Computer Physics Communications (CPC)** | Computational physics software & algorithms | **Computational Tooling First**: A high-performance neuro-symbolic solver for stiff string cosmology. | $1,520\times$ speedup over SciPy BVP, $7.1\times$ over C++ SUNDIALS; dual-tier AOT/IPC architecture ($12\text{--}18\,\text{ns}$ vs $45\,\text{ms}$); reproducible open-source Rust/Lean 4 suite. | Preempt physics over-reach: clarify that Lean 4 verifies discrete algebraic invariants, not continuous analysis from scratch. |
| **Journal of High Energy Physics (JHEP)** | Theoretical high-energy physics, string theory | **Duality-Locked Moduli Dynamics**: Resolving non-linear string cosmology without coordinate singularities. | Continuous Buscher symmetry co-simulation ($R \leftrightarrow \alpha'/R$); Swampland Distance Conjecture KK mass collapse via TDA barcodes; $w(t) \ge -1$ guaranteed by $\tau_{\text{im}} > 0$. | Address DFT section condition; clarify that 2D Langevin defect potential is an EFT proxy for D-brane boundary states. |
| **Journal of Cosmology and Astroparticle Physics (JCAP)** | Theoretical & numerical cosmology, dark energy, modified gravity | **Certified Cosmological Stability**: Moduli relaxation, Symmetron screening, and energy conditions. | Non-linear Symmetron BVP yielding exact screening $\Delta R/R \approx 1.24 \times 10^{-4}$ matching Cassini bounds; absence of phantom divide crossings; Mukhanov-Sasaki perturbations. | Clarify the connection between 4D Symmetron screening and the underlying F-theory $K3 \times T^2$ compactification. |
| **SciPost Physics** | High-rigor open access physics | **Mechanized Theoretical Physics**: Closed-loop numerical discovery to formal algebraic verification. | Unbroken closed loop ($\text{Simulation} \to \text{TDA} \to \text{Lean 4}$); transparent, reproducible code with zero `sorry` proofs; open science and artifact availability. | Ensure all definitions are mathematically watertight, especially the Gysin sequence disclaimer in Listing 9. |
| **NeurIPS / ICLR (AI for Science)** | AI/ML for scientific discovery | **Neuro-Symbolic Scientific Discovery**: Bridging stiff PDE solvers and formal interactive theorem provers. | First end-to-end framework linking numerical PDEs, simplicial TDA clustering, and formal interactive theorem provers; dual-tier execution architecture solving latency. | Emphasize why formal verification (Lean 4) provides guarantees beyond heuristic neural surrogates. |

---

## 2. Pre-Submission Verification Checklist

Before generating any submission package or uploading to arXiv / journal portals, execute the following technical audit:

1. **LaTeX Document Compilation**:
   ```bash
   cd "papers/T-dulaity alone"
   pdflatex -interaction=nonstopmode T_duality_Alone.tex && pdflatex -interaction=nonstopmode T_duality_Alone.tex
   ```
   *Requirement*: Exit code 0, 0 undefined citations, 0 undefined references, 0 broken labels.

2. **Lean 4 Proof Kernel Compilation**:
   ```bash
   lake build
   ```
   *Requirement*: 32/32 jobs complete successfully with zero `sorry` and zero fatal errors.

3. **Full Simulation Regression Test**:
   ```bash
   pytest tests/test_workshopcosmo.py
   ```
   *Requirement*: 27/27 tests pass (100% pass rate).

4. **Figure Integrity Audit**:
   - `figure_cosmo_simulations.png`: Verify 4-panel $2\times 2$ layout without orphan panel (e).
   - `figure_three_frontier_loops.pdf`: Verify SDC, Tachyon roll-down, and CDL bounce curves render sharply.
   - `figure_tda_mapper.png`: Verify 1-skeleton graph with 184 nodes and 493 edges.
   - `figure_serverless_spot_scaling.pdf` & `figure_invariant_conservation_5min.pdf`: Verify clean axes without legacy WRF/CAMB labels.

---

## 3. Tailored Cover Letter Templates

### Template A: High-Energy Theory (JHEP / JCAP)
```latex
Dear Editors of the Journal of High Energy Physics,

We submit our manuscript entitled "Topological T-Duality and Emergent Symmetries on K3 x T2: A Certified Neuro-Symbolic Approach to String Phenomenology and Cosmological Dynamics" for consideration as a Regular Article.

A longstanding challenge in string phenomenology and early-universe cosmology is that non-linear moduli trajectories frequently encounter unphysical coordinate singularities (such as the self-dual radius R = sqrt(alpha') or modular boundaries) when integrated using conventional numerical differential equation solvers. Furthermore, verifying that dynamical defect networks satisfy global string-theoretic consistency criteria (such as Ramond-Ramond tadpole cancellation) has historically relied on manual post-hoc bookkeeping.

In this work, we present LeanFlow: a high-throughput neuro-symbolic framework coupling production-grade stiff numerical integration (rusty-SUNDIALS) with automated discrete algebraic invariant verification in the Lean 4 proof assistant. Our key physical and computational contributions include:
1. Native Duality Engine: Continuous Buscher symmetry locks (R <-> alpha'/R) and SL(2, Z) modular folding that enable dual-frame co-simulation, precluding coordinate singularities without metric breakdown (min(y) = 0.289).
2. The Closed Verification Loop: Coupling 2D stochastic Langevin symmetry breaking to Topological Data Analysis (TDA) Mapper to extract defect equivalence classes, with Lean 4 kernel-verifying that the network identically satisfies Ramond-Ramond tadpole cancellation (sum Q_i = 0) with zero sorry axioms.
3. Non-Perturbative Frontiers: Simulating the Swampland Distance Conjecture (exponential mass tower collapse tracked via persistent homology barcodes), Sen's tachyon condensation (soliton roll-down with discrete Grothendieck group charge conservation), and Coleman-De Luccia vacuum decay (verifying holographic c-theorem irreversibility).
4. Symmetron Screening & Energy Conditions: Dynamically resolving the non-linear Poisson BVP for stellar screening (Delta R/R ~ 1.24e-4, satisfying Cassini bounds) and proving that target-space metric positivity tau_im > 0 strictly enforces w(t) >= -1, precluding phantom energy.

We believe this paper will be of strong interest to theorists working on string cosmology, the Swampland program, and computational methods in high-energy physics. All simulation code, Lean 4 proof scripts, and TDA pipelines are open-source and fully reproducible.

Thank you for your consideration.

Sincerely,
The Authors
```

### Template B: Computational Physics (Computer Physics Communications)
```latex
Dear Editors of Computer Physics Communications,

We submit our manuscript entitled "LeanFlow: A Certified Neuro-Symbolic Solver for Stiff Moduli Dynamics, Topological T-Duality, and Non-Linear String Cosmology" for consideration.

Integrating stiff, highly non-linear differential equations in quantum gravity and string cosmology presents extreme stability challenges. Explicit integrators suffer from timestep collapse, while standard implicit solvers can drift into unphysical regions of moduli space. Moreover, interactive theorem provers (ITPs) like Lean 4 have historically been confined to pure mathematics due to the prohibitive interpretation latency of running formal proof kernels inside numerical inner loops.

LeanFlow resolves this bottleneck through a dual-tier execution architecture:
- Inner Loop: Lean-proven mathematical invariants (such as upper half-plane metric positivity and the Weak Energy Condition bound w >= -1) are compiled Ahead-Of-Time (AOT) into zero-overhead Rust SIMD contract assertions (12-18 ns latency), preserving bare-metal BDF stepping performance (0.12 ms step latency).
- Outer Loop: Full Lean 4 kernel verification is executed asynchronously via Unix domain socket IPC (45 ms latency) at macroscopic checkpoints and TDA defect extraction milestones.

On identical non-linear cosmological equations, LeanFlow achieves a 1,520x speedup over Python/SciPy BVP baselines and a 7.1x speedup over standard C++ SUNDIALS/CVODE, while cutting compute costs by 99.1% on serverless spot GPU/TPU infrastructure.

We look forward to your review.
```
