# Publication & Scientific Framing Rules for String Cosmology & LeanFlow

These rules govern all writing, paper revisions, cover letters, abstract polishing, and referee rebuttals for the manuscript and related publications.

---

## 1. Primary Positioning Directive: Tooling & Verification First

* **Always position LeanFlow primarily as a computational and formal verification paradigm** for string theory, quantum gravity, and non-linear cosmology.
* **Frame physical models as rigorous benchmark realizations**: Present the $K3 \times T^2$ moduli flow, Symmetron screening, Swampland Distance Conjecture, and Kibble-Zurek defect simulations as explicit dynamical demonstrations showcasing how LeanFlow's neuro-symbolic triad ($\text{Simulation} \to \text{TDA} \to \text{Lean 4}$) prevents numerical drift, metric breakdown, and unphysical pathologies.
* **Never claim to replace analytical string theory**: Emphasize that LeanFlow mechanizes the bookkeeping of discrete topological invariants and enforces geometric guardrails, while the underlying physical field equations are integrated at bare-metal speeds.

---

## 2. Double Field Theory (DFT) & Generalized Geometry (Hull, Zwiebach, Hohm, Berman)

* **The "Algorithmic Section Condition" Pitch**: Explicitly frame the native duality engine and continuous Buscher symmetry co-simulation ($R \leftrightarrow \alpha'/R$, $\tau \leftrightarrow -1/\tau$) as an **operational numerical section condition**. 
* **Address Doubled Action Invariance**: 
  - Acknowledge that while DFT enforces $O(D, D; \mathbb{Z})$ invariance at the level of the continuous doubled action via coordinates $X^M = (x^i, \tilde{x}_i)$ subject to $\partial_M \partial^M \Psi = 0$, LeanFlow implements this symmetry dynamically as an invariant manifold projector on the effective moduli trajectory.
  - Highlight that this prevents coordinate singularities near self-dual string scales ($R \approx \sqrt{\alpha'}$) and modular boundaries without metric breakdown ($\min(y) = 0.289$).
  - Present the future extension to genuine doubled target-space metrics $\mathcal{H}_{MN}(X)$ as a natural next step for field-theoretic PDEs.

---

## 3. Topological T-Duality & Twisted $K$-Theory (Bouwknegt, Evslin, Mathai, Rosenberg)

* **Formalization Depth Guardrail**: Never claim to have formalized full differential geometry, Courant algebroids, or infinite-dimensional twisted $K$-theory from scratch in Lean 4.
* **Acknowledge the Gysin & Mathlib Roadmap**: Always acknowledge parenthetically that authentic topological T-duality from first principles requires mechanizing the Gysin exact sequence in cohomology, Courant algebroid automorphisms, and the degree-shifting isomorphism of twisted $K$-theory $K^{*+1}(E, H) \cong K^*(E^{\vee}, H^{\vee})$ in `mathlib`.
* **Discrete Invariant Gate Role**: Frame the current Lean 4 formalization as an automated **discrete algebraic invariant gate** that verifies the algebraic closure of the flux-topology exchange ($c_1(\widehat{E}) = \pi_* H$) and discrete Grothendieck group classes over modeled bundle structures during numerical evolution.

---

## 4. String Phenomenology & Swampland Program (Vafa, Valenzuela, Quevedo, McAllister)

* **EFT Proxy Distinction for Langevin Dynamics**: In the Kibble-Zurek multi-well simulation, explicitly acknowledge that the 2D Cartesian Landau-Ginzburg potential is an **effective field theory proxy** designed to test the numerical-to-symbolic closed loop.
* **Articulate the Genuine String Embedding**: Detail the theoretical requirements for embedding these proxy defects into authentic compactifications:
  - Domain walls $\to$ D-branes wrapping compact cycles with Callan-Harvey anomaly inflow ($\Delta F = H \wedge C$).
  - Cosmic strings $\to$ D1-strings coupled to $C_2$/$B_2$ with Green-Schwarz anomaly factorization ($I_8 = \frac{1}{2} X_4^2$).
  - Attractor vacua $\to$ moduli fixed points with orientifold tadpole cancellation ($\sum Q_{D7} + Q_{O7} = 0$).
* **Highlight Swampland Distance & Metric Positivity**: Strongly emphasize our genuine string results: TDA persistent homology tracking Kaluza-Klein mass tower collapse ($M \le M_0 e^{-\alpha \Delta d}$), and target-space metric positivity $\tau_{\text{im}} > 0$ strictly enforcing $w(t) \ge -1$ (ruling out phantom energy without unphysical pressure assumptions).

---

## 5. Computational String Theory & AI-for-Science (Primary Venue Audience)

* **Lead with the Closed Loop**: Champion the unbroken neuro-symbolic circuit: $\text{Numerical Integration (rusty-SUNDIALS)} \to \text{Topological Extraction (TDA Mapper)} \to \text{Formal Kernel Certification (Lean 4)}$.
* **Highlight the Dual-Tier Architecture**:
  - **Inner Loop**: AOT-compiled Rust SIMD contract assertions running in $\approx 12\text{--}18\,\text{ns}$, preserving the $0.12\,\text{ms}$ microsecond BDF step latency.
  - **Outer Loop**: Dynamic Lean 4 proof verification running asynchronously via non-blocking Unix socket IPC / C-FFI ($\approx 45\,\text{ms}$ batch latency) at macroscopic checkpoints.
* **Quantify Performance & Cost**: Emphasize the $1,520\times$ speedup over Python/SciPy BVP, $7.1\times$ speedup over C++ SUNDIALS via modular folding, and 99.1% cost reduction on serverless spot GPU/TPU infrastructure ($\text{min}=0$).
