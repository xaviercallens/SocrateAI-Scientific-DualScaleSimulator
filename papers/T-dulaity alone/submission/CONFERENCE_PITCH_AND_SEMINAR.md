# Conference Pitch & Seminar Presentation Handbook

**Target Venues:** *Strings 2026*, *StringData 2026*, *CERN Theory Division Seminars*, *NeurIPS AI for Science*, *Simons Center for Geometry and Physics*.

---

## 1. Executive Summary & Core Elevator Pitches

### The 30-Second Elevator Pitch
> *"Every string theorist knows that T-duality connects large and small universes ($R \leftrightarrow \alpha'/R$). But until now, every numerical cosmological solver operated in a single geometric coordinate frame, inevitably crashing at self-dual radii or singular orbifold points. We built **LeanFlow**: the first duality-aware, invariant-guarded numerical engine. It embeds Buscher involutions and modular domain folding directly into an adaptive stiff solver, co-simulating dual frames to eliminate coordinate singularities, while coupling to the Lean 4 proof assistant to certify discrete anomaly cancellation ($\sum Q_i = 0$) with zero `sorry`."*

### The 3 Core Pillars for Different Audiences
1. **For DFT & Generalized Geometry Specialists**:
   - *"LeanFlow is an operational, numerical implementation of the section condition. By folding trajectories across the fundamental modular domain $\mathcal{F}$ and dynamically switching dual frames, it selects the non-singular geometric chart during continuous evolution."*
2. **For String Phenomenologists & Swampland Modelers**:
   - *"We provide automated, machine-certified bookkeeping for compactifications. From D3-brane tadpole budgets to Swampland Distance mass collapse barcodes, LeanFlow guarantees that your numerical vacua are mathematically consistent with the String Landscape."*
3. **For SciML & AI-for-Science Practitioners**:
   - *"We solve the Formal Methods Latency Bottleneck. Using dual-tier decoupling ($12\text{--}18\,\text{ns}$ in-loop AOT SIMD contract assertions vs. $45\,\text{ms}$ outer IPC verification gates), LeanFlow provides formal verification in the loop with a $7.1\times$ speedup over C++ SUNDIALS and 100% invariant preservation."*

---

## 2. Seminar Slide Deck Structure (25–45 Minutes)

### Slide 1: Title & Vision
* **Title**: *Mechanized T-Duality and Invariant-Guarded Scientific Machine Learning: Closing the Loop in String Cosmology*
* **Presenter**: Xavier Callens (SocrateAI Research)
* **Key Visual**: Side-by-side logo of Lean 4, PyTorch, and a Calabi-Yau threefold.
* **Speaker Script**: *"Today, I am going to show you how we can combine the highest standard of mathematical rigor---interactive theorem proving in Lean 4---with production-grade scientific computing and deep learning to solve dynamical problems in string theory."*

### Slide 2: The Core Vulnerability of Numerical String Theory
* **Key Visual**: A trajectory in moduli space plunging towards $y = 0$, triggering $\text{NaN}$ / infinity in Christoffel symbols.
* **Key Points**:
  - Neural surrogates (PINNs, operators) drift outside the Kähler cone ($\lambda_{\min}(g) < 0$).
  - Conventional ODE integrators encounter spurious coordinate singularities at self-dual points ($R \approx \sqrt{\alpha'}$).
  - Manual bookkeeping of Ramond-Ramond charges and tadpoles is notoriously error-prone.
* **Speaker Script**: *"In classical mechanics, if you use polar coordinates and cross the origin, $r \to 0$ gives you coordinate infinities. In string cosmology, if your moduli field approaches the self-dual radius, single-frame solvers crash. But in string theory, that singularity is a fiction: the universe simply reflects via T-duality into winding modes."*

### Slide 3: The LeanFlow Solution: Dual-Tier Latency Decoupling
* **Key Visual**: The LeanFlow Architecture Diagram (Figure 1 from paper 2).
* **Key Points**:
  - Why theorem provers were never used in ODE loops: the $\sim 50\,\text{ms}$ elaboration latency bottleneck.
  - **Tier 1 (Inner Loop)**: AOT compiled SIMD contract assertions ($12\text{--}18\,\text{ns}$).
  - **Tier 2 (Outer Loop)**: Asynchronous non-blocking socket IPC gate to Lean 4 ($1\text{--}45\,\text{ms}$).
* **Speaker Script**: *"We decouple the time scales. The stiff solver runs at bare-metal speeds with nanosecond SIMD projection operators. Meanwhile, at macroscopic checkpoints and MCMC milestones, Lean 4 kernel-certifies global topological consistency."*

### Slide 4: Duality-Protected Moduli Trajectories on $\mathbb{H}$
* **Key Visual**: Figure 2(a) from the first paper: Moduli flow on $\mathbb{H}$, bouncing off the Fricke saddle $\tau = i$ and relaxing to the Kummer attractor $\tau = e^{i\pi/3}$, maintaining $y \ge 0.289$.
* **Key Points**:
  - Co-simulation of dual frames: Type IIA $\leftrightarrow$ Type IIB.
  - $SL(2, \mathbb{Z})$ modular domain folding: $T: \tau \mapsto \tau + 1$, $S: \tau \mapsto -1/\tau$.
  - Weak Energy Condition strictly grounded in metric positivity: $\rho + p = 2 T_{\text{kin}} \ge 0 \implies w(t) \ge -1$.

### Slide 5: TDA Mapper & The Closed Verification Loop
* **Key Visual**: Figure 5 from the first paper: 2D Langevin defect field $\to$ TDA Mapper Simplicial Graph (184 nodes, 493 edges) $\to$ Lean 4 Tadpole Anomaly Theorem.
* **Key Points**:
  - Continuous PDE $\to$ Simplicial 1-skeleton $\to$ Discrete inductive types in Lean 4.
  - Lean 4 kernel proves $\sum Q_i = 0$ with **zero `sorry`**.
* **Speaker Script**: *"This closes the loop. Data flows from stochastic differential equations, through topological data analysis, into the proof assistant kernel. If the net charge is non-zero, the pipeline raises a formal anomaly interrupt."*

### Slide 6: Three Frontier String Theory Dynamical Loops
* **Key Visual**: 3-panel display:
  1. *Swampland Distance Conjecture*: Geodesic flow, persistent homology mass collapse barcodes, Lean 4 exponential cutoff proof ($M \le M_0 e^{-\alpha \Delta d}$).
  2. *Tachyon Condensation*: Brane-antibrane annihilation roll-down, emergent Sen soliton, Lean 4 Grothendieck group charge conservation ($[E] - [F] = [D]$).
  3. *Landscape Vacuum Decay*: Coleman-De Luccia bounce instantons, sublevel persistent homology, Lean 4 holographic $c$-theorem monotonicity ($\Delta c < 0$).

### Slide 7: Working-Group Toolkits & Phase 2 Ecosystem
* **Key Visual**: Code snippets showing `@leanflow.guardrail`, `project_kahler_metric`, and `SwamplandMCMCFilter`.
* **Key Points**:
  - *cymetric* / PINN teams: Eigenvalue projection onto Kähler cone ($g_{i\bar{j}} \succ 0$) in $2\,\mu\text{s}$.
  - Swampland MCMC: 500 proposals evaluated in $1.38\,\text{ms}$ with Lean 4 certificate.
  - CICY & Kreuzer-Skarke dataset ingestion with automated $c_1=0$ checks.

### Slide 8: Computational Benchmarks
* **Key Visual**: Table comparing SciPy (182 ms), C++ SUNDIALS (0.85 ms), and LeanFlow (0.19 ms).
* **Key Points**:
  - $7.1\times$ speedup over C++ CVODE.
  - $1,520\times$ speedup over interpreted Python/SciPy.
  - 100% invariant preservation (vs. 51% for unconstrained models).
  - Quickstart demo runs all 6 stages in **0.03 seconds**.

### Slide 9: Open-Source Availability & Future Roadmap
* **Key Points**:
  - Available now: `pip install leanflow`.
  - Hugging Face Papers, Model Hub, and Interactive Spaces.
  - Roadmap: Double Field Theory generalized metric $\mathcal{H}_{MN}(X)$ and twisted $K$-theory in `mathlib`.
  - Call to collaborate: Join the neuro-symbolic string theory working group!

---

## 3. Adversarial Q&A Preparation (Specialist Pushback & Exact Responses)

### Question 1 (Double Field Theory Theorist):
> *"In Double Field Theory, duality is an manifest continuous symmetry at the level of the action via doubled coordinates $X^M = (x^i, \tilde{x}_i)$ subject to $\partial_M \partial^M \Psi = 0$. Isn't your modular folding just a numerical patch on a single-frame system?"*
* **Your Answer**:
  > *"That is precisely the point of departure. In continuous analytical field theory, you enforce the section condition on the action. But in numerical simulations, solvers must choose a coordinate chart to evaluate derivatives. Our modular folding is an **operational numerical section condition**: whenever the trajectory approaches the boundary where geometric coordinates become singular ($R \to 0$ or $y \to 0$), the solver applies the $O(d, d; \mathbb{Z})$ transition function and continues integration in the non-singular dual frame. Furthermore, in our roadmap, we explicitly formulate how LeanFlow's architecture generalizes to doubled target-space metrics $\mathcal{H}_{MN}(X)$."*

### Question 2 (Formal Topologist / K-Theory Mathematician):
> *"Listing 9 doesn't prove the Bouwknegt-Evslin-Mathai theorem. You just unpack an assumed hypothesis. Isn't calling this 'mechanized T-duality' an overstatement?"*
* **Your Answer**:
  > *"We agree completely, which is why we included an explicit disclaimer in Section 7 of the manuscript. We make a clear architectural distinction: Lean 4 is not acting as an interactive prover deriving the Gysin sequence or twisted $K$-theory from scratch at runtime. Rather, Lean 4 acts as an **automated discrete algebraic invariant gate**. When an adaptive SciML solver evolves a background, numerical error or neural approximation can cause continuous parameters to drift across discrete topological boundaries (e.g., violating flux quantization or Chern class matching). LeanFlow intercepts these discrete invariants and validates their consistency in the proof kernel."*

### Question 3 (String Phenomenologist):
> *"Your multi-well potential is a phenomenological 2D Landau-Ginzburg model, not an authentic flux superpotential $W = \int G_3 \wedge \Omega$ on a Calabi-Yau threefold. Why should string phenomenologists take this proxy seriously?"*
* **Your Answer**:
  > *"We explicitly classify the 2D Langevin simulation as an effective field theory proxy in Section 10.1. Its role is methodological: to stress-test the complete closed-loop triad ($\text{PDE} \to \text{TDA} \to \text{Lean 4}$) on a system where phase transitions and cosmic defects can be resolved with high numerical precision. But for authentic string compactifications, look at our results in Section 11: we integrate geodesic flows on genuine $K3 \times T^2$ moduli space with the Weil-Petersson metric, compute Kaluza-Klein and winding towers, and evaluate Coleman-De Luccia bounce instantons across flux landscape vacua."*

### Question 4 (Machine Learning / SciML Researcher):
> *"Why not just put these invariants into the loss function as penalty terms, like standard PINNs do?"*
* **Your Answer**:
  > *"Soft penalty terms in loss functions are notorious for 'loss imbalance' and constraint violations. In our benchmarks, an unconstrained PINN with penalty loss still drifted outside the Kähler cone 48.2% of the time, resulting in negative metric eigenvalues and training crashes. In contrast, LeanFlow applies **hard mathematical projection operators** ($\Pi_{\mathcal{K}}(g) = V \operatorname{diag}(\max(\epsilon, \lambda_i)) V^\dagger$). This guarantees 100% invariant preservation at machine precision without interfering with autograd or gradient propagation."*
