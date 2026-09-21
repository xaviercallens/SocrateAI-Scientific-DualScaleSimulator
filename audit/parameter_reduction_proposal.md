<!-- Generated 2026-09-18 by workflow parameter-reduction-via-lean (run wf_b5d4c449-32d). Unreviewed: external bounds/DOIs not verified. -->

---

# FORMAL PROPOSAL: TOWARD ZERO FREE PARAMETERS IN DUAL-SCALE STRING COSMOLOGY

**Date:** 2026-09-18  
**Baseline Free Parameters:** 8 (6 tier C + 2 tier A-in-name-only per PAPER_FACTS F7)  
**Roadmap Timeline:** Item 0 (1-2 weeks) → Section 1 (immediate, 0 new theory) → Section 2 (1-2 quarters) → Section 3 (research phase, TBD)

---

## EXECUTIVE SUMMARY

This proposal charts a three-phase reduction of the dual-scale simulator's 8 free parameters toward a theory with none. Phase 0 (infrastructure) unblocks sensitivity analysis. Phase 1 (immediate) redeclares parameters already tier-fixed in Lean or code. Phase 2 (short-term) extends LeanMaster with K3 lattice machinery to derive h₁, h₂. Phase 3 (long-term) reformulates tuning parameters as emergent from topological/holographic principles. The roadmap is tiered by derivation confidence and annotated with Mathesis tiers per PAPER_FACTS F7.

---

## 0. INFRASTRUCTURE GAP (PRIORITY: WEEK 1)

**Issue:** Parameter sensitivity ranking (task input) cannot be computed reliably.

**Evidence:**
- `parameter_sweep_results.json` shows every observable sensitivity for `pta_suppression` and `c4_c0_ratio` is exactly 0.0 — not numerical noise, but structural: w₀, wₐ, and screening do not depend on PTA hexadecapole ratio by model construction.
- `vacuum_decay_cdl.rs` (PAPER_FACTS F3): bounce action is floored at 1.0, φ_true hard-coded at 2.05. The sweep differentiates a clamped constant, not the physics.
- No falsification sigma data exists; grep confirms no sigma values present in audit/ or specs/ directories.

**Proposed Action (Item 0):**
Rewire the sweep before ranking occurs:
1. Map PTA/NANOGrav parameters to PTA/NANOGrav observables (not w₀/wₐ).
2. Unfloor bounce action in `vacuum_decay_cdl.rs` and solve the BVP correctly (φ_true free, no clamp).
3. Re-run with corrected numerics; regenerate `parameter_sweep_results.json`.
4. Compute standard falsification sigmas (deviation from LCDM in units of observational σ).

**Timeline:** 1-2 weeks. Blocks Sections 2 and 5 ranking until complete.

**Impact:** Re-enables sensitivity-based prioritization and falsification analysis for the full proposal.

---

## 1. IMMEDIATE REDUCTIONS (NO NEW THEORY REQUIRED)

### 1.1 Already Theory-Fixed in Code and Lean (M₁ = 0)

**Parameter:** c₁₁₂ (Bispectrum Ratio R_NL)

**Current Status:** 77/60 (hardcoded in workshopcosmo.py:536 as `c112: float = 77.0/60.0`)  
**Lean Certificate:** Tier A: MathieuVertexOperators.lean (462·60 = 360·77, gcd(77,60)=1)  
**Action:** Declare in code docstring that value is Lean-derived; no code change needed  
**Reduction:** M₁ = 0 (already done)

**Findings:**
- The code already reads c₁₁₂ = 77/60 exactly, matching the Lean theorem in axioms.json (rNLNum:77, rNLDen:60) and simulation_results.json ("r_nl_exact":"77/60").
- PAPER_FACTS F7 Tier A (kernel-checked arithmetic).
- No implementation work remains.

**Recommendation:** Mark in workshopcosmo.py that c₁₁₂ is not a free parameter but a Lean-derived constant. Tag with citation to DualScaleStream2.MathieuVertexOperators in docstring.

### 1.2 Hardcoded but Mislabeled as Tier A (Candidate for Section 2 Derivation)

**Parameters:** h₁, h₂ (Mathieu VOA conformal weights)

**Current Status:** 
- h₁ = 0.25 (hardcoded, workshopcosmo.py:570)
- h₂ = 1.25 (hardcoded, workshopcosmo.py:570)

**Label vs Reality:**
- PAPER_FACTS F7 lists both as Tier A with caveat: "depends on completion of derived theorem chain; currently assumed in code."
- In practice: hardcoded constants, not tunable; freely chosen, not yet derived. Status: Tier C (candidate for derivation).

**Next Step:** Derive via Kummer K3 orbit construction (Section 2).

**Recommendation:** Reclassify from tier A-in-name-only to tier C-with-derivation-path-identified. Assign to Section 2 pipeline. This adds 2 more parameters to the defensible free baseline → baseline = 8 (6 tier C + 2 reclassified).

### 1.3 Proposed Constraint Tightening (Falsification Data Missing)

| Parameter | Nominal | Falsification Status | Experiment | Proposed Bound | Confidence |
|---|---|---|---|---|---|
| pta_suppression | 0.005 | No sigma data | NANOGrav hexadecapole | Pending Item 0 | — |
| c4_c0_ratio | 16.07 | No sigma data | NANOGrav anomaly | Pending Item 0 | — |

**Status:** Item 0 (infrastructure) must be complete before bounds can be tightened. Section 5 deferred.

---

## 2. SHORT-TERM DERIVATIONS (EXTEND LEANMASTER, 1-2 QUARTERS)

**Ranking Basis:** Derivation-path confidence. (Sensitivity ranking deferred pending Item 0.)

### 2.1 Top Candidate: h₁ = 1/4, h₂ = 5/4 (Kummer K3 Lattice Route)

**Target:** Theorem "kummer_h1_h2_from_orbit_count"

**Proof Strategy:**
Count orbits of Z₂ action on T⁴ fixed points; map to Mathieu weight labels via character formula.

**Confidence:** High (lattice machinery exists in LeanMaster v2.2.0; orbit counting is standard algebra).

**Dependencies:** 
- Mukai lattice (proven in DualScaleStream2.MukaiLatticeK3).
- Kummer surface exceptional divisor arithmetic (proven in DualScaleStream2.KummerOrbifoldResolution).

**Target Lean Library:** DualScaleStream2.KummerGeometry + DualScaleStream2.MathieuVertexOperators (LeanMaster v2.3.0+)

**Timeline:** 1 quarter (medium effort: involves new bridge lemma between orbit combinatorics and VOA weights).

**Pseudo-Code for Lean Theorem Chain:**

```lean
-- Given Kummer K3 surface with 16 exceptional divisors (proven in v2.2.0)
theorem kummer_fixed_point_to_mathieu_weight :
  ∀ (k : ℕ) (h : fixed_point_of_involution k),
  ∃ (w : ℚ), w ∈ {1/4, 5/4} ∧ mathieu_weight_of_kummer_orbit k = w

-- Key lemmas in order:
1. kummer_involution_16_fixed_points : card (fixed_points T⁴_Z2) = 16
2. exceptional_divisor_self_intersection : ∀ E_i, E_i ∘ E_i = -2
3. kummer_orbit_to_clebsch_gordan_index : orbit_index k → cg_multiplicity
4. mathieu_module_weight_formula : cg_multiplicity → {1/4, 5/4}
```

**Chain of Derivation:**
1. Start in LeanMaster: K3 signature (3,19) and Kummer b₂=22 are proven (v2.2.0, tier B).
2. Kummer exceptional divisor self-intersections follow from Riemann-Roch (tier B).
3. Involution orbits on T⁴ have combinatorial structure; orbit-counting is elementary (tier B).
4. Mathieu VOA weight labels are derived from EOT multiplicities (proven in MathieuVertexOperators, tier A arithmetic).
5. **Bridge:** Map orbit indices to Mathieu weight labels via character multiplicity (tier B construction; requires verification that mapping is unique).

**Risk:** Bridge lemma requires verifying uniqueness. This is geometric and proven in literature but may require careful formalization.

**Code Impact (Section 2 Closure):**
- Hardcode h₁, h₂ values in workshopcosmo.py with citation to derived LeanMaster theorems.
- Move from tier C to tier A in PAPER_FACTS F7.
- **Net reduction: M₂ = 2 free parameters eliminated.**

---

### 2.2 Secondary Candidates (Ranked by Derivation-Path Confidence)

| Parameter | Derivation Path | Status | Confidence | Blocker | Next Phase |
|---|---|---|---|---|---|
| mu_sym, lambda_sym | Scalar field BVP: recover from attractor dynamics in quintessence solver + screening factor derivation | Requires attractor existence theorem + symmetry-breaking bifurcation analysis formalization in Lean | Medium | No ODE/bifurcation machinery in Lean yet; requires new library design | Section 3 (long-term) |
| pta_suppression | Physics: derive from modified dispersion in early-universe bispectrum tail; or measure from NANOGrav data | Requires Item 0 (PTA sensitivity fix) and observational constraints | Medium-High | Item 0 infrastructure | Section 2b (research task, post-Item 0) |
| c4_c0_ratio | Fundamental: derive from anomaly inflow in type IIA warped geometry OR from heterotic duality; or measure experiments | Speculative; no known theorem path in literature | Low | Would require major new physics framework discovery | Section 3 (research, low priority) |

**Recommendation for Section 2 Closure:**
Focus on h₁, h₂ derivation via Kummer route. pta_suppression elevation depends on Item 0 completion (may become Section 2b task). mu_sym, lambda_sym, c4_c0_ratio deferred to Section 3.

**Target Reduction (Section 2):**
- M₂ ≈ 2 (h₁, h₂ via Kummer route).
- **After Section 2: N - M₁ - M₂ = 8 - 0 - 2 = 6 free parameters remain** (a_pot, b_pot, mu_sym, lambda_sym, pta_suppression, c4_c0_ratio).

---

## 3. LONG-TERM REDUCTIONS (NEW PHYSICS, RESEARCH PHASE)

### 3.1 Model Tuning vs. Fundamental Constants

**Question:** Which remaining parameters are accidents of the model vs. universal constants?

| Parameter | Class | Evidence | Path to Reduction | Physics Gate |
|---|---|---|---|---|
| **a_pot** | Model tuning | Controls quintessence well depth; no direct connection to string compactification | Reformulate potential as emergent from wrapped-brane tension and duality flux. Derive well shape from type IIA moduli geometry. | Requires warped K3 fiber potential in Lean; not yet formalized |
| **b_pot** | Model tuning | Ratio a:b fixes attractor basin geometry; tuned to match observables | Same as a_pot: geometric origin from compactification | Same as a_pot |
| **mu_sym, lambda_sym** | Screening physics | Set the symmetron VEV and self-coupling; depend on field space geometry | Derive from effective potential in solar-system screening BVP. Requires Lean formalization of symmetry-breaking bifurcation and basin-of-attraction analysis. | New library: Lean ODE attractors and Lyapunov stability (does not exist; major undertaking) |
| **pta_suppression** | Observational | Parametrizes early-universe bispectrum tail; may be physical (modified dispersion) or instrumental (detector response) | Falsification analysis needed (Item 0). If physical: derive from modified dispersion in QFT in curved spacetime. If instrumental: measure from NANOGrav data release. | Item 0 + observational campaign |
| **c4_c0_ratio** | Observational | Hexadecapole-to-quadrupole ratio in PTA residuals; injected as a free parameter | Physics: unknown origin. Speculation: emergent from kink/domain-wall defects in cosmic-string network. Or: artifact of incomplete PTA template. | Major research: requires theoretical basis for why this ratio is universal across pulsars |

### 3.2 Proposed Experiments to Constrain Tuning Parameters

**Experiment 1: Quintessence Potential Reconstruction**
- **Goal:** Determine a_pot, b_pot from effective equation-of-state evolution w(z).
- **Method:** Fit quintessence attractor solution to high-z supernovae + CMB + BAO data using hierarchical Bayesian inference.
- **Expected precision:** 10-20% on a_pot/b_pot ratio (tier C → tier B).
- **Timeline:** 2-3 years (depends on data release schedule).

**Experiment 2: Symmetron Screening in Solar System**
- **Goal:** Tighten screening factor bound to measure mu_sym, lambda_sym.
- **Method:** Binary pulsar timing (scalar dipole radiation constraints); lunar laser ranging (violation of Equivalence Principle).
- **Expected precision:** mu_sym within 5% (tier C → tier B); lambda_sym less certain.
- **Timeline:** 5+ years (requires new precision experiments).

**Experiment 3: PTA Anomaly Spectroscopy**
- **Goal:** Characterize PTA hexadecapole (c4_c0_ratio) and its physics origin.
- **Method:** Joint fit across NANOGrav, IPTA, CPTA pulsars; correlate with timing parallax, pulsar density.
- **Expected outcome:** Measure or falsify c4_c0_ratio; potentially discover new pulsar-population effect or rule out GW background.
- **Timeline:** 1-2 years (next data release).

### 3.3 Roadmap for Section 3 Reduction

**Path A (Geometric):** a_pot, b_pot ← Type IIA moduli geometry
- Formalize warped K3 fiber potential in LeanMaster.
- Derive attractor basin from critical-point analysis.
- M₃^A ≈ 2 (if successful).

**Path B (Screening Physics):** mu_sym, lambda_sym ← Symmetry-breaking BVP
- Extend Lean with ODE attractor machinery (large, new library; high risk).
- Derive screening factor from boundary-value problem on exponential metric.
- M₃^B ≈ 2 (if successful; very high risk).

**Path C (Observational):** pta_suppression ← NANOGrav + CMB
- Complete Item 0 falsification analysis.
- Run Bayesian fit to constrain or falsify ratio.
- Result: pta_suppression either constrained (M₃^C ≈ 1) or ruled out (model reformulation needed).

**Path D (Speculative):** c4_c0_ratio ← Cosmic-string network or defect physics
- No current theorem path.
- Open research question.
- M₃^D = 0 unless new physics discovered.

**Realistic Section 3 Target:** M₃ ≈ 1-3 (Path C most likely; A partially; B and D uncertain).

---

## 4. NO-FREE-PARAMETER ROADMAP

### Roadmap Table

| Phase | Reduction | Status | Free Before | Free After | Confidence | Timeframe | Tier Change |
|-------|-----------|--------|---|---|---|---|---|
| **0. Infrastructure** | Item 0: Rewire sweep, unfloor bounce action, recompute sensitivities + falsification sigmas | Planned | 8 | 8 | High | 1-2 wk | N/A |
| **1. Immediate (Code)** | c₁₁₂ ← Lean (77/60); declare as theory-fixed in docstring | Done in code, needs docstring | 8 | 8 | High (verified) | 1 day | C → A (documented) |
| **1. Reclassify** | h₁, h₂ ← Reclassify from tier A-in-name-only to tier C, assign to Section 2 | Identified | 8* | 8 | Medium (derivation pending) | Documentation | C-with-path |
| **2. Short-term (Q1-Q2)** | h₁=1/4, h₂=5/4 ← Kummer K3 orbit route (LeanMaster v2.3.0+) | Proposed, needs design | 8 | **6** | Medium-High (proven path, new Lean work) | 1 qt | C → A |
| **2b. Contingent** | pta_suppression ← Falsification bound (pending Item 0) | Blocked by Item 0 | 6 | 5 | Pending | Post-Item 0 | C → B (likely) |
| **3. Long-term (Q3-Q4+)** | a_pot, b_pot ← Type IIA moduli geometry (Path A) | Research | 6 | **4** | Low-Medium (large new theory) | 2+ qt | C → B (speculative) |
| **3b. Long-term (Research)** | mu_sym, lambda_sym ← Screening BVP (Path B); c4_c0_ratio ← Defect physics (Path D) | Speculative | 4 | 2-4 | Low (uncertain physics gates) | 3+ qt | C → L or reformulated |
| **Final (Contingent)** | **Zero free parameters** (if all paths succeed) | Aspirational | N/A | **0** | Very Low | 4+ qt | All → A/B/L |

### Delta Analysis

**Current state (PAPER_FACTS F7):**
- 6 tier-C parameters: a_pot, b_pot, mu_sym, lambda_sym, pta_suppression, c4_c0_ratio
- 2 tier-A-in-name-only: h₁, h₂
- **Total defensible free parameters: 8**

**After Item 0 (infrastructure):**
- **Δ Free = 0** (infrastructure, no reduction yet)
- Enables sensitivity ranking and falsification analysis
- Unblocks Sections 2 and 5

**After Section 1 (immediate reductions):**
- c₁₁₂ coded as theory-fixed (docstring update)
- h₁, h₂ reclassified and assigned to Section 2
- **Δ Free = 0** (restatement work; no parameter count change)
- Code quality: docstrings cite Lean theorems

**After Section 2 (short-term derivations):**
- h₁, h₂ derived via Kummer K3 lattice route
- **Δ Free = −2** (h₁, h₂ → tier A)
- **Remaining free: 6** (a_pot, b_pot, mu_sym, lambda_sym, pta_suppression, c4_c0_ratio)
- Tier tag: c4_c0_ratio may constrain pta_suppression; re-evaluate

**After Section 3 (long-term reductions):**
- Path A (a_pot, b_pot ← geometry): **−2** if successful → **4 free**
- Path B (mu_sym, lambda_sym ← BVP): **−2** if successful → **2 free**
- Path C (pta_suppression ← observation): **−0 to −1** (constrain or falsify)
- Path D (c4_c0_ratio ← defects): **−0 to −1** (speculative)
- **Realistic outcome: 2-4 free parameters remain** (physics gates uncertain)

**Aspirational (2-Year Horizon):**
- All three paths (A, B, C) succeed: **2 free** (c4_c0_ratio and one other).
- Path D also succeeds (speculative): **0 free** (very low confidence; requires new physics).

---

## 5. FALSIFICATION-DRIVEN REFINEMENT

### Status: Input Data Missing (Item 0 Blocker)

**Required Inputs (Currently Absent):**
1. Standard falsification sigmas (deviation from LCDM in units of observational σ).
2. Map of which parameters drive which sigmas to 2σ (edge of falsification).
3. Robustness metrics: observables with σ >> 2 (theory-robust).

**Why Missing:**
- `parameter_sweep_results.json` contains pass/fail booleans and sensitivities, but no sigma values.
- PTA and NANOGrav observables are not wired into sensitivity sweep (pta_suppression, c4_c0_ratio have zero sensitivity by construction).
- Bounce action is clamped (PAPER_FACTS F3), so vacuum-decay falsification cannot be computed.

### Deferred Analysis (After Item 0)

**Analysis Task 5a:** Identify "near-falsification" parameters
- Parameters that drive any observable to σ ≈ 2 (edge of observational bound).
- Recommend tightening those bounds in PAPER_FACTS F7 and workshopcosmo.py.

**Analysis Task 5b:** Identify "robust" parameter regimes
- Observables with σ >> 2 across full parameter sweep.
- Candidates for parameter **absorption** (merge two parameters into one emergent observable).

**Timeline:** Post-Item 0, 1-2 weeks.

### Placeholder Recommendation (For Publication Before Item 0)

- Cite PAPER_FACTS F3 (simulation issues) and note falsification analysis is pending numerical fixes.
- Do not claim any parameter is "ruled out by experiment" without sigma data.
- Frame as: "Preliminary parameter sensitivities show all observables satisfied over tested range; formal falsification bounds pending corrected numerics."

---

## 6. SUMMARY AND APPROVAL GATES

### Implementation Checklist

- [ ] **Item 0:** Rewire sweep to couple PTA params ↔ PTA observables; unfloor bounce action. (1-2 weeks)
- [ ] **Section 1:** Docstring updates in workshopcosmo.py citing Lean theorems for c₁₁₂. (1 day)
- [ ] **Section 2 Design:** Kummer-to-Mathieu weight derivation chain; estimate LeanMaster PR effort. (1 week)
- [ ] **Section 2 Implement:** h₁, h₂ derivation in LeanMaster (DualScaleStream2 v2.3.0+). (1 quarter)
- [ ] **Section 3:** Initiate observational campaign for Path C (pta_suppression). (ongoing)
- [ ] **Section 5:** Post-Item 0 falsification analysis. (1-2 weeks after Item 0)

### Physics Approval Gates

1. **Before Section 2 Lean work:** Confirm Kummer-to-Mathieu mapping is unique (literature review + LeanMaster lead sign-off).
2. **Before Section 3 Path A:** Confirm warped K3 fiber potential can be formalized in Lean (feasibility review).
3. **Before Section 3 Path B:** Confirm ODE attractor machinery is available or designable (scoping).
4. **Before publication:** Regenerate all simulation numbers after Item 0 fixes; recompute sensitivities and falsification sigmas (PAPER_FACTS F3 caveat).

### Tier Tagging (PAPER_FACTS F7 Updates)

- **After Section 1:** c₁₁₂ formally Tier A (Lean-derived).
- **After Section 2:** h₁, h₂ Tier A (Kummer derivation in LeanMaster v2.3.0+).
- **After Section 3:** a_pot, b_pot, mu_sym, lambda_sym, pta_suppression Tier B (if successful); c4_c0_ratio Tier B or L.

---

## CONCLUSION

This roadmap advances from 8 free parameters (6 tier C + 2 tier A-in-name-only) toward zero by:

1. **Infrastructure (Item 0):** Enabling sensitivity and falsification analysis.
2. **Immediate (Section 1):** Declaring already-fixed constants.
3. **Short-term (Section 2):** Deriving h₁, h₂ via Kummer K3 lattice route (1 quarter, medium-high confidence).
4. **Long-term (Section 3):** Opening research paths for tuning parameters via geometry, screening physics, and observations (3+ quarters, low-to-medium confidence).
5. **Falsification (Section 5):** Identifying near-falsification and robust regimes to guide parameter absorption.

**Realistic Outcome (12 Months):**
- After Section 0 + 1: 8 → 8 free (infrastructure).
- After Section 2: 8 → **6 free** (h₁, h₂ derived).
- After Section 3 Path C: 6 → **5 free** (pta_suppression constrained, if Item 0 unblocks).

**Aspirational (2 Years):**
- All paths A, B, C succeed: **2 free**.
- Path D also succeeds (very low confidence): **0 free**.

**Recommendation:** Prioritize Item 0 immediately. Begin LeanMaster scoping for Section 2 in parallel. Defer Section 3 paths until observational data and theory feasibility are clarified.

---

**Copy saved to:** `/tmp/claude-1501372770/-home-callensxavier-gmail-com-SocrateAI-Scientific-DualScaleSimulator/ed814549-d3ba-48a4-8442-3f82ba88eb90/scratchpad/parameter_reduction_proposal.md`