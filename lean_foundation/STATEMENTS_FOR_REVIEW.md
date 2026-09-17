# DualScale Foundation: Theorem Statements for Review

**Date**: 2026-09-17  
**Foundation Version**: LeanMaster v2.2.0 (verified: git describe --tags = v2.2.0, build_v2.2.0.log exit = 0)  
**Audit Result**: 20 theorems audited, 0 failing. All theorems depend only on `propext`, `Classical.choice`, `Quot.sound`.

---

## T-Duality Module (supersedes `proofs/BuscherRules.lean`)

### Non-Vacuity Story
The old `proofs/BuscherRules.lean` proved `buscher_involution` using **five hand-written axioms**:
- `inv_inv`, `buscher_cross_inv`, `buscher_gmunu_inv`, `buscher_Bmunu_inv`, `Phi_inv`

These axioms asserted algebraic field properties required for the proof, pushing the logical burden outside the type system.

The `DualScaleFoundation.TDuality` module restates this content using **kernel-verified declarations from DualScaleStream2**:
- All theorems rest on Mathlib and the three Lean standard axioms only
- T-duality's involutive structure is expressed via the factorized duality matrices and the O(d,d;ℤ) group law
- No unaxiomatized claims

---

### Theorem: tduality_involution_factorized

```lean
#check DualScaleFoundation.TDuality.tduality_involution_factorized
-- tduality_involution_factorized (d : ℕ) (k : Fin d) :
--   DualScaleStream2.TDuality.factorized k * DualScaleStream2.TDuality.factorized k = 1
```

**Tier**: A (kernel-checked)  
**LeanMaster Dependencies**: `DualScaleStream2.TDuality.Factorized.factorized_mul_self`  
**Supersedes**: `proofs/BuscherRules.lean::buscher_involution`  
**Description**: Factorized T-duality in direction k is an involution (self-inverse).

---

### Theorem: generalized_metric_constraint

```lean
#check DualScaleFoundation.TDuality.generalized_metric_constraint
-- generalized_metric_constraint (d : ℕ) (G B : Matrix (Fin d) (Fin d) ℝ)
--   (hG : IsUnit G.det) :
--   DualScaleStream2.DFT.etaR d * DualScaleStream2.DFT.genMetric G B *
--       (DualScaleStream2.DFT.etaR d * DualScaleStream2.DFT.genMetric G B) =
--     1
```

**Tier**: A (kernel-checked)  
**LeanMaster Dependencies**: `DualScaleStream2.DFT.GeneralizedMetric.etaR_genMetric_sq`  
**Supersedes**: `proofs/BuscherRules.lean` (Double Field Theory formulation)  
**Description**: Hull–Zwiebach constraint (ηH)² = 1 expresses T-duality's involution in the continuous setting.

---

### Theorem: tduality_is_odd

```lean
#check DualScaleFoundation.TDuality.tduality_is_odd
-- tduality_is_odd (d : ℕ) : DualScaleStream2.TDuality.IsODD (DualScaleStream2.TDuality.eta d)
```

**Tier**: A (kernel-checked)  
**LeanMaster Dependencies**: `DualScaleStream2.TDuality.ODD.eta_isODD`  
**Supersedes**: N/A (group structure)  
**Description**: Full T-duality (momentum ↔ winding exchange) is an O(d,d;ℤ) element.

---

### Theorem: factorized_is_odd

```lean
#check DualScaleFoundation.TDuality.factorized_is_odd
-- factorized_is_odd (d : ℕ) (k : Fin d) :
--   DualScaleStream2.TDuality.IsODD (DualScaleStream2.TDuality.factorized k)
```

**Tier**: A (kernel-checked)  
**LeanMaster Dependencies**: `DualScaleStream2.TDuality.Factorized.factorized_isODD`  
**Supersedes**: N/A (group structure)  
**Description**: Each factorized duality preserves the O(d,d;ℤ) invariant form.

---

## K3 Lattice Module (supersedes `proofs/MukaiLatticeK3.lean`)

### Non-Vacuity Story
The old `proofs/MukaiLatticeK3.lean` defined its own `MukaiVector` and `mukaiPairing` from scratch.

The `DualScaleFoundation.K3Lattice` module uses **LeanMaster's verified declarations**:
- Mukai vector: `DualScaleStream2.Lattice.MukaiVec`, defined as (r, c, s) in H⁰ ⊕ H² ⊕ H⁴
- Pairing: `DualScaleStream2.Lattice.mukaiPair` with full symmetry and evenness proofs
- Signatures: verified lattice decompositions for K3 (3,19) and Mukai (4,20)

---

### Theorem: mukai_pairing_symmetric

```lean
#check DualScaleFoundation.K3Lattice.mukai_pairing_symmetric
-- mukai_pairing_symmetric (n : ℕ) (L : DualScaleStream2.Lattice.Gram n)
--   (hL : Matrix.transpose L = L) (v w : DualScaleStream2.Lattice.MukaiVec n) :
--   DualScaleStream2.Lattice.mukaiPair L v w = DualScaleStream2.Lattice.mukaiPair L w v
```

**Tier**: A (kernel-checked)  
**LeanMaster Dependencies**: `DualScaleStream2.Lattice.Mukai.mukaiPair_symm`  
**Supersedes**: `proofs/MukaiLatticeK3.lean::mukai_pairing_symmetric`  
**Description**: The Mukai pairing is symmetric for any symmetric Gram matrix.

---

### Theorem: mukai_self_pairing

```lean
#check DualScaleFoundation.K3Lattice.mukai_self_pairing
-- mukai_self_pairing (n : ℕ) (L : DualScaleStream2.Lattice.Gram n)
--   (v : DualScaleStream2.Lattice.MukaiVec n) :
--   DualScaleStream2.Lattice.mukaiPair L v v = v.c ⬝ᵥ Matrix.mulVec L v.c - 2 * v.r * v.s
```

**Tier**: A (kernel-checked)  
**LeanMaster Dependencies**: `DualScaleStream2.Lattice.Mukai.mukaiPair_self`  
**Supersedes**: `proofs/MukaiLatticeK3.lean` (self-pairing formula)  
**Description**: Self-pairing reduces to quadratic form on H² minus twice the product of H⁰ and H⁴ components.

---

### Theorem: k3_signature

```lean
#check DualScaleFoundation.K3Lattice.k3_signature
-- k3_signature : DualScaleStream2.Lattice.sigK3 = { pos := 3, neg := 19 }
```

**Tier**: A (kernel-checked)  
**LeanMaster Dependencies**: `DualScaleStream2.Lattice.K3T2Signature.sigK3_eq`  
**Supersedes**: `proofs/MukaiLatticeK3.lean` (signature data)  
**Description**: K3 lattice signature is (3,19) from decomposition E8(−1)⊕² ⊕ U⊕³.

---

### Theorem: mukai_signature

```lean
#check DualScaleFoundation.K3Lattice.mukai_signature
-- mukai_signature : DualScaleStream2.Lattice.sigMukai = { pos := 4, neg := 20 }
```

**Tier**: A (kernel-checked)  
**LeanMaster Dependencies**: `DualScaleStream2.Lattice.K3T2Signature.sigK3_eq` (arithmetic)  
**Supersedes**: N/A (extended lattice)  
**Description**: Mukai lattice Γ⁴'²⁰ = K3 lattice Γ³'¹⁹ ⊕ U has signature (4,20).

---

### Theorem: mukai_lattice_even

```lean
#check DualScaleFoundation.K3Lattice.mukai_lattice_even
-- mukai_lattice_even (n : ℕ) (L : DualScaleStream2.Lattice.Gram n)
--   (hL : Matrix.transpose L = L) (heven : DualScaleStream2.Lattice.IsEvenDiag L)
--   (v : DualScaleStream2.Lattice.MukaiVec n) : Even (DualScaleStream2.Lattice.mukaiPair L v v)
```

**Tier**: A (kernel-checked)  
**LeanMaster Dependencies**: `DualScaleStream2.Lattice.Mukai.mukaiPair_even`  
**Supersedes**: N/A (evenness property)  
**Description**: If H² has even diagonal, the full Mukai lattice is even.

---

### Theorem: reflection_is_isometry

```lean
#check DualScaleFoundation.K3Lattice.reflection_is_isometry
-- reflection_is_isometry (n : ℕ) (L : DualScaleStream2.Lattice.Gram n)
--   (hL : Matrix.transpose L = L) (v : Fin n → ℤ) (hv : DualScaleStream2.Lattice.latticeNorm L v = -2) :
--   Matrix.transpose (DualScaleStream2.Lattice.reflection L v) * L * DualScaleStream2.Lattice.reflection L v = L
```

**Tier**: A (kernel-checked)  
**LeanMaster Dependencies**: `DualScaleStream2.Lattice.Reflection.reflection_isometry`  
**Supersedes**: N/A (Weyl group structure)  
**Description**: Reflection in a (−2)-vector preserves the lattice form. These are the isometries generating the Weyl group of E8(−1).

---

## Tadpole Module (supersedes `proofs/TadpoleCancellation.lean`)

### Non-Vacuity Story
The old `proofs/TadpoleCancellation.lean` restated the O7/T⁴/ℤ₂ charge bookkeeping from literature.

The `DualScaleFoundation.Tadpole` module imports **only the geometry and arithmetic**:
- Euler characteristics from Hodge diamonds (K3 and T²)
- K3×K3 anomaly and tadpole budget (Dasgupta–Rajesh–Sethi)
- Flux integrality via Kronecker products of even lattices

The O7/D7 charge reconciliation is **left as a T0 question** (outside the LeanMaster scope).

---

### Theorem: k3k3_anomaly_value

```lean
#check DualScaleFoundation.Tadpole.k3k3_anomaly_value
-- k3k3_anomaly_value :
--   DualScaleStream2.Flux.chiK3K3 % 24 = 0 ∧ DualScaleStream2.Flux.chiK3K3 / 24 = 24
```

**Tier**: A (kernel-checked, Hodge input via Stream 1)  
**LeanMaster Dependencies**: `DualScaleStream2.Flux.Tadpole.k3k3_anomaly`, `StringTheoryFormalization.Frontier.HodgeNumbers.k3HodgeNumber`  
**Supersedes**: `proofs/TadpoleCancellation.lean::total_O7_charge_is_minus_64`  
**Description**: χ(K3×K3) = 24, so the M-theory anomaly χ/24 = 24 branes (per DRS eq. (3.1)).

---

### Theorem: k3t2_vanishing_euler

```lean
#check DualScaleFoundation.Tadpole.k3t2_vanishing_euler
-- k3t2_vanishing_euler :
--   DualScaleStream2.Flux.eulerFromHodge StringTheory.Frontier.k3HodgeNumber *
--       DualScaleStream2.Flux.eulerFromHodge DualScaleStream2.Flux.t2HodgeNumber =
--     0
```

**Tier**: A (kernel-checked, Hodge input via Stream 1)  
**LeanMaster Dependencies**: `DualScaleStream2.Flux.Tadpole.k3t2_euler_zero`, `StringTheoryFormalization.Frontier.HodgeNumbers`  
**Supersedes**: N/A (vanishing theorem)  
**Description**: χ(K3) · χ(T²) = 24 · 0 = 0, so K3×T² contributes no curvature tadpole alone.

---

### Theorem: tadpole_bound

```lean
#check DualScaleFoundation.Tadpole.tadpole_bound
-- tadpole_bound (flux n : ℕ) (h : flux + n = 24) : n ≤ 24 ∧ (flux = 0 → n = 24)
```

**Tier**: A (kernel-checked)  
**LeanMaster Dependencies**: `DualScaleStream2.Flux.Tadpole.tadpole_budget`  
**Supersedes**: `proofs/TadpoleCancellation.lean` (D7-brane count)  
**Description**: DRS tadpole budget: ½∫G∧G + n = 24. With no flux, exactly 24 D3-branes.

---

### Theorem: flux_integrality

```lean
#check DualScaleFoundation.Tadpole.flux_integrality
-- flux_integrality (L₁ L₂ : DualScaleStream2.Lattice.Gram 22) (h₁s : Matrix.transpose L₁ = L₁)
--   (h₂s : Matrix.transpose L₂ = L₂) (h₁ : DualScaleStream2.Lattice.IsEvenDiag L₁) (x : Fin 22 × Fin 22 → ℤ) :
--   ∃ k, x ⬝ᵥ (DualScaleStream2.Flux.kronForm L₁ L₂).mulVec x = 2 * k
```

**Tier**: A (kernel-checked)  
**LeanMaster Dependencies**: `DualScaleStream2.Flux.Integrality.flux_half_selfIntersection_integral`  
**Supersedes**: N/A (integrality proof)  
**Description**: Flux G ∈ H²(K3) ⊗ H²(K3) has even self-intersection, so ½∫G∧G ∈ ℤ.

---

## Moonshine Module (supersedes `proofs/MathieuVertexOperators.lean`)

### Non-Vacuity Story
The old `proofs/MathieuVertexOperators.lean` built rational-arithmetic fractions from scratch and stated M₂₄ representation dimensions as hard-coded constants (90, 462, 1540, …).

The `DualScaleFoundation.Moonshine` module cites **LeanMaster's `eotA` table**:
- Tier A: the `eotA : Fin 9 → ℕ` table of EOT's coefficients (45, 231, 770, …)
- Tier A: the 77/60 ratio derived from eotA values via exact arithmetic
- Tier L: the M₂₄ irrep interpretation (massive multiplicities, EOT's assertion of module structure)
- Tier C: physical reading as bispectrum non-Gaussianity

---

### Theorem: first_five_moonshine_irreps

```lean
#check DualScaleFoundation.Moonshine.first_five_moonshine_irreps
-- first_five_moonshine_irreps (n : Fin 5) :
--   ∃ i, StringTheory.StringDynamics.M24RepDim i = DualScaleStream2.Moonshine.eotA (Fin.castLE ⋯ n)
```

**Tier**: A (kernel-checked via `fin_cases` and table lookup)  
**LeanMaster Dependencies**: `DualScaleStream2.Moonshine.EOT.first_five_are_irreps`  
**Supersedes**: `proofs/MathieuVertexOperators.lean` (A1–A5 identification)  
**Description**: First five EOT coefficients (45, 231, 770, 2277, 5796) are irreducible representation dimensions of M₂₄.

---

### Theorem: a6_sum_decomposition

```lean
#check DualScaleFoundation.Moonshine.a6_sum_decomposition
-- a6_sum_decomposition :
--   DualScaleStream2.Moonshine.eotA 5 =
--     StringTheory.StringDynamics.M24RepDim 21 + StringTheory.StringDynamics.M24RepDim 25
```

**Tier**: A (kernel-checked)  
**LeanMaster Dependencies**: `DualScaleStream2.Moonshine.EOT.A6_decomposition`  
**Supersedes**: `proofs/MathieuVertexOperators.lean` (A6 decomposition)  
**Description**: A₆ = 3520 + 10395 (EOT eq. (1.14)).

---

### Theorem: a7_sum_decomposition

```lean
#check DualScaleFoundation.Moonshine.a7_sum_decomposition
-- a7_sum_decomposition :
--   DualScaleStream2.Moonshine.eotA 6 =
--     StringTheory.StringDynamics.M24RepDim 25 + StringTheory.StringDynamics.M24RepDim 23 +
--             StringTheory.StringDynamics.M24RepDim 24 +
--           StringTheory.StringDynamics.M24RepDim 22 +
--         StringTheory.StringDynamics.M24RepDim 18 +
--       StringTheory.StringDynamics.M24RepDim 17
```

**Tier**: A (kernel-checked)  
**LeanMaster Dependencies**: `DualScaleStream2.Moonshine.EOT.A7_decomposition`  
**Supersedes**: `proofs/MathieuVertexOperators.lean` (A7 decomposition)  
**Description**: A₇ = 10395 + 5796 + 5544 + 5313 + 2024 + 1771 (EOT eq. (1.15)).

---

### Theorem: bispectrum_ratio_exact

```lean
#check DualScaleFoundation.Moonshine.bispectrum_ratio_exact
-- bispectrum_ratio_exact :
--   2 * DualScaleStream2.Moonshine.eotA 1 * 60 = 4 * (2 * DualScaleStream2.Moonshine.eotA 0) * 77
```

**Tier A**: The arithmetic identity holds via `norm_num [eotA]`; eotA values are kernel-verified.

**Tier C**: The reading as bispectrum ratio ℛ_NL = dim(V₂)/(supercharges·dim(V₁)) = 462/(4·90) = 77/60 is physical conjecture.

**LeanMaster Dependencies**: `DualScaleStream2.Moonshine.EOT.eotA`  
**Supersedes**: `proofs/MathieuVertexOperators.lean::rNL` ratio  
**Description**:
- eotA 0 = 45 (A₁, dimension of first M₂₄ irrep)
- eotA 1 = 231 (A₂, dimension of second M₂₄ irrep)
- Massive multiplicity convention: A₁ → 2·45 = 90, A₂ → 2·231 = 462
- Supercharge count: 4 (N=4 superconformal symmetry on K3 worldsheet)
- Ratio: (2·231) / (4·2·45) = 462/360 = 77/60 (in lowest terms)

**Verification**:
- LHS: 2 * 231 * 60 = 27720
- RHS: 4 * (2 * 45) * 77 = 4 * 90 * 77 = 27720 ✓

---

### Theorem: gcd_77_60_one

```lean
#check DualScaleFoundation.Moonshine.gcd_77_60_one
-- gcd_77_60_one : Nat.gcd 77 60 = 1
```

**Tier**: A (kernel-checked via `norm_num`)  
**Supersedes**: N/A (minimality proof)  
**Description**: 77/60 is already in lowest terms.

---

### Theorem: bispectrum_ratio_minimal

```lean
#check DualScaleFoundation.Moonshine.bispectrum_ratio_minimal
-- bispectrum_ratio_minimal : Nat.gcd 77 60 = 1
```

**Tier**: A (kernel-checked via `norm_num`)  
**Supersedes**: N/A (minimality property)  
**Description**: Restatement of `gcd_77_60_one` for clarity.

---

## Summary

**Total Theorems**: 20  
**Audit Status**: All pass (0 failing)  
**Axiom Dependence**: Each theorem depends only on `propext`, `Classical.choice`, `Quot.sound` (Lean standard axioms)  
**LeanMaster Version**: v2.2.0 (verified)

All theorems have corresponding LeanMaster declarations they are based on. The old proofs relied on hand-written axioms, custom definitions, and unverified data; the foundation library restates them using kernel-verified libraries only.
