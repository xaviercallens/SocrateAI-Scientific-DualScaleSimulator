# DualScaleFoundation Theorem Statements for Review

**Date**: 2026-09-17  
**Version**: DualScaleFoundation downstream on LeanMaster v2.2.0  
**Total Theorems**: 17  
**Axiom Audit Result**: 17 theorems, 0 failing; all depend only on {propext, Classical.choice, Quot.sound}

---

## TDuality Module

### DualScaleFoundation.theta_shift_is_odd
```
theorem theta_shift_is_odd (d : ℕ) (Θ : Matrix (Fin d) (Fin d) ℤ) (hΘ : Θᵀ = -Θ) :
    IsODD (thetaShift Θ)
```
**Supersedes**: `proofs/TDuality.lean`  
**Depends on**: `DualScaleStream2.TDuality.thetaShift_isODD`  
**Tier**: A — Kernel-verified via LeanMaster  
**Status**: ✓ Proved

### DualScaleFoundation.basis_change_preserves_odd
```
theorem basis_change_preserves_odd (d : ℕ) (A B : Matrix (Fin d) (Fin d) ℤ) (h : Aᵀ * B = 1) :
    IsODD (basisChange A B)
```
**Supersedes**: `proofs/TDuality.lean`  
**Depends on**: `DualScaleStream2.TDuality.basisChange_isODD`  
**Tier**: A — Kernel-verified via LeanMaster  
**Status**: ✓ Proved

### DualScaleFoundation.odd_mul_property
```
theorem odd_mul_property (d : ℕ) (g h : Matrix (Charge d) (Charge d) ℤ)
    (hg : IsODD g) (hh : IsODD h) : IsODD (g * h)
```
**Supersedes**: `proofs/TDuality.lean`  
**Depends on**: `DualScaleStream2.TDuality.isODD_mul`  
**Tier**: A — Kernel-verified via LeanMaster  
**Status**: ✓ Proved

### DualScaleFoundation.tduality_metric_inversion
```
theorem tduality_metric_inversion (d : ℕ) (G : Matrix (Fin d) (Fin d) ℝ) (hG : IsUnit G.det) :
    etaR d * genMetric G 0 * etaR d = genMetric G⁻¹ 0
```
**Supersedes**: `proofs/TDuality.lean`  
**Depends on**: `DualScaleStream2.DFT.tduality_inverts_metric`  
**Tier**: A — Kernel-verified via LeanMaster  
**Status**: ✓ Proved

---

## K3Lattice Module

### DualScaleFoundation.k3_signature_is_three_nineteen
```
theorem k3_signature_is_three_nineteen : sigK3 = ⟨3, 19⟩
```
**Supersedes**: `proofs/K3Lattice.lean`  
**Depends on**: `DualScaleStream2.Lattice.K3T2Signature.sigK3_eq`  
**Tier**: A — Kernel-verified via LeanMaster  
**Status**: ✓ Proved

### DualScaleFoundation.mukai_signature_is_four_twenty
```
theorem mukai_signature_is_four_twenty : sigMukai = ⟨4, 20⟩
```
**Supersedes**: `proofs/K3Lattice.lean`  
**Depends on**: `DualScaleStream2.Lattice.K3T2Signature.sigMukai_eq`  
**Tier**: A — Kernel-verified via LeanMaster  
**Status**: ✓ Proved

### DualScaleFoundation.k3t2_signature
```
theorem k3t2_signature : sigK3T2 = ⟨6, 22⟩
```
**Supersedes**: `proofs/K3Lattice.lean`  
**Depends on**: `DualScaleStream2.Lattice.K3T2Signature.sigK3T2_eq`  
**Tier**: A — Kernel-verified via LeanMaster  
**Status**: ✓ Proved

### DualScaleFoundation.k3_rank
```
theorem k3_rank : sigK3.rank = 22
```
**Supersedes**: `proofs/K3Lattice.lean`  
**Depends on**: `DualScaleStream2.Lattice.K3T2Signature.rank_K3`  
**Tier**: A — Kernel-verified via LeanMaster  
**Status**: ✓ Proved

---

## Tadpole Module

### DualScaleFoundation.k3_cross_k3_euler_over_24
```
theorem k3_cross_k3_euler_over_24 : chiK3K3 % 24 = 0 ∧ chiK3K3 / 24 = 24
```
**Supersedes**: `proofs/Tadpole.lean`  
**Depends on**: `DualScaleStream2.Flux.Tadpole.k3k3_anomaly`  
**Tier**: A — Kernel-verified via LeanMaster  
**Status**: ✓ Proved

### DualScaleFoundation.k3t2_euler_is_zero
```
theorem k3t2_euler_is_zero : eulerFromHodge StringTheory.Frontier.k3HodgeNumber *
    eulerFromHodge t2HodgeNumber = 0
```
**Supersedes**: `proofs/Tadpole.lean`  
**Depends on**: `DualScaleStream2.Flux.Tadpole.k3t2_euler_zero`  
**Tier**: A — Kernel-verified via LeanMaster  
**Status**: ✓ Proved

### DualScaleFoundation.tadpole_conservation
```
theorem tadpole_conservation (flux n : ℕ) (h : flux + n = 24) : n ≤ 24 ∧ (flux = 0 → n = 24)
```
**Supersedes**: `proofs/Tadpole.lean`  
**Depends on**: `DualScaleStream2.Flux.Tadpole.tadpole_budget`  
**Tier**: A — Kernel-verified via LeanMaster  
**Status**: ✓ Proved

---

## Moonshine Module

### DualScaleFoundation.eot_values
```
theorem eot_values : eotA = ![45, 231, 770, 2277, 5796, 13915, 30843, 65550, 132825]
```
**Supersedes**: `proofs/Moonshine.lean`  
**Depends on**: `DualScaleStream2.Moonshine.EOT.eotA`  
**Tier**: A — Kernel-verified via LeanMaster  
**Status**: ✓ Proved

### DualScaleFoundation.eot_a1_value
```
theorem eot_a1_value : eotA 0 = 45
```
**Supersedes**: `proofs/Moonshine.lean`  
**Depends on**: `DualScaleStream2.Moonshine.EOT.eotA`  
**Tier**: A — Arithmetic (decided)  
**Status**: ✓ Proved

### DualScaleFoundation.eot_a4_value
```
theorem eot_a4_value : eotA 3 = 2277
```
**Supersedes**: `proofs/Moonshine.lean`  
**Depends on**: `DualScaleStream2.Moonshine.EOT.eotA`  
**Tier**: A — Arithmetic (decided)  
**Status**: ✓ Proved

### DualScaleFoundation.first_five_m24_irreps
```
theorem first_five_m24_irreps :
    ∀ n : Fin 5, ∃ i : Fin 26, StringTheory.StringDynamics.M24RepDim i = eotA (Fin.castLE (by norm_num) n)
```
**Supersedes**: `proofs/Moonshine.lean`  
**Depends on**: `DualScaleStream2.Moonshine.EOT.first_five_are_irreps`  
**Tier**: A — Kernel-verified via LeanMaster  
**Status**: ✓ Proved

### DualScaleFoundation.a6_not_m24_irrep
```
theorem a6_not_m24_irrep : ∀ i : Fin 26, StringTheory.StringDynamics.M24RepDim i ≠ eotA 5
```
**Supersedes**: `proofs/Moonshine.lean`  
**Depends on**: `DualScaleStream2.Moonshine.EOT.A6_not_irrep`  
**Tier**: A — Kernel-verified via LeanMaster  
**Status**: ✓ Proved

### DualScaleFoundation.a6_decomposition
```
theorem a6_decomposition :
    eotA 5 = StringTheory.StringDynamics.M24RepDim 21 + StringTheory.StringDynamics.M24RepDim 25
```
**Supersedes**: `proofs/Moonshine.lean`  
**Depends on**: `DualScaleStream2.Moonshine.EOT.A6_decomposition`  
**Tier**: A — Kernel-verified via LeanMaster  
**Status**: ✓ Proved

---

## Summary

- **Total Theorems**: 17
- **All Tier A** (Kernel-verified via LeanMaster v2.2.0)
- **Build**: ✓ 3664 jobs
- **Axiom Audit**: ✓ 17 theorems, 0 failing
- **Sorries**: 0
- **Statement Lock**: Pending review (T1 gate)

