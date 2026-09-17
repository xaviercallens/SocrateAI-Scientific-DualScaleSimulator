/-
Supersedes proofs/TadpoleCancellation.lean (see audit/lean_replacement_map.md)

Tadpole counting and flux integrality via DualScaleStream2.Flux.
NOTE: Does NOT restate the O7/T^4/Z2 charge bookkeeping — it is unresolved (T0 question).
-/

import DualScaleStream2.Flux.Tadpole
import DualScaleStream2.Flux.Integrality

namespace DualScaleFoundation.Tadpole

open DualScaleStream2.Flux DualScaleStream2.Lattice Matrix

/-! # K3 × K3 Anomaly

The M-theory anomaly χ/24 for K3 × K3 equals 24 (Dasgupta–Rajesh–Sethi).
This is the tadpole contribution from the geometry alone, before flux and branes.
-/

theorem k3k3_anomaly_value : chiK3K3 % 24 = 0 ∧ chiK3K3 / 24 = 24 :=
  k3k3_anomaly

/-! # K3 × T² Euler Characteristic

K3 × T² contributes zero to the Euler characteristic because χ(T²) = 0.
-/

theorem k3t2_vanishing_euler : eulerFromHodge StringTheory.Frontier.k3HodgeNumber *
    eulerFromHodge t2HodgeNumber = 0 :=
  k3t2_euler_zero

/-! # Tadpole Budget

The tadpole cancellation condition is: (flux contribution) + (D3-brane number) = 24.
With no flux, exactly 24 D3-branes are required.
-/

theorem tadpole_bound (flux n : ℕ) (h : flux + n = 24) : n ≤ 24 ∧ (flux = 0 → n = 24) :=
  tadpole_budget flux n h

/-! # Flux Contribution is Integral

When G is a flux in H²(K3) ⊗ H²(K3), its self-intersection ∫ G ∧ G is even,
so the flux contribution ½ ∫ G ∧ G is an integer.
This is proved using the Kronecker product structure and the fact that the K3 lattice is even.
-/

theorem flux_integrality (L₁ L₂ : Gram 22) (h₁s : L₁.transpose = L₁) (h₂s : L₂.transpose = L₂)
    (h₁ : IsEvenDiag L₁) (x : Fin 22 × Fin 22 → ℤ) :
    ∃ k : ℤ, x ⬝ᵥ (kronForm L₁ L₂ *ᵥ x) = 2 * k :=
  flux_half_selfIntersection_integral L₁ L₂ h₁s h₂s h₁ x

end DualScaleFoundation.Tadpole
