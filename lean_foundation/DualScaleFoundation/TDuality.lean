-- Supersedes proofs/BuscherRules.lean and proofs/LeanscratchDB/DoubleScaleT2.lean (see audit/lean_replacement_map.md)
-- T-duality foundations via kernel-verified LeanMaster declarations

import DualScaleStream2.TDuality.ODD
import DualScaleStream2.DFT.GeneralizedMetric
import DualScaleStream2.TDuality.Factorized

namespace DualScaleFoundation

open DualScaleStream2.TDuality DualScaleStream2.DFT Matrix

/-!
# T-Duality (Buscher Transformation)

Tier A: theorems directly cite kernel-verified LeanMaster declarations.
-/

/-- Theta-shift is ODD: kernel-verified by `thetaShift_isODD`. -/
theorem theta_shift_is_odd (d : ℕ) (Θ : Matrix (Fin d) (Fin d) ℤ) (hΘ : Θᵀ = -Θ) :
    IsODD (thetaShift Θ) :=
  thetaShift_isODD Θ hΘ

/-- Basis changes preserve ODD: kernel-verified by `basisChange_isODD`. -/
theorem basis_change_preserves_odd (d : ℕ) (A B : Matrix (Fin d) (Fin d) ℤ) (h : Aᵀ * B = 1) :
    IsODD (basisChange A B) :=
  basisChange_isODD A B h

/-- Product of ODD matrices is ODD: kernel-verified by `isODD_mul`. -/
theorem odd_mul_property (d : ℕ) (g h : Matrix (Charge d) (Charge d) ℤ)
    (hg : IsODD g) (hh : IsODD h) : IsODD (g * h) :=
  isODD_mul g h hg hh

/-- T-duality inverts the metric at B=0: kernel-verified by `tduality_inverts_metric`. -/
theorem tduality_metric_inversion (d : ℕ) (G : Matrix (Fin d) (Fin d) ℝ) (hG : IsUnit G.det) :
    etaR d * genMetric G 0 * etaR d = genMetric G⁻¹ 0 :=
  tduality_inverts_metric G hG

/-- Factorized T-duality along the `k`-th circle is an O(d,d;ℤ) element: `factorized_isODD`. -/
theorem factorized_tduality_is_odd (d : ℕ) (k : Fin d) : IsODD (factorized k) :=
  factorized_isODD k

/-- Factorized T-duality along one circle squares to the identity: `factorized_mul_self`.
This is the kernel-checked replacement for `buscher_involution` in proofs/BuscherRules.lean,
which rested on five project axioms that are jointly inconsistent. Scope: the integer O(d,d)
charge action, not the field-level Buscher rules for metric, B-field and dilaton. -/
theorem factorized_tduality_involution (d : ℕ) (k : Fin d) : factorized k * factorized k = 1 :=
  factorized_mul_self k

end DualScaleFoundation
