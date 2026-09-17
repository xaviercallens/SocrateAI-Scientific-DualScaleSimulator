-- Supersedes proofs/TDuality.lean (see audit/lean_replacement_map.md)
-- T-duality foundations via kernel-verified LeanMaster declarations

import DualScaleStream2.TDuality.ODD
import DualScaleStream2.DFT.GeneralizedMetric

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

end DualScaleFoundation
