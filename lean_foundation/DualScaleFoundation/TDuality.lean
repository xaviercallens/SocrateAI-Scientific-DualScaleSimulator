/-
Supersedes proofs/BuscherRules.lean (see audit/lean_replacement_map.md)

T-duality content via DualScaleStream2.TDuality.ODD and DFT.GeneralizedMetric.
-/

import DualScaleStream2.TDuality.ODD
import DualScaleStream2.TDuality.Factorized
import DualScaleStream2.DFT.GeneralizedMetric

namespace DualScaleFoundation.TDuality

open DualScaleStream2.TDuality DualScaleStream2.DFT

/-! # T-Duality Involution

T-duality applied twice returns the identity. In the LeanMaster library:
- The factorized T-duality in direction `k ∈ Fin d` is `factorized k : Matrix (Charge d) (Charge d) ℤ`
- It is an involution: `factorized k * factorized k = 1`
-/

theorem tduality_involution_factorized (d : ℕ) (k : Fin d) :
    factorized k * factorized k = (1 : Matrix (Charge d) (Charge d) ℤ) :=
  factorized_mul_self k

/-! # Generalized Metric Involution

In Double Field Theory, the generalized metric H satisfies the constraint (ηH)² = 1,
which expresses T-duality's involutive property in the continuous (real) setting.
-/

theorem generalized_metric_constraint (d : ℕ) (G B : Matrix (Fin d) (Fin d) ℝ)
    (hG : IsUnit G.det) :
    (etaR d * genMetric G B) * (etaR d * genMetric G B) = 1 :=
  etaR_genMetric_sq G B hG

/-! # Example: T-duality is an O(d,d;ℤ) group element

The full T-duality (swapping all momenta and windings) is itself in O(d,d;ℤ).
-/

theorem tduality_is_odd (d : ℕ) : IsODD (eta d) :=
  eta_isODD

example (d : ℕ) : IsODD (eta d) := eta_isODD

/-! # Structure-preserving property

All factorized T-dualities lie in O(d,d;ℤ).
-/

theorem factorized_is_odd (d : ℕ) (k : Fin d) : IsODD (factorized k) :=
  factorized_isODD k

end DualScaleFoundation.TDuality
