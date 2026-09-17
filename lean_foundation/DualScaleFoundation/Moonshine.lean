/-
Supersedes proofs/MathieuVertexOperators.lean (see audit/lean_replacement_map.md)

Moonshine coefficients and the bispectrum ratio via DualScaleStream2.Moonshine.EOT.
-/

import DualScaleStream2.Moonshine.EOT

namespace DualScaleFoundation.Moonshine

open DualScaleStream2.Moonshine StringTheory.StringDynamics

/-! # Mathieu Moonshine Coefficients

The elliptic genus of K3 carries an action of the sporadic group M₂₄.
The first five Fourier coefficients (normalized per EOT eq. 1.12) are irreducible
representation dimensions of M₂₄, with A₁ = 45, A₂ = 231, A₃ = 770, etc.
-/

theorem first_five_moonshine_irreps :
    ∀ n : Fin 5, ∃ i : Fin 26, M24RepDim i = eotA (Fin.castLE (by norm_num) n) :=
  first_five_are_irreps

/-! # Higher Coefficients as Sums

A₆ and A₇ decompose as sums of irreducible dimensions.
-/

theorem a6_sum_decomposition :
    eotA 5 = M24RepDim 21 + M24RepDim 25 :=
  A6_decomposition

theorem a7_sum_decomposition :
    eotA 6 = M24RepDim 25 + M24RepDim 23 + M24RepDim 24 + M24RepDim 22 +
      M24RepDim 18 + M24RepDim 17 :=
  A7_decomposition

/-! # Bispectrum Ratio (Tier A arithmetic; Tier C physical interpretation)

The vertex operator algebra's 3-point function Ward identities yield a ratio
of representation dimensions. With the massive-multiplicity convention (each
irrep appears twice: 45 → 90, 231 → 462), and 4 supercharges, the ratio is:

  ℛ_NL = dim(V₂) / (supercharges × dim(V₁)) = (2·A₂) / (4 × 2·A₁)

Cross-multiplying to avoid ℕ division:

  (2·A₂) × 60 = (4 × 2·A₁) × 77

This proves the exact arithmetic identity 77/60 in lowest terms.

[Tier A] The identity holds as integer arithmetic derived from eotA values.
[Tier C] Its reading as the bispectrum ratio ℛ_NL is conjectural.
-/

theorem bispectrum_ratio_exact :
    (2 * eotA 1) * 60 = (4 * (2 * eotA 0)) * 77 := by
  norm_num [eotA]

lemma gcd_77_60_one : Nat.gcd 77 60 = 1 := by
  norm_num

/-! # Verification

The ratio 77/60 is already in lowest terms.
-/

theorem bispectrum_ratio_minimal : Nat.gcd 77 60 = 1 :=
  gcd_77_60_one

end DualScaleFoundation.Moonshine
