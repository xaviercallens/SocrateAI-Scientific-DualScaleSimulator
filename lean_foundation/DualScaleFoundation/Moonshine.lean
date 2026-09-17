-- Supersedes proofs/MathieuVertexOperators.lean (see audit/lean_replacement_map.md)
-- Moonshine arithmetic via kernel-verified LeanMaster declarations

import DualScaleStream2.Moonshine.EOT

namespace DualScaleFoundation

open DualScaleStream2.Moonshine

/-!
# Moonshine Arithmetic and M24 Representation Theory

Tier A: arithmetic identities derived from kernel-verified EOT A_n values.
Tier C: physical interpretation (bispectrum ratio) is conjectural.
-/

/-- First nine EOT A_n values indexed by Fin 9. -/
theorem eot_values : eotA = ![45, 231, 770, 2277, 5796, 13915, 30843, 65550, 132825] := by
  rfl

/-- A₁ (index 0) has value 45: kernel-verified from `eotA` definition. -/
theorem eot_a1_value : eotA 0 = 45 := by
  decide

/-- A₄ (index 3) has value 2277: kernel-verified from `eotA` definition. -/
theorem eot_a4_value : eotA 3 = 2277 := by
  decide

/-- First five EOT values are irreducible M24 representations.
    Kernel-verified: `first_five_are_irreps` from Moonshine.EOT. -/
theorem first_five_m24_irreps :
    ∀ n : Fin 5, ∃ i : Fin 26, StringTheory.StringDynamics.M24RepDim i = eotA (Fin.castLE (by norm_num) n) :=
  first_five_are_irreps

/-- A₆ (index 5) is not a single irreducible M24 representation.
    Kernel-verified: `A6_not_irrep` from Moonshine.EOT. -/
theorem a6_not_m24_irrep : ∀ i : Fin 26, StringTheory.StringDynamics.M24RepDim i ≠ eotA 5 :=
  A6_not_irrep

/-- A₆ decomposes as sum of two M24 irreps (EOT eq. 1.14).
    Kernel-verified: `A6_decomposition` from Moonshine.EOT. -/
theorem a6_decomposition :
    eotA 5 = StringTheory.StringDynamics.M24RepDim 21 + StringTheory.StringDynamics.M24RepDim 25 :=
  A6_decomposition

/-- The ratio 77/60 of `proofs/MathieuVertexOperators.lean`, rebuilt from the EOT values instead of
retyped constants: massive multiplicities are `2 * eotA n`, so 462 = 2·A₂ and 90 = 2·A₁ and
462 / (4 · 90) = 77/60. Tier A as arithmetic only; reading it as a primordial bispectrum
ratio is a conjecture (tier C) with no derivation in this project. -/
theorem r_bps_from_eot : ((2 * eotA 1 : ℕ) : ℚ) / (4 * (2 * eotA 0 : ℕ)) = 77 / 60 := by
  simp [eotA]
  norm_num

end DualScaleFoundation
