/-
Supersedes proofs/MukaiLatticeK3.lean (see audit/lean_replacement_map.md)

Mukai lattice and K3 surface geometry via DualScaleStream2.Lattice.
-/

import DualScaleStream2.Lattice.Mukai
import DualScaleStream2.Lattice.K3T2Signature
import DualScaleStream2.Lattice.Reflection

namespace DualScaleFoundation.K3Lattice

open DualScaleStream2.Lattice Matrix

/-! # Mukai Pairing

The Mukai pairing on the cohomology lattice of a K3 surface is symmetric.
-/

theorem mukai_pairing_symmetric (n : ℕ) (L : Gram n) (hL : L.transpose = L)
    (v w : MukaiVec n) :
    mukaiPair L v w = mukaiPair L w v :=
  mukaiPair_symm L hL v w

theorem mukai_self_pairing (n : ℕ) (L : Gram n) (v : MukaiVec n) :
    mukaiPair L v v = v.c ⬝ᵥ (L *ᵥ v.c) - 2 * v.r * v.s := by
  unfold mukaiPair
  exact mukaiPair_self L v

/-! # K3 Signature

The K3 lattice decomposes as E8(−1)⊕2 ⊕ U⊕3, yielding signature (3,19).
-/

theorem k3_signature : sigK3 = ⟨3, 19⟩ :=
  sigK3_eq

theorem mukai_signature : sigMukai = ⟨4, 20⟩ := by
  unfold sigMukai sigK3 sigE8Neg sigU
  rfl

/-! # Mukai Lattice Evenness

The Mukai lattice is an even lattice whenever the K3 lattice (the H² component) is even.
-/

theorem mukai_lattice_even (n : ℕ) (L : Gram n) (hL : L.transpose = L) (heven : IsEvenDiag L)
    (v : MukaiVec n) :
    Even (mukaiPair L v v) :=
  mukaiPair_even L hL heven v

/-! # Reflection Isometries

Reflections in (−2)-vectors are isometries of the lattice.
-/

theorem reflection_is_isometry (n : ℕ) (L : Gram n) (hL : L.transpose = L)
    (v : Fin n → ℤ) (hv : latticeNorm L v = -2) :
    (reflection L v).transpose * L * reflection L v = L :=
  reflection_isometry L hL v hv

end DualScaleFoundation.K3Lattice
