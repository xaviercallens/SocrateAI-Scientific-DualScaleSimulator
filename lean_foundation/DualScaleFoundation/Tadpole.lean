-- Supersedes proofs/Tadpole.lean (see audit/lean_replacement_map.md)
-- Tadpole cancellation and flux via kernel-verified LeanMaster declarations

import DualScaleStream2.Flux.Tadpole

namespace DualScaleFoundation

open DualScaleStream2.Flux

/-!
# Tadpole Cancellation and Flux Constraints

Tier A: theorems cite kernel-verified LeanMaster declarations.

NOTE: O7/T⁴/Z2 charge bookkeeping is unresolved (T0); not restated here.
-/

/-- K3 × K3 Euler characteristic: χ(K3 × K3) = 24 × 24 = 576, so χ/24 = 24.
    Kernel-verified: `k3k3_anomaly` from Flux.Tadpole. -/
theorem k3_cross_k3_euler_over_24 : chiK3K3 % 24 = 0 ∧ chiK3K3 / 24 = 24 :=
  k3k3_anomaly

/-- K3 × T² has zero Euler characteristic (topological constraint).
    Kernel-verified: `k3t2_euler_zero` from Flux.Tadpole. -/
theorem k3t2_euler_is_zero : eulerFromHodge StringTheory.Frontier.k3HodgeNumber *
    eulerFromHodge t2HodgeNumber = 0 :=
  k3t2_euler_zero

/-- Tadpole budget constraint: when flux + n = 24, then n ≤ 24 and if flux=0 then n=24.
    Kernel-verified: `tadpole_budget` from Flux.Tadpole. -/
theorem tadpole_conservation (flux n : ℕ) (h : flux + n = 24) : n ≤ 24 ∧ (flux = 0 → n = 24) :=
  tadpole_budget flux n h

end DualScaleFoundation
