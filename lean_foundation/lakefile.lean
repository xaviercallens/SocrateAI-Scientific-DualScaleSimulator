import Lake
open Lake DSL

/-!
DualScale Foundation: downstream project building ON TOP of LeanMaster v2.2.0.
Restates content from proofs/*.lean using kernel-verified libraries.

Depends on PINNED LeanMaster v2.2.0 clone and reuses its built Mathlib.
-/

package «dualscale-foundation» where
  packagesDir := "/mnt/disks/disk-socrateai-local-1/leanmaster/lake/packages"
  leanOptions := #[⟨`maxHeartbeats, (1000000 : Nat)⟩]

-- Pinned v2.2.0 clone (verified v2.2.0 tag, EXIT 0 build log)
require «SocrateAI-Scientific-Agora-LeanMaster» from
  "/mnt/disks/disk-socrateai-local-1/leanmaster-v2.2.0"

@[default_target]
lean_lib «DualScaleFoundation» where
  roots := #[`DualScaleFoundation]
