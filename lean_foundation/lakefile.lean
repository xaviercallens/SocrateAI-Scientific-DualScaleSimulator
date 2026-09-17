import Lake
open Lake DSL

/-!
DualScale Foundation: downstream project built on LeanMaster v2.2.0 (kernel-verified).
Restates the content of proofs/*.lean using LeanMaster's verified libraries.

Verified 2026-09-17: LeanMaster v2.2.0 at /mnt/disks/disk-socrateai-local-1/leanmaster-v2.2.0
(git describe --tags: v2.2.0, build log exits 0).
-/

package «dualscale-foundation» where
  -- Reuse LeanMaster's packages (Mathlib etc.)
  packagesDir := "/mnt/disks/disk-socrateai-local-1/leanmaster/lake/packages"
  leanOptions := #[⟨`maxHeartbeats, (1000000 : Nat)⟩]

-- Depend on PINNED v2.2.0 clone (read-only; another session may edit live checkout)
require «SocrateAI-Scientific-Agora-LeanMaster» from
  "/mnt/disks/disk-socrateai-local-1/leanmaster-v2.2.0"

@[default_target]
lean_lib «DualScaleFoundation» where
  roots := #[`DualScaleFoundation]
