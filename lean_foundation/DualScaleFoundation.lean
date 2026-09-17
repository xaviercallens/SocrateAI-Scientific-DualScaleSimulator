-- DualScale Foundation Library
-- Downstream project building on LeanMaster v2.2.0 kernel-verified libraries

import DualScaleFoundation.TDuality
import DualScaleFoundation.K3Lattice
import DualScaleFoundation.Tadpole
import DualScaleFoundation.Moonshine

namespace DualScaleFoundation

/-!
# DualScale Foundation

A downstream Lean project that restates content from proofs/*.lean non-vacuously,
using only kernel-verified declarations from LeanMaster v2.2.0.

Each module restates specific string-theory content:
- TDuality: Buscher transformation and ODD symmetry
- K3Lattice: K3 surface and Mukai lattice structure
- Tadpole: Tadpole cancellation constraints
- Moonshine: Moonshine arithmetic and M24 representation theory

Verified: depends only on Lean's three standard axioms (propext, Classical.choice, Quot.sound).
-/

end DualScaleFoundation
