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

Axiom status is established by running LeanMaster's tools/axiom_audit.py on this library (see STATEMENTS_FOR_REVIEW.md), not by this comment.
-/

end DualScaleFoundation
