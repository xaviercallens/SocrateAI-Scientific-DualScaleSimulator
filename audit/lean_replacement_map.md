# Lean Replacement Map: proofs/*.lean → lean_foundation/DualScaleFoundation/

This document maps each old Lean proof file to its replacement in the verified foundation library.

## Mapping

| Old Proof File | New Foundation Module | LeanMaster Base | Key Changes |
|---|---|---|---|
| `proofs/BuscherRules.lean` | `lean_foundation/DualScaleFoundation/TDuality.lean` | `DualScaleStream2.TDuality.{ODD,Factorized}` + `DualScaleStream2.DFT.GeneralizedMetric` | Removed 5 hand-written axioms (`inv_inv`, `buscher_cross_inv`, etc.); use kernel-verified O(d,d;ℤ) matrices instead |
| `proofs/MukaiLatticeK3.lean` | `lean_foundation/DualScaleFoundation/K3Lattice.lean` | `DualScaleStream2.Lattice.{Mukai,K3T2Signature,Reflection}` | Replaced custom `MukaiVector`/`mukaiPairing` definitions with LeanMaster's verified versions; verified signature (3,19) for K3 |
| `proofs/TadpoleCancellation.lean` | `lean_foundation/DualScaleFoundation/Tadpole.lean` | `DualScaleStream2.Flux.{Tadpole,Integrality}` | Retained K3×K3 anomaly and tadpole budget (arithmetic); omitted O7/T⁴/ℤ₂ charge bookkeeping (T0 question) |
| `proofs/MathieuVertexOperators.lean` | `lean_foundation/DualScaleFoundation/Moonshine.lean` | `DualScaleStream2.Moonshine.EOT` | Replaced hard-coded representation dimensions (90, 462, …) with LeanMaster's `eotA` table; derived 77/60 ratio from `eotA` values |

---

## Coverage

**Files replaced**: 4 (out of 14 proofs/*.lean total)
- TDuality involutions and O(d,d;ℤ) group closure
- K3 lattice geometry and Mukai pairings
- Tadpole cancellation (geometry + flux; excludes charge accounting)
- Moonshine coefficients and bispectrum ratio

**Files not yet addressed**:
- `ExportAxioms.lean` (axiom enumeration; superseded by axiom audit)
- `DbraneBoundaryStates.lean` (boundary state formalism)
- `FourierMukai.lean` (derived category transforms)
- `GysinBEMSequence.lean` (cohomology sequences)
- `KummerOrbifoldResolution.lean` (Kummer surface)
- `KummerTDAAnomalyCertification.lean` (anomaly polynomial)
- `FluxVacuumDecayCTheorem.lean` (decay rates)
- `FTheoryCosmology.lean` (F-theory duality)
- `SwamplandDistanceConjecture.lean` (moduli geometry)
- `TachyonCondensationKTheory.lean` (open string dynamics)

These may be addressed in future iterations of the foundation bridge.

---

## Non-Vacuity Summary

### BuscherRules → TDuality
- **Old**: 5 axioms asserting field properties + TargetSpace structure
- **New**: O(d,d;ℤ) matrices, kernel-verified closure laws
- **Win**: Proofs via LeanMaster declarations; no new axioms beyond Lean's standard three

### MukaiLatticeK3 → K3Lattice
- **Old**: Custom MukaiVector, mukaiPairing; signature hardcoded
- **New**: LeanMaster's lattice hierarchy; verified E8(−1) ⊕ U decomposition
- **Win**: Signature (3,19) derived from component signatures, not asserted

### TadpoleCancellation → Tadpole
- **Old**: 16 fixed points + O7/D7 charge reconciliation
- **New**: K3×K3 anomaly (χ=24) and DRS budget; O7/D7 deferred
- **Win**: Anomaly is geometric (Hodge diamonds); charge balance is T0 question

### MathieuVertexOperators → Moonshine
- **Old**: Rational arithmetic library from scratch; eotA as magic numbers
- **New**: LeanMaster's eotA table + Burnside check against M24RepDim sum
- **Win**: 77/60 ratio proven from table values, not re-verified

---

## Axiom Audit Results

```
20 theorems audited, 0 failing.
Dependency set for all theorems: {propext, Classical.choice, Quot.sound}
```

**Control Test**:
- Before removal of `NegControl.lean`: 21 theorems, 1 failing (`sorryAx`)
- After removal: 20 theorems, 0 failing

Negative control theorem (false, with `sorry`) was correctly flagged by axiom audit as `FAIL` (contains `sorryAx`).

---

## Links

- **LeanMaster Documentation**: `/home/callensxavier_gmail_com/SocrateAI-Scientific-Agora-LeanMaster/docs/VERIFIED_FOUNDATION.md` (v2.2.0)
- **Verification Certificate**: LeanMaster build log `/mnt/disks/disk-socrateai-local-1/leanmaster-v2.2.0/build_v2.2.0.log` (EXIT 0)
- **Foundation Statements**: `lean_foundation/STATEMENTS_FOR_REVIEW.md` (this session)
