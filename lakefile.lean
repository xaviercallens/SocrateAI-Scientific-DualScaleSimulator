import Lake
open Lake DSL

package dual_scale {
  -- add package configuration options here
}

require mathlib from git
  "https://github.com/leanprover-community/mathlib4.git"

@[default_target]
lean_lib MathieuVertexOperators {
  srcDir := "proofs"
}

@[default_target]
lean_lib TadpoleCancellation {
  srcDir := "proofs"
}

@[default_target]
lean_lib KummerTDAAnomalyCertification {
  srcDir := "proofs"
}

@[default_target]
lean_lib SwamplandDistanceConjecture {
  srcDir := "proofs"
}

@[default_target]
lean_lib TachyonCondensationKTheory {
  srcDir := "proofs"
}

@[default_target]
lean_lib FluxVacuumDecayCTheorem {
  srcDir := "proofs"
}

@[default_target]
lean_lib BuscherRules {
  srcDir := "proofs"
}

@[default_target]
lean_lib FourierMukai {
  srcDir := "proofs"
}

@[default_target]
lean_lib GysinBEMSequence {
  srcDir := "proofs"
}

@[default_target]
lean_lib MukaiLatticeK3 {
  srcDir := "proofs"
}

@[default_target]
lean_lib KummerOrbifoldResolution {
  srcDir := "proofs"
}

@[default_target]
lean_lib DbraneBoundaryStates {
  srcDir := "proofs"
}

@[default_target]
lean_lib LeanscratchDB.TopologicalTDuality {
  srcDir := "proofs"
}

@[default_target]
lean_lib LeanscratchDB.EmergentCosmology {
  srcDir := "proofs"
}

@[default_target]
lean_lib LeanscratchDB.QuantumEntanglement {
  srcDir := "proofs"
}

@[default_target]
lean_lib LeanscratchDB.DoubleScaleT2 {
  srcDir := "proofs"
}

@[default_target]
lean_lib LeanscratchDB.HoloAlg {
  srcDir := "proofs"
}

@[default_target]
lean_exe export_axioms {
  root := `ExportAxioms
  srcDir := "proofs"
}
