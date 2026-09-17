# Paper facts sheet (single source of truth for manuscript rewrites)

Every fact below was established by a command run on 2026-09-17 in this repository or in the
LeanMaster clone at tag v2.2.0. Writers may use ONLY these facts plus text already in the paper and
files they open themselves. Anything else: leave the sentence's scope smaller, or tag it.

Tiers (Mathesis): A kernel-checked with only propext/Classical.choice/Quot.sound AND adequate statement ·
B exact arithmetic with a negative control · L literature (quoted source) · C conjecture/interpretation ·
X exploratory numerics.

## F1. Local Lean files in `proofs/` (superseded)
- `proofs/BuscherRules.lean` declares 5 global axioms (`inv_inv`, `buscher_cross_inv`, `buscher_gmunu_inv`,
  `buscher_Bmunu_inv`, `Phi_inv`). `inv_inv` instantiated at `Nat` with x = 2 proves `False`
  (checked with `lake env lean`; `#print axioms` = [propext, inv_inv]). Hence `buscher_involution`
  certifies nothing. Its dilaton shift also omits the factor 1/2, and its sign convention differs from
  `proofs/LeanscratchDB/DoubleScaleT2.lean`.
- LeanMaster axiom audit over the 18 local library modules: 83 theorems audited, 5 failing
  (`buscher_involution`; HoloAlg `topology_exchange`, `tadpole_cancellation`, `flop_charge_conservation`,
  `gravitino_anomaly_cancellation`, all resting on project axioms).
- The other 78 pass the axiom audit but many are true by construction of their definitions:
  SDC mass bound defined as 0 for d > 0; `dualTowerMass` has identical branches; c(N) = 100(N+1) and
  bounce action 42(n+1) chosen by hand; `is_minimum_action := true`; `rrCharge = rank + c1 + c2`;
  `classTadpoleAnomaly` returns 0 for every class; Callan–Harvey terms defined as −x and +x;
  von Neumann entropy defined as 0.0. These are NOT verification of the physics they name.
- Genuinely sound but modest: Mukai pairing symmetry/additivity and FM-isometry dimension preservation
  (`proofs/MukaiLatticeK3.lean`); Kummer b₂ = 22, χ = 24, signature (3,19) as arithmetic on entered
  constants (tier B); EOT multiplicity arithmetic and 462·60 = 360·77 (tier B).
- `proofs/FTheoryCosmology.lean` is not a build target and does not compile.
- Orientifold bookkeeping: files cite Gimon–Polchinski (hep-th/9601038) for the 16 fixed points of
  T⁴/ℤ₂ while assigning O7⁻ planes (charge −4) and 32 D7-branes (+2). The construction cited and the
  O-plane type assigned have not been reconciled against the source. Status: UNRESOLVED (tier C,
  under revision). Do not state the correct O-plane type.

## F2. Replacement: kernel-checked foundation (LeanMaster)
- LeanMaster clone at tag v2.2.0 (commit b27ce8b) built (3708 jobs). Its own gates re-run here:
  `DualScaleStream2` 99 theorems, 0 failing; `StringTheoryFormalization` 89 theorems, 0 failing;
  statement lock OK. Its certificate `docs/VERIFIED_FOUNDATION.md` is pinned to v2.1.0; per its
  author the tree differs from v2.1.0 only in comments.
- New downstream project `lean_foundation/` (library `DualScaleFoundation`, Lean v4.33.1, Mathlib v4.33.1):
  build 3665 jobs, 0 `sorry`, 20 theorems audited, 0 failing; negative control detected. Statement
  review (lock) pending with the project lead. Theorems and what they rest on:
  - T-duality: `theta_shift_is_odd` (`thetaShift_isODD`), `basis_change_preserves_odd`
    (`basisChange_isODD`), `odd_mul_property` (`isODD_mul`), `tduality_metric_inversion`
    (`tduality_inverts_metric`: ηR·H(G,0)·ηR = H(G⁻¹,0) for invertible G),
    `factorized_tduality_is_odd` (`factorized_isODD`), `factorized_tduality_involution`
    (`factorized_mul_self`: factorized k · factorized k = 1). Scope: integer O(d,d;ℤ) charge action and
    the generalized metric at B = 0 — NOT the field-level Buscher rules for metric, B-field and dilaton.
  - Lattices: `k3_signature_is_three_nineteen` (`sigK3_eq`), `mukai_signature_is_four_twenty`
    (`sigMukai_eq`), `k3t2_signature` (`sigK3T2_eq`, (6,22)), `k3_rank` (22).
  - Flux: `k3_cross_k3_euler_over_24` (χ(K3×K3)/24 = 24), `k3t2_euler_is_zero`, `tadpole_conservation`
    (flux + n = 24 ⇒ n ≤ 24 ∧ (flux = 0 → n = 24)).
  - Moonshine: `eot_values` (A_n = 45, 231, 770, 2277, 5796, 13915, 30843, 65550, 132825),
    `first_five_m24_irreps`, `a6_not_m24_irrep`, `a6_decomposition`,
    `r_bps_from_eot` (2·A₂ / (4·2·A₁) = 462/360 = 77/60 as a rational identity).
    77/60 is tier A as arithmetic; reading it as a primordial bispectrum / f_NL ratio is tier C with no
    derivation in this project.
- One-direction scope limits (from LeanMaster): `thetaShift_isODD`, `basisChange_comm_thetaShift` are
  one direction; generation of O(d,d;ℤ) and classification of even unimodular lattices are tier L.
- Citation form (use verbatim, filling <Module>.<name>): "Kernel-checked in Lean 4 (v4.33.1, Mathlib
  v4.33.1) in DualScaleStream2.<Module>.<name>, repository xaviercallens/SocrateAI-Scientific-Agora-LeanMaster,
  release v2.2.0; depends only on Lean's three standard axioms." Downstream restatements live in
  `lean_foundation/` of this repository.
- Forbidden: "zero axioms", "100% verified", "proves string theory", "formalizes string theory",
  "certified" for anything below tier A.

## F3. Simulations (Rust `rust_simulator/`, Python)
- `vacuum_decay_cdl.rs` is a flat-space Coleman bounce (no gravity), not Coleman–De Luccia. Its
  shooting never detected overshoot (clamp at φ_false), φ_true was hard-coded 2.05 while the potential's
  minimum with the default parameters is ≈ 1.925, and the action was floored at 1.0. The central charge
  c(N) = 100(N+1) is a model assumption, not a derived c-theorem. Fixes are on branch; numbers will be
  regenerated.
- `swampland_geodesic.rs`: explicit Taylor step (not symplectic/Velocity-Verlet); bound check had 1.5×
  slack; α = 1/√2 hard-coded, not measured. Fix (RK4, strict bound, α_effective) on branch.
- Committed summaries disagree: bounce action 137.15 / bubble radius 2.82 (`vacuum_decay_cdl_summary.json`)
  vs 316.0 / 1.65 (`frontier_loops_summary.json`).
- TDA Mapper: committed `tda_mapper_skeleton.json` = 187 nodes, 557 edges, β₁(graph cycles) = 376;
  `simulation_results.json` = 186 / 528 / 348; re-running pytest regenerates different values (Mapper not
  seeded). β₁ here is the cycle rank of the Mapper graph, not persistent homology. Its
  `tadpole_cancellation_certified` field is a stored flag, not a kernel result.
- Rust simulator had 0 tests before the fix branches.
- rusty-SUNDIALS: loaded only from an absent sibling repository; no SUNDIALS run exists in this repo.
- All simulation outputs are tier X.

## F4. Benchmarks and costs
- No benchmark harness, timing log, or execution record exists in the repository for: 1,520× vs SciPy,
  >900× vs SciPy, 7.1× vs C++ SUNDIALS/CVODE, 12–18 ns per step, 45 ms IPC gate, < 3 µs per sample,
  $0.0163 / $0.0293 / $1.542 costs, 99.1% cost reduction. Required wording: "pending hardware
  verification"; remove the numbers from abstracts and conclusions.

## F5. Publication record
- Zenodo record 22683565 (DOI 10.5281/zenodo.22683565, concept 22683564) is published with the
  original claims. The corrected manuscript will be issued as a new version of that concept record.
  Papers should include a short "Correction notice" stating what was withdrawn and why (F1, F3, F4).

## F6. Numbers not yet final
Simulation numbers will be regenerated after the numerics fixes merge. Wherever a simulation number is
kept, wrap it as `\numreconcile{<value>}{<quantity>}` (macro defined in the preamble) so the post-merge
pass can find and replace it. Do not invent new values.
