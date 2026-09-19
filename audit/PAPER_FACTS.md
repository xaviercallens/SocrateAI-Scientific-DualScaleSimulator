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
- Orientifold bookkeeping: RESOLVED 2026-09-18 at literature level (tier L). Sources pinned on branch
  `loop/k3t2-rigidity`, `audit/k3t2_rigidity_v2/sources/` (sha256 in SHA256SUMS).
  - Type IIB on K3 × T²/ℤ₂ (Tripathy–Trivedi, hep-th/0301139, §2.2, text lines 160–176):
    - "The ℤ₂ orientifold symmetry has 4 fixed points on the T², an O7-plane is located at each of
      these fixed points. To cancel the resulting 7-brane charge 16 D7-branes need to be added."
    - O7s and D7s wrap K3. That induces 2 units of D3 charge per O7 and 1 per D7: 4·2 + 16·1 = 24.
    - Their eq. (2.3) is **½ N_flux + N_D3 = 24**, where N_flux = (2π)⁻⁴(α′)⁻² ∫ H₃∧F₃ is taken over the
      covering space (footnote: in F-theory, 24 (p,q) 7-branes each carry one unit).
    - In their symmetric configuration each O7 carries 4 D7s (lines 312–314).
    - Hence −4 per O7 in D7 units. That value is an inference from 16/4, not a quote.
  - The 16 fixed points of T⁴/ℤ₂ are the orbifold limit of the **K3 factor**. The IIB orientifold does
    not act on K3, so **no O7-plane sits at a Kummer point**.
  - In the Type I / Gimon–Polchinski description (T-dual along T²) those points carry O5-planes and
    half five-branes. BLPSSW (hep-th/9605184, §1, text lines 227–229) say "there must be one instanton
    hidden at each fixed point", and after blow-up "eight 5-branes on a smooth K3" remain (p. 10).
    16 + 8 = 24 is our inference.
  - Therefore the assignment "O7⁻ (−4) at each of the 16 Kummer points + 32 D7 (+2)" in
    `proofs/KummerTDAAnomalyCertification.lean`, `proofs/TadpoleCancellation.lean` and the old paper
    text describes no construction. Its arithmetic 64 − 64 = 0 is true but has no physical content.
    **Withdrawn.**
  - LeanMaster `DualScaleStream2.Flux.Tadpole.tadpole_budget` (flux + n = 24, restated here as
    `lean_foundation` `tadpole_conservation`) has the same form as eq. (2.3) with flux := ½ N_flux.
    - LeanMaster sources it to Dasgupta–Rajesh–Sethi (M-theory on K3 × K3, ½∫G∧G + n = 24 = χ(K3×K3)/24).
    - That setting is dual to IIB on K3 × T²/ℤ₂, so two independent literature routes give the same 24.
  - LeanMaster had the same 16-O7 assignment and a D3 target "= χ(K3)/24 = 1". This session
    reported both (`audit/LEANMASTER_NOTE_orientifold.md`). **Fixed in LeanMaster v3.21.0**, checked
    here read-only:
    - `TadpoleCancellation.lean` now has 4 O7 at −4, 16 D7 at +1, and 4·2 + 16·1 = χ(K3×K3)/24 = 24.
      `d3_tadpole_target_is_one` and `total_O7_charge_is_minus_64` are removed.
    - `KummerTadpole.lean` docstrings: units are stated and the citation now points to TT and Sen.
      LeanMaster found that GP (3.12) is a Chan–Paton projection, not a charge normalisation.
    - The TT eq. (2.3) citation is in `StringTheoryFormalization/StringDynamics/TadpoleConstraint.lean`.
    - LeanMaster's review record: `docs/reviews/2026-09-18_dualscalesimulator_orientifold_note.md`.

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

## F7. Parameter classification (constant source accountability)
The manuscript should not claim "no free parameters" without specifying the source tier for each constant.
Below is the full classification (Mathesis tier A/B/L/C/X methodology):

**Tier A (Lean kernel-derived, adequate statement):**
- `c112 = 77/60`: Kernel-checked arithmetic (462·60 = 360·77, gcd(77,60)=1) in MathieuVertexOperators.lean
- `h1 = 0.25, h2 = 1.25`: Candidate Lean-derivable via Kummer K3 geometry (depends on completion of derived theorem chain; currently assumed in code)

**Tier B (exact arithmetic with negative control):**
- Kummer surface invariants: b₂=22, χ=24, signature (3,19) via Mukai pairing (proven in MukaiLatticeK3.lean)
- EOT multiplicities A_n for n=1..9 via explicit enumeration in MathieuVertexOperators.lean

**Tier L (literature fixed):**
- `RHO_M0 = 0.315`: Planck 2018 ΩM (https://doi.org/10.1051/0004-6361/201833910)
- `RHO_R0 = 9.2e-5`: Planck 2018 ΩR
- `DESI_w0_mu = -0.827, DESI_wa_mu = -0.75`: DESI 2024 DR1 BAO+CMB+SNe fit (hard-coded at workshopcosmo.py:438–439; cite release DOI)

**Tier C (conjecture/physical assumption, no derivation):**
- `a_pot = 1.0, b_pot = 0.01`: Potential well depth (ansatz for screening model)
- `mu_sym = 1.0, lambda_sym = 1.0`: Symmetron coupling (model scaling; no derivation from fundamental theory)
- `pta_suppression = 0.005`: PTA l=4 hexadecapole suppression (workshopcosmo.py:700; comment cites 10^-3, code uses 5×10^-3; needs justification)
- `NANOGrav_c4_c0_ratio = 16.07`: Injected into simulation (workshopcosmo.py:691); described as "injected" in docstring

**Tier X (exploratory/numerical, regenerable, not claimed as fixed):**
- All Langevin SDE step-size and noise coefficients (numeric integration tuning, not physical constants)
- TDA Mapper clustering resolution and filtration thresholds
- Simulation output (κummer_langevin_*, tda_mapper_*, vacuum_decay_*, swampland_geodesic_* JSON/CSV files)

**How to apply:** Before any manuscript paragraph claims "no free parameters" or "all parameters determined by theory," enumerate the constants involved, list their tier, and confirm that tier A+B+L cover all required terms. Any tier C or X constant means the statement should be qualified: "N of M parameters are theory-fixed; the remaining M-N are determined by [external data / physical assumption / model tuning]."
