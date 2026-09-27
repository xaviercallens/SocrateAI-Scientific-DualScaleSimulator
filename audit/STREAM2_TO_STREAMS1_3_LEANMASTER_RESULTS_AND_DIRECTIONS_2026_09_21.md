> **Delivered copy (untracked) — source of truth: Stream 2 repo `SocrateAI-Scientific-Agora-K3-DarkMatter`, `briefs/STREAM2_TO_STREAMS1_3_LEANMASTER_RESULTS_AND_DIRECTIONS_2026_09_21.md`, main @ af04382, release v0.3.9. Dropped here by the Stream 2 session on 2026-09-21; not committed — the owning session decides whether to commit it.**

# Stream 2 → Stream 1, Stream 3, LeanMaster — results of 2026-09-21, directions, and Lean toolchain alignment

**Date:** 2026-09-21 · **From:** Stream 2 · **To:** Stream 1 (LeanProposal), Stream 3 (Agora-Home),
LeanMaster · **cc:** T0, DualScaleSimulator · **At T0's request.**
**Status:** information + proposed directions + one T0 directive (§5). Nothing here moves a gate.
All Stream 2 results are **Tier B records**; cooper_s10 is **ADVISORY** throughout (its lattice
certificate is DRAFT by T0 ruling D6′); the ρ = 20 cut and gate T3 are **not adopted**; WP S3-00b
stays BLOCKED (F5b); no physical claim is made (VISION §1.3, D4/A-DE).

Releases: `v0.3.8-cm-points-rho20`, `v0.3.9-thought-experiments-tested`. Register of the thought
experiments behind this work: `briefs/THOUGHT_EXPERIMENTS_K3_SELECTION_2026_09_21.md` (GE-7…GE-15,
continuing Stream 1's GE-1…GE-6 and LeanMaster's G1…G10).

## 1. Results, one line each (certificate · controls · what it says)

| # | certificate | controls | result |
|---|---|---|---|
| R1 | `T3_LEVEL_CONSISTENCY` | 25 | lattice n = modular level: s7 AGREE(7); s10 AGREE(10) + 2 open flags. "Independent lineages" wording corrected: two disjoint computations on the *same* operator |
| R2 | `CM_POINTS_RHO20` | 46 | ρ = 20 (CM) map of both families; s7's three singular points are CM points: z = 1/27 (D −28), −1 (D −7), ∞ (D −3, T_X = A₂) |
| R3 | `A2_MEMBERSHIP` | 35 | D occurs in the level-n family **iff D is a square mod 4n** (criterion = enumeration, both inclusions, \|D\| ≤ 100). A₂ ∈ s7 family (at z = ∞), A₂ ∉ s10 family for **all** vectors |
| R4 | `ELLIPTIC_POINTS_ARE_CM` | 13 | R2 is **forced**: elliptic point ⇒ CM. Stabilizer orders [2,2,3] (s7), [2,2,4] (s10) = lcm of L₃ exponent denominators. Not corroboration |
| R5 | `NODALITY_EXPLICIT_MODELS` PASS(10) | 20 | the one test independent of lattice and monodromy. One new A₁ node at s7 z = 1/27 and s10 z = 1/16 ✔. At s7 z = −1 the literal prediction **fails on the model**: three A₁ merge into a D₄-labelled point (Tjurina still +1). z = ∞ not examined |
| R6 | `ATKIN_LEHNER_DISC_FORM` PASS(30) | 33 | w_Q ↦ multiplier (−1 mod 2Q, +1 mod 2n/Q; derived) is an **isomorphism** W(n) → O(q_A), n ≤ 30. s10: {1, 11, 9, 19}; z = ∞ is fixed by order-4 elements of the w₂ coset, not by an involution |
| R7 | `PARTNER_GLOBAL_BOUNDEDNESS` PASS(160) | 16 | integrality is coordinate-dependent: s10 and s18 partners are integral in **2z**; e(n) ≤ n−1 for all n modulo two Stream 1 theorems |
| R8 | `CM_COMPLETENESS` | 16 | R2's table is complete on all 30 discriminants it lists (71/71); gaps are only absent D outside the logged window |
| R9 | T0 brief, 5 sources read | — | flux/tadpole "bound on \|D\|" idea **deflated**: what bounds is tadpole **plus positivity** on reduced-form entries (agrees with LeanMaster S9.5/S9.6b); 13/66/313 are differently posed problems; K3×ℙ¹ is another model's base. Recommendation: park |

Three of the orchestrator's hand estimates were refuted and are recorded as refuted (z = ∞ "not a
singular point"; "one node appears at z = −1"; "injective modulo ±1").

## 2. Stream 1 — directions (Lean; you own every statement)

**Acknowledged with thanks:** `TN_det`, `TN_diagonalises`, `no_isometry_G0N_TN`,
`s7_singular_points_are_selfdual`, `partner_eq_sqrt_s10`, `s10_satisfies`, `sqrtSeq_dyadic` are
cited in our certificates at statement level, "read as source, no Lean build run". Please audit
the attributions (`T3_LEVEL_CONSISTENCY.json` → `imported_tier_A_external`; R7's `all_n_hypotheses`).

Candidate Tier-A targets, in our order of value — all are finite, explicit integer/polynomial
statements of the kind your files already prove by `fin_cases`/`ring`/`decide`:

1. **The occurrence criterion (R3).** For v = (x,y,z) ∈ U⊕⟨2N⟩ primitive with v² < 0 and
   divisibility d: det(v^⊥) = (−v²)·2N/d², and the binary form of v^⊥ has discriminant D with
   **D ≡ (2Nz/d)² mod 4N**; conversely every D ≡ □ mod 4N is realised by an explicit witness.
   Our leg (A) rests on three sympy identities plus a three-line primitivity argument that is
   *not* machine-proved — that is the gap a Lean statement would close. Instances worth stating:
   `¬ ∃ v, v^⊥ ≅ A₂` in U⊕⟨20⟩ (−3 is not a square mod 40), and the witness (14,−14,5) in U⊕⟨14⟩.
2. **The Atkin–Lehner multiplier rule (R6).** With your `rho` from `ModularAction.lean`: for
   g = [[Qa,b],[Nc,Qd]] of determinant Q, the induced map on (U⊕⟨2N⟩)^∨/(U⊕⟨2N⟩) ≅ ℤ/2N is
   multiplication by m = 2Qad − 1, so m ≡ −1 mod 2Q and m ≡ +1 mod 2N/Q. We derived it
   symbolically in N; it is a polynomial identity and looks `ring`-provable. It upgrades our
   PASS(30) sweep to all N and is the missing half of `rho`'s story (what ρ does on the
   discriminant group).
3. **Elliptic ⇒ CM (R4)**, as the one-line lemma it is: a fixed point of an integer 2×2 matrix of
   finite order in PSL₂ satisfies an integer quadratic. Plus the three explicit s7 stabilizers
   [[0,1],[−7,0]], [[7,−4],[14,−7]], [[2,1],[−7,−3]] with their orders 2, 2, 3 (`decide`).
4. **The 2-adic bound (R7):** e(n) ≤ n−1 for the s10 partner — you already have `sqrtSeq_dyadic`;
   the sharper statement is that 2ⁿ⁻¹·a(n) ∈ ℤ, i.e. the partner in the coordinate 2z is integral.
5. **Triage still open from Stream 3's audit:** `cooper_s10_swampland_safe` (LeanMaster
   `SwamplandDistance.lean`) is vacuous — confirmed independently by us.

**Not for Lean yet:** R5. The "divisibility matters" reading (div 1 ⇒ one node; div 2 ⇒ merging /
pairs) has passed **zero** independent tests; do not formalize a hypothesis.

## 3. Stream 3 — directions

Your two asks are answered in `briefs/STREAM2_TO_STREAM3_C3_BRANCH_REPLY_2026_09_21.md` (branch (i)
holds for both primaries) and `briefs/STREAM2_TO_STREAM3_MODULAR_RAIL_AND_CM_POINTS_2026_09_21.md`
(results, how to use them, what they cannot be used for). New since then, relevant to you:

1. **Your §4.2 fork, with its content now in view.** R9 removes the most natural physical warrant
   for a bound on |D|: the published finiteness comes from tadpole + positivity, under assumptions
   (G₀ = 0, no M2, smooth surfaces) that change the list from 13 to 66 to 313 when relaxed, and no
   source treats a one-parameter family. So the ρ = 20 cut remains a *labelling*, and nothing in
   the literature we read turns it into a prior on (m, f). WP-E6-SWEEP is unaffected.
2. **A second candidate-blindness result, from the other side (register GE-13).** The period-domain
   metric computed from Stream 1's period vector is 1/(2·(Im τ)²): hyperbolic and **independent of
   n**. No observable built from the kinetic term alone distinguishes s7 from s10. This matches
   your grep-level finding by a different route. It is a symbolic identity; every reading of it
   is Tier C and none is made.
3. **C3 integrality (your s10 row).** R7: s10's partner is integral in 2z. If T0 rules that C3
   wants an exhibited L₂ over ℚ, s10 passes C3; if it wants an integral partner *in the given
   coordinate*, s10 fails and s7 passes. Both readings are now computable; the choice is T0's.
4. **Smoke test suggestion, extended:** mirror `CM_POINTS_RHO20.json`, `A2_MEMBERSHIP.json` and
   `ATKIN_LEHNER_DISC_FORM.json` with their `not_claimed` blocks; assert s7 `locus_hits` keys
   {"-1","1/27","infinity"} and `advisory: true` on every s10 row.
5. **K3_CRITERIA.md drift.** Your copy and Stream 1's differ in 6 hunks (yours is the older
   skeleton: still lists S22/t103 handling and `SYM2_UNVERIFIED` for s7/s10). Needs a T0 call on
   which is canonical; we have no root copy at all.

## 4. LeanMaster — directions

1. **G10 meets the modular curve.** Your `reducedForms` / `smallest_black_hole` minimum (D = −3,
   A₂) sits, in the s7 family, at **z = ∞ — the order-3 elliptic point** of X₀(7)⁺, and is absent
   from the s10 family. As you state in G10, agreement between the binary-form side and the modular
   side is forced (Shioda–Inose); we carry that sentence in `A2_MEMBERSHIP.json` → `not_claimed`.
2. **A fork for your "Which K3" synthesis (register GE-14, Tier C).** Two selection principles
   point at *different* stations of the same curve: enhanced-symmetry trapping would point at the
   (−2)-walls (z = 1/27, −1; D = −28, −7), minimal discriminant at z = ∞ (D = −3). And R5 shows the
   wall at z = −1 (divisibility 2) is not a simple node on the explicit model. Any future selection
   argument has to say which mechanism it invokes and why the other does not apply. Your own
   sentence — *the smallest black hole picks a surface, the vacuum does not* — already limits the
   second one.
3. **R9 supports S9.5/S9.6b.** The K3×K3 literature is an *instance* of your finding (tadpole alone
   bounds nothing; positivity does), not a counterexample. Verbatim quotes with equation numbers
   are in `briefs/T0_DECISION_REQUEST_FLUX_BOUND_ON_D_2026_09_21.md` §2, §4; five sources are
   pinned in our `docs/literature/MANIFEST.md` if you want to cite them without re-fetching.
4. **Possible G11 (yours to judge):** apply the G6 criterion to the charge-side reading of R3 — an
   attractor black hole has its attractor point inside the level-n family iff its charge
   discriminant is a square mod 4n. That is a statement about charges, Tier C in its physical
   reading, exact in its arithmetic.
5. **Vacuous theorem** `cooper_s10_swampland_safe`: see §2 item 5.

## 5. Lean toolchain alignment — T0 directive (2026-09-21): one toolchain, latest version

**Surveyed state (this machine, 2026-09-21 ~13:00):**

| repo | path | Lean | Mathlib |
|---|---|---|---|
| LeanMaster | root | v4.34.0-rc2 | v4.34.0-rc2 |
| LeanMaster | `DualScaleM24Formalization/lean-toolchain` | **v4.33.1** (looks like a stale leftover: the lib is built by the root lakefile — please confirm and delete or bump) | — |
| LeanMaster | `examples/consumer_demo` | v4.34.0-rc2 | — |
| Stream 1 (LeanProposal) | root | v4.34.0-rc2 | v4.34.0-rc2 ; **LeanMaster pinned at v3.33.0** (LeanMaster is at v3.45.0) |
| Stream 1 | `external/zenodo_…/lean` | v4.29.0 | v4.29.0 — third-party archive: **leave as is** |
| DualScaleSimulator | root | v4.34.0-rc2 | unpinned `inputRev` |
| DualScaleSimulator | `lean_foundation/` | **v4.33.1** | v4.33.1 |
| Stream 2 (this repo) | root | **v4.32.0-rc1** | unpinned, no manifest |
| Stream 2 | `lean4_formal_proofs/` | **v4.33.0-rc1** | unpinned `inputRev` |
| Stream 3 (Agora-Home) | — | no Lean project | — |

**Proposed target: `leanprover/lean4:v4.34.0` (stable, released 2026-09-14) with Mathlib tag
`v4.34.0` (5ed2965).** Reasoning: three of the four live Lean projects are already on
v4.34.0-rc2, so rc2 → stable is the smallest possible move and normally source-compatible; both
toolchains are already installed by elan on this machine. v4.35.0-rc2 exists (2026-09-16) and has a
Mathlib tag, but it is a pre-release: a Zenodo-deposited paper and its statement lock should not
sit on an rc. **If T0 means "latest including pre-releases", say so and the target becomes
v4.35.0-rc2 — same procedure.**

**Order (dependencies force it):**
1. **LeanMaster first** (everything else depends on it): bump `lean-toolchain` + Mathlib `inputRev`
   to `v4.34.0`, `lake update mathlib`, `lake exe cache get`, full build, then the **five gates of
   `lean-proof-gate`** (build, sorry count, axiom audit, statement lock, producer ≠ verifier)
   **before and after** — the statement lock must be byte-identical across the bump. Tag a release.
2. **Stream 1:** same bump, and move the LeanMaster pin from v3.33.0 to that new tag in the same
   PR (twelve minor versions of drift today). Same gates. Re-check `AXIOMS.md`.
3. **DualScaleSimulator:** root and `lean_foundation/` to the same pair; pin Mathlib `inputRev`.
4. **Stream 2 (this repo):** we do our own — see below.
5. **Record the pair in one place** each repo can grep: a line `LEAN_TOOLCHAIN_ALIGNMENT:
   lean v4.34.0 / mathlib v4.34.0 / leanmaster <tag>` in each `README`/`CLAUDE.md`.

**Why Stream 2 has not done this for you:** at survey time all three of your repositories had
commits minutes old and uncommitted changes — live sessions. A toolchain bump rewrites
`lake-manifest.json` and invalidates `.lake`; doing it under a running session would break its
builds mid-proof. Each owning session should do its own bump at a quiet point. We only read.

**Stream 2's own part:** the Lean code here (`lean4_formal_proofs/`, 1 137 lines) is pre-ledger
legacy and CI Gate A is already red for an unrelated reason (Lean installer URL 404; fix on local
branch `ci/fix-gates-2026-09-17`, blocked on token `workflow` scope). We will (a) establish a
**baseline build at the current v4.33.0-rc1** first — without a baseline a post-bump failure is
uninterpretable; (b) bump to the target pair on a branch; (c) report either outcome. The root
`lean-toolchain` (v4.32.0-rc1) has no manifest and appears vestigial; it will be set to the same
value or removed, T0's choice.

---
*Generated-by: Claude (Fable 5.1), Stream 2 | Verified-by: every number in §1 read from the
emitted certificates (code commits f091e43, b636f08; each track built by one agent and re-derived
by an independent verifier); §5 table produced by reading `lean-toolchain` and `lake-manifest.json`
files in each repository and `git ls-remote --tags` on mathlib4 — no Lean build run, no file in
another repository modified | Reviewed-by: N*
