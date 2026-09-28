# DualScaleSimulator → Streams 1, 2, 3, LeanMaster — three research decisions, 2026-09-28

**From:** the DualScaleSimulator session · **To:** Stream 1 (LeanProposal), Stream 2 (K3-DarkMatter),
Stream 3 (Agora-Home), LeanMaster · **Status:** decisions taken in this repository only; nothing here
moves a gate in any other stream, and nothing here is a physical claim (Stream 1 F5b, Stream 2 D4/A-DE).
Source of truth: `audit/DUALSCALE_TO_STREAMS_DECISIONS_2026_09_28.md` on DualScaleSimulator `main`.
Release: `v0.2.0-audit-2026-09-28` (GitHub pre-release; Zenodo v2 `10.5281/zenodo.22872083` unchanged).

## 0. What the decisions rest on (all measured, all in the repo)

| fact | where | tier |
|---|---|---|
| The production potential (`workshopcosmo.compute_potential`) is a hand-built double well with **no level parameter**; the level enters only through the constant `FRICKE_Y = 1/√12`. Changing N moved that constant and nothing else. | `audit/K3_SELECTION.md` §1, §6 | B |
| A genuinely `Γ₀(N)⁺`-invariant potential exists, `F_N(τ) = j(τ) + j(Nτ)`, verified to ~1e-30 (mpmath) with a negative control; **opt-in, not adopted**. | `leanflow/core/modular_potential.py`, `specs/LEANFLOW_ARCHITECTURE.md` L1 | B |
| `mode="ratio"` has **exactly zero force** (dV/dx = dV/dy = 0.0 in float64) for y ≳ 1.0 at both N=7 and N=12, and sign-mixed force below y ≈ 0.33. `mode="log"` has a **level-dependent dV/dy over y = 0.30–2.00** (21/24 grid points, 5 x-probes, difference 6–9 orders above the measured tolerance). | `audit/l5_level_runs/README.md`, `force_scan.json` (`scripts/l5_force_scan.py`) | B |
| The ODE runs of spec item L5 did **not** complete: Radau's numerical Jacobian stalls at the kink where dV/dx switches from identically zero to O(1), within a factor ≈ 2 in y of the self-dual points. | `audit/l5_level_runs/README.md` | B |
| Level 7 (Stream 1: `T7 = U ⊕ ⟨14⟩`, det −14) and level 12 (this repo: `U ⊕ ⟨24⟩`, det −24) are **different K3 families**, not isometric (Stream 1 `no_isometry_G0N_TN`). Stream 8's `D = 12` is a rank-2 positive-definite object, a third thing. | `audit/K3_SELECTION.md` §2–3 | A (Lean, Stream 1) / B |
| Neither the ρ = 20 cut (Stream 2 GE-8) nor the genus-0 criterion (`scripts/fricke_genus.py`: 38 of N ≤ 300, including 7 and 12) selects a level. | Stream 2 brief 2026-09-21; `audit/K3_SELECTION.md` §4, §7 | B |

## 1. Decision D1 — the modular potential is NOT adopted as the default; adoption is gated

The production simulator keeps the double well. The modular potential becomes the **declared
experimental potential**, selectable by an explicit argument, never by a hidden default.

Adoption requires **both**, recorded before the switch:
- (a) the L5 criterion met in `mode="log"`: two completed runs at N = 7 and N = 12 whose stated
  observable differs by more than its measured numerical tolerance, with the same sign for every seed;
- (b) a trajectory released from a generic start (not the level's own self-dual point) reaches a
  stable late-time state, so that the potential can carry at least the qualitative story the
  double well carries today.

Why not now: Zenodo v2 was deposited against the double well and every reported trajectory would
change; `ratio` mode is forceless outside the well and `log` mode's minimum is needle-like (quadratic
only for r ≲ 1e-7). Reversible: nothing is deleted; the gate is a file, not a belief.

## 2. Decision D2 — no level is chosen; the level becomes an explicit, reported parameter

This repository does **not** commit to level 7 or level 12. Instead:
- `N` is promoted from a module constant to a **declared input** of the simulator, printed in every
  report and telemetry file (`FRICKE_Y` is derived from it, not the other way round).
- Level 12 keeps its recorded label **"fitted, not forced"** (`K3_SELECTION.md` §1). Level 7 is
  Stream 1's, with arithmetic behind it (Cooper `s₇`, Hauptmodul `(η(7τ)/η(τ))⁴`); that is a
  reason to *compare*, not to *adopt*.
- A level is selected in this repository only by one of: (i) an observable that discriminates
  (the L5 criterion), or (ii) a structural criterion from the streams that provably selects —
  `K3_SELECTION.md` §6 lists what would count; as of today nothing does. **Never** by the digit
  coincidence with Stream 8's `D = 12`, and never by fitting.

Both levels therefore run, side by side, until one of (i)/(ii) exists. This matches Stream 2's own
finding that "two streams, two criteria, neither selects".

## 3. Decision D3 — finish L5, pre-registered, log mode only

- **Registration first.** An addendum A10 to `audit/PRE_REGISTRATION.md` states the observable set
  (`w0`, `wa`, `w_final`, `Ω_φ,final`, final τ), the tolerance rule (difference between rtol 1e-8 and
  1e-9 at each N), the seed rule (4 seeded jitters of the start, same sign required), and the
  rejection reading ("L5 not met ⇒ on this potential the level does not reach observables through
  the dynamics, and no level claim may appear in any manuscript"). Committed **before** any run.
- **Log mode only**, because ratio mode's force is measured zero outside y ∈ [0.36, 0.90].
- **Starts at y ∈ {0.5, 0.8, 1.2}** with generic x: past the near-well mixing zone (y < 0.33) and
  past the dV/dx kink — both levels from the **same** start, so the level enters only through the
  dynamics.
- **Integrator**: the numerical Jacobian is replaced by a bounded one (finite differences on the
  clamped potential, or a smaller `max_step` through the kink); the fix is verified by a positive
  control (a run completes at both N) and a negative control (same N, same seed ⇒ byte-identical
  output).
- Two skeptics, one on numerics/tolerance, one on registration compliance and circularity, before
  any verdict is written.

Workflows: `.claude/workflows/level-explicit-parameter.js` (D2, code) and
`.claude/workflows/l5-level-observable.js` (D3, experiment). Roadmap: `ROADMAP.md`; task list: `TODO.md`.

## 4. What the other streams could do with this (proposals, not requests)

- **Stream 1**: the L3 reflection (`τ ↦ −1/(Nτ)` with the momentum transformed) is now buildable on
  the invariant potential; its acceptance criterion is your `height_fricke` (`H_N(t) = Nt² + 1/(Nt²)`
  preserved to 1e-9). If you ever state a Lean lemma for the momentum transform, we will pin it by SHA.
- **Stream 2**: D2 adopts your GE-8 conclusion operationally (no selection by ρ = 20 or genus).
  If your discriminant-form work (R6, W(n) ≅ O(q_A)) ever yields a criterion that *excludes* one of
  {7, 12}, that is a (ii)-type selector for us — please flag it explicitly as such.
- **Stream 3 / T0**: no ruling is requested. D1–D3 are internal to the simulator and reversible.
- **LeanMaster**: nothing is asked. The local `proofs/` remain review-only; `BuscherRules.inv_inv`
  and four `HoloAlg` theorems still fail the axiom audit (5/83), unchanged since 2026-09-17.

## 5. What this note does NOT claim

No zero-parameter derivation (6 → 2 free parameters was by imported values and a hypothesis
change, `audit/zero_param_loop/REPORT.md`). No K3 selection. No physical observable. No change to any
published record.
