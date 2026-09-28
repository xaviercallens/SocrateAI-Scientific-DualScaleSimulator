# DualScaleSimulator — roadmap

Last updated 2026-09-28. Decisions and their evidence: `audit/DUALSCALE_TO_STREAMS_DECISIONS_2026_09_28.md`.
Task list: `TODO.md`. Session re-orientation: `scripts/session_status.sh`.
(`specs/roadmap.md` is the September strategy note that started the AXE 1–3 work; it is history, not this roadmap.)

Tiers follow Mathesis: **A** kernel-verified and adequate · **B** exact/numeric with negative control ·
**L** literature · **C** conjecture · **X** exploratory. Nothing in this repository is above B except
what it imports from LeanMaster / Stream 1 and cites by SHA.

## Where the programme stands

| track | state | evidence |
|---|---|---|
| Six-axis simulation pipeline (quintessence, symmetron, tadpole Lean, moonshine, NANOGrav, Kummer-Langevin + TDA) | runs end-to-end, telemetry byte-reproducible, CI green | `TELEMETRY_PROVENANCE.json`, CI run 36473488754 |
| Local `proofs/` | review-only; 78/83 theorems standard-axiom, 5 fail (`BuscherRules.inv_inv`, 4 × `HoloAlg`); replacement is external LeanMaster | `audit/lean_axiom_report.json`, `audit/lean_replacement_map.md` |
| Free parameters | 6 → 2 (`mu_sym`, `c4_pta_product`), by imported Ω_Λ and a ψ = φ/φ₀ rewrite — **not** a derivation; 2 → 0 closed by derivation, tested only as the hypothesis change M0 | `audit/zero_param_loop/REPORT.md`, `audit/reverse_zero/` |
| M0 vs data | 0.86σ from ΛCDM without CMB priors, **2.16σ disfavoured with** them; DR2 CPL 1.63σ is a pipeline action, not a P1 verdict | `audit/reverse_zero_r2/` |
| K3 × T² blind re-derivation | v3: 68 agree / 0 disagree / 261 not computed of 332 LeanMaster targets; rigidity 3 of 35 | `audit/k3t2_rigidity_v3/REPORT.md` |
| TDA instrument | alpha path validated on real Re6Zr vortex maps; lower-star path needs a value-distribution-preserving null; three pipeline defects repaired with a regression guard; CMB lenses null | `audit/tda_validation/`, `audit/FLUID_TO_COSMOLOGY_BRIDGE.md` |
| Which K3 | level 12 is **fitted, not forced**; level 7 is Stream 1's; no criterion selects | `audit/K3_SELECTION.md` |
| Modular potential (L1) | built, invariant to ~1e-30, opt-in; log mode shows a level-dependent force over y = 0.30–2.00; L5 ODE runs incomplete | `specs/LEANFLOW_ARCHITECTURE.md`, `audit/l5_level_runs/` |
| Papers / deposit | Zenodo v2 `10.5281/zenodo.22872083` (2026-09-21); GitHub pre-release `v0.2.0-audit-2026-09-28` | `project_zenodo_deposit` memory, `audit/PAPER_FACTS.md` |

## The three decisions of 2026-09-28

- **D1 — modular potential: gated, not adopted.** Adoption needs the L5 criterion met in log mode
  *and* a generic-start trajectory reaching a stable late-time state.
- **D2 — no level chosen; N becomes an explicit reported parameter.** Both levels run; selection only
  by an observable that discriminates or a structural criterion that provably selects.
- **D3 — finish L5, pre-registered (A10), log mode, starts at y ∈ {0.5, 0.8, 1.2}, bounded Jacobian,
  two skeptics.**

## Milestones (in order; each has a checkable exit)

1. **M1 — level explicit** (`.claude/workflows/level-explicit-parameter.js`). Exit: `workshopcosmo`
   takes `level` and `potential` as arguments; `FRICKE_Y` is derived; a test shows the modular
   potential's output changes with `level` and the double well's does not (documenting the defect,
   not hiding it); default behaviour byte-identical to `TELEMETRY_PROVENANCE.json`.
2. **M2 — L5 verdict** (`.claude/workflows/l5-level-observable.js`). Exit: A10 committed before data;
   ≥ 1 completed run at each N from each start; skeptic reports; `audit/l5_level_runs/REPORT.md`
   with one of {met, not met, inconclusive} and the numbers.
3. **M3 — adoption gate decision** (manual, user). Reads M2. If met + attractor: L2 (Γ₀(N)⁺ fold) and
   L3 (Fricke reflection) become buildable; if not met: the manuscript sentence "the level has no
   observable consequence on this potential" is written and L2/L3 stay closed.
4. **M4 — LeanFlow spec closure** L4 (bridge honest or deleted), L6 (performance story retired or
   logged), L7 (declared-but-unused dependencies dropped). Independent of M1–M3; mechanical.
5. **M5 — next Zenodo version** only when M2 changes a manuscript claim. Always a NEW version under
   concept `22683564`; regenerate and hash-check the bundle zip first.

## Standing rules that shape every milestone

- Nothing is adopted or claimed without its acceptance criterion passing in a committed file.
- `proofs/` is not edited; LeanMaster is read-only from here; every Lean citation carries a SHA.
- No zero-parameter claim, no K3-selection claim, no physical-observable claim (Stream 1 F5b).
- Haiku-first workflows, evidence-bound; producer never verifies itself; negative controls everywhere.
- Push, `reset --hard`, `rm`, and any Zenodo action need explicit confirmation.
