# TODO — DualScaleSimulator

Kept short on purpose; the reasoning is in `ROADMAP.md` and the decisions note. Tick by commit SHA.

## Now (M1, M2)

- [ ] **Run `level-explicit-parameter`** — `Workflow({name: "level-explicit-parameter"})`. Branch
      `loop/level-explicit`. Produces `level=` / `potential=` arguments in `workshopcosmo.py`, a
      derived `FRICKE_Y`, a level-sensitivity test, byte-identical default telemetry.
- [ ] **Write and commit A10** (the L5 registration) — the workflow's first phase does this and
      refuses to run data if the commit is missing. Review the addendum text before the run phase.
- [ ] **Run `l5-level-observable`** — `Workflow({name: "l5-level-observable"})`. Needs ≥ 12 GB free
      (check `free -g`; the 2026-09-27 run was killed by the memory guard while ollama held 15 GB).
      Expected wall time 1–3 h on 3 processes.
- [ ] **Save the L5 report** — subagents cannot write `.md`; the orchestrator writes
      `audit/l5_level_runs/REPORT.md` from the returned `report` field.

## Next (M3, M4)

- [ ] **M3 decision** (user): adopt or not, from the L5 report. Record the answer in
      `specs/LEANFLOW_ARCHITECTURE.md` L1 and in `audit/K3_SELECTION.md` §5.
- [ ] **L4** — `leanflow/bridge/lean_ipc.py`: either point it at a real lake target with an exit-code
      gate (S1-F14 is fixed, but the charges checked still sum to zero by construction), or delete it.
- [ ] **L6** — delete the benchmark table / GPU narrative from the manuscript sources, or commit a
      timing log with hardware, command and timestamp.
- [ ] **L7** — drop `rayon` from `rust_simulator/Cargo.toml`; sweep other declared-unused deps.
- [ ] **`workshopcosmo.main()`** prints "100% CONSILIENCE" unconditionally; make it read the axiom
      audit (5/83 non-standard) and the run results.

## Later

- [ ] `reverse-to-zero-v3` with the staged r3 data (Planck plik-lite, DESI DR1 LSS with randoms,
      SMICA) once the X2 mock gate has power.
- [ ] Bound `c4_pta_product` unconditionally (X3 gave a CURN-conditional interval only).
- [ ] TopoDB: decide whether the 557-run corpus moves to the LeanFlow MCP repo.
- [ ] Next Zenodo version only if a manuscript claim changes (M5).

## Done this cycle (for the record)

- [x] 2026-09-28 six `loop/*` branches merged and pushed (`7a43244`); release `v0.2.0-audit-2026-09-28`.
- [x] 2026-09-27 L5 force scan: log mode level-dependent over y = 0.30–2.00 (`1fe526f`).
- [x] 2026-09-27 `scripts/session_status.sh`; Stream 2 brief filed (`a6daf7b`).
