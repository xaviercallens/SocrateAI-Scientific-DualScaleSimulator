# LeanFlow — architecture as built, and a specification for what to build next

**Written 2026-09-21.** Supersedes the architecture prose in `T_duality_Alone.tex` §`sec:leanflow`,
which described a system that does not exist; four of its claims were withdrawn on 2026-09-21
(`audit/STREAM1_BRIDGE.md` S1-F10, S1-F13).

This document has one rule. **The as-built column states only what is verifiable in this repository
today**, with the finding that established it. **Every intended item carries an acceptance criterion a
test could check.** An intended item without such a criterion is not a specification, it is the
aspirational prose the audit just removed.

---

## 1. As built — the verified negative inventory

Each row was measured; see the named finding in `audit/STREAM1_BRIDGE.md` for the command.

| Component | Claimed (withdrawn) | Actually present | Finding |
|---|---|---|---|
| Stiff integrator | "adaptive BDF", `rusty-SUNDIALS` CVODE/IDA | **SciPy `solve_ivp` Radau** and `solve_bvp`. `rusty-SUNDIALS` is **not bundled and not used**; no SUNDIALS in `Cargo.toml` | S1-F13 |
| Inner-loop guards | "AOT-compiled Rust contract assertions" | **None.** All **41** `assert` in `rust_simulator/src/` are inside `#[cfg(test)]`, compiled out of release | S1-F10 |
| Vectorisation | "SIMD-vectorized boundary guards" | **None.** `simd`, `std::simd`, `packed_simd`, `#[repr(simd)]`: no match | S1-F10 |
| Parallelism | — | `rayon` declared in `Cargo.toml`, **never used** (no `par_iter`, no `rayon::`) | S1-F13 |
| Python↔Rust | "Rust/FFI" | **No FFI.** No `extern "C"`, `#[no_mangle]`, `ctypes`, `cffi`, `pyo3`, `libc::`. Driven by `subprocess.run` on a built binary | S1-F13 |
| Python↔Lean | "non-blocking Unix domain socket IPC / C-FFI, JSON-RPC" | `leanflow/bridge/lean_ipc.py`: **blocking `subprocess.run`**, 15 s timeout | S1-F10 |
| Lean verification result | `lean_verified: True` + a named certificate | Was reported **without Lean ever running** (target file absent, fallback returned `True`). Fixed: exit code 0 or nothing | S1-F14 |
| Lean→Rust constants | "Lean 4 Consilience Bridge … exact certified bounds" | `scripts/exact_math_middleware.py` wrote `config.rs` to `../rusty-SUNDIALS/`, **outside the repo, a project that does not exist**. `rust_simulator/` has no `config.rs` and no `FRICKE_Y` | S1-F13 |
| Invariant projections | "embedded in the integrator step" | **Python, applied to the solution array after integration** (`leanflow/core/projections.py`, `solver.py`) | S1-F2 |
| Deployment | "Operating on serverless spot GPU/TPU" | No deployment, no timing harness, no execution log | S1-F13 |

### What *is* real, and works

- **SciPy Radau / `solve_bvp`** integration in `workshopcosmo.py`, plus standalone explicit Rust
  integrators (RK4, semi-implicit Euler) in `rust_simulator/`, run as batch subprocesses.
- **Deterministic, reproducible telemetry.** Seeded end to end; byte-identical across three
  independent reruns; `TELEMETRY_PROVENANCE.json` carries commands, seeds, configs, toolchain, sha256
  (S1-F11).
- **A kernel-anchored `Γ₀(N)⁺` layer.** `leanflow/core/gamma0n_plus.py` transcribes Stream 1's exact
  integer representation with 34 tests pinning it to the Lean statements.
- **Lean checking of discrete invariants**, out of band, via the CLI bridge.
- 136 Python tests, 21 Rust tests, CI green on `lean`, `rust`, `consistency`, `python`.

### Known-broken or misleading, still open

| Item | State |
|---|---|
| `metric_positivity_projection` default `min_y = FRICKE_Y` | Advertised as a positivity guard; actually imposes the level-12 self-dual point — a modelling choice (S1-F3) |
| Clamping instead of reflecting | A trajectory crossing the self-dual point should map by `τ ↦ −1/(Nτ)`, which preserves the dual-scale height; clamping discards the momentum (S1-F3) |
| `modular_domain_fold` | `S`-step disabled by default because `S` is not a symmetry of the implemented potential. **No fold is currently licensed beyond `T`** (S1-F2) |
| `Solver` state layout | Now requires explicit `modulus_indices`; skips with a warning on an unknown layout rather than guessing (S1-F8) |
| The potential itself | A hand-built double well with no modular symmetry; its two stationary points are placed by construction (S1-F2, `audit/K3_SELECTION.md`) |

---

## 2. Specification for what to build next

Ordered by dependency. Each item is written so that "done" is decidable by a test. **None of these is a
physics claim**, and none is registered as a prediction (`PRE_REGISTRATION.md` A5).

### L1 — A potential that is actually modular *(blocks L2, L3, L5)* — **BUILT 2026-09-21, opt-in; not yet adopted**

**Problem.** Every fold, projection and duality claim presupposes an invariance the implemented
potential does not have. Until that is fixed, the "native duality engine" cannot exist, because there
is no duality to be native to.

**Build.** Replace the hand-built double well with a genuine `Γ₀(N)⁺`-invariant function of `τ` — the
natural candidate being a rational function of the Hauptmodul of `X₀(N)⁺`, which exists because that
curve has genus 0 (`scripts/fricke_genus.py`).

**Acceptance criterion.** With `V` the new potential and `N` the declared level:
```python
ok_T, _ = is_invariant_under(V, lambda z: z + 1,            atol=1e-9)   # already true
ok_F, _ = is_invariant_under(V, lambda z: -1/(N*z),         atol=1e-9)   # MUST become true
```
Both `True`, using `leanflow.core.gamma0n_plus.is_invariant_under`. Plus a negative control: a
deliberately non-invariant perturbation must fail the same check.

**STATUS: criterion met.** `leanflow/core/modular_potential.py`, 20 tests in
`tests/test_modular_potential.py`.

| | `T : τ↦τ+1` | Fricke `τ↦−1/(Nτ)` |
|---|---|---|
| new `modular_potential`, `mode="ratio"` | True, `2.2e-15` | **True, `6.7e-16`** |
| new `modular_potential`, `mode="log"` | True, `0.0` | **True, `4.5e-13`** |
| old `compute_potential` (still the default) | True, `3.6e-15` | **False, `32.8`** |

**The construction.** `j` is `SL(2,ℤ)`-invariant, so `W_N : τ ↦ −1/(Nτ)` sends `j(τ) ↦ j(Nτ)` and
`j(Nτ) ↦ j(τ)` — it *swaps* them. And for `g = [[a,b],[Nc,d]] ∈ Γ₀(N)`, `N·gτ = g'(Nτ)` with
`g' = [[a,bN],[c,d]]` of determinant `ad−Nbc = 1`. Hence

```
F_N(τ) = j(τ) + j(Nτ)   is Γ₀(N)⁺-invariant,
```

and any function of it is too. Verified to `~1e-30` at 30 digits, with a negative control that
`S : τ↦−1/τ` (the `N=1` map) does *not* fix it.

**Three things measured and disclosed rather than smoothed over.**
- The first form tried, `V = u/(1+u)`, was invariant but **useless**: `j(Nτ) ~ e^{2πNy}` saturates it
  to `1.0` at `τ=i`, `0.3+0.9i` and `0.1+5i` alike. Hence the two modes.
- `mode="ratio"` has a clean minimum — measured `V ~ r²` with fitted slope **2.000 over four decades**
  and gradient `(0, −5e-9)` — but saturates to `1` away from the well, so the plateau is near
  forceless. `mode="log"` spreads over the domain (`0 → 0.97`) but its minimum is **needle-like**:
  quadratic only for `r ≲ 1e-7`.
- **The invariance is exact but not usable at float64.** `F_N` reaches `1e55` on this domain, so
  computing the transformed point in double precision costs ~13 digits: the same check reads `1e-30`
  in `mpmath` and `1e-13` from a float64 `τ`. A test records this.

**Not adopted, deliberately.** No default changed, no telemetry regenerated, the manuscript untouched
— the Zenodo draft is under review. Switching the simulation over would change every reported
trajectory and is a separate decision. What this unblocks is L2, L3 and L5, which were all waiting on
an invariance that did not exist.

### L2 — Fold by the right group, and only after checking *(depends on L1)*

**Build.** Replace the `SL(2,ℤ)` fold with a `Γ₀(N)⁺` fold built from the exact integer matrices in
`gamma0n_plus.py` (`rho`, `rho_AL`, `fricke_matrix`), and refuse to fold at all unless the target
function passes `is_invariant_under`.

**Acceptance criteria.** (a) The level-`N` self-dual point `i/√N` is a **fixed point** of the fold
(today the `SL(2,ℤ)` fold moves it to `i√12`). (b) `V(fold(τ)) == V(τ)` to `1e-9` over a random sample.
(c) Folding an already-folded point is idempotent. (d) A potential that fails the invariance check
causes the fold to raise, not to run.

### L3 — Reflection instead of clamping *(depends on L1)*

**Build.** Where a trajectory crosses the self-dual locus, map it by the Fricke involution
`τ ↦ −1/(Nτ)` with the momentum transformed accordingly, instead of pinning `y` to `FRICKE_Y`.

**Acceptance criterion.** The dual-scale height `H_N(t) = Nt² + 1/(Nt²)` is preserved across the
reflection to `1e-9` — this is Stream 1's `height_fricke`, already executable as
`gamma0n_plus.height`. And: no trajectory point equals `FRICKE_Y` exactly, which today is the
signature of a clamp. Energy drift across the reflection bounded and reported.

### L4 — Make the Lean bridge honest, or delete it

**Problem (the reporting half is now FIXED, S1-F14).** The bridge referenced
`proofs/KummerLangevinTDA.lean`, which does not exist, and its missing-file and exception branches both
set `lean_verified = True` — so Lean was never invoked and the result still reported success with a
named Lean certificate. Fixed 2026-09-21: no path reports verification without a Lean exit code of 0,
and two tests that asserted the false claim were corrected. **Still open:** the bridge is a blocking
CLI call, not IPC; and the charges it checks are constructed to sum to zero, so even a real Lean build
would only be checking `sum([1,-1,0,…]) = 0`.

**Build.** Either (a) a real asynchronous bridge with a persistent Lean process, or (b) — the cheaper
and more honest option — keep the subprocess call, **remove the fallback that reports success on a
missing file**, and rename the module so it does not say IPC.

**Acceptance criterion.** A missing or non-compiling Lean file yields `verified = False` with the
compiler output attached. A test asserts this by pointing the bridge at a deliberately broken file.
**No path returns `True` without a Lean exit code of 0.**

### L5 — Give the level observable consequences *(depends on L1; prerequisite for any K3 claim)* — **precondition met; criterion NOT met**

**Problem.** `audit/K3_SELECTION.md` §6: changing `N` currently moves `FRICKE_Y` and nothing else, so
**no observable depends on the level and no level can be falsified.**

**Build.** Derive at least one reported quantity from `N` through the dynamics rather than through the
constant alone.

**Acceptance criterion.** Two runs at different `N` produce a stated observable differing by more than
its numerical tolerance, and the difference survives a seed change. Until then, no level-selection
claim may be made in any manuscript — and `audit/K3_SELECTION.md` records level 12 as *fitted*.

**Status.** The *precondition* is met: with L1's potential the level moves the landscape everywhere —
mean `|V₇ − V₁₂| = 0.068`, max `0.099` over 14 random domain points, with no probe where it vanishes,
against exactly `0.0` for the potential currently integrated, which takes no level argument
(`test_level_changes_the_landscape_unlike_the_old_potential`). **The criterion itself is NOT met**,
because it asks for two *runs*, and that needs the adoption decision in L1. The claim stays closed.

### L6 — Retire or implement the performance story

**Build.** Either delete the benchmark table and the GPU/TPU narrative, or produce a timing harness
that writes a committed log.

**Acceptance criterion.** Every speedup figure in any manuscript resolves to a committed log file with
the hardware, command and timestamp recorded, or the figure is removed. Today: `1520×` and `7.1×` have
no log and are already marked "target, pending hardware verification".

### L7 — Use or drop what is declared

**Acceptance criterion.** `rayon` is either used (with a benchmark showing it helps) or removed from
`Cargo.toml`. `exact_math_middleware.py` either writes somewhere the Rust build actually reads, with a
test proving the constant reaches the binary, or is documented as Python-only. (Partly done: it now
reports instead of writing into the void.)

---

## 3. Rules for anyone editing this document

1. **The as-built column changes only when the code changes**, and each change cites the test or
   command that shows it.
2. **No component is described in the present tense before its acceptance criterion passes.** That is
   the specific failure S1-F10 and S1-F13 recorded: "embeds", "operating on", "compiled ahead-of-time"
   were all written about things that did not exist.
3. **Intended design goes in §2 with a criterion**, never in §1.
4. Tier every claim that touches physics, and remember Stream 1's program-wide statement: no exact
   physical observable exists anywhere in this programme.
