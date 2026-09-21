# Stream 1 bridge — what the newest Lean 4 formalization licenses, and what it refutes here

**Written 2026-09-21.** This is an audit note. **No prediction is registered by this document**
(cf. `PRE_REGISTRATION.md` A5: a statistic, null, seeds and decision rule must be filed before the
relevant data are loaded). Nothing here is a physical claim.

## Pins

Quotations below were taken with `grep`/`sed` from the working trees at these commits. They were
**not** taken from any agent's summary, and neither Lean project was rebuilt for this note
(their own recorded gate results are cited *as theirs*).

| Project | Path | Commit | Tag |
|---|---|---|---|
| **Stream 1** — DualScaleTopologicalUniverseModel-LeanProposal | `~/SocrateAI-DualScaleTopologicalUniverseModel-LeanProposal` | `bb74acb56f386a97e433f94eb0b2632ed03bc4ca` | Zenodo concept 10.5281/zenodo.22853239 |
| **LeanMaster** — Agora-LeanMaster | `~/SocrateAI-Scientific-Agora-LeanMaster` | `ede49f06800cde8177867c02bb08d44f7be275c5` | v3.44.0 |
| This repo | `~/SocrateAI-Scientific-DualScaleSimulator` | `f0747ee` + working tree | — |

Caveats on the pins, stated because they matter:
- Stream 1's working tree at `bb74acb` is **dirty** (`M Agora/Geometry.lean`, `M README.md`,
  `?? Agora/Geometry/EmbeddingAssembly.lean`). Every file quoted below — `Agora/Geometry/SelfDual.lean`
  (blob `989c3e33`), `Agora/Geometry/SymSquareForms.lean`, `Agora/Geometry/MnLattice.lean`,
  `Agora/Geometry/ModularAction.lean` — is **clean at `bb74acb`**, checked with `git status --porcelain`.
- LeanMaster's tree carries one untracked file (`DualScaleDyons/FrickeRepair.lean`), not quoted here.
- LeanMaster went v3.42 → v3.43 → v3.44 on 2026-09-20 and **v3.44 amended a v3.43 claim**. Cite the SHA,
  never "at HEAD".

## Note on which "Stream 1"

Two different things carry that name. LeanMaster's internal Stream 1 is its 2026-09-15
`StringTheoryFormalization` Mathlib integration (`docs/STREAM1_COMPLETION_REPORT.md`). **That is not
what this note is about.** Here "Stream 1" means the separate repository above, the one LeanMaster's
two most recent commits (`f1a40dd`, `ede49f0`, both 2026-09-20) reconcile with. It is the most recent
Lean 4 formalization in the programme, it reports **zero `sorry`**, and it and LeanMaster
independently re-verified each other's theorems rather than accepting them on report.

---

## Part 1 — What Stream 1 proves that bears on this repository

Verbatim from `Agora/Geometry/SelfDual.lean` at `bb74acb`:

```lean
/-- **The self-dual locus is exactly where the root becomes orthogonal to the
    period.** -/
theorem root_orthogonal_iff_selfdual (N τ : K) :
    bN N rootK (period N τ) = 0 ↔ N * τ ^ 2 = -1

/-- Height on the imaginary axis `τ = it`: `H_N(t) = Nt² + 1/(Nt²)`. -/
noncomputable def height (N t : ℝ) : ℝ := N * t ^ 2 + (N * t ^ 2)⁻¹

/-- `H_N ≥ 2` — this IS LeanMaster's circle dual-scale bound at `R = Nt²`. -/
theorem height_ge_two (N t : ℝ) (hN : 0 < N) (ht : 0 < t) : 2 ≤ height N t :=
  DualScaleStream2.DualScale.circle_effective_scale_ge_two (N * t ^ 2) (by positivity)

/-- The height is Fricke-invariant: `H_N(1/(Nt)) = H_N(t)`. -/
theorem height_fricke (N t : ℝ) (hN : 0 < N) (ht : 0 < t) :
    height N (1 / (N * t)) = height N t

/-- Equality holds exactly on the self-dual locus `Nt² = 1` (i.e. `Nτ² = −1`). -/
theorem height_eq_two_iff (N t : ℝ) (hN : 0 < N) (ht : 0 < t) :
    height N t = 2 ↔ N * t ^ 2 = 1
```

Three consequences that this repository can use, and that it could not state precisely before:

1. **The self-dual / Fricke point of level `N` is `τ = i/√N`** — `N·τ² = −1`. The relevant involution
   at level `N` is the **Fricke** involution `t ↦ 1/(Nt)`, under which the height is invariant
   (`height_fricke`). It is `S : τ ↦ −1/τ` **only when `N = 1`**.
2. The dual-scale bound of the whole programme (`tr G + tr G⁻¹ ≥ 2d`, here the circle case) **is** the
   statement that the height is `≥ 2`, with equality **exactly** on the self-dual locus. Stream 1 gets
   this by literally invoking LeanMaster's theorem, so the two projects agree by construction, not by
   coincidence.
3. Also proved, on the `h`-line: `zOf_fricke` (`z(1/(49h)) = z(h)`), `fricke_fixed_iff`
   (`h = ±1/7` are exactly the fixed points), and `s7_singular_points_are_selfdual` — the finite
   singular locus `{−1, 1/27}` of the `s₇` operators **is** the image of the Fricke fixed points.

Two further results, quoted because they are binding on any future work in `proofs/`:

```lean
-- Agora/Geometry/SymSquareForms.lean
theorem sym2_contravariant (M M' : Matrix (Fin 2) (Fin 2) K)   -- sym2 (M · M') = sym2 M' · sym2 M
theorem sym2_not_covariant                                      -- negative control
theorem G0N_det_ne_TN_det (N : ℤ) (hN : 1 ≤ N) : (G0N N).det ≠ (TN N).det
theorem no_isometry_G0N_TN (N : ℤ) (hN : 1 ≤ N) (P : Matrix (Fin 3) (Fin 3) ℤ) ...
```

- **`Sym²` is an ANTI-homomorphism.** It acts by substitution into the form, `Q ↦ Q ∘ M`, and
  substitution reverses composition. The covariant form is **false**, with an explicit negative
  control on `M = ![![1,1],![0,1]]`, `M' = ![![1,0],![1,1]]`. LeanMaster (`ede49f0`) re-derived this
  in sympy before adopting it and calls it "a trap in our own file, caught by Stream 1".
- **There are two different rank-3 lattices, both signature (2,1), and they are NOT isometric:**
  `⟨1⟩ + U(2N)` with `det = −4N²`, and `U + ⟨2N⟩` with `det = −2N`. Conflating them under one name is
  the error LeanMaster amended in v3.44.0. An isometry has determinant `±1` and so preserves the Gram
  determinant, which is the whole proof.

## Part 2 — What Stream 1 explicitly FORBIDS

Quoted verbatim from Stream 1's `README.md` at `bb74acb` (lines 77–81):

> **Blocked, program-wide (Tier C).** No exact physical observable exists anywhere in this program
> (F5b). The Sym² relation supplies **no physical coupling** (VISION §1.3). In particular the
> coincidence that one integer matrix is both the Fricke involution and the Narain T-duality
> generator is a fact about a lattice isometry and **not** a physical identification. Stream 3
> independently found that T-duality is homologically invisible, which is consistent.

This reaches the same place, from the mathematics side, as this repository's own
`PRE_REGISTRATION.md` A7 ruling that "zero free parameters by derivation" is **closed as not
reachable today**. It is recorded here so that a future session cannot re-open `Sym²` as a route to
parameter reduction: **the Sym² relation is not a coupling and yields no free-parameter reduction.**
Stream 1 also records the retraction "ρ = 4 / T = 18; any Kodaira reading of the L₂/L₃ exponents" —
never cite those.

---

## Part 3 — Findings against this repository

Every measurement below was reproduced in-session; the commands are given so they can be re-run.

### S1-F1 — The paper and the code disagree about which point the Fricke point is. The paper is wrong. *(defect, manuscript)*

`papers/T-dulaity alone/T_duality_Alone.tex` says, twice, that the Fricke point is `τ = i`:

- line 50: "… between the Fricke ($\tau = i$) and Orbifold ($\tau = e^{i\pi/3}$) attractors."
- line 383: "… the Fricke involution point $\tau = i$ and the Kummer orbifold fixed point
  $\tau = e^{i\pi/3}$ are recognized as dual geometric attractors …"

The code uses `FRICKE_Y = 1.0 / math.sqrt(12.0)`, i.e. `τ = i/√12`, in five places
(`workshopcosmo.py:43`, `leanflow/core/projections.py:24`, `leanflow/core/solver.py`,
`parameter_sweep_v2.py:78`, and generated into Rust by `scripts/exact_math_middleware.py:31`).

By `root_orthogonal_iff_selfdual`, `τ = i` is the self-dual point at **`N = 1`** and `τ = i/√12` is the
self-dual point at **`N = 12`**. Measured: `1/FRICKE_Y² = 11.999999999999996`. The two statements are
incompatible.

**The code is right and the paper is wrong.** The implemented potential is stationary at the code's
point and not at `τ = i`:

```
tau = i/sqrt12   V= 1.00000000e+00   dV/dx=-0.000e+00   dV/dy= 0.000e+00   <- stationary
tau = i          V= 1.50598306e+00   dV/dx=-0.000e+00   dV/dy= 1.423e+00   <- NOT stationary
tau = e^{i pi/3} V= 1.00000000e-02   dV/dx=-5.091e-16   dV/dy= 4.329e-33   <- stationary
```

Level 12 is independently what the repository intends elsewhere: `specs/spec phase 1.md` carries an
"EXTREMAL LEVEL-12 ETA-QUOTIENT" sub-article throughout.

The likely origin of the error: `τ = i` and `τ = e^{iπ/3}` are the two elliptic points of `SL(2,ℤ)`,
a natural-looking pair, so the prose describes an `SL(2,ℤ)` picture the code does not implement.

### S1-F2 — The modular fold uses the wrong group, and demonstrably destroys the model's own critical point. *(defect, numerics)*

`leanflow/core/projections.py:56` `modular_domain_fold` folds with the `SL(2,ℤ)` generators
`S : τ ↦ −1/τ` and `T : τ ↦ τ ± 1` into `F = {|τ| ≥ 1, |Re τ| ≤ 1/2}`. But `|i/√12| = 0.2887 < 1`, so
the model's own Fricke saddle **is not in `F`** and the fold moves it:

```
Fricke saddle tau=i/sqrt12   (0.000000,0.288675) -> (-0.000000,3.464102)  folds=1
orbifold tau=e^{i pi/3}      (0.500000,0.866025) -> ( 0.500000,0.866025)  folds=0
```

The saddle is teleported to `i√12`. Per `height_fricke`, the involution that fixes a level-`N`
self-dual point is the **Fricke** involution `τ ↦ −1/(Nτ)`, not `S`.

**Worse, and independent of the level: `S` is not a symmetry of the implemented potential at all.**
`compute_potential` is a hand-built double well,

```python
v = a_pot*cos2 + b_pot*sin2 + cos2*(y - FRICKE_Y)**2 + sin2*(y - ORBIFOLD_Y)**2
```

which is periodic in `x` but has no `S`-invariance. Measured over 8 random points in `F`:

```
T-invariance |V(tau+1) - V(tau)| = 6.66e-16          <- holds, to machine precision
max |V(S tau) - V(tau)|          = 8.77              <- fails; ratios V(S tau)/V(tau) from 0.10 to 5.37
```

`solver.py:138-146` applies this fold **to the solution array after integration**, rewriting
`y_out[0,i]` and `y_out[1,i]`. Since `V(Sτ) ≠ V(τ)`, the fold silently moves reported trajectory
points to physically inequivalent points of the potential. Any trajectory that leaves `F` and is
folded back is corrupted, not symmetrised.

The paper's line 383 describes this as a feature — "Rather than interpreting large excursions as
unphysical runaways, the solver folds the trajectory back into `F`" — so the manuscript claim rests
on an invariance the potential does not have.

### S1-F3 — The `FRICKE_Y` floor inside `modular_domain_fold` is unreachable dead code. *(defect, minor)*

`projections.py:91` ends the fold with `y = max(FRICKE_Y, y)`. After a successful fold `y ≥ √3/2 ≈
0.8660 > 0.2887`, so the clamp can never bind. Measured over 200 000 random inputs: minimum folded
`y = 0.866605917`, clamp bound **0** times.

The separate `metric_positivity_projection` (default `min_y = FRICKE_Y`, applied at
`solver.py:150-151`) **does** bind. Its docstring says it enforces "`tau_im = y > 0`" and
"Precludes coordinate runaways and division-by-zero singularities", but the floor it actually
imposes is the level-12 self-dual point — a modelling choice, not a positivity guard. Note also that
clamping is not the involution: a trajectory crossing the self-dual point should be mapped by
`τ ↦ −1/(Nτ)` (which preserves the height, `height_fricke`), not pinned, which discards the momentum.

### S1-F4 — `t_duality_calculus.py` asserts a Sym² result that Stream 1 refutes as stated. *(defect, paper-support code)*

`papers/T-dulaity alone/t_duality_calculus.py:40-46`:

```python
E_0 = sp.Rational(-425, 6)
c_eff = 1 - 24 * E_0
if c_eff == 1701:
    print("Symmetric square modular invariant verified for L_3 = Sym^2 L_2.")
```

Two independent problems:

1. **Non-sequitur.** `1 − 24·(−425/6) = 1701` is an arithmetic identity about two hard-coded rationals.
   It verifies nothing about `L₂`, `L₃` or `Sym²`. The test can never fail.
2. **The claimed relation is wrong as stated.** Stream 1 proves the operator relation with a
   **non-trivial prefactor**: `L₃ = P₂ · Sym²(L₂)` where `P₂ = 1 − 26z − 27z²` (`s7_P2_eval`,
   `s7_P2_discriminant`). `P₂ ≢ 1` — it is exactly the polynomial whose roots `{−1, 1/27}` are the
   images of the Fricke fixed points. "`L_3 = Sym^2 L_2`" drops the prefactor.

The same file closes with `"Calculus verification complete. Hallucination bounds are strictly zero."`

The manuscript itself was **already corrected** by the tier-honesty pass — `T_duality_Alone.tex:377`
now says only that `c_eff = 1 − 24 E_0 = 1701` "is an exact arithmetic identity" and admits no solver
record of it. The script is a survivor of that cleanup. Separately, no source is recorded anywhere in
this repository for `E_0 = −425/6`; it is labelled "Fractional Pole" with no citation.

### S1-F8 — The solver projected the WRONG state components. *(defect, numerics, found while measuring S1-F2)*

`solver.py` applied both projections at state indices `(0, 1)`, and
`projections.py`'s docstring said "Assumes y is at index 1 in `[x, y, vx, vy]` or `[x, y]`". But the
canonical cosmology state of `workshopcosmo.cosmology_rhs` is documented in its own docstring as

```python
y_vec = [a, x, y, u, v] where u = dx/dt, v = dy/dt
```

so index 0 is the **scale factor** and `Im τ` is at index **2**. Applied to that state the solver

- fed `(a, x)` to `modular_domain_fold` as if it were `(Re τ, Im τ)`, T-folding the scale factor by
  `round(a)`. Measured on a `t ≤ 200` run: **`modular_folds_count = 1 984 560 473`** — the count is
  large precisely because `a` grows;
- clamped `y_out[1]`, i.e. **`Re τ`**, up to `FRICKE_Y = 1/√12`. `Re τ` is legitimately `0` at the
  Fricke point and `0.5` at the orbifold point, so the floor is meaningless there and actively wrong;
- never touched the real modulus at index 2.

**Fix.** `Solver` takes `modulus_indices=(re, im)`. Left unset it uses `(0, 1)` only for the
2- and 4-component layouts the docstrings name, and otherwise **skips the projections with a
`RuntimeWarning`** rather than guessing. With `modulus_indices=(1, 2)` the same run gives
`Re τ ∈ [−0.5000, 0.5000]`, `Im τ ∈ [0.2894, 0.9672]` — spanning the Fricke saddle `0.2887` up
toward the orbifold point `0.8660` — and leaves the scale factor alone (`max a = 1.51e4`).

### Blast radius of the S1-F2 and S1-F8 changes — measured

**No committed simulation artifact is affected.** Every tracked output
(`swampland_geodesic_telemetry.csv`, `kummer_langevin_pointcloud.csv`, `parameter_sweep_results.json`,
`tda_mapper_skeleton.json`, `vacuum_decay_cdl_summary.json`, and the `*_summary.json` files) is written
by `workshopcosmo.py`, `parameter_sweep*.py`, `scripts/tda_mapper.py` or the Rust simulator. Checked in
a clean interpreter:

```
$ python3 -c "import workshopcosmo, sys; print([k for k in sys.modules if k.startswith('leanflow')])"
[]        # workshopcosmo imports NO leanflow module; it calls scipy.solve_ivp directly
```

`leanflow.core.solver.Solver` — the only code path that folds — is used only inside the `leanflow`
package, its tests, `scripts/generate_colab_notebook.py` and the Hugging Face deployment copy.

**On that path the change is large, not cosmetic.** Replaying the old `apply_S=True` behaviour over
the trajectory the solver actually produces:

| horizon | points | points the old S-fold would move | max \|Δx\| | max \|Δy\| |
|---|---|---|---|---|
| `t ≤ 5`   | 155  | 155 (100%) | 6.36e-01 | 3.175 |
| `t ≤ 50`  | 1419 | 1419 (100%) | 7.16e-01 | 3.175 |
| `t ≤ 200` | 2831 | 2831 (100%) | 7.16e-01 | 3.175 |

Every point was being moved, by up to `3.175` in `y`. Anyone re-running the LeanFlow solver path
should expect different numbers from before 2026-09-21, and the old numbers should not be trusted.

### S1-F9 — Other copies of the manuscript carried the same error; one is stale in general. *(fixed here; a re-sync decision is left to the user)*

`τ = i` did not appear only in `papers/T-dulaity alone/T_duality_Alone.tex`. Grepping the tracked tree
(excluding `.claude/worktrees/`, which is untracked scratch) found it in **ten** places across four live files, all corrected to `τ = i/√12` on
2026-09-21. The first pass on the main manuscript caught only 2 of its 5 sites, so the grep was
re-run across the whole tracked tree rather than one glob:

| file | sites |
|---|---|
| `papers/T-dulaity alone/T_duality_Alone.tex` (3 sites missed by the first pass) | 3 |
| `zenodo_bundle/T_duality_Alone.tex` (tracked) | 5 |
| `papers/T-dulaity alone/submission/COVER_LETTER_SCIPOST.tex` | 1 |
| `papers/T-dulaity alone/submission/COVER_LETTER_JHEP.tex` | 1 |

**Resolved 2026-09-21 on the user's instruction.** `zenodo_bundle/T_duality_Alone.tex` was a **stale
pre-tier-honesty snapshot**: it differed from the corrected paper in more than the `τ = i` error, and
still carried `\min(y) = 0.289`, a number the tier-honesty pass withdrew as untraceable to any
committed output (the cover letters already recorded that withdrawal; this copy did not). The bundle's
`.tex` and `.pdf` have now been **re-synced from `papers/T-dulaity alone/`** and are byte-identical to
the corrected manuscript. **No Zenodo deposit was made** — that remains a separate decision.

### S1-F10 — The "Dual-Tier Execution Architecture" paragraph did not match the code. *(defect, manuscript; the last open item from the recovery notes)*

This paragraph (`T_duality_Alone.tex` §"Dual-Tier Execution Architecture") had been flagged as
unverified since 2026-09-18. Audited 2026-09-21 against the code; four of its claims are false and are
now withdrawn in a dated note in the manuscript.

| Claim in the paragraph | What the repository contains |
|---|---|
| "compiled ahead-of-time into Rust contract assertions" in the inner loop | **None.** All **41** `assert` statements in `rust_simulator/src/` are inside `#[cfg(test)]` modules — measured 0 before the test module in each of the five files — so they are unit tests, compiled out of release builds |
| "SIMD-vectorized boundary guards" | **Absent.** `simd`, `std::simd`, `packed_simd`, `#[repr(simd)]` return no match anywhere in `rust_simulator/src/` |
| `τ_im > 0` and `w(t) ≥ −1` enforced "during implicit time-stepping" in Rust | Implemented in **Python** (`leanflow/core/projections.py`) and applied to the solution array **after** integration. Defective in three ways: S1-F2, S1-F3, S1-F8 |
| "non-blocking Unix domain socket IPC / C-FFI using structured JSON-RPC payloads" | `leanflow/bridge/lean_ipc.py` makes a **blocking** `subprocess.run` CLI call with a 15 s timeout. No socket, no C-FFI, no JSON-RPC anywhere in the repository |

Also recorded: the Rust crate is not an inner loop of the Python solver at all — `workshopcosmo.py`
invokes a built binary through `subprocess.run` as a separate batch process. The paragraph's own
closing sentence already conceded that the socket-IPC gate "has not been exercised end-to-end in this
repository"; the corrections above go further, because the mechanism described does not exist.

### S1-F5 — Convention now binding on any future `proofs/` work. *(no change needed today)*

`Sym²` is **contravariant**. Any future Lean or Python code in this repository that composes a `Sym²`
lift must use `sym2(M·M') = sym2(M')·sym2(M)`. A grep found no place where this repository currently
composes two `Sym²` lifts, so nothing is broken today — this is recorded to prevent the trap being
walked into, as it was in LeanMaster.

### S1-F6 — The two rank-3 lattices must not be conflated. *(no occurrence found)*

Grepped for a rank-3 lattice / transcendental-lattice construction in this repository: none found.
Recorded so the distinction (`det −4N²` vs `det −2N`) is available before any is written.

### S1-F7 — Untracked scratch Lean files that would pass every gate while proving nothing. *(disposition is the user's call)*

`test_sym2.lean` and `test_sym2_axioms.lean` (untracked, at the repository root) contain:

```lean
def dimSym2A1 : Nat := (dimA1 * (dimA1 + 1)) / 2
theorem sym2_A1_dimension_is_4095 : dimSym2A1 = 4095 := by ...
```

This is a `Nat` arithmetic identity with no connection to `Sym²` of a quadratic form or of a lattice,
in exactly the area Stream 1 has just corrected. Stream 1's own README warns: "⚠️ **A vacuous theorem
would pass every one of those gates.** It compiles, evades the `sorry` grep, reports the three
standard axioms and locks cleanly … Only reading the *statement* catches it."

These files are **not committed** and their origin is unrecorded (see the session memory), together
with ~20 other `check_*.lean` / `proofs/Test*.lean` / `proofs/Check*.lean` scratch files. **No action
taken**: deleting them or moving them to a scratch directory is the user's decision.

---

## Part 4 — What was checked and found sound

- `ORBIFOLD_X, ORBIFOLD_Y = 0.5, √3/2` is `τ = e^{iπ/3}` exactly, matching the manuscript, and the
  potential is stationary there (`dV/dx = −5.09e-16`, `dV/dy = 4.33e-33`).
- The `T : τ ↦ τ ± 1` half of the fold is a genuine symmetry of the potential (machine precision).
- `T_duality_Alone.tex:377` (the `c_eff` claim) is already honest after the tier-honesty pass.
- No `Sym²` composition and no rank-3 lattice construction exists in this repository to be wrong.

## Part 5 — What this does NOT license

- Nothing here is evidence for or against the dual-scale model. Stream 1 states that **no exact
  physical observable exists anywhere in the programme** and that **`Sym²` supplies no physical
  coupling**; the Fricke/Narain matrix coincidence "is a fact about a lattice isometry and **not** a
  physical identification".
- The agreement in Part 1 item 2 between Stream 1's `height_ge_two` and LeanMaster's circle bound is
  **not independent corroboration**: Stream 1 obtains it by invoking LeanMaster's theorem directly.
- No threshold, statistic or decision rule is registered by this note (`PRE_REGISTRATION.md` A5).

## Part 6 — Dolgachev Thm 7.1: one degree of freedom, and an exact integer Fricke isometry

Added 2026-09-21, prompted by the user. Both halves check out; the tiers differ and that matters.

### 6.1 The Fricke involution is an exact integer matrix isometry — kernel-proved, Tier A

`Agora/Geometry/ModularAction.lean` at `bb74acb`, verbatim:

```lean
/-- `ρ(W)` for the Atkin–Lehner element `W = [[Na, b], [Nc, Nd]]/√N`. The `√N`
    cancels in the symmetric square, so this is again an INTEGER matrix — which
    is why Atkin–Lehner involutions are lattice isometries at all. -/
def rhoAL (N a b c d : K) : Matrix (Fin 3) (Fin 3) K :=
  !![N * d ^ 2, -(c ^ 2), 2 * N * c * d;
     -(b ^ 2), N * a ^ 2, -(2 * N * a * b);
     b * d, -(a * c), N * a * d + b * c]

theorem rhoAL_isometry (N a b c d : K) (h : (N * a * d - b * c) ^ 2 = 1) :
    (rhoAL N a b c d)ᵀ * TNR N * rhoAL N a b c d = TNR N

theorem rhoAL_fricke (N : K) : rhoAL N 0 (-1) 1 0 = !![0, -1, 0; -1, 0, 0; 0, 0, -1]

theorem rhoAL_fricke_eq_neg_swap (N : ℤ) : rhoAL N 0 (-1) 1 0 = -swap

/-- The Atkin–Lehner elements are also determinant-cubes, hence also land in
    `SO(2,1)` when `Nad − bc = 1`. So the WHOLE of `Γ₀(N)⁺`, Fricke included,
    acts by orientation-preserving isometries of `U ⊕ ⟨2N⟩`. -/
theorem rhoAL_det (N a b c d : K) : (rhoAL N a b c d).det = (N * a * d - b * c) ^ 3
```

The Fricke matrix `!![0,-1,0; -1,0,0; 0,0,-1]` carries **no `N` and no `√N`** — the `√N` of the
Atkin–Lehner element cancels in the symmetric square. So the involution is an exact integer isometry
of `U ⊕ ⟨2N⟩` **for every level `N` at once**, of determinant `+1`. Verified here for
`N ∈ {1, 2, 7, 12, 49}` (`tests/test_stream1_bridge.py`): `FᵀT_N F = T_N`, `F² = I`, `det F = 1`.

On the period line it is the Möbius map `τ ↦ −1/(Nτ)` (`rhoAL_mulVec_period_field` at
`(0,−1,1,0)`), which fixes `τ = i/√N` — the self-dual point of `root_orthogonal_iff_selfdual`. At
`N = 12` that is exactly this repository's `FRICKE_Y`.

### 6.2 `X₀(N)⁺` — one degree of freedom. **Literature (Tier L), not kernel-proved**

From Stream 1's `paper/sections/02-preliminaries.tex`, verbatim:

> for the rank-19 lattices `M_n = U ⊕ E₈² ⊕ ⟨−2n⟩`, Dolgachev computes
> `(M_n)^⊥ = U ⊕ ⟨2n⟩` (\cite{Dolgachev1996}, §7) and proves that the coarse moduli space of
> `M_n`-polarized K3 surfaces is the Fricke modular curve `H/Γ₀(n)+` (\cite{Dolgachev1996}, Thm 7.1).

A modular **curve** is one complex dimension. So the moduli space of `M_n`-polarized K3 surfaces
carries **exactly one modulus** — the degree of freedom that survives after the rank-19 polarization
is fixed and the Fricke involution is quotiented out. That is the "one degree of freedom `X₀(N)⁺`".

**Tier discipline, stated because it is easy to lose.** Stream 1 labels the surrounding statement
explicitly, `ModularAction.lean:47`:

> NO physics. That `Γ₀(N)⁺ ≅ O⁺(U⊕⟨2N⟩)/±1` is Dolgachev's theorem (literature, not proved here):
> we prove the inclusion `⊇` constructively, which is the direction the applications need.

So:

| Claim | Tier | Status |
|---|---|---|
| Fricke = exact integer matrix isometry of `U ⊕ ⟨2N⟩`, `det = +1`, all `N` | **A** | kernel-proved, 0 sorry, re-verified here |
| self-dual locus is `Nτ² = −1`, height `≥ 2` with equality exactly there | **A** | kernel-proved, re-verified here |
| `Γ₀(N)⁺ ≅ O⁺(U⊕⟨2N⟩)/±1`; moduli space is `H/Γ₀(N)⁺` (Dolgachev Thm 7.1) | **L** | literature; Stream 1 proves only `⊇` |
| the Hauptmodul identification of this period line with the `s₇` family | **L/B** | Stream 1: "NOT kernel-proved … PASS(40) exact plus literature" |
| that this one modulus **is** a cosmological degree of freedom of our universe | **C** | no support anywhere; see below |

### 6.3 What this does and does not say about the parameter count

It is a real sharpening of `PRE_REGISTRATION.md` A7 and worth recording, **as mathematics**:

- A7 closed "zero free parameters by derivation" as **not reachable**. Dolgachev Thm 7.1 says why
  more precisely than "not reachable": for an `M_n`-polarized K3, the mathematics does not have zero
  moduli to offer. It has **one**, and `X₀(N)⁺` is its moduli space. Chasing `2 → 0` was chasing a
  number the geometry does not contain; `1` is the floor on that side.
- The residual modulus has a canonical coordinate (a Hauptmodul of `X₀(N)⁺`) and a canonical
  involution acting on it (Fricke, `τ ↦ −1/(Nτ)`), both now executable here
  (`leanflow/core/gamma0n_plus.py`).

**It does not say the universe has one degree of freedom.** That step needs an identification of this
modulus with a physical field, and Stream 1 states program-wide that there is none
(README, quoted in Part 2): "No exact physical observable exists anywhere in this program (F5b)",
"the Sym² relation supplies **no physical coupling**", and the Fricke/Narain matrix coincidence "is a
fact about a lattice isometry and **not** a physical identification". A7's `mu_sym` has no unit
bridge in any harness. So:

> **`X₀(N)⁺` gives one degree of freedom to the moduli space of `M_n`-polarized K3 surfaces. Reading
> it as the free parameter of a cosmology is Tier C, and no prediction is registered from it here.**

### 6.4 Consequence for the code — already applied

Finding S1-F2 said the `SL(2,ℤ)` fold uses the wrong group. Part 6.1 makes the right one explicit:
at level `N` the group is `Γ₀(N)⁺` and the involution is `τ ↦ −1/(Nτ)`, with the integer matrices
now transcribed in `leanflow/core/gamma0n_plus.py`.

**But swapping `Γ₀(12)⁺` in for `SL(2,ℤ)` would not have been a fix.** The implemented potential is
not invariant under it either — measured `max |V(Fricke τ) − V(τ)| = 35.81`, against `8.77` for `S`
and `6.7e-16` for `T`. `compute_potential` is a hand-built double well,
`a·cos²(πx) + b·sin²(πx) + cos²(πx)(y−y_F)² + sin²(πx)(y−y_O)²`, with the two stationary points
placed by construction; it is periodic in `x` and has no modular symmetry at all. The fix applied was
therefore to make the `S`-step **opt-in and off by default**, leaving only `T`, which is a genuine
symmetry. `gamma0n_plus.is_invariant_under` is provided so the check is run before any future fold.

## Reproducing the measurements

```bash
cd ~/SocrateAI-Scientific-DualScaleSimulator
python3 - <<'PY'
import warnings; warnings.filterwarnings('ignore')
import math, random, workshopcosmo as wc
from leanflow.core.projections import modular_domain_fold, FRICKE_Y
print("N = 1/t^2 =", 1.0/FRICKE_Y**2)                       # 11.999999999999996
for lab,x,y in [("i/sqrt12",0.0,FRICKE_Y),("i",0.0,1.0),("rho",0.5,math.sqrt(3)/2)]:
    print(lab, wc.compute_potential(x,y))                   # stationary at i/sqrt12 and rho, not at i
print("fold(Fricke) ->", modular_domain_fold(0.0, FRICKE_Y))# -> (~0, 3.464102) = i*sqrt12
random.seed(1); worst=0.0
for _ in range(8):
    x=random.uniform(-0.5,0.5); y=10**random.uniform(-0.6,0.8); n=x*x+y*y
    worst=max(worst, abs(wc.compute_potential(-x/n,y/n)[0]-wc.compute_potential(x,y)[0]))
print("max |V(S tau)-V(tau)| =", worst)                     # 8.77 -> S is NOT a symmetry of V
print("|V(tau+1)-V(tau)| =",
      abs(wc.compute_potential(1.17,0.9)[0]-wc.compute_potential(0.17,0.9)[0]))  # 6.66e-16
PY
```
