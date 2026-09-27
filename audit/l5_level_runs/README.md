# L5 level experiment: partial run, 2026-09-27

Script: `scripts/l5_level_experiment.py`. It is a measurement of spec item L5 (does the level N
reach a reported observable through the dynamics?). It adopts nothing: `workshopcosmo.py`
still uses the hand-built double well.

**Status: INCOMPLETE. L5 is neither met nor refuted.** Claude Code's low-memory guard killed the
run after 20 of 64 runs. Other processes on the VM held about 24 GB; this job was not the cause.
The partial log is `partial_run_2026-09-27.txt`.

What the 20 runs show:

1. **Ratio mode from a generic start cannot test L5.** At the common start tau = 0.2 + 1.2i,
   V = 1.0 with gradient exactly 0 at both N=7 and N=12. This is the flat region disclosed in
   `modular_potential`. With no force the field never moves, so the level cannot reach any
   observable from there. All 8 runs failed within about 5 s: Radau's numerical Jacobian kept
   enlarging its probe step against a right-hand side that never changed, until mpmath
   overflowed. Earlier probes at x ~ 3e12 were a separate problem, fixed exactly by reducing
   x mod 1, since T is in Gamma_0(N).
2. **Log mode from the common start** failed the same way (Jacobian overflow) after 200-260 s.
   All 3 runs that finished were at N=7, rtol 1e-10. None completed.
3. **N=7, ratio, self-dual start** (4 seeds x 2 tolerances, all OK): the field sits in its own well
   at tau = i/sqrt(7) (Im tau = 0.37796). w0 is between -0.22141 and -0.22178 across seeds; the
   rtol 1e-8 vs 1e-10 difference is about 1e-9 or less. No N=12 counterpart finished. In any case
   this start places the level through the initial condition, which is the "constant only" route
   that L5 excludes.

Cost: an rtol 1e-10 run takes 23-29 min on one core, so the full 64-run design takes hours.

Before a rerun: guard the Jacobian probes (clamp Im tau to a bounded range, or pass an explicit
Jacobian), use rtol 1e-9 as the tolerance reference, drop ratio mode from the common start
(already answered: no force), and run serially or on 2-3 processes while memory is short.

## 2026-09-27, second attempt: found why, not yet fixed

Fixed one crash (`math.fmod` on `inf` when a Jacobian probe pushed a coordinate to infinity) by
clamping the evaluation domain to `y in [1e-8, 60]`; `Y_CAP=60` is justified because both modes
are documented to reach their `y -> infinity` asymptote well before that (measured below).
`scripts/l5_level_experiment.py` `_bounded`/`Y_CAP`/`CLAMP_EVENTS`.

That crash gone, log mode from the common start (`0.2, 1.2`) still would not finish — not a
crash this time, a stall (`RuntimeWarning: overflow encountered in multiply` from
`scipy...common.py:345/367` repeating without bound). Traced the real cause:

**`dV/dx` is exactly `0.0` in float64 whenever `|j(N*tau)| >> |ref|`,** for *both* modes, not an
approximation — measured directly on `gamma0_plus_invariant`, not just on `V`. `F_N(tau) =
j(tau) + j(N*tau)` and `|j(N*tau)| ~ exp(2*pi*N*y)` depends only on `y`; its `x`-dependence is a
pure phase. Once `|F_N|` exceeds `|ref|` by more than float64's ~1e-16 relative precision, that
phase cannot survive the float() cast, and `V`, a function of `|F_N - ref|`, loses **all**
`x`-dependence to machine precision. Measured `|F_N|/|ref|` and the `x`-spread of `|F_N|` across
`x in {0.05, 0.2, 0.35}`, both `N`:

| y | N=7 relative x-spread | N=7 \|F\|/\|ref\| | N=12 relative x-spread | N=12 \|F\|/\|ref\| |
|---|---|---|---|---|
| 1.2 | 0 | 2.5e15 | 0 | 3.5e29 |
| 0.8 | 1.3e-12 | 5.8e7 | 0 | 2.8e16 |
| 0.5 | 8.0e-5 | 107 | 3.0e-12 | 4.2e6 |
| 0.4 | 6.4e-2 | 1.24 | 3.7e-7 | 2210 |

So the crossover from "x is dynamically inert" to "x matters" happens within a factor of ~2 in y,
right where the self-dual points sit (Im = 0.378 at N=7, 0.289 at N=12). A Jacobian column that
goes from identically zero to O(1) over that short a range is a textbook stiff kink, and is almost
certainly what stalls Radau's adaptive numerical-Jacobian step search — not a bug in this script,
a genuine feature of the potential's steepness (`j(N*tau)` is a modular function; its q-expansion
is inherently exponential in y).

**This is itself indirect evidence for L5, without a completed ODE run:** at the *same* point
(0.2, 1.2), `dV/dy` already differs between levels: 0.1116 (N=7) vs 0.0747 (N=12), a 33% difference,
in a regime where `dV/dx` is identically zero at both. That is a level-dependent force from the
dynamics, at a point neither level's self-dual constant was used to reach — but it is one point,
not an integrated trajectory, and does not by itself meet the stated L5 criterion (two *runs*).

**Not done, and not attempted further this session:** completing an ODE run through the stiff
transition (needs either a bounded/analytic Jacobian, a max_step small enough to resolve the
crossover, or accepting non-convergence there as a reportable outcome), and the full 64-run design.
Whether to keep debugging the integrator or to test level-dependence more cheaply — by sampling
`dV/dy(y; N)` directly across a grid, no ODE at all — was left to the user rather than decided
here. **User chose the grid scan.**

## 2026-09-27, third pass: grid scan (`scripts/l5_force_scan.py`, no ODE)

Not the L5 acceptance criterion (that needs completed runs); a cheaper, prior question: is there
a level-dependent force anywhere along y, robust to what shouldn't matter? Scanned `dV/dy(x, y; N)`
for `N in {7, 12}`, both modes, `y` log-spaced over `[0.20, 2.00]` (24 points), at 5 arbitrary `x`
probes standing in for "seed". Tolerance at each point is measured, not assumed: the larger of the
spread from halving the finite-difference step and from raising `dps` 25->30, each isolated with
the other held fixed. A cell is flagged only if **every** x-probe exceeds tolerance **and** the
sign of `dV/dy(N=12) - dV/dy(N=7)` agrees across all 5. Full grid: `audit/l5_level_runs/force_scan.json`.

**Result: 28 of 48 (mode, y) cells pass.** Two very different regimes:

- **`mode="log"`: level-dependent from y=0.30 to y=2.00, essentially the whole scanned range**
  (21 of 24 y-values; the 3 misses are y=0.20-0.27, inside the wells' immediate mixing zone,
  where the sign varies by which x-probe — itself informative, not a null result). The difference
  decays smoothly with y (0.70 at y=0.30 down to 0.017 at y=2.00) but stays 6-9 orders of
  magnitude above its shrinking tolerance throughout. This is the mode the earlier single-point
  check (0.2, 1.2) already flagged, now shown to hold over a 7x range in y, not just one point.
- **`mode="ratio"`: level-dependent only in a band, y=0.36-0.90** (8 of 24), matching the
  documented saturation: beyond y~1.0 the difference is machine-zero at both levels (confirms the
  `dV/dx`-plateau finding above extends to `dV/dy` in this mode) and below y~0.33 the sign is
  x-dependent, same caveat as log mode's near-well band.

**Reading this against L5:** it is evidence *for* the level reaching an observable through the
dynamics — robust to x, to two independent numerical-precision checks, and (in log mode) to a
wide range of y — but it is a force-field measurement, not a completed trajectory. It does not
by itself satisfy the stated criterion, and it does not tell you whether an actual field, released
somewhere, spends enough time in the y=0.3-2.0 band for the difference to show up in `w0`/`wa`
rather than washing out. It does say where to aim any further integrator work: log mode, `y > 0.4`
or so, safely past the near-well mixing zone and the `dV/dx`-plateau's own stiff transition.
