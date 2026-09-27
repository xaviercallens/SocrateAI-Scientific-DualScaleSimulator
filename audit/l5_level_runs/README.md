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
