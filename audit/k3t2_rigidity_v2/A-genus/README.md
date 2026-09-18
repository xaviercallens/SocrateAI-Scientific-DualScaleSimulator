# Track A v2 -- K3 elliptic genus / Mathieu moonshine / twining

## Exact reproduction command

```
cd audit/k3t2_rigidity_v2/A-genus
/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python run_track_a_v2.py
```

Runs from a clean checkout of this directory alone (only imports
`series2d.py`, `thetas.py`, `appell.py`, all three present here, copied
byte-for-byte from `audit/k3t2_rigidity/genus-moonshine/` -- v1's blind
track, reuse permitted by the task). No arguments. Writes `results.json`
and `exports.json` in the current directory (overwriting). Takes ~2s.

## What's in results.json

- `edge_row_truncation_bug`: a defect found in this session (not listed by
  the task) -- `compute_ZK3_and_phi01(cutoff_q, N)` alone is inexact at its
  own top reported q-order; fixed by `compute_ZK3_and_phi01_safe` (computes
  with margin +2, reports only the safe range). Evidence: raw vs
  margin-corrected (n=10, l=+-5) coefficient, cross-checked against the
  n=4, n=6 rows at the same discriminant D=15 (all three must agree by the
  discriminant-dependence theorem; the raw one didn't, the corrected one
  does).
- `H_tau_appell_lerch.A_n_1_to_10`: A_1..A_10 of H(tau), Tier B (exact,
  truncation-stable at cutoff+2 and cutoff+4, two-slice-agreement and
  multiply-back verified).
- `rigidity_scan_a_mu_coefficient`: N in 20..28 selected by z-independence
  of the two y-power slices -- structural, classified RIGID.
- `rigidity_scan_b_overall_factor`: both the Euler-characteristic selector
  (v1's, circular, NORMALISATION) and the task-suggested (q^0,y^1)-
  coefficient selector (also circular on inspection -- the target value 2
  is the RHS of the equation being solved for k) are reported and BOTH
  classified NORMALISATION, honestly, with the reasoning stated in-line.
  No literal-free structural condition pinning k was found.
- `twining_2A_3A_5A_7A`: full twining for all four classes using the
  task-supplied Tier-L chi(g)/F_g inputs (Cheng-Duncan-Harvey Table 3, from
  memory -- NOT re-derived). Structural checks -- integrality and
  A_n^(g) === A_n mod ord(g) for n<=8 -- are done on the **raw** (undivided)
  H_g Fourier coefficients, not on A_n^(g)=H_g[i]/2: dividing by 2 first
  produced spurious "non-integer" results at class 7A that are NOT a
  defect (see below), and testing congruence on raw coefficients sidesteps
  that entirely. **All four classes pass integrality and the mod-ord(g)
  congruence at every reported n<=8** (spot-checked by hand for 7A: raw
  untwined - raw 7A at n=1,4,6,7,8 is 91, 4550, 27832, 61684, 131103, all
  divisible by 7). `A_n_g` is still reported per class (halved), and at
  7A this includes genuine half-integers at n=1,8 (-1/2, -3/2) -- these
  are `Re(chi_45(7A))`-type values: K_1=45+45bar decomposes at 7A into a
  conjugate pair of irrational M24 characters ((-1+-i*sqrt(7))/2), so
  halving the (integer) raw trace can legitimately land on a half-integer.
  Flagged in `A_n_g_non_integer_note`, not treated as a bug.
  The "27720 lock" identity A_2^(g)*60=4*A_1^(g)*77 is computed (not
  skipped) for all four classes and **fails for all of them** -- expected,
  since it is a relation among untwined dimensions, not preserved by
  traces.
- `twining_truncation_stability_cutoff_plus2`: the entire twining
  computation (independent H(tau) extraction, Lambda_N, 1/eta^3) rerun at
  H_cutoff_q+2; raw H_g[n], n=0..10, unchanged for all four classes.
- `negative_control_wrong_class_pairing`: chi=6 (3A's) deliberately paired
  with F_2A (wrong pairing) -- fails integrality and congruence, confirming
  the structural checks have discriminating power (not vacuously true).
- `lambda_N_cross_checks_match_direct_log_derivative`: Lambda_N closed form
  independently re-derived by direct series log-differentiation for every
  N in {2,3,5,7} used; all match.

## Revision note

An earlier pass of this script computed the congruence/integrality checks
on A_n^(g)=H_g[i]/2 (post-division) and reported 7A as failing. An advisor
review caught that this was a misdiagnosis (the irrational-character
halving above, not a bad tier-L input), and separately caught that
series2d's zero-dropping convention meant several checks (congruence at
vanishing coefficients, the identity check for 3A/5A/7A, sign pattern at
zero) were silently skipped rather than computed. Both are fixed in the
current script (raw-coefficient checks + explicit zero-fill); this file
and results.json reflect the corrected run.

## Could not do

See `could_not_do` in results.json: no literal-free structural condition
was found for rigidity scan (b) (the overall factor k of Z=k*phi_{0,1});
both selectors tried are reported as NORMALISATION rather than mislabelled
RIGID. Everything else in the task's Track A scope was attempted and
reported, including the negative/failing results (7A's non-integer A_n,
the 27720 lock failing under every twining).
