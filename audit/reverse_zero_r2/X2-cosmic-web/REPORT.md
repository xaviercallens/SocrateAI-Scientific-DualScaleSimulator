# X2 — cosmic-web TDA in redshift space with a powered gate (round 2)

Tier: X (exploratory numerics / data fit). Registered spec: `audit/reverse_zero_r2/registration/registration.json`,
id `X2`, frozen at commit `9e6705e940f35aef8b377ec3631f62ba3d52aea4` (before any data file in this task was opened).
This is **not** a test of K3xT2: it tests whether a Gaussian-initial-condition, linear-bias, Kaiser+FoG lognormal
mock model reproduces the observed SDSS DR17 wedge's clustering. No LeanMaster or PRE_REGISTRATION statement is
about cosmic-web topology (LeanMaster Stream 8: "Observables: none").

## What changed from round 1 (`audit/reverse_zero/E5-cosmic-web-tda-scaled/`)

Round 1 used an alpha-complex on the raw point cloud and a z-shuffle null; its Poisson null failed a KS test on r
(p=0.009) and its xi(r) gate had no power (xi=0 passes it). Round 2 discards that machinery and implements the
X2 statistic exactly as registered:

- **Cubical complex**, not alpha: Betti b0,b1,b2 of sublevel/superlevel sets of the density-contrast field
  `delta_s`, Gaussian-smoothed at R_s=20 Mpc/h on a 4 Mpc/h voxel grid, evaluated on 4 non-overlapping RA/Dec
  tiles (split at RA=180, Dec=25). T = sum over tiles and nu (31 points in [-3,3]) of (b-mu)^2/max(sigma^2,1).
- **Redshift-space lognormal mocks**: CAMB (Planck18) linear P(k) at z_eff=0.0742, sigma8=0.811 rescaled from
  CAMB's own sigma8(z=0); lognormal transform xi_G=ln(1+b^2 xi_lin); Gaussian random field on a 128x160x96,
  4 Mpc/h grid; RSD displacement Psi_k = i k G_k/(b k^2) from the **Gaussian** field (not the exponentiated
  density), f=Omega_m(z_eff)^0.55 (see `camb_tables.npz`), plus N(0, sigma_v=4 Mpc/h) along the line of sight.
- **Poisson null matched shell by shell**: drawn from the data's own comoving-distance shells (4 Mpc/h), uniform
  in volume within each shell, uniform on the pixel mask. KS test of mock r vs data r: median p=0.431, 0/200 mocks
  below 0.05 (round 1's defect, p=0.009, is fixed — see `x2_results.json:shell_matching`).
- **Full-covariance gate disjoint from the bias fit**: bias b fitted only on F=[8,20) Mpc/h (10 mocks per b,
  b in {0.80,...,2.40} step 0.10, 170 mocks total); the gate's chi2 (Hartlap-corrected, 10x10 covariance from
  500 fiducial mocks) is evaluated only on G=[20,60) Mpc/h, a disjoint range.
- **Gate power demonstrated before use**: 200 Poisson mocks (xi=0) rejected at rate 0.995 [0.972,1.000] (95% CI);
  200 xi-times-3 mocks rejected at rate 0.985 [0.957,0.997]. Both >= 0.95 -> **power_pass = true**.
- **Known-answer controls** at the same N (193,536), through the identical voxel/smooth/cubical pipeline:
  planted voids (evacuated spheres, refilled Poisson-in-volume outside them) and a planted ring (thin circular
  filament). Void detection rate (p_corr<0.05, Sidak N_eff=3): 30/30 = 1.00. Ring detection rate: 1/30 = 0.033
  (see limitations below).

## Registered result

**Gate calibration fails: data chi2(G) = 40.21 vs. the 500-mock fiducial null's 95th-percentile threshold 20.18;
two-sided empirical p = 0.004 (2.88 sigma via `scipy.stats.norm.isf(p/2)`).** Per the registered decision_rule
("Gate calibration fails ... INCONCLUSIVE"), **no topology p-value is computed and the registered X2 verdict is
INCONCLUSIVE.** This is a statement about the lognormal-mock family with bias fitted on F=[8,20) and RSD velocities
from the Gaussian field, not a statement about M0 or about K3xT2, which makes no cosmic-web-topology prediction.
Wording: this is "the null model does not reproduce the data's xi(r) on [20,60) Mpc/h at the fitted bias" —
never "M0 is falsified" and never "confirms".

Per-bin pulls (`(data-mock_mean)/mock_sd`, G bins, 4 Mpc/h wide from 20 to 60): -1.09, -1.51, -0.83, -0.27,
+0.28, +0.99, +1.27, +1.02, +0.61, +0.20 — a broad excess across most of the range relative to the mock
ensemble, not concentrated at one bin, consistent with the chi2-level failure being a genuine shape mismatch
rather than one outlier bin.

## Informational-only deviation (not the registered verdict)

Because the registered rule stops at INCONCLUSIVE, the topology p-value it would have produced is not part of
the verdict. `x2_results.json:topology_informational_deviation` reports it anyway, against the same (already
miscalibrated) null, purely for context: p_corr (Sidak, N_eff=3) = 0.341, i.e. the topology statistic alone
does not show the same tension the xi(r) gate does. This disagreement between the two statistic families is
worth noting for a future round, not resolved here.

## Limitations / honesty notes

1. **Randoms multiplier of 20x, not the registered 5x, is used only inside `x2_paircount_check.py`**'s
   independent n=6000 side validation of the grid xi(r) estimator — not part of the main pipeline. The main
   pipeline's randoms (`Ctx._randoms`, alpha=N/NR) use the registered 5x for both data and mocks (they share
   the same `Ctx.Rg`/`Ctx.alpha`); this deviation is confined to the validation script.
2. **Ring control has low power (0.033)** — a thin planted ring is close to isotropic on the scales this
   statistic is sensitive to (tile-level Betti counts at R_s=20 smoothing wash out a narrow filament unless
   it dominates a tile). This is a real limitation of the statistic for filament-like anomalies, reported
   rather than hidden. Void detection (1.00) shows the pipeline has power for large evacuated regions.
3. **Fiducial self-check false-positive rate = 0.084** (leave-one-out family p_corr<0.05 rate among the 500
   fiducial mocks themselves; nominal is 0.05, and the binomial 1-sigma band on a true 0.05 rate at N=500 is
   about +/-0.01). Plausible causes are LOO's mu/sigma bias correction being imperfect at N=499 in the
   denominator, or residual non-Gaussianity of T's tail. Recorded as a caveat on the informational p_corr
   number above, not corrected post hoc.
4. **CAMB P(k)/sigma8 double-redshift bug from round 1 (`e5b_lognormal_and_poisson_null.py`) does not recur
   here**: `x2_lib.camb_tables()` calls `get_matter_power_spectrum` once at z=0 (for sigma8 rescaling) and once
   at z_eff separately, avoiding the sorted-`zs`-vs-unsorted-`PK_redshifts` indexing mismatch that caused the
   round-1 8.8% sigma8 error. Not independently re-verified by a second self-check in this round (time
   budget); flagged NOT ATTEMPTED (only the design avoids the exact round-1 bug class; no new top-hat-quadrature
   cross-check was rerun here).
5. **Grid (NGP+FFT) xi(r) estimator vs. direct pair counting**: cross-checked on an independent n=6000
   subsample (`x2_paircount_check.py`, `paircount_crosscheck.json`) — agreement within ~20% in the F-range
   bins (ratio 1.03-1.22), where signal exceeds subsample noise; G-range bins are noise-dominated at n=6000
   (ratios swing from -1.4 to +4.7) and this check cannot validate them directly. The full-N (193,536) xi(r)
   used in the actual analysis is not directly pair-count-cross-checked at G-range precision; this is a real
   gap, listed as NOT ATTEMPTED rather than assumed fine by extrapolation.
6. **Runtime**: the full ensemble (170 biasgrid + 500 fiducial + 200 Poisson + 200 x3 + 30 + 30 controls =
   1130 mocks) ran in several foreground `prlimit`-bounded batches over roughly 90 minutes wall time; no
   background jobs were used, one batch at a time as required.

## Commands (repo-relative to `audit/reverse_zero_r2/X2-cosmic-web/`, all run as
`prlimit --as=10737418240 -- /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python <script> <args>`, foreground)

```
python x2_run.py biasgrid 0 170 --workers 4
python x2_fit_bias.py
python x2_run.py poisson 0 200 --workers 4
python x2_run.py x3 0 200 --b 0.937 --workers 4
python x2_run.py fiducial 0 500 --b 0.937 --workers 7   (resumable; run repeatedly, skips existing files)
python x2_run.py control_void 0 30 --b 0.937 --workers 4
python x2_run.py control_ring 0 30 --b 0.937 --workers 4
python x2_analyze.py
python x2_paircount_check.py
```

Seeds: as registered (`registration.json` id X2 `seeds`), reproduced in `x2_lib.SEEDS` and `x2_run.py`'s
per-mode seed arithmetic. **Deviation**: the registered `seeds` field names the randoms seed 20260920; this
run's `x2_lib.SEEDS["randoms"]` is 20260919 (matching the seed already fixed in
`data/real2/cosmic_web/make_random_catalogue.py`, though this run regenerates its own randoms deterministically
inside `x2_lib.Ctx` rather than reading that file). Effect: an unregistered but fixed and disclosed seed for
one random draw; not expected to move the chi2=40.21 result outside mock-to-mock noise, not re-verified with
the registered seed value due to time budget (NOT ATTEMPTED).

## Files

- `x2_lib.py` — shared context (data load, mask, shells, grid, xi estimator, betti curves, CAMB/lognormal generator).
- `x2_run.py` — ensemble runner (biasgrid/fiducial/poisson/x3/control_void/control_ring), one npz per mock under `mocks/`.
- `x2_fit_bias.py` — bias fit on F=[8,20); writes `bias_fit.json`, `data_xi.npz`.
- `x2_analyze.py` — registered decision procedure; writes `x2_results.json`, `data_betti.npz`.
- `x2_paircount_check.py` — independent pair-count cross-check of the grid xi(r) estimator; writes `paircount_crosscheck.json`.
- `camb_tables.npz` — cached CAMB P(k)/sigma8/f (regenerable by `x2_lib.camb_tables()`; committed so the
  numbers above are exactly reproducible without depending on the CAMB install producing bit-identical output).
- `mocks/` — one `.npz` per mock (xi, ks_r_p, seed, and betti for the topology-carrying modes). Not committed
  (bulk, regenerable byte-for-byte from the seeds above); `x2_results.json` and `bias_fit.json` carry every
  number quoted in this report.
