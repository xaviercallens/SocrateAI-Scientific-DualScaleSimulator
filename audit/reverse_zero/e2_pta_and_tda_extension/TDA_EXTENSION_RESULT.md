# TDA extension: SDSS DR17 cosmic web vs. matched radial-shuffle null

**Tier: X** (exploratory, descriptive, un-preregistered).

## Data

- Real: `sdss_dr17_cosmic_web_galaxies` -- 193,536 SDSS DR17 spectroscopic
  galaxies, ra∈[140,220], dec∈[0,50], 0.02≤z≤0.12 (sha256-verified, exact
  row-count match to an independent `SELECT COUNT(*)`, per manifest).
- Null: `sdss_dr17_random_shuffled_z` -- the SAME 193,536 rows, SAME
  observed ra/dec, with z (and paired zErr) permuted, seed 20260919
  (`data/real2/cosmic_web/make_random_catalogue.py`). This is a
  **methodological upgrade** over the Poisson/box nulls used in rounds
  1-3: it preserves the real angular selection function exactly and
  destroys only the radial (line-of-sight) structure, so any measured
  difference isolates genuine 3D clustering, not a survey-footprint
  artifact.

Neither dataset had been run through TDA in this worktree before this
round (`audit/zero_param_loop/round1-3` used the Pantheon+ SN comoving
cloud, not SDSS).

## Method

1. Comoving Cartesian conversion (flat ΛCDM, Ωm=0.31115, the frozen
   LeanMaster value) via `comoving_distance_flat_lcdm`, imported
   **unmodified** from `audit/zero_param_loop/round2/tda_gudhi.py`
   (reuse, not reimplementation).
2. 193,536 points is not Rips-tractable at H1/H2. Per advisor review,
   used `gudhi.AlphaComplex` (exact 3D Delaunay-based persistence, no
   `max_edge_length` fudge) on seeded subsamples of N=1200 points,
   rescaled to unit diameter.
3. Repeated over 8 independent subsample seeds
   (20260919..20260926) for BOTH real and null, through the identical
   pipeline. Reported the distribution (mean, std, all 8 values), not a
   single draw.

Command: `.venv-tda/bin/python audit/reverse_zero/e2_pta_and_tda_extension/sdss_cosmic_web_tda.py`
Full output: `audit/reverse_zero/e2_pta_and_tda_extension/sdss_cosmic_web_tda_result.json`

## Result

`betti_numbers()` (evaluated at the end of the filtration, i.e. on the
full Alpha complex) is `(1, 0, 0)` for every draw of both real and null
-- expected and uninformative, because a full 3D Alpha complex spans a
contractible region generically. The informative statistic is **total
persistence** (sum of bar lengths, i.e. how much topological structure
exists across scales before dying), summed per homology dimension:

| statistic | real (mean ± std, n=8) | null (mean ± std, n=8) | Welch t | p |
|---|---|---|---|---|
| total persistence, H1 (loops) | 0.3734 ± 0.0177 | 0.4240 ± 0.0221 | -4.74 | 3.6e-4 |
| total persistence, H2 (voids) | 0.0792 ± 0.0040 | 0.0959 ± 0.0057 | -6.30 | 3.3e-5 |

The real catalogue shows systematically **lower** total persistence in
both H1 and H2 than the radial-shuffle null, at this N_SUB=1200,
unit-rescaled filtration, across all 8 seeds.

## Confound check: is this just a rescaling artifact?

Before trusting the sign above, checked whether it is explained by
real and null subsamples having a different realized bounding diameter
(unit-diameter rescaling would then mechanically inflate bar lengths
for whichever side has the smaller raw diameter, producing this exact
sign with zero topological content). Per draw, recorded the raw
(pre-rescale) diameter and the raw mean nearest-neighbor distance:

| statistic | real mean (n=8) | null mean (n=8) | Welch t | p |
|---|---|---|---|---|
| raw diameter (pre-rescale) | 0.14272 | 0.14177 | 0.58 | 0.57 |
| raw mean NN distance | 0.003684 | 0.004157 | -15.1 | 2.6e-9 |
| mean NN distance / diameter | 0.02583 | 0.02933 | -9.0 | 3.3e-7 |

**Raw diameters are statistically indistinguishable (p=0.57)** -- the
specific mechanical confound (differing bounding-box size under
unit-diameter rescaling) is excluded. The real catalogue's mean
nearest-neighbor distance IS significantly smaller than the null's,
even relative to diameter. That is consistent with genuine clustering
(filamentary structure packs points more tightly on average than a
radially-scrambled version at the same footprint and overall extent),
not a separate artifact to correct for -- it is plausibly part of the
same signal that produces the persistence difference, not an
independent confound. This does not promote the result to confirmed;
it removes the one specific mechanical alternative explanation that was
checked.

## What this is -- and is not

This is a reproducible signal in this specific pipeline (exact command
+ seeds above; anyone can rerun it), with the one checked mechanical
confound (diameter mismatch under rescaling) excluded. What it is
**not**:

- **Not a pre-registered test.** `PRE_REGISTRATION.md` contains no
  threshold for any Betti number or persistence statistic, so this
  result is descriptive, not a pass/fail verdict.
- **Not evidence for or against the LeanMaster Golay-octad/Kummer 8+16
  split.** Per the FRAMING RULE and the E3/E4 verdicts (Stream 8), that
  split is Tier C with no registered cosmological prediction. No
  measured Betti number here "matches" or "supports" 8+16 -- that
  sentence is exactly what the framing rule exists to forbid.
- **Not a constraint on `mu_sym` or `c4_pta_product`.**
  `parameter_effect = none`: this result does not touch either free
  parameter of the simulator.
- **Not fully independent draws.** The 8 seeds subsample from the same
  193,536-row parent catalogue without replacement per draw, so
  different seeds' subsamples can share points; n=8 is small. The
  p-values above should be read as "worth pre-registering and
  re-testing with a cleaner design," not as a confirmed result.

The one honest, defensible statement: **a real ΛCDM comoving embedding
of the observed SDSS DR17 footprint has measurably less low-dimensional
persistent topological structure (fewer/shorter-lived loops and voids
at this scale) than the same footprint with redshifts radially
scrambled**, in this specific pipeline. Plausible reading: real
large-scale structure (filaments, walls) creates a more "connected"
point distribution along the line of sight than a radially randomized
version, which produces more transient small-scale loops/voids by
chance -- but this is a plausibility read, not something established
here; it would need a dedicated test to confirm.
