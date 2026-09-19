# Improvement proposals (grounded in this round's results)

Per the FRAMING RULE: none of these claim K3×T2 predicts or fixes
`mu_sym`, `c4_pta_product`, or Ω_Λ. Each proposal is either (a) a concrete
data-procurement step that could bound one of the 2 remaining free
parameters, or (b) a pre-registration step that would convert this
round's exploratory TDA signal into something falsifiable next round.

## 1. Close the PTA sector (`c4_pta_product`) -- the only route that can
   move this parameter

E2 this round established that no angular HD correlation table exists
publicly (checked, ABSENT, URL trail in the manifest and in
`E2_RESULT.md`). The exact-arithmetic sensitivity scan
(`e2_pta_result.json`) gives a concrete spec: a binned Γ(θ) table with
covariance, precise to ~0.6% of the HD curve's dynamic range per bin,
would be enough to detect `|c4_pta_product| ≳ 0.006` (the 1% threshold);
cruder tables only bound it above ~0.03-0.06.

**Proposal:** track NANOGrav 15yr / IPTA DR3 for a released HD table
with covariance (not just a figure). If the team wants a number sooner,
running `enterprise_extensions`'s binned-correlation pipeline on the
already-downloadable ~1 GB raw NANOGrav 15yr dataset (~18 min per the
manifest) is a defined, boundable side project -- it is a *derived*
computation from a public raw dataset, not a literal answer, so it
should be scoped and pre-registered (expected output format, expected
runtime, exact repo/commit of `enterprise_extensions`) *before* running
it, so its result cannot be read as tuned after the fact.

## 2. Close the screening sector (`mu_sym`) -- currently untouched this
   round

`PRE_REGISTRATION.md:64` lists lunar laser ranging, MICROSCOPE-type
equivalence-principle tests, atom interferometry, and binary pulsars as
the deciding experiments for `mu_sym`. None of these appear in this
round's USABLE DATA list (which was cosmology/PTA/cosmic-web focused).

**Proposal:** the next data-fetch round should target one of these
machine-readable products specifically -- e.g. the published
MICROSCOPE Eötvös-parameter bound (a single number + uncertainty, likely
already citable) or a public lunar-laser-ranging equivalence-principle
constraint -- and run the same style of exact-arithmetic sensitivity
scan on whatever local formula ties `mu_sym` to a screening observable
(analogous to `compute_screening_observable` in
`scripts/param_loop_sim.py`), before claiming any bound.

## 3. Pre-register a TDA statistic and threshold before the next data
   pass

This round's TDA extension found a statistically suggestive (p<1e-3,
Tier X, n=8, not fully independent draws) difference in total H1/H2
persistence between the real SDSS DR17 comoving cloud and its
radial-shuffle null. Right now this is purely descriptive because
`PRE_REGISTRATION.md` has no TDA threshold at all.

**Proposal, concretely:** before the next SDSS/CMB TDA round, freeze in
`PRE_REGISTRATION.md`:
  - the exact statistic (e.g. "total H1 persistence of an AlphaComplex
    on an N=1200 seeded subsample of a flat-ΛCDM, Ωm=0.31115 comoving
    embedding, rescaled to unit diameter"),
  - the exact null (the shuffled-z catalogue, same construction as
    `make_random_catalogue.py`),
  - the exact seeds and N_SEEDS (bump to N_SEEDS≥30 with
    non-overlapping subsamples -- i.e. partition the 193,536 rows into
    disjoint blocks of 1200 rather than resampling with replacement
    across seeds -- to get genuinely independent draws before trusting
    a p-value),
  - a decision rule (e.g. "real total persistence more than 2σ below the
    null distribution's mean, replicated on both H1 and H2, counts as
    X"),
  - a pre-specified confound check: this round verified raw
    (pre-rescale) bounding diameter is statistically indistinguishable
    between real and null draws (p=0.57), which excludes the specific
    "unit-diameter rescaling inflates whichever side has smaller raw
    extent" artifact; the next round should keep recording per-draw
    diameter and mean-NN distance and re-check this before trusting any
    persistence comparison.
  This converts "we saw a suggestive difference" into a test that can
  actually be failed next round, which is the gap E4's own verdict
  names explicitly: *"This would be a new hypothesis (H'), to be frozen
  before comparison. No such prediction is registered in
  PRE_REGISTRATION.md."*

## 4. Extend the TDA to the fetched CMB pair (not yet run this round)

`wmap_ilc_9yr_v5.fits` + `wmap_kq85_mask` (both fetched, sha256-verified,
gitignored) have not yet been run through any TDA pipeline in this
worktree. A natural next step, in the same spirit as #3, is a masked
HEALPix (via `healpy`, already in the venv) persistent-homology or
Euler-characteristic-curve analysis of the WMAP ILC map, with a matched
Gaussian-random-field null generated from the same power spectrum
(`camb`, already in the venv) -- again pre-registered (statistic +
threshold) before running, per #3's lesson.

## 5. Framing discipline for any future octad/Kummer comparison

If a future round wants to actually test the E3/E4 Golay-octad/Kummer
8+16 structure against cosmic-web or CMB topology, the LeanMaster
verdict itself specifies the missing step: *"a TDA prediction would be:
'given cosmological data, the topology of the large-scale structure
encodes this octad partition.' This is frozen as Tier C and requires a
fresh comparison."* Concretely: state the specific topological signature
the octad/complement split would imply (e.g. a specific Betti-number
ratio, or a specific count of independent structures) as a numerical
prediction, write it into `PRE_REGISTRATION.md` BEFORE running the
comparison, then run it once. Until that numerical prediction is
written down first, no TDA result (including this round's) can be
read as evidence either way -- which is why this round's result was
kept strictly descriptive.
