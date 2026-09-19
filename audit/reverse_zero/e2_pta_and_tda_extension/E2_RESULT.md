> **Correction notice (2026-09-19):** statements in this file were corrected after the skeptic review. See `audit/reverse_zero/ERRATA.md` and REPORT.md §6. Wherever they differ, REPORT.md is authoritative.

# E2 -- is the PTA sector GR? (NO-DATA branch)

**Tier: X/B** (exact arithmetic on a fixed pre-existing formula; no fit,
no literature bound invented).

## 1. Data check (re-run this session)

```
cd /mnt/disks/disk-socrateai-local-1/dualscale-wt-reverse
find data -iname '*hellings*' -o -iname '*angular*'
```
Output: empty.

```
ls -la data/real/pulsar_timing/
```
Output: only `nanograv_kde_freespectrum.zip` (a ceffyl free-spectrum
posterior KDE -- a function of GW **frequency**, not of pulsar-pair
**angular separation**). It is orthogonal to `c4_pta_product`, which
lives entirely in the angular (Legendre-P4) part of Gamma(theta). It is
NOT used below as an angular constraint.

Both checks confirm, independently of the run manifest, that
`nanograv_15yr_hd_angular_correlation` is **ABSENT**. The manifest
already records the fetch trail:

> "3+ sources checked (nanograv/15yr_stochastic_analysis figure_1,
> PTArcade, Zenodo searches). Only posterior chains and a pipeline
> notebook exist that COMPUTES the binned correlation from a ~GB raw
> dataset (~18min enterprise_extensions run); not reconstructed per
> ground rules."

So the "fit Gamma(theta) to it" branch of E2 does not apply this round.

## 2. Pre-registered rule (quoted verbatim)

> PRE_REGISTRATION.md:65 -- "`c4_pta_product` (ℓ = 4 term in the PTA
> angular correlation) | bound it; a measured correlation consistent
> with pure Hellings–Downs would exclude large values | NANOGrav, IPTA,
> then SKA"

The registered decisive product is explicitly "a measured correlation"
-- exactly the thing confirmed absent in §1.

## 3. What would decide it (concrete procurement spec)

A binned angular two-point correlation table with:
- bin centers in angular separation θ (degrees),
- ρ(θ) or Γ(θ) per bin,
- the **bin-to-bin covariance matrix** (HD estimators from a fixed
  pulsar array are correlated across angular bins; diagonal errors
  alone understate χ² and would bias any c4 interval).

Named candidates (checked, none available with covariance as of this
session):
1. NANOGrav 15yr stochastic-analysis release's own HD-figure data, if
   ever released as a machine-readable table+covariance (currently a
   plotted figure only).
2. IPTA DR3/DR2 combined angular-correlation product, if released with
   covariance.
3. Running the `enterprise_extensions` binned-correlation pipeline on
   the full ~1 GB raw NANOGrav 15yr dataset ourselves (~18 min per the
   manifest) -- **not done this round**, because it is a multi-step
   derived computation, not a fetch of an existing public product. This
   is the concrete follow-up if the team wants a first angular
   constraint before NANOGrav/IPTA publish one.

A figure-only release (points + error bars, no covariance) would be
usable only for a rough diagonal-χ² look, not a defensible 68/95%
interval on `c4_pta_product`.

## 4. Sensitivity scan (exact arithmetic, no data, no fit)

Using `scripts/param_loop_sim.compute_pta_observable` (imported
unmodified) at `pta_suppression=1.0` (so `c4_c0_ratio == c4_pta_product`
directly), `gamma_theta` is exactly linear in `c4_pta_product` (the
Legendre-P4 term `l4_response` does not depend on it). Fitted slope:
`max_deviation_from_hd / c4_pta_product = 1.000` [CORRECTED from 6.171; the script JSON and skeptic_statistics/check_e2_slope.json give 1.000] (spread across the
scanned grid: see `e2_pta_result.json`, effectively zero -- confirms
linearity to floating-point precision).

Thresholds for `|c4_pta_product|` producing a given fractional
perturbation of the HD curve's own dynamic range across the harness's
15 angular bins:

| perturbation of HD curve | c4_pta_product |
|---|---|
| 1% | 0.00648 |
| 5% | 0.03241 |
| 10% | 0.06482 |

Reading: a future angular table would need per-bin precision better than
roughly these fractions of the HD curve's swing to say anything about
`c4_pta_product` at that scale; anything cruder only re-confirms
"consistent with pure HD" without narrowing the parameter.

Command: `.venv-tda/bin/python audit/reverse_zero/e2_pta_and_tda_extension/e2_pta_sensitivity_scan.py`
Full numbers: `audit/reverse_zero/e2_pta_and_tda_extension/e2_pta_result.json`

## 5. Verdict / parameter_effect

**parameter_effect = none.** Free-parameter count stays at 2 (`mu_sym`,
`c4_pta_product`). No angular-correlation dataset exists to fit this
round, so `c4_pta_product` is neither constrained nor eliminated by
data. Setting `c4_pta_product = 0` inside the M0 zero-parameter
hypothesis is a **hypothesis change** -- the PTA sector set to its GR
value (pure Hellings–Downs) -- **tier L for GR itself**, not a
derivation from K3×T2 and not a data-driven reduction. This matches the
FRAMING RULE and every LeanMaster Stream 6-8 verdict reviewed this
round (E2/E3/E4): `mu_sym` and `c4_pta_product` are unchanged by all of
them; K3×T2 geometry (Golay octad / Kummer 8+16 split, moduli-trapping
at (ω,ω)) yields no PTA observable.
