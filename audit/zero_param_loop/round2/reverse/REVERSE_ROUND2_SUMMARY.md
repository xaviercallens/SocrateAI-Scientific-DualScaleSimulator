# Round-2 reverse loop summary (reduced -> TDA -> experiment -> hypothesis)

Commands, in order, each writing the file named:

```
python3 reduced_model.py                # selftest only, no output file
python3 sweep_reduced.py                # -> sweep_reduced.csv (225 rows, seed 43)
python3 tda_reduced.py                  # -> tda_reduced_report.json, tda_persistence_reduced_{2d,mu_only,c4_only}.json
python3 refit_reduced.py                # -> refit_reduced_report.json
```
(run with `.venv-tda/bin/python`, gudhi 3.13.0, numpy/scipy/pandas)

## 1. Reduced model

Both round-2 forward proposals applied via harness params
(`reduced_model.py`): **a_pot, b_pot deleted** -> analytic flat LCDM
(Omega_m = 0.31115, H0 = 67.6611 km/s/Mpc, from LeanMaster
`DarkEnergyScale.lean:84` / `SelfDualCutoff.lean:100`); **lambda_sym
deleted** -> held fixed at 1.0, not swept (forward round's exact-algebra
result, not re-derived here). Remaining free parameters: **mu_sym,
c4_pta_product**.

## 2. TDA (`tda_reduced.py`, seed 43)

All 225 swept points were numerically stable (mu_sym restricted to
[0.2, 3.0] to avoid the known BVP-divergence region; 129 joint + 48+48
single-parameter positive-control points).

| sweep | swept params | predicted Betti | observed Betti | predicted dim | observed dim (PCA var.) |
|---|---|---|---|---|---|
| 2d (joint) | mu_sym, c4_pta_product | [1,0] | **[1,0]** | 2 | **2** |
| mu_only (control) | mu_sym | [1,0] | **[1,0]** | 1 | **1** |
| c4_only (control) | c4_pta_product | [1,0] | **[1,0]** | 1 | **1** |

`topology_match = True` per the definition fixed in `PREREGISTRATION.md`.

**Disclosure (sequencing, not just results):** `PREREGISTRATION.md` said
"effective dimension" without disambiguating between two conventions
already in use elsewhere in this loop (the SVD-singular-value threshold
`jacobian.py` uses vs. the PCA-cumulative-variance-fraction k
`tda_gudhi.py` uses). Both were computed; **the PCA-variance metric was
selected as primary only after both had been run and the SVD-threshold
metric was seen to give the unstable (3,2,1) reading** (matching the
forward round's own already-reported instability of that exact metric).
Likewise, the log10 transform of `pta_max_deviation_from_hd` below was
adopted **after** the untransformed run produced a `Betti_0=9` artifact
for the `c4_only` control, not decided in advance. Both choices are
principled and documented, but neither was locked in before the data was
seen -- stated here rather than left implicit, since that is exactly what
pre-registration is meant to prevent.

Two findings worth recording honestly rather than smoothing over:

- A **stricter** dimensionality metric (singular value > 1e-3 of the
  largest, the convention `../jacobian.py` uses) gives (3, 2, 1) instead
  of (2, 1, 1) -- inflated by mild curvature of the mu_sym -> (log
  screening_suppression_factor, log phi_center_ratio) map and by residual
  correlation in the joint sweep. This is the SAME instability the
  forward round already reported for this exact metric ([2,3,2,2,3]
  across 5 points for a model of fixed true dimension); it is reported as
  a secondary, non-deciding number, not silently substituted for the
  primary PCA-variance-fraction metric.
- The `c4_only` control's observable (`pta_max_deviation_from_hd`) is
  EXACTLY linear in `c4_pta_product`, which is drawn log-uniformly over 8
  decades. Building the Rips complex on the raw (linear) values gives a
  heavy-tailed point cloud that fragments into 9 components
  (`Betti_0=9`) under a median-based edge-length threshold -- a pure
  sampling-density artifact, directly analogous to the c4_c0_ratio /
  pta_suppression sampling-density caveat the round-1 forward loop
  already flagged for the same reason. Log-transforming all three
  observable columns (matching the log-uniform sampling density) removes
  the artifact; this transform is applied uniformly to all three sweeps,
  not selectively to force a match.

**Dark-energy sub-vector alone** (Omega_m, H0, and every dark-energy
column in `sweep_reduced.csv`): confirmed constant across all 225 points
-- max pairwise distance = 0.0, exactly the pre-registered single-point
prediction.

**Bottleneck vs. forward round's real-data diagram**: reported in
`tda_reduced_report.json -> bottleneck_vs_forward_real_data_DESCRIPTIVE_ONLY`,
carrying the forward round's own caveat verbatim (different, independently
rescaled embeddings; not on equal footing). Per `PREREGISTRATION.md`, this
number does NOT set `topology_match`.

**Note on round-1's own bottleneck numbers** (task asked to check
`round1/reverse/reduced_tda.py` for "two bottleneck distances identical to
all digits"): re-derived directly with GUDHI in this session --
`model_vs_real_sn_H1` and `model_reduced_vs_model_round1_H1` are both
exactly `0.027069824244136675`, and this is **not a copy-paste error**.
Both equal exactly half the persistence of the round-1 reduced model's
single most-persistent H1 bar, `(0.30690863919317485 -
0.2527689907049015)/2 = 0.027069824244136675`: that bar is so much more
persistent than anything in either comparison target (the real-SN diagram
and round-1's own full-model diagram) that it is left unmatched (matched
to the diagonal) in BOTH comparisons, saturating the bottleneck cost at
the same value regardless of what the rest of either target diagram looks
like. The metric was saturated by one dominant bar and is therefore
UNINFORMATIVE about the actual difference between those two (very
different) target diagrams -- a caveat that applies to any bottleneck
number in this loop dominated by a single outlier bar, and is recorded
here in the round-2 ledger rather than left only in this session's
transcript.

## 3. Experiment (`refit_reduced.py`)

Independent re-implementation (own DESI-BAO/Pantheon+ loading and chi2,
not an import of `../decisive_experiment.py`):

| quantity | value | cross-check |
|---|---|---|
| chi2_reduced (frozen LCDM, 0 data-facing params) | **702.9066083231721** | vs forward's `decisive_experiment_report.json` 702.9066083231721 -- diff = 0.0 |
| chi2_lcdm (Om fitted, k=1) | **702.6753004830802** | vs forward's 702.6753004830802 -- diff = 0.0 |
| chi2_full (a_pot,b_pot best of 513-pt sweep, k=2) | 703.8946371287691 | REUSED from `../sweep_chi2.csv`, not recomputed here |

AIC (k_reduced=0, k_lcdm=1, k_full=2; mu_sym/lambda_sym/c4_pta_product
excluded from k for all three models -- no dataset fits them):

- `delta_aic_vs_lcdm` = AIC_reduced - AIC_lcdm = **-1.7686921599080279** (favors reduced)
- `delta_aic_vs_full` = AIC_reduced - AIC_full = **-4.988028805596969** (favors reduced)
- `delta_chi2_reduced_minus_full` = **-0.9880288055969686** (reduced model's chi2 is already LOWER than the 2-parameter quintessence sector's own best sweep point)
- `fit_degraded` = **False** under both the 2-removed-parameter threshold (4) and the conservative 3-removed threshold (6), since the delta is negative (an improvement, not a degradation).

**Load-bearing caveat**: `chi2_reduced` fits DESI+Pantheon+ with 0
parameters, not 2. `mu_sym` and `c4_pta_product` are free in the model but
enter no term that either dataset constrains (no fifth-force/screening or
angular-PTA dataset exists in `data/real/`; both ABSENT in
`data/real/MANIFEST.json`). Calling this a "2-parameter model that fits
the data" would be false.

## 4. Hypothesis (theory_feedback)

**Ladder count this round: 3 removed (a_pot, b_pot, lambda_sym), 2 remain
free (mu_sym, c4_pta_product); 5 -> 2.** The reduced model fits
DESI+Pantheon+ with **zero** data-facing parameters, not two -- and does
so at least as well as the 2-parameter quintessence sector's own best fit
(delta_chi2 = -0.99, delta_AIC = -4.99 vs. the full model, -1.77 vs.
fitted LCDM; both favor the reduced model). That is the single strongest
result of this round.

The dark-energy sector reaches zero data-facing free parameters this
round, but **by importing two Planck-2018 literature constants through
Lean** (`omegaLambda`, `hubbleRadius_m` -- tier L physics, kernel-fixed
tier as Lean `def`s), not by deriving them from dual-scale or moonshine
structure. The literal dual-scale hypothesis (l_micro ~ l_Planck paired
with l_macro ~ 1/H0) is independently FALSIFIED at cosmological scales by
~10^30 (`CKNInstance.lean`, `planckEnergy_gt_ckn_horizon_cutoff`), so this
round's "success" is best read as: *the quintessence potential this
harness implemented adds fitting freedom the fetched data does not
reward, and a 0-parameter flat-LCDM alternative (imported, not derived)
already does at least as well* -- not as evidence for any specific
dual-scale mechanism generating dark energy.

The two parameters that remain free after this round (`mu_sym`,
`c4_pta_product`) are algebraically real (nonzero Jacobian loadings on
their respective observable blocks, per the forward round) but
empirically invisible to every dataset this repo has fetched. The theory
cannot be pushed further toward zero parameters along this path without
either (a) fetching a real fifth-force/short-range-gravity dataset to
constrain `mu_sym`, or (b) fetching a real angular-correlation
Hellings-Downs product (not the frequency-only NANOGrav KDE already
verified) to constrain `c4_pta_product` -- or (c) accepting the S3/S4
"scope, not derivation" downgrade already flagged in
`../reduce_proposals.json` and reporting the observable vector as
correspondingly smaller, honestly, rather than as a further parameter
reduction.

LeanMaster's advance to 587 theorems this round is genuinely
zero-parameter mathematics (Hurwitz class numbers for immortal K3xT2
dyons, Goettsche numbers p24(k), chi(K3)=24 -- all passing independent
checks per the harness's supplied context) but touches NO observable in
this cosmology harness: it lives in a black-hole-microstate-counting
sector this loop has no bridge to. It reduces nothing here this round.
The "27720 lock" (77/60) downgrade -- now shown to fail the M24 twining
test at all 25 non-identity classes -- is correctly not used as evidence
anywhere in this round's proposals, consistent with round 1's decision to
already drop it.

**Net message for the next round**: the zero-free-parameter program is
advancing on two disconnected fronts (a Lean-verified, genuinely
zero-parameter black-hole-microstate sector, and a data-fitted,
literature-imported cosmology sector), and closing the gap between them
-- deriving `omegaLambda`/`H0`, or `mu_sym`/`c4_pta_product`, from the
same K3xT2/moonshine structure that fixes the dyon counts -- remains the
open problem, not a solved one.
