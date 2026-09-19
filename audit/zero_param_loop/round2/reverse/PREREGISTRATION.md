# Round-2 reverse loop: pre-registration

Committed BEFORE `sweep_reduced.py`, `tda_reduced.py` or `refit_reduced.py`
are run, so the "state the prediction before looking at data" claim in the
task is checkable against git history, not asserted after the fact.

## Reduced model

`reduced_model.py` applies both round-2 forward proposals
(`../reduce_proposals.json`): `a_pot,b_pot` deleted (frozen flat LCDM,
Omega_m = 0.31115, H0 = 67.6611 km/s/Mpc, both from LeanMaster
DarkEnergyScale.lean:84 / SelfDualCutoff.lean:100); `lambda_sym` deleted
(held fixed at 1.0, not swept, per the forward round's exact-algebra
result -- not re-derived here). Remaining swept parameters: `mu_sym`,
`c4_pta_product`.

## Predictions (stated now, checked in `tda_reduced.py` / `refit_reduced.py`)

1. **Dark-energy sub-vector alone**, evaluated at any number of distinct
   (mu_sym, c4_pta_product) points: every row is byte-identical (it is a
   pure function of Z_GRID/Z_POINTS and two frozen Lean constants, with no
   dependence on either swept parameter). Predicted persistence diagram:
   a single point in observable space -> all pairwise distances 0 ->
   Betti_0 = 1, Betti_1 = 0, and the Vietoris-Rips filtration is trivial
   at every scale > 0.

2. **Full reduced observable vector** (dark-energy + screening + PTA
   columns together), swept over BOTH mu_sym and c4_pta_product within a
   numerically-stable range: the dark-energy columns will have exactly
   zero variance and must be dropped before embedding (z-scoring a
   constant column divides by zero); the effective dimension of the
   remaining (screening, PTA) columns is predicted to be **2**, and the
   cloud is predicted to be **connected with no 1-dimensional hole**
   (Betti_0 = 1, Betti_1 = 0), because it is the smooth image of a
   connected, simply-connected 2D rectangle in (log mu_sym,
   log c4_pta_product) under a continuous map, PROVIDED the swept range
   avoids the mu_sym divergence region (mu_sym >= ~2-3 sends the symmetron
   BVP fallback to NaN per round-1/round-2 forward findings) -- points
   that fail `numerically_stable` are excluded, and their count is
   reported, not silently imputed.

3. **Positive controls** (pre-registered so "matches the prediction" is
   falsifiable, not just declared): (a) sweep `c4_pta_product` alone with
   `mu_sym` fixed at nominal (1.0) -> predicted effective dimension 1,
   Betti [1, 0]; (b) sweep `mu_sym` alone (within the stable range) with
   `c4_pta_product` fixed at nominal (0.08035) -> predicted effective
   dimension 1, Betti [1, 0]. If either control fails to show dimension
   1, the "dimension = number of free parameters" reading of the 2D sweep
   is not trustworthy either.

4. **Bottleneck comparison with the forward round's real-data diagram**
   (`../tda_persistence_real_3d_sn_comoving.json`): reported as a
   descriptive number only, carrying the forward round's own caveat
   verbatim (model cloud and real cloud are built from different
   embeddings -- PCA of standardized model observables vs a native 3D
   comoving construction -- each independently rescaled, so raw Betti_0
   counts and bottleneck magnitudes are not on equal footing). This
   number does NOT set `topology_match`; `topology_match` is defined
   below as computed-vs-predicted, not model-vs-real-universe.

## `topology_match` definition (fixed now)

`topology_match` = True iff the OBSERVED Betti numbers of the full
2-parameter reduced sweep (item 2) equal the PREDICTED Betti numbers
[1, 0] stated in item 2, AND both positive controls (item 3) show
effective dimension 1. It does not depend on the bottleneck distance to
the real-data cloud (item 4), which is reported separately.

## chi2 / AIC bookkeeping (fixed now)

- `chi2_reduced`: the frozen-flat-LCDM chi2 against DESI 2024 BAO +
  Pantheon+ SH0ES, RECOMPUTED independently in `refit_reduced.py` (own
  script, not a re-import of `../decisive_experiment.py`), then compared
  against the forward round's `decisive_experiment_report.json ->
  lcdm_frozen_lean.chi2_total` (702.9066083231721) as a cross-check.
  `mu_sym` and `c4_pta_product` contribute 0 to this number: no
  fifth-force or angular-PTA dataset exists in `data/real/` to fit them
  against (both ABSENT in `data/real/MANIFEST.json`). This model therefore
  fits DESI+Pantheon+ with **0** parameters, not 2; that distinction is
  carried through to `theory_feedback` explicitly.
- `chi2_lcdm`: LCDM with Omega_m FITTED (1 free parameter), also
  recomputed independently in `refit_reduced.py`, cross-checked against
  `decisive_experiment_report.json -> lcdm_om_fitted.chi2_total`
  (702.6753004830802).
- `chi2_full`: the full (pre-round-2) model's best point, taken from the
  forward round's own 513-point Sobol sweep over (a_pot, b_pot)
  (`../sweep_chi2.csv`, best chi2_total = 703.8946371287691) -- reused,
  not recomputed, because reproducing the full quintessence ODE optimizer
  is out of scope for a reverse-loop cross-check; this reuse is stated
  explicitly, not disguised as an independent recomputation.
- AIC parameter counts: k=0 for the reduced model (nothing is fit to
  data), k=1 for fitted-Omega_m LCDM, k=2 for the full quintessence sector
  (a_pot, b_pot only -- mu_sym, lambda_sym, c4_pta_product are not fit to
  any of these datasets in either model and are excluded from k for all
  three, so the AIC comparison isolates the dark-energy-sector deletion).
  Sign convention: `delta_aic_vs_X = AIC_reduced - AIC_X`; negative favors
  the reduced model.
- `fit_degraded` = True iff `chi2_reduced - chi2_full > 2 * n_removed`,
  with `n_removed` = 2 (a_pot, b_pot; the only two chi2-relevant deletions
  this round -- lambda_sym and the mu_sym/c4_pta_product retention do not
  change chi2_full or chi2_reduced against DESI+Pantheon+ at all, since
  neither model fits mu_sym/c4_pta_product to that data either). Both the
  2-removed-parameter threshold (4) and the more conservative
  3-removed-parameter threshold (6, counting lambda_sym too) are reported;
  the sign of `chi2_reduced - chi2_full` is expected to be negative (the
  forward round already found frozen LCDM chi2 slightly BETTER than the
  quintessence sector's own best sweep point), which would make
  `fit_degraded` False under either threshold.

## Known caveats carried forward unchanged

- `vacuum_decay_cdl`, `tadpole_cancellation_certified`, local
  `proofs/BuscherRules.lean`: not used as evidence anywhere below.
- `c4_pta_product` and `mu_sym` remain free but untestable by any dataset
  fetched into this repo; their "removal" is not attempted this round
  (S3/S4 from the round-1 staged plan are still open).
- The literal dual-scale hypothesis (l_micro ~ l_Planck paired with
  l_macro ~ 1/H0) is FALSIFIED at cosmological scales by ~10^30
  (CKNInstance.lean, `planckEnergy_gt_ckn_horizon_cutoff`); the
  omegaLambda/hubbleRadius_m constants used above are imported as Planck
  2018 literature values (tier L), not derived from dual-scale physics.
  Identifying them as THIS model's dark-energy sector remains tier C.
