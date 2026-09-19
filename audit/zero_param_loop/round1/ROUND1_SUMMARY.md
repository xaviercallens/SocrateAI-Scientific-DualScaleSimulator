# Zero-param loop -- Round 1 (FORWARD: hypothesis -> experiment -> TDA -> reduce)

Worktree: `/mnt/disks/disk-socrateai-local-1/dualscale-wt-loop` (branch `loop/zero-param`).
Python: `/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python` (gudhi 3.13.0, seed 42 everywhere).
All files below are under `audit/zero_param_loop/round1/`.

## 0. Starting point

6 free parameters: `a_pot, b_pot, mu_sym, lambda_sym, pta_suppression, c4_c0_ratio`.
Harness: `scripts/param_loop_sim.py` (`evaluate_point`), unchanged this round.

## 1. EXPERIMENT: wide-range sweep

Command: `python param_sweep.py` (writes `sweep.csv`).
Method: Latin hypercube (`scipy.stats.qmc.LatinHypercube`, seed 42), 420 points + the nominal
point (421 total), each parameter spanning 4 decades in log10 (a_pot,mu_sym,lambda_sym:
[1e-2,1e2]; b_pot,pta_suppression: [1e-4,1e0]; c4_c0_ratio: [1e-1,1e3]). Evaluated with
`multiprocessing.Pool(6)` (`_worker_init` overrides `param_loop_sim`'s `Z_GRID` to
`linspace(0,2.45,80)` and `Z_POINTS` to DESI's 7 z's, so no `np.interp` boundary-clamping
happens for DESI z=2.330). Per-point `signal.alarm(45)` timeout guard.

Result: **421/421 points completed (0 errors/timeouts)**; **105/421 (25%) numerically_unstable**
(the `mu_sym`-divergence zone where the symmetron relaxation fallback silently returns
non-finite `screening_suppression_factor`, per round-0's Check-B finding -- caught here, not
misread as "unchanged"). 316 numerically-stable points used downstream.

### Data used (real, verified)

- `data/real/dark_energy/desi_2024_bao_all.txt` (12 rows, verified round 0) **+ its covariance**,
  fetched THIS round from the same trusted repo (`CobayaSampler/bao_data`,
  `desi_2024_gaussian_bao_ALL_GCcomb_cov.txt`), sha256 `bbafa907...` -- see
  `data_provenance.json`. Row/column block structure verified by direct inspection
  (`python -c "...nonzero pattern..."`) to match the 12-row mean file exactly:
  index 0 (z=0.295, DV alone), 1-2 (z=0.510, DM+DH), 3-4 (z=0.706), 5-6 (z=0.930),
  7-8 (z=1.317), 9 (z=1.491, DV alone), 10-11 (z=2.330).
- `data/real/dark_energy/pantheon_plus_sh0es.dat`: Hubble-flow subset (`0.01 < zHD <= 2.4`,
  finite `MU_SH0ES`/`MU_SH0ES_ERR_DIAG`), **N=1590 SNe**, diagonal errors only (the
  1701x1701 stat+sys covariance was not fetched -- ground rule "diagonal if covariance too
  large").
- **PTA: ABSENT.** `unzip -l data/real/pulsar_timing/nanograv_kde_freespectrum.zip` and a
  `find ... -iname '*angsep*' -o -iname '*pair*' -o -iname '*corr*'` over the unzipped tree
  show only `density.npy/freqs.npy/log10rhogrid.npy/bandwidths.npy` per analysis folder -- a
  per-frequency free-spectrum KDE, **not** an inter-pulsar angular-separation dataset. No real
  `Gamma(theta)` likelihood exists in the verified data to test `pta.gamma_theta` against.
  (A raw `NANOGrav15yr_PulsarTiming_v2.1.0` .tim/.par tarball was found on disk at
  `/mnt/disks/disk-socrateai-local-1/NANOGrav15yr_PulsarTiming_v2.1.0/` during this round, but
  turning that into an angular correlation requires running a full timing pipeline -- out of
  scope this round; flagged for round 2.)

### chi2 (`chi2_and_jacobian.py` -> `chi2_report.json`, `sweep_chi2.csv`)

One nuisance per sector, both handled exactly (analytic linear least squares, not a bounded
numerical search -- a first version used a bounded 1D search over `ln(s)` that silently
clipped at its boundary for the LCDM baseline because the model's and LCDM's native distance
scales differ by orders of magnitude; caught by comparing against the analytic solution and
fixed): BAO nuisance `s = r_s*H0/c` (data ~ base/s, quadratic in `u=1/s`, solved exactly);
SN nuisance additive offset (absorbs `5log10(c/H0) - M_B`, solved exactly as a weighted mean
residual).

**Bug found and fixed in-round**: the harness's `distance_modulus_shape_grid` (`=5*log10(D_L*H0)`)
diverges like `log(z)` near z=0; linearly interpolating that *already-logged* quantity on the
80-point grid at Pantheon+'s many `z~0.01-0.03` SNe gave `mu_model(z=0.01) ~ -43` instead of the
correct `~-7`, inflating `chi2_sn` to ~1.6e6. Fixed by storing the smooth
`D_M_times_H0_grid` instead and computing `5*log10((1+z)*D_M_H0(z))` **after** interpolating
the smooth (not logged) quantity at the exact query z's (`param_sweep.py`,
`chi2_and_jacobian.py::mu_shape_from_DM_H0`). Re-ran the sweep once (`sweep.csv` regenerated,
same 421/0/105 counts) after this fix.

| | chi2_bao (dof) | chi2_sn (dof) | chi2_total (dof) | chi2/dof |
|---|---|---|---|---|
| **Model best-fit** (idx=282, a_pot=3.471, b_pot=4.28e-4) | 16.672 (9) | 687.210 (1587) | 703.882 (1596) | 0.4409 |
| **Flat LCDM baseline** (Om=0.2939 fit) | 12.741 (10) | 684.600 (1588) | 697.340 (1598) | 0.4364 |

**Headline finding (as required: report honestly if the model is LCDM-like everywhere)**:
**it is NOT.** `w0_cpl_latetime` ranges **[-1.0737, -0.0776]** and `wa_cpl_latetime` ranges
**[-0.0006, 1.1003]** over the whole stable sweep (`chi2_report.json.headline_finding`) --
nowhere close to uniformly pinned at LCDM's (-1, 0). The model spans a wide range of
dark-energy behavior across its (a_pot,b_pot) sweep, most of which is a poor-to-mediocre fit,
but the sweep's own chi2-minimum reaches **chi2/dof=0.441, statistically indistinguishable
from LCDM's 0.436** (Delta-chi2 ~ 6.5 for 2 fewer dof) -- at THAT specific point, despite
w0=-0.694, wa=0.955 being far from LCDM's values. This is an important, non-obvious result:
the model is not "secretly LCDM" (the w0/wa summary statistic varies hugely), but there exists
a non-LCDM corner of (a_pot,b_pot) space that fits BAO+SN about as well as LCDM does.

## 2. RESPONSE RANK (`jacobian_response_rank.py` -> `jacobian_report.json`)

Central differences in `ln(param)`, step 1e-2 (not 1e-6: the ODE/BVP solves have finite
tolerance). Decade-spanning observables (`screening_suppression_factor`, `phi_center_ratio`)
log10-transformed before standardization. 11-component observable vector (see file docstring).
Only numerically-stable points used; 316 of them recomputed (parallel, 6 workers) to get
sweep mean/std for standardization. Probed at the chi2-best-fit point (idx=282) and 5 random
stable points (seed 42).

**Two machine-precision-exact null directions at ALL 6 probed points, independent of where in
parameter space the point sits:**

1. `ln(pta_suppression) - ln(c4_c0_ratio)` (loadings exactly `+-0.7071`/`-+0.7071`):
   singular value `0.0` to `3.5e-18` at every point. This is the **exact algebraic
   degeneracy** already found in round 0's Check C (max diff 1.110e-16) -- now confirmed
   independently by a completely different method (SVD of a numerical Jacobian) at 6
   different points in the 6-dimensional space. **This is this round's tier-B reduction: 6 -> 5.**
2. `ln(lambda_sym)` alone (loading exactly `1.0`): singular value `3.2e-14` down to `1.1e-19`
   at every point. Stronger than round 0's borderline 1.2e-6-relative finding, but **still
   insensitivity, not a derivation** (ground rule: "insensitivity alone is NOT derivation").
   Kept at **0.5-evidence / convention**, in the count.

**a_pot, b_pot are NOT universally degenerate** -- unlike the two directions above, their
near-null singular values vary by point and are never machine-precision zero (1e-3 to 1e-14,
point-dependent): `effective_dimension_at_1e-3_of_max` is **3 at the best-fit point**
(singular values `[4.65, 2.58, 1.72, 3.7e-4, 3.2e-14, 0.0]`) but drops as low as **1** at one
random point (idx=99, where BOTH a_pot's and b_pot's own directions become locally
unconstrained, `1.0e-5` and `9.4e-14`) and rises to **4** at another (idx=46). This is a
genuine nonlinear, point-dependent weak-identification effect in the dark-energy sector, not
an algebraic degeneracy -- reported as the honest, non-uniform picture the data show, per the
advisor's explicit instruction to report the full spectrum rather than a single threshold count.

`mu_sym` never appears in a null direction at any probed point (own resolvable singular
direction, e.g. SV=1.7176 at best-fit) -- confirms round 0's Check D independently.

## 3. TDA (real GUDHI, `import gudhi`; `tda_gudhi.py` -> `tda_report.json` + diagrams)

**Interpretation caveat (stated up front, per advisor guidance)**: Betti numbers of a smooth
image of a parameter box are NOT a parameter count -- a contractible blob gives beta0=1,
beta1=beta2=0 regardless of the box's true dimension. TDA here is used for cloud SHAPE and for
the mandated controls; the dimension count is the Jacobian/PCA result above.

- **Model side**: standardized 20-dim observable vector (w0,wa; DM*H0 and H/H0 at DESI's 7
  z's; log-ssf, log-pcr; pta max-dev, log-c4pta) over the 316 stable points, PCA'd to
  **4 dims at 97.3% cumulative variance** (`pca_cumulative_variance_fraction`). **Correction
  (post-advisor review): this is NOT an independent cross-check of the Jacobian's eff_dim.**
  The Jacobian is a *local* derivative rank at one point (how many directions move the
  observables there); the PCA is a *global* variance decomposition of 316 points spread over a
  4-decade box (how many linear axes carry the cloud's spread). These measure different things
  and their landing at similar numbers (4 vs. 3 at best-fit) is not evidence either is right --
  a strongly nonlinear 2-parameter model can produce a high-rank PCA, and a near-degenerate
  6-parameter model confined to a thin curved sheet can produce a low-rank PCA at high local
  Jacobian rank. Reported here only as its own finding (dark energy + screening observables
  compress to a 4-dim linear subspace over this sweep box), not as confirmation of anything
  else. **The effective-dimension number this round is the Jacobian's: 3, at the best-fit
  point, per Section 2.** Rips complex (rescaled to max pairwise distance=1, `sparse=0.2`):
  betti_numbers `[1,0]` (one connected, contractible blob, as expected).
- **Real side (primary)**: Pantheon+ SN residual-vs-best-fit-flat-LCDM embedding
  `(zHD, c, x1, residual)`, 4D, N=316 (subsampled to match model cloud size, seed 42). Also
  betti `[1,0]`.
- **Real side (extra, qualitative)**: 2MRS RA/Dec angular patch (N=316 of 50000, embedded on
  the unit sphere `x,y,z`). **NOT** the comoving-coordinate analysis the task's first TDA
  option names -- `data/real/cosmic_web/2mrs_sample.tsv` was verified (again) to carry only
  `2MASS-ID, RAJ2000, DEJ2000, Kmag`, no redshift/velocity column, so no comoving 3D embedding
  is available from the verified set without inventing a Kmag-to-distance relation (forbidden).
- **Controls**: noisy circle (unit radius, sigma=0.05, N=316) -- **assertion PASSED**: one
  dominant H1 bar, top persistence 0.644 vs second-highest 0.0102, ratio **62.9**
  (`circle_control_assertion.ASSERTION_one_dominant_H1_bar: true`). Poisson/uniform null cloud,
  same N=316 and (post-rescale) dimension as the model cloud (4, after PCA).
- **Bottleneck distances** (all clouds rescaled to max pairwise distance=1 first; still
  labeled qualitative): `model_vs_real_sn_H1 = 0.0195`, `null_vs_real_sn_H1 = 0.0399`; same
  pattern against 2MRS (`model_vs_real_2mrs_H1 = 0.0291` vs `null_vs_real_2mrs_H1 = 0.0403`).
  **Caveat added post-advisor review**: the null cloud is drawn uniformly in a box; this
  ordering (model closer to real than null is) is reported **with a uniform-box null only,
  NOT tested against a null matched to the model cloud's own per-axis (PCA) variances** -- if
  the ordering is actually a bounding-box artifact of the uniform null being "more spread out"
  than either the model or real sheet, a variance-matched null could flip it. Not claiming the
  model resembles the real topology more than chance until that control is run (round 2).

## 4. REDUCE -- proposals (`reduce_proposals.json`)

**This round's result: 6 -> 5 (tier B), not 6 -> 4.** One clean, doubly-confirmed
`absorb_combination`:

- `pta_suppression`, `c4_c0_ratio` -> single parameter `c4_pta_product` (nominal 0.08035).
  Tier B. Evidence: round-0 Check C (1.11e-16) + round-1 Jacobian (exact null direction, all
  6 probed points).

Kept at 0.5-evidence (NOT counted as removed, per ground rule):

- `lambda_sym` -> convention (fixed at 1.0), NOT derived. Same two evidence sources, but
  insensitivity is explicitly not derivation.

Two Tier-C candidates identified but explicitly NOT adopted this round (both have caveats
from their own source and neither is supported by this round's data):

- `(a_pot, b_pot)` -> LeanMaster `dark_energy_w0=-1, dark_energy_wa=0`
  (`DualScaleValidation/Observables.lean:97-98`, verified present in LeanMaster this round).
  The sweep's own chi2-best-fit point lands at (w0,wa)=(-0.694, 0.955), far from (-1,0), and
  fits the data statistically indistinguishably from a point that DOES sit near (-1,0) --
  the identification is not preferred by the data, and the Lean constant's own docstring flags
  tension with DESI DR2 at 3.1sigma.
- `c4_pta_product` -> LeanMaster `planckGmuBound=1.5e-7`
  (`DualScaleCosmology/CosmicString.lean:83`, verified present). Off by ~5.4e5x from the
  nominal value, and the constant's own caveat says it applies to cosmic-string tension, not
  GW suppression; no real PTA angular dataset exists this round to test it either way.

## Files written this round

- `param_sweep.py`, `sweep.csv` (421 rows)
- `chi2_and_jacobian.py`, `chi2_report.json`, `sweep_chi2.csv`
- `jacobian_response_rank.py`, `jacobian_report.json`
- `tda_gudhi.py`, `tda_report.json`, `tda_persistence_{model,real_sn,real_2mrs,null,circle}.{json,png}`
- `reduce_proposals.json`
- `data_provenance.json` (this round's one new fetch: DESI BAO covariance)
- this file

## What did NOT happen this round (explicit, per ground rules)

- No use of `vacuum_decay_cdl`, `tadpole_cancellation_certified`, or `proofs/BuscherRules.lean`
  as evidence anywhere above.
- No edit to the shared `data/real/MANIFEST.json` (another session also writes to this
  worktree -- commit `6f56bf4`); the new covariance fetch is recorded separately in
  `data_provenance.json`.
- No touch to `proofs/`, no push, no Zenodo deposit, no edit outside this worktree (LeanMaster
  was read-only, via `grep`/`find` only).
- Reverse loop (5-param model re-swept/re-TDA'd/re-fit) is NOT done this round -- listed as
  round-2's first priority in `reduce_proposals.json.next_round_priorities`.
