# X3: c4_pta_product bound from NANOGrav 15-yr data

Tier: X (exploratory numerics / data fit — a pipeline result, not a kernel-checked or literature
statement). This bounds `c4_pta_product` as defined in
`scripts/param_loop_sim.py:compute_pta_observable` (`gamma_theta = HD(theta) + c4_c0_ratio *
pta_suppression * l4_response`, so the fitted ratio `b_hat/a_hat` below IS `c4_pta_product`, the
product, not `c4_c0_ratio` alone). LeanMaster makes no statement about `c4_pta_product` (Stream 8:
"Observables: none"; verified again at LeanMaster's true current HEAD, see
`x3_leanmaster_check.json`). No claim here that K3xT2 predicts or fixes `c4_pta_product`.

**Note on the round-2 task's supplied "LEANMASTER AT HEAD" block**: it claims LeanMaster HEAD is
`eb791e7` (v3.25.0) with no commits since. That was stale at the time of this run — LeanMaster's
true current HEAD is `64f905f` (v3.28.0), 4 commits ahead of `eb791e7`. `x3_leanmaster_check.json`
re-runs the prohibited-term search at the true HEAD, using full-tree `git grep` rather than
`git show <commit> | grep` (which only searches that one commit's diff, not the whole tree, and so
is not by itself a reliable absence check). The conclusion (no mention of `mu_sym`,
`c4_pta_product`, "unit-bearing", or a TDA prediction) is unchanged, but independently reconfirmed
rather than re-quoted from the stale block.

## Environment

Separate venv, built exactly as instructed (did not touch `.venv-tda`):
`/mnt/disks/disk-socrateai-local-1/venv-pta` — `enterprise-pulsar==3.5.0`,
`enterprise_extensions==3.0.3`, `pint-pulsar==1.1.7`, `numpy==1.26.4`, `scipy==1.15.3`. Install
succeeded on the first attempt (`pip install "numpy<2" scipy enterprise-pulsar
enterprise_extensions pint-pulsar h5py`); no ABSENT/stop condition was hit.

## Data

`/mnt/disks/disk-socrateai-local-1/NANOGrav15yr/NANOGrav15yr_PulsarTiming_v2.1.0/narrowband`
(the registration's declared path
`data/real2/pulsar_timing/NANOGrav15yr_PulsarTiming_v2.1.0` is empty in this worktree; the real
tarball lives outside it). sha256 of the tarball:
`91476bf20d4c8baa9f5ad39c6e114b581478a5a7e269876489cafc8f62c5708f`, identical to
`/mnt/disks/disk-socrateai-local-1/zenodo_16051178/NANOGrav15yr_PulsarTiming_v2.1.0.tar.gz` (the
round-1 sha256-verified copy) — verified again here with
`sha256sum NANOGrav15yr_PulsarTiming_v2.1.0.tar.gz zenodo_16051178/....tar.gz`.

## Pulsar selection (data-blind: tim-file MJDs only)

`x3_select.py` -> `selection.json`. From the 76 narrowband entries: drop 6 per-telescope split
duplicates (`B1937+21ao/gbt`, `J1600-3053gbt`, `J1643-1224gbt`, `J1713+0747ao/gbt`,
`J1903+0327ao`, `J1909-3744gbt` — same physical pulsar as a combined entry already selected) and
1 short-span pulsar (`J0614-3329`, span < 3 yr) -> 67 selected, common `Tspan` = 16.03 yr. No
pulsar subset was named in the registration beyond "14 frequencies, gamma=13/3, fixed noise", so
this criterion is this script's own, stated explicitly.

Of the 67, **66 were used**: `J1713+0747` (59395 TOAs across 33 backend `.tim` files, the
largest/most heavily-observed pulsar in the release) failed to load as an `enterprise.Pulsar`
object in 3 foreground attempts of ~10 minutes each (prlimit 10 GB) — not a data-quality issue,
a computational-budget one. This is recorded as a deviation, not silently dropped.

## Per-pulsar model and noise (`x3_persr.py`, cached to `cache/*.npz`)

Per pulsar: EFAC + ecorr (basis) + t2equad per backend, intrinsic red noise (30 Fourier
components, own `Tspan`), and a **fixed** common process (CURN) at `gamma=13/3`, `log10_A=-14.62`,
14 frequencies, common `Tspan` across all pulsars (registration: "14 frequencies, gamma=13/3,
fixed noise"). **`log10_A=-14.62` is NOT attributed to the NANOGrav 15-yr paper** — fetching/
citing that paper's exact CURN amplitude was NOT ATTEMPTED this round; -14.62 is used only as a
representative fixed value. A full 66-pulsar sensitivity rerun at `log10_A=-14.0` (`x3_persr.py
--log10a -14.0 --outdir cache_alt_amp`, then `x3_combine.py --alt`) shows `c4_pta_product` is
**not** invariant to this choice — see "Sensitivity" below. White
and red noise parameters are fixed at the **median of the released MCMC noise chain**
(`noise/*.pars.txt` + `noise/*.chain_1.txt`, 25% burn-in dropped), not sampled — this is the "fixed
noise" the registration specifies. Full (non-marginalized) timing-model design matrix (SVD), as
`enterprise_extensions.OptimalStatistic` requires. The `X_a`, `Z_a` (Sigma-marginalized Fourier
residual/covariance) ingredients used by the optimal statistic are saved per pulsar; PINT
`Pulsar()` construction is cached separately (pickled) so re-running the analysis never reloads
the raw `.tim` files.

## Combination (`x3_combine.py` -> `x3_pta_result.json`)

- 2145 pairs among 66 pulsars; `min(xi) > 0` asserted (no duplicate-position pairs — the split
  duplicates removed in selection were exactly the risk here).
- 7 equal-pair-count angular bins (306-307 pairs each).
- GLS fit `rho_bin = a*HD(theta) + b*P4(theta)` (`HD`/`P4` match `hd_orf` /
  `compute_pta_observable`'s `l4_response` formulas exactly). `c4_pta_product = b_hat/a_hat`.
- **Deviation** from the registered "7x7 covariance that includes shared-pulsar correlations":
  `enterprise_extensions==3.0.3` has no analytic pair-covariance function (`grep -rn "pair_cov"
  .../frequentist/*.py` -> no match). Substituted a delete-one-pulsar jackknife covariance (66
  jackknife samples), which removes all pairs sharing a pulsar together and so captures
  shared-pulsar correlation empirically.
- Sky-scramble gate (1000 scrambles, seeds `8000000+k`): pulsar sky positions randomly permuted,
  `rho_ij`/`sig_ij` kept attached to their real pulsar-index pair. **Two candidate statistics were
  computed and both are reported** (not one chosen after seeing which passed): the raw scrambled
  bin vector against `C`, dof=7 (`rho_bin_scrambled^T C^-1 rho_bin_scrambled`, mean chi2/7 =
  **0.931**, passes [0.8,1.2]), and the post-GLS-fit residual, dof=5 (mean chi2/5 = **0.872**,
  which is *not* comparable to the [0.8,1.2] band since fitting 2 of 7 parameters away already
  predicts ~5/7=0.71 on the dof=7 denominator for a perfectly calibrated `C`). The registered
  wording "chi2/7" is read as naming the raw dof=7 statistic, which is the one used for the gate
  decision. **Deviation**: the same frozen jackknife covariance `C` is reused for every scramble (a
  fresh 66-pulsar jackknife per scramble x 1000 was infeasible in the stated time budget); the gate
  therefore checks calibration against the real-data covariance, not a scramble-specific one, and —
  because `C` is the jackknife covariance of the same `rho_ij` being rebinned — this is closer to an
  internal-consistency check of the jackknife scale than an independent calibration test.
- Monopole/dipole negative controls are fit **jointly with HD+P4** (3-parameter GLS), because a
  single-regressor fit is confounded by the real HD signal (monopole/dipole are not orthogonal to
  HD over 7 bins); the unmarginalized single-regressor numbers are reported too, for comparison
  only.
- Verified: a hand-written 3-pulsar re-implementation check
  (`x3_verify_os.py` -> `x3_verify_os_result.json`) reproduces
  `enterprise_extensions.OptimalStatistic.compute_os`'s `xi`, `rho`, `sig` to machine precision
  (0.0 max relative difference), and confirms no `"<param> not set!"` warnings were silently
  swallowed by earlier output filtering.

## Result

- `c4_pta_product = -0.187` (**conditional on `log10_A_CURN=-14.62`** — see Sensitivity below)
- 95% Fieller interval: **[-0.841, 0.101]** — contains 0. The HD amplitude itself is only
  `a_hat_snr_sigma = 2.42 sigma` from zero, so this bound is marginal: a modestly weaker HD
  detection here would flip the registered rule to "no bound" instead of an interval.
- 68% Fieller interval: [-0.351, -0.063] — does not contain 0.
- Monopole negative control (marginalized over HD+P4): amplitude z = 1.05 sigma (not significant).
- Dipole negative control (marginalized over HD+P4): amplitude z = 0.60 sigma (not significant).
  (Unmarginalized single-regressor z: monopole 1.11 sigma, dipole 1.83 sigma — confounded by HD,
  reported for comparison only, see above.)

Per the registered decision rule: **M0 (c4 = 0) is consistent with the data at 95%** (0 lies
inside the 95% interval), but this is a marginal result (2.4 sigma HD detection with 66/67
pulsars, one high-precision pulsar excluded — see Sensitivity and the J1713+0747 deviation). The
68% interval excluding 0 is reported as-is and is not, by the registered rule, a rejection — no
claim of "confirms" or "rejects" is made. This is a data fit (Tier X); nothing here is a K3xT2
prediction, and LeanMaster makes no statement about `c4_pta_product` (re-verified at LeanMaster's
true current HEAD — see below).

## Sensitivity (log10_A_CURN) — result is NOT amplitude-independent

A full 66-pulsar rerun at `log10_A=-14.0` instead of `-14.62` (`x3_persr.py --log10a -14.0
--outdir cache_alt_amp` then `x3_combine.py --alt` -> `x3_pta_result_alt_amp.json`) gives:

| | log10_A=-14.62 (base) | log10_A=-14.0 (alt) |
|---|---|---|
| `c4_pta_product` | -0.187 | -0.324 |
| 95% Fieller interval | [-0.841, 0.101] (bounded) | **unbounded** |
| sky-scramble gate (chi2/7) | 0.931 (pass) | 0.424 (**fail**) |

`c4_pta_product` is **not invariant** to the assumed fixed CURN amplitude — contrary to a naive
expectation that `c4=b_hat/a_hat` would cancel a common rescaling (the amplitude enters the
per-pulsar Sigma-marginalization, not just an overall normalization of `rho`). The base-run bound
should be read as conditional on `log10_A_CURN=-14.62`, which was itself **NOT ATTEMPTED** to
verify against the NANOGrav 15-yr paper this round.

## Commands (from repo root, `dualscale-wt-reverse`)

```
cd audit/reverse_zero_r2/X3-pta && /mnt/disks/disk-socrateai-local-1/venv-pta/bin/python x3_select.py
cd audit/reverse_zero_r2/X3-pta && /mnt/disks/disk-socrateai-local-1/venv-pta/bin/python x3_persr.py --workers 6
cd audit/reverse_zero_r2/X3-pta && /mnt/disks/disk-socrateai-local-1/venv-pta/bin/python x3_combine.py
cd audit/reverse_zero_r2/X3-pta && /mnt/disks/disk-socrateai-local-1/venv-pta/bin/python x3_verify_os.py
cd audit/reverse_zero_r2/X3-pta && /mnt/disks/disk-socrateai-local-1/venv-pta/bin/python x3_persr.py --workers 6 --log10a -14.0 --outdir cache_alt_amp
cd audit/reverse_zero_r2/X3-pta && /mnt/disks/disk-socrateai-local-1/venv-pta/bin/python x3_combine.py --alt
```
All run as `prlimit --as=10737418240 -- <python> <script> <args>` in the foreground, per ground
rule 6. Seeds: sky scrambles `8000000+k` (k<1000, `numpy.random.default_rng`); no other random
draws are used (noise parameters are chain medians, not sampled).
