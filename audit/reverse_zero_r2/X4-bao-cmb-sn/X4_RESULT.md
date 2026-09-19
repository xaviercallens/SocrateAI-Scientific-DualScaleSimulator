# X4 -- DESI DR2 BAO + compressed CMB distance priors + Pantheon+ full covariance

Registration: `audit/reverse_zero_r2/registration/registration.json`, test `X4`, frozen at commit
`9e6705e940f35aef8b377ec3631f62ba3d52aea4` (verified an ancestor of HEAD; the committed file's
sha256 equals the blob at that commit -- see `registration.commit_is_ancestor_of_HEAD` and
`registration.file_sha256_equals_commit_blob` in `x4_results.json`).

**Status for PRE_REGISTRATION P1 (the cosmological constant): INFORMATIONAL.** PRE_REGISTRATION.md's
falsification rule for P1 names DESI DR3 or final BAO + CMB + one SN compilation, Euclid, or
Rubin/LSST supernovae (quoted verbatim and grep-asserted in `x4_results.json.quotes_verbatim`).
DESI **DR2** BAO + compressed Planck-2018 distance priors + Pantheon+ is none of those, so **no P1
verdict follows from X4**, whatever the number. This is the same reading round 2's E1/A1 already
established (ERRATA.md, PRE_REGISTRATION.md addendum A1) -- X4 only adds a CMB likelihood on top.
LeanMaster Stream 6's own P1 (extra-dimension radius) is a **different** P1 and is untouched by X4.

## Scripts (unified, ONE SN cut)

- `x4_fetch_prior.py` -- fetches arXiv:1808.05724 (abstract + e-print), extracts the TeX source,
  records sha256 in `data/fetch_log.json`.
- `x4_core.py` -- unified likelihood: DESI DR2 ALL_GCcomb BAO (13x13 cov), the CMB distance-prior
  chi2, and Pantheon+ with the full 1701x1701 STAT+SYS covariance, sliced by **one** boolean mask
  function (`sn_cut_mask`) for both the 1580-SN primary cut (`zHD in (0.01,2.4]`, `IS_CALIBRATOR==0`)
  and the 1590-SN disclosed variant (same range, no `IS_CALIBRATOR` filter). This replaces and
  unifies `audit/reverse_zero/E1-desi-dr2/e1_desi_dr2.py` and
  `audit/reverse_zero/E3-M0/m0_model.py`, whose separate, undocumented SN cuts were flagged in
  `ERRATA.md`'s last row ("1590 vs 1580 SNe... Different cuts").
- `x4_camb_check.py` -- validates a CAMB `r_drag(omega_b, omega_m)` interpolation table (40 seeded
  points, seed 9100000, max relative error 1.4e-5) and cross-checks the published CMB prior table
  against CAMB 2.0.4 at the parameters the table itself implies.
- `x4_bao_cmb_sn.py` -- `fit` (20 seeded minimiser starts, seed `9000000+i`, per model), `polish`
  (CAMB-exact CMB, a sensitivity variant), `regress` (no-CMB regression mode, reproduces the E1/E3
  convention with `u`/offset profiled), `assemble` (writes `x4_results.json`, does all verbatim-quote
  and provenance assertions).

Commands actually run (`<venv>` = `/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python`,
all under `prlimit --as=10737418240 --`, foreground, from this directory):
```
<venv> x4_fetch_prior.py
<venv> x4_camb_check.py
<venv> x4_bao_cmb_sn.py regress --cut 1580
<venv> x4_bao_cmb_sn.py regress --cut 1590
<venv> x4_bao_cmb_sn.py fit --cut {1580,1590} --model {m0,lcdm,cpl} --start-lo I --start-hi J   (covering starts 0..19, run in several chunks to respect timeouts)
<venv> x4_bao_cmb_sn.py polish --cut {1580,1590} --model {m0,lcdm,cpl}
<venv> x4_bao_cmb_sn.py assemble
```

## 1. CMB distance-prior source (quoted, checked)

**Chen, Lu; Huang, Qing-Guo; Wang, Ke, "Distance Priors from Planck Final Release," arXiv:1808.05724
(JCAP 02 (2019) 028, doi:10.1088/1475-7516/2019/02/028).** Fetched from arXiv (abstract page +
e-print), sha256-recorded in `data/fetch_log.json`.

Base flat-LCDM (Planck TT,TE,EE+lowE) row of the paper's Table (`\label{distance priors}`, TeX
lines 140-142), used as the CMB data vector and (with its correlation matrix) covariance:

| Quantity | Value | TeX source line |
|---|---|---|
| R | 1.7502 +/- 0.0046 | 140 |
| l_A | 301.471 (+0.089/-0.090) | 141 |
| omega_b h^2 | 0.02236 +/- 0.00015 | 142 |

The **inverse** covariance (before normalization) used directly in the chi2, copied verbatim from
the paper's own CosmoMC snippet (`\verb|~/cosmomc/data/Distance_invcov.txt|`, TeX lines 303-306):
```
     94392.3971   -1360.4913   1664517.2916
    -1360.4913   161.4349   3671.6180
     1664517.2916   3671.6180   79719182.5162
```
and the matching central values from `distance.ini` (TeX lines 289-291): `r=1.750235`,
`la=301.4707`, `omegabh2=0.02235976`. All of these are grep-asserted verbatim against the fetched
TeX in `x4_results.json.paper.quotes_verbatim` (every entry `matched: true`).

**Covariance cross-check:** inverting the published inverse-covariance gives sigma_R=0.00462,
sigma_lA=0.0889, sigma_wb=0.000149 and correlations (0.463, -0.660, -0.327), matching the table's
0.0046, 0.0895, 0.00015 and (0.46, -0.66, -0.33) to within rounding
(`x4_camb_check.json.covariance_check`).

**Recomputation check against CAMB 2.0.4:** solving the paper's own R(z*), l_A(z*) formulas for
(Omega_m, h) at omega_b=0.02235976 that reproduce the table's (R, l_A) gives Omega_m=0.3163,
h=0.6744 (omega_m=0.1438). Running CAMB at that point gives a theta*-based l_A offset from the
table by -0.61 sigma_lA and R by +0.08 sigma_R (`x4_camb_check.json.table_implied_point`); most of
this is the well-known difference between the paper's Hu-Sugiyama z* fit and CAMB's exact
recombination z*, not an error in the table (z*_HS=1092.0 vs CAMB z*=1090.0). **A one-line code/prose
discrepancy inside the paper itself is disclosed**: the g_1 coefficient is 0.0738 in the prose (eq.
after "Hu:1995en", TeX line 200) but 0.0783 in the paper's own Fortran likelihood code (TeX line
435); this script uses the code value (0.0783), the one that actually generates the table.

## 2. Regression against E1/E3 (no CMB; sanity check on the unification)

`regress` mode (`Omega_r=0`, `u`/SN-offset profiled, no CMB) reproduces round-2's numbers closely
(`x4_results.json.regression_no_cmb`):

| Quantity | This script (1580 / 1590) | E1/E3 JSON (round 2) |
|---|---|---|
| M0 vs fitted Lambda, Delta chi2 (1 dof) | 0.746 / 0.788 | E3: 0.746 (1580) |
| CPL vs Lambda, Delta chi2 (2 dof), E1-style unbounded start | 4.735 (1580) / 4.547 (1590) | E1: 4.547 (1590); ERRATA: "4.547 and 4.735" |

The match confirms the unified likelihood reduces to the same numbers as the two separate round-2
scripts once CMB is turned off, given the same SN cut.

## 3. Primary result: BAO + CMB distance priors + SN, Omega_m and h free (paper's own R, l_A formulas)

20 seeded Nelder-Mead starts per model (seed `9000000+i`). For M0 and Lambda **all 20/20** starts
land within 1e-3 of the reported minimum. For CPL, **19/20** converge to the reported minimum;
start index 3 in both cuts lands on a distinct, higher local minimum (chi2=1526.7 at 1580,
1544.3 at 1590, vs. 1397.3/1413.9 for the other 19) and is excluded as non-converged
(`x4_results.json.fits.<cut>.models.cpl.non_converged_start_indices` = `[3]`). No parameter sits at
a search-bound at the reported best fit for any model (`at_bound` is empty in all six cases).

| Cut | Model | chi2_min | Omega_m | h | omega_b | w0, wa | chi2 (bao / sn / cmb) |
|---|---|---|---|---|---|---|---|
| 1580 | M0 (frozen Omega_L=0.68885) | 1408.955 | 0.31106 | 0.6781 | 0.02239 | -1, 0 | 20.05 / 1388.39 / 0.51 |
| 1580 | Lambda (Omega_m free) | 1404.280 | 0.30277 | 0.6843 | 0.02251 | -1, 0 | 11.92 / 1389.69 / 2.67 |
| 1580 | CPL (Omega_m, w0, wa free) | 1397.281 | 0.31148 | 0.6771 | 0.02243 | -0.858, -0.510 | 9.32 / 1387.27 / 0.69 |
| 1590 | M0 | 1425.438 | 0.31106 | 0.6781 | 0.02239 | -1, 0 | 20.05 / 1404.88 / 0.51 |
| 1590 | Lambda | 1420.712 | 0.30273 | 0.6843 | 0.02251 | -1, 0 | 11.90 / 1406.13 / 2.68 |
| 1590 | CPL | 1413.944 | 0.31119 | 0.6774 | 0.02244 | -0.861, -0.502 | 9.36 / 1403.89 / 0.70 |

M0's Omega_m = 0.31106, not the literal 1-0.68885 = 0.31115: `om_m0(h)` solves Omega_m
self-consistently against the same Omega_r(Omega_m, h) radiation prescription used for every other
model, at flat Omega_Lambda = 0.68885 fixed; the 0.00009 shift is the radiation density at that h.

The chi2_cmb contribution (3 data points) stays O(0.5-2.7) across all six fits -- it is not driving
the comparisons through some large, unaccounted internal tension; the CAMB-exact sensitivity check
below gives the same qualitative answer, confirming this is genuine (small) added information from
CMB, not a definitional artifact.

**CPL vs Lambda (registered primary, N_eff=1, uncorrected), 2 dof:**
- 1580: Delta chi2 = 6.999, p (Wilks) = 0.0302, **sigma = 2.167** (scipy `chi2.sf` then
  `norm.isf(p/2)`). Below both the 3-sigma (11.83) and 5-sigma (28.74) thresholds. Fitted point
  (w0=-0.858, wa=-0.510) **is** in the w0>-1, wa<0 quadrant that the falsification rule names, but
  the magnitude is far short of the 5-sigma bar.
- 1590: Delta chi2 = 6.768, sigma = 2.121. Same qualitative reading.

**M0 vs fitted Lambda (reported alongside, uncorrected, not a rejection test), 1 dof:**
- 1580: Delta chi2 = 4.675, p = 0.0306, **sigma = 2.162**. `p < 0.05`, so by the registered wording
  this is **not** "not rejected" -- it reads as "disfavoured relative to fitted Lambda at that DR2 +
  compressed-CMB combination" (registration's own phrasing for this case), still **not a P1
  verdict**.
- 1590: Delta chi2 = 4.726, sigma = 2.174. Same reading.

This is a genuine change from round 2's CMB-free result (E3: Delta chi2 = 0.746, "not rejected,
0.86 sigma"): adding even three compressed CMB numbers moves M0 from consistent to mildly
disfavoured, because R and l_A break the Omega_m-h degeneracy that BAO+SN alone leaves open. It is
still nowhere near a 3-sigma (9.00) or 5-sigma (25.00) threshold, and neither PRE_REGISTRATION nor
this script's registration named this M0-vs-Lambda comparison as a rejection rule.

**Sensitivity variant (CAMB-exact R, l_A, r_drag in place of the paper's own closed-form
formulas):** re-optimizing from the paper-formula best fit with CAMB 2.0.4 supplying R, l_A and
r_drag directly gives essentially the same numbers (`x4_results.json.fits.<cut>.sensitivity_camb_exact`):
CPL-vs-Lambda sigma = 2.153 (1580) / 2.107 (1590); M0-vs-Lambda sigma = 2.247 (1580) / 2.259 (1590).
The qualitative reading is unchanged between the paper's own distance-prior formulas and a direct
CAMB evaluation.

## 4. Is the published ~3.1 sigma approached?

**No.** PRE_REGISTRATION.md R4 records "Disfavored at 3.1sigma" for DESI DR2 BAO + CMB + SNe, quoting
LeanMaster's `DualScaleValidation/Observables.lean` module docstring, itself citing arXiv:2503.14738
and arXiv:2503.14743 (both grep-asserted verbatim in `x4_results.json.quotes_verbatim`; LeanMaster
Stream 8 elsewhere records "Observables: none" for its own K3xT2 predictions, also grep-asserted).
This script's CPL-vs-Lambda sigma is **2.12-2.17** (2.11-2.15 with the CAMB-exact sensitivity
variant) -- about 1 sigma short of 3.1, and using only 3 compressed CMB numbers (R, l_A, omega_b)
rather than the full Planck power-spectrum likelihood plus DESI's actual DR2 pipeline. The shortfall
is recorded, not adjusted toward 3.1 (`x4_results.json.fits.<cut>.published_3.1sigma_comparison`).
This is expected: a compressed 3-parameter CMB summary carries less constraining power than the full
likelihood, and this script's own CPL model differs in nuisance handling (e.g., the disclosed r_d*h
profiled variant vs. this run's non-profiled r_d(omega_b,omega_m)) from whatever exact combination
the published 3.1-sigma result used.

## 5. Framing-rule and ground-rule compliance

- M0 (Omega_Lambda=0.68885, GR tensor sector, no symmetron sector) is stated as a **hypothesis
  change, not a K3xT2 derivation**, throughout the script's docstrings and this report. Nothing here
  claims K3xT2 predicts or fixes mu_sym, c4, Omega_Lambda or any TDA outcome.
- Every quote from PRE_REGISTRATION.md and LeanMaster is asserted with a fixed-string `grep -F`
  command recorded in `x4_results.json` (`quotes_verbatim`, all `matched: true`); qualifiers are kept
  (e.g. "as recorded in LeanMaster `DualScaleValidation/Observables.lean` module docstring").
- Two different P1s: this section only discusses PRE_REGISTRATION's P1 (the cosmological constant);
  LeanMaster Stream 6's P1 (extra-dimension radius) is untouched by X4 and not conflated with it.
- DR2 results are informational / a pipeline action, never a P1 verdict (Sec. "Status for P1" above).
- Wording: "not rejected", "disfavoured relative to fitted Lambda", "consistent with" are used; not
  "confirms", "cleared decisively" or "the data reward". All sigma values are `scipy.stats.chi2.sf`
  then `scipy.stats.norm.isf(p/2)`, never eyeballed.
- Registration: `audit/reverse_zero_r2/registration/registration.json` test `X4` was committed
  (9e6705e) before any data file here was opened by this script; `x4_results.json.registration`
  records the file-hash-equals-commit-blob and ancestor-of-HEAD checks that confirm this.

## 6. Deviations from the registered rule (disclosed)

- The registered rule asks for "r_d computed from (omega_b, omega_m) by one stated fitting formula,
  NOT profiled when the CMB priors are included": this script computes r_drag with a
  CAMB-2.0.4-built interpolation table (validated to 1.4e-5 relative error against direct CAMB
  calls, `x4_camb_check.json.rd_spline_vs_camb`) rather than a hand-coded literature fitting
  formula (e.g. Eisenstein-Hu 1998), because CAMB is available in this environment and is more
  accurate than any closed-form fit; this is a strengthening, not a weakening, of the registered
  method, and is disclosed here rather than edited into `registration.json`.
- z* uses the paper's own Hu-Sugiyama fit with its **code** coefficient (0.0783), not its **prose**
  coefficient (0.0738); both are quoted verbatim in `x4_results.json.paper.verbatim_lines` so the
  discrepancy in the source paper itself is visible, and the code value is used because it is the
  one that actually reproduces the paper's own published table (Sec. 1).
- Radiation density Omega_r uses the paper's own prescription (Omega_r = Omega_m/(1+z_eq),
  z_eq=2.5e4 Omega_m h^2 (T/2.7K)^-4) rather than a standard photon+3-neutrino calculation, matching
  what the registered `cmb_observables` field describes ("massless neutrinos (approximation
  stated)").
- **The CPL minimiser rejects w0+wa >= 0 with a smooth penalty** (a DESI-style prior carried over
  from round 2's `e1_desi_dr2.py` convention), which the registration's own CPL spec ("free:
  Omega_m, h, omega_b, w0, wa") does not impose. **Effect: inactive at every reported best fit** --
  w0+wa = -1.368 (1580) and -1.363 (1590), both well inside the unconstrained region -- so none of
  the reported chi2/sigma numbers are changed by it. Recorded in
  `x4_results.json.registration.deviations`.
- The sample-size assert added to `x4_core.py`'s `load_sn` (`assert int(m.sum()) == cut`) was added
  to the source file after this run's `fit`/`polish` calls had already completed; it is a no-op
  guard (both cuts already gave exactly 1580/1590 by construction of `sn_cut_mask`) added for future
  reruns, and `assemble` does not reload the SN data, so it changes nothing reported here.
