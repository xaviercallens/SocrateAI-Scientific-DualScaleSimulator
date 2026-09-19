> **Correction notice (2026-09-19):** statements in this file were corrected after the skeptic review. See `audit/reverse_zero/ERRATA.md` and REPORT.md §6. Wherever they differ, REPORT.md is authoritative.

# E3 -- M0 (zero-parameter hypothesis) vs M2 vs fitted LCDM, and the LeanMaster P2/kappa=1 note

Script: `audit/reverse_zero/E3-M0/m0_model.py`
Command: `/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python audit/reverse_zero/E3-M0/m0_model.py` (run from the worktree root)
Output: `audit/reverse_zero/E3-M0/m0_result.json`
Tier: **X** (data fit, no new Lean statement). Seeds: none needed (deterministic closed-form/least-squares chi2, same convention as round2/round3).

**Data-match disclosure:** the computed task asked to run "on the SAME data as E1-DR1." `audit/reverse_zero/E1-desi-dr2/` exists in this worktree but is **empty** (no output was produced there this loop). This script therefore could not match E1's data and instead uses the newest verified datasets in the manifest: DESI **DR2** ALL_GCcomb + the Pantheon+SH0ES **full STAT+SYS covariance** (Sec. 2). If E1 is later run on different data, its numbers must not be silently compared against this round's without checking they used the same dataset.

## 1. M0 defined, and its parameter count verified

M0 = flat LCDM with Omega_Lambda = 0.68885 **frozen** (LeanMaster `DarkEnergyScale.lean:84`), GR tensor sector (`c4_pta_product = 0`), symmetron sector **removed** (`mu_sym` has no value -- the sector is deleted, not set to a number). H0 (equivalently the BAO amplitude nuisance `u_star` and the SN additive offset) is a fitted **nuisance**, not a theory parameter, by the same convention rounds 2-3 already used and disclosed (marginalizing it analytically leaves it unable to distinguish any Omega_m, so "freezing H0" changes nothing in this fit).

Every number, listed with source (`m0_result.json` -> `parameters.M0`):

| Number | Value | Source |
|---|---|---|
| Omega_Lambda | 0.68885 | LeanMaster `DarkEnergyScale.lean:84` (Planck 2018; tier L physics) |
| Omega_m | 0.31115 | = 1 - Omega_Lambda |
| c4_pta_product | 0 | GR value (Hellings-Downs exactly) |
| mu_sym | -- (no value) | sector removed |
| H0 / u_star / SN offset | fitted per dataset | nuisance, marginalized in closed form |

**M0 has 0 free theory parameters.**

## 2. Data (recomputed, not reused, since both datasets changed from round 2/3)

- DESI **DR2** ALL_GCcomb consensus BAO (13 points, 13x13 covariance) -- `data/real2/dark_energy/desi_dr2/desi_gaussian_bao_ALL_GCcomb_{mean,cov}.txt`, fetched_verified.
- Pantheon+SH0ES **full STAT+SYS 1701x1701 covariance**, sliced (data vector and covariance with the same boolean mask, so row/column correspondence is exact) to a cosmology-only cut `zHD>0.01 & zHD<=2.4 & IS_CALIBRATOR==0` -> **1580 SNe** (advisor-review fix: the raw 1701-row file includes 77 Cepheid-host calibrator SNe whose `MU_SH0ES` is tied to `CEPH_DIST`, a distance-ladder anchor, not a Hubble-flow distance modulus -- fitting them with the same shape+offset cosmological model is a likelihood error, not a data upgrade). Rounds 2-3 used the same z-cut without the `IS_CALIBRATOR` filter and diagonal errors only, landing on 1590 SNe. This is the "full Pantheon+ covariance" pipeline action `PRE_REGISTRATION.md`'s P1 calls for.
- **Primary N_total = 1593** (13 BAO + 1580 SN). A disclosed, non-primary variant using the full 1701-row sample (including calibrators) is also reported for comparison (`variant_full_1701_including_calibrators` in `m0_result.json`); the headline verdict is unchanged between the two (Sec. 3).

## 3. Results (primary: 1580-SN cosmology-only cut)

| Model | k (theory params) | chi2 | AIC | BIC |
|---|---|---|---|---|
| **M0** (frozen Omega_m, no symmetron, c4=0) | 0 | 1401.071 | 1401.071 | 1401.071 |
| M2 (M0 + mu_sym, c4 reinstated free) | 2 | 1401.071 (identical) | 1405.071 | 1415.818 |
| LCDM fitted (Omega_m free) | 1 | 1400.325 | 1402.325 | 1407.699 |

- Delta chi2 (M0 - LCDM fitted) = **0.746** (1 dof). Threshold table (`PRE_REGISTRATION.md`:35; this comparison was NOT pre-registered): "1 dof: Delta chi2 = 9.00 is 3sigma". **NOT REJECTED, and not close** -- M0's frozen Omega_m sits well inside the data's own preferred value (fitted Omega_m = 0.30422, close to the frozen 0.31115). Read this as "consistent at 0.86 sigma equivalent (chi2.sf(0.746, 1) = 0.388)", not primarily via AIC/BIC (see caveat below).
- Delta AIC (M0 - LCDM fitted) = **-1.25**, Delta BIC = **-6.63**: on parsimony grounds alone M0 costs nothing in fit quality and drops a parameter. **Caveat (do not over-read):** a negative Delta BIC here is a restatement that the frozen value happens to land close to the data's fitted optimum -- it would flip sign (favor the fitted model) had the frozen Omega_Lambda missed by only a little more. The Delta chi2 = 0.746 (1 dof) figure is the primary, load-bearing number; AIC/BIC are subordinate color, not independent confirmation.
- Delta AIC (M2 - M0) = **+4.00**, Delta BIC (M2 - M0) = **+14.75** (= 2 ln N, purely definitional -- see note): reinstating mu_sym and c4 as free parameters buys **zero** improvement in chi2 (they enter no term in this data's likelihood -- confirmed identical to 15 decimal places) and is penalized on both criteria for that reason alone. **This reproduces round-3's `chi2_independent.json` "probes_note" finding on the new DR2 + full-covariance data: chi2 is flat across the whole (mu_sym, c4) sweep because no dataset in hand tests either parameter.**
- Negative controls (same z-cut, 1580 SNe): Omega_Lambda = 0.5 gives a large chi2 penalty vs fitted; 0.72 and 0.65 give smaller but still clearly non-zero penalties (`m0_result.json` -> `negative_controls`) -- the frozen value is not merely "not rejected", it sits close to the data's own minimum.
- **Disclosed full-1701 variant** (including calibrators, N_total=1714): M0 chi2 = 1772.870, fitted LCDM Omega_m = 0.31036, chi2 = 1772.861, Delta chi2 = 0.0094. The qualitative verdict (M0 not rejected, close to the data's own optimum) is stable across both variants; the exact Delta chi2 differs because the calibrator SNe change what the SN chi2 minimum looks like.

## 4. Screening: does M0's "no fifth force" prediction have a bound?

- **Fifth-force data status: ABSENT.** The Eot-Wash tarball (`arXiv:2002.11761`, `data/real2/fifth_force/2002.11761.tar.gz`, fetched and sha256-verified) contains only the LaTeX source and 9 figure PDFs -- no machine-readable alpha-lambda table (grep for `tabular`/`table` environments returns nothing). The exclusion curve exists only as a plotted figure and was not digitised, per ground rules.
- **Unit bridge status: NO BRIDGE, cannot map.** `mu_sym` in `workshopcosmo.run_symmetron_screening_simulation` / `scripts/param_loop_sim.py` is a dimensionless BVP parameter with no documented conversion to a physical mass, length or coupling scale. Even with a machine-readable Eot-Wash table in hand, there is today no formula that turns an alpha-lambda exclusion point into a bound on `mu_sym` in the harness's own units.
- **Does any data prefer the symmetron sector?** No. Confirmed again on this round's DR2 + full-covariance fit (Sec. 3, M2 vs M0 identical chi2): no dataset in this session's manifest constrains `mu_sym` at all. Its removal in M0 is free.

## 5. LeanMaster P2 (self-dual length) and the [21,105] um bracket

Quoted verbatim from LEANMASTER_CONSTRAINTS:

> **P6.2 (kappa=1):** "Within the programme's own T-duality R -> alpha'/R with alpha' = s^2, the self-dual radius is the unique positive fixed point R = s (`p62_selfdual_fixed_point`) ... With the standard KK normalisation m_n = n/R, this gives kappa = 1: P1 is the programme's prediction, not one convention among several."
>
> **P1 (excluded):** "P1 is excluded. R exceeds each of the three bounds. By the pre-registered decision rule (T1 fails), P1 is excluded. The failure is not a rounding effect: R > 1.5 x 30 um and R > 1.06 x 44 um."

What this means for `PRE_REGISTRATION.md`'s registered **[21.0, 105.3] um** bracket (that file is **not edited** here, per ground rules):

- The bracket exists because of an **independent, unrelated** O(1) ambiguity: the reduced-vs-non-reduced Planck-mass convention, factor (8*pi)^(1/4) = 2.239. LeanMaster's kappa=1 result does not touch or shrink that bracket -- kappa=1 is about T-duality's OWN normalisation (`m_n = n/R` vs some other KK convention), a *different* O(1) choice from the Planck-mass one.
- What kappa=1 DOES remove is the escape route "maybe kappa is a free tuning knob you could dial down within T-duality to save P1". It is not: self-duality forces it to exactly 1, so P1's central value (47.008 um) is *the* T-duality prediction, not a convention-dependent guess -- and LeanMaster separately proves that exact value fails all three cited bounds.
- The bracket's lower half (per `PRE_REGISTRATION.md`'s own language, "only the lower part of the bracket survives") still formally survives as a Planck-mass-convention artifact, but T-duality's own preferred point (kappa=1, 47.008 um) sits **outside** that surviving lower half. There is no rescue for P1 inside the programme's own T-duality.

## 6. Framing-rule compliance

M0 is a hypothesis change (Omega_Lambda imported from Planck, c4 set to the GR value, symmetron sector deleted), not a K3xT2 derivation. Nothing in this script or report claims K3xT2 predicts or fixes mu_sym, c4, or Omega_Lambda.
