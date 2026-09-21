# Pre-registered falsification criteria

**Written 2026-09-18, before the future data named below.** The commit hash of this file is
the registration. A public timestamp exists only once it is pushed to origin.

Sources:
- simulator branch `loop/zero-param`: commits 8b8911c (lambda_sym deleted) and 1d6639f (round-3 skeptic);
- LeanMaster at `v3.17.0-2-gf4da238`, read-only.

Column "checked": **S** = verified by a command in this session (file/line or computation quoted);
**M** = from memory, not checked in this session (instrument dates and sensitivities).

## 0. What this covers, and what it does not

| Part | Status | Falsifiable by experiment? |
|---|---|---|
| K3×T² mathematics: χ(K3) = 24, Göttsche numbers p₂₄(k), lattice signatures (3,19), (4,20), (6,22), dyon counts via class numbers | Kernel-checked in LeanMaster (tier A). 8 statements were independently re-derived by a blind Python/GUDHI route (`loop/k3t2-rigidity`, comparator 9feff72; skeptic audit still running). | **No.** These are theorems. Their physical meaning (e.g. dyon microstate counting) is tier L and concerns black-hole physics that no planned instrument probes directly. |
| Cosmology sectors of the simulator | ΛCDM with imported constants, plus two parameters no dataset constrains (`mu_sym`, `c4_pta_product`) | **Partly.** Rows below. |

A falsification of any row below hits the cosmological postulate named in that row.
It does not hit the K3×T² mathematics: none of these numbers is derived from K3×T².

## 1. Tests already run (reported, not pre-registered)

| # | Claim | Data | Result | Checked |
|---|---|---|---|---|
| R1 | H₀ = c / `hubbleRadius_m` = 67.661 km/s/Mpc (`SelfDualCutoff.lean:100`), imposed on the SH0ES-calibrated supernova magnitudes (SN offset frozen at the value implied by H₀ = 67.66; left free, the SNe prefer 73.4) | DESI 2024 DR1 BAO (12×12 cov) + Pantheon+SH0ES (1590 SNe, diagonal errors) | **Rejected: Δχ² = +1204.5** (1907.4 vs 702.9; `loop/zero-param` commit 1d6639f, `audit/zero_param_loop/round3/skeptic/chi2_independent.json`). This is the Hubble tension. H₀ cannot be frozen to this constant against SH0ES-calibrated supernovae. | S |
| R2 | Flat ΛCDM with Ω_Λ = 0.68885 frozen (`DarkEnergyScale.lean:84`) vs Ω_Λ fitted | same | **Not rejected**: Δχ² = 0.231 (1 dof). The frozen value lies inside the pipeline's 1σ interval Ω_m ∈ [0.305, 0.329]. | S |
| R3 | w = −1 exactly (frozen Λ) vs CPL (w₀, wₐ) | same, our pipeline | CPL preferred by Δχ² = 6.09 (2 dof) = **2.0σ**. Uses diagonal SN errors only (χ²_SN/N = 0.43), so this is weaker than a full-covariance analysis. | S |
| R4 | w = −1 exactly | DESI DR2 BAO + CMB + SNe (arXiv:2503.14738, 2503.14743), as recorded in LeanMaster `DualScaleValidation/Observables.lean` module docstring | **Disfavored at 3.1σ** (w₀ > −1, wₐ < 0 quadrant preferred). Published result, **not re-analysed** by this pipeline. | S (that LeanMaster records it); M (the analysis itself) |

## 2. Pre-registered predictions and rejection rules

Thresholds (computed with `scipy.stats.chi2.isf`):
- 1 dof: Δχ² = 9.00 is 3σ, and 25.00 is 5σ;
- 2 dof: Δχ² = 11.83 is 3σ, and 28.74 is 5σ.

### P1. Dark energy is a cosmological constant (w₀ = −1, wₐ = 0), with Ω_Λ = 0.68885
- **Current status.** Already disfavored at 3.1σ by published DESI DR2 combinations (R4). Our own DR1-based pipeline gives 2.0σ (R3). This is not an open question waiting for data.
- **Falsified if** any of these, analysed with full published covariances, prefers the w₀ > −1, wₐ < 0 quadrant over w = −1 at **≥ 5σ (Δχ² ≥ 28.74, 2 dof)**:
  - DESI DR3 or final BAO combined with CMB and any one major SN compilation;
  - Euclid;
  - Rubin/LSST supernovae.
- The 5σ bar is the conventional discovery level. It is **not** a retreat from the 3.1σ already observed, which is recorded above as "disfavored".
- **Also falsified if** a ΛCDM fit to those data puts Ω_Λ away from 0.68885 by **≥ 5σ** of the measurement error.
- **Pipeline action, to be run with these thresholds fixed now:**
  - re-run `decisive_experiment.py` with DESI DR2 BAO and the full Pantheon+ covariance;
  - report the CPL-vs-Λ Δχ² as it comes out.
- Instruments and dates (M): DESI DR3 about 2026–2027; Euclid cosmology releases from about 2026; Rubin/LSST supernovae from about 2027–2030.

### P2. Self-dual length s = √(ℓ_P · c/H₀) ("dark dimension" reading, LeanMaster Stream 3)
- **Value.** s = 47.008 µm (S), from ℓ_P = 1.616255×10⁻³⁵ m and c/H₀ = 1.3672×10²⁶ m (`SelfDualCutoff.lean:96,100`; `selfDual_length_sq_bracket`).
- **The O(1) factor is not derived.** The reduced vs non-reduced Planck-mass convention alone contributes (8π)^{1/4} = 2.239 (S). So the registered prediction is the **bracket [21.0, 105.3] µm**, not 47 µm.
- **What is already known.**
  - LeanMaster `selfDual_above_torsion_radius_bound` proves s ≥ 1.5 × 30 µm against the Eöt-Wash 30 µm toroidal-radius bound. The gravitational-strength Yukawa range bound is < 38.6 µm (Lee et al., arXiv:2002.11761, as pinned in LeanMaster).
  - The central value 47 µm is therefore already above the bound. Only the lower part of the bracket survives.
  - LeanMaster `selfDual_lambda_above_mvv_range`: λ_s = 0.5355 (S) is at least 5× above MVV's central range.
- **Falsified if** a short-range gravity experiment excludes a gravitational-strength (α = 1) Yukawa deviation for every range ≥ **21.0 µm** at 95% CL.
- Instruments (M): Eöt-Wash-type torsion balances and other sub-mm gravity experiments.

### P3. Parameters that are bounded, not predicted (cannot be falsified as they stand)
| Parameter | What data would do | Instruments (M) |
|---|---|---|
| `mu_sym` (symmetron mass scale; λ_sym was deleted, 8b8911c) | bound it | lunar laser ranging, MICROSCOPE-type equivalence-principle tests, atom interferometry, binary pulsars |
| `c4_pta_product` (ℓ = 4 term in the PTA angular correlation) | bound it; a measured correlation consistent with pure Hellings–Downs would exclude large values | NANOGrav, IPTA, then SKA |
| `planckGmuBound` = 1.5×10⁻⁷ (`CosmicString.lean:83`) | nothing: it is an imported upper bound, not a prediction | LISA, PTAs |

These become falsifiable only when the theory predicts their values.

## 3. Numbers in the corpus that are NOT registered, and why

| Number | Source | Reason it is not registered |
|---|---|---|
| r = 1/252 ≈ 0.00397 | `Observables.lean`, `tensor_to_scalar_ratio_reduction` (proves only 110·252 = 27720·1) | It is derived from 27720, the "BPS lock". LeanMaster `ratio_fails_at_every_class` shows that lock fails the M₂₄ twining test at all 25 non-identity classes, so the formula has no group-theoretic basis. LeanMaster's own docstring calls its window "illustrative" and says the master theorem is "not a 'zero free parameters' theorem". **Cannot be pre-registered until its provenance is replaced by a derivation.** If one is found, the intended rule would be: falsified if r < 0.002 at 95% CL (LiteBIRD / CMB-S4, σ(r) ~ 10⁻³, M). |
| n_s = 0.965 | `Observables.lean`, `spectral_index_scaled` | This is Planck's measured value used as a benchmark, not a prediction. |
| δ_CP = 283° | `Observables.lean`, `neutrino_cp_phase_deg` | The window is marked "illustrative, not a NuFIT confidence interval". It is not derived in this repository. |

## 4. What would make the theory falsifiable in a stronger sense

All registered numbers are postulates or imported constants. None is derived from the K3×T² mathematics that LeanMaster verifies. A zero-parameter theory whose failure would matter needs **one** of these to be derived before the data arrive:
- the O(1) factor in the dark-energy scale (P2), which would turn the [21, 105] µm bracket into a number;
- a replacement for the 27720-based formula for r.

Until then, a failure of P1 or P2 refutes a cosmological postulate. It does not refute K3×T².

## Addendum 2026-09-19 (appended; original rules above unchanged)
Source: simulator branch `loop/reverse-zero` (E1 commit 4eba061, E3 commit 0c590a5, skeptic commits 240e421/1d42c16/e5cf853, ledger 5b6302e, report 13756dc, errata b77d11a); LeanMaster read-only at `eb791e7` (v3.25.0).

**A1. P1 pipeline action (lines 46-48) executed.** DESI DR2 BAO (ALL_GCcomb, 13×13 cov) + Pantheon+SH0ES full STAT+SYS covariance, 1590 SNe, no CMB, script `audit/reverse_zero/E1-desi-dr2/e1_desi_dr2.py` (not `decisive_experiment.py`). CPL vs Λ, Ω_m free in both: Δχ² = 4.547 (2 dof, 1.63σ), w₀ = −0.896, wₐ = −0.186. ΛCDM fit Ω_Λ = 0.6960 (0.89σ from 0.68885). The falsification rules at lines 40-45 name DESI DR3/final BAO + CMB + SN, Euclid, or Rubin; they are **not applicable** to this dataset, and no verdict on P1 is recorded from it. The published DESI DR2 + CMB + SN 3.1σ (R4) is not reproduced (no CMB likelihood used).

**A2. R3 superseded in accuracy (R3 kept as the historical record).** Same DR1 BAO and same 1590-SN cut: diagonal SN errors give Δχ² = 6.091 (reproduces R3), full STAT+SYS covariance gives Δχ² = 3.191 (1.27σ). R3's 2.0σ is an artefact of neglecting SN covariance.

**A3. Reported, not pre-registered: M0.** M0 = frozen flat ΛCDM with Ω_Λ = 0.68885 (imported Planck18 value), symmetron sector removed, c4_pta_product = 0 (GR). It is a hypothesis change, not a derivation from K3×T². It carries two profiled nuisances (r_d·h, SN offset); H0 is not frozen (R1). DESI DR2 + Pantheon+ full cov (1580 SNe): Δχ²(M0 − fitted ΛCDM) = 0.746 (1 dof, 0.86σ).

**A4. P2.** No new short-range gravity data analysed; the rule at line 58 is not triggered. LeanMaster Stream 6 (STREAM6_EXPERIMENT_PLAN.md:74-83 at `eb791e7`) excludes R = s ≈ 47 µm as an extra-dimension radius and gives κ = 1 under the Tier C identification α′ = s²; this does not remove the (8π)^{1/4} convention spread of line 53.

**A5. No TDA threshold is registered in this file.** The TDA thresholds used on 2026-09-19 (ξ(r) gate RMS z < 3, LOO 95th percentile, Bonferroni over 6) were fixed in script source, not here. Any future TDA or defect test enters this file, with statistic, null, seeds and decision rule, before its data are loaded.

## Addendum 2026-09-19b (appended; A1-A5 and the original rules above unchanged)
Source: branch `loop/reverse-zero`, round 2 (registration `audit/reverse_zero_r2/registration/registration.json` frozen at commit 9e6705e **before any data was loaded**; report `bc29b3f`); LeanMaster read-only at `eb791e7` (v3.25.0); TDA validation on branches `loop/tda-validation`, `loop/tda-simple`.

**A6. Round-2 registered tests X1-X4.**
- **X1 (CMB TDA, spectrum-matched null).** The spectrum gate passed this time (max |z| = 0.85 against a 2.81 threshold; round 1's null was mismatched and produced no usable p-values). WMAP ILC, primary: the Gaussian isotropic null is **not rejected**, p_corr = 0.472 (Sidak, N_eff = 2). Planck PR3 SMICA, conditional and not pooled: **anomaly at p_corr = 0.004**, causes not separable (chance, foregrounds, mask edge, component-separation residuals; none isolated). Neither result is a statement about K3 x T2.
- **X2 (cosmic web, redshift space).** Gate power demonstrated (rejection rates 0.995 and 0.985, required >= 0.95), then calibration **failed**: data chi2 = 40.21 against a 20.18 threshold, p = 0.004 (2.88 sigma). **INCONCLUSIVE** by the registered rule. This is a statement about the lognormal + linear-bias + finger-of-god mock family, not about M0.
- **X3 (NANOGrav 15-yr, `c4_pta_product`).** The first analysis of real pulsar timing data in this programme (66 pulsars, 2145 pairs, 7 angular bins). **No unconditional bound.** At one assumed common-noise amplitude (log10 A = -14.62) the Fieller interval is [-0.841, 0.101], which contains the GR value 0; at log10 A = -14.0 the calibration gate fails and no bound follows. The P3 row for `c4_pta_product` is therefore **unchanged**: still bounded-in-principle, not bounded in fact.
- **X4 (DESI DR2 BAO + compressed CMB + Pantheon+).** **INFORMATIONAL for P1**, since DR2 is not one of the datasets P1's rule names. CPL versus Lambda: 2.17 sigma (1580 SNe) and 2.12 sigma (1590 SNe), about 1 sigma short of the published 3.1 sigma of R4 because the CMB enters only through three distance priors. M0 against fitted Lambda, reported alongside: **2.16 sigma**, up from 0.86 sigma (A3) once the CMB priors are added.

**A7. Verdicts recorded on 2026-09-19.**
- **P1 (dark energy is a cosmological constant with Omega_Lambda = 0.68885): NOT FALSIFIED, and not confirmed.** No registered falsification rule has fired: line 40-45 names DESI DR3 or final BAO with CMB and SN, Euclid, or Rubin, and requires 5 sigma. What exists today: 1.63 sigma on DR2 BAO + SN (A1), 2.17 sigma once compressed CMB priors are added (A6), against a published 3.1 sigma that this pipeline does not reproduce. **Omega_Lambda = 0.68885 remains an imported Planck 2018 value (tier L), not a derived one**, so P1 is a constrained phenomenological fit, not an ab-initio prediction.
- **P2 (self-dual length as the radius of a large extra dimension): the physical reading is REFUTED; the registered rule is NOT YET TRIGGERED.** Both halves matter and are kept separate:
  - LeanMaster Stream 6 records that `R = s ~ 47.0 um` fails all three bounds (Eot-Wash toroidal radius < 30 um, Yukawa range < 38.6 um, neutron-star heating < 44 um) and that the programme's own T-duality fixes `kappa = 1` under the tier-C identification `alpha' = s^2`, so no O(1) factor can rescue it. Its conclusion, quoted: "the dual-scale hypothesis as an extra-dimension or UV/IR statement does not survive the existing data".
  - The rule written at line 58 of this file falsifies P2 only if a gravitational-strength Yukawa deviation is excluded **for every range >= 21.0 um**. The bounds in hand exclude ranges above about 30 um, so the 21-30 um window of the registered bracket is not yet excluded. **No new short-range gravity data were analysed in rounds 1 or 2.**
  - Recorded position: the extra-dimension interpretation is abandoned as refuted by the programme's own verdict; the registered rule stays open and can only be closed by data covering 21-30 um.
- **"Zero free parameters by derivation": CLOSED AS NOT REACHABLE TODAY.** The 2 -> 0 step cannot be made by derivation from the mathematics this programme has verified, for reasons that are recorded, not conjectured: Stream 6 and Stream 7 exclude every physical reading of the dual-scale identification tried so far; Stream 8 states "Observables: none" for its results (lines 112, 153, 272, 374 at eb791e7); and `mu_sym` still has no unit bridge in any harness. A further obstruction is recorded in Stream 8 and is **kernel-checked arithmetic with a tier-C physical reading**: on the full `Gamma_6,22` moduli space of K3 x T2 the maximal gauge enhancement is `SO(44)` with 924 roots, while any point that factorises into a K3 point and a T2 point carries at most 760 + 6 = 766 roots (`product_points_not_maximal`). The reading that global string-gas trapping therefore does not preserve the K3 x T2 factorisation, and so pulls away from the UV attractor `tau = omega`, is **tier C (LeanMaster's own label), not a theorem about our universe**.
- **Zero-parameter status achieved by hypothesis change (M0) remains reported, not pre-registered** (A3), and X4 now disfavours it at 2.16 sigma.

**A8. Defects found in the TDA code, and what they affect.** Disclosed here because two of them touch results reported in earlier addenda.
- `cmb_tda.coarse_stats` counted bins that are empty across the whole simulation ensemble toward the degrees of freedom, and a 1e-8 ridge before `pinv` gave those bins enormous weight. On 100 test maps whose p-values should be uniform, the b1 p-values failed uniformity at KS p = 1.1e-4. **Round-1 CMB chi2 p-values that used this function are withdrawn.**
- `cmb_tda.build_topology` fills every 4-clique of adjacent pixels, so the full-sky complex has b2 = 49147 instead of 1. b0 and b1 are unaffected, but the curve labelled Euler characteristic is b0 - b1.
- **X1 (A6) is not affected by either defect**: it computes empirical rank p-values, p = (1 + #{T_null >= T_data}) / (N + 1), with the variance floored, and uses only b0 and b1 curves, with the same complex for data and simulations.
- Both defects are being repaired on branch `loop/tda-simple` with a regression suite that fails on the old code.

**A9. Reorientation of the TDA search (no prediction is registered by this note).** Searches for smooth continuous defects (cosmic strings) returned nothing at N = 400, N = 25000 and on the CMB. Future TDA work will instead target **residual discrete angular anisotropy** matching the finite symmetry groups that Stream 8 singles out: the order-192 Kummer group `(Z_2)^4 |x A_4` (Frame shapes `1^24, 1^8 2^8, 1^6 3^6, 1^4 2^2 4^4`, verified at `VERIFIED_FOUNDATION.md:194`), the Hurwitz lattice `D_4` and its automorphism group of order 1152, and the attractive Kummer surface with `D = 12` at `tau = omega`. **This is an exploratory lens (tier X), not a prediction**: LeanMaster derives no CMB observable from any of it. Per A5, no threshold is registered here; any such test enters this file with its statistic, null, seeds and decision rule before its data are loaded.

## Addendum 2026-09-21 (appended; A1-A9 and the original rules above unchanged)

**A10. Stream 1 bridge — no prediction is registered by this addendum.** Per A5, nothing enters this
file as a prediction without its statistic, null, seeds and decision rule stated before its data are
loaded. A10 is bookkeeping: it records what the newest Lean 4 formalization in the programme licenses
and forbids, so that a later session does not re-open a route that is already closed. Full evidence
and reproduction commands: `audit/STREAM1_BRIDGE.md`.

- **Pins.** "Stream 1" here is the repository `SocrateAI-DualScaleTopologicalUniverseModel-LeanProposal`
  at commit `bb74acb56f386a97e433f94eb0b2632ed03bc4ca` (Zenodo concept 10.5281/zenodo.22853239),
  **not** LeanMaster's internal 2026-09-15 Stream 1. LeanMaster is pinned at
  `ede49f06800cde8177867c02bb08d44f7be275c5` (v3.44.0). LeanMaster released v3.42 → v3.43 → v3.44 in
  one day and v3.44 amended a v3.43 claim, so these are cited by SHA and not "at HEAD".

- **A10.1 The zero-parameter closure of A7 is sharpened, not reopened.** Dolgachev (1996) Thm 7.1 —
  literature, quoted through Stream 1's `paper/sections/02-preliminaries.tex`, and **not**
  kernel-proved there — states that the coarse moduli space of `M_n`-polarized K3 surfaces is the
  Fricke modular curve `H/Γ₀(n)+`. A modular curve is one complex dimension. So the geometry does not
  contain zero moduli to be derived: it contains **one**. A7 said `2 → 0` is not reachable; A10
  records *why* on the mathematics side, and that **1** is the floor there. Tier L for the moduli
  statement; Tier A for the lattice facts it rests on in Stream 1 (`root_orthogonal_iff_selfdual`,
  `height_ge_two`/`height_fricke`/`height_eq_two_iff`, `rhoAL_isometry`, `rhoAL_det`).

- **A10.2 Reading that modulus as a cosmological free parameter is Tier C and is NOT registered.**
  Stream 1's README states program-wide: "No exact physical observable exists anywhere in this program
  (F5b)"; "The Sym² relation supplies **no physical coupling**"; and the coincidence that one integer
  matrix is both the Fricke involution and the Narain T-duality generator "is a fact about a lattice
  isometry and **not** a physical identification". `mu_sym` still has no unit bridge in any harness
  (A7). **Sym² is therefore closed as a route to parameter reduction**, and this line is written so
  that a future session does not re-open it.

- **A10.3 Corrections to previously reported material.** Two claims in `T_duality_Alone.tex` were
  wrong and are withdrawn there with a dated note: (i) the Fricke point of this model is
  `τ = i/√12` (level 12), not `τ = i` (level 1) — the implemented potential is stationary at the
  former and not the latter; (ii) the `SL(2,ℤ)` fold was described as keeping trajectories in the
  physical domain, but `S : τ ↦ −1/τ` is not a symmetry of the implemented potential
  (`max |V(Sτ) − V(τ)| = 8.77`, against `6.7e-16` for `T`), so it moved reported trajectory points to
  physically inequivalent points. The `S`-step is now off by default. **This changes solver output**, and the blast radius was
  measured rather than left open: **no committed simulation artifact is affected**, because every
  tracked output is written by `workshopcosmo.py`, `parameter_sweep*.py`, `scripts/tda_mapper.py` or
  the Rust simulator, none of which import `leanflow` (verified in a clean interpreter). On the
  `leanflow` solver path itself the change is large — the old `S`-fold moved **100% of trajectory
  points**, by up to `3.175` in `y` — so numbers produced on that path before 2026-09-21 should not be
  trusted. (iii) A third defect, **S1-F8**, was found while making that measurement: the solver applied
  both projections at state indices `(0, 1)`, but the cosmology state is `[a, x, y, u, v]`, so it was
  folding the **scale factor** as `Re τ` (reporting 1.98e9 folds on a `t ≤ 200` run) and clamping
  `Re τ` — legitimately `0` at the Fricke point and `0.5` at the orbifold point — to `1/√12`, while
  never touching the real modulus at index 2. The solver now takes explicit `modulus_indices` and
  skips with a warning rather than guessing an unknown layout.

- **A10.4 A defect disclosed in paper-support code.** `papers/T-dulaity alone/t_duality_calculus.py`
  printed "Symmetric square modular invariant verified for `L_3 = Sym^2 L_2`" conditioned on
  `c_eff == 1701`, an arithmetic identity in the hard-coded `E_0 = −425/6` alone that could never
  fail. Stream 1 proves the relation **with a prefactor**, `L₃ = P₂·Sym²(L₂)`, `P₂ = 1 − 26z − 27z²`,
  so the printed form was also wrong as stated. Corrected. No source is recorded anywhere in this
  repository for `E_0 = −425/6`; it is now marked unsourced.

- **A10.5 Conventions now binding on future work** (no current occurrence in this repository, recorded
  before one is written): `Sym²` is **contravariant**, `sym2(MM') = sym2(M')sym2(M)`, with Stream 1's
  negative control `sym2_not_covariant`; and the two rank-3 lattices of signature (2,1) —
  `⟨1⟩ + U(2N)` with `det = −4N²` and `U + ⟨2N⟩` with `det = −2N` — are **not isometric**
  (`no_isometry_G0N_TN`) and must not be named interchangeably.
