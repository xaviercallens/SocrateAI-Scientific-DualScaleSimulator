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
