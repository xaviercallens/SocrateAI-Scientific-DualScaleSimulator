# Future observational targets: what the reverse-to-zero campaign leaves to whoever comes next

**Sealed 2026-09-20.** Branch `loop/reverse-zero`. This memo is the transmission document of the campaign. It records one methodological lesson that cost a full run to learn, the two concrete blockers that stop the cosmological side today, and the state every future run should start from.

Everything here is tier X (numerics) unless marked otherwise. Nothing in this memo is derived from K3 × T².

---

## 1. The lesson of 3-D TDA: use the COUNT, not the spacing

The quantum-fluid validation succeeded on real scanning-tunnelling data of vortex cores in an a-Re₆Zr film (2-D), carried by the **spacing statistic**: the interquartile spread of H₀ death radii divided by their median. It recovered the published lattice constant (36.0 nm against 34.6 nm) and vortex counts (ratio 1.005–1.106).

That statistic was transferred faithfully to cosmology — the transferred code reproduces the Re₆Zr values on the saved point clouds **to machine zero on all 11 fields** — and then measured against injected signals in both geometries. The measured outcome:

| Geometry | Spacing (IQR/median) | Count | Orientational (ψ₆ / Q₆) |
|---|---|---|---|
| **Sphere** (CMB, nside 128) | first to fire: 0.97 power at A = 4σ, 1000 cores | 0.67 at the same cell | 0.09 at the same cell |
| **3-D** (DESI, volume-limited) | **never reaches 95 %** (best 0.17) | **saturates at A = 1σ_δ** | **saturates at A = 1σ_δ** |

**The ordering inverts with dimension.** In three dimensions the statistic that carried the Re₆Zr success collapses, while the count and the orientational order saturate at the lowest amplitude tested.

**Instruction for any future search for a 3-D topological-defect imprint in galaxy data: build the test on the COUNT of detected cores against a clustering-matched null, with an orientational statistic as the second leg. Do not make the spacing spread the primary statistic.** The reason is structural, not a tuning accident: the spread measures regularity, so a denser injected population pushes the cloud back toward the null's own regularity. The same non-monotonicity in amplitude was measured independently for the crystallinity statistic, which *falls before it rises*.

**Second-order lesson, same origin.** The spacing spread measures regular spacing, not orientational order: on both the Re₆Zr and the simulated condensate data it stayed small while the orientational parameter ψ₆ collapsed. Never read it as evidence of order without an explicit orientational check.

---

## 2. The two blockers on the cosmological side

### (a) DESI-footprint mocks — required to lift Test B out of INCONCLUSIVE

Test B (DESI DR1 BGS_BRIGHT-21.5, volume-limited) ran its detector and found **832 density peaks in 1.13 × 10⁹ (Mpc/h)³**, with spacing spread 0.3010 and Q₆ 0.3137. **No p-value was quoted, and none should be, until a clustering-matched null exists.**

- Why the official randoms are **not** a substitute: they contain no clustering, so a peak detector applied to them gives systematically more peaks and a smaller spread (1235 ± 23 peaks against the data's 832). The gap measures galaxy clustering, not defects. Substituting them would manufacture a detection.
- What is missing, precisely: **lognormal (or N-body) mock catalogues on the DESI footprint**, with the survey's own n(z), angular mask and redshift-space distortions. The existing mock code (`audit/reverse_zero/E5-cosmic-web-tda-scaled/x2_lib.py`) is wired to the SDSS DR17 window and was not ported; that port is the single piece of work standing between this campaign and a DESI verdict.
- A second requirement rides on the same work: round 2's X2 found the SDSS-window lognormal family **mis-calibrated at 2.88σ** against the data's own two-point function. A DESI mock family must pass that same calibration gate before any topology p-value is read from it. If it fails, the correct answer is again INCONCLUSIVE.
- Also note the measured bias in the injection study: DESI sensitivity was calibrated against a **clustering-free** baseline, so the quoted "95 % at A = 1σ_δ" **overstates** sensitivity against the real clustered field.

### (b) The southern hemisphere — required for replication

Only the DESI northern cap (NGC) was analysed. The southern cap (SGC) was never run, so **the replication leg of Test B was never exercisable**. Any future candidate signal in one cap must be reproduced in the other, on data that were not used to choose the detector or its threshold. The data are already staged: `BGS_BRIGHT-21.5_SGC_clustering.dat.fits` and its randoms, in `/mnt/disks/disk-socrateai-local-1/dualscale-data-r3/desi_dr1_lss/`, with sha256 in `audit/data_r3_manifest.json`.

---

## 3. What was measured, so nobody re-derives it

**CMB, Test A: NULL, with a weak and quantified sensitivity.** WMAP 9-yr ILC (primary): p = 0.535 (spacing), 0.523 (count), 0.507 (orientational); Planck SMICA (secondary, the same sky, not an independent test): 0.144 / 0.950 / 0.371. Bonferroni threshold 0.0083 over the registered family of 6. The variance gate passed on both maps.

- **Smallest injected population detected at 95 %: A = 4 σ_T with 1000 full-sky cores** (~397 inside the eroded mask).
- **The 0.73σ figure is not the sensitivity limit**; it is its cause. A 1 σ_T point-like spot survives the detector's own smoothing at only **0.73 σ of the smoothed field**, against a 1.0 σ threshold (loss factor 0.500; σ(T_smoothed)/σ_T = 0.688). Anything below ≈ 1.4 σ_T therefore cannot become a candidate at all. A future search that wants weak cores must change the detector, not the statistic.
- False-positive row at zero amplitude: 0.07 / 0.01 / 0.04 per statistic, 0.12 union, against a nominal 0.05 / 0.143.
- **So the null excludes a strong, sparse population and is close to vacuous against a weak or dense one.** Do not cite it as a constraint without this sentence.

**Sky-symmetry lens: NULL, and one structural correction.** The order-192 Kummer group `(ℤ₂)⁴ ⋊ A₄` **does not act on the sphere** — the `(ℤ₂)⁴` factor is 2-torsion translation on the abelian surface. Only `A₄` acts, of order 12, through the binary-tetrahedral double cover; the Frame class `4B` has no spherical counterpart. Results: WMAP A₄ p = 0.627, SMICA A₄ p = 0.577. Sensitivity needs the symmetric component to carry **tens of percent of the total ℓ ≤ 64 power**, and the specificity is asymmetric (an injected C₁₂ leaks into the A₄ statistic at 0.37, but not the reverse).

**The hard limit on internal-space topology, measured on branch `loop/tda-k3t2`.** Sampled persistent homology recovered flat T² and T³ but **no 4-dimensional space at all** within budget: not T⁴, not T⁴/Z₂, not K3 (whose b₂ plateaued at 27, never 22), not even S² × S² whose b₂ is 2. Required sample size grows ≈ 16× per dimension. A naive persistence-ratio rule fired on 9 of 12 density-matched nulls. **No point-sample TDA claim about K3 × T² topology in cosmological data is supportable**, and a non-null signal would not single out K3 × T².

---

## 4. Tooling state for the next run

- **Use the fixed TDA library**: `audit/tda_validation/tda_fixed/` on branch `loop/tda-simple` (source commit `d8175f1`). Full-sky Betti is now (1, 0, 1); dead and atomic bins are dropped with a recorded reason; degrees of freedom are the retained rank. **Prefer the rank p-value**: the χ² branch remains mildly anti-conservative (0.058 against 0.05).
- **Do not import** `audit/reverse_zero/E5-cmb-tda/cmb_tda.py` for new work. It carries both defects and is annotated as the historical record.
- Round-1 CMB χ² p-values that used the old statistic stay **withdrawn**. X1 and the vorticity run are unaffected (rank p-values, b₀ and b₁ only).
- Every run in this campaign committed its registration **alone, before reading data**. Keep that: it is what made the INCONCLUSIVE verdicts possible to state without embarrassment.

---

## 5. Ranked next steps

1. **Port the mock generator to the DESI footprint** and pass the 2.88σ calibration gate. Without it there is no galaxy verdict, only absolute numbers.
2. **Run the southern cap** and treat it strictly as replication.
3. **Rebuild the count leg against an independent expectation.** The Re₆Zr count test worked because B·A/Φ₀ gave an expected number from outside the data. Cosmology has no such number here: no theory in this programme predicts a core density. Until one exists, the count test compares against a null, not against a prediction — a weaker test, and it should be labelled as such.
4. **If a detector for weak cores is wanted**, fix the 0.73σ smoothing loss first; no statistic can recover what the detector discards.
