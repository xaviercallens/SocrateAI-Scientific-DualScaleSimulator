# From quantum fluids to cosmology: what transfers, what does not, and one testable conjecture

**Date:** 2026-09-19. **Branch:** `loop/reverse-zero`. **Tier of everything proposed here: C (conjecture).** Nothing in this note is derived from K3 × T², and nothing in it has been tested against cosmological data yet.

## 0. A correction of framing, stated first

The quantum-fluid work of 2026-09-19 (`loop/tda-validation`, `audit/tda_validation/quantum_fluid/report.json`) **validated an instrument, not a theory**. What succeeded is this project's alpha-complex TDA path, which recovered published vortex physics from real scanning-tunnelling data on an a-Re₆Zr film:

| Quantity | Published | Recovered by our pipeline |
|---|---|---|
| Vortex lattice constant at 20 kOe | 34.6 nm (1.075·√(Φ₀/B)) | **36.0 nm** (ratio 1.043) |
| Vortex count | B·A/Φ₀ | ratio **1.005–1.106** (the 3 kOe excess of +10.6% matches the paper's stated 10–15%) |
| Re-entrant order–disorder sequence across field | ordered near 20 kOe | H₀-spread minimal at 20 kOe (0.069), larger at 3 and 70 kOe |

**That is a statement about our software.** No dual-scale or T-duality prediction was tested by it, and no superconductor result says anything about the expansion of the universe. The sentence "the dual-scale proposal triumphed in quantum fluids" is not supported by the evidence in this repository, and this note does not adopt it.

Two further results from the same validation constrain what may be claimed here:
- On rough fields the **lower-star path failed its own negative control**: site-shuffled fields reproduced the XY-model "transition", so that signal came from the value distribution, not from topology.
- The **alpha-path H₀ spread measures regular spacing, not orientational order**: it stayed small while the orientational order parameter ψ₆ collapsed.

## 1. P2 as a scale error: the hypothesis, and its cost

**The claim to be examined.** P2 measured the wrong object: a self-dual length of ≈ 47 µm was tested as the radius of a *global* extra dimension, i.e. a property of the whole bath, when it could instead be a property of *defect cores* — the size of a vortex-like object in which T-duality is realised — which no torsion-balance experiment integrates over.

**What supports it.**
- The falsification is unambiguous and is the programme's own: LeanMaster Stream 6 records that R = s ≈ 47.0 µm fails all three bounds (Eöt-Wash toroidal radius < 30 µm, Yukawa range < 38.6 µm, neutron-star heating < 44 µm), and that the programme's T-duality fixes κ = 1 under the tier-C identification α′ = s², so no O(1) factor can rescue it. Verbatim: *"the dual-scale hypothesis as an extra-dimension or UV/IR statement does not survive the existing data"*.
- A bound on a homogeneous extra dimension does not automatically bound a rare, localised core. Short-range gravity experiments constrain the Yukawa deviation of the *mean* field between macroscopic plates.

**What it costs, and this must be paid before the idea is used.**
1. **It is a reinterpretation after a failure.** P2 was pre-registered as a global length and refuted. Re-describing the same number as a core size is not a new prediction; it is a rescue. It becomes science only when it is frozen as a new hypothesis (P1′ in Stream 6's language) with an observable and a threshold written down before any comparison.
2. **The registered rule stays formally open.** `PRE_REGISTRATION.md:58` falsifies P2 only if a gravitational-strength Yukawa deviation is excluded for *every* range ≥ 21.0 µm; the bounds in hand reach to ≈ 30 µm, so the 21–30 µm window is untested. The interpretation is abandoned, the rule is not yet closed.
3. **A core-size reading owes an abundance.** A gas of 47 µm objects is not automatically invisible: it must be shown what number density and coupling evade Eöt-Wash, MICROSCOPE, neutron-star heating and fifth-force searches. Until that calculation exists, "it hides in the cores" is not a defence.

## 2. The conjecture: dark energy as topological friction

**Statement (tier C).** Let the early universe contain a gas of micro-defects whose cores realise the self-dual scale. Expansion stretches this network. The work done against the network's tension, dissipated as the defect gas is stretched and reconnects, appears in the stress tensor as an effective **bulk viscosity**, i.e. an isotropic negative pressure that mimics dark energy.

**Why the analogy is attractive.** In superfluid turbulence a tangle of quantised vortices dissipates energy through reconnection and Kelvin-wave cascades, and the tangle's coarse-grained effect on the fluid is described by mutual friction. The mathematics of a defect network coupled to a background is standard in cosmology too: strings and walls obey a velocity-dependent one-scale model with a friction term.

**Why it is only a conjecture here.** Three gaps, each of which is enough to block a claim:
1. **No derivation.** Nothing in this programme derives the defect gas, its tension, its density or its reconnection rate from K3 × T². LeanMaster records "Observables: none" for its Stream 8 results, and no unit-bearing quantity follows from the verified mathematics.
2. **Bulk viscosity is severely constrained.** An effective w ≈ −1 from bulk viscosity is generally accompanied by entropy production and by a scale dependence that shows up in the CMB and in structure growth. A conjecture that only reproduces w ≈ −1 today, with free parameters, is weaker than ΛCDM, which reproduces it with one constant.
3. **Known defect networks are already bounded.** Scaling string networks are limited to Gμ ≲ 10⁻⁷ (`planckGmuBound = 1.5e-7` in LeanMaster, a literature bound). Any new network must say why it evades these limits.

**What would make it testable.** The conjecture must produce at least one number with units, frozen before comparison. Candidate observables, in the order that this project could reach them:
- an equation-of-state history w(z) with a stated shape and one fewer free parameter than CPL, confronted with DESI DR2 + Pantheon+ by the existing X4 pipeline;
- a bulk-viscosity-induced deviation of the growth rate fσ₈(z), which DESI and eBOSS already measure;
- a defect-gas contribution to the stochastic gravitational-wave background at NANOGrav frequencies, which the X3 pipeline could constrain with the same machinery that produced the c4 interval;
- a **spatial statistic of defect cores**, which is what §3 makes operational.

Until one of those exists with a threshold written down in `PRE_REGISTRATION.md`, "topological friction" remains a name for an intuition.

## 3. The new lens: search for a discrete defect-core population, not for smooth strings

**Why the old lens failed.** E4 (N = 400) and E5 (N = 25 000) looked for large smooth structures in the cosmic web and found nothing that survived its own controls; round 2's X2 could not even calibrate its mocks (2.88σ) and returned INCONCLUSIVE. On the sky, X1 does not reject Gaussian isotropy for WMAP (p = 0.47).

**Two reasons, and only one of them is physics.**
1. *Possible physics:* if the infrared dynamics smooths the geometry (Stream 8's tier-C reading around the SO(44) point), a coherent network of long smooth defects need not survive.
2. *Method, and this one is measured:* sampled persistent homology **cannot recover 4-dimensional topology at feasible sample sizes at all**. On branch `loop/tda-k3t2`, PH recovered flat T² and T³ but failed on T⁴, on T⁴/Z₂, on K3 (where b₂ plateaued at 27, never 22) and even on S²×S² whose b₂ is 2, with the required sample size growing ≈ 16× per dimension. A naive persistence-ratio rule fired on 9 of 12 density-matched nulls.

**The transfer that is legitimate.** What worked on the Re₆Zr data was not "topology" in the abstract: it was a specific, reproducible pipeline applied to a **point cloud of defect cores** — detect cores, build the alpha complex, and read two statistics: the spread of H₀ death radii (IQR/median, a regularity-of-spacing measure) and the count against an independent expectation. Those statistics are dimension-appropriate (2-D or 3-D point sets), they have published answer keys, and they passed on real data.

**So the cosmological search is reconfigured as follows** (to be pre-registered before it is run):
- **CMB (nside 128, WMAP ILC and Planck SMICA):** extract a point cloud of candidate defect cores from the temperature field, using a detector with a stated rule (e.g. local extrema above a threshold in the wavelet- or gradient-filtered map, or Kaiser–Stebbins step candidates), then apply the *same* alpha-path statistics with the *same* hyperparameters used for Re₆Zr. Null: the spectrum-matched Gaussian simulations already built and gate-validated by X1 (max |z| = 0.85 against 2.81).
- **Galaxy survey (DESI DR1 BGS_BRIGHT-21.5, which is volume-limited, plus its official randoms):** the same statistics on a point cloud of density peaks or void centres, with the survey's own random catalogue and the CAMB lognormal mocks as nulls.
- **Mandatory controls, from the validation work:** a site-shuffled field control for anything lower-star; an explicit orientational statistic alongside the H₀ spread, since the spread alone misses orientational disorder; an absolute persistence floor, since ratios alone are meaningless; and a null-injection row at zero amplitude to measure the false-positive rate.
- **Honest prior.** Nothing in the verified mathematics predicts such a population, its density or its scale. A null result is the expected outcome and is worth recording; a non-null result would be a signal in need of an astrophysical explanation first, and would not single out K3 × T².

## 4. Status of the TDA code these searches depend on

Both defects found by the known-answer suite are now repaired on `loop/tda-simple` (`audit/tda_validation/tda_fixed/`), with a regression guard that fails on the originals (19 failures) and passes on the fix (30 tests):
- `build_topology`: the full-sky complex now gives **Betti (1, 0, 1)** at nside 32 and 64 instead of b₂ = 49 147; a disk mask gives a disc; every edge lies in at most two triangles; the Euler key is now the true Euler characteristic.
- `coarse_stats`: dead and atomic bins are dropped with a recorded reason, the degrees of freedom are the retained rank, and the pooled rank p-value is the recommended statistic. Residual limitation, disclosed: with an estimated covariance the χ² branch remains mildly anti-conservative (0.058 instead of 0.05).
- Round-1 CMB χ² p-values that used the old function stay withdrawn. X1's results are unaffected, because X1 uses rank p-values and only b₀ and b₁.

## 5. First result of the discrete-symmetry lens: NULL, with a measured sensitivity

Run on `loop/tda-validation` (`audit/tda_validation/crystallography/`, commits `33e4bca`..`810ff51`), pre-registered before the maps were read.

**A correction that the implementation forced, and that matters for every future statement.** The order-192 Kummer group `(ℤ₂)⁴ ⋊ A₄` **does not act on the sphere**. The `(ℤ₂)⁴` factor is 2-torsion translation on the abelian surface, not a rotation. Only `A₄` can act, and it is reached from the 24 Hurwitz units (the binary tetrahedral group 2T) through the double cover q ↦ (v ↦ q v q̄), whose image has **order 12, not 24**. So the sky statistic tests `A₄`, with `C₁₂` as the control group. The Frame-shape class `4B` has no counterpart on the sphere at all: ±i, ±j, ±k map to π rotations of order 2, and the order-4 structure lives in the double cover.

**Result on the real maps** (200 Gaussian nulls per map from each map's own masked spectrum, frozen 192-point orientation grid, rank p-values, Bonferroni over the family of 4):

| Map | Group | Z | rank p |
|---|---|---|---|
| WMAP 9-yr ILC (primary) | A₄ | 0.761 | **0.627** |
| WMAP 9-yr ILC | C₁₂ | 1.675 | 0.224 |
| Planck PR3 SMICA | A₄ | 0.858 | **0.577** |
| Planck PR3 SMICA | C₁₂ | 1.647 | 0.214 |

Nothing approaches the 0.0125 threshold. The two maps are the same sky and are not independent tests.

**How to read this null, honestly.** The measured sensitivity is poor: the statistic is a power *fraction*, so a coherent component interferes linearly before it dominates quadratically, and crystallinity **falls before it rises** as the injected amplitude grows. The 95% detection point needs the symmetric component to carry **tens of percent of the total ℓ ≤ 64 power** (narrow ℓ ∈ {3,4,6} pattern: 0.90 power at a 0.24 fraction; broadband: only 0.54 at 0.30). Therefore:
- the null **excludes only a very large discrete-symmetry component**, and is close to vacuous against a weak one;
- the specificity is **asymmetric**: an injected A₄ signal never fires the C₁₂ control, but an injected C₁₂ signal leaks into the A₄ statistic (power 0.37 at the largest amplitude), so a hypothetical A₄ detection could not by itself be attributed to A₄;
- and it falsifies nothing about K3 × T², because nothing connects the two. LeanMaster still records, verbatim, *"(iii) Observables: none — N = 4, non-chiral."*

## 6. Next step already defined: the vortex-core statistic, transferred from Re₆Zr

The lens above tests a *global* symmetry of the field. The quantum-fluid validation succeeded on something else entirely: a **point cloud of defect cores**, with two statistics that have published answer keys (spacing regularity through the H₀-death-radius spread, and core count against an independent expectation). That pipeline, with its hyperparameters, is being transferred to the CMB and to DESI as a search for a discrete population of defect-core candidates, with the nulls and mandatory controls listed in §3. It is registered before it is run, and a null result is the expected outcome.
