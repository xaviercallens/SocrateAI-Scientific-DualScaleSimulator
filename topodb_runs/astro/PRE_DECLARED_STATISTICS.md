# Pre-declared statistics, nulls and rules — astro persistent-homology campaign

**Committed before any dataset in this campaign was opened.** The prior campaign
(`docs/FUTURE_OBSERVATIONAL_TARGETS.md`, §4) credits exactly this practice for
making its INCONCLUSIVE verdicts statable without embarrassment. Everything below
is tier **X** (numerics). Nothing here is "proved".

## 0. Betti-extraction rule

`common.gap_betti`: sort finite persistences descending, count the bars above the
largest multiplicative gap. Parameter-free.

**Scope limit, declared in advance:** this rule is used **only** in
`step0_known_answers.py`, where the true topology is known and `expected` is
pre-stated. It is **not** a detector. The prior campaign measured that a naive
persistence-ratio rule *fired on 9 of 12 density-matched nulls*. On real data the
primary record is the COUNT statistic against an explicit null, with `expected`
left NULL in the `betti` table.

## 1. Three-dimensional point clouds (DESI, SDSS) — GUDHI AlphaComplex

Filtration value is the squared circumradius; every statistic is quoted after
`sqrt`, i.e. in the length units of the input coordinates (Mpc/h).

Primary statistic (the memo's §1 lesson: **use the COUNT, not the spacing spread**):

| name | definition |
|---|---|
| `n_h1_bars_over_5mpc` | number of H1 bars with persistence > 5 Mpc/h |
| `n_h2_bars_over_5mpc` | number of H2 bars with persistence > 5 Mpc/h (voids) |
| `total_persistence_h1` | sum of H1 persistences (Mpc/h) |
| `total_persistence_h2` | sum of H2 persistences (Mpc/h) |
| `max_persistence_h2` | longest H2 bar (Mpc/h) |
| `h0_death_iqr_over_median` | the Re₆Zr *spacing* statistic — recorded as a **secondary**, because §1 of the memo measured that it never reaches 95 % power in 3-D |

Null: **N disjoint subsamples of the official random catalogue**, of the same size,
drawn from the same survey window and n(z).

**The randoms are an UNCLUSTERED Poisson null.** They contain no galaxy
clustering. A p-value against them therefore tests *"does this point set differ
topologically from an unclustered field with the same selection function?"* — a
question whose answer is known in advance to be yes, because galaxies cluster. It
does **not** test for any anomaly with respect to ΛCDM. The memo states that
substituting the randoms for a clustering-matched null "would manufacture a
detection"; this campaign therefore records the null's identity in the
`null_model` string itself and writes the finding as a statement about clustering
detectability, never as a topological anomaly. A clustering-matched (lognormal or
N-body) DESI-footprint mock family does not exist on this disk and is **NOT
ATTEMPTED** here.

Distances use flat ΛCDM with Ω_m = 0.3137, h = 0.6736 (DESI DR1 fiducial /
AbacusSummit base). The cosmology sets the unit of length and nothing else; no
result below depends on it.

## 2. Sky maps (Planck SMICA, WMAP ILC, Haslam, IRIS/COBE) — FIXED complex

`lib/cmb_topology.build_topology_fixed` + `betti_curves_from_topology`, lower-star
filtration on the field normalised by σ of the unmasked pixels (this is what makes
K_CMB and MJy/sr maps comparable at all). ν grid: 41 points on [−3, 3].

| name | definition |
|---|---|
| `b0_at_nu_1` | b0 of the sublevel set at ν = +1 |
| `b1_at_nu_0` | b1 of the sublevel set at ν = 0 |
| `b0_curve_integral` | Σ b0(ν) Δν |
| `b1_curve_integral` | Σ b1(ν) Δν |
| `euler_char_at_nu_0` | true Euler characteristic #V−#E+#F at ν = 0 |
| `coarse_rank_p_b0`, `coarse_rank_p_b1` | `lib/stats.coarse_stats_fixed` rank p over 8 coarse bins |

Null: Gaussian realisations synthesised from **each map's own masked pseudo-C_ℓ,
divided by f_sky** to undo the mask's power suppression, then masked with the same
mask. n stated per run. **Rank p-values only** — the χ² branch of the fixed
library remains mildly anti-conservative (0.058 against 0.05).

`null_calibration` control on every such run: one null realisation is held out and
ranked against the remaining n−1; its p must not be extreme. Without this control
the p-value is uninterpretable, because a mis-normalised null makes the data look
anomalous for free.

**Foreground maps (Haslam 408 MHz, IRIS/COBE) are a POSITIVE control.** They are
strongly non-Gaussian, so a Gaussian null must be rejected at the rank floor. If it
is not, the CMB nulls are vacuous and must be read as "the pipeline cannot fire",
not as "no signal".

## 3. CAMELS HI maps — GUDHI CubicalComplex

Sublevel cubical persistence on log10(HI + ε) of each 256×256 map, ε = the map's
own positive minimum. Per-map statistics (the **full list, fixed here before any
correlation was computed**):

1. `n_h0_bars`  2. `n_h1_bars`  3. `total_persistence_h0`
4. `total_persistence_h1`  5. `max_persistence_h1`  6. `mean_persistence_h1`

Correlation: Spearman ρ of each statistic against Ω_m and against σ_8.
**6 statistics × 2 parameters = 12 tests; multiplicity = Bonferroni/12,
threshold 0.05/12 = 0.00417.** Permutation null: labels permuted, n stated.

Effective-N check, run **before** any correlation: if the 750 test maps carry
fewer than 750 distinct label rows, maps share a simulation and a per-map
permutation is anti-conservative; the permutation is then done at the
simulation-group level and N_eff is reported.

Known-answer control: the same Spearman on **shuffled labels** must show no
correlation surviving the threshold.

## 4. What a p-value may say here

`n_null = 20` gives a rank floor of 1/21 = 0.048. A p at the floor is
**resolution-limited, not significant**, and is recorded as such in `p_method`.
No p-value is stored without its `null_model` and `n_null` (TopoDB refuses it).
Where no null exists, the number is stored **without** a p-value.
