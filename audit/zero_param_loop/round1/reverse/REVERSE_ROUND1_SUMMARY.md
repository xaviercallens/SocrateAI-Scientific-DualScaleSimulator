# Reverse loop, round 1: reduced (5-param) model -> TDA -> experiment -> hypothesis

Every number below is reproducible from a script in this directory. Commands
and output files are named per step.

## 1. Reduced model

Applied ONLY the round-1 tier-B accepted reduction (`../reduce_proposals.json`
`accepted_this_round_tier_B`): `(pta_suppression, c4_c0_ratio)` collapse to
`c4_pta_product = c4_c0_ratio * pta_suppression`, since no observable this
harness computes depends on anything but that product (round-0 Check C
negative control, max diff 1.11e-16; round-1 Jacobian: exact SV=0.0 null
direction at all 6 probed points). `lambda_sym` is **kept** in the free
count (0.5-evidence / "fixed by convention" per ground rule (c) — insensitivity
is not derivation).

**Free parameters still remaining: 5 — `a_pot, b_pot, mu_sym, lambda_sym, c4_pta_product`**
(down from 6; not 4, not 2 — the task's target ladder is a multi-round plan,
not a claim reached this round).

- `reduced_model.py` — wrapper; splits `c4_pta_product` arbitrarily as
  `sqrt(P), sqrt(P))`, verified inconsequential:
  `python reduced_model.py` -> `ALL_PASS=True` (diff vs full model and vs an
  alternate split both < 1.4e-17).

## 2. TDA — pre-registered prediction, then measurement

`tda_prediction.json` was written (file mtime 17:22) **before**
`reduced_sweep.py` (mtime 17:28) or `reduced_tda.py` (mtime 20:13) ran —
verifiable by file mtimes even though all of `reverse/` lands in a single
git commit after the fact. It records both the task's literal instruction
("dimension = remaining free params; zero params -> single point") and a
corrected prediction carried from round 1's own `interpretation_caveat`:
Betti numbers of a smooth image of a still-continuous (>0-parameter) domain
are `[1,0]` regardless of how many parameters remain, and the number that
actually tracks parameter *content* is the Jacobian-SVD effective dimension,
predicted to stay at **3** (unchanged from round 1's `best_fit` probe),
because the collapsed direction was already an exact `SV=0.0` null direction
in the 6-param model.

Command: `python reduced_sweep.py` (seed 43, 521 points, 0 errors, 388
numerically stable) -> `reduced_sweep.csv`.
Command: `python reduced_jacobian.py` -> `reduced_jacobian_report.json`:
`best_fit` `effective_dimension_at_1e-3_of_max = 3` (singular values
`[4.177, 2.543, 0.751, 5.0e-4, 1.9e-14]`) — **prediction confirmed at the
best-fit probe** (the same probe round 1's headline number came from). The
5 random numerically-stable points give effective dimensions
`[2, 3, 4, 3, 2]` — not uniformly 3; reported here rather than only citing
the on-message number.
Command: `python reduced_tda.py` (real `import gudhi` 3.13, `RipsComplex` +
`create_simplex_tree` + `persistence` + `bottleneck_distance`, subsampled to
N=316 to match round 1's N exactly) -> `reduced_tda_report.json`:
`betti_numbers_observed = [1, 0]`, `betti_matches_corrected_prediction = true`,
`effective_dimension_matches_prediction = true`.
Circle control still asserts one dominant H1 bar (ratio 34.3) — pipeline
sane.

Bottleneck distances vs round 1's own model diagram: H0=0.039, H1=0.027 —
small, consistent with the pre-registered "seed + density difference, not a
topology change" prediction (no qualitative shape change).

**Open flag from round 1, resolved this round (retraction):** round 1 found
`null_vs_real_sn_H1 (0.0399) > model_vs_real_sn_H1 (0.0195)` and explicitly
flagged it as untested against a variance-matched null / possible
bounding-box artifact. This round built that control
(`null_variance_matched`, per-axis variance matched to the model cloud's own
post-PCA spectrum instead of a uniform box) and found
`null_variance_matched_vs_real_sn_H1 = 0.0203` vs
`model_vs_real_sn_H1 = 0.0271` — **the ordering flips**.
`reduced_tda_report.json["variance_matched_null_check"]["ordering_..._survives_variance_matching"] = false`.
**Round 1's implicit suggestion that the model cloud sits topologically
closer to real SN data than a null cloud does was a bounding-box artifact
and is retracted.** No model-vs-real topological preference is supported by
either round's data.

## 3. Experiment — chi2/AIC/BIC vs full model and vs LCDM

`reduced_chi2.py`:
- **Byte-level equivalence proof** (not an assumption): the reduction lives
  entirely in the PTA sector (`workshopcosmo.py:766`), which BAO/SN chi2
  never touches. Evaluating the reduced model at forward's own best-fit
  `(a_pot, b_pot)` reproduces forward's `w0_cpl_latetime`/`wa_cpl_latetime`
  to `<4e-10` and `chi2_total` to `<1e-6`
  (`equivalence_check.byte_level_equivalence_confirmed = true`) once the
  script's `Z_GRID`/`Z_POINTS` are set to match round 1's `param_sweep.py`
  worker init exactly — an earlier version of this script, without that
  override, gave a spurious 5.1-point chi2 discrepancy from grid-boundary
  clipping at DESI's z=2.330 row; caught and fixed, not silently absorbed.
- `chi2_reduced = 703.8819056575242` vs `chi2_full = 703.8819057544594`
  (forward-reported) -> `delta_chi2_reduced_minus_full = -9.7e-8`.
  `fit_degraded = false` (rule: delta > 2 per removed param; threshold here
  is 2). **This has zero discriminating power**: it shows the reduction
  didn't break an unrelated fit, not that the reduction is supported by
  data (PTA chi2 remains ABSENT — no real angular-correlation dataset this
  round either).
- Independent LCDM recompute (own implementation, same public DESI+Pantheon+
  files, not importing round 1's code): `chi2_total = 697.3402981903151`,
  consistent with round 1's reported value to `abs_diff = 0.0`.
- AIC/BIC, **k_fit convention (primary)** — k = fit-relevant scanned params
  + 2 identical nuisances (s_BAO, offset_SN) for all three models, so only
  real data content survives in the deltas: `k_fit_full=k_fit_reduced=4`,
  `k_fit_lcdm=3`. `delta_aic(reduced-full) = -9.7e-8` (~0, as expected).
  `delta_aic(reduced-lcdm) = +8.54`; `delta_bic(reduced-lcdm) = +13.74` —
  **LCDM is still preferred over the reduced model by essentially the same
  margin round 1 found for the full model.** The reduction changes this
  verdict in neither direction.
- A nominal-k bookkeeping variant is reported separately and labeled
  explicitly as **not data content** (it trivially gives `delta_aic=-2` for
  any 1-parameter drop, by construction).

## 4. Lean-pinned nested test (tier-C, secondary, NOT counted as a reduction)

`lean_pinned_test.py` searches whether LeanMaster's `dark_energy_w0:=-1,
dark_energy_wa:=0` (`DualScaleValidation/Observables.lean:97-98`, re-verified
present at those exact lines this round via
`grep -n 'dark_energy_w0\|dark_energy_wa' .../Observables.lean`) is reachable
in the model's `(a_pot,b_pot)` image, using multiple starts checked against
the searched bounds.

**Reachable, almost exactly**: at `(a_pot=0.0985, b_pot=0.1111)`,
`w0_achieved=-1.0`, `wa_achieved=-1.15e-17` (residual `1.3e-34`) — the model
has a wide basin where the CPL fit lands essentially exactly on `(-1,0)`.
But at that point, `chi2_total = 6483.82` vs the free best fit's `703.88` —
**`delta_chi2 = +5779.9` for 2 fewer fitted parameters** (threshold for
"not degraded" is 4). `supported_by_data = false`. This is a stronger,
more direct falsification than round 1's (round 1 only showed the free
best-fit point doesn't land on `(-1,0)`; this round shows the point that
DOES land on `(-1,0)` is catastrophically disfavored by BAO+SN).

## 5. Hypothesis feedback

The parameter removed this round (the `pta_suppression`/`c4_c0_ratio` split)
lives entirely in a sector with **no real data to test it against** (PTA
chi2 remains ABSENT both rounds) — removing it costs nothing
(`delta_chi2~0`) and proves nothing about the theory's viability; it is a
bookkeeping simplification, not a validated derivation.

The parameters the data *can* see — `a_pot, b_pot` — are exactly the ones a
LeanMaster theorem wants pinned at `(w0,wa)=(-1,0)`. This round shows that
target is reachable (a wide basin of `(a_pot,b_pot)` hits it almost exactly),
but BAO+SN pull the best fit far away from that basin: the data-preferred
point sits at `(w0,wa)=(-0.694,+0.955)` (`chi2=703.88`), while the
Lean-pinned point costs `+5780` in chi2. **Pinning `a_pot,b_pot` to the Lean
constant is not merely "not adopted" (round 1's caution) — it is actively
disfavored by the same public data this loop already fetched.** Meanwhile
flat LCDM beats even the *free* best fit by `~8.5` AIC units. Reaching
zero free parameters therefore cannot proceed by (a) deleting parameters
that live in unconstrained sectors (free, but contentless — this round's
mechanism) or (b) pinning `a_pot,b_pot` to the quoted Lean constant (costly,
now measured directly rather than argued qualitatively).

Effective dimension (Jacobian SVD) is **3**, not 5 and not 6: of the 5
nominal free parameters, `mu_sym` (screening) and the `a_pot,b_pot` pair
(dark energy, but only 2 independent directions of it — confirmed by round
1's `sv_0,sv_1` loadings) carry essentially all the resolvable signal;
`lambda_sym` and `c4_pta_product` both ride consistently near-null
directions (`lambda_sym` SV~1e-14–3e-14 both rounds; `c4_pta_product`
SV~5e-4 this round). The theory is over-parameterized relative to its own
observable content by roughly 2 parameters beyond what any of this round's
real data can see — but "invisible to current data" and "derivable from
first principles" are different claims, and zero-parameter status needs the
second, not the first. The `mu_sym`/screening sector is the one direction
this loop has NOT yet tried to match against a real dataset (no real
fifth-force/screening bound was fetched either round; still ABSENT) — that
is the highest-value next real-data target, since it is the one
observably-alive sector without a Lean-constant or public-data check yet.

## Files

- `tda_prediction.json` — pre-registration, written before any reduced-model
  sweep/TDA ran.
- `reduced_model.py`, `reduced_sweep.py` -> `reduced_sweep.csv`.
- `reduced_chi2.py` -> `reduced_chi2_report.json`.
- `reduced_jacobian.py` -> `reduced_jacobian_report.json`.
- `reduced_tda.py` -> `reduced_tda_report.json`,
  `reduced_tda_persistence_{model,real_sn,circle,null_uniform_box,null_variance_matched}.{json,png}`.
- `lean_pinned_test.py` -> `lean_pinned_report.json`.
