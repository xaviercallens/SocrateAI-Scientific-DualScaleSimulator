# Reverse pass: LeanMaster (K3×T2) statements → predictions → blind experiment

Branch `loop/k3t2-rigidity`, worktree
`/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2`. Python:
`/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python`
(exact `fractions.Fraction`/`int` arithmetic throughout; `sympy.Rational` for
one item). Tier B at best — nothing here is "proved". LeanMaster read at
`7be7626` (`DualScaleDyons/DMVV.lean`, `DualScaleDyons/ImmortalHigher.lean`,
read-only, no LeanMaster files touched).

## Method

Five of LeanMaster's own theorems are stated only "through `q^N`" or for a
small finite list of `m`/`k`. Each is treated as a **hypothesis**: build an
*independent* Python engine for the same defining formulas (a fresh
`(n,l)`-Laurent-series engine, `nlseries.py`/`dmvv_product.py`, not a copy of
Lean's code or of the sibling `dyons/series.py` `P,Y`-power engine), feed it
our own independently-computed theta data
(`dyons/theta_forms_cache.json`, θ-functions to `QMAX=22` ⇒ `c(D)` for
`D ≤ 88`, already cross-checked AGREE against Lean in the forward pass), and
**first reproduce every already-Tier‑A Lean result** before trusting the
engine to extrapolate beyond Lean's checked range. Every extrapolation below
passed that validation step exactly (bit-for-bit match on the Lean-checked
sub-range, see each script's `validation_*` block).

| # | Lean theorem(s) extended | Lean's checked range | New range tested | Result |
|---|---|---|---|---|
| P1 | `p24_values`, `goettsche` (`DMVV.lean`) | `k ≤ 5` | `k = 6,7,8` | **HOLDS** |
| P2 | `dmz_516_q1`/`dmz_516_q2` (`DMVV.lean`), `m=4` row | `q^1` only | `q^2, q^3, q^4` | **HOLDS** |
| P3 | `polar_part_removes_pole`/`polar_coefficient_pinned` (`DMVV.lean`), `immortal_exact_m23` (`ImmortalHigher.lean`) | `m ∈ {1,2,3}` | `m = 4` | **HOLDS** (with caveat, see below) |
| P4 | same theorems | `m ∈ {1,2,3}` | `m = 5` | **HOLDS** (same caveat) |
| P5 | blind Track C `lattices-duality/04_dual_scale_bound.py` (not a located Lean statement — see note) | `d = 1..6` | `d = 7..10` | **HOLDS** |

No counterexample was found. All five predictions are confirmations, not
failures — see the honesty caveats below for exactly how strong each
confirmation is.

## P1 — Goettsche numbers `p24(6), p24(7), p24(8)`

Script: `p1_p2_dmz516_m4.py` → `p1_p2_results.json` (`P1_goettsche_k6_7_8`).
Run: `python p1_p2_dmz516_m4.py`.

```
p24(6) = 1073720
p24(7) = 5930496
p24(8) = 30178575
```

**Honesty correction**: this is *not* two independent numerical routes. The
`q^0,y=1` restriction of the DMVV product provably (not just numerically)
collapses at `Q=0` to `∏_r(1-p^r)^{-(c(0)+2c(-1))}` — the *same* convolution
as the direct `∏(1-q^k)^{-24}` computation, once the exponent is fixed at
`c(0)+2c(-1)=24`. So this is **one route plus a derivation that the exponent
is 24**, cross-checked by a **negative control**: exponent 23 or 25 give a
visibly different table (`negative_control_exponent_perturbation` in the
JSON), confirming `c(0)=20, c(-1)=2` (Lean's `c_first`) are load-bearing.

**Proposed Lean statement** (reachable at Lean's current `zK3 := ellipticGenus 9`,
since `dmvvDs 8 1` needs `D ≤ 32 ≤ 36`, no truncation bump required):
```lean
theorem p24_next : (List.range 3).map (fun i => p24 (6 + i)) =
    [1073720, 5930496, 30178575] := by decide

theorem goettsche_ext :
    (dmvv 8 1).map (fun Gk => Gk.map LP.eval1) =
      (List.range 9).map fun k => [p24 k, 0] := by decide +kernel
```

## P2 — DMZ (5.16), `m = 4` row, beyond Lean's `q^1`-only coverage

Script: `p1_p2_dmz516_m4.py` → `p1_p2_results.json` (`P2_dmz516_m4`).

Lean's docstring quotes (Tier L, from DMZ arXiv:1208.4074, via
`DualScaleDyons/DMVV.lean`):
`72 Δψ₄ = 51 A⁻¹B⁵ + 155 E₄AB³ + 93 E₆A²B² + 102 E₄²A³B + 31 E₄E₆A⁴`.
`dmz_516_q1` checks the cleared form (`72 G₅ = 51B⁵+155E₄A²B³+93E₆A³B²+102E₄²A⁴B+31E₄E₆A⁵`)
through `q^1`; `dmz_516_q2` explicitly **excludes** `k=5` (its range is
`k ≤ 4`). We computed `G₅` from an independently re-implemented DMVV product
(validated against Lean's `dmz_516_q1`/`dmz_516_q2` for `k ≤ 4` at `q^1,q^2`,
and against the earlier blind Track's hand-derived `G₂` —
`dyons/part2_psi_m_results.json`'s `sample_G2_terms` — bit-for-bit) and
against `A,B,E4,E6` from `dyons/theta_forms_cache.json`.

Result: the identity **holds through `q^1, q^2, q^3, q^4`**. Negative
control: perturbing any one of the five coefficients `{51,155,93,102,31}` by
`+1` breaks the identity at 39 coefficients (mirrors
`dyons/part2_psi_m.py`'s `c(D)+1` control).

**Provenance, stated plainly**: the RHS *ansatz* itself is Tier L (quoted
from the paper); what is Tier B here is that our independently-computed
`G₅` (from the K3-elliptic-genus/DMVV product alone, not from the ansatz)
matches it through four more orders in `q` than Lean has checked, and does
so non-tautologically (negative control fails).

**Proposed Lean statement** (⚠ prerequisite: `dmvv K Q` reads `c(D)` up to
`D = 4KQ`; Lean's `zK3 := ellipticGenus 9` only reaches `D ≤ 36`, so `K=5,Q=2`
needs `D ≤ 40` ⇒ **`zK3` must first be raised to `ellipticGenus N`, `N ≥ 10`**,
and the `dmvv_reachable`-style guard re-checked at `(5,2)` before this
`decide +kernel` is attempted; `Q=3,4` need `N ≥ 15, 20` respectively — cost
of `decide +kernel` on `dmvv 5 4` will be substantially higher than on
`dmvv 4 2`):
```lean
theorem dmz_516_q2_m4 :
    Ser.eqB (Ser.scale 72 ((dmvv 5 2).getD 5 []))
      (Ser.add (Ser.add (Ser.scale 51 (mulAll [phiB 2,phiB 2,phiB 2,phiB 2,phiB 2] 2))
                        (Ser.scale 155 (mulAll [e4 2,phiA 2,phiA 2,phiB 2,phiB 2,phiB 2] 2)))
        (Ser.add (Ser.add (Ser.scale 93 (mulAll [e6 2,phiA 2,phiA 2,phiA 2,phiB 2,phiB 2] 2))
                          (Ser.scale 102 (mulAll [e4 2,e4 2,phiA 2,phiA 2,phiA 2,phiA 2,phiB 2] 2)))
                (Ser.scale 31 (mulAll [e4 2,e6 2,phiA 2,phiA 2,phiA 2,phiA 2,phiA 2] 2))))
      = true := by decide +kernel
```

## P3, P4 — Immortal-dyon polar cancellation at `m = 4, 5`

Script: `p3_p4_polar_m4_m5.py` → `p3_p4_results.json`.

`polarDefect m coef := G_{m+1} - coef · A · A₂,ₘ` (same object under three
names: `polarDefect`/`polar_part_removes_pole`/`polar_coefficient_pinned` in
`DMVV.lean`, `immortalNumM`/`immortal_exact_m23` in `ImmortalHigher.lean`).
`vanish2` requires `eval1 = 0` (the pole coefficient) **and** `d1 = 0` (its
`y`-derivative at `y=1`). Lean checks `m ∈ {1,2,3}`; "Not proved: `m ≥ 4`".

We reproduced Lean's `m=1,2,3` pattern exactly (`validation_m1_2_3`, all
`matches_lean_pattern: true`) with an independently re-implemented `A₂,ₘ`
(DMZ (9.55) strip expansion) and `rPart` (both fresh code from the defining
formulas), then extended:

* `m = 4`: `vanish2` holds at `coef = p24(5) = 176256`; fails at `176255` and
  `176257` (negative control), exactly as predicted.
* `m = 5`: `vanish2` holds at `coef = p24(6) = 1073720` (a value from P1, not
  from Lean); fails at `1073719` and `1073721`.

**Honesty caveat (carried over from Lean's own docstring, not discovered by
us)**: `d1` was **identically zero** in every case tested here, including
both negative controls — `A`, `A₂,ₘ` and the `dmvv` product are all
`y ↔ y⁻¹`-symmetric, so `d1` is an odd functional of an even series and never
discriminates anything in this test. What actually distinguishes the
positive from the negative controls is `eval1` alone: the single first-order
statement `G_{m+1}(q^n, y=1) = p24(m+1)` for `n ≤ Q`. This is **not** an
independently-verified double-pole cancellation at `m=4,5` — it is the
weaker, still non-trivial, pole-coefficient identity. (Lean's own text under
`polar_part_removes_pole` says the same about `m ≤ 3`.)

**Proposed Lean statement** (⚠ prerequisite: current `polarDefect`/
`immortalNumM` are hard-coded to `dmvv 4 2`, i.e. `K=4`, which only reaches
`m+1 ≤ 4`; extending to `m=5` needs `dmvv 6 2`, reading `c(D)` up to
`D = 48`, so **`zK3` must be raised to `ellipticGenus N`, `N ≥ 12`**):
```lean
theorem polar_part_removes_pole_ext :
    [4, 5].all (fun m => vanish2 (polarDefect m (p24 (m + 1)))) = true := by
  decide +kernel

theorem polar_coefficient_pinned_ext :
    [4, 5].all (fun m => !vanish2 (polarDefect m (p24 (m + 1) + 1)) &&
      !vanish2 (polarDefect m (p24 (m + 1) - 1))) = true := by decide +kernel
```
(Both statements would need `polarDefect`/`aTimesA2m`/`immortalNumM` in
`DMVV.lean`/`ImmortalHigher.lean` re-parametrized from the hard-coded `4` to
a `K` argument, plus the reachability bump above, before they typecheck as
written; not attempted here, since we do not touch LeanMaster.)

We deliberately did **not** extend the Hecke-class-number (`H|V_m`) ansatz
itself (`dmz_911_table`/`dmz_912_table`/`psiOpt3`/`psiOpt4`,
`immortal_m2`/`immortal_m3`) to `m=4,5`: Lean's own file states the general
`m`-formula is only Tier L ("Not proved: ... the general statement (9.13),
DMZ §10"), and it is given explicitly only for `m=2,3`; guessing the degree-6
(`m=4`) or degree-7 (`m=5`) polynomial ansatz in `A,B,E4,E6` without the
paper's general formula would be an unfounded guess, not a derivation from a
defining formula, so we left it out rather than risk reporting a
speculative match/mismatch as if it were principled.

## P5 — Dual-scale trace bound, `d = 7..10`

Script: `p5_trace_bound_d7_10.py` → `p5_trace_bound_results.json`.

Extends `lattices-duality/04_dual_scale_bound.py` (Track C of the forward
pass; 2000 exact-`sympy.Rational` SPD samples per dimension, `d=1..6`) to
`d=7,8,9,10`. The bound `tr(G) + tr(G⁻¹) - 2d ≥ 0` holds on all 8000 new
samples (`all_f_nonnegative: true` for every `d`), with exact positive
minima observed (e.g. `d=7`: `42128270630/271025753`); the perturbation scan
around `G=I` again shows `f > 0` for every `ε ≠ 0` sampled.

Algebraically this is `Σᵢ (xᵢ-1)²/xᵢ ≥ 0` for the `d` eigenvalues, so
`d`-independence is structural, not a coincidence — this extension is
expected to hold for *every* `d`, and did. We did not locate a LeanMaster
theorem literally named for this bound in a quick read-only search of
`StringTheoryFoundation/`/`DoubleFieldTheory/`; if LeanMaster does state it
for general `d` (as the `string-theory-foundation` skill's description
suggests), this Track-C numerical extension is redundant with an already-
general Lean theorem, and we report that uncertainty honestly rather than
claim a specific Lean range was extended.

## Summary

Two genuinely new numerical facts were produced (`p24(6)=1073720`,
`p24(7)=5930496`, `p24(8)=30178575`, and the confirmation that DMZ (5.16)'s
`m=4` row and the polar-cancellation pattern extend at least to `m=5` and
`q^4`); zero counterexamples were found. Every extrapolation was gated by
first reproducing the corresponding Lean-checked sub-range bit-for-bit with
the same independent engine, and every "holds" claim above carries the
honesty caveats stated in its section (P1: one route, not two; P3/P4: `d1`
never discriminates, only `eval1` does). Two Lean-side prerequisites were
identified and are *not* satisfied by the current LeanMaster tag
(`zK3 := ellipticGenus 9` must be raised to `N ≥ 10..20` before the proposed
`decide +kernel` statements above would even typecheck against real data)
— reported for LeanMaster's own maintainers, not attempted here since this
worktree does not touch LeanMaster or `proofs/`.
