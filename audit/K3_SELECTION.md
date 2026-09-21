# Which K3? — what this repository has actually committed to, and on what evidence

**Written 2026-09-21.** Audit note. **No prediction is registered by this document** (`PRE_REGISTRATION.md`
A5). Nothing here is a physical claim: Stream 1 records program-wide that no exact physical observable
exists anywhere in the programme. This is about **which mathematical family**, not about our universe.

Pins: Stream 1 `bb74acb56f386a97e433f94eb0b2632ed03bc4ca`; LeanMaster `ede49f06800cde8177867c02bb08d44f7be275c5`
(v3.44.0). Evidence for the code facts is in `STREAM1_BRIDGE.md`.

---

## 1. The finding: this repository already committed to a level, and nobody chose it

`FRICKE_Y = 1.0 / math.sqrt(12.0)` is the Fricke point of the implemented potential — verified, not
assumed: the potential is stationary there (`dV/dx = dV/dy = 0`) and not at `τ = i` (`dV/dy = 1.423`),
see S1-F1. By Stream 1's `root_orthogonal_iff_selfdual`, the self-dual locus is `N·τ² = −1`, so
`τ = i/√N`. Measured, `1/FRICKE_Y² = 11.999999999999996`.

**So the model is at level `N = 12`.** Under Dolgachev's framework that means the `M₁₂`-polarized
family, with transcendental lattice

```
T(X) = U ⊕ ⟨24⟩,   rank 3,   signature (2,1),   det = −24
```

and coarse moduli space the Fricke modular curve `X₀(12)⁺` — **one modulus** (Dolgachev 1996 Thm 7.1,
Tier L literature; see `PRE_REGISTRATION.md` A10.1).

**Nobody derived this.** `FRICKE_Y` is a hard-coded constant, and the potential it sits in is a
hand-built double well,

```python
v = a_pot*cos2 + b_pot*sin2 + cos2*(y - FRICKE_Y)**2 + sin2*(y - ORBIFOLD_Y)**2
```

whose two stationary points are placed there **by construction**. It has no modular symmetry at all —
measured: invariant under `T : τ ↦ τ+1` to `6.7e-16`, and *not* under `S` (`8.77`) or under the
level-12 Fricke involution (`35.81`), S1-F2.

Stream 8 states the standard this programme holds itself to (`docs/STREAM8_WHICH_K3.md` §0, verbatim):

> **which mechanism of the theory chooses it**, since a hand-picked point is a free parameter like
> `α'` was. A choice must be forced, not fitted.

**By that standard, level 12 in this repository is fitted, not forced.** That is the finding. It does
not make the level wrong; it makes it an input, and it should be counted as one.

## 2. The programme is running two different levels, and no document reconciles them

| | level | transcendental lattice | det | moduli space |
|---|---|---|---|---|
| **Stream 1** (`Agora/Geometry/MnLattice.lean`, `T7 := TN 7`, `T7_det = −14`) | `N = 7` | `U ⊕ ⟨14⟩` | `−14` | `X₀(7)⁺` |
| **This repository** (from `FRICKE_Y`) | `N = 12` | `U ⊕ ⟨24⟩` | `−24` | `X₀(12)⁺` |

Stream 1's level 7 is not arbitrary — it is tied to Cooper's `s₇` sequence, the Hauptmodul
`h = (η(7τ)/η(τ))⁴`, the rational map `z = h/(1+13h+49h²)`, the Fricke fixed points `h = ±1/7`, and the
prefactor `P₂ = 1 − 26z − 27z²` whose roots `{−1, 1/27}` are the images of those fixed points
(`zOf_fricke`, `fricke_fixed_iff`, `s7_P2_discriminant`, `s7_singular_points_are_selfdual`). There is a
whole arithmetic story at 7.

At 12 there is a constant in a double well.

**Open, and recorded as open:** either the level-12 constant should be replaced by the level the
mathematics actually supports, or a reason for 12 should be produced. Neither is done here. Note also
that the lattices themselves are not isometric — `det −14 ≠ det −24` — and by Stream 1's
`no_isometry_G0N_TN` argument (an isometry has determinant `±1` and so preserves the Gram determinant)
no change of basis relates them. **These are different K3 families, not two descriptions of one.**

## 3. The `12` of Stream 8 is NOT this `12`

This is the coincidence to refuse, and it is exactly the kind of conflation LeanMaster amended in
v3.44.0. Memory of the K3×T² work records "the attractive Kummer surface with `D = 12` at `τ = ω`".
That `D` and this `N` are different invariants of different objects:

| | Stream 8's `D` | this note's `N` |
|---|---|---|
| Picard number | `ρ = 20` (attractive / singular) | `ρ = 19` |
| transcendental lattice | **rank 2**, positive definite, `T(X) = [[2a,b],[b,2c]]` | **rank 3**, signature (2,1), `U ⊕ ⟨2N⟩` |
| invariant | discriminant `D = b² − 4ac < 0` | level `N`; lattice determinant `−2N` |
| moduli | a **point** (rigid — no modulus) | a **curve** `X₀(N)⁺` — one modulus |
| at the value 12 | a specific attractive surface | `T(X) = U ⊕ ⟨24⟩`, determinant **−24** |

The determinants settle it: level `N = 12` gives `−24`, not `12`. **There is no relation between the
two numbers, and the agreement of the digits is a coincidence.** An `M_n`-polarized family does contain
attractive (`ρ = 20`) points where the transcendental rank drops 3 → 2, so the two pictures meet — but
meeting is not equality, and no computation here identifies which attractive surfaces sit in `X₀(12)⁺`.

## 4. Appendix — a criterion that does NOT discriminate

Dolgachev's moduli space is a Hauptmodul-bearing curve only when its genus is 0, so one might hope
"genus 0" selects a level. It does not: **both levels in play pass.** Computed by
`scripts/fricke_genus.py` (`--self-test`):

```
   N   g(X0(N))   fix(w_N)   g(X0(N)+)
   7          0          2           0      <- Stream 1
  12          0          2           0      <- this repository
  37          2          2           1      <- first level that FAILS
  71          6         14           0
```

Of `N = 1..300`, **38** levels have `g(X₀(N)⁺) = 0`. The criterion excludes most levels but separates
neither 7 nor 12, so it is recorded as a **negative control**, not a selection principle.

The script's self-test is what makes it quotable. It checks `g(X₀(N))` against textbook values
(including `g(X₀(37)) = 2`, `g(X₀(49)) = 1`, `g(X₀(121)) = 6`); class numbers against 18 known values
(`h(−163) = 1`, `h(−71) = 7`, `h(−44) = 3` — the last requires counting only **primitive** forms, since
`(2,2,6)` has content 2, an error caught by the test); that `g(X₀(N)⁺)` is a non-negative integer for
every `N ≤ 300`, which Riemann–Hurwitz would break otherwise; and — the check the author did not put in
by hand — that the **primes** with `g(X₀(p)⁺) = 0` come out as exactly

```
{2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71}
```

Ogg's supersingular primes, the primes dividing the order of the Monster.

## 5. Status, by tier

| Claim | Tier | Basis |
|---|---|---|
| `FRICKE_Y = 1/√12`, so the implemented model sits at level `N = 12` | **A** (about this code) | measured; `1/FRICKE_Y² = 12.000`, potential stationary there and not at `τ = i` |
| level `N` ⟹ `T(X) = U ⊕ ⟨2N⟩`, det `−2N`; self-dual locus `Nτ² = −1` | **A** | Stream 1, kernel-proved, re-verified in `tests/test_stream1_bridge.py` |
| moduli space of `Mₙ`-polarized K3s is `X₀(n)⁺`, hence one modulus | **L** | Dolgachev 1996 Thm 7.1, literature; Stream 1 proves only the `⊇` inclusion |
| level 12 is **fitted, not forced** | **A** (about this code) | it is a hard-coded constant in a potential with no modular symmetry |
| `g(X₀(N)⁺) = 0` for `N = 7` and `N = 12` | **A** | `scripts/fricke_genus.py`, self-tested against Ogg's primes |
| that any of this selects the K3 of our universe | **C** | nothing supports it; Stream 1 forbids the physical reading |

## 6. What would actually count as forcing a level

Recorded so the next attempt has a bar to clear, **not** registered as a prediction:

1. A level that is **output** by a mechanism already in the programme rather than typed into a
   constant — e.g. a tadpole or flux-quantisation condition that admits `N` only for certain values.
2. A check that the chosen `N` is consistent across the programme: Stream 1's arithmetic is at 7, so
   an argument for 12 must say what happens to `s₇`, `P₂` and the Hauptmodul at 12.
3. A statement of what the choice would forbid. Under the construction this repository actually
   integrates, changing `N` moves `FRICKE_Y` and nothing else observable, because the potential is a
   hand-built well — so **no observable depends on the level at all**, and no level can be falsified.
   Fixing that is prior to arguing for any particular `N`.

   **Update 2026-09-21: the obstruction is now liftable, though not yet lifted.**
   `leanflow/core/modular_potential.py` (spec item L1) is a genuinely `Γ₀(N)⁺`-invariant potential
   built from `F_N(τ) = j(τ) + j(Nτ)`, and it *does* read the level. The potential currently
   integrated cannot: `compute_potential` has no level parameter at all, so `N` reaches it only
   through a module constant.

   **How much of that dependence is real, measured rather than assumed.** Over 14 random domain
   points:

   | | mean `\|V₇ − V₁₂\|` at the same `τ` | after rescaling so `12·Im τ' = 7·Im τ` |
   |---|---|---|
   | `mode="log"` | `0.068` | `0.016` — **76% removed** |
   | `mode="ratio"` | `0.043` | `0.041` — essentially unchanged |

   In `log` mode `V` is dominated by `j(Nτ) ~ e^{2πNy}`, so it is close to a function of the product
   `N·y`: **most of its apparent level-dependence is a stretched imaginary axis, not structure.**
   About a quarter survives. In `ratio` mode the rescaling removes almost nothing, so that dependence
   *is* structural — it sits in where the wells are and what shape they have. Read the `log` number as
   the weaker of the two.

   So the level can now in principle have consequences. It does not yet, because the modular potential
   is **opt-in and not adopted** — adopting it would change every reported trajectory, a decision
   deliberately left open while the Zenodo draft is under review. Until then §1 stands unchanged:
   level 12 is fitted, and nothing in this repository could falsify it.

---

## 7. Addendum 2026-09-21b — reconciled with Stream 2's independent register

Stream 2 (`SocrateAI-Scientific-Agora-K3-DarkMatter`, commit `da84a90`) published
`briefs/THOUGHT_EXPERIMENTS_K3_SELECTION_2026_09_21.md`, a GE-7…GE-15 register asking the same
question — "which K3?", seen from the modular curve. Read as a source; nothing of theirs was rebuilt.

### 7.1 They reached the same negative conclusion by a different route

Their **GE-8** (*"la poussière d'étoiles"*) shows that **`ρ = 20` alone selects nothing**: CM points are
dense on the modular curve, so requiring maximal Picard number cuts out no piece of it. What makes any
list finite is a **bound on `|D|`** — and that bound, not the cut, does the work.

This is the same shape as §4 of this note: my genus-0 criterion excludes most levels and separates
neither 7 nor 12, so it is a negative control. **Two independent criteria, from two streams, each
selects nothing.** Their commit message states it plainly: *"nothing selects, nothing executed."*

Their **GE-9** is worth recording against a temptation: the three singular points of the `s₇` operator
are `ρ = 20` points (`D = −28, −7, −3`), which looks like a discovery and **is forced** — a fixed point
of a finite-order integer matrix satisfies an integer quadratic equation, hence is automatically CM.
Consistent with Stream 1's `s7_singular_points_are_selfdual`. One of their three was *not* predicted by
the prior hand estimate, and they record that clause as **REFUTED** rather than quietly absorbing it.

### 7.2 What this repository can contribute: nobody had run the criterion at level 12

Their certificate `data/certificates/A2_MEMBERSHIP.json` (checker `f091e43`) gives the exact occurrence
criterion, **Tier E** on their side:

> `D` occurs in the level-`n` family **iff** `D` is a square modulo `4n`.

They computed the `n = 7` and `n = 10` columns; **their certificate contains no `n = 12` column**, and
this repository is at `n = 12` (§1). (That is what was checked — not that no other stream has computed
it somewhere else.) `scripts/level_membership.py` transcribes their criterion and applies it;
`--self-test` **reproduces both of their columns for all 39 discriminants in their certificate**, so
the transcription is validated against their own data before being extended.

```
  level n= 7: 27 of 50 discriminants |D| <= 100
  level n=10: 22 of 50
  level n=12: 16 of 50            <- this repository; the most restrictive of the three

  Stream 8's two "most attractive" K3s (the self-dual torus points):
    tau = omega (D = -3): n7=True   n10=False  n12=False
    tau = i     (D = -4): n7=False  n10=True   n12=False
```

### 7.3 The finding: level 12 contains neither surface that E2 names

LeanMaster Stream 8's **E2** proposes, as its *candidate answer* to "which K3":

> **E2 — "The self-dual points are the most attractive K3s."** … Their forms `(1,0,1)`, `(1,1,1)` are
> the two most attractive K3s, `D = 4, 3` … *Candidate answer:* the dual-scale K3 is the attractive K3
> over the self-dual torus.

**At level 12, neither `D = −3` nor `D = −4` is admitted.** So the level this repository is implicitly
committed to **contains neither of the two surfaces that Stream 8's candidate answer names.** Level 7
admits `D = −3` and level 10 admits `D = −4`; level 12 admits neither.

**What this is: arithmetic.** Which CM discriminants lie on `X₀(12)⁺`, computed and validated against
Stream 2's certificate.

**What it is not: a contradiction between two answers.** §3 of this note separates the two objects and
that separation applies here too. E2 asks *which attractive (`ρ = 20`, rigid) surface*; a level asks
*which `ρ = 19` polarized family*, a one-dimensional moduli space. A family not containing a particular
rigid point is not a competing proposition to a claim about that point — the two are not answers to the
same question, and it would be the same conflation §3 refuses to say they cannot both hold.

So this note does **not** say level 12 is wrong, does not rank 12 against `D = 3, 4`, and does not put
any level in or out of contention — nothing selects, as both this note (§4) and Stream 2's GE-8
establish independently.

What it does is **constrain what a level-12 story would have to give up**: it could not also claim the
attractive surfaces over the self-dual torus points, because they are not in that family. Levels 7 and
10 each keep one of them. That is a cost, stated, and it is the first statement in this repository
about the level that could have come out otherwise.

The separate argument that level 12 is **fitted rather than forced** stands on §1 alone — a hard-coded
constant in a potential with no modular symmetry — and does not need this computation to support it.

### 7.4 Recorded, not registered

Per A5 no prediction is registered. Tiers: the criterion is Stream 2's Tier E, transcribed and
validated here; "level 12 admits neither `D = −3` nor `D = −4`" is Tier A about the arithmetic;
E2 is Stream 8's own **Tier C** candidate, not a theorem; and any reading of any of this as a statement
about our universe remains Tier C and forbidden by Stream 1's program-wide note.

**Open question handed back, not answered here:** Stream 2's GE-8 asks *who bounds `|D|`?* and routes it
to their GE-13. This repository has nothing to offer on that until spec item L5 gives the level an
observable consequence at all (`specs/LEANFLOW_ARCHITECTURE.md`).
