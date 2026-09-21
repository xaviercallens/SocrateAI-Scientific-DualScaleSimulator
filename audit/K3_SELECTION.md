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
3. A statement of what the choice would forbid. Under the current construction, changing `N` moves
   `FRICKE_Y` and nothing else observable, because the potential is a hand-built well — so **no
   observable currently depends on the level at all**, and no level can be falsified. Fixing that is
   prior to arguing for any particular `N`.
