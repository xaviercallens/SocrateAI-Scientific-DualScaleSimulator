# Reverse pass v2: LeanMaster (K3×T²) theorems → new predictions → blind computation

Branch `loop/k3t2-rigidity`, worktree
`/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2`. Python:
`/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python`
(exact `fractions.Fraction`/`int` arithmetic throughout). Tier B at best; nothing here is
"proved". LeanMaster read at commits `4c4ad233`…`41d65194` (it advanced twice during this
session; the theorem statements quoted below — `DualScaleMoonshine/{QSeries,Twining,
TwiningAll,Characters,CharactersAll,Decompositions}.lean`, `DualScaleDyons/{DMVV,Immortal,
ImmortalHigher}.lean` — did not change across those commits), read-only, no LeanMaster
files touched.

**Scope note**: the task asked for 5–8 theorems; this pass completed **3 items in depth**
(plus one item shares its engine with a second Lean theorem, listed separately below) —
each required building and validating an independent computation engine from the primary-
source formulas (Cheng–Duncan–Harvey, Dabholkar–Murthy–Zagier) rather than transliterating
Lean code, and each carries a working negative control. Depth was prioritised over count.
Candidates considered but **not** attempted here: umbral exactness beyond Lean's range
(would reuse the same theta engine but was not reached), and the Tripathy–Trivedi tadpole
cross-check (already resolved in `comparison.json`'s three orientifold DISAGREE rows by
LeanMaster's own `adc85e7`/`v3.21.0`; re-testing a settled row was judged the lowest-value
use of remaining budget).

Excluded (already done by v1, per the task): `p24` `k=6..8`, DMZ (5.16) `m=4`, polar
cancellation `m=4,5`, the dual-scale trace bound `d=7..10`.

## Method

For each item: (1) read the Lean theorem(s) and their stated range; (2) build an
independent Python engine from the same primary-source defining formulas; (3) **validate**
the engine by bit-for-bit reproduction of Lean's own checked range (a gate, not the
deliverable); (4) extend to a new range or a new Hecke-like-operator hypothesis; (5) run a
negative control that must fail. Every script states its exact run command and the exact
regeneration command for any input data file it reads.

| # | Lean theorem(s) extended | Lean's checked range | New content tested | Result |
|---|---|---|---|---|
| 1 | `c_depends_only_on_D` (`DMVV.lean`) | `n ≤ 9` (`zK3 := ellipticGenus 9`) | `n = 10..20` | **HOLDS** |
| 2 | `Twining.twined_div24`, `TwiningAll.twinedAll_div`, `Decompositions.moonshine_modules_decompose` (`DualScaleMoonshine`) | `n ≤ 9` (all three) | `n = 10..20` | **HOLDS** (new fact: CDH's own Table 48 stops at `n=9`) |
| 3 | `ImmortalHigher.lean`'s `m=2,3` pattern (`dmz_911_verified`/`dmz_913_verified_m3`) — **Lean's own gap**: "Not proved: m ≥ 4; the general statement (9.13), DMZ §10" | `m ∈ {2,3}` (both prime) | `m = 4` (**not** prime — the genuinely new structural question), with the correct Hecke-like operator identified from DMZ §4.4/§9.2 | **HOLDS**, with two working negative controls |

No counterexample was found in any item; all three "holds" are the outcome of a
genuinely blind computation against a target built from Tier-L literature formulas
containing no numerical literal from LeanMaster's proof, checked against negative controls
that fail.

## Item 1 — K3 elliptic genus index-1 Jacobi property, `q^{10}..q^{20}`

Script: `item1_index1_q20.py` → `item1_index1_q20_results.json`. Run:
`python item1_index1_q20.py`. Reads `../B-dyons/theta_forms_cache.json` (QMAX=40,
regenerate with `cd ../B-dyons && python theta_forms.py 40`), an independent theta-
function engine (`A-genus/thetas.py`, `B-dyons/theta_forms.py`) never derived from Lean.

Lean's `c_depends_only_on_D` checks, through `q⁹`, that every coefficient of the computed
`Z_K3 = 2φ₀,₁` depends only on `D = 4n − l²` (the defining index-1 property of a weak
Jacobi form). We re-derived this collapse independently from `B_series` for `n = 10..20`:
**zero conflicts**, 22 new `D` values exercised beyond Lean's range. Negative control:
perturbing one `(n,l)` coefficient at `n=15` by `+1` creates a same-`D` disagreement with
its partner `(n,l')` — the check is not vacuous.

**Proposed Lean statement** (directly reachable: `c_depends_only_on_D`'s own proof is
`by decide` at `N=9`; nothing in its statement or proof depends on `N≤9` specifically,
only the truncation argument does):
```lean
theorem c_depends_only_on_D_ext :
    (List.range 21).all fun n => (List.range 45).all fun i =>
      let l : ℤ := (i : ℤ) - 22
      (zK3.getD n LP.zero).coef l == cK3 (4 * n - l * l) := by decide
```
(⚠ prerequisite, same as v1 found for other `dmvv`-based extensions: `zK3 := ellipticGenus 9`
only reaches `n ≤ 9`; this statement needs `zK3 := ellipticGenus 20` first.)

## Item 2 — Twined-series divisibility and the Mathieu-moonshine module decomposition, `n = 10..20`

Script: `moonshine_engine.py` (the shared independent engine, re-implementing
`QSeries.lean`/`Twining.lean`/`TwiningAll.lean`/`CharactersAll.lean`/`Decompositions.lean`'s
formulas fresh from Cheng–Duncan–Harvey, general truncation `N`) plus
`item2_moonshine_q20.py` → `item2_moonshine_q20_results.json`. Run:
`python item2_moonshine_q20.py`.

**Validation gate** (all pass): the fresh engine reproduces, bit-for-bit at `N=9`, all four
`Twining.lean` tables (`2A,3A,5A,7AB`), all sixteen `TwiningAll.lean` tables
(`2B,3B,4A,4B,4C,6A,6B,8A,10A,11A,12A,12B,14AB,15AB,21AB,23AB`), `twinedAll_div`, and
`Decompositions.moonshine_modules_decompose` against CDH's printed Table 48 for `n=0..9`.
(Two implementation bugs were found and fixed during this validation gate — a `1-q^0`
constant-series edge case and an eta-quotient shift/truncation off-by-`s` slice — exactly
the kind of error this gate exists to catch before trusting any extrapolation.)

**Extension to `n=10..20`** (26 classes × 11 levels = 286 new structural conditions,
sharing 16 divisibility checks per level with `twinedAll_div`): every one of
`⟨T_n,χ_i⟩`'s omega-parts (`ω₇,ω₁₅,ω₂₃`) vanishes, the rational part is divisible by
`|M₂₄| = 244823040`, and the resulting multiplicity is `≥ 0`, for every class and every
`n=10..20`. Since Cheng–Duncan–Harvey's own Table 48 stops at `n=9`, **this is a genuinely
new fact** (the first 11 new graded pieces of the Mathieu-moonshine module are consistent
`M₂₄`-representations), not a check against a known table. All sixteen
`TwiningAll`-style divisibility-by-`24·D` conditions (plus the four `Twining`-style
divisibility-by-24 conditions) also continue to hold at every new order.

Negative controls (both have teeth): (a) perturbing the `F_2A` combination's `-16` to `-15`
breaks divisibility-by-24 at several `n∈[10,20]`; (b) perturbing one entry of the
character table (`χ₂` at `1A`, `+1`) breaks integrality/non-negativity at several
`n∈[10,20]`.

**Proposed Lean statement**:
```lean
theorem moonshine_modules_decompose_ext :
    (List.range 21).all (fun n => (List.range 26).all fun i =>
      let v := moonshineInner n i
      v.2.1 == 0 && v.2.2.1 == 0 && v.2.2.2 == 0 && 244823040 ∣ v.1 && 0 ≤ v.1 / 244823040) = true := by
  decide +kernel
```
(No `ellipticGenus`/`dmvv` reachability prerequisite here — unlike `dmvv`-based statements,
`twined24`/`etaQ`/`inner4` are not capped by any fixed truncation constant in Lean; `N=20`
just needs to be substituted for the `9` hard-coded in each definition's call site.)

## Item 3 — Immortal dyons at `m = 4`: the correct Hecke-like operator (non-prime `m`)

Script: `item3_immortal_m4_hecke.py` → `item3_immortal_m4_hecke_results.json`. Run:
`python item3_immortal_m4_hecke.py`. Reads `../B-dyons/theta_forms_cache.json` (QMAX=40)
and `../B-dyons/hurwitz_class_numbers_results.json` (DMAX=40); regenerate with
`cd ../B-dyons && python theta_forms.py 40 && python hurwitz_class_numbers.py`.

**The question the task posed**: LeanMaster's `ImmortalHigher.lean` proves the DMZ
(9.11)/(9.13) identity `(ψ₀,ₘ₊₁^opt − 12^{m+1}·A·A₂,ₘ)/A = 12^{m+1}·Φ₂,ₘ^opt` for `m=2,3`
(both **prime**), using the naive two-term operator `H|Vₘ`
(`hurwitzVSer`: `12H(∆)+m·12H(∆/m²)·[m²∣∆]`) and its docstring explicitly flags
`m ≥ 4` as "Not proved". **What should the operator be for `m=4` (`=2²`, not prime)?**
DMZ §4.4 eq. (4.37) defines the true Hecke-like `V_{k,t}`:
`c(φ|V_{k,t};n,r) = Σ_{d∣gcd(n,r,t)} d^{k−1}·c(φ; nt/d², r/d)`. For prime `m` the divisors
of `gcd(n,r,m)` are only `{1,m}` — exactly Lean's two-term formula. **For `m=4` the
divisors can also be `{1,2,4}`**, i.e. there is a middle `d=2` term the naive "prime
pattern" misses entirely. DMZ's own text (papers/foundations/1208_4074.txt, ~l.3851)
states this explicitly and gives the fix: `Φ₂,₄^opt = −H∣V₄ + 2·H∣U₂`, with `U₂` the
index-raising operator `c(n,r) ↦ c(n,r/2)` (eq. 4.36). Table 2 of the same paper gives the
companion weak Jacobi form `ψ₀,₅^opt = B⁵ − 10E₄A²B³ + 20E₆A³B² − 15E₄²A⁴B + 4E₄E₆A⁵`.

**What was computed**: `G₅ = [p⁵]` of the DMVV product (from the K3 elliptic genus alone,
independent of the `ψ₀,₅^opt`/`A₂,₄` ansatz; validated first against the already-Tier-A
`p₂₄(5) = 176256` at `q⁰,y=1`), and `(ψ₀,₅^opt − 12⁵·A·A₂,₄)/A` (built from `A,B,E₄,E₆` via
the theta-product engine and the generic `A₂,ₘ` construction of `part3_polar.py`,
independent of any Hecke-operator choice), compared against the **full-divisor-sum**
`Φ₂,₄^opt`.

**Result: HOLDS**, exactly, over `n=0..5`, `|l|≤15` (68 keys, zero mismatches, over the
window `n ≤ QCHK−1` with `QCHK=6`, margin `QMAX=40 ≥ 6×QCHK`).

**Two negative controls, both with teeth**:
* the naive two-term "prime pattern" (`d∈{1,4}` only, skipping `d=2`) **fails** — 10
  examples found in the tested window where the `d=2` term is genuinely nonzero and the
  two formulas differ (e.g. `(n,r)=(2,0)`: full `= 5`, naive `= 1`).
* DMZ's own alternative eq. (9.14), `Φ₂,ₘ^opt = −2·H∣V^{(1)}_{2,m}` for prime powers,
  expands at `m=4` (via (4.38), `V^{(1)}_{2,4}=V_{2,4}−2U₂`) to `−2H∣V4+4H∣U2` — **twice**
  the `m=4`-specific formula quoted two paragraphs earlier in the same source. This is a
  genuine internal tension in the cited text (not a transcription error on our part: both
  formulas are quoted verbatim); we tested both, honestly, rather than picking one by
  fiat, and report that the computation sides with the `m=4`-specific formula
  (`Φ₂,₄^opt = −H∣V4+2H∣U2`), not the (9.14) factor-two variant, which also fails to
  match.

(Sanity: comparing the *raw* `G₅/A`, rather than the `ψ₀,₅^opt`-based finite part, against
the same target correctly **fails** — `G₅/A` is `∆ψ₄`, not the "optimal" mock Jacobi form;
matching Φ₂,₄^opt requires the specific weak-Jacobi-form subtraction Table 2 encodes, not
just any finite part. This is the expected, honest non-match, analogous to why Lean's own
`immortal_m2`/`immortal_m3` — the raw-`∆ψ_m` identities — carry *extra* `E₄AB^k`/`E₆A²B^k`
terms that the `ψ^opt`-based `dmz_911_verified`/`dmz_913_verified_m3` do not.)

**Proposed Lean statement** (⚠ prerequisites: `hurwitzVSer` is hard-coded to the two-term
prime pattern and would need a genuine divisor-sum generalisation `hurwitzV4Full`;
`immortalNumM`/`dmvv 4 2` reach only `m+1≤4`, so `m=4` needs `dmvv 6 2`, i.e.
`zK3 := ellipticGenus N`, `N≥12`, per v1's already-identified reachability gap — not
attempted here since this worktree does not touch LeanMaster):
```lean
def hurwitzV4Full (n r : ℤ) : ℤ :=  -- 12*(H|V4), full divisor sum d | gcd(n,r,4)
  (divisorsGcd n r 4).foldl (fun acc d => acc + d * h12 ((4*n - r^2) / (d*d))) 0
def hurwitzU2 (n r : ℤ) : ℤ := if r % 2 = 0 then h12 (4*n - (r/2)^2) else 0  -- times 12 already via h12
theorem dmz_913_hecke_m4 :
    Ser.eqB (divByA (Ser.sub psiOpt5 (Ser.scale 248832 (aTimesA2m 4 2))) 2)
      (Ser.sub (Ser.scale (-1) (hurwitzV4FullSer 4 2)) (Ser.scale (-2) (hurwitzU2Ser 2))) = true := by
  decide +kernel  -- needs the two def's above and psiOpt5 first
theorem dmz_913_hecke_m4_naive_fails :
    ¬ Ser.eqB (divByA (Ser.sub psiOpt5 (Ser.scale 248832 (aTimesA2m 4 2))) 2)
      (Ser.sub (Ser.scale (-1) (hurwitzVSer 4 2 true)) (Ser.scale (-2) (hurwitzU2Ser 2))) := by
  decide +kernel
```

## Summary

Three items, each requiring a fresh engine validated against Lean's own checked range
before extrapolating: (1) the K3 elliptic genus's defining index-1 property extends
cleanly to `q^{20}`; (2) the Mathieu-moonshine module's graded pieces at `n=10..20` are
genuinely new, structurally-forced facts (integer, non-negative, `M₂₄`-consistent
multiplicities) with two working negative controls sharing the engine with
`twinedAll_div`'s extension; (3) the specific, literature-identified generalisation of the
Hecke-like operator to the first non-prime case `m=4` — `Φ₂,₄^opt = −H∣V4+2H∣U2`, the full
divisor-sum operator, not the naive prime-pattern extrapolation Lean's own docstring
flagged as open — holds exactly against an independently-computed target, with the naive
generalisation and one of the paper's own alternative formulas both failing as required.
Two Lean-side reachability prerequisites (raising `zK3 := ellipticGenus 9` to `N≥12..20`,
and generalising `hurwitzVSer` to a true divisor sum) were identified and are reported for
LeanMaster's maintainers, not attempted here since this worktree does not touch LeanMaster.
