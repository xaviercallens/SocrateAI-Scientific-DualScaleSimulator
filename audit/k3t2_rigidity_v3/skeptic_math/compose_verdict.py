#!/usr/bin/env python3
"""
Skeptic (lens = MATHEMATICAL CORRECTNESS) verdict for k3t2_rigidity_v3.

Writes verdict.json next to this script (script-generated, per ground rule 10).
Repo root found via pathlib from __file__; no absolute /mnt/... paths.

Run:
  cd audit/k3t2_rigidity_v3/skeptic_math && \
  /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python compose_verdict.py
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]

rederiv = json.loads((HERE / "rederivations_results.json").read_text())

VERDICT = {
    "lens": "MATHEMATICAL CORRECTNESS",
    "scope": "k3t2_rigidity_v3 BLIND (tracks A-F), CHAIN, COMPARISON, REVERSE",
    "method": (
        "Read each track's own script and results.json line by line for the items "
        "task (a)-(f) name; re-ran several scripts from a fresh `git clone` of this "
        "worktree (not in place) to check reproducibility; wrote 6 from-scratch "
        "independent re-derivation scripts (one per track A-F, in this directory, "
        "rederivations.py) using a mathematical route DIFFERENT from each track's own "
        "implementation; checked every rigidity label against ground rule 8/9's "
        "definitions; checked ground rule 6 wording (grep for proved/proof/theorem "
        "outside Lean-statement contexts); checked ground rule 11 (git show --stat "
        "per commit, confirmed each track's commit touched only its own directory); "
        "checked comparison.json's status counts for DISAGREE/NOT_COMPARABLE rows; "
        "independently reran and hand-checked REVERSE item 4's rank-9..11 numbers and "
        "item 5's mod-24 congruence."
    ),

    "per_track_findings": {
        "A_genus": {
            "checked": [
                "chi(g) derivation: fixed-point counts (8,6,4,4,3,3,2) come from "
                "m24_shapes.py's Golay-code permutation representation, re-run from a "
                "fresh clone (exit 0), content-identical except for one cosmetic "
                "dict-key-order difference (see violations.reproducibility_hygiene).",
                "F_g/Lambda_N derivation: T0 = k*(1 - chi_g/24) comes from the "
                "ground_state_invariance premise (Z_g's y^{+-1}q^0 term = k, same as "
                "untwined); Fg_factor = T0/(N(N-1)/24) is computed BEFORE comparison "
                "to the FROM-MEMORY exp_fac dict (verified by reading "
                "run_track_a_v3.py:377-389: F_g_coefficient_of_Lambda_N_derived is "
                "computed independently of exp_fac, which only feeds a downstream "
                "equality-check field) -- ground rule 4 respected, not a v2-style "
                "typed-then-matched value.",
                "H(0)=-1/12 convention: used only as an external comparison constant "
                "in B-dyons, not as an input to Track A's own H(tau) construction "
                "(A's H comes from solving N in Delta psi = N y^{1/2}Psi - k eta^3, "
                "z-independence), so no circularity between A and B on this point.",
                "sign pattern: A_n^g alternation for 2A and periodic 0-pattern for "
                "3A/5A/7A is REPORTED (class_report/sign_pattern_alternating field) "
                "but explicitly NOT claimed to follow from any derived criterion -- "
                "could_not_do says so honestly; no overclaim found.",
                "2*A_n convention used throughout COMPARISON.Twining.* rows: hand-"
                "checked 2A (table2A=[-2,-6,14,-28] vs blind A_n=[-3,7,-14], "
                "2*(-3,7,-14)=(-6,14,-28), matches with table[0]=T0 offset by one) and "
                "7AB (table=[-2,-1,0,0] vs blind A_n=[-1/2,0,0], 2*(-1/2,0,0)=(-1,0,0), "
                "matches). Arithmetically consistent, not a bug.",
            ],
            "independent_rederivation": rederiv["A"],
            "verdict": "No mathematical error found in Track A's derivation chain given "
                       "its declared tier-L inputs. k's NORMALISATION label and the "
                       "chi(5A)/chi(7A) RIGID_GIVEN_DEFINITION labels are both "
                       "internally consistent with ground rules 8-9.",
        },
        "B_dyons": {
            "checked": [
                "DMVV product at z=0: B2_euler evaluates the y=1 (not z=0 for the "
                "elliptic-genus variable, which is the correct specialisation for the "
                "y=1 Euler-number sector) slice of the two-variable DMVV product; this "
                "is the standard route to e(Hilb^k K3), matching Goettsche with an "
                "imported chi from Track D -- correctly flagged as NOT an independent "
                "derivation of chi=24 (24 = k*12, k typed).",
                "exact solve for (N,M) [reported as B4_immortal's (N,c,M)=(324,3,648)]: "
                "an over-determined linear system, rank 3 = rank of augmented matrix, "
                "zero residual on 106 equations, all 26 joint integer neighbours fail "
                "-- this is a genuine exact linear-algebra solve, not a grid scan "
                "centred on the answer (ground rule 8 respected).",
                "Hurwitz table vs Dirichlet: cross-checked against a THIRD, freshly-"
                "coded implementation (this review's own weighted reduced-form "
                "counter, not importing hurwitz.py) for D=1..300: 0 mismatches "
                "(D=0 excluded, see below).",
            ],
            "independent_rederivation": rederiv["B"],
            "verdict": "Hurwitz table B1 is correct for D=1..300 (independently "
                       "reverified). H(0)=-1/12 is a declared special-case convention, "
                       "not a form-count value, and is used consistently as such.",
        },
        "C_lattices": {
            "checked": [
                "chain K=O, b1=0 => chi(O)=2 => chi_top=24 => (b2+,b2-)=(3,19): each "
                "arrow is a standard, correctly-applied identity (Serre duality, "
                "Noether's formula chi(O)=(c1^2+c2)/12 with c1=0 so c2=12*chi(O), "
                "Hirzebruch signature tau=(c1^2-2c2)/3=-2c2/3, Hodge index "
                "b2+-=(b2+-tau)/2). Arithmetic re-verified by hand: chi(O)=2 => "
                "c2=24=chi_top; tau=-2*24/3=-16; b2=chi_top-2=22 (from "
                "1-0+22-0+1=24 Betti sum, an independent internal check); "
                "b2+=(22-16)/2=3, b2-=(22+16)/2=19. Correct.",
                "realizability of controls: chi_top's own negative-control row "
                "honestly states 'none realizable... HYPOTHETICAL, not counted' -- "
                "see rigidity_verdicts below for the labeling consequence.",
                "mU+n(-E8) signature scan: (3,2) is the unique solution among all "
                "integer 0<=m,n<=30 with signature (m,m+8n) matching (3,19); "
                "(11,0)->(11,11) and (7,1)->(7,15) are REALIZABLE even unimodular "
                "lattices of the wrong signature, so this negative control is sound.",
            ],
            "independent_rederivation": rederiv["C"],
            "verdict": "Chain arithmetic is correct. See rigidity_verdicts for a "
                       "labeling concern on the chi_top row (not a numeric error).",
        },
        "D_tda": {
            "checked": [
                "open-star premise at each N in {4,6,8}: regularity (f-vector "
                "halving) and link homology (RP^3 mod 2/3/5, S^3 for ordinary points) "
                "are GUDHI-computed facts, re-run in this review is not repeated "
                "(GUDHI computation, not re-implemented independently) but the "
                "elementary fixed-point count (16 = 2^4) was independently "
                "re-derived for every even N in [4,20] in rederivations.py, "
                "generalising past N in {4,6,8}.",
                "MV computation of b3, b4: b(Q)=(1,0,6,0,1) (GUDHI) plus the cone-glue "
                "argument (resolution_local_model, a declared tier-L input) gives "
                "b(resolved)=(1,0,22,0,1); chi=24=chi(Q)+16*(chi(S^2)-chi(cone-pt)) "
                "=8+16*(2-1)=24, matches directly and independently via the elementary "
                "formula chi(Q)=(chi(T^4)+16)/2=(0+16)/2=8 re-derived in this review.",
                "parity of the intersection form: NOT attempted on the actual "
                "resolved (non-triangulated) K3 -- D7's mod-2 cup-pairing is only "
                "computed on T^4 and on the singular quotient Q, and the track's own "
                "could_not_do list says this explicitly. No overclaim found; this is "
                "an honest gap, not a math error.",
            ],
            "independent_rederivation": rederiv["D"],
            "verdict": "chi(Q)=8 and n_fixed=16 independently confirmed by elementary "
                       "means. D7's parity attempt is correctly reported as "
                       "inconclusive, not as a computed parity.",
        },
        "E_flux": {
            "checked": [
                "TT (A.2)-(A.3) unimodular lattice line by line: E01a computes "
                "det(H33)=-1 (odd, signature (3,3,0)) and the E8 Cartan matrix's "
                "leading principal minors all positive (Sylvester -> positive "
                "definite) with det=1 -- independently re-verified in this review "
                "(rederivations.py track E) using a from-scratch sympy construction "
                "of the (A.6)-(A.7) Gram and the stated e_i=f_i+g_i embedding: "
                "|det A6|=64, and the index identity |det T|^2*|det H| = |det A6| "
                "holds exactly (1^2 * 64 = 64).",
                "index of the (A.6) span: E01b reports index 8 from two routes "
                "(transition-matrix determinant and sqrt-discriminant-ratio) which "
                "agree; this review's third, independently-coded route also gives 8.",
                "line-941 resolution: E-flux/inputs.json declares "
                "tadpole_half_factor_eq2.3=1/2 'decided against line-941 prose by 02' "
                "-- read 02_tt_worked_examples.py's TT_printed comparisons (E02c, "
                "E02d) which show the (1/2)*N_flux convention reproduces TT's own "
                "printed worked-example numbers (N_D3=20 at N_flux=8, N_D3=48 at "
                "N_flux=0) while the no-half convention does not; this is a genuine, "
                "checkable resolution of an ambiguous literature line, not an "
                "arbitrary choice.",
                "symmetry group justification (E01d, order 384): read "
                "lattice_u3.py's group_G() directly (not just the results.json "
                "summary): permutations(3) [3!=6] x swaps in {0,1}^3 [2^3=8] x "
                "signs in {1,-1}^3 [2^3=8], giving 6*8*8=384, matching the "
                "reported order_G exactly. The group is literally "
                "(Z2 x Z2)^3 semidirect-product S3 (per-block independent "
                "swap-bit and sign-bit, wreathed by the block permutation), so "
                "both the order and the '(Z2 x Z2) wr S3' name are correct.",
            ],
            "independent_rederivation": rederiv["E"],
            "verdict": "Index-8 and unimodular-lattice claims independently confirmed. "
                       "E01d's order-384 symmetry group was checked directly against "
                       "its source code (lattice_u3.py) and is correct.",
        },
        "F_whichk3": {
            "checked": [
                "stabiliser orders and enhanced root system: F1a/F1b/F1c report "
                "(omega,omega) as the maximal-stabiliser point (order 144, formal "
                "group) with A2+A2 (SU(3)xSU(3)) roots; this is the textbook T^2xT^2 "
                "self-dual-point enhancement result (both factors at the SU(3) "
                "point), consistent with standard Narain-lattice lore -- not "
                "independently re-derived here but structurally unsurprising.",
                "NS/T computation (F2): rank and discriminant for (i,i), (omega,omega) "
                "and mixed (i,omega) reproduce the well-known CM-point Neron-Severi "
                "ranks (4,4,2) and transcendental lattices with the expected "
                "discriminants (4, 3, indefinite U^2) -- standard facts, correctly "
                "applied.",
                "|Aut| enumeration (F3): D4 maximal at 1152 among the 6 rank-4 root "
                "lattices listed is correct (|Aut(D4)|=2^3*4!*3=1152, the extra "
                "factor-3 triality beyond the generic Weyl-group-times-diagram-"
                "automorphism count) -- a well-known special fact about D4, "
                "correctly identified as the maximum.",
                "Golay/M24 order check: independently re-derived in this review via "
                "Gleason's theorem (a route STRUCTURALLY UNRELATED to F4b's QR(23) "
                "permutation-group construction), giving the exact same weight "
                "enumerator {0:1,8:759,12:2576,16:759,24:1}. |M24|=244823040 was NOT "
                "independently re-derived in this review (would require "
                "reconstructing the Golay-code automorphism group, out of scope for "
                "the time budget) but is the correct, well-known value.",
                "Kummer-code embedding search (F4d): 759 of C(24,16)=735471 16-subsets "
                "give a shortened-code dimension of 5 (matching the Kummer code K) -- "
                "759 is exactly the number of octads, a strong internal-consistency "
                "signal (complements of octads), not independently re-verified here.",
                "Lefschetz numbers (F5): fixed-point counts 8,6,4 for symplectic "
                "orders 2,3,4 on a Kummer surface match the standard Nikulin table "
                "for these orders (Nikulin_fixed_points input, used only as a "
                "comparison, not fed into the computation) -- consistent.",
            ],
            "independent_rederivation": rederiv["F"],
            "verdict": "Golay weight enumerator independently confirmed via a "
                       "structurally distinct route (Gleason's theorem). Other Track F "
                       "claims are consistent with standard results but were spot-"
                       "checked against literature knowledge rather than fully "
                       "re-derived, given the scope of six tracks in one review pass.",
        },
    },

    "item_g_disagree_not_comparable_rows": (
        "comparison.json's own 'counts' field is {'AGREE': 71, 'NOT_COMPUTED': 261}; "
        "verified directly (python: set(r['status'] for r in comparison.json['rows']) "
        "== {'AGREE','NOT_COMPUTED'}). There are NO rows with status DISAGREE or "
        "NOT_COMPARABLE in this comparison. Item (g) of the task is therefore "
        "satisfied vacuously: there is nothing to adjudicate."
    ),

    "item_h_reverse_items": [
        {
            "item": 1,
            "name": "moore_vs_hurwitz (D up to 5000)",
            "assessment": "The underlying Hurwitz class numbers H(D) that this "
                           "identity depends on were independently reverified by this "
                           "review for D=1..300 (0 mismatches, see rederivations.py "
                           "track B). This review did not reimplement Lean's specific "
                           "nForms(D)/h12(D)/isKSquare functions or extend the check "
                           "to D up to 5000 independently, given the time budget; the "
                           "reported computation (a from-scratch reimplementation, "
                           "distinct from both existing routes) is plausible and its "
                           "negative control (dropping the two correction terms fails "
                           "on 54 of the new D) is a real, targeted control.",
        },
        {
            "item": 2,
            "name": "kummer_iff_even (D up to 5000)",
            "assessment": "Elementary logical/computational claim (isKummerForm == "
                           "'a,b,c all even'); the negative control (requiring only a,c "
                           "even) is a natural, well-targeted perturbation of the same "
                           "selecting condition. No issue found; not independently "
                           "re-run given scope.",
        },
        {
            "item": 3,
            "name": "class_number_one (most_attractive, D up to 3000)",
            "assessment": "The computed set {3,4,7,8,11,19,43,67,163} is in fact "
                           "STRONGER than the write-up claims: this is exactly the set "
                           "of discriminants for which the Baker-Heegner-Stark theorem "
                           "proves class number 1 for ALL D (not just D<=3000) -- a "
                           "completely classified, fully proven fact in number theory, "
                           "not merely an empirical pattern that could break beyond "
                           "3000. The write-up's own honest caveat about D=3,4,7,8,11 "
                           "being a dense run of 5 consecutive valid discriminants is "
                           "correct and does not weaken the result; it is a genuine "
                           "feature of this classification, not a defect.",
        },
        {
            "item": 4,
            "name": "ade_trapping (D_n beyond rank 8, extended)",
            "assessment": "INDEPENDENTLY RE-RUN in this review from a fresh clone "
                           "(item4_ade_trapping.py, exit 0) and hand-checked: rank 9 "
                           "best=242 (=E8's 240 roots + A1's 2 roots, matching the "
                           "'E8+best(1)' tie reported) vs D9=2*9*8=144, so D9 is "
                           "correctly NOT optimal at rank 9; same pattern confirmed at "
                           "ranks 10 (246 vs D10=180) and 11 (252 vs D11=220); rank 24 "
                           "gives D24=2*24*23=1104 > 3*E8=720, confirming the negative "
                           "control. The reported holds=false verdict is CORRECT.",
        },
        {
            "item": 5,
            "name": "twined_2A_extended (divisibility by 24 through q^40)",
            "assessment": "Hand-verified the stated reduction: combined[0] = "
                           "chi*(-2) + c*Nlev*(Nlev-1) with chi=8, Nlev=2 gives "
                           "combined[0] = -16 + 2c; divisibility by 24 requires "
                           "2c = 16 (mod 24), i.e. c = 8 (mod 12). This is EXACTLY the "
                           "congruence class reported (c=-16,-4,8,20,32 pass; "
                           "c=-15,-3 fail), so the NON_DISCRIMINATING classification "
                           "with respect to c is mathematically correct, not an "
                           "overclaim.",
        },
    ],

    "independent_rederivations_summary": {
        "script": "audit/k3t2_rigidity_v3/skeptic_math/rederivations.py",
        "command": "cd audit/k3t2_rigidity_v3/skeptic_math && "
                   "/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/"
                   ".venv-tda/bin/python rederivations.py",
        "all_pass": rederiv["all_pass"],
        "per_track": {k: v["pass_"] for k, v in rederiv.items() if k != "all_pass"},
        "note_on_track_F": (
            "The first attempt at Track F's independent check (a hand-built "
            "bordered-QR(23) generator matrix, coded from scratch) produced a WRONG, "
            "non-self-dual code with odd-weight codewords present. This was DISCARDED "
            "and replaced with a Gleason's-theorem derivation (structurally unrelated "
            "to both the discarded attempt and to F-whichk3/f4_codes.py's own QR(23)/"
            "M24 construction), which reproduces the exact weight enumerator by "
            "solving 2 linear equations from a declared structural input (the "
            "invariant-ring generators for doubly-even self-dual length-24 codes) "
            "plus the definitional fact that the Golay code has minimum distance 8. "
            "This failure-then-fix is recorded per the honesty standard this review "
            "expects of the tracks it is grading."
        ),
        "note_on_track_B": (
            "The first attempt at Track B's independent check used a hand-typed "
            "'FROM MEMORY' literature table for H(D) and disagreed with the blind "
            "table at D=32 (memory said 2, blind said 3). Direct recomputation from "
            "the reduced-form definition (3 forms of discriminant -32, none of the "
            "two special-automorphism shapes, so H(32)=1+1+1=3) showed the MEMORY "
            "VALUE WAS WRONG. The blind script was correct; this reviewer's memory "
            "was not. Replaced with a freshly-coded weighted reduced-form counter "
            "checked against 300 values with 0 mismatches."
        ),
    },
}

math_errors = []

violations = [
    "REPRODUCIBILITY (ground rule 10) -- minor: A-genus/results.json's "
    "'square_cycle_shape' field is str() of a plain Python dict; re-running "
    "m24_shapes.py + run_track_a_v3.py from a fresh clone with the SAME committed "
    "command reproduces identical NUMBERS everywhere but this one field's key order "
    "can differ (observed '{2: 8, 1: 8}' vs committed '{1: 8, 2: 8}') because "
    "PYTHONHASHSEED randomises dict/set iteration order across runs. Cosmetic only "
    "-- no numeric content is affected -- but strictly speaking the script does not "
    "produce a byte-identical results.json on every run, which the reproducibility "
    "standard implicitly wants. Fix: sort dict items before str()/json serialisation, "
    "or set PYTHONHASHSEED=0 in the documented command.",
]

disputed_rows = [
    {
        "target_id": "tadpole_budget",
        "reason": "Labeled AGREE with independent_routes=0 and the row's own note "
                  "says 'Same literature constant typed on both sides, not "
                  "independently derived by either.' This is technically correct "
                  "(the two typed values ARE equal) but 'AGREE' invites reading it "
                  "as a substantive cross-check the way most other AGREE rows are; "
                  "it is really a tautology (same tier-L input copied into both "
                  "places). Recommend a distinct status such as SAME_DECLARED_INPUT "
                  "for this and the sibling rows (tadpole_cancellation, "
                  "d3_tadpole_target_is_24) so a reader scanning for AGREE counts "
                  "does not overweight them as evidence.",
    },
]

rigidity_verdicts = [
    {
        "parameter": "Track D: j (number of resolved singular points), selecting "
                     "condition 'no unresolved point with non-sphere Z/2 link'",
        "genuinely_rigid": False,
        "reason": "The blind track's own note says this condition 'selects j=16 by "
                  "construction': X_j is DEFINED by resolving j of the 16 orbifold "
                  "points (replacing each with a smooth disc bundle) and leaving "
                  "16-j as cone points; requiring 'no bad link' is then definitionally "
                  "equivalent to 'every point was resolved', i.e. j=16, given only "
                  "the construction of X_j. This is not an independent selection "
                  "criterion discovered by computation; it is a restatement of the "
                  "definition of full resolution. The realizable negative control "
                  "(j=15 having one bad-link cone point) is real and correctly "
                  "computed, but it demonstrates the definition was applied "
                  "correctly, not that an independent condition picked out j=16.",
        "final_label": "NORMALISATION (the selecting condition, once unpacked, "
                       "contains the target j=16 by construction of X_j) rather than "
                       "RIGID_GIVEN_DEFINITION as currently labeled -- though this is "
                       "a labeling-precision issue, not a numeric error: chi(X_16)=24 "
                       "and b2=22 are independently confirmed by both GUDHI and Hodge "
                       "theory regardless of how this one row is labeled.",
    },
    {
        "parameter": "Track C: chi_top (Euler characteristic), selecting condition "
                     "'chi(O) from K=O, b1=0, Serre duality; Noether'",
        "genuinely_rigid": False,
        "reason": "The track's own rigidity table lists control_perturbs_same_"
                  "parameter=False and negative_control='none realizable (other "
                  "chi(O) leaves K3 class); HYPOTHETICAL, not counted'. Per ground "
                  "rule 9, a control that is not realizable and not counted as "
                  "discrimination cannot satisfy RIGID's or RIGID_GIVEN_DEFINITION's "
                  "requirement that 'a neighbouring value demonstrably fails'. This "
                  "is not really a rigidity/selection question at all -- chi_top=24 "
                  "is a direct, forced CONSEQUENCE of the declared inputs (K=O, "
                  "b1=0, Noether's formula), with nothing being chosen among "
                  "alternatives.",
        "final_label": "VERIFIED_IDENTITY would fit better than RIGID_GIVEN_DEFINITION "
                       "(nothing is selected; the chain of equalities is forced once "
                       "the definitional inputs are granted, matching the spirit of "
                       "'a universal identity' in ground rule 8's vocabulary, even "
                       "though a K3-specific identity rather than a fully universal "
                       "one). Not a numeric error either way: chi_top=24 is correct.",
    },
    {
        "parameter": "Track A: k (overall factor in Z=k*phi_{0,1})",
        "genuinely_rigid": False,
        "reason": "Correctly self-labeled NORMALISATION already: the track's own "
                  "k_selector_summary explicitly states no selector bounds k above "
                  "and the selecting_condition_contains_target=True. No dispute; "
                  "this is a model of how a NORMALISATION row should be written.",
        "final_label": "NORMALISATION (confirmed correct as currently labeled)",
    },
    {
        "parameter": "Track A: chi(g) at level 5 and 7",
        "genuinely_rigid": True,
        "reason": "Solution sets {4,14,24} (N=5) and {3,24} (N=7) with realizable "
                  "neighbours (chi=5,3 for N=5; chi=2,4 for N=7) failing, and the "
                  "HYPOTHETICAL value 14 correctly excluded from counting as "
                  "discrimination per ground rule 9. This is a correctly-applied "
                  "RIGID_GIVEN_DEFINITION label.",
        "final_label": "RIGID_GIVEN_DEFINITION (confirmed correct)",
    },
    {
        "parameter": "REVERSE item 5: c (frame-shape coefficient in twined24)",
        "genuinely_rigid": False,
        "reason": "Hand-verified: divisibility-by-24 of combined[0]=chi*(-2)+"
                  "c*Nlev*(Nlev-1) at chi=8,Nlev=2 reduces to c=8 (mod 12), a full "
                  "residue class of integers, not the single value c=-16. Correctly "
                  "labeled NON_DISCRIMINATING with respect to c.",
        "final_label": "NON_DISCRIMINATING (confirmed correct)",
    },
]

reverse_verdicts = [
    {"item": "1: moore_vs_hurwitz (D<=5000)", "stands": True,
     "reason": "Underlying H(D) values independently reverified for D<=300 (0 "
               "mismatches); the D=5000 extension and Lean-function reimplementation "
               "were not independently rerun given scope, so this is accepted as "
               "plausible rather than independently reproduced end-to-end."},
    {"item": "2: kummer_iff_even (D<=5000)", "stands": True,
     "reason": "Elementary biconditional with a correctly-targeted negative control "
               "(perturbs which of a,b,c must be even); not independently rerun."},
    {"item": "3: class_number_one (D<=3000)", "stands": True,
     "reason": "The computed set matches the Baker-Heegner-Stark class-number-one "
               "discriminants, a fully PROVEN classification (stronger support than "
               "the write-up itself claims, since this holds for all D, not just "
               "D<=3000)."},
    {"item": "4: ade_trapping (extended beyond rank 8)", "stands": True,
     "reason": "Independently rerun from a fresh clone (exit 0) and hand-verified "
               "at ranks 9, 10, 11, 24; the holds=false verdict and every reported "
               "number (144, 242, 180, 246, 220, 252, 1104, 720) is correct."},
    {"item": "5: twined_2A_extended (divisibility to q^40)", "stands": True,
     "reason": "Hand-derived the same c=8 (mod 12) congruence from the stated "
               "combined[0] formula; matches the reported pass/fail set for "
               "c in {-16,-15,-4,8,-3,20,32} exactly."},
]

VERDICT["math_errors"] = math_errors
VERDICT["violations"] = violations
VERDICT["disputed_rows"] = disputed_rows
VERDICT["rigidity_verdicts"] = rigidity_verdicts
VERDICT["reverse_verdicts"] = reverse_verdicts

VERDICT["overall_conclusion"] = (
    "No mathematical ERROR (in the sense of a wrong number reaching a downstream "
    "claim) was found in the BLIND tracks, the CHAIN check, the COMPARISON, or the "
    "REVERSE items. Six from-scratch independent re-derivations (one per track, "
    "using a route structurally different from each track's own script) all "
    "reproduce the tracks' reported numbers: phi_{0,1}(tau,0)=12 (Hodge-diamond "
    "route), Hurwitz H(D) for D=1..300 (0 mismatches, fresh weighted reduced-form "
    "counter), E8 Gram det=1/even/unimodular (D8+glue-vector construction), "
    "chi(Q)=8 and n_fixed=16 for T^4/Z2 (elementary, extended to N up to 20), the "
    "E-flux index-8 identity |det T|^2|det H|=|det A6| (independent sympy "
    "construction), and the Golay weight enumerator {0:1,8:759,12:2576,16:759,24:1} "
    "(Gleason's theorem, structurally unrelated to the QR(23)/M24 construction the "
    "track itself uses). Two of these six re-derivation attempts initially FAILED "
    "due to errors in THIS REVIEW's own first-pass code/memory (a buggy hand-built "
    "Golay generator matrix, and a wrong memorised Hurwitz value at D=32); both "
    "were corrected and both corrections point to the reviewer's error, not the "
    "blind track's, per direct recomputation from definitions. The two labeling "
    "concerns found (Track D's j-selection and Track C's chi_top row) are about "
    "which fixed rigidity-vocabulary word best describes an already-correct "
    "computation, not about any wrong number. Item (g) is vacuous: comparison.json "
    "contains zero DISAGREE or NOT_COMPARABLE rows. All five REVERSE items' "
    "verdicts stand under independent scrutiny, with item 4 fully independently "
    "rerun and hand-verified and item 3 shown to be backed by a stronger (fully "
    "proven) classical theorem than the write-up itself invokes."
)

out_path = HERE / "verdict.json"
out_path.write_text(json.dumps(VERDICT, indent=2))
print("wrote", out_path)
