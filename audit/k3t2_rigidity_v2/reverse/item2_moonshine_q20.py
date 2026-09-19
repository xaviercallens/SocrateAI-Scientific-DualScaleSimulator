"""
Reverse-pass item 2 (and item 4, `twinedAll_div`/`twined_div24`, sharing this engine):
extend LeanMaster's `Decompositions.moonshine_modules_decompose` and
`TwiningAll.twinedAll_div`/`Twining.twined_div24` (DualScaleMoonshine) from Lean's checked
range n = 0..9 to n = 10..20, using the independent engine in moonshine_engine.py.

PREDICTION (fixed before `computed` below is touched):
  * `twined_div24_q20`: every coefficient of the 26 computed twined series (`24*D*H_g`
    normalised) stays exactly divisible by `24*D` at n = 10..20 (D=1 for the classes using
    `twined24` directly). This is expected because divisibility by 24*D is forced by the
    modular-weight-2 structure of `24*D*H_g = (D*chi_g*numer + 24*D*F_g)/eta^3`; Lean's own
    proof by `decide` is order-by-order, so nothing in it special-cases n<=9.
  * `moonshine_modules_decompose_q10_20`: <T_n, chi_i> (n=10..20, i=0..25) has vanishing
    omega_7/omega_15/omega_23 parts, is divisible by |M24|=244823040, and the resulting
    integer multiplicity is >= 0 for every class i. This is the SAME structural condition
    Lean's `decide +kernel` checks for n=0..9 (26 conditions per n, none built from a known
    "correct" multiplicity -- the quotient by 244823040 is arithmetic, not a target). Since
    CDH's own printed Table 48 stops at n=9, there is no literature value to compare
    against for n=10..20: a genuine new fact, not a check against a known answer.

VALIDATION GATE (run first): reproduce bit-for-bit, with THIS independent engine, Lean's
own theorems at N=9 -- twined_2A/3A/5A/7AB (Twining.lean), twined_2B..23AB and
twinedAll_div (TwiningAll.lean), and moonshine_modules_decompose (n=0..9 against CDH Table
48, Decompositions.lean) -- before trusting the engine to extrapolate.

Negative control: perturb the coefficient "-16" in `twined24(N, 8, 2, -16)` (the F_2A
combination) by +1 and show `twined_div24`-style divisibility breaks at n=10..20 (not just
at n<=9); also perturb one CHAR_TAB entry by +1 and show a level in n=10..20 loses
integrality/non-negativity, i.e. the checks are not vacuous at the new orders either.

Run: /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python item2_moonshine_q20.py
Writes item2_moonshine_q20_results.json.
"""
import json
import moonshine_engine as E

N_LEAN = 9
N_TARGET = 20

# ---------------------------------------------------------------------------
# VALIDATION GATE: reproduce Lean's own N=9 theorems with this fresh engine.
# ---------------------------------------------------------------------------
table2A = [-2, -6, 14, -28, 42, -56, 86, -138, 188, -238]
table3A = [-2, 0, -6, 10, 0, -18, 20, 0, -30, 42]
table5A = [-2, 0, 2, 0, -6, 2, 0, 6, 0, -10]
table7AB = [-2, -1, 0, 0, 4, 0, -2, 2, -3, 0]

val_2A = E.twined24(9, 8, 2, -16) == [24 * x for x in table2A]
val_3A = E.twined24(9, 6, 3, -6) == [24 * x for x in table3A]
val_5A = E.twined24(9, 4, 5, -2) == [24 * x for x in table5A]
val_7AB = E.twined24(9, 3, 7, -1) == [24 * x for x in table7AB]

tables16 = {
    "2B": [-2, 10, -18, 20, -38, 72, -90, 118, -180, 258],
    "3B": [-2, 6, 0, -14, 12, 0, -16, 30, 0, -42],
    "4A": [-2, -6, -2, 4, -6, -8, 6, 6, -4, -14],
    "4B": [-2, 2, -2, -4, 2, 8, -2, -10, 4, 10],
    "4C": [-2, 2, 6, -4, -6, 0, 6, -2, -12, 10],
    "6A": [-2, 0, 2, 2, 0, -2, -4, 0, 2, 2],
    "6B": [-2, -2, 0, 2, 4, 0, 0, -2, 0, 6],
    "8A": [-2, -2, -2, 0, -2, 0, 2, -2, 0, -2],
    "10A": [-2, 0, 2, 0, 2, 2, 0, -2, 0, -2],
    "11A": [-2, 2, 0, 0, 0, -2, 0, -2, 2, 0],
    "12A": [-2, 0, -2, -2, 0, -2, 0, 0, 2, -2],
    "12B": [-2, 2, 0, 2, 0, 0, 0, -2, 0, -2],
    "14AB": [-2, 1, 0, 0, 0, 0, 2, 2, -1, 0],
    "15AB": [-2, 0, -1, 0, 0, 2, 0, 0, 0, 2],
    "21AB": [-2, -1, 0, 0, -2, 0, -2, 2, 0, 0],
    "23AB": [-2, -2, 2, -1, 0, 0, 0, 0, 0, 0],
}
d9 = E.data_all(9)
val16 = {}
for name, tbl in tables16.items():
    computed = E.twined_d(9, d9[name])
    expected = [24 * d9[name][0] * x for x in tbl]
    val16[name] = (computed == expected)

val_div24_at9 = all(all(x % (24 * d9[name][0]) == 0 for x in E.twined_d(9, d9[name])) for name in tables16)

table48 = [
 [-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
 [0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
 [0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
 [0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0],
 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0],
 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2],
 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,2,0,0,0,2,2,2,2],
 [0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,2,0,0,2,2,2,4,2,2,6],
 [0,0,0,0,0,0,0,0,2,2,2,0,0,2,2,2,0,2,2,2,4,4,4,8,8,10],
]

series9 = E.computed_series(9)
decompose9_ok = True
mismatches9 = []
for n in range(10):
    for i in range(26):
        v = E.moonshine_inner(9, n, i, series9)
        want = (E.ORDER_M24 * table48[n][i], 0, 0, 0)
        if v != want:
            decompose9_ok = False
            mismatches9.append({"n": n, "i": i, "got": v, "want": want})

validation = {
    "twined_2A": val_2A, "twined_3A": val_3A, "twined_5A": val_5A, "twined_7AB": val_7AB,
    "twined_all16_at_q9": val16,
    "twinedAll_div_at_q9": val_div24_at9,
    "moonshine_modules_decompose_at_n0_9": decompose9_ok,
    "mismatches_n0_9": mismatches9[:5],
    "all_pass": all([val_2A, val_3A, val_5A, val_7AB, val_div24_at9, decompose9_ok] + list(val16.values())),
}

# ---------------------------------------------------------------------------
# EXTENSION: N = 20.
# ---------------------------------------------------------------------------
series20 = E.computed_series(N_TARGET)
d20 = E.data_all(N_TARGET)

div24_new = {}
for name in tables16:
    coeffs = E.twined_d(N_TARGET, d20[name])
    D = d20[name][0]
    bad = [(n, coeffs[n]) for n in range(N_LEAN + 1, N_TARGET + 1) if coeffs[n] % (24 * D) != 0]
    div24_new[name] = {"holds_n10_20": len(bad) == 0, "violations": bad}
# also the four twined24-direct classes (1A,2A,3A,5A,7AB): divisibility by 24
for name, args in [("1A", (24, 1, 0)), ("2A", (8, 2, -16)), ("3A", (6, 3, -6)),
                    ("5A", (4, 5, -2)), ("7AB", (3, 7, -1))]:
    coeffs = E.twined24(N_TARGET, *args)
    bad = [(n, coeffs[n]) for n in range(N_LEAN + 1, N_TARGET + 1) if coeffs[n] % 24 != 0]
    div24_new[name] = {"holds_n10_20": len(bad) == 0, "violations": bad}

div24_all_hold = all(v["holds_n10_20"] for v in div24_new.values())

decompose_new = {}
all_ok = True
new_multiplicities = {}
for n in range(N_LEAN + 1, N_TARGET + 1):
    row_ok = True
    row_mult = []
    for i in range(26):
        v = E.moonshine_inner(N_TARGET, n, i, series20)
        rational, w7, w15, w23 = v
        omega_zero = (w7 == 0 and w15 == 0 and w23 == 0)
        divides = (rational % E.ORDER_M24 == 0)
        mult = rational // E.ORDER_M24 if divides else None
        nonneg = (mult is not None and mult >= 0)
        ok = omega_zero and divides and nonneg
        row_ok = row_ok and ok
        row_mult.append(mult)
        if not ok:
            all_ok = False
    decompose_new[n] = row_ok
    new_multiplicities[n] = row_mult

# ---------------------------------------------------------------------------
# NEGATIVE CONTROLS at the new orders.
# ---------------------------------------------------------------------------
# (a) perturb F_2A's coefficient -16 -> -15; divisibility by 24 should break somewhere in 10..20
coeffs_pert = E.twined24(N_TARGET, 8, 2, -15)
bad_pert = [(n, coeffs_pert[n]) for n in range(N_LEAN + 1, N_TARGET + 1) if coeffs_pert[n] % 24 != 0]
neg_control_a = {"perturbed": "twined24(N,8,2,-15) instead of -16",
                  "violations_n10_20": bad_pert, "control_has_teeth": len(bad_pert) > 0}

# (b) perturb CHAR_TAB[1] (the 2A column's row, chi_2 dimension 23 entry) by +1 at class j=0
CHAR_TAB_PERT = [list(row) for row in E.CHAR_TAB]
CHAR_TAB_PERT[1] = list(CHAR_TAB_PERT[1])
CHAR_TAB_PERT[1][0] = (CHAR_TAB_PERT[1][0][0] + 1, CHAR_TAB_PERT[1][0][1])
broke_at = []
for n in range(N_LEAN + 1, N_TARGET + 1):
    row = [(series20[j][n] if n < len(series20[j]) else 0, 0) for j in range(26)]
    v = E.inner4(row, CHAR_TAB_PERT[1])
    rational, w7, w15, w23 = v
    ok = (w7 == 0 and w15 == 0 and w23 == 0 and rational % E.ORDER_M24 == 0 and rational // E.ORDER_M24 >= 0)
    if not ok:
        broke_at.append(n)
neg_control_b = {"perturbed": "CHAR_TAB[1][0] (chi_2 at 1A) +1",
                  "levels_n10_20_where_integrality_or_nonneg_breaks": broke_at,
                  "control_has_teeth": len(broke_at) > 0}

result = {
    "prediction": {
        "twined_div24_q20": "every coefficient of the 26 computed twined series stays "
            "divisible by 24*D at n=10..20 (structural: forced by the weight-2 modular "
            "combination, not by any target value)",
        "moonshine_modules_decompose_q10_20": "<T_n,chi_i> has vanishing omega-parts, is "
            "divisible by |M24|=244823040, and the quotient is >=0, for n=10..20, i=0..25 "
            "-- a genuine new fact since CDH Table 48 stops at n=9",
    },
    "validation_gate_lean_n9": validation,
    "lean_theorems_extended": [
        "Twining.twined_div24 (DualScaleMoonshine)",
        "TwiningAll.twinedAll_div (DualScaleMoonshine)",
        "Decompositions.moonshine_modules_decompose (DualScaleMoonshine)",
    ],
    "computed": {
        "twined_divisibility_q10_20_per_class": div24_new,
        "twined_divisibility_all_hold_q10_20": div24_all_hold,
        "moonshine_modules_decompose_holds_per_n": decompose_new,
        "moonshine_modules_decompose_all_hold_q10_20": all_ok,
        "multiplicities_n10_20": new_multiplicities,
    },
    "negative_controls": {"perturb_F2A_coefficient": neg_control_a,
                           "perturb_character_table_entry": neg_control_b},
    "tier": "B",
    "command": "cd audit/k3t2_rigidity_v2/reverse && "
               "/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python item2_moonshine_q20.py",
}

with open("item2_moonshine_q20_results.json", "w") as f:
    json.dump(result, f, indent=1)

print(json.dumps({
    "validation_all_pass": validation["all_pass"],
    "div24_all_hold_q10_20": div24_all_hold,
    "decompose_all_hold_q10_20": all_ok,
    "neg_control_a_has_teeth": neg_control_a["control_has_teeth"],
    "neg_control_b_has_teeth": neg_control_b["control_has_teeth"],
}, indent=1))
