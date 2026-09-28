#!/usr/bin/env python3
"""E-flux v3, script 06: tier-L quotation arithmetic (kept OUT of results), unfixed-moduli quotations, TT section 6
(null-vector shifts) check, Euler-characteristic bookkeeping using Track D's exports.

Run: cd audit/k3t2_rigidity_v3/E-flux/scripts && <venv python> 06_literature_and_moduli.py
(no arguments; needs results/00_track_d_poll.json written by poll_track_d.py)
"""
import json, itertools
import sympy as sp
from sympy import Rational as R
from common import *
from lattice_u3 import H33
from susy import gram_of, sig_rank
from inputs_decl import TADPOLE_TOTAL

SCRIPT = "06_literature_and_moduli.py"
CMD = rel_command(SCRIPT)
H = sp.Matrix(H33)

def main():
    poll = json.load(open(RES / "00_track_d_poll.json"))
    n_sing = poll["T4Z2_fixed_points_from_track_D"]
    n_O7 = poll["T2Z2_fixed_points_computed_here_tierB"]
    dexp = json.load(open(ROOT / poll["track_d_export_used"]))
    lit = []
    # (L) TT eq (2.3) bookkeeping: each O7 induces 2 units of D3 charge, each of the 16 D7 one unit  [quoted numbers: 2, 16, 1]
    q_o7 = quote("plane is located at each of these fixed points")
    q_d7 = quote("16 D7-branes need to be added")
    q_charge = quote("2 units of D3 brane charge")
    q_tot = quote("a total of 24 units of three brane charge")
    D3_per_O7, N_D7, D3_per_D7 = 2, 16, 1   # QUOTED from TT lines above (tier L), not computed
    val = n_O7 * D3_per_O7 + N_D7 * D3_per_D7
    lit.append({"id": "L06a", "claim": "TT: total induced D3 charge = (#O7)*2 + (#D7)*1", "tier": "L",
                "arithmetic": f"{n_O7}*{D3_per_O7} + {N_D7}*{D3_per_D7} = {val}", "value": val,
                "inputs": {"#O7": f"{n_O7} (computed here, tier B: half-lattice fixed points of x->-x on T2; from 00_track_d_poll.json)", "D3 per O7 = 2": "quoted TT", "#D7 = 16": "quoted TT", "D3 per D7 = 1": "quoted TT"},
                "expected_from_source": TADPOLE_TOTAL, "source_lines": [q_o7, q_d7, q_charge, q_tot], "agrees": val == TADPOLE_TOTAL,
                "independence_note": "this is quotation arithmetic: 2, 16 and 1 are quoted; only the O7 count is computed. It is not an independent derivation of 24."})
    lit.append({"id": "L06b", "claim": "TT footnote 1: F-theory 24 (p,q) 7-branes with 1 unit of D3 charge each", "tier": "L", "arithmetic": "24*1 = 24", "value": 24 * 1,
                "expected_from_source": TADPOLE_TOTAL, "source_lines": [quote("there are 24 (p, q) 7-branes each of which acquires one unit of 3-brane charge")], "agrees": 24 * 1 == TADPOLE_TOTAL,
                "independence_note": "24 and 1 are both quoted; the 24 is the same input, not a second derivation"})
    bl = BL_TXT.read_text(encoding="utf-8").split("\n")
    qb1 = quote("n1 + n5 = 24", bl); qb2 = quote("Since 24 are expected", bl)
    lit.append({"id": "L06c", "claim": "BLPSSW hep-th/9605184: 24 expected, 8 five-branes in the GP model, 24-8 = 16 = number of T4/Z2 orbifold singularities (hidden instantons)", "tier": "L",
                "arithmetic": f"24 - 8 = {24 - 8}; #T4/Z2 singular points from Track D = {n_sing}", "value": 24 - 8, "expected_from_source": n_sing,
                "inputs": {"24": "quoted BLPSSW n1+n5=24", "8": "quoted BLPSSW eight five-branes", "#singular": f"{n_sing} read from {poll['track_d_export_used']} (also computed here tier B: {poll['T4Z2_fixed_points_computed_here_tierB']})"},
                "source_lines": [qb1, qb2], "agrees": (24 - 8) == n_sing,
                "note": "BLPSSW themselves say 'one might intuitively think that one instanton is hiding at each singularity': the equality 16(hidden)+8=24 is an inference from their subtraction and singularity count, not a printed sentence"})
    chi_T4 = 0
    chiK3 = R(chi_T4 - n_sing, 2) + n_sing * 2
    lit.append({"id": "L06d", "claim": "Euler characteristic of the resolved T4/Z2 from the singular-point count: (chi(T4) - n)/2 + n*chi(P1) (standard formula; FROM MEMORY as a formula)", "tier": "L",
                "arithmetic": f"({chi_T4}-{n_sing})/2 + {n_sing}*2 = {chiK3}", "value": int(chiK3), "expected_from_source": f"Track D export chi = {dexp.get('chi')}", "agrees": int(chiK3) == dexp.get("chi"),
                "note": "cross-track consistency only; no relation of this 24 to the TT tadpole 24 is claimed or tested here. The F-theory relation chi(CY4)/24 = 24 for K3 x K3 (Sethi-Vafa-Witten, FROM MEMORY) is a tier-L statement not checked in this track"})
    b2 = dexp.get("b2")
    lit.append({"id": "L06e", "claim": "TT unfixed Kahler directions of K3: b2 - dim V_flux = 22-4 = 18 (dim V=4, section 3.4) and 22-2 = 20 (dim V=2, section 4.1)", "tier": "L",
                "arithmetic": f"{b2}-4 = {b2 - 4}; {b2}-2 = {b2 - 2}", "value": [b2 - 4, b2 - 2],
                "inputs": {"b2": f"{b2} read from Track D export (rank of H^2), dim V_flux quoted"},
                "expected_from_source": [18, 20], "source_lines": [quote("leaves an 18 dimensional subspace"), quote("allowing for all twenty")],
                "agrees": [b2 - 4, b2 - 2] == [18, 20]})

    # unfixed moduli quotations (line numbers by grep)
    moduli = {
        "abstract": quote("All the complex structure moduli and some of the"),
        "summary_3.4": quote("Unlike the complex structure"),
        "3.4_18": quote("18 dimensional subspace of the K3"),
        "T2_volume_3.4": quote("of the T 2 , unfixed"),
        "overall_volume_2.1": quote("modulus which is left unfixed is the overall volume"),
        "open_string_2.2": quote("there are moduli that arise from the open string sector"),
        "D3_positions_2.4": quote("Turning on flux does not freeze these fields"),
        "D7_positions_2.4": quote("generically the 7-brane moduli"),
        "conclusion_volume": quote("the most serious limitation of these"),
        "4.1_twenty": quote("allowing for all twenty"),
    }
    res = [{"id": "E06a", "quantity": "moduli that remain unfixed, quoted from TT with line numbers (lines found by grep on the pinned text)", "computed": {"quotes": moduli,
            "summary": "ALL complex structure moduli of K3 x T2 and the dilaton-axion are fixed in the first branch (c != 0); of the closed-string Kahler moduli the volume, the T2 Kahler modulus and (22 - dim V_flux) K3 Kahler directions "
                       "(18 for dim V_flux=4, 20 for dim V_flux=2) remain unfixed; open-string moduli (D3 positions, 7-brane gauge fields, 7-brane positions, expected but not shown to be lifted) are not analysed. "
                       "Second branch: see E05c (N=2, only partial fixing)."},
            "shared_inputs": ["TT_H33_gram_A3"], "script": SCRIPT, "command": CMD}]

    # TT section 6: null-vector shifts of the (4.15) example, exact
    a = [1, 1, 0, 0, 0, 0]; b = [0, 0, 1, 1, 0, 0]     # (4.15) lambda_x = e1, lambda_x' = e2 in fg coordinates (e1=f1+g1,e2=f2+g2)
    nullv = [0, 0, 0, 0, 1, 0]                        # f3 : null, orthogonal to a and b
    dot = lambda u, v: sum(u[i] * H[i, j] * v[j] for i in range(6) for j in range(6))
    assert dot(nullv, nullv) == 0 and dot(nullv, a) == 0 and dot(nullv, b) == 0
    rows = []
    for n in range(0, 6):
        ax = [2 * x + 2 * n * y for x, y in zip(a, nullv)]; bx = [2 * x for x in b]; ay = [-x for x in bx]; by = [2 * x for x in a]
        Nf = -dot(bx, ay) + dot(by, ax)
        Gm = gram_of(H, [ax, bx, ay, by]); rk = sp.Matrix([ax, bx, ay, by]).rank()
        sg = sig_rank(Gm)
        rows.append({"n": n, "N_flux": int(Nf), "alpha_x": ax, "coordinate_rank": int(rk), "gram_signature_(pos,neg,null)": sg, "null_vector_in_V_flux": bool(sg[2] > 4 - rk)})
    res.append({"id": "E06b", "quantity": "TT section 6 check: alpha_x -> alpha_x + 2 n v with v null and orthogonal to V_flux: N_flux constant, V_flux acquires a null vector (no spacelike Kahler form)",
                "computed": {"rows": rows, "N_flux_constant": len({r["N_flux"] for r in rows}) == 1, "null_vector_appears_for_n_ne_0": all(r["null_vector_in_V_flux"] for r in rows[1:]) and not rows[0]["null_vector_in_V_flux"],
                             "meaning": "an infinite family of inequivalent flux data with the SAME D3 charge exists, but only n=0 is an allowed vacuum: this is why the O(Gamma)-orbit count of admissible fluxes can be finite while window counts grow"},
                "shared_inputs": ["TT_H33_gram_A3", "sec41_conditions_4.8_4.10", "N_flux_formula_2.8", "flux_quantisation_even_coefficients"], "script": SCRIPT, "command": CMD})
    write_json("06_literature_and_moduli.json", {"script": SCRIPT, "results": res, "rigidity": [], "literature_checks": lit, "could_not_do": []})
    print(json.dumps(lit, indent=1, default=str)[:3500])

if __name__ == "__main__":
    main()
