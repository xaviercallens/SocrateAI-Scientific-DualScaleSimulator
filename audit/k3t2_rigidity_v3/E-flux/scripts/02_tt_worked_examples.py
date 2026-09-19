#!/usr/bin/env python3
"""E-flux v3, script 02: TT's own worked examples, recomputed exactly, decide the (2.3) 1/2 vs line ~941.

Run: cd audit/k3t2_rigidity_v3/E-flux/scripts && <venv python> 02_tt_worked_examples.py
(no arguments)

Every example is computed in (A.6) e-coordinates AND in the unimodular U^3 presentation through the embedding
e_i=f_i+g_i, e_{i+3}=f_i-g_i (same numbers: it is an isometry).  TT's printed values are 'expected' fields,
filled after computing, each with the line where it is printed (found by grep, not typed).

Tadpole: N_D3 = 24 - N_flux/2 (2.3: the 24 and the 1/2 are declared tier-L inputs).  The alternative reading
(no 1/2) is computed as a counterfactual to see which reading reproduces TT's printed numbers.
"""
import itertools
import sympy as sp
from sympy import I, Rational as R, sqrt
from common import *
from lattice_u3 import H33
from susy import *

SCRIPT = "02_tt_worked_examples.py"
CMD = rel_command(SCRIPT)
H = sp.Matrix(H33)
A6 = sp.diag(2, 2, 2, -2, -2, -2)
T = sp.zeros(6, 6)
for i in range(3):
    T[i, 2 * i] = 1; T[i, 2 * i + 1] = 1
    T[i + 3, 2 * i] = 1; T[i + 3, 2 * i + 1] = -1

def e2u(c):
    """vector given by coefficients c_j of e_j (j=1..6) -> (f,g)-coordinates"""
    return list(sp.Matrix([c]) * T)

def nflux_28(ax, bx, ay, by, Hm):
    """(2.8): N_flux = -bx.ay + by.ax"""
    return -cdot(bx, ay, Hm) + cdot(by, ax, Hm)

def ND3(Nf, half=True):
    return 24 - (Nf * R(1, 2) if half else Nf)

def main():
    res = []
    lit = []
    ex = {}
    # TT printed values (expected), with the lines where they are printed
    ex["4.16_ND3"] = (16, quote("ND3 = 24 − 8 = 16"))
    ex["4.32_Nflux"] = (32, quote("Nf lux = αx .βy − βx .αy = 32"))
    ex["4.32_D3_added"] = (8, quote("As a result 8 D3 branes need to be added"))
    ex["5.3_D3_added"] = (20, quote("(2.3), so that 20"))  # 20 D3 printed at end of this line
    ex["5.3_halfNflux"] = (4, quote("three brane charge is Nf lux /2 = 4"))
    ex["941_remark"] = quote("choices of δαx , δβx which give rise to a vacuum where")
    ex["2.3"] = quote("Nf lux + ND3 = 24")
    ex["4.13"] = quote("Nf lux = 2α2x")
    ex["4.14"] = quote("α2x + ND3 = 24")
    ex["4.15"] = quote("αx = 2e1 , βx = 2e2")
    ex["4.31_alpha_yy"] = quote("For these flux vectors αxx = 16 and αyy = 8")

    # ---------- (4.15): section 4.1 example
    ax = e2u([2, 0, 0, 0, 0, 0]); bx = e2u([0, 2, 0, 0, 0, 0]); ay = [-x for x in bx]; by = ax
    a2 = cdot(ax, ax, H); Nf = nflux_28(ax, bx, ay, by, H)
    cond = tt_conditions(H, ax, bx, ay, by, I, I)
    Gm = gram_of(H, [ax, bx, ay, by]); sg = sig_rank(gram_of(H, [ax, bx]))
    res.append({"id": "E02a", "quantity": "TT (4.15)-(4.16): alpha_x=2e1, beta_x=2e2 with (4.8), phi=tau=i",
        "computed": {"alpha_x_sq": int(a2), "beta_x_sq": int(cdot(bx, bx, H)), "alpha_x_dot_beta_x": int(cdot(ax, bx, H)),
                     "N_flux_from_2.8": int(Nf), "N_flux_equals_2_alpha_sq_(4.13)": bool(Nf == 2 * a2),
                     "V_flux_signature_(pos,neg,null)": sg,
                     "N_D3_with_half_(2.3)": int(ND3(Nf)), "N_D3_without_half_counterfactual": int(ND3(Nf, False)),
                     "TT_printed_N_D3": ex["4.16_ND3"][0], "printed_at": ex["4.16_ND3"][1],
                     "reproduces_TT_with_half": bool(ND3(Nf) == ex["4.16_ND3"][0]),
                     "reproduces_TT_without_half": bool(ND3(Nf, False) == ex["4.16_ND3"][0]),
                     "susy_conditions_3.16abc_at_phi=tau=i": {k: str(cond[k]) for k in ("3.16a", "3.16b", "3.16c")},
                     "3.31_value_positive": str(cond["3.15_3.31_value"])},
        "shared_inputs": ["TT_A6_gram_A6_A7", "embedding_A6_into_U3", "tadpole_total_24_eq2.3", "tadpole_half_factor_eq2.3", "N_flux_formula_2.8", "sec41_conditions_4.8_4.10"],
        "script": SCRIPT, "command": CMD})
    # ---------- (4.31): section 4.2 example
    ax = e2u([2, -2, 0, 0, 0, 0]); ay = e2u([2, 2, 0, 2, 0, 0]); bx = e2u([0, -4, 0, 0, 0, 0]); by = e2u([4, 0, 0, 2, 2, 0])
    axx = cdot(ax, ax, H); ayy = cdot(ay, ay, H)
    c424 = [cdot(ax, ay, H), cdot(bx, by, H), cdot(ax, by, H) + cdot(ay, bx, H)]
    c425 = [(cdot(bx, bx, H), 2 * cdot(ax, ax, H), 2 * cdot(ax, bx, H)), (cdot(by, by, H), 2 * cdot(ay, ay, H), 2 * cdot(ay, by, H))]
    Nf = nflux_28(ax, bx, ay, by, H)
    Nf_alt = cdot(ax, by, H) - cdot(bx, ay, H)   # TT (4.32) written as alpha_x.beta_y - beta_x.alpha_y
    tau = I * sqrt(R(ayy, axx)); phis = [R(1, 2) * (1 + I), R(1, 2) * (1 - I)]
    conds = [tt_conditions(H, ax, bx, ay, by, ph, tau) for ph in phis]
    sig4 = sig_rank(gram_of(H, [ax, bx, ay, by]))
    # orbifold criterion with phi=(1+i)/2 (Im phi > 0)
    Om = conds[0]["Gzbar"]
    orb, Bsat, nker = orbifold_hits(H, [ax, bx, ay, by], Om, 6)
    res.append({"id": "E02b", "quantity": "TT (4.31)-(4.32): section 4.2 example, (2+,2-)",
        "computed": {"alpha_xx": int(axx), "alpha_yy": int(ayy), "(4.24)_three_products": [int(x) for x in c424],
                     "(4.25)_triples_equal": [[int(t) for t in tr] for tr in c425],
                     "N_flux": int(Nf), "N_flux_TT_form_ax.by-bx.ay": int(Nf_alt),
                     "N_D3_with_half": int(ND3(Nf)), "N_D3_without_half_counterfactual": int(ND3(Nf, False)),
                     "TT_printed_N_flux": ex["4.32_Nflux"][0], "TT_printed_D3_added": ex["4.32_D3_added"][0],
                     "printed_at": [ex["4.32_Nflux"][1], ex["4.32_D3_added"][1]],
                     "reproduces_TT_with_half": bool(ND3(Nf) == ex["4.32_D3_added"][0] and Nf == ex["4.32_Nflux"][0]),
                     "reproduces_TT_without_half": bool(ND3(Nf, False) == ex["4.32_D3_added"][0]),
                     "tau_(3.29-4.30)": str(tau), "phi_choices": [str(p) for p in phis],
                     "susy_3.16abc_for_phi_(1+i)/2": {k: str(conds[0][k]) for k in ("3.16a", "3.16b", "3.16c")},
                     "susy_3.16abc_for_phi_(1-i)/2": {k: str(conds[1][k]) for k in ("3.16a", "3.16b", "3.16c")},
                     "3.31_value": str(conds[0]["3.15_3.31_value"]),
                     "V_flux_signature": sig4,
                     "orbifold_test_(lattice vector in V_flux orthogonal to Omega)": {"exists": bool(orb), "dim_of_solution_space": int(nker), "saturation_rank": len(Bsat)}},
        "shared_inputs": ["TT_A6_gram_A6_A7", "embedding_A6_into_U3", "tadpole_total_24_eq2.3", "tadpole_half_factor_eq2.3", "N_flux_formula_2.8", "sec42_conditions_4.24_4.25", "orbifold_criterion_TT_3.3"],
        "script": SCRIPT, "command": CMD})
    # ---------- (5.8): second-branch example
    ax0 = [0] * 6; bx0 = e2u([0, 0, 0, 2, 0, 0]); ay0 = e2u([0, 0, 0, 2, 0, 0]); by0 = [0] * 6
    Nf = nflux_28(ax0, bx0, ay0, by0, H)
    res.append({"id": "E02c", "quantity": "TT section 5.3, (5.8): alpha_x=0, beta_x=2e4, alpha_y=2e4, beta_y=0 (second branch)",
        "computed": {"N_flux": int(Nf), "N_flux_over_2": str(Nf * R(1, 2)), "e4_norm": int(cdot(e2u([0, 0, 0, 1, 0, 0]), e2u([0, 0, 0, 1, 0, 0]), H)),
                     "N_D3_with_half": int(ND3(Nf)), "N_D3_without_half_counterfactual": int(ND3(Nf, False)),
                     "TT_printed_half_N_flux": ex["5.3_halfNflux"][0], "TT_printed_D3": ex["5.3_D3_added"][0],
                     "printed_at": [ex["5.3_halfNflux"][1], ex["5.3_D3_added"][1]],
                     "reproduces_TT_with_half": bool(ND3(Nf) == ex["5.3_D3_added"][0] and Nf * R(1, 2) == ex["5.3_halfNflux"][0]),
                     "reproduces_TT_without_half": bool(ND3(Nf, False) == ex["5.3_D3_added"][0])},
        "shared_inputs": ["TT_A6_gram_A6_A7", "embedding_A6_into_U3", "tadpole_total_24_eq2.3", "tadpole_half_factor_eq2.3", "N_flux_formula_2.8", "branch2_conditions_5.1_5.3"],
        "script": SCRIPT, "command": CMD})

    # ---------- line ~941: N_flux needed for N_D3 = 0, and explicit constructions
    Nf_needed = sp.solve(sp.Eq(24 - sp.Symbol("N") * R(1, 2), 0), sp.Symbol("N"))[0]
    Nf_needed_nohalf = sp.solve(sp.Eq(24 - sp.Symbol("N"), 0), sp.Symbol("N"))[0]
    # explicit section-4.1 flux with alpha_x^2 = 24 (lambda^2 = 6): search small lambda in U3
    found = None
    rng = range(-2, 3)
    vecs6 = [v for v in itertools.product(rng, repeat=6) if sum(2 * v[2 * i] * v[2 * i + 1] for i in range(3)) == 6]
    for a in vecs6:
        for b in vecs6:
            if sum(a[2 * i] * b[2 * i + 1] + a[2 * i + 1] * b[2 * i] for i in range(3)) == 0:
                found = (a, b); break
        if found: break
    a, b = found
    axf = [2 * x for x in a]; bxf = [2 * x for x in b]; ayf = [-x for x in bxf]; byf = axf
    Nf41 = nflux_28(axf, bxf, ayf, byf, H)
    # the delta construction of section 4.1 (4.17)-(4.23) from (4.15): enumerate timelike delta_ax, delta_bx (even coefficients), orthogonal
    # to alpha_x, beta_x and to each other; report the achievable N_D3 values
    ax = e2u([2, 0, 0, 0, 0, 0]); bx = e2u([0, 2, 0, 0, 0, 0])
    Nf0 = 16
    dvals = {}
    box = range(-1, 2)  # lambda-window in e-coordinates for the delta vectors, delta = 2*lambda
    lam = [c for c in itertools.product(box, repeat=6) if c[0] == 0 and c[1] == 0]   # orthogonal to e1,e2 (diagonal metric)
    def n_e(c):
        return 2 * sum(s * x * x for s, x in zip((1, 1, 1, -1, -1, -1), c))
    timelike = [c for c in lam if n_e(c) < 0]
    for c1 in timelike:
        for c2 in timelike:
            dot12 = 2 * sum(s * x * y for s, x, y in zip((1, 1, 1, -1, -1, -1), c1, c2))
            dN = -(4 * n_e(c2)) - (4 * n_e(c1))       # -(delta bx)^2 - (delta ax)^2 with delta = 2*lambda
            dvals.setdefault(int(Nf0 + dN), 0)
            dvals[int(Nf0 + dN)] += 1
    ND3_by_Nflux = {k: int(ND3(k)) for k in sorted(dvals)}
    res.append({"id": "E02d", "quantity": "line ~941 remark 'N_flux = 24 and no D3-branes': which N_flux gives N_D3 = 0 under (2.3)?",
        "computed": {"N_flux_for_N_D3=0_with_half_(2.3)": str(Nf_needed), "N_flux_for_N_D3=0_without_half": str(Nf_needed_nohalf),
                     "remark_printed_at": ex["941_remark"], "eq_2.3_printed_at": ex["2.3"],
                     "explicit_section41_flux_in_U3": {"lambda_x": a, "lambda_y": b, "alpha_x_sq": int(cdot(axf, axf, H)), "N_flux": int(Nf41), "N_D3_with_half": int(ND3(Nf41))},
                     "delta_construction_from_4.15_(4.20-4.23)_window_[-1,1]_e-coords": {"N_flux_values_reached_with_counts": {str(k): v for k, v in sorted(dvals.items())}, "N_D3_for_each_(with_half)": {str(k): v for k, v in ND3_by_Nflux.items()}},
                     "N_flux_24_gives_N_D3_(with_half)": int(ND3(24)),
                     "tally_of_TT_worked_examples_reproduced": {"with_half": [r["id"] for r in res[:3] if r["computed"]["reproduces_TT_with_half"]], "without_half": [r["id"] for r in res[:3] if r["computed"]["reproduces_TT_without_half"]]},
                     "verdict": "All numbered equations and all three printed worked examples (4.16), (4.32), (5.3) reproduce with (2.3) as printed (1/2 N_flux + N_D3 = 24). "
                                "N_D3 = 0 needs N_flux = 48 (attained here, e.g. alpha_x^2 = 24); N_flux = 24 leaves N_D3 = 12. The line-941 prose reads consistently only "
                                "if N_flux there means N_flux/2 (a 'downstairs' normalisation, cf. the covering-manifold remark after (2.4)); TT use N_flux/2 = 4 explicitly in 5.3. "
                                "Reported as a convention slip in a prose aside, not as an error in any numbered formula."},
        "shared_inputs": ["tadpole_total_24_eq2.3", "tadpole_half_factor_eq2.3", "N_flux_formula_2.8", "sec41_conditions_4.8_4.10"], "script": SCRIPT, "command": CMD})
    write_json("02_tt_worked_examples.json", {"script": SCRIPT, "results": res, "rigidity": [
        {"parameter_inserted": "the factor 1/2 in eq (2.3) (versus no factor)", "selecting_condition": "agreement with TT's printed D3 counts in three independent worked examples (4.16), (4.32), (5.3)",
         "condition_uses_true_value": False, "solution_set": "with 1/2: all three reproduce; without 1/2: none reproduce",
         "classification": "CONDITIONAL_ON_INPUT", "input_it_depends_on": "TT's printed numbers (4.16),(4.32),(5.3) as literature inputs; the 1/2 is tier L",
         "negative_control": "the same three examples with N_D3 = 24 - N_flux (no 1/2) fail all three", "control_perturbs_same_parameter": True}],
        "literature_checks": lit, "could_not_do": []})
    print(json.dumps(res, indent=1, default=str)[:7000])

if __name__ == "__main__":
    main()
