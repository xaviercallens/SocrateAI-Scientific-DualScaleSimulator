"""
Assembles results.json and exports.json for Track B v3 from the per-script JSON outputs. Nothing is hand-authored:
every value is read from a *_results.json written by the scripts named in the `script`/`command` fields.

Run (after the five scripts): cd audit/k3t2_rigidity_v3/B-dyons && <venv-python> build_results.py
(no arguments)
"""
import json
from fractions import Fraction as Fr
from common import HERE, REPO

PYV = "/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python"
CD = "cd audit/k3t2_rigidity_v3/B-dyons && " + PYV + " "


def J(name):
    with open(HERE / name) as f:
        return json.load(f)


if __name__ == "__main__":
    inputs = J("inputs.json")
    h = J("hurwitz_results.json")
    p1 = J("part1_euler_results.json")
    p2 = J("part2_psi_m_results.json")
    p3 = J("part3_immortal_results.json")
    mg = J("theta_margin_check_results.json")

    s3 = p3["main_solution"]
    N, c, M = s3["solution_N_c_M"]
    g = p2["G_solutions_exact"]
    exp2 = p2["expected_after_computing"]
    solved_equals_expected = {
        "G_2": g["G_2"]["solution"] == exp2["G_2 = (9/4) B^2 + (3/4) E4 A^2  [4 Dpsi_1 = 9B^2/A + 3E4A]"],
        "G_3": g["G_3"]["solution"] == [str(Fr(x)) for x in exp2["G_3 = (50/27) B^3 + (48/27) E4 A^2 B + (10/27) E6 A^3"]],
        "immortal_c_equals_3": c == "3",
    }
    results = [
        {"id": "B1_hurwitz", "quantity": "Hurwitz class numbers H(D), D<=400: reduced-form counting vs Dirichlet L-function route, Kronecker-Hurwitz relation n<=100",
         "computed": f"{h['num_D_values']} values of D; counting == Dirichlet for all: {not h['counting_vs_dirichlet_mismatch_D']}; "
                     f"{h['num_fundamental_discriminants_checked_with_pure_class_number_formula']} fundamental discriminants use the pure formula h=-(1/|D|)sum chi(a)a; "
                     f"KH failures: {h['kronecker_hurwitz_failures']}; wrong-weight control breaks {h['negative_control_wrong_weight_1_for_a_eq_c_b_eq_0']}",
         "shared_inputs": ["Hurwitz_H0", "Kronecker_Hurwitz_relation"], "script": "hurwitz.py", "command": CD + "hurwitz.py 400"},
        {"id": "B2_euler", "quantity": "e(Hilb^k(K3)), k=0..8, from the DMVV product evaluated at y=1 (own expansion)",
         "computed": p1["euler_numbers_from_DMVV_at_y1_q0"],
         "extra": {"equals_Goettsche_with_chi_from_TrackD": p1["product_equals_goettsche"],
                   "chi_TrackD": p1["trackD_exports_source"]["chi"], "TrackD_source_file": p1["trackD_exports_source"]["path"],
                   "q_positive_coeffs_vanish": p1["q^j (j>=1) coefficients of p^k, k<=KMAX (must be empty)"] == {},
                   "k_solved_from_chi_D": p1["k_solved_exactly_from_chi_D = chi_D/cB_sum"],
                   "declared_k_equals_solved_k (a consistency check between typed k=2 and Track D chi; NOT an independent derivation of 24: 24 = k x 12)": p1["declared_k_equals_solved_k"]},
         "shared_inputs": p1["shared_inputs"], "script": "part1_euler.py", "command": CD + "part1_euler.py 40 8 4"},
        {"id": "B3_G_solutions", "quantity": "coefficients of G_1..G_4 in weak-Jacobi monomial basis, solved exactly (Delta psi_m = G_{m+1}/A, m=0..3)",
         "computed": {k: {"monomials": v["monomials"], "solution": v["solution"], "rank_M": v["rank_M"], "rank_augmented": v["rank_augmented"],
                          "num_equations": v["num_equations"], "residual_rows": v["num_nonzero_residual_equations"],
                          "all_neighbours_fail": v["all_neighbours_fail"]} for k, v in g.items()},
         "extra": {"solved_equals_expected_after_computing (DMZ 5.16, FROM MEMORY)": solved_equals_expected},
         "shared_inputs": p2["shared_inputs"], "script": "part2_psi_m.py", "command": CD + "part2_psi_m.py 40 7 2"},
        {"id": "B4_immortal", "quantity": "(N,c,M) in Delta psi_1 = N A_{2,1} + c E4 A - M Hhat, exact linear solve",
         "computed": {"N": N, "c": c, "M": M, "num_rows": s3["num_rows"], "rank_M": s3["rank_M"], "rank_augmented": s3["rank_augmented"],
                      "residual_rows": s3["num_nonzero_residual_rows"], "all_26_integer_neighbours_fail": s3["all_26_neighbours_fail"],
                      "k_scan_solutions": {k: v["solution_N_c_M"] for k, v in p3["k_scan_1_to_6"].items()}},
         "shared_inputs": p3["shared_inputs"], "script": "part3_immortal.py", "command": CD + "part3_immortal.py 40 8 20 60"},
        {"id": "B5_margin", "quantity": "theta-cache truncation stability, QMAX 40 vs 60", "computed": mg["agree"],
         "shared_inputs": ["Jacobi_theta_definitions"], "script": "theta_margin_check.py", "command": CD + "theta_margin_check.py 40 60 8"},
    ]
    rigidity = [
        {"parameter_inserted": "overall factor k in Z_K3 = k*phi_{0,1}", "selecting_condition": "k*(cB(0)+2cB(-1)) = chi from Track D (contains the target chi)",
         "condition_uses_true_value": True, "solution_set": "{2} (exact solve k = chi/12)", "classification": "NORMALISATION",
         "input_it_depends_on": "Track D chi; typed k=2 is a declared input",
         "negative_control": "k in 1..6 (all integers) each change e(Hilb^k); values other than 2 are HYPOTHETICAL (no realizable object), so not counted as discrimination",
         "control_perturbs_same_parameter": True},
        {"parameter_inserted": "coefficients x of G_2,G_3,G_4 in the weak-Jacobi monomial basis",
         "selecting_condition": "linear identity for all (n,l), n<=QCHK (target-free)", "condition_uses_true_value": False,
         "solution_set": "unique (rank = #unknowns, zero residual) for k=2", "classification": "RIGID_GIVEN_DEFINITION",
         "input_it_depends_on": "K3_elliptic_genus_factor_k, weak_Jacobi_ring_monomials",
         "negative_control": "every neighbour x_i +- 1/d fails (all_neighbours_fail true for G_1..G_4)", "control_perturbs_same_parameter": True},
        {"parameter_inserted": "(N,c,M) immortal decomposition", "selecting_condition": "exact identity of series for all (n,l) in window (target-free)",
         "condition_uses_true_value": False, "solution_set": f"unique ({N},{c},{M}) for k=2; k-dependent for other k (scan)",
         "classification": "RIGID_GIVEN_DEFINITION", "input_it_depends_on": "immortal_ansatz_form, K3_elliptic_genus_factor_k",
         "negative_control": "all 26 joint integer neighbours fail", "control_perturbs_same_parameter": True},
        {"parameter_inserted": "k (again) tested against solvability of the G_2 and immortal systems", "selecting_condition": "solvability of the linear systems",
         "condition_uses_true_value": False, "solution_set": "every k in 1..6 admits a solution", "classification": "NON_DISCRIMINATING",
         "input_it_depends_on": "-", "negative_control": "k=1..6, all consistent", "control_perturbs_same_parameter": True},
        {"parameter_inserted": "Hurwitz values H(D)", "selecting_condition": "two independent computations agree; Kronecker-Hurwitz identity",
         "condition_uses_true_value": False, "solution_set": "counted == Dirichlet for all D<=400", "classification": "VERIFIED_IDENTITY",
         "input_it_depends_on": "Hurwitz_H0", "negative_control": "wrong weight breaks both", "control_perturbs_same_parameter": True},
    ]
    could_not_do = [
        "Delta*psi_3 (G_4) is solved as coefficients in a monomial basis only; no literature line was available to compare (no DMZ formula for m=3 quoted), and the mock/polar decomposition for m=3 was not attempted.",
        "The immortal decomposition for m>=2 (needs A_{2,m} plus extra mock pieces) was not done.",
        "Track D exports.json was read from the v3 D-tda directory as found at run time (a parallel track, status: depends on declared local-model inputs); rerun part1_euler.py if it changes.",
        "k_norm=2 is a typed input; the check chi=k*12=24 against Track D is a cross-check, not an independent derivation of 24.",
        "Truncation margins were checked only by QMAX 40 vs 60 agreement (theta_margin_check.py); no rigorous bound on the margin.",
    ]
    out = {"track": "B-dyons (v3)", "branch": "loop/k3t2-rigidity", "inputs_file": "audit/k3t2_rigidity_v3/B-dyons/inputs.json",
           "inputs": inputs, "results": results, "rigidity": rigidity, "could_not_do": could_not_do}
    with open(HERE / "results.json", "w") as f:
        json.dump(out, f, indent=1)
    exports = {
        "chi_from_own_series": p1["chi_from_own_series = k*cB_sum"],
        "euler_numbers_Hilb_k_K3_k0_to_8": p1["euler_numbers_from_DMVV_at_y1_q0"],
        "G_solutions": {k: v["solution"] for k, v in g.items()},
        "immortal_N_c_M": [N, c, M],
        "hurwitz_table_path": "audit/k3t2_rigidity_v3/B-dyons/hurwitz_results.json",
        "depends_on_inputs": ["K3_elliptic_genus_factor_k", "Track_D_chi_K3"],
        "source": "build_results.py from part*_results.json",
    }
    with open(HERE / "exports.json", "w") as f:
        json.dump(exports, f, indent=1)
    print(json.dumps({"solved_equals_expected": solved_equals_expected, "exports": exports}, indent=1)[:1500])
