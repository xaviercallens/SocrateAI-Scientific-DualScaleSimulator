#!/usr/bin/env python3
"""E-flux v3, script 01: TT's unimodular presentation vs the (A.6) diagonal basis; group G; evenness.

Run: cd audit/k3t2_rigidity_v3/E-flux/scripts && <venv python> 01_lattice_index.py
(no arguments)

Defect fixed (v2 #1): v2 refused the U-basis claiming e1 would be null.  That confuses the VECTOR e1 of (A.6)
with the first basis vector of U.  Here e_i = f_i + g_i, e_{i+3} = f_i - g_i, embedded in U^3, reproduce (A.6)
exactly; |det| of the 6x6 transition matrix is the index of span(A.6) in Gamma_{3,3}.
"""
import itertools
import numpy as np
import sympy as sp
from fractions import Fraction
from common import *
from lattice_u3 import *

SCRIPT = "01_lattice_index.py"
CMD = rel_command(SCRIPT)

def signature_exact(G):
    """Exact signature (p, q, z) of a symmetric integer matrix from its characteristic polynomial roots' signs
    (via sympy's exact real root isolation)."""
    M = sp.Matrix(G)
    x = sp.symbols("x")
    cp = sp.Poly(M.charpoly(x).as_expr(), x)
    roots = sp.roots(cp, multiple=True) if cp.degree() <= 6 else None
    if roots is None or len(roots) != cp.degree():
        # fallback: count real roots by sign using Sturm sequences
        p = sp.real_roots(cp)
        pos = sum(1 for r in p if r.evalf(50) > 0); neg = sum(1 for r in p if r.evalf(50) < 0)
        zero = sum(1 for r in p if r == 0)
        return pos, neg, zero
    pos = sum(1 for r in roots if sp.N(r, 50) > 0); neg = sum(1 for r in roots if sp.N(r, 50) < 0)
    zero = sum(1 for r in roots if r == 0)
    return pos, neg, zero

def main():
    res, rig, lit = [], [], []
    H = sp.Matrix(H33)
    E8 = sp.Matrix(parse_E8())
    # ---- unimodularity / evenness / signature of Gamma_{3,3} and full Gamma_{3,19}
    det_H = H.det()
    even_diag_H = all(H33[i][i] % 2 == 0 for i in range(6))
    sig_H = signature_exact(H33)
    E8_minors = [E8[:k, :k].det() for k in range(1, 9)]
    full = sp.diag(H, -E8, -E8)
    det_full = full.det()
    even_diag_full = all(full[i, i] % 2 == 0 for i in range(full.shape[0]))
    # signature of full: H33 has (3,3); -E8 twice is negative definite (leading minors of E8 all > 0)
    sig_full = (sig_H[0], sig_H[1] + 16, 0)
    res.append({
        "id": "E01a", "quantity": "Gamma_{3,3} (A.3) and full Gamma_{3,19} (A.2): det, evenness, signature",
        "computed": {"det_H33": int(det_H), "H33_even_diagonal": even_diag_H, "H33_signature": sig_H,
                     "E8_leading_principal_minors": [int(m) for m in E8_minors],
                     "E8_positive_definite_by_Sylvester": all(m > 0 for m in E8_minors),
                     "det_full_A2": int(det_full), "full_even_diagonal": even_diag_full, "full_signature": sig_full,
                     "note": "A symmetric integer matrix with even diagonal has even q(x)=x^T G x for all integer x; "
                             "so evenness of the FULL lattice is verified by exact computation (tier B) from the parsed Gram."},
        "shared_inputs": ["TT_H33_gram_A3", "TT_E8_cartan_A4"], "script": SCRIPT, "command": CMD})

    # ---- (A.6) diagonal basis and the embedding into U^3
    A6 = sp.diag(2, 2, 2, -2, -2, -2)
    T = sp.zeros(6, 6)   # rows: e_j in (f1,g1,f2,g2,f3,g3) coordinates
    for i in range(3):
        T[i, 2 * i] = 1;     T[i, 2 * i + 1] = 1      # e_{i+1} = f_i + g_i
        T[i + 3, 2 * i] = 1; T[i + 3, 2 * i + 1] = -1  # e_{i+4} = f_i - g_i
    gram_check = T * H * T.T
    index_witness = abs(T.det())
    index_from_discs = sp.sqrt(abs(A6.det()) / abs(H.det()))
    res.append({
        "id": "E01b", "quantity": "index of span(TT (A.6) basis) in unimodular Gamma_{3,3}",
        "computed": {"T_H_Ttranspose_equals_A6": bool(gram_check == A6), "abs_det_A6": int(abs(A6.det())),
                     "index_from_transition_matrix_det": int(index_witness),
                     "index_from_sqrt_disc_ratio": str(index_from_discs),
                     "two_routes_agree": bool(index_witness == index_from_discs),
                     "transition_matrix_rows_e1_to_e6_in_fg_coords": T.tolist(),
                     "consequence": "index > 1 => TT's (A.6) forms span a proper SUBLATTICE D of Gamma_{3,3}; "
                                    "any count restricted to fluxes 2*D is a count of a SUBFAMILY of the unimodular count. "
                                    "For ANY isometric embedding into a unimodular lattice |det T|^2 * |det H| = |det A6| forces the same index."},
        "shared_inputs": ["TT_H33_gram_A3", "TT_A6_gram_A6_A7", "embedding_A6_into_U3"], "script": SCRIPT, "command": CMD})
    # D is not isometric to U(2)^3 : norm/2 parity
    D_has_norm_2_mod_4 = True  # e1 has norm 2
    U2_all_norms_div4 = all((2 * 2 * (x * y + z * w + u * v)) % 4 == 0 for x, y, z, w, u, v in itertools.product(range(-2, 3), repeat=6))
    res.append({
        "id": "E01c", "quantity": "is span(A.6) isometric to U(2)^3 (the pullback lattice)? parity-of-norm/2 test",
        "computed": {"D_contains_vector_of_norm_2": True, "norm_of_e1": int(A6[0, 0]),
                     "U(2)^3_all_norms_divisible_by_4_on_window": U2_all_norms_div4,
                     "D_isometric_to_U2_cubed": False,
                     "note": "D = 2*I_{3,3} (odd unimodular scaled), U(2)^3 = 2*U^3 (even unimodular scaled): different Z-lattices "
                             "of the same rank/signature/determinant -64. So TT's (A.6) integral basis is NOT U(2)^3; we make no claim on which is "
                             "geometrically correct, only that the two differ as Z-lattices (tier B)."},
        "shared_inputs": ["TT_A6_gram_A6_A7"], "script": SCRIPT, "command": CMD})

    # ---- group G
    Gm = group_G()
    n_iso = sum(is_isometry(M) for M in Gm)
    dets = [det_on_positive_plane(M) for M in Gm]
    Gplus = [M for M, d in zip(Gm, dets) if d == 1]
    preserves_D = all(in_D(tuple(M @ np.array(v))) for M in Gm for v in [(1, 1, 0, 0, 0, 0), (1, -1, 0, 0, 0, 0), (0, 0, 2, 0, 0, 0), (0, 0, 0, 2, 0, 0)])
    # G-invariance of D checked on generators of D (rows of T) for every M
    Trows = [tuple(int(t) for t in T.row(i)) for i in range(6)]
    preserves_D = all(in_D(tuple(int(z) for z in (M @ np.array(v)))) for M in Gm for v in Trows)
    distinct = len({M.tobytes() for M in Gm})
    res.append({
        "id": "E01d", "quantity": "finite symmetry group G (order, isometry check, orientation of positive 3-block)",
        "computed": {"name": "(Z2 x Z2) wr S3 = block permutations x per-block swap(x<->y) x per-block sign, acting on U+U+U",
                     "order_G": len(Gm), "all_distinct": distinct == len(Gm), "n_isometries_g^T H g == H": int(n_iso),
                     "orientation_det_on_positive_plane_counts": {str(k): sum(1 for d in dets if d == k) for k in (1, -1)},
                     "order_G_plus_(det_pos=+1)": len(Gplus),
                     "G_preserves_the_sublattice_D_(A.6)": preserves_D,
                     "justification": "(a) every element is an integer matrix with g^T H g = H, i.e. an element of O(Gamma_{3,3}) extended by the identity "
                                      "on the rest of Gamma_{3,19}, hence in O(Gamma_{3,19}); (b) the SUSY/tadpole conditions used (4.8)-(4.10), (2.8), "
                                      "(3.10)-(3.15), (3.31) are polynomial in inner products of the flux vectors, so they are preserved when "
                                      "all four flux vectors are moved by an isometry (checked symbolically in 03/04); (c) TT quotient the moduli by O+(3,19) "
                                      "(2.13),(2.16),(5.4): Omega is an ORIENTED spacelike plane and orientation reversal is Omega<->Omegabar; "
                                      "G+ = elements with det +1 on the positive 3-block (every element preserves the positive/negative splitting, asserted).",
                     "not_quotiented": "S-duality (H3<->F3), T-duality/SL(2,Z) on T2 and the full infinite O(Gamma) are NOT quotiented; "
                                       "counts 'modulo G' are counts of G-orbits inside the stated window, nothing more."},
        "shared_inputs": ["TT_H33_gram_A3", "group_G_definition"], "script": SCRIPT, "command": CMD})

    # ---- v2's coordinate group (perm within blocks of the e_i, sign flips): how many elements preserve Gamma?
    Tinv = T.inv()
    pos = [0, 1, 2]; neg = [3, 4, 5]
    n_total = 0; n_pres = 0; bad_example = None
    for pp in itertools.permutations(pos):
        for pn in itertools.permutations(neg):
            for fl in itertools.product((1, -1), repeat=6):
                n_total += 1
                Pm = sp.zeros(6, 6)  # acts on e-coordinates (row-vector convention) : e_j -> fl_j * e_{perm(j)}
                perm = list(pp) + list(pn)
                for j in range(6):
                    Pm[j, perm[j]] = fl[j]
                # in fg coordinates: row-vector map  v -> v * (Tinv * Pm * T) ... rows of T are e_j
                Mfg = Tinv * Pm * T
                if all(x.q == 1 for x in Mfg):
                    n_pres += 1
                elif bad_example is None:
                    bad_example = {"perm": perm, "signs": list(fl), "matrix_in_fg_coords": [[str(x) for x in row] for row in Mfg.tolist()]}
    res.append({
        "id": "E01e", "quantity": "v2's group (perm within pos/neg block of the e_i, sign flips; order 2304): how many elements map Gamma_{3,3} to itself",
        "computed": {"order_v2_group": n_total, "n_elements_integral_in_U3_coordinates": n_pres, "example_not_in_O(Gamma)": bad_example,
                     "consequence": "v2's group is NOT a subgroup of O(Gamma_{3,3}); only the subgroup counted here preserves the unimodular lattice. "
                                    "v3 uses G (order in E01d), verified to be a subgroup of O(Gamma)."},
        "shared_inputs": ["TT_H33_gram_A3", "TT_A6_gram_A6_A7", "embedding_A6_into_U3"], "script": SCRIPT, "command": CMD})

    # ---- divisibility alpha^2 in 8Z : unimodular even lattice, and odd-lattice control
    # exact identity: for lambda in Gamma, (2 lambda)^2 = 4 q(lambda), q(lambda) even
    box = range(-2, 3)
    vals_even = set(); vals_odd = set()
    Hodd = np.diag([1, 1, 1, -1, -1, -1])
    Hn = np.array(H33, dtype=np.int64)
    n_checked = 0
    for v in itertools.product(box, repeat=6):
        va = np.array(v)
        qe = int(va @ Hn @ va); qo = int(va @ Hodd @ va)
        vals_even.add((4 * qe) % 8); vals_odd.add((4 * qo) % 8); n_checked += 1
    res.append({
        "id": "E01e2", "quantity": "alpha_x^2 mod 8 for alpha_x = 2*lambda, lambda in Gamma: unimodular EVEN lattice vs ODD unimodular control",
        "computed": {"lambdas_checked_in_window_[-2,2]^6": n_checked,
                     "residues_alpha_sq_mod_8_even_lattice_U3": sorted(vals_even),
                     "residues_alpha_sq_mod_8_odd_control_I33": sorted(vals_odd),
                     "structural_reason_even_case": "even Gram diagonal => q(lambda) in 2Z => (2 lambda)^2 = 4 q in 8Z (all of Gamma_{3,19}, E01a)",
                     "control_meaning": "the odd lattice <1>^3+<-1>^3 is unimodular with the SAME signature; it differs only in evenness, the property the divisibility depends on. "
                                        "It is realizable as a lattice but is NOT H^2(K3,Z); it is a control for the logic, HYPOTHETICAL as a physical spectrum."},
        "shared_inputs": ["TT_H33_gram_A3", "flux_quantisation_even_coefficients", "TT_E8_cartan_A4"], "script": SCRIPT, "command": CMD})
    rig.append({
        "parameter_inserted": "the factor 8 in 'alpha_x^2 in 8Z' (hence the admissible set {8,16,24} of alpha_x^2 <= 24)",
        "selecting_condition": "Gamma even (Gram diagonal even) and flux quantisation alpha = 2*lambda, lambda in Gamma; the target 8 does not appear in the condition",
        "condition_uses_true_value": False,
        "solution_set": "alpha_x^2 in 8Z for every lambda (even lattice); odd unimodular control gives alpha_x^2 in 4Z (residue 4 mod 8 occurs)",
        "classification": "RIGID_GIVEN_DEFINITION",
        "input_it_depends_on": "flux_quantisation_even_coefficients (TT eq 2.5) and evenness of H^2(K3,Z) (from A.2, checked in E01a)",
        "negative_control": "odd unimodular lattice I_{3,3}: residues mod 8 of alpha^2 include 4 (E01e2) -- perturbs the SAME property (evenness) the 8 depends on",
        "control_perturbs_same_parameter": True,
        "note": "The admissible set {alpha_x^2 <= 24} = {8,16,24} is the same in the unimodular lattice and in the subfamily 2D: the index changes counts, not the divisibility."})
    alt = sp.diag(2, 2, 2, -2, -2, -4)
    rig.append({"parameter_inserted": "the index 8 of span(A.6) in Gamma_{3,3}", "selecting_condition": "|det T|^2 * |det H33| = |det A6| for any isometric embedding of the A.6 Gram into a unimodular lattice (identity; nothing is chosen)",
                "condition_uses_true_value": False, "solution_set": "index = sqrt(|det A6|/|det H33|) = %s" % str(index_from_discs), "classification": "VERIFIED_IDENTITY",
                "input_it_depends_on": "TT_A6_gram_A6_A7 (the determinant 64) -- a typed literature Gram", "negative_control": "perturbing the Gram to diag(2,2,2,-2,-2,-4) (det %d) gives sqrt(%d) = %s: not an integer, so no finite-index embedding into a unimodular lattice exists" % (int(abs(alt.det())), int(abs(alt.det())), str(sp.sqrt(abs(alt.det())))),
                "control_perturbs_same_parameter": True})
    write_json("01_lattice_index.json", {"script": SCRIPT, "results": res, "rigidity": rig, "literature_checks": lit, "could_not_do": [],
        "provenance": {"tt_txt_sha256": sha256(TT_TXT), "tt_txt": rel(TT_TXT)}})
    print(json.dumps(res, indent=1, default=str)[:6000])

if __name__ == "__main__":
    main()
