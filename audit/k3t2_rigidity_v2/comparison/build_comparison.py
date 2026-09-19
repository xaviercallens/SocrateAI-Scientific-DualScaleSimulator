"""Build audit/k3t2_rigidity_v2/comparison.json: one row for EVERY sealed-v2 target.

Run (no arguments; paths resolved relative to this file):
  /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python \
      audit/k3t2_rigidity_v2/comparison/third_computation_orientifold.py
  /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python \
      audit/k3t2_rigidity_v2/comparison/build_comparison.py
The sealed file is read from /mnt/disks/disk-socrateai-local-1/k3t2-sealed-v2/targets_sealed_v2.json
(outside the worktree, never committed).  If it is absent the script aborts.

Design
------
* A skeleton row is emitted for every entry of sealed['tracks'] (271), plus
  kummer_betti_euler_elsewhere (27), orientifold_statements (12) and
  out_of_scope (2) -- 312 rows.  target_id = "<group>#<index>:<lean name>".
* Every row defaults to NOT_COMPUTED with a reason; the MAP table below
  overrides rows where some blind track computed something comparable.
* Every comparable row carries a check: the blind value is LOADED from the
  blind track's JSON (never typed here), the Lean value is typed from the
  sealed statement (Lean values are allowed in this file), and the status is
  decided by the comparison actually evaluating -- a failing comparison
  becomes DISAGREE automatically.  NOT_COMPARABLE rows carry an explicit
  conversion function and record whether the values match under it.
* independent_routes counting rule (applied uniformly):
    - Tracks A and B both build Jacobi theta / phi_{0,1} series: ONE route
      ("theta") for any quantity both compute (c(D), chi from the genus).
    - Track D (GUDHI simplicial homology of a resolved Kummer) is a separate route.
    - Track C exact Gram-matrix lattice computations are a separate route; C's
      chained derivation (topology -> Hodge -> signature) consumes D's numbers,
      so it adds a route ONLY for the signature (3,19) (it is independent of
      C's Gram-matrix diagonalisation), never for b2 or chi.
    - Track E literature bookkeeping (quoted tier-L inputs + arithmetic) is one route.
    - The comparison-stage third computation (Kunneth on D's Betti numbers,
      torus fixed points) counts as one route and is labelled as such.
    - A value merely DERIVED in this script from a blind track's output
      counts under that track's route, not as a new one.
"""
import json
import os
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
V2 = os.path.dirname(HERE)
SEALED = "/mnt/disks/disk-socrateai-local-1/k3t2-sealed-v2/targets_sealed_v2.json"


def load(rel):
    with open(os.path.join(V2, rel)) as f:
        return json.load(f)


sealed = json.load(open(SEALED))
A = load("A-genus/results.json")
B1 = load("B-dyons/part1_euler_numbers_results.json")
B2 = load("B-dyons/part2_psi_m_results.json")
B3 = load("B-dyons/part3_polar_results.json")
BH = load("B-dyons/hurwitz_class_numbers_results.json")
BT = load("B-dyons/theta_forms_cache.json")
C = load("C-lattices/results.json")["items"]
CE8 = C["01_e8_result.json"]
CLAT = C["02_k3_mukai_gamma_result.json"]["results"]
CTD = C["03_symbolic_tduality_proofs_result.json"]
CTAD = C["04_tadpole_budget_result.json"]
CCH = C["06_chained_rigidity_result.json"]
D = load("D-tda/results.json")
D02 = load("D-tda/02_invariance_and_quotient_results.json")
D05 = load("D-tda/05_rigidity_scan_results.json")
D01 = load("D-tda/01_controls_results.json")
DK = load("D-tda/04_kunneth_results.json")
E = load("E-flux/results.json")
EB = load("E-flux/tadpole_bookkeeping_results.json")
T3 = load("comparison/third_computation_orientifold_results.json")

# ---------------------------------------------------------------- blind accessors
TW = A["twining_2A_3A_5A_7A"]


def A_n(n):
    return Fr(A["H_tau_appell_lerch"]["A_n_1_to_10"][str(n)])


def raw_untwined(n):
    """coefficient n of q^{1/8}H (H = 2q^{-1/8}(-1+sum A_n q^n)), n=0..10, from Track A."""
    if n == 0:
        return Fr(A["rigidity_scan_a_mu_coefficient"]["scan"][4]["polar_term_from_j0_slice"])
    return 2 * A_n(n)


def raw_twined(g, n):
    """coefficient n of q^{1/8}H_g, n=0..10, from Track A (raw table n<=8, 2*A_n^(g) beyond)."""
    t = TW[g]
    if n == 0:
        return Fr(t["H_g_polar_term_q_neg_1_8"])
    c = t["congruence_An_g_equiv_An_mod_ord_g_n_le_8"]
    if str(n) in c:
        return Fr(c[str(n)]["raw_H_g"])
    return 2 * Fr(t["A_n_g"][str(n)])


def cD(D_):
    return Fr(A["c_of_D"][str(D_)])


def cD_B(D_):
    return Fr(BT["cD_table"][str(D_)])


def H_B(D_):
    if D_ < 0:
        return Fr(0)
    if D_ == 0:
        return Fr(-1, 12)
    if D_ % 4 in (1, 2):
        return Fr(0)
    return Fr(BH["H_table"][str(D_)]["H"])


def forms_B(D_):
    return [tuple(f[:3]) for f in BH["H_table"][str(D_)]["reduced_forms"]]


def sig(name):
    s = CLAT[name]["signature_method_A_congruence_diagonalization"]
    s2 = CLAT[name]["signature_method_B_sturm_charpoly"]
    assert (s["p"], s["q"]) == (s2["p"], s2["q"])
    return (s["p"], s["q"])


def k3_betti():
    r = D["resolved_K3_betti_numbers"]
    b6 = [r["N=6"][f"b{k}"] for k in range(5)]
    b8 = [r["N=8"][f"b{k}"] for k in range(5)]
    assert b6 == b8
    return b6


def chi(b):
    return sum((-1) ** k * v for k, v in enumerate(b))


def all_dims(pred):
    return all(pred(CTD["by_dimension"][str(d)]) for d in (1, 2, 3, 4))


CHI_THETA = int(A["ZK3_tau0_const"])           # Track A elliptic genus
CHI_THETA_B = int(B1["chi_computed_c0_plus_2cm1"])
assert CHI_THETA == CHI_THETA_B
CHI_D = chi(k3_betti())                          # Track D GUDHI

# ---------------------------------------------------------------- mapping
# Each entry: dict(status_rule, lean, blind (callable -> value), conv (optional), src, routes, coverage, note)
# kind:
#   "eq"   : AGREE if blind() == lean else DISAGREE
#   "conv" : NOT_COMPARABLE; matches_under_conversion = (conv(blind()) == lean)
#   "true" : AGREE if blind() is True (a boolean blind check that the Lean claim's content holds) else DISAGREE
#   "fixed": status given explicitly (OUT_OF_SCOPE, NOT_COMPUTED with a reason)
M = {}


def eq(key, lean, blind, src, routes=1, coverage="full", note=""):
    M[key] = dict(kind="eq", lean=lean, blind=blind, src=src, routes=routes, coverage=coverage, note=note)


def tr(key, lean, blind, src, routes=1, coverage="full", note=""):
    M[key] = dict(kind="true", lean=lean, blind=blind, src=src, routes=routes, coverage=coverage, note=note)


def cv(key, lean, blind, conv, conv_text, src, routes=1, coverage="full", note=""):
    M[key] = dict(kind="conv", lean=lean, blind=blind, conv=conv, conv_text=conv_text, src=src,
                  routes=routes, coverage=coverage, note=note)


def fx(key, status, note, lean=None):
    M[key] = dict(kind="fixed", status=status, note=note, lean=lean)


def dis(key, lean, blind, src, routes, note):
    """expected-DISAGREE row: still decided by evaluation (AGREE if they happen to match)."""
    M[key] = dict(kind="eq", lean=lean, blind=blind, src=src, routes=routes, coverage="full", note=note)


# ======================= A-genus =======================
TABLES = {  # Lean CDH Table 20 columns, typed from the sealed statements' definitions (Twining.lean)
    "2A": [-2, -6, 14, -28, 42, -56, 86, -138, 188, -238],
    "3A": [-2, 0, -6, 10, 0, -18, 20, 0, -30, 42],
    "5A": [-2, 0, 2, 0, -6, 2, 0, 6, 0, -10],
    "7A": [-2, -1, 0, 0, 4, 0, -2, 2, -3, 0],
}
EOTA = [45, 231, 770, 2277, 5796, 13915, 30843, 65550, 132825]  # Lean eotA

eq("A-genus#7", [90, 462, 1540, 4554, 11592, 27830, 61686],
   lambda: [int(raw_untwined(n)) for n in range(1, 8)],
   "A-genus/results.json H_tau_appell_lerch.A_n_1_to_10 (x2) / twining raw_H_untwined", coverage="partial",
   note="RHS values (hComputed coefficients 1..7 = 2A_n) agree; the LHS irrep decomposition eotMult.chi1A is NOT computed by any track.")
eq("A-genus#8", TABLES["2A"][1:8], lambda: [int(raw_twined("2A", n)) for n in range(1, 8)],
   "A-genus/results.json twining_2A_3A_5A_7A.2A raw_H_g", coverage="partial",
   note="RHS (computed twined 2A series, levels 1..7) agrees; group-trace LHS (chi2A.eotMult) not computed. Blind uses the same tier-L CDH chi(g), F_g inputs as Lean.")
eq("A-genus#9", TABLES["3A"][1:8], lambda: [int(raw_twined("3A", n)) for n in range(1, 8)],
   "A-genus/results.json twining_2A_3A_5A_7A.3A raw_H_g", coverage="partial",
   note="RHS agrees; group-trace LHS not computed; same tier-L chi(g), F_g inputs as Lean.")
tr("A-genus#10", "3A traces (= 3A twined series) != 2A twined series, levels 1..7",
   lambda: [raw_twined("3A", n) for n in range(1, 8)] != [raw_twined("2A", n) for n in range(1, 8)],
   "A-genus/results.json twining 2A vs 3A (derived here)", note="derived in this script from Track A's two series; they differ at level 1 (0 vs -6).")
eq("A-genus#29", [90, 462, True, 27720],
   lambda: [int(raw_untwined(1)), int(raw_untwined(2)), raw_untwined(2) * 60 == 4 * raw_untwined(1) * 77, int(raw_untwined(2) * 60)],
   "A-genus/results.json identity_A2_60_eq_4_A1_77_untwined + A_n",
   note="Track A reports the lock in halved units (A_2*60=13860=4*A_1*77); Lean's coeff untwined = 2A_n; the lock is homogeneous so it holds in both; values given in Lean's (doubled) units.")
for i, g in zip((30, 31, 32, 33), ("2A", "3A", "5A", "7A")):
    tr(f"A-genus#{i}", f"ratio relation A2(g)*(4*A1(1A)) = A2(1A)*(4*A1(g)) FAILS at {g}",
       (lambda g=g: raw_twined(g, 2) * 4 * raw_untwined(1) != raw_untwined(2) * 4 * raw_twined(g, 1)),
       f"A-genus/results.json twining {g} (derived here)", note="derived in this script from Track A's computed coefficients.")
tr("A-genus#34", "literal lock A2*60 = 4*A1*77 fails at 2A",
   lambda: TW["2A"]["identity_A2_60_eq_4_A1_77_under_twining"]["holds"] is False,
   "A-genus/results.json twining_2A_3A_5A_7A.2A.identity_A2_60_eq_4_A1_77_under_twining (420 vs -924)")
eq("A-genus#59", [-2] + [2 * a for a in EOTA], lambda: [int(raw_untwined(n)) for n in range(0, 10)],
   "A-genus/results.json H_tau_appell_lerch.A_n_1_to_10 + polar term",
   note="Same normalisation H = 2q^{-1/8}(-1 + sum A_n q^n) on both sides. Lean computes h from (-2E2+48F2)/eta^3; blind extracts H as the Appell-Lerch finite part of the elliptic genus: different formulas, one blind route. Blind also has A_10=260568 (beyond Lean's range).")
eq("A-genus#60", [24, 0, 0, 0, 0, 0, 0, 0, 0, 0],
   lambda: [int(A["ZK3_tau_0"].get(str(n), "0")) for n in range(10)],
   "A-genus/results.json ZK3_tau_0 (only nonzero entries stored; q^1..q^10 are exactly 0)", routes=1,
   note="Track B's c(0)+2c(-1)=24 is the same theta route (q^0 only).")
eq("A-genus#61", [[2, 20, 2], [20, -128, 216, -128, 20]],
   lambda: [[int(cD(-1)), int(cD(0)), int(cD(-1))], [int(cD(0)), int(cD(3)), int(cD(4)), int(cD(3)), int(cD(0))]],
   "A-genus/results.json c_of_D (q^0: c(-1),c(0),c(-1); q^1: c(0),c(3),c(4),c(3),c(0))",
   note="Track B cD_table gives the same c(D) (same theta route).")
tr("A-genus#62", "(star) holds with chi=24 and h = computed H",
   lambda: A["rigidity_scan_a_mu_coefficient"]["scan"][4]["slices_agree"] is True and A["H_tau_appell_lerch"]["multiply_back_check_passes"] is True,
   "A-genus/results.json H_tau_appell_lerch (Z*eta^3 = 24*y^{1/2}*Psi - H*theta_1^2) + rigidity_scan_a N=24",
   note="Blind form of the same polar/finite decomposition (CDH 2.22-2.25) with polar multiplicity 24 and H whose coefficients agree with Lean's h (row A-genus#59); cleared differently (theta_1^2 vs (y-1)B1^2), so this is agreement of the identity, not of an intermediate series.")
tr("A-genus#63", "(star) fails at chi=23 and chi=25",
   lambda: [r["slices_agree"] for r in A["rigidity_scan_a_mu_coefficient"]["scan"] if r["N"] in (23, 25)] == [False, False],
   "A-genus/results.json rigidity_scan_a_mu_coefficient (N=20..28; only N=24 has agreeing slices)",
   note="Blind classification RIGID (structural two-slice z-independence selector, target-free). Covers a superset (20..28) of Lean's {23,25}.")
tr("A-genus#68", "24*H_g divisible by 24 (H_g integral) for 1A,2A,3A,5A,7AB",
   lambda: all(TW[g]["H_g_all_reported_coeffs_integral"] for g in TW) and all(raw_untwined(n).denominator == 1 for n in range(11)),
   "A-genus/results.json twining_*.H_g_all_reported_coeffs_integral + untwined raw coefficients")
for i, g in zip((70, 71, 72, 73), ("2A", "3A", "5A", "7A")):
    eq(f"A-genus#{i}", TABLES[g], (lambda g=g: [int(raw_twined(g, n)) for n in range(10)]),
       f"A-genus/results.json twining_2A_3A_5A_7A.{g} (raw_H_g n<=8; 2*A_n_g for n=9)",
       note="Same tier-L inputs (chi(g), F_g = c*Lambda_N from CDH) on both sides, so this checks the arithmetic of the twining formula and Track A's independent H(tau), not the literature inputs. Lean's column is CDH Table 20; blind never read Table 20.")
tr("A-genus#74", "twined24 9 8 1 0 (chi=8, F=0) != table2A*24",
   lambda: [8 * raw_untwined(n) for n in range(10)] != [24 * raw_twined("2A", n) for n in range(10)],
   "A-genus/results.json raw H and raw H_2A (derived here: 8*H vs 24*H_2A)", note="derived in this script; differ already at n=0 (-16 vs -48).")
tr("A-genus#75", "3A series != 2A table", lambda: [raw_twined("3A", n) for n in range(10)] != [raw_twined("2A", n) for n in range(10)],
   "A-genus/results.json (derived here)", note="derived in this script.")
for i, note in ((101, "A_1..A_5 = 45,231,770,2277,5796 agree with eotA (row A-genus#59); the M24 irrep-dimension side is not computed."),
                (102, "A_6 = 13915 agrees with eotA 5; the M24 decomposition 3520+10395 is not computed."),
                (103, "A_7 = 30843 agrees with eotA 6; the M24 decomposition is not computed."),
                (104, "no M24 dimension table in any blind track.")):
    fx(f"A-genus#{i}", "NOT_COMPUTED", note)

# ======================= B-dyons =======================
tr("B-dyons#0", "coefficients of Z_K3 depend only on 4n-l^2 (n<10, |l|<=12)",
   lambda: A["discriminant_dependence_holds"] is True and A["discriminant_mismatches"] == {} and BT["checks"]["B_index1_structure_holds"] is True,
   "A-genus/results.json discriminant_dependence_holds (q<=10, all l) ; B-dyons theta_forms_cache checks.B_index1_structure_holds",
   note="A and B are the same theta route (counted once). Blind found and fixed an edge-row truncation bug in its own code before this held.")
eq("B-dyons#1", [2, 20, -128, 216], lambda: [int(cD(-1)), int(cD(0)), int(cD(3)), int(cD(4))],
   "A-genus/results.json c_of_D (== B-dyons theta_forms_cache cD_table)")
fx("B-dyons#2", "NOT_COMPUTED", "Lean-implementation truncation guard (which c(D) the Lean product reads); blind analogue is B's QMAX-margin study, not the same statement.")
fx("B-dyons#3", "NOT_COMPUTED", "Lean-implementation truncation guard; no blind analogue.")
fx("B-dyons#4", "NOT_COMPUTED", "Lean-implementation exactness of integer divisions; blind uses Fractions throughout.")
eq("B-dyons#5", [1, 24, 324, 3200, 25650, 176256], lambda: [int(B1["Gk_z0_k0_8_DMVV_restriction"][str(k)]) for k in range(6)],
   "B-dyons/part1_euler_numbers_results.json Gk_z0_k0_8_DMVV_restriction (k=0..5)",
   note="Blind values are the DMVV product restricted to z=0, q^0 -- equal to p24(k) values; blind also has k=6,7,8.")
eq("B-dyons#6", [1, 24, 324, 3200, 25650, 176256], lambda: [int(B1["Gk_z0_k0_8_DMVV_restriction"][str(k)]) for k in range(6)],
   "B-dyons/part1_euler_numbers_results.json Gk_z0 (q^0 only)", coverage="partial",
   note="Only the q^0 coefficient is computed blind; Lean also proves the higher q-coefficients of G_k(tau,0) vanish (tau-independence) -- NOT computed blind. Blind itself labels the Goettsche comparison self-consistency (kappa = NORMALISATION).")
tr("B-dyons#7", "all coefficients of Z_K3 even", lambda: all(Fr(v).denominator == 1 and int(Fr(v)) % 2 == 0 for v in A["c_of_D"].values()),
   "A-genus/results.json c_of_D (derived here: every c(D), D<=40, is even)", note="derived in this script from the blind c(D) table (Z_K3 coefficients are c(4n-l^2)).")
tr("B-dyons#8", "DMZ (5.16) identities hold (k<=5, through q^1)",
   lambda: B2["m_1"]["holds_to_order_q^QCHK"] is True and B2["m_2"]["holds_to_order_q^QCHK"] is True and B2["G1_equals_2B_exactly"] is True,
   "B-dyons/part2_psi_m_results.json m_1, m_2 (to q^10), G1=2B", coverage="partial",
   note="Blind checks the m=1 and m=2 lines (k=2,3) to q^10 with a joint (c3,c4) negative control; m=0 is tautological; the m=3 and m=4 lines (k=4,5) were NOT computed blind (B could_not_do).")
tr("B-dyons#9", "DMZ (5.16) identities hold (k<=4, through q^2)",
   lambda: B2["m_1"]["holds_to_order_q^QCHK"] is True and B2["m_2"]["holds_to_order_q^QCHK"] is True,
   "B-dyons/part2_psi_m_results.json", coverage="partial", note="as B-dyons#8: m=1,2 lines only; m=3 line not computed blind.")
fx("B-dyons#10", "NOT_COMPUTED", "blind built A_{2,1} from the Fourier form only (with its s<=-1 bug fix); the equality of the direct strip expansion with (9.55) for m=1,2,3 was not computed.")
eq("B-dyons#11", 324, lambda: int(B3["tail_slope_structural_N_check"]["N_from_tail_slope_ratio"]),
   "B-dyons/part3_polar_results.json tail_slope_structural_N_check + N_M_joint_grid (m=1)", coverage="partial",
   note="Lean: polar coefficient p24(m+1) removes the pole for m=1,2,3. Blind: m=1 only, N=324=p24(2) from a target-free tail-slope and the unique joint (N,M) grid solution; m=2,3 not computed.")
tr("B-dyons#12", "p24(2)+-1 = 325, 323 fail (m=1)",
   lambda: B3["N_M_joint_grid"]["solutions_found (N,M) matching on drop_top2_keep_l20_no_D_filter"] == [[324, 648]] and B3["N_M_joint_grid"]["N_window"][0] <= 323 and B3["N_M_joint_grid"]["N_window"][1] >= 325,
   "B-dyons/part3_polar_results.json N_M_joint_grid (N in 310..339 unique 324)", coverage="partial", note="m=1 only.")
cv("B-dyons#13", [-1, 4, 6, 12, 12, 12, 16, 24], lambda: [H_B(d) for d in (0, 3, 4, 7, 8, 11, 12, 15)],
   lambda v: [int(12 * x) for x in v], "Lean h12(D) = 12*H(D); blind tabulates H(D) itself",
   "B-dyons/hurwitz_class_numbers_results.json H_table (reduced-form counting, Kronecker-Hurwitz validated)",
   note="blind H(0..15) = -1/12,1/3,1/2,1,1,1,4/3,2.")
tr("B-dyons#14", "G2 - 324*A*A21 vanishes to 2nd order at y=1 (division by A exact)",
   lambda: B3["N_M_joint_grid"]["solutions_found (N,M) matching on drop_top2_keep_l20_no_D_filter"] == [[324, 648]] and B3["invA0_times_A0_correct_near_low_l"] is True,
   "B-dyons/part3_polar_results.json", coverage="partial",
   note="implied, not checked directly: blind finds G2/A - 324*A21 equal to the holomorphic-in-z series 3E4A - 648*H through q^8, so G2 - 324*A*A21 = A*(finite) is divisible by A.")
cv("B-dyons#15", ["N=324 (p24 2)", "coefficient of 12H: 54"],
   lambda: [int(B3["tail_slope_structural_N_check"]["N_from_tail_slope_ratio"]), B3["N_M_joint_grid"]["solutions_found (N,M) matching on drop_top2_keep_l20_no_D_filter"][0][1]],
   lambda v: [f"N={v[0]} (p24 2)", f"coefficient of 12H: {Fr(v[1], 12)}"],
   "Lean hurwitzSer = 12*H, coefficient 54 on it; blind Hhat = H, coefficient M=648 = 54*12",
   "B-dyons/part3_polar_results.json (joint N,M grid, target-free)",
   note="blind: G2/A - 324*A21 = 3E4A - 648*Hhat through q^8 (Lean through q^3); joint (c3,c4) perturbation survives only at (0,0). The coefficient 3 of E4*A is part of the blind ANSATZ (DMZ-derived), not scanned, so it is not compared.")
tr("B-dyons#16", "without the class-number term the identity fails",
   lambda: B3["negative_control_M0"]["fails_as_expected"] is True, "B-dyons/part3_polar_results.json negative_control_M0")
fx("B-dyons#17", "NOT_COMPUTED", "DMZ Example 5 tables (E4*A and (B^2/A)^F coefficients) not computed by any track.")


def dmz_table(Ds, k):
    def f(D_):
        extra = Fr(0)
        if D_ >= 0 and D_ % (k * k) == 0:
            extra = k * H_B(D_ // (k * k))
        return -(H_B(D_) + extra)
    return [f(d) for d in Ds]


cv("B-dyons#18", [0, 0, 3, -6, -12, -12, -24, -24, -30], lambda: dmz_table([-4, -1, 0, 4, 7, 8, 12, 15, 16], 2),
   lambda v: [int(12 * x) for x in v], "Lean uses h12 = 12H; combination -(H(D) + 2H(D/4)) evaluated on blind H then x12",
   "B-dyons/hurwitz_class_numbers_results.json (combination formula taken from the Lean statement, evaluated here)",
   note="only the class-number TABLE is blind; the combination is Lean's. Checks Lean's table values, not DMZ (9.11) itself.")
cv("B-dyons#19", [0, 0, 0, 4, -4, -12, -12, -16, -24], lambda: dmz_table([-9, -4, -1, 0, 3, 8, 11, 12, 15], 3),
   lambda v: [int(12 * x) for x in v], "as B-dyons#18 with -(H(D) + 3H(D/9))",
   "B-dyons/hurwitz_class_numbers_results.json (combination from Lean statement)", note="as B-dyons#18.")
for i in range(20, 26):
    fx(f"B-dyons#{i}", "NOT_COMPUTED", "m=2,3 immortal / optimal-form identities (DMZ 9.11-9.13) not computed by any track (B covers m=1 only).")
for i in range(26, 33):
    fx(f"B-dyons#{i}", "NOT_COMPUTED", "self-dual T^2 charge-lattice root counting not computed by any track.")
for i in (33, 34, 35):
    fx(f"B-dyons#{i}", "NOT_COMPUTED", "M24 frame shapes / power maps not computed by any track.")
eq("B-dyons#36", [1, 24, 324, 3200, 25650], lambda: [int(B1["Gk_z0_k0_8_DMVV_restriction"][str(k)]) for k in range(5)],
   "B-dyons/part1_euler_numbers_results.json Gk_z0", coverage="partial",
   note="both sides of the Lean equation equal these values; blind reproduces the values (untwined Hilbert-scheme Euler numbers) but not Lean's frame-shape construction twinedHilb.")
for i in (37, 38):
    fx(f"B-dyons#{i}", "NOT_COMPUTED", "twined Goettsche character decomposition not computed by any track.")
for i in range(39, 51):
    fx(f"B-dyons#{i}", "NOT_COMPUTED", "twisted (M24-twined) DMVV / twined elliptic genera not computed by any track (blind twining is only of H(tau), Track A).")
fx("B-dyons#51", "NOT_COMPUTED", "Moore's form-count identity to D<=400 not computed (blind H table is D<=40 and uses a single counting route).")
fx("B-dyons#52", "NOT_COMPUTED", "Kummer-form criterion not computed by any track.")
eq("B-dyons#53", [[(1, 0, 3), (2, 2, 2)], True], lambda: [forms_B(12), (2, 0, 2) in forms_B(16)],
   "B-dyons/hurwitz_class_numbers_results.json H_table[12], H_table[16] reduced_forms", coverage="partial",
   note="reduced forms at D=12 and D=16 agree (same reduction convention); the first conjunct (Kummer-form filter over D<=16 = [12,16]) is NOT computed blind.")
eq("B-dyons#54", [[(1, 1, 1)], [(1, 0, 1)]], lambda: [forms_B(3), forms_B(4)],
   "B-dyons/hurwitz_class_numbers_results.json H_table[3], H_table[4]", coverage="partial",
   note="reduced forms agree; isKummerForm conjuncts not computed.")
fx("B-dyons#55", "NOT_COMPUTED", "D=64 is outside the blind Hurwitz table (D<=40); Kummer form not computed.")
for i in range(56, 62):
    fx(f"B-dyons#{i}", "NOT_COMPUTED", "Kummer code / Golay code / octad computations not done by any track.")
M["B-dyons#62"] = dict(kind="eq",
   lean=[8, 16, 24, 24, 22],
   blind=lambda: [sum(D02["N=6"]["quotient_T4_mod_Z2"]["betti_by_field"]["Z3"]["betti"][0::2]),
                  D02["N=6"]["num_fixed_vertices_computed"],
                  CHI_D,
                  D05["N=6"]["chi_U_computed"] + D05["N=6"]["num_fixed_points_computed"] * D05["N=6"]["chi_S2_computed"],
                  k3_betti()[2]],
   src="D-tda: 02 quotient Betti (1,0,6,0,1) even sum; 02 fixed vertices; results chi; 05 chi_U + 16*chi_S2; results b2",
   routes=1, coverage="partial",
   note="conjuncts 1-5 computed by Track D GUDHI (T^4/Z2 even Betti sum 8, 16 fixed points, chi 24, chi(U)+16*chi(S^2) = -8+32 = 24 i.e. Lean's (0-16)/2+16*2, b2 = 6+16 = 22); conjunct 6 (Kummer discriminant |(-4)^3| = 2^6) NOT computed and not included in the compared list; KummerE3.lean was untracked when the sealing scan started and committed in v3.20.0 before it ended (sealed _meta.drift_note).")

# ======================= C-lattices =======================
fx("C-lattices#0", "NOT_COMPUTED", "general lemma (any Gram with integer left inverse is unimodular); blind checks instances (see C-lattices#5).")
fx("C-lattices#1", "NOT_COMPUTED", "general lemma; blind checks evenness of specific lattices.")
tr("C-lattices#2", "cartanE8 has an integer inverse (cartanE8Inv * cartanE8 = 1)", lambda: CE8["determinant"]["value"] == "1",
   "C-lattices/01_e8_result.json determinant = 1", note="implied: det = 1 => integer inverse exists; Lean's explicit inverse matrix not compared entrywise.")
fx("C-lattices#3", "NOT_COMPUTED", "symmetry of the Lean matrix not checked by any track (blind Cartan matrix is symmetric by construction; labelling may differ).")
tr("C-lattices#4", "cartanE8 has even diagonal", lambda: CE8["is_even"] is True, "C-lattices/01_e8_result.json is_even")
tr("C-lattices#5", "cartanE8 unimodular", lambda: CE8["determinant"]["value"] == "1", "C-lattices/01_e8_result.json determinant")
fx("C-lattices#6", "NOT_COMPUTED", "symmetry of e8Neg not checked by any track.")
tr("C-lattices#7", "-E8 even diagonal", lambda: CE8["is_even"] is True and CLAT["K3_3U_2negE8"]["is_even"] is True,
   "C-lattices/01 is_even; 02 K3 = 3U+2(-E8) is_even", note="derived: diagonal of -E8 is -2.")
tr("C-lattices#8", "-E8 unimodular", lambda: Fr(CE8["determinant"]["value"]) * (-1) ** 8 in (1, -1),
   "C-lattices/01_e8_result.json determinant (det(-E8) = (-1)^8 det E8)", note="derived here.")
fx("C-lattices#9", "NOT_COMPUTED", "explicit LDL factor matrices not computed (blind has leading principal minors 2,4,6,5,4,3,2,1 in Bourbaki labelling; the Lean labelling was not checked, so D entries are not compared).")
eq("C-lattices#10", 1, lambda: int(CE8["determinant"]["value"]), "C-lattices/01_e8_result.json determinant")
tr("C-lattices#11", "E8 Cartan positive definite", lambda: CE8["is_positive_definite"] is True and all(int(m) > 0 for m in CE8["leading_principal_minors"]),
   "C-lattices/01_e8_result.json leading_principal_minors")
fx("C-lattices#12", "NOT_COMPUTED", "symmetry of U not separately checked.")


def mn_row(m, n):
    for r in CCH["lattice_mn_scan"]["full_scan_results"]:
        if (r["m"], r["n"]) == (m, n):
            return r


tr("C-lattices#13", "U even", lambda: mn_row(11, 0)["is_even"] is True,
   "C-lattices/06 lattice_mn_scan (m,n)=(11,0): 11U even", note="derived: 11U even => U even.")
fx("C-lattices#14", "NOT_COMPUTED", "U*U = 1 not computed.")
tr("C-lattices#15", "U unimodular", lambda: mn_row(11, 0)["is_unimodular"] is True and mn_row(11, 0)["determinant"] in ("-1", "1"),
   "C-lattices/06 lattice_mn_scan (11,0): det(11U) = -1", note="derived: det(U)^11 = -1 => det U = -1.")
fx("C-lattices#16", "NOT_COMPUTED", "explicit congruence of U to diag(2,-2) not computed.")
fx("C-lattices#17", "NOT_COMPUTED", "Lean-internal identification of two definitions (hyperbolicU vs Narain gram); no blind analogue.")
eq("C-lattices#18", (3, 19), lambda: sig("K3_3U_2negE8"), "C-lattices/02 K3 signature (congruence diag == Sturm) ; 06 chained derivation from D topology gives (3,19) too",
   routes=2, note="two blind routes: C Gram-matrix signature (2 exact methods) and C's chained topology->Hodge->signature derivation (consumes D's GUDHI Betti numbers + Calabi-Yau/Kahler/simply-connected inputs).")
eq("C-lattices#19", (4, 20), lambda: sig("Mukai_K3_plus_U"), "C-lattices/02 Mukai signature")
eq("C-lattices#20", (6, 22), lambda: sig("Gamma_6_22"), "C-lattices/02 Gamma^{6,22} signature")
eq("C-lattices#21", 22, lambda: CLAT["K3_3U_2negE8"]["rank"], "C-lattices/02 K3 rank ; D-tda b2 = 22", routes=2,
   note="C lattice rank and D GUDHI b2 are independent.")
eq("C-lattices#22", 28, lambda: CLAT["Gamma_6_22"]["rank"], "C-lattices/02 Gamma^{6,22} rank")
eq("C-lattices#23", [0, 0, 0], lambda: [(p - q) % 8 for p, q in (sig("K3_3U_2negE8"), sig("Mukai_K3_plus_U"), sig("Gamma_6_22"))],
   "C-lattices/02 signatures (derived here)")
eq("C-lattices#24", (3, 19), lambda: (int(CCH["derivation_chain"]["8_b2_plus_input_kahler"]), int(CCH["derivation_chain"]["11_b2_minus"])),
   "C-lattices/06 derivation_chain (b2+ from Kahler/Hodge, b2- = b2 - b2+)", note="b2+ uses the Kahler formula 2h^{2,0}+1 as a named input; b2 from D.")
for i in (25, 26, 27):
    fx(f"C-lattices#{i}", "NOT_COMPUTED", "general Mukai-pairing lemma; no blind analogue.")
for i in (28, 29):
    fx(f"C-lattices#{i}", "NOT_COMPUTED", "-U properties not computed.")
fx("C-lattices#30", "NOT_COMPUTED", "Mukai vector of the structure sheaf not computed.")
for i in (31, 32, 33):
    fx(f"C-lattices#{i}", "NOT_COMPUTED", "general reflection lemma; blind uses Weyl reflections to enumerate 240 roots but does not check the isometry statement.")
tr("C-lattices#34", "simple roots of -E8 have norm -2", lambda: all(CE8["cartan_matrix"][i][i] == 2 for i in range(8)),
   "C-lattices/01_e8_result.json cartan_matrix diagonal (derived: -E8 diagonal = -2)")
fx("C-lattices#35", "NOT_COMPUTED", "Weyl-reflection isometry of -E8 not checked (only used to generate roots).")

# ======================= D-tda (Lean T-duality/DFT + K3Topology) =======================
for i in (0, 1, 2):
    fx(f"D-tda#{i}", "NOT_COMPUTED", "coordinate projections: used implicitly in Track C's K_i, not checked as statements.")
tr("D-tda#3", "factorized T-duality K_k squares to 1", lambda: all_dims(lambda r: r["A_factorized_duality"]["all_squares_to_identity"]),
   "C-lattices/03_symbolic_tduality_proofs_result.json A_factorized_duality", coverage="partial",
   note="scope: blind symbolic check for d = 1..4 only; Lean is for all d.")
tr("D-tda#4", "K_k in O(d,d;Z)", lambda: all_dims(lambda r: r["A_factorized_duality"]["all_preserve_eta"]),
   "C-lattices/03 A_factorized_duality.all_preserve_eta", coverage="partial", note="d = 1..4 only.")
for i in (5, 6, 7, 8, 9):
    fx(f"D-tda#{i}", "NOT_COMPUTED", "not computed by any track (commutation / product to eta / charge norm).")
for i in (10, 11, 12, 13, 14):
    fx(f"D-tda#{i}", "NOT_COMPUTED", "mirror / tau-shift / theta-shift O(d,d) elements not computed by any track.")
tr("D-tda#15", "basisChange(A, A^-T) in O(d,d)", lambda: all_dims(lambda r: r["B_basis_change_preserves_eta"]["holds"]),
   "C-lattices/03 B_basis_change_preserves_eta (generic symbolic invertible A)", coverage="partial",
   note="d = 1..4 only; Lean hypothesis A^T B = 1 is the same as B = A^-T.")
for i in range(16, 35):
    fx(f"D-tda#{i}", "NOT_COMPUTED", "O(d,d) group / spectrum / B-shift / generalized-metric lemma not computed by any track.")
tr("D-tda#35", "eta H(G,0) eta = H(G^-1,0)", lambda: all_dims(lambda r: r["C_eta_H_eta_eq_Hinv"]["holds"]),
   "C-lattices/03 C_eta_H_eta_eq_Hinv", coverage="partial", note="d = 1..4 only (generic symbolic G).")
for i in (36, 37, 38, 39, 40, 41, 42, 43, 44, 45):
    fx(f"D-tda#{i}", "NOT_COMPUTED", "not computed by any track (mass form / section condition / definitional dualScale identities; v1's radius check was not redone in v2).")
tr("D-tda#46", "tr G + tr G^-1 >= 2d for G positive definite",
   lambda: all_dims(lambda r: r["D1_minimizer_identity_in_eigenvalues"]["difference_simplifies_to_zero"]) and CTD["D2_fully_symbolic_d2_spd_check"]["f_entries_minus_f_eigs_simplifies_to_zero"],
   "C-lattices/03 D1 (eigenvalue identity sum(x+1/x-2) = sum (x-1)^2/x, d=1..4) + D2 (entry-level tie-back, d=2)", coverage="partial",
   note="d=1..4 at eigenvalue level; entry level only d=2; relies on the standard spectral facts named in C's implicit_structural_inputs_for_D1_general_d.")
tr("D-tda#47", "dualScale(1) = 2d", lambda: all_dims(lambda r: r["D1_minimizer_identity_in_eigenvalues"]["difference_simplifies_to_zero"]),
   "C-lattices/03 D1 (equality case x_i = 1)", coverage="partial", note="derived: at x_i = 1 every (x_i-1)^2/x_i term vanishes; d=1..4.")
tr("D-tda#48", "dualScale G = 2d iff G = 1", lambda: all_dims(lambda r: r["D1_minimizer_identity_in_eigenvalues"]["difference_simplifies_to_zero"]),
   "C-lattices/03 D1 structural_nonnegativity (equality iff x_i = 1)", coverage="partial",
   note="d=1..4, eigenvalue level (x_i = 1 for all i and G symmetric => G = 1).")
fx("D-tda#49", "NOT_COMPUTED", "d=1 closed form (R+1/R)^2-2 not computed.")
tr("D-tda#50", "R + 1/R >= 2 for R > 0", lambda: CTD["by_dimension"]["1"]["D1_minimizer_identity_in_eigenvalues"]["difference_simplifies_to_zero"],
   "C-lattices/03 D1 at d=1 (x + 1/x - 2 = (x-1)^2/x)")
eq("D-tda#51", 24, lambda: chi(k3_betti()), "D-tda/results.json resolved K3 Betti (1,0,22,0,1) at N=6,8, Z/3,Z/5; A-genus Z_K3(tau,0) = 24",
   routes=2, note="Lean statement: K3BettiSum 1 0 22 0 1 = 24; blind Betti vector is exactly (1,0,22,0,1) (b1,b2 computed; b0,b3,b4 stated). chi also from the elliptic genus (theta route).")
eq("D-tda#52", -16, lambda: sig("K3_3U_2negE8")[0] - sig("K3_3U_2negE8")[1], "C-lattices/02 K3 signature (3,19) (derived: 3-19)", routes=2,
   note="C Gram signature and C chained derivation both give (3,19).")
eq("D-tda#53", 22, lambda: int(CCH["derivation_chain"]["7_h20_eq_h02_hodge_symmetry"]) + (k3_betti()[2] - 2) + int(CCH["derivation_chain"]["6_h02_from_hodge_decomposition"]),
   "D-tda b2 = 22 ; C-lattices/06 h^{2,0} = h^{0,2} = 1", coverage="partial",
   note="h^{1,1} = 20 is not independently computed (it is b2 - 2 here), so this row checks b2 = 22 and h^{2,0} = h^{0,2} = 1.")
eq("D-tda#54", [22, -16], lambda: [CLAT["K3_3U_2negE8"]["rank"], sig("K3_3U_2negE8")[0] - sig("K3_3U_2negE8")[1]],
   "C-lattices/02 K3 = 3U+2(-E8): rank 22, signature (3,19)", routes=1)
fx("D-tda#55", "NOT_COMPUTED", "holonomy group dimensions not computed.")
cv("D-tda#56", 2, lambda: int(CCH["derivation_chain"]["4_chi_O_noether"]), lambda v: v,
   "Lean: Dirac index -sigma/8; blind: holomorphic Euler characteristic chi(O) = (c1^2+c2)/12 (Noether). Equal for c1 = 0 surfaces (Todd = A-hat, tier L).",
   "C-lattices/06 derivation_chain.4_chi_O_noether", note="different quantity, same value under the tier-L identification.")
fx("D-tda#57", "NOT_COMPUTED", "parallel chiral spinor count not computed.")

# ======================= E-flux =======================
for i in (0, 1, 2):
    fx(f"E-flux#{i}", "NOT_COMPUTED", "general Kronecker-form lemma; no blind analogue.")
fx("E-flux#3", "NOT_COMPUTED", "general integrality lemma; Track E's instance (alpha_x^2 = 8*sum(+-n_i^2), TT basis A.6-A.7) is consistent with it but is not the general statement.")
eq("E-flux#4", 0, lambda: chi(D01["cubical"]["T2"]["betti"]), "D-tda/01_controls_results.json cubical T2 Betti (1,2,1)")
eq("E-flux#5", 24, lambda: CHI_D, "D-tda/results.json chi(resolved K3); A-genus Z_K3(tau,0)", routes=2,
   note="Lean computes it from the Hodge diamond; blind from GUDHI Betti numbers and from the elliptic genus.")
eq("E-flux#6", [0, 24], lambda: [T3["d3_tadpole"]["chi_K3xK3_mod_24"], int(Fr(T3["d3_tadpole"]["target_chi_K3xK3_over_24"]))],
   "comparison/third_computation_orientifold_results.json (Kunneth on D-tda's blind K3 Betti numbers)", routes=1,
   note="computed at the COMPARISON stage (not by a blind track) from Track D's Betti numbers: chi(K3xK3) = 576.")
eq("E-flux#7", 0, lambda: DK["N=6_Z3"]["chi_K3xT2_from_multiplicativity_chiK3_times_chiT2"], "D-tda/04_kunneth_results.json")
tr("E-flux#8", "flux + n = 24 => n <= 24 and (flux = 0 => n = 24)",
   lambda: CTAD["all_n_le_budget"] is True and CTAD["flux_zero_implies_n_eq_budget"] is True and CTAD["budget_used_in_enumeration"] == 24,
   "C-lattices/04_tadpole_budget_result.json (budget read from D-tda chi)", note="blind classifies this NORMALISATION (trivial arithmetic once the budget 24 is given).")

# ======================= elsewhere =======================
EW = "elsewhere"
eq(f"{EW}#0", 24, lambda: sum(k3_betti()[0::2]), "D-tda/results.json b0+b2+b4", note="Lean defs 1, 22, 1, 24.")
fx(f"{EW}#1", "NOT_COMPUTED", "Picard number of the singular Kummer not computed.")
fx(f"{EW}#2", "NOT_COMPUTED", "transcendental rank not computed.")
fx(f"{EW}#3", "NOT_COMPUTED", "Picard decomposition 16 + 4 not computed (Track D's 16 + 6 is the b2 decomposition, a different quantity).")
eq(f"{EW}#4", 24, lambda: sum(k3_betti()[0::2]), "D-tda/results.json", coverage="partial",
   note="only the Betti conjunct is computed; Picard/transcendental/cycle conjuncts not computed.")
eq(f"{EW}#5", 16, lambda: D02["N=6"]["num_fixed_vertices_computed"], "D-tda/02 fixed vertices",
   note="Lean statement is the type-cardinality tautology card(Fin 16) = 16; blind computes the 16 singular points.")
fx(f"{EW}#6", "NOT_COMPUTED", "intersection form not built by any track (D could_not_do).")
fx(f"{EW}#7", "NOT_COMPUTED", "exceptional-curve self-intersection not computed.")
fx(f"{EW}#8", "NOT_COMPUTED", "moduli-space geodesics not computed by any track.")
eq(f"{EW}#9", 24, lambda: CHI_D, "D-tda chi; A-genus chi", routes=2)
eq(f"{EW}#10", 20, lambda: k3_betti()[2] - 2 * int(CCH["derivation_chain"]["7_h20_eq_h02_hodge_symmetry"]),
   "D-tda b2 - 2*h^{2,0} (C-lattices/06)", coverage="partial", note="value h^{1,1} = 20 derived; the decomposition 16 + 4 is not computed.")
eq(f"{EW}#11", [22, 22], lambda: [sum(sig("K3_3U_2negE8")), k3_betti()[2]], "C-lattices/02 signature sum ; D-tda b2", routes=2)
eq(f"{EW}#12", 22, lambda: sum(sig("K3_3U_2negE8")), "C-lattices/02 (3,19)", routes=2, note="plus D b2 = 22.")
eq(f"{EW}#13", 16, lambda: D02["N=6"]["num_fixed_vertices_computed"], "D-tda/02 fixed vertices (also comparison third computation: 16 on T^4)", routes=1)
eq(f"{EW}#14", 0, lambda: chi(D01["cubical"]["T2"]["betti"]), "D-tda/01 T2")
eq(f"{EW}#15", 24, lambda: CHI_D, "D-tda; A-genus", routes=2)
fx(f"{EW}#16", "OUT_OF_SCOPE", "simulator swampland-safety flag (sealed track = UNCERTAIN); not K3xT2 mathematics.")
eq(f"{EW}#17", 24, lambda: CHI_D, "D-tda; A-genus", routes=2)
fx(f"{EW}#18", "NOT_COMPUTED", "exceptional intersection matrix not computed.")
eq(f"{EW}#19", 0, lambda: chi(DK["N=6_Z3"]["betti_K3xT2_via_kunneth"]), "D-tda/04 Kunneth (1,2,23,44,23,2,1)")
fx(f"{EW}#20", "NOT_COMPUTED", "Picard rank not computed.")
fx(f"{EW}#21", "NOT_COMPUTED", "Golay code not computed.")
eq(f"{EW}#22", 0, lambda: DK["N=6_Z3"]["chi_K3xT2_from_multiplicativity_chiK3_times_chiT2"], "D-tda/04")
fx(f"{EW}#23", "NOT_COMPUTED", "elliptic fibration singular fibres not computed.")
eq(f"{EW}#24", 24, lambda: CHI_D, "D-tda", routes=2)
eq(f"{EW}#25", 20, lambda: k3_betti()[2] - 2 * int(CCH["derivation_chain"]["7_h20_eq_h02_hodge_symmetry"]), "D-tda b2; C h^{2,0}")
eq(f"{EW}#26", 44, lambda: DK["N=6_Z3"]["betti_K3xT2_via_kunneth"][3], "D-tda/04 Kunneth b3(K3xT2)",
   note="Lean: 2*b2(K3) + b3(K3) = b3(K3xT2) = 44.")

# ======================= orientifold (sealed v3.20.0 values) =======================
OR = "orientifold"
tr(f"{OR}#0", "flux + n = 24 => ...", lambda: CTAD["all_n_le_budget"] and CTAD["flux_zero_implies_n_eq_budget"],
   "C-lattices/04", note="DUPLICATE of E-flux#8 (same file:line DualScaleStream2/Flux/Tadpole.lean:116).")
tr(f"{OR}#1", "net 7-brane charge 16*4 + 4*(-16) = 0", lambda: T3["seven_brane_charge"]["net"] == 0 and T3["seven_brane_charge"]["n_O7_computed"] == 4,
   "comparison/third_computation (4 O7 at T^2 fixed points, 16 D7 needed) ; E-flux TT quote lines 160-163", routes=1,
   note="counts (16 D7, 4 O7) agree with TT; per-object charges (+4, -16) are 4x the D7-unit charges (+1, -4); net 0 agrees in any unit.")
fx(f"{OR}#2", "OUT_OF_SCOPE", "master contract bundles the tadpole (see orientifold#1, AGREE) with simulator constants (sdc_mass_bound, bounce_action_numerator).")
eq(f"{OR}#3", 24, lambda: EB["part1"]["computed"]["orientifold_bookkeeping_total"], "E-flux/tadpole_bookkeeping_results.json part1 (TT eq 2.3)",
   routes=2, note="the D3 budget 24 agrees (E: 4*2+16*1 and the F-theory count; comparison third computation: chi(K3xK3)/24 = 24). Lean's flux term H3*F3 vs TT's (1/2)N_flux is a normalisation difference not compared.")
tr(f"{OR}#4", "total D7 + total O7 = 0", lambda: T3["seven_brane_charge"]["net"] == 0, "comparison/third_computation",
   note="net 0 agrees; the ingredient totals (+-64 in sealed TadpoleCancellation) DISAGREE, see orientifold#8, #9.")
tr(f"{OR}#5", "net_tadpole_charge = 0", lambda: T3["seven_brane_charge"]["net"] == 0, "comparison/third_computation")
cv(f"{OR}#6", 64, lambda: T3["seven_brane_charge"]["total_D7_charge_D7_units"], lambda v: 4 * v,
   "sealed KummerTadpole uses D7 charge 4 per brane; blind/third computation uses D7 units (charge 1). Lean/4 = 16 D7 units.",
   "comparison/third_computation_orientifold_results.json seven_brane_charge", note="counts agree (16 D7); matches under the x4 unit conversion.")
cv(f"{OR}#7", -64, lambda: T3["seven_brane_charge"]["total_O7_charge_D7_units"], lambda v: 4 * v,
   "sealed KummerTadpole O7 charge -16 per plane = 4 x (-4 D7 units)", "comparison/third_computation", note="4 O7 planes agree; matches under x4.")
dis(f"{OR}#8", -64, lambda: T3["seven_brane_charge"]["total_O7_charge_D7_units"], "comparison/third_computation_orientifold_results.json (T^2 fixed points = 4, Q(O7^-) = -4)",
    routes=1, note="DISAGREE. Sealed TadpoleCancellation.lean (v3.20.0) puts an O7 at each of the 16 fixed points of T^4/Z2 with Q = -4 D7 units (-64). The O7 planes of K3xT2/Z2 sit at the 4 fixed points of T^2/Z2 (TT lines 160-163, and the third computation's fixed-point count), giving -16. No unit conversion reconciles it: the plane COUNT differs (16 vs 4) and the sealed file's charge ratio Q(O7)/Q(D7) = -4/2 = -2 differs from the required -4 (a rescaling preserves ratios). The negative control in the third computation reproduces -64 exactly from the T^4 count. Post-seal LeanMaster (adc85e7, v3.21.0) changed this to numO7Planes = 4, total_O7_charge_is_minus_16.")
dis(f"{OR}#9", 64, lambda: T3["seven_brane_charge"]["total_D7_charge_D7_units"], "comparison/third_computation (16 D7 needed, unit charge)",
    routes=1, note="DISAGREE. Sealed: 32 D7 x charge 2 = 64. Required: 16 D7 x 1 = 16 (TT; third computation). Discriminator is the charge RATIO: sealed TadpoleCancellation has Q(O7)/Q(D7) = -4/2 = -2, whereas TT and the third computation have -4/1 = -4 (and sealed KummerTadpole has -16/4 = -4, which is why orientifold#6/#7 are NOT_COMPARABLE-matching). A unit rescaling preserves the ratio, so no rescaling reconciles -2 with -4. Post-seal LeanMaster (adc85e7) changed to 16 D7 x 1 = 16.")
tr(f"{OR}#10", "totalD7 + totalO7 = 0", lambda: T3["seven_brane_charge"]["net"] == 0, "comparison/third_computation",
   note="net 0 agrees; its two summands DISAGREE (orientifold#8, #9).")
dis(f"{OR}#11", 1, lambda: int(Fr(T3["d3_tadpole"]["target_chi_K3xK3_over_24"])),
    "E-flux/tadpole_bookkeeping (TT eq 2.3: 24) ; comparison/third_computation chi(K3xK3)/24 = 576/24 = 24 from D's Betti numbers",
    routes=2, note="DISAGREE. Sealed d3TadpoleTarget = chi(K3)/24 = 1. Blind: 24 (E: 4*2+16*1 = 24 and 24 (p,q) 7-branes x 1). Third computation: chi(K3xK3)/24 = 24 via Kunneth on Track D's blind K3 Betti numbers; chi(K3)/24 = 1 is reproduced only by using K3 instead of the F-theory 4-fold K3xK3. The sealed corpus is internally inconsistent here: sealed DualScaleStream2/Flux/Tadpole.lean:99 (k3k3_anomaly, chiK3K3/24 = 24; row E-flux#6 AGREE) and :116 (tadpole_budget, flux+n = 24; E-flux#8 AGREE) use 24, this file uses 1; the blind value sides with the former. Post-seal LeanMaster (adc85e7) changed to d3_tadpole_target_is_24.")

# ======================= out_of_scope =======================
fx("out_of_scope#0", "OUT_OF_SCOPE", "simulator TDA Mapper graph constant (V=187, E=557), not K3xT2 mathematics.")
fx("out_of_scope#1", "OUT_OF_SCOPE", "simulator TDA Mapper graph constant (V=187, E=557), not K3xT2 mathematics.")

# ---------------------------------------------------------------- default reasons
DEFAULT_REASON = {
    "A-genus": "M24 character-table / decomposition / HMN / umbral / other-class twining content: Track A took chi(g) as tier-L input for 2A,3A,5A,7A only and built no character table.",
    "B-dyons": "not computed by Track B.",
    "C-lattices": "not computed by Track C.",
    "D-tda": "not computed by any track.",
    "E-flux": "not computed by Track E.",
    "elsewhere": "not computed by any track.",
    "orientifold": "not computed by any track.",
}


def compact(s, n=400):
    s = " ".join(s.split())
    return s if len(s) <= n else s[: n - 3] + "..."


def jsonable(v):
    if isinstance(v, Fr):
        return str(v)
    if isinstance(v, tuple):
        return [jsonable(x) for x in v]
    if isinstance(v, list):
        return [jsonable(x) for x in v]
    return v


rows = []
groups = [(t, sealed["tracks"][t]) for t in sealed["tracks"]] + [
    ("elsewhere", sealed["kummer_betti_euler_elsewhere"]),
    ("orientifold", sealed["orientifold_statements"]),
    ("out_of_scope", sealed["out_of_scope"])]
used = set()
for gname, entries in groups:
    for i, e in enumerate(entries):
        key = f"{gname}#{i}"
        tid = f"{key}:{e['name']}"
        row = {"target_id": tid, "lean_file_line": f"{e['file']}:{e['line']}",
               "lean_statement": compact(e["statement"]), "committed_at_seal": e.get("committed")}
        m = M.get(key)
        if m is None:
            row.update(status="NOT_COMPUTED", independent_routes=0, lean_value=compact(e["statement"], 200),
                       blind_value=None, blind_source=None, coverage=None,
                       note=DEFAULT_REASON.get(gname, "not computed."))
        else:
            used.add(key)
            if m["kind"] == "fixed":
                row.update(status=m["status"], independent_routes=0, lean_value=compact(e["statement"], 200),
                           blind_value=None, blind_source=None, coverage=None, note=m["note"])
            else:
                bv = m["blind"]()
                lv = m["lean"]
                row["lean_value"] = json.dumps(jsonable(lv)) if not isinstance(lv, str) else lv
                row["blind_source"] = m["src"]
                row["independent_routes"] = m["routes"]
                row["coverage"] = m["coverage"]
                if m["kind"] == "eq":
                    ok = jsonable(bv) == jsonable(lv)
                    row["blind_value"] = json.dumps(jsonable(bv))
                    row["status"] = "AGREE" if ok else "DISAGREE"
                elif m["kind"] == "true":
                    row["blind_value"] = f"check evaluates to {bool(bv)}"
                    row["status"] = "AGREE" if bv is True else "DISAGREE"
                elif m["kind"] == "conv":
                    conv = m["conv"](bv)
                    ok = jsonable(conv) == jsonable(lv)
                    row["blind_value"] = json.dumps(jsonable(bv))
                    row["status"] = "NOT_COMPARABLE"
                    row["conversion"] = m["conv_text"]
                    row["blind_value_converted"] = json.dumps(jsonable(conv))
                    row["matches_under_conversion"] = ok
                row["note"] = m["note"]
                if row["status"] == "DISAGREE":
                    row["investigation"] = "third minimal computation: comparison/third_computation_orientifold.py" if gname == "orientifold" else "UNINVESTIGATED"
        rows.append(row)

unused = set(M) - used
assert not unused, f"mapping keys with no sealed target: {unused}"

counts = {}
for r in rows:
    counts[r["status"]] = counts.get(r["status"], 0) + 1
counts["TOTAL"] = len(rows)
counts["TOTAL_tracks_only"] = sum(len(v) for v in sealed["tracks"].values())
counts["AGREE_full_coverage"] = sum(1 for r in rows if r["status"] == "AGREE" and r.get("coverage") == "full")
counts["AGREE_partial_coverage"] = sum(1 for r in rows if r["status"] == "AGREE" and r.get("coverage") == "partial")
counts["NOT_COMPARABLE_matching_under_conversion"] = sum(1 for r in rows if r["status"] == "NOT_COMPARABLE" and r.get("matches_under_conversion"))
counts["TOTAL_distinct_targets"] = len(rows) - 1  # orientifold#0 == E-flux#8
track_groups = set(sealed["tracks"])
counts["DISAGREE_within_271_track_targets"] = sum(1 for r in rows if r["status"] == "DISAGREE" and r["target_id"].split("#")[0] in track_groups)
counts["DISAGREE_in_orientifold_group"] = sum(1 for r in rows if r["status"] == "DISAGREE" and r["target_id"].startswith("orientifold#"))
counts["rows_with_2plus_independent_routes"] = sum(1 for r in rows if (r.get("independent_routes") or 0) >= 2)
by_group = {}
for r in rows:
    g = r["target_id"].split("#")[0]
    by_group.setdefault(g, {})
    by_group[g][r["status"]] = by_group[g].get(r["status"], 0) + 1

BLIND_ONLY = [
    "A: A_10 = 260568 (Lean eotA stops at A_9) and c(D) for D = 19..40 (Lean c_first only D <= 4; c_depends_only_on_D covers q^9)",
    "A: twined H_g at q^10 for 2A,3A,5A,7A (Lean tables stop at q^9); congruence A_n^(g) = A_n mod ord(g) (raw coefficients, n <= 8)",
    "A: rigidity scan of the mu-term multiplicity N in 20..28 (Lean decomposition_pins_chi tests only 23, 25)",
    "B: DMZ (5.16) m=1, m=2 identities to q^10 (Lean dmz_516_q1/_q2 go to q^1/q^2) with a joint (c3,c4) negative control",
    "B: immortal m=1 identity through q^8 (Lean immortal_m1 through q^3) and the joint (N,M) = (324,648) uniqueness over a 900-pair grid",
    "B: Hurwitz class-number TABLE values H(D) for D = 16..40 with reduced forms, Kronecker-Hurwitz validated n=1..10 (Lean asserts individual values only for D <= 15; B-dyons#51 moore_vs_hurwitz asserts a relation over D < 401, not the values)",
    "B: Goettsche Euler numbers G_k(z=0) for k = 6,7,8 (1073720, 5930496, 30178575)",
    "C: 240 E8 roots by Weyl-reflection closure; Sturm-sequence second signature method; chained topology -> (m,n)=(3,2) selection with chi perturbation control",
    "D: GUDHI Betti numbers of resolved Kummer K3 (1,0,22,0,1) at N=6,8 over Z/3,Z/5; RP^3 link mod-2 torsion jump; translation-quotient control; N=4 premise failure",
    "D: Kunneth Betti vector of K3 x T^2 (1,2,23,44,23,2,1) as a full vector (Lean has only b3 = 44 and chi = 0 elsewhere)",
    "E: TT section-4.1 flux-vacuum pair counts (rank-4: 80/656/1616/3280 with tadpole, 80/1136/5600/18592 without; orbit counts 3/13/30/57; rank-6: 6960/339600) and the admissible set N_D3 in {16,8,0} from 8 | alpha_x^2",
    "E: TT worked-example checks (4.15)-(4.16), (4.31)-(4.32), symbolic N_flux = 2 alpha_x^2, and the flagged factor-2 slip at TT line 941",
    "E: BLPSSW 24 - 8 = 16 hidden-instanton count vs 16 orbifold singularities (tier L + arithmetic)",
    "E: unfixed moduli of TT's example (20 of 22 K3 Kahler directions + T2 Kahler modulus)",
    "comparison-stage: chi(K3 x K3) = 576 via Kunneth on D's blind Betti numbers; torus fixed points 4 (T^2), 16 (T^4)",
]

out = {
    "_meta": {
        "purpose": "Row-per-target comparison of the blind K3xT2 v2 results against the sealed LeanMaster target list.",
        "sealed_file": SEALED,
        "sealed_leanmaster_version": sealed["_meta"]["leanmaster_git_describe_tags"],
        "sealed_leanmaster_head": sealed["_meta"]["leanmaster_git_rev_parse_HEAD"],
        "headline": "All 271 track targets (plus 41 extra sealed statements) are rowed. 0 DISAGREE among the 271 track targets. 3 DISAGREE, all in the orientifold group (sealed v3.20.0 TadpoleCancellation.lean: O7 total -64, D7 total 64, D3 target 1). The sealed corpus itself is inconsistent on the D3 budget (Flux/Tadpole.lean uses 24); the blind value 24 sides with it. All three DISAGREEs were changed to the blind values by LeanMaster commit adc85e7 (v3.21.0), which the blind tracks never saw.",
        "leanmaster_drift_note": "LeanMaster HEAD at comparison time is 52ce329 (v3.21.0-1). Commit adc85e7 ('Tadpole correction: K3 x T2/Z2 has 4 O7-planes and 16 D7-branes; D3 budget is chi(K3 x K3)/24 = 24') changed exactly the orientifold rows found DISAGREE here (orientifold#8, #9, #11) in the direction of the blind values. lean_value in every row is the SEALED (v3.20.0) value; post-seal values appear only inside notes, labelled.",
        "commands": [
            "/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python audit/k3t2_rigidity_v2/comparison/third_computation_orientifold.py",
            "/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python audit/k3t2_rigidity_v2/comparison/build_comparison.py",
        ],
        "blind_inputs_read": ["A-genus/results.json", "B-dyons/part1_euler_numbers_results.json", "B-dyons/part2_psi_m_results.json",
                              "B-dyons/part3_polar_results.json", "B-dyons/hurwitz_class_numbers_results.json", "B-dyons/theta_forms_cache.json",
                              "C-lattices/results.json", "D-tda/results.json", "D-tda/01_controls_results.json",
                              "D-tda/02_invariance_and_quotient_results.json", "D-tda/04_kunneth_results.json",
                              "D-tda/05_rigidity_scan_results.json", "E-flux/results.json", "E-flux/tadpole_bookkeeping_results.json",
                              "comparison/third_computation_orientifold_results.json"],
        "status_definitions": {
            "AGREE": "a blind computation (or an arithmetic derivation in build_comparison.py from a blind output, labelled 'derived') gives the Lean value in the same convention. coverage='full' if the blind computation covers the whole numeric content of the statement over Lean's range; 'partial' if it covers a stated part (the note says which part is not covered).",
            "DISAGREE": "the evaluated comparison fails; every DISAGREE row names its investigation (third minimal computation).",
            "NOT_COMPUTED": "no blind track computed the statement's content; the note says what is missing.",
            "NOT_COMPARABLE": "the blind quantity is in a different convention/normalisation or is a different quantity with the same value; 'conversion' states the map and 'matches_under_conversion' whether the values agree after it.",
            "OUT_OF_SCOPE": "simulator constants or simulator contracts, not K3xT2 mathematics.",
        },
        "independent_routes_rule": "see module docstring of build_comparison.py: A and B theta series = one route; D GUDHI = one; C Gram-matrix = one (C's chained derivation adds a route only for the signature); E literature bookkeeping = one; comparison-stage third computation = one; derived-in-script values count under their source track. NOT_COMPUTED / OUT_OF_SCOPE rows have 0.",
        "excluded_false_positives_not_rowed": [x["name"] + " @" + x["file"] + ":" + str(x["line"]) for x in sealed["excluded_false_positives"]],
        "duplicates": ["orientifold#0 == E-flux#8 (DualScaleStream2/Flux/Tadpole.lean:116); both rowed, count it once when summing."],
        "tier_policy": "blind numbers are tier B (exact arithmetic with negative controls) or tier L (literature inputs), never 'proved'. AGREE means agreement of values, not a proof of the Lean statement's physical reading.",
    },
    "counts": counts,
    "counts_by_group": by_group,
    "blind_only": BLIND_ONLY,
    "rows": rows,
}
with open(os.path.join(V2, "comparison.json"), "w") as f:
    json.dump(out, f, indent=1, ensure_ascii=False)
print(json.dumps(counts, indent=1))
print(json.dumps(by_group, indent=1))
for r in rows:
    if r["status"] == "DISAGREE" or (r["status"] == "NOT_COMPARABLE" and not r.get("matches_under_conversion")):
        print("CHECK:", r["target_id"], r["status"], r.get("lean_value"), r.get("blind_value"), r.get("blind_value_converted"))
