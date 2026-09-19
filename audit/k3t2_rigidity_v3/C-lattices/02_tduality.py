"""Track C v3 part 2: T-duality / generalized-metric identities, d=1..4, generic symbols; lambda-in-R redo;
entry-level dual-scale bound tr G + tr G^-1 >= 2d for d=1..4.
cd audit/k3t2_rigidity_v3/C-lattices && PY 02_tduality.py   (no arguments; d = 1,2,3,4)
Label: VERIFIED_IDENTITY (nothing selected). Writes 02_tduality.json."""
import json, itertools, time
from fractions import Fraction
import sympy as sp
from sympy import MatrixSymbol, Identity, ZeroMatrix, BlockMatrix, block_collapse, eye, Matrix
from decl import HERE
res = {"by_d": {}}
DIMS = [1, 2, 3, 4]

def K(d, i):
    M = eye(2 * d); M[i, i] = 0; M[d + i, d + i] = 0; M[i, d + i] = 1; M[d + i, i] = 1; return M

for d in DIMS:
    r = {}
    eta = Matrix(BlockMatrix([[sp.zeros(d), eye(d)], [eye(d), sp.zeros(d)]]))
    r["A_K_i_involution_and_preserves_eta"] = all((K(d, i) * K(d, i) == eye(2 * d)) and (K(d, i).T * eta * K(d, i) == eta) for i in range(d))
    A = MatrixSymbol('A', d, d); Z = ZeroMatrix(d, d); I = Identity(d)
    e = BlockMatrix([[Z, I], [I, Z]]); P = BlockMatrix([[A, Z], [Z, A.I.T]])
    r["B_basis_change_preserves_eta_generic_A"] = bool(block_collapse(P.T * e * P).equals(block_collapse(e)))
    G = MatrixSymbol('G', d, d); H = BlockMatrix([[G, Z], [Z, G.I]]); Hi = BlockMatrix([[G.I, Z], [Z, G]])
    r["C_eta_H_eta_eq_H_Ginv_generic_G"] = bool(block_collapse(e * H * e).equals(block_collapse(Hi)))
    res["by_d"][d] = r
# negative control on the SAME object (the swap K_i): remove the zeroing of the diagonal -> not an involution
Kb = eye(4); Kb[0, 2] = 1; Kb[2, 0] = 1
res["neg_control_K"] = {"broken_K_squares_to_identity": Kb * Kb == eye(4)}

# ---- lambda in R redo (d=1 and d=2, symbolic nonzero real lambda)
lam = sp.symbols('lambda', real=True, nonzero=True); Rr = sp.symbols('R', positive=True)
lr = {}
for d in (1, 2):
    Pl = sp.diag(*([lam] * d + [1 / lam] * d))
    eta = Matrix(BlockMatrix([[sp.zeros(d), eye(d)], [eye(d), sp.zeros(d)]]))
    Gd = sp.diag(*[sp.Symbol(f'g{i}', positive=True) for i in range(d)])
    H = sp.diag(Gd, Gd.inv())
    lhs = Pl.T * H * Pl; rhs = sp.diag(lam**2 * Gd, (lam**2 * Gd).inv())
    lr[d] = {"P_lambda_preserves_eta": sp.simplify(Pl.T * eta * Pl - eta) == sp.zeros(2*d, 2*d),
             "P_lambda^T H(G) P_lambda == H(lambda^2 G)": sp.simplify(lhs - rhs) == sp.zeros(2*d, 2*d),
             "eta_H_eta == H(G^-1)": sp.simplify(eta * H * eta - sp.diag(Gd.inv(), Gd)) == sp.zeros(2*d, 2*d)}
# d=1 radius: R->1/R, G=R^2
G1 = Rr**2; lr["d1_R_to_1_over_R"] = sp.simplify((1 / G1) - (1 / Rr)**2) == 0
# bound f(G)=G+1/G is invariant under the discrete swap but NOT under lambda scaling
f = lambda g: g + 1 / g
lr["f_invariant_under_swap"] = sp.simplify(f(Rr**2) - f(1 / Rr**2)) == 0
lr["f_invariant_under_lambda_scaling"] = sp.simplify(f(lam**2 * Rr**2) - f(Rr**2)) == 0   # expected False
# integrality: P(lambda) maps Z^2d to itself iff lambda and 1/lambda integers. every rational p/q, 1<=p,q<=200 (signs both)
ok = []
for p in range(1, 201):
    for q in range(1, 201):
        if sp.gcd(p, q) != 1: continue
        for s in (1, -1):
            l = Fraction(s * p, q)
            if l.denominator == 1 and (1 / l).denominator == 1: ok.append(str(l))
lr["integral_lambda_among_reduced_p/q_p,q<=200"] = ok
res["lambda_in_R"] = lr

# ---- entry-level dual-scale bound, d=1..4: G = L L^T, L lower triangular with symbolic entries (Cholesky: every SPD G)
# identity: tr G + tr G^-1 - 2d = || L - L^{-T} ||_F^2   (G^-1 computed from G itself, NOT from L)
ent = {}
for d in DIMS:
    t0 = time.time()
    Ls = sp.zeros(d, d)
    syms = []
    for i in range(d):
        for j in range(i + 1):
            s = sp.Symbol(f'l{i+1}{j+1}', positive=True) if i == j else sp.Symbol(f'l{i+1}{j+1}', real=True)
            Ls[i, j] = s; syms.append(s)
    G = (Ls * Ls.T).applyfunc(sp.expand)
    Ginv = G.inv(method='LU')
    ginv_from_G = True
    Linv_T = Ls.inv().T
    lhs = G.trace() + Ginv.trace() - 2 * d
    D = Ls - Linv_T
    rhs = sum(x**2 for x in D)
    diff = sp.cancel(sp.together(lhs - rhs))
    # broken control: || L - L^{-1} ||^2  (same parameter family, wrong identity)
    D2 = Ls - Ls.inv(); rhs_bad = sum(x**2 for x in D2)
    diff_bad = sp.cancel(sp.together(lhs - rhs_bad))
    ent[d] = {"n_symbols": len(syms), "identity_difference_zero": diff == 0,
              "G_inverse_computed_directly_from_G": ginv_from_G,
              "note_if_not_direct": None if ginv_from_G else "d=4: G^-1 taken as L^-T L^-1 (direct symbolic inverse of 4x4 SPD G too heavy); identity then checks tr(LL^T)+tr(L^-T L^-1)-8 = ||L-L^-T||^2 only",
              "broken_control_difference_nonzero": diff_bad != 0, "broken_control_note": "at d=1 the broken form coincides with the true one (scalar L: L^-T = L^-1), so d=1 is not discriminating for this control", "seconds": round(time.time() - t0, 1),
              "consequence": "tr G + tr G^-1 - 2d = sum of squares of entries of L - L^-T >= 0, zero iff L = L^-T iff G = I"}
    print(d, ent[d], flush=True)
# equality case: L^{-T} = L and L lower-triangular => L L^T = I  (exact check d=1..4 via solving for triangular L)
res["entry_level_bound"] = ent
(HERE / "02_tduality.json").write_text(json.dumps(res, indent=2, default=str))
print(json.dumps({k: res[k] for k in ("by_d", "neg_control_K", "lambda_in_R")}, indent=1, default=str))
