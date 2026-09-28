#!/usr/bin/env python3
"""E-flux v3, script 05: TT section 5 (second branch, G_zbar = 0, N=2 supersymmetry).

Usage: <venv python> 05_second_branch.py NW KMAX
   NW   : entry bound of the fg-coordinate window for the dim V_flux = 2 enumeration (committed run: 1)
   KMAX : scan bound for the exact dim V_flux = 1 solve (j,k up to KMAX) (committed run: 24)
Committed run: cd audit/k3t2_rigidity_v3/E-flux/scripts && <venv python> 05_second_branch.py 1 24

Conditions (TT 5.1): (alpha_x - phi beta_x) tau = alpha_y - phi beta_y with Im phi, Im tau != 0; V_flux spanned by
timelike (negative norm) vectors, dim <= 2.  Convention: Im phi > 0 (phi = C0 + i e^-Phi) and Im tau > 0.

dim V_flux = 1 (EXACT solve, no grid): flux = 2*lambda0*(a,b,c,d), lambda0 primitive, lambda0^2 = -2j (j>=1).
   tau = (c - d phi)/(a - b phi) gives Im tau = (bc - ad) Im phi/|a - b phi|^2, so a solution with Im phi, Im tau > 0 exists iff k := bc - ad > 0,
   and N_flux = -beta_x.alpha_y + beta_y.alpha_x = 8 j k.  N_D3 = 24 - N_flux/2 >= 0 <=> j k <= 6 (the 24 is the declared input).
dim V_flux = 2 (window enumeration): alpha_y = p alpha_x + q beta_x, beta_y = r alpha_x + s beta_x, and (5.1) reduces to
   r phi^2 + (s-p) phi - q = 0, tau = p - r phi; Im phi > 0 and Im tau > 0 iff Delta=(s-p)^2+4rq < 0 and r < 0.
Counts are inside the window only and 'mod G+' means orbits of the finite group G+ (script 01), not physically distinct vacua.
"""
import sys, time, itertools, math
import numpy as np
import sympy as sp
from sympy import Rational as R
from common import *
from lattice_u3 import H33_np, group_G, orient_plus
from inputs_decl import TADPOLE_TOTAL, HALF

SCRIPT = "05_second_branch.py"
H = sp.Matrix(H33_np.tolist())
BASE = 16; OFF = 8

def enc(V):
    return ((V + OFF) * (BASE ** np.arange(6))).sum(axis=1)

def nflux_lam(b, v, w, a, Hn):
    """N_flux = -beta_x.alpha_y + beta_y.alpha_x with alpha = 2 lambda"""
    return 4 * (-(b @ Hn @ v) + (w @ Hn @ a))

def dim1_exact(KMAX):
    rows = []
    for j in range(1, KMAX + 1):
        for k in range(1, KMAX + 1):
            Nf = 8 * j * k
            ND3 = TADPOLE_TOTAL - R(Nf, HALF)
            if ND3 >= 0:
                # realisation: lambda0 = f1 - j g1  (norm -2j), (a,b,c,d) = (0,1,k,0)
                l0 = [1, -j, 0, 0, 0, 0]
                n0 = sum(l0[i] * H[i, jj] * l0[jj] for i in range(6) for jj in range(6))
                ax = [0] * 6; bx = [2 * x for x in l0]; ay = [2 * k * x for x in l0]; by = [0] * 6
                Nf_real = -sum(bx[i] * H[i, jj] * ay[jj] for i in range(6) for jj in range(6)) + sum(by[i] * H[i, jj] * ax[jj] for i in range(6) for jj in range(6))
                # complex moduli from (5.1): phi tau = -k ; take phi = i t (t>0): tau = -k/(i t) = i k / t
                rows.append({"j": j, "k": k, "lambda0_norm": int(n0), "N_flux_scan": Nf, "N_flux_realised_from_lattice_Gram": int(Nf_real), "N_D3": str(ND3),
                             "phi_tau_relation": f"phi*tau = -{k}", "consistent": bool(int(Nf_real) == Nf and n0 == -2 * j)})
    return rows

def dim2_window(NW, Gplus):
    r = np.arange(-NW, NW + 1)
    g = np.meshgrid(*([r] * 6), indexing="ij")
    V = np.stack([x.ravel() for x in g], axis=1).astype(np.int64)
    Hn = H33_np
    n = len(V)
    keys = enc(V); order = np.argsort(keys); ks = keys[order]
    def index_of(W):
        k = enc(W); p = np.clip(np.searchsorted(ks, k), 0, n - 1)
        return np.where(ks[p] == k, order[p], -1)
    q = np.einsum("ij,jk,ik->i", V, Hn, V)
    imgs = np.array([index_of(V @ M.T) for M in Gplus]); assert (imgs >= 0).all()
    lab = imgs.min(axis=0)
    NEG = np.nonzero(q < 0)[0]
    reps = [i for i in NEG if lab[i] == i]
    HV = V @ Hn
    total = 0; orbit_total = 0; nf_hist = {}; examples = []
    all_reps = []
    for a in reps:
        stab = [k for k in range(len(Gplus)) if imgs[k][a] == a]
        osz = len(Gplus) // len(stab)
        Tl = []
        for b in NEG:
            D2 = q[a] * q[b] - (HV[a] @ V[b]) ** 2
            if D2 <= 0:      # need independent with negative-definite span:  qa<0 and det>0
                continue
            # coefficients of every window vector in basis (a,b): pick pivot 2x2 minor
            A = V[a]; B = V[b]
            piv = None
            for i in range(6):
                for j in range(i + 1, 6):
                    if A[i] * B[j] - A[j] * B[i] != 0:
                        piv = (i, j); break
                if piv: break
            i, j = piv
            Dm = A[i] * B[j] - A[j] * B[i]
            P = V[:, i] * B[j] - V[:, j] * B[i]          # p*Dm
            Q = A[i] * V[:, j] - A[j] * V[:, i]           # q*Dm
            rec = P[:, None] * A[None, :] + Q[:, None] * B[None, :]
            inspan = np.all(rec == V * Dm, axis=1)
            Mem = np.nonzero(inspan)[0]
            Pm = P[Mem]; Qm = Q[Mem]
            # pairs (v,w) in Mem^2: alpha_y=v (p=Pm/Dm,q=Qm/Dm), beta_y=w (r=Pm/Dm, s=Qm/Dm)
            Pv = Pm[:, None]; Qv = Qm[:, None]; Pw = Pm[None, :]; Qw = Qm[None, :]
            r_neg = (Pw * Dm) < 0
            disc = ((Qw - Pv) ** 2 + 4 * Pw * Qv) < 0             # Delta*Dm^2 < 0
            ok = r_neg & disc
            if not ok.any():
                continue
            vi, wi = np.nonzero(ok)
            v = Mem[vi]; w = Mem[wi]
            Nf = 4 * (-(HV[b] @ V[v].T) + (HV[a] @ V[w].T))
            keep = (Nf > 0) & (Nf <= 2 * TADPOLE_TOTAL) & (Nf % 2 == 0)
            for x, y, nfv in zip(v[keep], w[keep], Nf[keep]):
                Tl.append((int(b), int(x), int(y), int(nfv)))
        if not Tl:
            continue
        T = np.array(Tl, dtype=np.int64)
        total += osz * len(T)
        for nfv in np.unique(T[:, 3]):
            nf_hist[int(nfv)] = nf_hist.get(int(nfv), 0) + osz * int((T[:, 3] == nfv).sum())
        can = None
        for k in stab:
            im = imgs[k]
            code = (im[T[:, 0]] * n + im[T[:, 1]]) * n + im[T[:, 2]]
            can = code if can is None else np.minimum(can, code)
        uc = np.unique(can)
        orbit_total += len(uc)
        for c in uc[:2]:
            w_ = c % n; v_ = (c // n) % n; b_ = c // (n * n)
            examples.append({"lambda_ax": V[a].tolist(), "lambda_bx": V[b_].tolist(), "lambda_ay": V[v_].tolist(), "lambda_by": V[w_].tolist()})
        all_reps += [(int(a), int(c // (n * n)), int((c // n) % n), int(c % n)) for c in uc]
    # exact check of (5.1) solvability for a sample of representatives
    ex_check = []
    rs = np.random.RandomState(3)
    sel = all_reps if len(all_reps) <= 60 else [all_reps[i] for i in rs.choice(len(all_reps), 60, replace=False)]
    n_ok = 0
    for a, b, v, w in sel:
        A, B, Vy, Wy = [2 * V[i] for i in (a, b, v, w)]
        # solve alpha_y = p alpha_x + q beta_x etc exactly with sympy
        M = sp.Matrix([list(map(int, A)), list(map(int, B))]).T
        sol_v = M.gauss_jordan_solve(sp.Matrix(list(map(int, Vy))))[0]; sol_w = M.gauss_jordan_solve(sp.Matrix(list(map(int, Wy))))[0]
        p, qq = sol_v; rr, ss = sol_w
        phi_sym = sp.symbols("phi")
        roots = sp.solve(rr * phi_sym ** 2 + (ss - p) * phi_sym - qq, phi_sym)
        good = False
        for ph in roots:
            if sp.im(ph) > 0:
                tau = p - rr * ph
                lhs = [(A[i] - ph * B[i]) * tau for i in range(6)]; rhs = [Vy[i] - ph * Wy[i] for i in range(6)]
                if all(sp.simplify(lhs[i] - rhs[i]) == 0 for i in range(6)) and sp.im(tau) > 0:
                    good = True
        n_ok += int(good)
    return {"window": f"fg entries |.|<={NW}", "n_window_vectors": int(n), "negative_norm_lambda_orbit_reps_for_alpha_x": len(reps),
            "tuples_(dim V_flux=2, second branch, tadpole 0<N_flux<=48)": int(total), "orbits_mod_Gplus": int(orbit_total),
            "tuples_by_N_flux": {str(k): v for k, v in sorted(nf_hist.items())}, "N_D3_by_N_flux": {str(k): str(TADPOLE_TOTAL - R(k, HALF)) for k in sorted(nf_hist)},
            "exact_check_(5.1)_solution_with_Im_phi>0,Im_tau>0": {"checked": len(sel), "ok": n_ok}, "example_orbit_reps": examples[:4]}

def main():
    NW = int(sys.argv[1]); KMAX = int(sys.argv[2])
    CMD = rel_command(SCRIPT, " ".join(sys.argv[1:]))
    t0 = time.time()
    Gplus = [M for M in group_G() if orient_plus(M)]
    d1 = dim1_exact(KMAX)
    d2 = dim2_window(NW, Gplus)
    # negative control: without tadpole, every (j,k) is admissible (unbounded)
    n_without = KMAX * KMAX
    quotes = {"5.4_5.5": quote("which is 36 dimensional"), "5.4_5.5b": quote("which is 52 dimensional"), "5.3_38": quote("which is 38 dimensional and the moduli space of"), "5.3_55": quote("identifications), which is 55 dimensional"),
              "5.3_phi_tau": quote("φτ = −1."), "second_branch_all_orbifold": quote("branch correspond to orbifold singularities. From §2.3"),
              "N2": quote("essary and sufficient conditions for N = 2 supersymmetry"), "T2_volume": quote("T 2 is not constrained at all")}
    res = [{"id": "E05a", "quantity": "second branch, dim V_flux = 1: exact solve of the tadpole for (j,k) with N_flux = 8jk",
            "computed": {"allowed_(j,k)_from_scan_up_to_KMAX": [(r["j"], r["k"]) for r in d1], "n_allowed": len(d1), "N_flux_values": sorted({r["N_flux_scan"] for r in d1}),
                         "N_D3_values": sorted({r["N_D3"] for r in d1}, key=int), "rows": d1, "all_realisations_consistent": all(r["consistent"] for r in d1),
                         "KMAX": KMAX, "scan_pairs_without_tadpole_(negative_control)": n_without,
                         "note": "The scan grid is (j,k) in 1..KMAX x 1..KMAX stepping through every integer; the allowed set is read from the scan. "
                                 "Each (j,k) is only a CLASS label: flux = 2 lambda0 (a,b,c,d) with bc-ad = k has infinitely many integer matrices; "
                                 "identification of matrices by SL(2,Z) x SL(2,Z) is NOT asserted here."},
            "shared_inputs": ["TT_H33_gram_A3", "flux_quantisation_even_coefficients", "branch2_conditions_5.1_5.3", "N_flux_formula_2.8", "tadpole_total_24_eq2.3", "tadpole_half_factor_eq2.3"], "script": SCRIPT, "command": CMD},
           {"id": "E05b", "quantity": "second branch, dim V_flux = 2: enumeration in a window, mod G+",
            "computed": d2, "shared_inputs": ["TT_H33_gram_A3", "flux_quantisation_even_coefficients", "branch2_conditions_5.1_5.3", "N_flux_formula_2.8", "tadpole_total_24_eq2.3", "tadpole_half_factor_eq2.3", "group_G_definition", "window_bounds"], "script": SCRIPT, "command": CMD},
           {"id": "E05c", "quantity": "unfixed moduli in the second branch (TT quotations)", "computed": {"quotes": quotes,
             "summary": "second branch: dim V_flux<=2; only ONE complex relation between phi and tau (dim 1: phi*tau = -k, dim 2: tau = p - r phi with r phi^2+(s-p)phi-q=0 fixes both); T2 Kahler modulus unconstrained; "
                        "K3 complex structure only partially fixed (Grassmannian O(3,17) or O(3,18) type, 36 or 38 dimensional); N=2 supersymmetry; all such solutions sit at K3 orbifold points (TT)."},
            "shared_inputs": ["branch2_conditions_5.1_5.3"], "script": SCRIPT, "command": CMD}]
    rig = [{"parameter_inserted": "the set of (j,k) with N_flux = 8jk <= 48, i.e. jk <= 6 (second branch, dim 1)", "selecting_condition": "N_D3 = 24 - N_flux/2 >= 0 with N_flux = 8jk (from the lattice Gram) and k>0 for Im phi, Im tau>0",
            "condition_uses_true_value": True, "solution_set": "the (j,k) with jk<=6 (14 pairs; read from E05a)", "classification": "CONDITIONAL_ON_INPUT", "input_it_depends_on": "tadpole_total_24_eq2.3",
            "negative_control": "without the tadpole every one of the KMAX^2 scanned (j,k) is admissible", "control_perturbs_same_parameter": True}]
    write_json("05_second_branch.json", {"script": SCRIPT, "results": res, "rigidity": rig, "literature_checks": [],
        "could_not_do": ["second branch dim V_flux = 2 enumerated only in the stated window; dim V_flux = 1 solved exactly as classes (j,k) but not counted modulo a duality group"], "wall_seconds": round(time.time() - t0, 1)})
    print(json.dumps(res[1]["computed"], default=str)[:2500])
    print("dim1 allowed", [(r["j"], r["k"]) for r in d1])

if __name__ == "__main__":
    main()
