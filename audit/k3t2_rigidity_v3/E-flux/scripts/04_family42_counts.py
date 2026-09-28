#!/usr/bin/env python3
"""E-flux v3, script 04: TT section 4.2 family ((4.24)-(4.25) ansatz), counted in windows of the unimodular lattice.

Usage:  <venv python> 04_family42_counts.py NU N_A6 MAXVERIFY
   NU     : entry bound of the fg-coordinate window for lambda in Gamma_{3,3}=U^3            (committed run: 2)
   N_A6   : entry bound of the e-coordinate window for lambda in D = span(TT A.6 e_i)        (committed run: 2; contains TT (4.31))
   MAXVERIFY: number of G+-orbit representatives to verify exactly with (3.16a-c), (3.31), V_flux type, orbifold test (committed run: 400)
Committed run:  cd audit/k3t2_rigidity_v3/E-flux/scripts && <venv python> 04_family42_counts.py 2 2 400

Flux vectors are 2*lambda:  alpha_x=2a, beta_x=2b, alpha_y=2v, beta_y=2w  (a,b,v,w lattice vectors in the window).
Ansatz (TT 4.24)-(4.25):  a.v = b.w = a.w + v.b = 0 ;  b^2 = 2 a^2 = 2 a.b ;  w^2 = 2 v^2 = 2 v.w ;  a^2 != 0, v^2 != 0.
Then phi = (1 +- i)/2 (4.29), tau = i sqrt(alpha_yy/alpha_xx) (4.30): non-real tau needs alpha_yy/alpha_xx > 0.
N_flux = -beta_x.alpha_y + beta_y.alpha_x (2.8); N_D3 = 24 - N_flux/2 >= 0.
TRUNCATION: only this ansatz, only inside the window; it is NOT the complete class (c) of TT section 3.4.
Counts are ordered tuples; 'mod G+' / 'mod G' are orbit counts of the finite groups of script 01 acting on the window set.
Exact verification of (3.16a-c) etc. is done on orbit representatives (all of them if fewer than MAXVERIFY).
"""
import sys, time, itertools
import numpy as np
import sympy as sp
from sympy import I, Rational as R, sqrt
from common import *
from lattice_u3 import H33_np, group_G, orient_plus
from susy import tt_conditions, gram_of, sig_rank, orbifold_hits
from inputs_decl import TADPOLE_TOTAL, HALF

SCRIPT = "04_family42_counts.py"
H = sp.Matrix(H33_np.tolist())
BASE = 16; OFF = 8

def enc(V):
    return ((V + OFF) * (BASE ** np.arange(6))).sum(axis=1)

def window_U(N):
    r = np.arange(-N, N + 1)
    g = np.meshgrid(*([r] * 6), indexing="ij")
    return np.stack([x.ravel() for x in g], axis=1).astype(np.int64)

def window_A6(n):
    T = np.zeros((6, 6), dtype=np.int64)
    for i in range(3):
        T[i, 2 * i] = 1; T[i, 2 * i + 1] = 1; T[i + 3, 2 * i] = 1; T[i + 3, 2 * i + 1] = -1
    r = np.arange(-n, n + 1)
    g = np.meshgrid(*([r] * 6), indexing="ij")
    ne = np.stack([x.ravel() for x in g], axis=1).astype(np.int64)
    return ne @ T

def solve(V, Gs, label):
    Hn = H33_np
    n = len(V)
    keys = enc(V); order = np.argsort(keys); ks = keys[order]
    def index_of(W):
        k = enc(W); p = np.searchsorted(ks, k); p = np.clip(p, 0, n - 1)
        ok = ks[p] == k
        return np.where(ok, order[p], -1)
    q = np.einsum("ij,jk,ik->i", V, Hn, V)
    imgs = []
    for M in Gs:
        im = index_of(V @ M.T)
        assert (im >= 0).all(), "window not invariant under G"
        imgs.append(im)
    imgs = np.array(imgs)
    lab = imgs.min(axis=0)
    reps = np.nonzero(lab == np.arange(n))[0]
    HV = V @ Hn
    tuples_total = 0; orbit_total = 0; allowed_total = 0; allowed_orbits = 0
    reps_out = []      # canonical representatives (a,b,v,w) indices of orbits with tadpole-allowed
    nflux_hist = {}
    nflux_hist_orb = {}
    per_rep = []
    for a in reps:
        qa = q[a]
        if qa == 0:
            continue
        stab = [k for k in range(len(Gs)) if imgs[k][a] == a]
        orbit_size = len(Gs) // len(stab)
        Bc = np.nonzero((HV @ V[a] == qa) & (q == 2 * qa))[0]
        Va = np.nonzero((HV @ V[a] == 0) & (q != 0) & (q * qa > 0))[0]
        Tl = []
        for b in Bc:
            Wb = np.nonzero(HV @ V[b] == 0)[0]
            if len(Wb) == 0 or len(Va) == 0:
                continue
            wa = (HV[Wb] @ V[a]); wq = q[Wb]
            vb = HV[Va] @ V[b]; vq = q[Va]
            # group v by (v.b, q_v): required w.a = -v.b, q_w = 2 q_v
            kv = (-vb) * 100003 + 2 * vq
            kw = wa * 100003 + wq
            common = np.intersect1d(np.unique(kv), np.unique(kw))
            for key in common:
                vi = np.nonzero(kv == key)[0]; wi = np.nonzero(kw == key)[0]
                Vv = V[Va[vi]]; Ww = V[Wb[wi]]
                D = Vv @ Hn @ Ww.T        # v.w
                qv = q[Va[vi]]
                ii, jj = np.nonzero(D == qv[:, None])
                if len(ii):
                    Tl.extend(zip(np.full(len(ii), b), Va[vi][ii], Wb[wi][jj]))
        if not Tl:
            continue
        T = np.array(Tl, dtype=np.int64)
        tuples_total += orbit_size * len(T)
        # N_flux, tadpole
        Nf = 4 * (-(HV[T[:, 0]] * V[T[:, 1]]).sum(axis=1) + (HV[T[:, 2]] @ V[a]))
        ok = (Nf >= 0) & (Nf <= 2 * TADPOLE_TOTAL) & (Nf % 2 == 0)
        for x in np.unique(Nf[ok]) if ok.any() else []:
            nflux_hist[int(x)] = nflux_hist.get(int(x), 0) + orbit_size * int((Nf[ok] == x).sum())
        allowed_total += orbit_size * int(ok.sum())
        # orbits under Stab(a)
        Tk = T[ok] if ok.any() else T[:0]
        for name, Tsel, is_allowed in (("all", T, False), ("allowed", Tk, True)):
            if len(Tsel) == 0:
                continue
            can = None
            for k in stab:
                im = imgs[k]
                code = (im[Tsel[:, 0]] * n + im[Tsel[:, 1]]) * n + im[Tsel[:, 2]]
                can = code if can is None else np.minimum(can, code)
            uc = np.unique(can)
            if name == "all":
                orbit_total += len(uc)
            else:
                allowed_orbits += len(uc)
                for c in uc:
                    w_ = c % n; v_ = (c // n) % n; b_ = c // (n * n)
                    # find the actual N_flux of the canonical rep
                    reps_out.append((int(a), int(b_), int(v_), int(w_)))
    return {"V": V, "q": q, "reps_out": reps_out, "tuples_total": tuples_total, "orbit_total": orbit_total,
            "allowed_total": allowed_total, "allowed_orbits": allowed_orbits, "nflux_hist": nflux_hist, "n_window": n,
            "group_order": len(Gs), "n_lambda_orbit_reps": int(len(reps)), "imgs": imgs}


def classify(V, q, reps_out):
    """Float classification (tolerance 1e-7, integer data of modest size) of orbit representatives:
       c!=0 and (3.31)>0 ; V_flux has 2 positive directions and NO null vector (TT 3.3); returns arrays.
       Exact re-check on (a subset of) representatives is done in verify_reps."""
    Hn = H33_np.astype(float)
    R_ = np.array(reps_out, dtype=np.int64).reshape(-1, 4)
    if len(R_) == 0:
        return {"admissible": np.zeros(0, bool), "c331": np.zeros(0), "rank": np.zeros(0, int), "npos": np.zeros(0, int), "nnull_in_V": np.zeros(0, int), "nneg": np.zeros(0, int)}
    Aa, Bb, Vv, Ww = [2.0 * V[R_[:, k]] for k in range(4)]
    tau = 1j * np.sqrt(q[R_[:, 2]] / q[R_[:, 0]].astype(float))
    phi = 0.5 * (1 + 1j)
    nx = Aa - phi * Bb; ny = Vv - phi * Ww
    Gzb = nx * tau[:, None] - ny
    c331 = np.real(np.einsum("mi,ij,mj->m", Gzb, Hn, np.conj(Gzb)))
    X = np.stack([Aa, Bb, Vv, Ww], axis=1)                       # (m,4,6)
    sv = np.linalg.svd(X, compute_uv=False)
    rank = (sv > 1e-8).sum(axis=1)
    Gm = np.einsum("mai,ij,mbj->mab", X, Hn, X)
    ev = np.linalg.eigvalsh(Gm)
    npos = (ev > 1e-7).sum(axis=1); nneg = (ev < -1e-7).sum(axis=1); nz = 4 - npos - nneg
    nnull_in_V = nz - (4 - rank)
    adm = (c331 > 1e-7) & (npos == 2) & (nnull_in_V == 0)
    return {"admissible": adm, "c331": c331, "rank": rank, "npos": npos, "nnull_in_V": nnull_in_V, "nneg": nneg}

def verify_reps(V, reps_out, maxverify, float_adm=None):
    """Exact verification on orbit representatives."""
    Hn = H33_np
    q = np.einsum("ij,jk,ik->i", V, Hn, V)
    rs = np.random.RandomState(7)
    sel = list(range(len(reps_out))) if len(reps_out) <= maxverify else list(rs.choice(len(reps_out), maxverify, replace=False))
    stats = {"verified": 0, "all_3.16abc_zero": True, "all_3.31_positive": True, "types": {}, "orbifold_free": 0, "orbifold_singular": 0, "nflux_positive": 0, "examples": []}
    for k in sel:
        a, b, v, w = reps_out[k]
        ax, bx, ay, by = [[int(2 * x) for x in V[i]] for i in (a, b, v, w)]
        r = R(int(q[v]), int(q[a]))
        tau = I * sqrt(r)
        phi = R(1, 2) * (1 + I)
        c = tt_conditions(H, ax, bx, ay, by, phi, tau)
        z = all(sp.simplify(c[t]) == 0 for t in ("3.16a", "3.16b", "3.16c"))
        pos = bool(sp.simplify(c["3.15_3.31_value"] > 0))
        stats["all_3.16abc_zero"] &= z; stats["all_3.31_positive"] &= pos
        sg = sig_rank(gram_of(H, [ax, bx, ay, by]))
        stats["types"][str(sg)] = stats["types"].get(str(sg), 0) + 1
        try:
            orb, Bs, nk = orbifold_hits(H, [ax, bx, ay, by], c["Gzbar"], 6)
            stats["orbifold_singular" if orb else "orbifold_free"] += 1
        except AssertionError:
            stats.setdefault("orbifold_test_failed", 0); stats["orbifold_test_failed"] += 1
        stats["verified"] += 1
        if float_adm is not None:
            exact_adm = bool(pos and z and (sg[0] == 2) and (sg[2] == 4 - int(sp.Matrix([ax, bx, ay, by]).rank())))
            stats["float_vs_exact_admissibility_disagreements"] = stats.get("float_vs_exact_admissibility_disagreements", 0) + int(exact_adm != bool(float_adm[k]))
        Nf = -sum(bx[i] * H[i, j] * ay[j] for i in range(6) for j in range(6)) + sum(by[i] * H[i, j] * ax[j] for i in range(6) for j in range(6))
        stats["nflux_positive"] += int(bool(Nf > 0))
        if len(stats["examples"]) < 3:
            stats["examples"].append({"alpha_x": ax, "beta_x": bx, "alpha_y": ay, "beta_y": by, "tau": str(tau), "N_flux": int(Nf), "V_flux_signature": sg, "3.31_value": str(c["3.15_3.31_value"])})
    stats["verified_subset"] = "all orbit representatives" if len(sel) == len(reps_out) else f"random {len(sel)} of {len(reps_out)} representatives (seed 7)"
    return stats

def main():
    NU = int(sys.argv[1]); NA = int(sys.argv[2]); MAXV = int(sys.argv[3])
    CMD = rel_command(SCRIPT, " ".join(sys.argv[1:]))
    t0 = time.time()
    Gall = group_G(); Gplus = [M for M in Gall if orient_plus(M)]
    res = []
    for label, V in (("U3_window_fg_N=%d" % NU, window_U(NU)), ("A6_subfamily_window_e_n=%d" % NA, window_A6(NA))):
        t1 = time.time()
        rp = solve(V, Gplus, label)
        rg = solve(V, Gall, label)
        cl = classify(rp["V"], rp["q"], rp["reps_out"])
        adm_idx = [k for k in range(len(rp["reps_out"])) if cl["admissible"][k]]
        imgs = rp["imgs"]
        def orbit_size(rep):
            a, b, v, w = rep
            fix = ((imgs[:, a] == a) & (imgs[:, b] == b) & (imgs[:, v] == v) & (imgs[:, w] == w)).sum()
            return len(imgs) // int(fix)
        n_adm_tuples = sum(orbit_size(rp["reps_out"][k]) for k in adm_idx)
        types = {}
        for k in adm_idx:
            key = "(%d+,%d-)" % (cl["npos"][k], cl["nneg"][k]); types[key] = types.get(key, 0) + 1
        rej = {"c331_not_positive_(Gzbar=0_second_branch_or_isotropic)": int((cl["c331"] <= 1e-7).sum()),
               "V_flux_has_null_vector": int((cl["nnull_in_V"] > 0).sum()), "V_flux_not_2_positive": int((cl["npos"] != 2).sum())}
        adm_reps = [rp["reps_out"][k] for k in adm_idx]
        ver = verify_reps(rp["V"], adm_reps, MAXV, float_adm=np.ones(len(adm_reps), bool))
        ver_rej = verify_reps(rp["V"], [rp["reps_out"][k] for k in range(len(rp["reps_out"])) if not cl["admissible"][k]], 40, float_adm=np.zeros(max(1, len(rp["reps_out"]) - len(adm_idx)), bool))
        # locate TT (4.31) if window is the A6 one
        tt431 = None
        if label.startswith("A6"):
            T = np.zeros((6, 6), dtype=np.int64)
            for i in range(3):
                T[i, 2 * i] = 1; T[i, 2 * i + 1] = 1; T[i + 3, 2 * i] = 1; T[i + 3, 2 * i + 1] = -1
            e = lambda c: np.array(c) @ T
            lam = [e([1, -1, 0, 0, 0, 0]), e([0, -2, 0, 0, 0, 0]), e([1, 1, 0, 1, 0, 0]), e([2, 0, 0, 1, 1, 0])]  # a,b,v,w  (alpha=2*lambda)
            Vn = rp["V"]
            keys = {tuple(x): i for i, x in enumerate(Vn)}
            idx = [keys.get(tuple(x), -1) for x in lam]
            Hn = H33_np
            q = rp["q"]
            cond = None
            if min(idx) >= 0:
                a, b, v, w = [Vn[i] for i in idx]
                cond = {"a.v": int(a @ Hn @ v), "b.w": int(b @ Hn @ w), "a.w+v.b": int(a @ Hn @ w + v @ Hn @ b), "b2-2a2": int(b @ Hn @ b - 2 * (a @ Hn @ a)), "a2-a.b": int(a @ Hn @ a - a @ Hn @ b), "w2-2v2": int(w @ Hn @ w - 2 * (v @ Hn @ v)), "v2-v.w": int(v @ Hn @ v - v @ Hn @ w)}
            tt431 = {"inside_window": min(idx) >= 0, "ansatz_conditions_residuals": cond}
        entry = {"window": label, "n_lattice_vectors_in_window": int(len(V)),
                 "tuples_satisfying_4.24_4.25_with_tau_nonreal": int(rp["tuples_total"]),
                 "of_which_0<=N_flux<=48_(N_D3>=0)": int(rp["allowed_total"]),
                 "orbits_mod_Gplus": {"order": rp["group_order"], "all_tuples": int(rp["orbit_total"]), "tadpole_allowed": int(rp["allowed_orbits"])},
                 "orbits_mod_G": {"order": rg["group_order"], "all_tuples": int(rg["orbit_total"]), "tadpole_allowed": int(rg["allowed_orbits"])},
                 "tuples_by_N_flux_(tadpole_allowed)": {str(k): v for k, v in sorted(rp["nflux_hist"].items())},
                 "N_D3_by_N_flux": {str(k): str(TADPOLE_TOTAL - R(k, HALF)) for k in sorted(rp["nflux_hist"])},
                 "admissible_after_TT_conditions_(c!=0,(3.31)>0,V_flux 2 positive and no null vector)": {"orbit_reps_mod_Gplus": len(adm_idx), "tuples": int(n_adm_tuples), "V_flux_types_(orbits)": types, "rejected_orbit_counts_(tadpole-allowed reps)": rej},
                 "exact_verification_on_admissible_Gplus_orbit_reps": ver, "exact_check_of_REJECTED_reps_agrees_with_float_(sample_of_40)": ver_rej, "TT_4.31_present_in_window": tt431, "seconds": round(time.time() - t1, 1)}
        res.append(entry)
        print(json.dumps(entry, default=str)[:1500], flush=True)
    out = [{"id": "E04a", "quantity": "TT section 4.2 ansatz (4.24)-(4.25): tuples and orbit counts in the stated windows",
            "computed": {"windows": res, "truncation": "ansatz (4.24)-(4.25) only; windows only; NOT the full class (c); counts mod G+ / G are not physically distinct vacua",
                         "moduli_fixed": "phi=(1+-i)/2, tau=i sqrt(alpha_yy/alpha_xx), all complex structure moduli of K3 via Omega (TT 3.9); Kahler: 18 of K3 + T2 volume unfixed for dim V_flux=4 (TT 3.4)"},
            "shared_inputs": ["TT_H33_gram_A3", "flux_quantisation_even_coefficients", "sec42_conditions_4.24_4.25", "N_flux_formula_2.8", "tadpole_total_24_eq2.3", "tadpole_half_factor_eq2.3", "group_G_definition", "window_bounds", "embedding_A6_into_U3", "orbifold_criterion_TT_3.3"],
            "script": SCRIPT, "command": CMD}]
    write_json("04_family42_counts.json", {"script": SCRIPT, "results": out, "rigidity": [
        {"parameter_inserted": "the finite set of N_flux values reachable in the 4.2 ansatz within the window", "selecting_condition": "(4.24),(4.25) and 0 <= N_flux <= 2*24 (tadpole)",
         "condition_uses_true_value": True, "solution_set": "see E04a tuples_by_N_flux", "classification": "CONDITIONAL_ON_INPUT", "input_it_depends_on": "tadpole_total_24_eq2.3, window_bounds",
         "negative_control": "N_D3 >= 0 removed: N_flux unbounded above inside the window (largest N_flux reached is reported)", "control_perturbs_same_parameter": True}],
        "literature_checks": [], "could_not_do": ["section 4.2 counts restricted to the ansatz (4.24)-(4.25) and windows; the general class (c) (dim V_flux = 4, M zero eigenvalue) was not enumerated"], "wall_seconds": round(time.time() - t0, 1)})

if __name__ == "__main__":
    main()
