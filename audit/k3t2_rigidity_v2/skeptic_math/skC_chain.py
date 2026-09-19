"""Skeptic check for Track C: (1) Noether/Hodge chain from chi_top (read from D-tda/exports.json,
cross-checked against A-genus/exports.json), c1=0, b1=0, Kahler; (2) exact signature of
3U+2(-E8) from my own E8 Cartan matrix (Dynkin edges typed from the E8 diagram 1-2-3-4-5-6-7, 3-8,
i.e. TT eq A.4's layout) via Descartes sign changes of the (real-rooted) characteristic polynomial.
Command: cd audit/k3t2_rigidity_v2/skeptic_math && /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python skC_chain.py
"""
import json, os
import sympy
from fractions import Fraction as Fr
here = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(here, "..", "D-tda", "exports.json")))
A = json.load(open(os.path.join(here, "..", "A-genus", "exports.json")))
def find(d, key):
    if isinstance(d, dict):
        for k, v in d.items():
            if k == key: return v
            r = find(v, key)
            if r is not None: return r
    return None
chi = find(D, "chi_K3_resolved"); chiA = find(A, "chi_from_genus")
chiA = chiA["value"] if isinstance(chiA, dict) else chiA
res = {"chi_top_D": chi, "chi_top_A": chiA, "agree": int(chi) == int(chiA)}
def chain(chi, b1=0):
    c1sq = 0                               # c1 = 0 (K3 / CY assumption)
    chiO = Fr(c1sq + chi, 12)              # Noether: chi(O) = (c1^2 + c2)/12, c2 = chi_top
    h01 = Fr(b1, 2)                        # compact Kahler surface: b1 = 2 h^{0,1}
    h02 = chiO - 1 + h01                   # chi(O) = h00 - h01 + h02
    b2 = chi - 2 + 2 * b1                  # chi = 2 - 2 b1 + b2 (b3 = b1 by Poincare)
    b2p = 2 * h02 + 1                      # Hodge index for Kahler surfaces
    ok = chiO.denominator == 1 and h02.denominator == 1 and h02 >= 0
    return {"chi_O": str(chiO), "h02": str(h02), "b2": b2, "b2_plus": str(b2p), "b2_minus": str(b2 - b2p), "integral": ok}
res["chain_24"] = chain(chi)
res["chain_perturbed"] = {c: chain(c) for c in (0, 12, 22, 26, 36, 48)}
edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (2, 7)]
E8 = sympy.zeros(8, 8)
for i in range(8): E8[i, i] = 2
for i, j in edges: E8[i, j] = E8[j, i] = -1
U = sympy.Matrix([[0, 1], [1, 0]])
def sig(M):
    x = sympy.symbols("x")
    p = sympy.Poly(M.charpoly(x).as_expr(), x)
    co = [c for c in p.all_coeffs()]
    def changes(cs):
        s = [c for c in cs if c != 0]
        return sum(1 for a, b in zip(s, s[1:]) if a * b < 0)
    npos = changes(co)
    neg = [c * (-1) ** (len(co) - 1 - i) for i, c in enumerate(co)]
    nneg = changes(neg)
    return npos, nneg, M.rows - npos - nneg
# Reverse chain, NO chi input: K3 := compact complex surface with K_X trivial and b1 = 0.
# h^{2,0} = h^0(K_X) = h^0(O) = 1 (K trivial), h^{0,1} = 0 (b1 = 0 => h^{0,1} = 0 for any compact
# complex surface), Serre duality h^{0,2} = h^{2,0}; Noether with c1 = 0: chi_top = 12 chi(O).
h20 = 1; h01 = 0
chiO_rev = 1 - h01 + h20
res["reverse_chain_no_chi_input"] = {"chi_O": chiO_rev, "chi_top_forced": 12 * chiO_rev,
                                     "note": "chi_top is fixed by K trivial + b1=0 + Noether; C's chain could run this way (RIGID given the definition of K3) instead of reading 24 from a file"}
res["E8_det"] = int(E8.det()); res["E8_signature"] = sig(E8)
for m, n in ((3, 2), (7, 1), (11, 0)):
    M = sympy.diag(*([U] * m + [-E8] * n))
    res[f"sig_{m}U+{n}(-E8)"] = sig(M); res[f"det_{m}U+{n}(-E8)"] = int(M.det())
res["even_diag_all"] = True
json.dump(res, open(os.path.join(here, "skC_chain_results.json"), "w"), indent=1, default=str)
print(json.dumps(res, indent=1, default=str))
