"""Skeptic independent GUDHI computation for Track D (own triangulation code, no D-tda imports).
Kuhn/Freudenthal triangulation of T^4 = (Z/N)^4 (4-simplices v, v+e_s1, v+e_s1+e_s2, ...),
quotient by x -> -x (vertex canonical rep = lexicographic min of {x, -x}), then
  X  = T^4/+-  (full quotient),  U = X minus open stars of the 16 fixed vertices.
Checks: simpliciality of the quotient (every image simplex has 5 distinct vertices and the
f-vector halves as expected), Betti numbers of X and U over Z/3, and the Mayer-Vietoris
bookkeeping b2(K3) = b2(U) + 16*b2(S^2), chi(K3) = chi(U) + 16*chi(S^2) - 16*chi(RP^3).
Command: cd audit/k3t2_rigidity_v2/skeptic_math && <venv python> skD_kummer.py 6
"""
import sys, json, os, itertools
import gudhi
N = int(sys.argv[1]) if len(sys.argv) > 1 else 6
d = 4
def canon(x):
    y = tuple((-c) % N for c in x)
    return min(x, y)
def enc(x):
    r = 0
    for c in x: r = r * N + c
    return r
tops = set(); ntop_cover = 0; bad = 0
for v in itertools.product(range(N), repeat=d):
    for perm in itertools.permutations(range(d)):
        cur = list(v); simp = [tuple(cur)]
        for i in perm:
            cur[i] = (cur[i] + 1) % N; simp.append(tuple(cur))
        ntop_cover += 1
        img = frozenset(enc(canon(p)) for p in simp)
        if len(img) != d + 1: bad += 1
        tops.add(img)
fixed = {enc(p) for p in itertools.product((0, N // 2), repeat=d)}
X = gudhi.SimplexTree(); Ucx = gudhi.SimplexTree()
for s in tops:
    X.insert(sorted(s))
    rest = sorted(set(s) - fixed)
    if rest: Ucx.insert(rest)
def betti(st):
    st.compute_persistence(homology_coeff_field=3, persistence_dim_max=True)
    b = st.betti_numbers(); return b + [0] * (5 - len(b))
def fvec(st):
    f = [0] * 5
    for s, _ in st.get_simplices(): f[len(s) - 1] += 1
    return f
fX, fU = fvec(X), fvec(Ucx)
bX, bU = betti(X), betti(Ucx)
chi = lambda f: sum((-1) ** i * x for i, x in enumerate(f))
b2K3 = bU[2] + 16 * 1
chiK3 = chi(fU) + 16 * 2 - 16 * 0
out = {"N": N, "cover_top_simplices": ntop_cover, "quotient_top_simplices": len(tops),
       "degenerate_images": bad, "fixed_vertices": len(fixed),
       "f_X": fX, "f_U": fU, "betti_X_Z3": bX, "betti_U_Z3": bU,
       "chi_X_fvec": chi(fX), "chi_U_fvec": chi(fU), "chi_U_betti": sum((-1) ** i * x for i, x in enumerate(bU)),
       "b2_K3_MV": b2K3, "chi_K3_MV": chiK3}
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), f"skD_kummer_N{N}_results.json"), "w"), indent=1)
print(json.dumps(out))
