"""Track F part 4: Kummer code K, extended Golay code G24, M24 = <PSL(2,23), delta>, cycle shapes, and 16-subsets of coordinates.
Run: cd audit/k3t2_rigidity_v3/F-whichk3 && <venv python> f4_codes.py 40000
argument: number of random M24 elements sampled for cycle shapes (committed run: 40000).
Everything exact (bit masks over GF(2), permutations, sympy Schreier-Sims).
"""
import sys, itertools, math, random, json
from collections import Counter
from sympy.combinatorics import Permutation, PermutationGroup
from common import *

NSAMP = int(sys.argv[1]) if len(sys.argv) > 1 else 40000
random.seed(20260919)   # sampling seed (does not influence exact results)

def weight(x): return bin(x).count("1")
def span(gens):
    words = {0}
    for g in gens:
        words |= {w ^ g for w in words}
    return words
def rank2(vs):
    basis = []
    for v in vs:
        for b in basis: v = min(v, v ^ b)
        if v: basis.append(v)
    return len(basis)
def wenum(words): return dict(sorted(Counter(weight(w) for w in words).items()))

# ---------- (a) Kummer code: affine hyperplanes of F_2^4
pts = list(range(16))
hyper = []
for a in range(1, 16):
    for c in (0, 1):
        hyper.append(sum(1 << p for p in pts if (weight(a & p) & 1) == c))
Kwords = span(hyper)
resK = {"n_hyperplanes": len(hyper), "dim": rank2(hyper), "n_words": len(Kwords), "weight_enumerator": wenum(Kwords)}
print("K:", resK)

# ---------- (b) extended Golay code from quadratic residues mod 23
p = 23
Q = {(x * x) % p for x in range(1, p)}
Nn = set(range(1, p)) - Q
INF = 23
def shift_word(S, k): return sum(1 << ((s + k) % p) for s in S)
ones23 = (1 << p) - 1
variants = {"Q": Q, "N": Nn, "Q+0": Q | {0}, "N+0": Nn | {0}}
golay_variants = {}
for name, S in variants.items():
    gens = [shift_word(S, k) for k in range(p)] + [ones23]
    r = rank2(gens)
    ws = span(gens) if r <= 13 else set()
    # extend by overall parity bit at INF
    ext = {w | ((weight(w) & 1) << INF) for w in ws}
    golay_variants[name] = {"dim23": r, "ext": ext}
    print(name, "dim", r)

# PSL(2,23) generators on P^1(F_23) with point 23 = infinity
inv = {x: pow(x, -1, p) for x in range(1, p)}
def perm_from(f): return tuple(f(x) for x in range(24))
t_ = perm_from(lambda x: INF if x == INF else (x + 1) % p)
m_ = perm_from(lambda x: INF if x == INF else (4 * x) % p)          # 4 = 2^2 generates the squares
s_ = perm_from(lambda x: 0 if x == INF else (INF if x == 0 else (-inv[x]) % p))
def apply(perm, w):
    out = 0
    for i in range(24):
        if (w >> i) & 1: out |= 1 << perm[i]
    return out
def preserves(perm, code, basis=None):
    return all(apply(perm, w) in code for w in (basis if basis is not None else code))

chosen, golay_info = None, {}
for name, d in golay_variants.items():
    ext = d["ext"]
    if not ext: continue
    ws = wenum(ext)
    inv_ok = all(preserves(g, ext) for g in (t_, m_, s_)) if len(ext) == 4096 else False
    golay_info[name] = {"n_words": len(ext), "weight_enumerator": ws, "PSL_gens_preserve": inv_ok}
    print("variant", name, golay_info[name])
    if chosen is None and len(ext) == 4096 and inv_ok and ws.get(8) == 759:
        chosen = name
assert chosen is not None
G24 = golay_variants[chosen]["ext"]
octads = sorted(w for w in G24 if weight(w) == 8)
resG = {"construction_chosen": chosen, "variants_tested": golay_info, "n_words": len(G24), "weight_enumerator": wenum(G24), "n_octads": len(octads)}
gen12 = []
b = []
for w in sorted(G24):
    if rank2(b + [w]) > len(b): b.append(w)
G24basis = b
print("G24 dim", len(G24basis), wenum(G24))

# ---------- (c) M24
SP = lambda x: Permutation(list(x))
PSL = PermutationGroup([SP(t_), SP(m_), SP(s_)])
psl_order = PSL.order()
# search delta in the family x -> c1 x^k (x in Q), c2 x^k (x in N), 0,inf fixed
found, outside = [], []
for k in range(1, 23):
    for c1 in range(1, 23):
        for c2 in range(1, 23):
            img = [0] * 24
            img[0] = 0; img[INF] = INF
            for x in range(1, p):
                img[x] = (c1 if x in Q else c2) * pow(x, k, p) % p
            if len(set(img)) != 24: continue
            pm = tuple(img)
            if preserves(pm, G24, G24basis):
                found.append((k, c1, c2))
                if not PSL.contains(SP(pm)): outside.append((k, c1, c2, pm))
print("family members preserving G24:", len(found), "outside PSL(2,23):", len(outside))
k, c1, c2, delta = outside[0]
M24 = PermutationGroup([SP(t_), SP(m_), SP(s_), SP(delta)])
m24_order = M24.order()
print("|PSL(2,23)| =", psl_order, " |<PSL,delta>| =", m24_order)
gens_all = [t_, m_, s_, delta]
assert all(preserves(g, G24, G24basis) for g in gens_all)
# transitivity on octads
orb = {octads[0]}; stack = [octads[0]]
while stack:
    o = stack.pop()
    for g in gens_all:
        n_ = apply(g, o)
        if n_ not in orb: orb.add(n_); stack.append(n_)
basic_orbit_lengths = [len(o) for o in M24.basic_orbits]

def cyc_type(perm):
    seen, ct = [False] * 24, []
    for i in range(24):
        if not seen[i]:
            l, j = 0, i
            while not seen[j]: seen[j] = True; j = perm[j]; l += 1
            ct.append(l)
    return tuple(sorted(ct))
def lcm(a, b): return a * b // math.gcd(a, b)
def order_of_ct(ct):
    o = 1
    for c in ct: o = lcm(o, c)
    return o
def mul(a, b): return tuple(b[a[i]] for i in range(24))   # a then b
regs = [random.choice(gens_all) for _ in range(20)]
regs = list(gens_all) + [gens_all[i % 4] for i in range(16)]
acc = tuple(range(24))
def pr_step():
    global acc
    i, j = random.sample(range(len(regs)), 2)
    regs[i] = mul(regs[i], regs[j]) if random.random() < .5 else mul(regs[j], regs[i])
    acc = mul(acc, regs[i])
    return acc
for _ in range(500): pr_step()
shapes = {}
for _ in range(NSAMP):
    g = pr_step()
    ct = cyc_type(g)
    o = order_of_ct(ct)
    shapes.setdefault(o, Counter())[ct] += 1
def fmt_ct(ct):
    c = Counter(ct); return " ".join(f"{l}^{c[l]}" for l in sorted(c))
req = [2, 3, 4, 5, 6, 7, 8, 11, 12, 14, 15, 21, 23]
shape_tab = {str(o): {fmt_ct(ct): n for ct, n in sorted(shapes[o].items())} for o in sorted(shapes)}
n_types = sum(len(v) for v in shapes.values())
print("orders seen:", sorted(shapes), " distinct cycle types:", n_types)
for o in req: print(o, shape_tab.get(str(o)))
resM = {"psl_order": psl_order, "delta_family": {"parametrisation": "x->c1 x^k on QR, c2 x^k on NQR, 0 and inf fixed; k in 1..22, c1,c2 in 1..22",
        "n_members_preserving_G24": len(found), "n_outside_PSL": len(outside), "delta_chosen(k,c1,c2)": (k, c1, c2)},
        "order_of_generated_group": m24_order, "generators_preserve_G24": True, "octad_orbit_size": len(orb),
        "schreier_sims_basic_orbit_lengths": basic_orbit_lengths,
        "sample_size": NSAMP, "shapes_by_order_counts": shape_tab, "n_distinct_cycle_types": n_types, "required_orders_present": {str(o): (o in shapes) for o in req}}

# ---------- (d) 16-subsets
import numpy as np
oc = np.array(octads, dtype=np.int64)
Ts = list(itertools.combinations(range(24), 8))
Tmask = np.array([sum(1 << i for i in c) for c in Ts], dtype=np.int64)
cnt_in = np.zeros(len(Tmask), dtype=np.int64)
CH = 4000
for s0 in range(0, len(Tmask), CH):
    blk = Tmask[s0:s0 + CH]
    cnt_in[s0:s0 + CH] = ((blk[:, None] & oc[None, :]) == 0).sum(axis=1)
dist_octads_inside = dict(sorted(Counter(cnt_in.tolist()).items()))
# dimension of shortened code {c in G24: supp c subset S}, S = complement of T: 12 - rank of restriction to T
dims = Counter(); good = []
for T, cin in zip(Tmask.tolist(), cnt_in.tolist()):
    r = rank2([g & T for g in G24basis])
    d = 12 - r
    dims[d] += 1
    if d >= 5: good.append(T)
FULL = (1 << 24) - 1
n_equiv = 0; fails = 0; wenums = Counter()
for T in good:
    S = FULL ^ T
    Sidx = [i for i in range(24) if (S >> i) & 1]
    C = [w for w in G24 if (w & T) == 0]
    wenums[str(wenum(C))] += 1
    if len(C) != 32 or S not in C: fails += 1; continue
    base = [S]
    for w in sorted(C):
        if rank2(base + [w]) > len(base): base.append(w)
    vec = {i: tuple((base[j] >> i) & 1 for j in range(1, 5)) for i in Sidx}
    if len(set(vec.values())) != 16: fails += 1; continue
    # map the 30 weight-8 words to affine hyperplanes of F_2^4 under the bijection point -> vec
    hs = set()
    for a in range(1, 16):
        for c in (0, 1):
            m = 0
            for i in Sidx:
                v = sum(vec[i][j] << j for j in range(4))
                if (weight(a & v) & 1) == c: m |= 1 << i
            hs.add(m)
    w8 = {w for w in C if weight(w) == 8}
    if hs == w8 and len(w8) == 30: n_equiv += 1
    else: fails += 1
resD = {"n_16_subsets_total": len(Ts), "distribution_of_number_of_octads_inside_S": dist_octads_inside,
        "dim_of_shortened_code_distribution": dict(sorted(dims.items())), "n_S_with_dim>=5": len(good),
        "weight_enumerators_of_those": dict(wenums), "n_S_equivalent_to_K": n_equiv, "n_failures": fails,
        "test": "shortened code has dim 5, contains all-ones on S, its 16 columns (1,v_i) are distinct in F_2^4 (bijection S->F_2^4), and its 30 weight-8 words equal the 30 affine hyperplanes",
        "S_is_complement_of_an_octad": all((FULL ^ T) is not None and (T in set(octads)) for T in good)}
print(resD)

add_inputs([
 {"name": "QR_golay_construction", "value": "extended QR code mod 23 built from shifts of an indicator of Q, N, Q+0, N+0 with all-ones; the variant with dim 12, 759 octads and PSL(2,23)-invariance is used (tested, not typed)", "tier": "B", "why": "task item 4(b)"},
 {"name": "delta_search_family", "value": "x->c1 x^k (QR), c2 x^k (NQR), 0/inf fixed", "tier": "B", "why": "search space for the extra M24 generator; the solution is returned by the scan"},
 {"name": "Aut_G24_is_M24", "value": "Aut(G24)=M24 of order 244823040 (used only to say the generated group IS the automorphism group; we compute its order and check it preserves G24, not that it is the full automorphism group)", "tier": "L", "why": "FROM MEMORY"},
 {"name": "sampling_seed", "value": 20260919, "tier": "B", "why": "cycle-shape sampling only"},
])
add_results([
 {"id": "F4a_kummer_code", "quantity": "dim and weight enumerator of the code spanned by affine hyperplanes of F_2^4", "computed": resK, "shared_inputs": [], "script": RELDIR + "/f4_codes.py", "command": cmd("f4_codes.py", str(NSAMP)), "tier": "B"},
 {"id": "F4b_golay", "quantity": "extended Golay code from QR mod 23: weight enumerator, octads", "computed": resG, "shared_inputs": ["QR_golay_construction"], "script": RELDIR + "/f4_codes.py", "command": cmd("f4_codes.py", str(NSAMP)), "tier": "B"},
 {"id": "F4c_M24", "quantity": "<PSL(2,23), delta> order (Schreier-Sims) and cycle shapes of sampled elements", "computed": resM, "shared_inputs": ["delta_search_family", "Aut_G24_is_M24", "sampling_seed"], "script": RELDIR + "/f4_codes.py", "command": cmd("f4_codes.py", str(NSAMP)), "tier": "B (that this group is the full Aut(G24) is tier L)"},
 {"id": "F4d_16_subsets", "quantity": "16-subsets S of the 24 coordinates whose shortened Golay code is equivalent to the Kummer code", "computed": resD, "shared_inputs": ["QR_golay_construction"], "script": RELDIR + "/f4_codes.py", "command": cmd("f4_codes.py", str(NSAMP)), "tier": "B"},
])
add_rigidity([
 {"id": "F4d_subset", "parameter_inserted": "choice of 16 of 24 coordinates", "selecting_condition": "shortened code of G24 on S is equivalent to the Kummer code (dim 5, columns distinct, 30 weight-8 words = affine hyperplanes)",
  "condition_uses_true_value": False, "solution_set": f"{n_equiv} subsets of {len(Ts)}", "classification": "RIGID_GIVEN_DEFINITION",
  "input_it_depends_on": "QR_golay_construction", "negative_control": f"the other {len(Ts) - n_equiv} 16-subsets have shortened code of dimension 4 (see dim distribution); all are realizable subsets", "control_perturbs_same_parameter": True},
])
add_exports({"M24_cycle_shapes": {str(o): sorted(fmt for fmt in shape_tab[str(o)]) for o in req if str(o) in shape_tab},
             "M24_cycle_shapes_all_orders": {o: sorted(v) for o, v in shape_tab.items()},
             "M24_order": m24_order, "G24_n_octads": len(octads), "K_weight_enumerator": resK["weight_enumerator"], "F4d_n_subsets_equiv_to_K": n_equiv})
# stash the objects for part 5 (small)
json.dump({"octads": octads, "gens": gens_all, "G24basis": G24basis}, open(HERE / "_m24_data.json", "w"))
print("done")
