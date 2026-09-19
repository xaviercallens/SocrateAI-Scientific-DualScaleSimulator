"""Skeptic (math lens) check of the PREMISES of Track D's Mayer-Vietoris step, at a given grid N.
Own quotient construction (same Kuhn/Freudenthal recipe as skD_kummer.py, no D-tda imports).

What the MV step needs (topological input, tier L): K3 = U  u  16 x D(O(-2)), glued along 16 copies of
RP^3 = boundary of the disk bundle D(O(-2)) ~ S^2, where U ~ (T^4/+-) minus 16 points.
Premises checked here, per N:
 P1 open-star disjointness: no quotient simplex contains two fixed vertices. This is what makes
    U := X minus the open stars of the 16 fixed vertices a deformation retract of X minus 16 points.
    (Closed-star disjointness, which Track D tested, is sufficient but NOT necessary; reported too.)
 P2 each fixed-vertex link has the homology of RP^3: Z/2 Betti (1,1,1,1), Z/3 Betti (1,0,0,1).
    Over Z/3 the gluing locus therefore has H_1 = H_2 = 0, which makes the MV sequence split in degree 2:
      0 = H_2(dA) -> H_2(U) + H_2(D)^16 -> H_2(K3) -> H_1(dA) = 0   =>  b2(K3) = b2(U) + 16.
 P3 the rank of the connecting-side map H_3(dA) -> H_3(U) (dA = 16 links) is COMPUTED, not assumed,
    from the MV sequence of X = U u 16 cones (cones acyclic):
      0 -> H_4(X) -> H_3(dA) -> H_3(U) -> H_3(X) -> H_2(dA) = 0
    so rank = 16 - b4(X), and b3(U) = rank + b3(X) is a consistency check. The SAME map appears in the
    K3 gluing (H_3(D) = H_4(D) = 0), hence b3(K3) = b3(U) - rank and b4(K3) = 16 - rank.
All over Z/3 (and Z/5 as a cross-field check). Exact integer arithmetic; GUDHI for Betti numbers.
Command: cd audit/k3t2_rigidity_v2/skeptic_math && \
  /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python skD_premises.py 4
  (and the same with argument 6)
"""
import sys, json, os, itertools
import gudhi

N = int(sys.argv[1]) if len(sys.argv) > 1 else 6
d = 4


def canon(x):
    y = tuple((-c) % N for c in x)
    return min(x, y)


tops = set()
for v in itertools.product(range(N), repeat=d):
    for perm in itertools.permutations(range(d)):
        cur = list(v); simp = [tuple(cur)]
        for i in perm:
            cur[i] = (cur[i] + 1) % N; simp.append(tuple(cur))
        img = frozenset(canon(p) for p in simp)
        assert len(img) == d + 1
        tops.add(img)
fixed = [canon(p) for p in itertools.product((0, N // 2), repeat=d)]
fset = set(fixed)
assert len(fset) == 16

# P1
two_fixed = sum(1 for s in tops if len(s & fset) >= 2)
closed_star = {p: set().union(*[s for s in tops if p in s]) for p in fixed}
closed_overlap = sum(1 for p, q in itertools.combinations(fixed, 2) if closed_star[p] & closed_star[q])


def betti(simplices, field, dim):
    st = gudhi.SimplexTree()
    for s in simplices:
        st.insert(sorted(s))
    st.compute_persistence(homology_coeff_field=field, persistence_dim_max=True)
    b = st.betti_numbers()
    return b + [0] * (dim + 1 - len(b))


# encode vertices as ints for gudhi
idx = {}
def enc(s):
    return [idx.setdefault(v, len(idx)) for v in s]

# P2: links of fixed vertices
link_betti = {}
for p in fixed:
    lk = [enc(s - {p}) for s in tops if p in s]
    link_betti[str(p)] = {"Z2": betti(lk, 2, 3), "Z3": betti(lk, 3, 3)}
P2 = all(v["Z2"] == [1, 1, 1, 1] and v["Z3"] == [1, 0, 0, 1] for v in link_betti.values())

# P3: X and U
Xs = [enc(s) for s in tops]
Us = [enc(s - fset) for s in tops if s - fset]
res_fields = {}
for F in (3, 5):
    bX = betti(Xs, F, 4); bU = betti(Us, F, 4)
    rank = 16 - bX[4]
    res_fields[f"Z{F}"] = {
        "betti_X": bX, "betti_U": bU,
        "rank_H3(links)->H3(U)": rank,
        "consistency_b3U_eq_rank_plus_b3X": bU[3] == rank + bX[3],
        "K3_betti_from_MV": [bU[0], bU[1], bU[2] + 16, bU[3] - rank, 16 - rank],
    }
out = {
    "N": N, "command": f"cd audit/k3t2_rigidity_v2/skeptic_math && /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python skD_premises.py {N}",
    "quotient_top_simplices": len(tops),
    "P1_simplices_containing_two_fixed_vertices": two_fixed,
    "P1_open_stars_disjoint": two_fixed == 0,
    "closed_star_overlapping_pairs_(Track_D_criterion,_sufficient_not_necessary)": closed_overlap,
    "P2_all_16_links_have_RP3_homology_Z2_Z3": P2,
    "P2_link_betti_sample": link_betti[str(fixed[0])],
    "P3": res_fields,
    "reading": "MV b2(K3)=b2(U)+16 is the exact output of the Z/3 MV sequence (H_1, H_2 of the gluing locus vanish mod 3); "
               "the only nontrivial map, H_3(links)->H_3(U), has its rank computed from X's Betti numbers. "
               "Premises P1,P2 are what is needed; closed-star disjointness is not.",
}
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), f"skD_premises_N{N}_results.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
