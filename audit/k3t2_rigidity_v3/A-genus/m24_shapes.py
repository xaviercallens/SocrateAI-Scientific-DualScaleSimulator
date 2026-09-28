#!/usr/bin/env python3
"""
Track A v3 helper: build M24 as a permutation group on 24 points from the
extended binary Golay code (quadratic residues mod 23) and read cycle shapes
of elements of prime-power orders. Nothing is typed from a character table.

Steps (all computed):
  1. Golay code G24 = span of the 23 cyclic shifts of QR-indicator words with the
     point at infinity appended, plus the all-ones word. Check dim 12, weight
     distribution.
  2. Generators alpha: x->x+1, beta: x->2x, gamma: x->-1/x, delta: x->x^3/9 on
     NON-residues, 9x^3 on residues (the first tried assignment, opposite, failed the code-preservation check and was swapped) (points 0..22, infinity=23). Each generator
     is checked to preserve the code (image of every basis word is a codeword).
  3. Group order via Schreier-Sims (sympy); compare (expected field) with 244823040.
  4. Sample random elements; for each element of order o take powers to get
     prime-order and order-4 elements; record cycle shapes and fixed point counts,
     with class sizes |G|/|C(g)| (centralizer order by sympy) for one
     representative per (order, shape).

Command (from repo root the script lives at audit/k3t2_rigidity_v3/A-genus):
  cd audit/k3t2_rigidity_v3/A-genus && \
  /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python m24_shapes.py --samples 4000 --seed 1
Writes m24_shapes.json next to the script.
"""
import argparse, json, random
from pathlib import Path
from collections import Counter
from sympy.combinatorics import Permutation, PermutationGroup

HERE = Path(__file__).resolve().parent
P = 23
INF = 23


def golay():
    qr = sorted({(x * x) % P for x in range(1, P)})
    basis = []
    for s in range(P):
        w = 0
        for r in qr:
            w |= 1 << ((r + s) % P)
        w |= 1 << INF  # infinity appended
        basis.append(w)
    basis.append((1 << 24) - 1)
    # row-reduce
    red = []
    for w in basis:
        for r in red:
            w = min(w, w ^ r)
        if w:
            red.append(w)
    return qr, red


def span_all(red):
    words = [0]
    for r in red:
        words += [w ^ r for w in words]
    return words


def inv_mod(a):
    return pow(a, P - 2, P)


def make_gens(qr):
    qrs = set(qr)
    def alpha(x): return INF if x == INF else (x + 1) % P
    def beta(x): return INF if x == INF else (2 * x) % P
    def gamma(x):
        if x == INF: return 0
        if x == 0: return INF
        return (-inv_mod(x)) % P
    def delta(x):
        if x in (0, INF): return x
        if x in qrs: return (9 * pow(x, 3, P)) % P
        return (pow(x, 3, P) * inv_mod(9)) % P
    return [alpha, beta, gamma, delta]


def permute_word(w, f):
    out = 0
    for i in range(24):
        if w >> i & 1:
            out |= 1 << f(i)
    return out


def cyc_shape(perm):
    seen = set(); shape = Counter()
    for i in range(24):
        if i in seen: continue
        j = i; L = 0
        while j not in seen:
            seen.add(j); j = perm[j]; L += 1
        shape[L] += 1
    return dict(sorted(shape.items()))


def shape_str(sh):
    return " ".join(f"{a}^{m}" for a, m in sh.items())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=4000)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    qr, red = golay()
    words = span_all(red)
    wd = Counter(bin(w).count("1") for w in words)
    gens = make_gens(qr)
    gens_perm = [[g(i) for i in range(24)] for g in gens]
    for gp in gens_perm:
        assert sorted(gp) == list(range(24))
    wordset = set(words)
    preserved = [all(permute_word(w, g) in wordset for w in red) for g in gens]
    G = PermutationGroup([Permutation(gp) for gp in gens_perm])
    order = G.order()
    rng = random.Random(a.seed)
    found = {}   # (order, shape_str) -> dict
    counts = Counter()
    els = [G.random() for _ in range(a.samples)]
    for e in els:
        o = e.order()
        divs = [d for d in range(1, o + 1) if o % d == 0]
        for d in divs:
            h = e ** (o // d)  # element of order d
            arr = h.array_form
            sh = cyc_shape(arr)
            key = (d, shape_str(sh))
            counts[key] += 1
            if key not in found:
                found[key] = {"order": d, "shape": {str(k): v for k, v in sh.items()},
                              "fixed_points": sh.get(1, 0), "perm": arr}
    classes = []
    for key, v in sorted(found.items()):
        if v["order"] == 1: continue
        h = Permutation(v["perm"])
        C = G.centralizer(h)
        v = dict(v)
        v["centralizer_order"] = C.order()
        v["class_size_shape_representative"] = order // C.order()
        v["shape_str"] = key[1]
        v["sample_frequency"] = counts[key]
        classes.append(v)
    out = {
        "golay": {"dim": len(red), "n_words": len(words),
                  "weight_distribution": {str(k): v for k, v in sorted(wd.items())}},
        "generators_preserve_code": preserved,
        "group_order_computed": int(order),
        "group_order_expected_field_source": "244823040 = |M24| (literature, FROM MEMORY)",
        "group_order_matches_expected": int(order) == 244823040,
        "n_random_samples": a.samples,
        "note": "class_size is size of the conjugacy class of the representative; several M24 classes can share one cycle shape (e.g. 7A/7B), in which case this is one class only",
        "classes": classes,
    }
    (HERE / "m24_shapes.json").write_text(json.dumps(out, indent=1))
    print(json.dumps({k: out[k] for k in ("golay", "generators_preserve_code", "group_order_computed", "group_order_matches_expected")}))
    for c in classes:
        print(c["order"], c["shape_str"], "fix", c["fixed_points"], "C", c["centralizer_order"], "size", c["class_size_shape_representative"], "freq", c["sample_frequency"])

main()
