"""
Track C (v2), item (1): E8 Cartan matrix -- determinant, evenness, positive
definiteness (exact leading principal minors), and root count (exact
enumeration of norm-2 vectors), built entirely from the Dynkin diagram
adjacency (no imported/typed-in 8x8 grid, no R^8 embedding assumed).

Reused with only path changes from v1 (audit/k3t2_rigidity/lattices-duality/
01_e8.py, written blind); logic unchanged since it already had no defect
flagged for this item.

Method for root count: the Cartan matrix C IS the Gram matrix of the simple
roots (alpha_i, alpha_j) in the simply-laced normalization (alpha_i,alpha_i)=2.
A root beta = sum_j c_j alpha_j has coefficient vector c in Z^8. The simple
reflection s_i acts on coefficient vectors by
    c -> c - (C c)_i * e_i
(derived from <beta, alpha_i^vee> = (C c)_i for simply-laced C). Starting
from the 8 simple roots (c = e_i) and closing under all 8 simple reflections
(BFS) generates the *entire* root system exactly, in integer arithmetic,
directly from C -- no separate R^8 realization is used or assumed.
Each generated root is checked to satisfy c^T C c == 2 as an independent
sanity check.

Run: /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python \
     /mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity_v2/C-lattices/01_e8.py
"""
import json
import sys
from collections import deque

HERE = "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity_v2/C-lattices"
sys.path.insert(0, HERE)
from lattice_common import leading_principal_minors, det_fraction_matrix

N = 8

# Bourbaki-labeled E8 Dynkin diagram adjacency (built from the diagram, not
# copied as a pre-made matrix):
#   1 - 3 - 4 - 5 - 6 - 7 - 8
#           |
#           2
EDGES = [(1, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (2, 4)]

def build_cartan(n, edges):
    C = [[0] * n for _ in range(n)]
    for i in range(n):
        C[i][i] = 2
    for (a, b) in edges:
        i, j = a - 1, b - 1
        C[i][j] = -1
        C[j][i] = -1
    return C

C = build_cartan(N, EDGES)

# --- determinant, evenness, leading principal minors (positive definiteness) ---
det_C = det_fraction_matrix(C)
is_even = all(C[i][i] % 2 == 0 for i in range(N))
minors = leading_principal_minors(C)
is_pos_def = all(m > 0 for m in minors)

# --- exact root enumeration by Weyl-reflection closure on coefficient vectors ---
def apply_C(c):
    return [sum(C[i][j] * c[j] for j in range(N)) for i in range(N)]

def norm2(c):
    Cc = apply_C(c)
    return sum(c[i] * Cc[i] for i in range(N))

start_roots = [tuple(1 if k == i else 0 for k in range(N)) for i in range(N)]
seen = set(start_roots)
queue = deque(start_roots)
bad = []
while queue:
    c = queue.popleft()
    Cc = apply_C(list(c))
    if norm2(list(c)) != 2:
        bad.append(c)
    for i in range(N):
        coeff = Cc[i]
        if coeff == 0:
            continue
        new_c = list(c)
        new_c[i] = new_c[i] - coeff
        new_c = tuple(new_c)
        if new_c not in seen:
            seen.add(new_c)
            queue.append(new_c)

num_roots = len(seen)
all_norm2 = all(norm2(list(c)) == 2 for c in seen)

# negative control: verify a coefficient vector that is NOT in the closure
# (e.g. 2*e_1, norm^2 = 2^2*2 = 8 != 2) correctly fails the norm-2 test and
# is (correctly) absent from the generated root set.
neg_control_vec = tuple(2 if k == 0 else 0 for k in range(N))
neg_control_norm2 = norm2(list(neg_control_vec))
neg_control_in_roots = neg_control_vec in seen
neg_control_passes = (neg_control_norm2 != 2) and (not neg_control_in_roots)

result = {
    "track": "C", "version": "v2",
    "item": "1_e8_cartan",
    "method": "Cartan matrix built from Dynkin-diagram adjacency (Bourbaki labeling); "
              "det/minors via sympy exact rationals; roots via exact Weyl-reflection "
              "closure c -> c - (Cc)_i e_i on integer coefficient vectors (no R^8 embedding).",
    "cartan_matrix": C,
    "determinant": {"value": str(det_C), "exact": True},
    "is_even": is_even,
    "leading_principal_minors": [str(m) for m in minors],
    "is_positive_definite": is_pos_def,
    "root_enumeration": {
        "num_roots_found": num_roots,
        "all_have_norm2_exactly_2": all_norm2,
        "num_bad_norm2_during_bfs": len(bad),
    },
    "negative_control": {
        "description": "coefficient vector 2*e_1: has norm^2 = 8 != 2, and was checked to be "
                        "absent from the reflection-closure set generated above",
        "vector": list(neg_control_vec),
        "norm2": neg_control_norm2,
        "in_generated_root_set": neg_control_in_roots,
        "control_passes": neg_control_passes,
    },
    "expected": {
        "determinant": 1,
        "num_roots": 240,
        "source": "standard E8 root-system facts (determinant of E8 Cartan matrix = 1, "
                   "240 roots), from memory, unverified prior to this computation -- "
                   "filled in AFTER computing the fields above",
    },
}

out_path = HERE + "/01_e8_result.json"
with open(out_path, "w") as f:
    json.dump(result, f, indent=2)

print(json.dumps(result, indent=2))
