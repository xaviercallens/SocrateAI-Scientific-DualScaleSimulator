"""
Track C, item (2): K3 lattice = 3U + 2(-E8); Mukai lattice = K3 + U;
Gamma^{6,22} = 6U + 2(-E8).

For each: assemble the FULL explicit Gram matrix (block diagonal, built from
the same E8 Cartan matrix as item 1, and the hyperbolic plane U), then
compute rank, determinant, evenness, unimodularity, and signature (p,q)
by TWO independent exact methods on the assembled matrix itself (not by
assuming signature is additive over blocks):
  (A) symmetric congruence diagonalization over Q (exact Fractions), and
  (B) exact characteristic-polynomial sign counting (sympy Sturm-based
      count_roots on QQ), which is genuinely independent machinery.
Agreement of (A) and (B) is the "two different, same finding" check the
task asks for at the arithmetic level.
"""
import json
import sys

sys.path.insert(0, "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity/lattices-duality")
from lattice_common import (
    block_diag, U, det_fraction_matrix, signature_congruence, signature_sturm,
)

# Rebuild E8 Cartan matrix exactly as in 01_e8.py (self-contained; not imported
# as a pre-typed constant from that script's output -- rebuilt from the same
# Dynkin-diagram adjacency definition).
def build_e8_cartan():
    n = 8
    edges = [(1, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (2, 4)]
    C = [[0] * n for _ in range(n)]
    for i in range(n):
        C[i][i] = 2
    for (a, b) in edges:
        i, j = a - 1, b - 1
        C[i][j] = -1
        C[j][i] = -1
    return C

E8 = build_e8_cartan()
NEG_E8 = [[-x for x in row] for row in E8]


def analyze(name, M, expected_note):
    n = len(M)
    det_M = det_fraction_matrix(M)
    is_even = all(M[i][i] % 2 == 0 for i in range(n))
    p1, n1, z1 = signature_congruence(M)
    p2, n2, z2 = signature_sturm(M)
    methods_agree = (p1, n1, z1) == (p2, n2, z2)
    is_unimodular = bool(abs(det_M) == 1)
    return {
        "name": name,
        "rank": n,
        "determinant": str(det_M),
        "is_even": is_even,
        "is_unimodular": is_unimodular,
        "signature_method_A_congruence_diagonalization": {"p": p1, "q": n1, "zero": z1},
        "signature_method_B_sturm_charpoly": {"p": p2, "q": n2, "zero": z2},
        "methods_agree": methods_agree,
        "expected": expected_note,
    }


results = {}

# K3 lattice = 3U + 2(-E8): rank 3*2 + 2*8 = 22
K3 = block_diag(U, U, U, NEG_E8, NEG_E8)
results["K3_3U_2negE8"] = analyze(
    "K3 = 3U + 2(-E8)", K3,
    {"signature": [3, 19], "det": 1, "even": True,
     "source": "standard K3 lattice fact, from memory, unverified prior to computation"},
)

# Mukai lattice = K3 + U : rank 24
MUKAI = block_diag(U, U, U, U, NEG_E8, NEG_E8)
results["Mukai_K3_plus_U"] = analyze(
    "Mukai = K3 + U", MUKAI,
    {"signature": [4, 20], "det": 1, "even": True,
     "source": "standard Mukai lattice fact, from memory, unverified prior to computation"},
)

# Gamma^{6,22} = 6U + 2(-E8): rank 6*2 + 2*8 = 28
GAMMA_6_22 = block_diag(U, U, U, U, U, U, NEG_E8, NEG_E8)
results["Gamma_6_22"] = analyze(
    "Gamma^{6,22} = 6U + 2(-E8)", GAMMA_6_22,
    {"signature": [6, 22], "det": 1, "even": True,
     "source": ("task statement / standard K3xT2 heterotic-type-II charge lattice fact, "
                "from memory, unverified prior to computation")},
)

out = {
    "track": "C",
    "item": "2_k3_mukai_gamma622",
    "method": "Explicit block-diagonal Gram matrices from U=[[0,1],[1,0]] and the E8 Cartan "
               "matrix (rebuilt from Dynkin adjacency); rank/det/evenness read off the "
               "assembled matrix; signature computed by two independent exact methods "
               "(congruence diagonalization over Q, and Sturm-sequence sign counting on "
               "the exact characteristic polynomial) on the FULL matrix, not assumed additive.",
    "results": results,
}

out_path = "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity/lattices-duality/02_k3_mukai_gamma_result.json"
with open(out_path, "w") as f:
    json.dump(out, f, indent=2)

print(json.dumps(out, indent=2))
