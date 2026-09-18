"""
Track C, rigidity parameter (a): "for which (m,n) with rank 22 is the lattice
m*U + n*(-E8) even unimodular of signature (3,19)?"

This is the RIGIDITY TEST for the K3 lattice decomposition: (m,n) are
treated as free integer parameters. The scan is generated MECHANICALLY:
every (m,n) with m,n >= 0 and 2m + 8n = 22 (the rank constraint), full stop
-- no candidate is hand-picked, and the true value (m,n)=(3,2) is not
special-cased or excluded from the mechanical enumeration.

For each candidate: assemble the FULL m*U + n*(-E8) Gram matrix, and
COMPUTE (not assume) rank, determinant, evenness, and signature (p,q) by the
same two independent exact methods as item 2 (congruence diagonalization and
Sturm real-root counting with multiplicity). Then check p - q == 0 (mod 8)
computationally on the actually-observed (p,q) -- this is the even-
unimodular-lattice classification constraint (Milnor & Husemoller), verified
here rather than assumed.

"Zero free parameters" for (m,n) is supported only if exactly one candidate
in the mechanically generated solution set satisfies ALL of: rank 22, even,
unimodular, signature (3,19) [equivalently p-q=-16, consistent with p-q=0
mod 8 the moment p+q=22 forces p-q even but not automatically a multiple of
8 -- checked explicitly below], AND at least one nearby (m,n) in the same
scan fails (negative control).
"""
import json
import sys

sys.path.insert(0, "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity/lattices-duality")
from lattice_common import block_diag, U, det_fraction_matrix, signature_congruence, signature_sturm


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
TARGET_RANK = 22

# --- mechanical enumeration of the scan set: every (m,n), m,n>=0, 2m+8n=22 ---
scan_set = []
for n in range(0, TARGET_RANK // 8 + 1):
    remaining = TARGET_RANK - 8 * n
    if remaining < 0 or remaining % 2 != 0:
        continue
    m = remaining // 2
    scan_set.append((m, n))

results_list = []
for (m, n) in scan_set:
    blocks = [U] * m + [NEG_E8] * n
    if blocks:
        M = block_diag(*blocks)
    else:
        M = []
    rank = len(M)
    det_M = det_fraction_matrix(M) if rank > 0 else 1
    is_even = all(M[i][i] % 2 == 0 for i in range(rank)) if rank > 0 else True
    is_unimodular = bool(abs(det_M) == 1)
    if rank > 0:
        pA, qA, zA = signature_congruence(M)
        pB, qB, zB = signature_sturm(M)
    else:
        pA, qA, zA = 0, 0, 0
        pB, qB, zB = 0, 0, 0
    methods_agree = (pA, qA, zA) == (pB, qB, zB)
    p, q = pA, qA
    p_minus_q_mod_8 = (p - q) % 8
    is_signature_3_19 = (p == 3 and q == 19)
    # INTRINSIC conditions: everything derivable from the Gram matrix alone
    # (rank, evenness, unimodularity, no radical, and the p-q=0 mod 8 law for
    # even unimodular lattices). These do NOT reference the target (3,19).
    intrinsic_conditions_hold = (rank == TARGET_RANK and is_even and is_unimodular
                                  and zA == 0 and p_minus_q_mod_8 == 0)
    # EXTERNAL-INPUT condition: matching K3's actual signature (3,19). (3,19)
    # is NOT derived anywhere in this script -- it is the target signature
    # already computed independently in item 2 from K3's Hodge numbers
    # (b0..b4 = 1,0,22,0,1, giving b2=22=h^{1,1}+2, and the Hodge index
    # theorem giving signature (3,19)); that derivation is external to this
    # (m,n) scan and is NOT re-derived here. Treating (3,19) as given is
    # therefore an external-input condition, not an intrinsic one.
    matches_target_signature = is_signature_3_19
    all_conditions_hold = intrinsic_conditions_hold and matches_target_signature
    results_list.append({
        "m": m, "n": n, "rank": rank,
        "determinant": str(det_M),
        "is_even": is_even,
        "is_unimodular": is_unimodular,
        "signature_p_q_zero": [p, q, zA],
        "methods_agree": methods_agree,
        "p_minus_q_mod_8": p_minus_q_mod_8,
        "intrinsic_conditions_hold": bool(intrinsic_conditions_hold),
        "is_signature_3_19": is_signature_3_19,
        "all_consistency_conditions_hold": bool(all_conditions_hold),
    })

intrinsic_solutions = [r for r in results_list if r["intrinsic_conditions_hold"]]
intrinsic_non_solutions = [r for r in results_list if not r["intrinsic_conditions_hold"]]
solutions = [r for r in results_list if r["all_consistency_conditions_hold"]]
non_solutions = [r for r in results_list if not r["all_consistency_conditions_hold"]]

# negative control check: is there at least one nearby (m,n) in the SAME
# rank-22 scan that fails the FULL condition set? (must be true for a
# "zero free parameters" claim of any kind)
negative_control_present = len(non_solutions) >= 1

is_rigid_single_solution = (len(solutions) == 1)
is_rigid_from_intrinsic_conditions_alone = (len(intrinsic_solutions) == 1)

results = {
    "track": "C",
    "item": "rigidity_a_lattice_mn_enum",
    "rigidity_test": {
        "parameter_inserted": "(m, n) = number of U copies, number of (-E8) copies",
        "scanned": f"all (m,n) with m,n>=0 and 2m+8n={TARGET_RANK} (mechanically generated, "
                   f"{len(scan_set)} candidates: {scan_set})",
        "conditions_checked_intrinsic": [
            "rank == 22 (by construction of the scan, always true)",
            "even (all diagonal Gram entries even)",
            "unimodular (|det| == 1)",
            "signature has zero-dimensional radical (no null vectors, z==0)",
            "(p - q) mod 8 == 0 (even-unimodular-lattice constraint, computed not assumed)",
        ],
        "condition_checked_external_input": "signature equals (3,19) -- NOT derived in this script; this is the "
            "K3 signature computed independently in item 2 from K3's given Betti/Hodge numbers (Hodge index "
            "theorem), taken here as an external target, not re-derived from (m,n) alone",
        "solution_set_from_intrinsic_conditions_alone": [(r["m"], r["n"]) for r in intrinsic_solutions],
        "intrinsic_conditions_fix_mn_uniquely": bool(is_rigid_from_intrinsic_conditions_alone),
        "solution_set_after_also_requiring_target_signature": [(r["m"], r["n"]) for r in solutions],
        "negative_control": {
            "present": negative_control_present,
            "failing_candidates_under_full_condition_set": [
                (r["m"], r["n"], r["signature_p_q_zero"], r["p_minus_q_mod_8"]) for r in non_solutions
            ],
        },
        "zero_free_parameters_supported_by_intrinsic_lattice_conditions_alone": bool(
            is_rigid_from_intrinsic_conditions_alone
        ),
        "zero_free_parameters_supported_conditional_on_target_signature_3_19": bool(
            is_rigid_single_solution and negative_control_present
        ),
        "honest_summary": (
            "Rank-22 + even + unimodular + no-radical + (p-q mod 8 == 0) does NOT fix (m,n): "
            f"{len(intrinsic_solutions)} candidates satisfy all of these intrinsic conditions "
            f"({[(r['m'], r['n']) for r in intrinsic_solutions]}), each with a DIFFERENT signature "
            f"({[r['signature_p_q_zero'][:2] for r in intrinsic_solutions]}). Only after ALSO requiring "
            "the target signature (3,19) -- an input taken from K3's known Hodge numbers, not derived "
            "here from (m,n) -- does the solution set collapse to the single point (m,n)=(3,2)."
        ),
    },
    "full_scan_results": results_list,
    "expected": {
        "solution": [3, 2],
        "source": "standard K3 = 3U + 2(-E8) decomposition (task statement / memory), "
                  "unverified prior to this computation -- filled in AFTER the scan above",
    },
}

out_path = "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity/lattices-duality/06_rigidity_a_result.json"
with open(out_path, "w") as f:
    json.dump(results, f, indent=2)

print(json.dumps(results, indent=2))
