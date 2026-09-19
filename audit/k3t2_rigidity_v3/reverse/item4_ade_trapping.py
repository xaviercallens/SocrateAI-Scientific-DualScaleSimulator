"""
Item 4 -- reverse-tests DualScaleDyons.KummerD4.trapping_rank_table beyond its Lean range.

Lean (KummerD4.lean, `decide`) computes, for rank d = 1..8, the largest number of roots of a
semisimple simply-laced (ADE) system of rank <= d: bestTable 8 = [0,2,6,12,24,40,72,126,240], and
notes D4 (24 roots) uniquely beats A4 (20) and every split at rank 4.

This script reimplements the same unbounded-knapsack DP independently (component root-count
formulas A_n = n(n+1), D_n = 2n(n-1) computed structurally; E6/E7/E8 = 72/126/240 declared as a
Tier-L input FROM MEMORY, see inputs.json) and extends it from rank 8 to rank_max_ade_trapping.
It reports, at every rank beyond 8, whether D_n alone achieves the maximum (a RIGID pick) or ties
with some other combination (perturbing the SAME parameter the Lean claim is about: which
component-combination is optimal at a given rank).

IMPORTANT STRUCTURAL CAVEAT (see inputs.json's extended_adeSimple_rank_range): Lean's own
`adeSimple` list only contains A1..A8, D4..D8, E6, E7, E8 -- nothing of rank > 8. To ask the
question "is D_n still optimal beyond rank 8" at all, this script's component list must be
enlarged to A1..A32 and D4..D32 first. So this is not a pure range extension of an existing
Lean `decide` call (which item1/item2/item5 are); it also extends a structural input Lean's
`adeSimple` does not contain. `lean_reachable_now` is therefore False for this item.

Run:
    cd audit/k3t2_rigidity_v3/reverse && \
    /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python item4_ade_trapping.py
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
LEAN_RANGE_MAX = 8
LEAN_BEST_TABLE = [0, 2, 6, 12, 24, 40, 72, 126, 240]


def components(R: int, E: dict) -> list[tuple[int, int, str]]:
    comps = [(n, n * (n + 1), f"A{n}") for n in range(1, R + 1)]
    comps += [(n, 2 * n * (n - 1), f"D{n}") for n in range(4, R + 1)]
    comps += [(6, E["E6"], "E6"), (7, E["E7"], "E7"), (8, E["E8"], "E8")]
    return [c for c in comps if c[1] > 0]


def best_table_with_witness(R: int, E: dict):
    """b[k] = max roots of a semisimple ADE system of rank <= k (unbounded knapsack over
    simple components, exactly mirroring Lean's bestTable). Also returns, for each k, one
    component whose use achieves b[k] (None if b[k] = b[k-1], i.e. padding rank is optimal)."""
    comps = components(R, E)
    b = [0] * (R + 1)
    last_used = [None] * (R + 1)
    for k in range(1, R + 1):
        best = b[k - 1]
        used = None
        for (rank, roots, name) in comps:
            if rank <= k:
                cand = roots + b[k - rank]
                if cand > best:
                    best = cand
                    used = name
        b[k] = best
        last_used[k] = used
    return b, last_used


def d_n_roots(n: int) -> int:
    return 2 * n * (n - 1)


def main() -> None:
    inputs = json.loads((HERE / "inputs.json").read_text())
    R = next(i["value"] for i in inputs if i["name"] == "rank_max_ade_trapping")
    E = next(i["value"] for i in inputs if i["name"] == "E6_E7_E8_root_counts_from_memory")

    b, last_used = best_table_with_witness(R, E)

    regression_ok = b[:9] == LEAN_BEST_TABLE

    beyond = []
    for n in range(LEAN_RANGE_MAX + 1, R + 1):
        dn = d_n_roots(n)
        is_dn_optimal = (dn == b[n])
        # is D_n the UNIQUE way to reach b[n], or does some other combo also reach it (tie)?
        # brute-force: try every component alone plus every pair sum against b[n] restricted
        # to components of rank <= n; report whichever alternative (not "D{n}") also hits b[n].
        comps = components(n, E)
        alt_hits = []
        for (r1, roots1, name1) in comps:
            if name1 == f"D{n}":
                continue
            if r1 == n and roots1 == b[n]:
                alt_hits.append(name1)
            else:
                rem = n - r1
                if rem >= 0:
                    # does some second component (or padding) make up b[n]?
                    if roots1 + best_table_with_witness(rem, E)[0][rem] == b[n] and rem > 0:
                        alt_hits.append(f"{name1}+best({rem})")
        beyond.append({
            "rank": n,
            "best_roots": b[n],
            "D_n_roots": dn,
            "D_n_is_optimal": is_dn_optimal,
            "tie_with_other_combo": len(alt_hits) > 0,
            "alternative_optimal_combos": sorted(set(alt_hits))[:5],
        })

    ties = [row for row in beyond if row["tie_with_other_combo"]]
    strict_d_wins = [row for row in beyond if row["D_n_is_optimal"] and not row["tie_with_other_combo"]]
    d_loses = [row for row in beyond if not row["D_n_is_optimal"]]

    # Negative control: perturb the SAME thing the Lean claim is about (which combination is
    # optimal at a given rank) by checking rank 24 specifically against 3*E8 (a Niemeier-type
    # combination people might naively expect to dominate at rank 24).
    three_e8 = 3 * E["E8"]
    rank24_control = None
    if R >= 24:
        rank24_control = {
            "rank": 24,
            "D24_roots": d_n_roots(24),
            "3xE8_roots": three_e8,
            "best_at_24": b[24],
            "D24_beats_3xE8": d_n_roots(24) > three_e8,
        }

    result = {
        "theorem": "DualScaleDyons.KummerD4.trapping_rank_table",
        "lean_range": [1, LEAN_RANGE_MAX],
        "scanned_range": [1, R],
        "regression_matches_lean_table_ranks_0_to_8": regression_ok,
        "best_table_full": b,
        "beyond_lean_range": beyond,
        "num_ranks_beyond_8": len(beyond),
        "num_ranks_where_Dn_strictly_optimal": len(strict_d_wins),
        "num_ranks_where_Dn_ties_another_combo": len(ties),
        "num_ranks_where_Dn_is_not_optimal": len(d_loses),
        "tie_ranks": [row["rank"] for row in ties],
        "negative_control_rank24_vs_3xE8": rank24_control,
    }
    (HERE / "item4_results.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "beyond_lean_range"}, indent=2))
    print("first few beyond-range rows:", json.dumps(beyond[:6], indent=2))


if __name__ == "__main__":
    main()
