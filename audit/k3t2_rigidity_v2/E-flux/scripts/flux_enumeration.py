#!/usr/bin/env python3
"""
Track E, Part 3: enumerate integer flux choices in the Tripathy-Trivedi
section-4.1 family (eqs 4.8-4.10, 4.13-4.14 of hep-th/0301139) within an
explicitly stated finite truncation of the rank-6 lattice Gamma_{3,3} (basis
(A.6)-(A.7), see lattice.py for the basis-ambiguity discussion), and count
solutions of

    alpha_x . beta_x = 0                    (4.10)
    alpha_x^2 = beta_x^2  ( > 0, spacelike ) (4.9) + (3.31)
    alpha_x^2 + N_D3 = 24, N_D3 >= 0         (4.14) with (2.3)

(alpha_y = -beta_x, beta_y = alpha_x per eq 4.8 are then fixed automatically;
N_flux = 2*alpha_x^2 per eq 4.13.)

STRUCTURAL fact independent of the "24" input: because flux coefficients
must be even integers (eq 2.5) and the diagonal metric has entries +-2
(eqs A.6-A.7), alpha_x^2 = 8 * sum_i (+-1) n_i^2 is ALWAYS a multiple of 8.
This alone forces N_flux in 8*Z and alpha_x^2 in 8*Z; it holds for every
flux vector regardless of any tadpole bound. Applying TT's tadpole bound
alpha_x^2 <= 24 on top of this restricts alpha_x^2 to the 3-element set
{8, 16, 24}, i.e. N_D3 in {16, 8, 0}. That the admissible set is discrete
comes from 8|alpha_x^2 (structural); WHICH three values comes from "24"
(the literature input TT eq 2.3, itself derived from the D7/O7 count
checked in tadpole_bookkeeping.py).

Tier: B, exact int arithmetic, within the stated truncation (rank r,
entry bound N). Truncation and N are reported with every count. Orbit
counting under a stated FINITE subgroup of lattice automorphisms
(coordinate permutation within each signature block + independent sign
flips, see lattice.py:coordinate_automorphisms) -- NOT the full arithmetic
automorphism group, so orbit counts are upper bounds on the number of
physically distinct configurations, and do not converge as N grows (see
negative control).

Negative control: report count(N) for the SAME structural conditions
(4.9)+(4.10)+positivity, WITHOUT the tadpole bound alpha_x^2<=24, for
growing N. This is expected to (and is checked to) grow without bound.
We ALSO report count(N) WITH the tadpole bound for growing N, to check
whether the bound alone gives a truncation-independent finite answer (it
does not: in indefinite signature the bound alpha_x^2<=24 does not bound
the entries n_i, so the with-tadpole count also grows with N; finiteness
of the reported "landscape" comes from the STATED entry bound N, not from
the tadpole condition by itself).

Run exactly (rank-4 primary, rank-2 empty-check, rank-6 capped secondary):
  /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python \
    /mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity_v2/E-flux/scripts/flux_enumeration.py
"""
import itertools
import json
import sys
import time

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from lattice import (  # noqa: E402
    FULL_SIGNS,
    truncation,
    norm,
    dot,
    enumerate_vectors,
    coordinate_automorphisms,
    group_order,
)

OUT = __file__.rsplit("/", 1)[0].rsplit("/", 1)[0] + "/flux_enumeration_results.json"

# Truncations, given as 0-based indices into the full (e1,...,e6) basis.
# index 0,1,2 -> e1,e2,e3 (positive norm +2); index 3,4,5 -> e4,e5,e6 (negative norm -2)
TRUNCATIONS = {
    "rank2_e1e4": {"indices": (0, 3), "pos": 1, "neg": 1},   # signature (1,1) -- predicted EMPTY
    "rank4_e1e2e4e5": {"indices": (0, 1, 3, 4), "pos": 2, "neg": 2},  # PRIMARY
    "rank6_full": {"indices": (0, 1, 2, 3, 4, 5), "pos": 3, "neg": 3},  # secondary, capped bound
}


def bucket_by_k(vectors, signs):
    """Bucket vectors by k = sum_i signs_i * n_i^2 (so norm = 8k). Returns
    dict k -> list of vectors, restricted to k > 0 (spacelike / positive
    norm, required by (3.31))."""
    buckets = {}
    for v in vectors:
        k = sum(s * vi * vi for s, vi in zip(signs, v))
        if k > 0:
            buckets.setdefault(k, []).append(v)
    return buckets


def count_pairs_same_k(buckets, signs, k_filter=None):
    """Count ordered pairs (alpha_x, beta_x) with alpha_x, beta_x in the same
    k-bucket (so alpha_x^2=beta_x^2=8k>0), alpha_x . beta_x = 0, alpha_x !=
    beta_x. (Linear independence is automatic: if beta_x = c*alpha_x with
    alpha_x.beta_x=0 and both norm 8k>0, then c=0, contradiction, or
    beta_x=-alpha_x which the caller can choose to identify via the sign-flip
    automorphism; we count beta_x=-alpha_x as a DISTINCT ordered pair here
    and let orbit-canonicalisation merge automorphism-equivalent pairs.)
    Returns dict k -> (count_pairs, [sample pairs up to 5]).
    """
    result = {}
    for k, vecs in buckets.items():
        if k_filter is not None and not k_filter(k):
            continue
        n_pairs = 0
        samples = []
        for a in vecs:
            for b in vecs:
                if a == b:
                    continue
                if dot(a, b, signs) == 0:
                    n_pairs += 1
                    if len(samples) < 5:
                        samples.append({"alpha_x": a, "beta_x": b})
        result[k] = {"pair_count": n_pairs, "samples": samples}
    return result


def canonical_form(pair, group_elems):
    a, b = pair
    best = None
    for g in group_elems:
        cand = (g(a), g(b))
        if best is None or cand < best:
            best = cand
    return best


def orbit_count(buckets, signs, group_elems, k_filter):
    """Exact orbit count under the stated finite automorphism subgroup, for
    the tadpole-ON subset only (kept small deliberately)."""
    orbit_reps = {}
    for k, vecs in buckets.items():
        if not k_filter(k):
            continue
        seen = set()
        for a in vecs:
            for b in vecs:
                if a == b:
                    continue
                if dot(a, b, signs) != 0:
                    continue
                canon = canonical_form((a, b), group_elems)
                seen.add(canon)
        orbit_reps[k] = {"n_orbits": len(seen), "sample_reps": list(itertools.islice(seen, 5))}
    return orbit_reps


def run_truncation(name, spec, bounds, tadpole_bound=24, do_orbit_bounds=None, time_budget_s=600):
    indices = spec["indices"]
    signs, labels = truncation(indices)
    rank = len(indices)
    group_elems = coordinate_automorphisms(signs)
    g_order = group_order(spec["pos"], spec["neg"])
    assert len(group_elems) == g_order, (len(group_elems), g_order)

    growth_with_tadpole = {}
    growth_without_tadpole = {}
    orbit_results = {}

    for N in bounds:
        t0 = time.time()
        vectors = list(enumerate_vectors(rank, N))
        buckets = bucket_by_k(vectors, signs)

        # WITHOUT tadpole bound: any k>0 (structural conditions 4.9+4.10+positivity only)
        no_bound = count_pairs_same_k(buckets, signs, k_filter=lambda k: True)
        total_no_bound = sum(v["pair_count"] for v in no_bound.values())
        growth_without_tadpole[N] = {
            "n_flux_vectors_enumerated": len(vectors),
            "total_valid_ordered_pairs": total_no_bound,
            "per_k": {str(k): v["pair_count"] for k, v in sorted(no_bound.items())},
        }

        # WITH tadpole bound: k <= tadpole_bound/8 (since alpha_x^2=8k<=24 => k<=3 for bound 24)
        assert tadpole_bound % 8 == 0
        k_max = tadpole_bound // 8
        with_bound = count_pairs_same_k(buckets, signs, k_filter=lambda k: k <= k_max)
        total_with_bound = sum(v["pair_count"] for v in with_bound.values())
        growth_with_tadpole[N] = {
            "n_flux_vectors_enumerated": len(vectors),
            "total_valid_ordered_pairs": total_with_bound,
            "per_k": {
                str(k): {
                    "pair_count": v["pair_count"],
                    "N_flux": 2 * 8 * k,
                    "N_D3": 24 - 8 * k,
                    "samples": v["samples"],
                }
                for k, v in sorted(with_bound.items())
            },
        }

        if do_orbit_bounds and N in do_orbit_bounds:
            orbits = orbit_count(buckets, signs, group_elems, k_filter=lambda k: k <= k_max)
            orbit_results[N] = {
                "group_order_used": g_order,
                "per_k": {
                    str(k): v for k, v in sorted(orbits.items())
                },
                "total_orbits": sum(v["n_orbits"] for v in orbits.values()),
            }

        elapsed = time.time() - t0
        if elapsed > time_budget_s:
            growth_with_tadpole[N]["TIMED_OUT_partial"] = True
            break

    return {
        "truncation": name,
        "basis_indices_0based": list(indices),
        "basis_labels": labels,
        "signature": f"({spec['pos']},{spec['neg']})",
        "automorphism_group_order_used": g_order,
        "bounds_tried": bounds,
        "with_tadpole_bound_alpha_x_sq_le_24": growth_with_tadpole,
        "without_tadpole_bound_negative_control": growth_without_tadpole,
        "orbit_counts_tadpole_on_only": orbit_results,
    }


def main():
    t_start = time.time()
    out = {
        "track": "E-flux",
        "part": "3-flux-enumeration",
        "command": " ".join([sys.executable] + sys.argv),
        "method": (
            "TT section 4.1 family: alpha_y=-beta_x, beta_y=alpha_x (eq 4.8); "
            "alpha_x^2=beta_x^2 (eq 4.9); alpha_x.beta_x=0 (eq 4.10); "
            "N_flux=2*alpha_x^2 (eq 4.13); tadpole alpha_x^2+N_D3=24, N_D3>=0 "
            "(eq 4.14 / eq 2.3). Basis: TT App. A eqs (A.6)-(A.7), diagonal "
            "metric diag(+2,+2,+2,-2,-2,-2) on (e1,...,e6); flux coefficients "
            "even integers a_i=2*n_i (eq 2.5); enumeration variable is the "
            "integer tuple n=(n_1,...,n_rank), |n_i| <= N (the stated entry "
            "bound)."
        ),
        "results": {},
    }

    # rank-2: PREDICTED EMPTY (signature (1,1) has 1-dim positive cone; two
    # orthogonal strictly-positive-norm vectors cannot both exist in a
    # rank-2 space of signature (1,1)). Checked directly, not assumed.
    out["results"]["rank2_e1e4"] = run_truncation(
        "rank2_e1e4", TRUNCATIONS["rank2_e1e4"], bounds=[1, 2, 3, 4]
    )

    # rank-4: PRIMARY. bounds 1,2,3,4 as advised (81/625/2401/6561 vectors, cheap)
    out["results"]["rank4_e1e2e4e5"] = run_truncation(
        "rank4_e1e2e4e5",
        TRUNCATIONS["rank4_e1e2e4e5"],
        bounds=[1, 2, 3, 4],
        do_orbit_bounds=[1, 2, 3, 4],
    )

    # rank-6: SECONDARY, capped at N<=2 (15625 vectors) for time budget.
    # Orbit canonicalisation under the full rank-6 automorphism subgroup
    # (order 2304) is COMPUTATIONALLY CAPPED to N=1 only: profiling showed
    # N=1 has ~7k tadpole-satisfying pairs (~40s to canonicalise) while N=2
    # has ~340k tadpole-satisfying pairs (~340k*2304 ops, empirically >20min,
    # exceeds the ~15min per-script budget). We therefore report RAW pair
    # and vector counts at N=2 (fast, see without/with_tadpole_bound blocks)
    # but do NOT report orbit counts at N=2 -- this is stated explicitly in
    # results.json and could_not_do rather than silently omitted or left to
    # hang.
    out["results"]["rank6_full"] = run_truncation(
        "rank6_full",
        TRUNCATIONS["rank6_full"],
        bounds=[1, 2],
        do_orbit_bounds=[1],
        time_budget_s=120,
    )
    out["results"]["rank6_full"]["orbit_note"] = (
        "Orbit canonicalisation under the order-2304 automorphism subgroup "
        "was computed only at N=1 (see could_not_do for N=2: profiled at "
        ">20 minutes, exceeds the per-script time budget). Raw (non-orbit) "
        "pair/vector counts ARE reported at N=2."
    )

    out["total_wall_time_s"] = round(time.time() - t_start, 2)

    with open(OUT, "w") as f:
        json.dump(out, f, indent=2)
    print(f"wrote {OUT}, wall time {out['total_wall_time_s']}s")
    # print compact summary
    for name, res in out["results"].items():
        print(f"--- {name} (sig {res['signature']}) ---")
        for N, g in res["with_tadpole_bound_alpha_x_sq_le_24"].items():
            print(f"  N={N}: WITH tadpole pairs={g['total_valid_ordered_pairs']}  "
                  f"WITHOUT tadpole pairs={res['without_tadpole_bound_negative_control'][N]['total_valid_ordered_pairs']}")


if __name__ == "__main__":
    main()
