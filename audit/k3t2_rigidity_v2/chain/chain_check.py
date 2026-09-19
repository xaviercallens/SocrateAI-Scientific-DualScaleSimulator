#!/usr/bin/env python3
"""
Chain-check for the K3 x T2 rigidity v2 tracks.

Verifies numerically, from the tracks' own exports/results JSON files (no
literal values typed in this script), whether the claimed chain

  D (GUDHI TDA)  -> chi, b2
  D  -> C (Noether + Hodge index)      -> signature -> lattice 3U + 2(-E8)
  A (elliptic genus)                    -> Z(tau,0) must equal D's chi
  B (DMVV)                              -> Euler numbers of Hilb^k must follow
                                            from A's c(0)+2c(-1) and D's chi

is actually numerically consistent, i.e. that the different tracks'
independently-computed numbers agree with each other, and that D's exported
chi is the SAME NUMBER that C's script actually consumed (not just a file
that happens to sit at the expected path).

Every number used below is READ from a file under this same
audit/k3t2_rigidity_v2/ directory; nothing is hardcoded except structural
integers that define the U and E8 lattices themselves (rank(U)=2,
signature(U)=(1,1); rank(E8)=8, signature(-E8)=(0,8)), which are lattice-
theory definitions, not track outputs.

Run with:
  /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python \
      /mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity_v2/chain/chain_check.py
"""
import json
import os
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)  # audit/k3t2_rigidity_v2/


def load(*relpath):
    p = os.path.join(ROOT, *relpath)
    with open(p) as f:
        return json.load(f), p


def as_int(x):
    return int(Fraction(str(x)))


def main():
    links = []
    definitional = []

    # ---- load raw exports (no literals) ----
    A_exp, A_path = load("A-genus", "exports.json")
    B_res, B_path = load("B-dyons", "results.json")
    B_part1, B_part1_path = load("B-dyons", "part1_euler_numbers_results.json")
    C_exp, C_path = load("C-lattices", "exports.json")
    D_exp, D_path = load("D-tda", "exports.json")

    chi_A = as_int(A_exp["chi_from_genus"])
    chi_D = as_int(D_exp["chi_K3_resolved"])
    b2_D = as_int(D_exp["b2_K3_resolved"])
    b1_D = as_int(D_exp["b1"])

    chi_B_keyresult = as_int(B_res["key_results"]["chi_K3_computed"])
    c0_B = as_int(B_part1["c0"])
    cm1_B = as_int(B_part1["cm1"])
    chi_B_from_c0_cm1 = c0_B + 2 * cm1_B
    chi_B_reported = as_int(B_part1["chi_computed_c0_plus_2cm1"])

    chi_top_C = as_int(C_exp["chi_top_used"])
    chi_top_source_file_C = C_exp["chi_top_source_file"]
    p_C, q_C = C_exp["derived_signature_p_q"]
    p_C, q_C = as_int(p_C), as_int(q_C)
    m_sel, n_sel = C_exp["selected_mn"][0]
    m_sel, n_sel = as_int(m_sel), as_int(n_sel)

    # ================================================================
    # Link 1: Track A's independently computed elliptic-genus chi must
    # equal Track D's independently computed (GUDHI) chi. Neither
    # value was produced FROM the other -- this is a genuine
    # cross-track agreement check, not a tautology.
    # ================================================================
    ok1 = (chi_A == chi_D)
    links.append({
        "link": "A.chi_from_genus == D.chi_K3_resolved",
        "A_value": chi_A, "A_source": A_path,
        "D_value": chi_D, "D_source": D_path,
        "consistent": ok1,
    })

    # ================================================================
    # Link 2: Track B's own theta-series computation of c(0)+2c(-1)
    # (part1) must equal the chi_K3_computed it reports in results.json
    # (internal self-consistency of Track B), and both must equal
    # Track D's independently computed chi (cross-track agreement).
    # ================================================================
    ok2a = (chi_B_from_c0_cm1 == chi_B_reported)
    ok2b = (chi_B_reported == chi_B_keyresult)
    ok2c = (chi_B_keyresult == chi_D)
    links.append({
        "link": "B.c0 + 2*B.cm1 == B.chi_computed_c0_plus_2cm1 (internal)",
        "c0": c0_B, "cm1": cm1_B, "computed": chi_B_from_c0_cm1,
        "reported": chi_B_reported,
        "source": B_part1_path,
        "consistent": ok2a,
    })
    links.append({
        "link": "B.part1.chi_computed == B.results.key_results.chi_K3_computed (internal)",
        "part1_value": chi_B_reported, "part1_source": B_part1_path,
        "results_value": chi_B_keyresult, "results_source": B_path,
        "consistent": ok2b,
    })
    links.append({
        "link": "B.chi_K3_computed == D.chi_K3_resolved (cross-track)",
        "B_value": chi_B_keyresult, "B_source": B_path,
        "D_value": chi_D, "D_source": D_path,
        "consistent": ok2c,
    })

    # ================================================================
    # Link 3: Track C claims it read chi_top from D-tda/exports.json.
    # Verify (a) the file it names is actually the file D wrote, and
    # (b) the number C used (chi_top_used) equals what is actually
    # sitting in that file right now, i.e. the provenance pointer is
    # not stale / was not silently overridden.
    # ================================================================
    named_path_matches = (chi_top_source_file_C == "D-tda/exports.json")
    ok3 = named_path_matches and (chi_top_C == chi_D)
    links.append({
        "link": "C.chi_top_used == D.chi_K3_resolved (provenance pointer verified)",
        "C_chi_top_used": chi_top_C,
        "C_declared_source_file": chi_top_source_file_C,
        "D_actual_value": chi_D, "D_actual_source": D_path,
        "named_path_matches_actual_D_file": named_path_matches,
        "consistent": ok3,
    })

    # ================================================================
    # Link 4: rank check. C derives signature (p,q) from
    # chi_top (via Noether's formula) and b2. p+q must equal D's own
    # b2 (rank of H^2), independently exported by D.
    # ================================================================
    ok4 = (p_C + q_C == b2_D)
    links.append({
        "link": "C.(p+q) == D.b2_K3_resolved (rank agreement)",
        "p": p_C, "q": q_C, "p_plus_q": p_C + q_C,
        "D_b2": b2_D, "D_source": D_path,
        "consistent": ok4,
    })

    # ================================================================
    # Link 5: the selected lattice m*U + n*(-E8) must actually
    # reproduce the derived signature (p_C, q_C) that Track C computed,
    # using only the DEFINITIONAL signatures of U=(1,1) and
    # -E8=(0,8) (these two pairs are lattice-theory definitions, not
    # a track's fitted output -- flagged as definitional below).
    # ================================================================
    U_sig = (1, 1)       # definitional: hyperbolic plane U has signature (1,1)
    negE8_sig = (0, 8)   # definitional: -E8 (negative-definite E8) has signature (0,8)
    definitional.append({"quantity": "signature(U)", "value": U_sig,
                          "reason": "definition of the hyperbolic lattice U"})
    definitional.append({"quantity": "signature(-E8)", "value": negE8_sig,
                          "reason": "definition of the negative-definite E8 root lattice"})

    p_from_mn = m_sel * U_sig[0] + n_sel * negE8_sig[0]
    q_from_mn = m_sel * U_sig[1] + n_sel * negE8_sig[1]
    rank_from_mn = m_sel * 2 + n_sel * 8  # rank(U)=2, rank(E8)=8 -- also definitional
    definitional.append({"quantity": "rank(U), rank(E8)", "value": [2, 8],
                          "reason": "standard rank of the hyperbolic plane / E8 root lattice"})

    ok5a = (p_from_mn == p_C and q_from_mn == q_C)
    ok5b = (rank_from_mn == b2_D)
    links.append({
        "link": "selected (m,n)=(%d,%d): m*U + n*(-E8) signature == C.derived_signature_p_q" % (m_sel, n_sel),
        "computed_signature": [p_from_mn, q_from_mn],
        "C_derived_signature": [p_C, q_C],
        "consistent": ok5a,
    })
    links.append({
        "link": "selected (m,n) rank == D.b2_K3_resolved",
        "computed_rank": rank_from_mn, "D_b2": b2_D,
        "consistent": ok5b,
    })

    # ================================================================
    # Link 6: even-unimodular necessary condition p - q == 0 mod 8
    # (structural, definitional arithmetic fact about even unimodular
    # lattices -- not a track output, so also flagged as a
    # definitional input used to accept the chain, not a free result).
    # ================================================================
    ok6 = ((p_C - q_C) % 8 == 0)
    definitional.append({"quantity": "even-unimodular signature rule (p-q) mod 8 == 0",
                          "reason": "classical lattice-theory fact (Milnor-Husemoller), used as an acceptance criterion, not computed by any track"})
    links.append({
        "link": "(C.p - C.q) mod 8 == 0 (even-unimodular necessary condition)",
        "p_minus_q": p_C - q_C, "mod8": (p_C - q_C) % 8,
        "consistent": ok6,
    })

    # ================================================================
    # Additional definitional inputs already flagged upstream by the
    # tracks themselves (surfaced here, not re-derived, so a reader
    # sees them alongside the numeric links).
    # ================================================================
    definitional.append({"quantity": "c1(K3) = 0", "reason": "Calabi-Yau condition, assumed by Track C, not computed",
                          "source": C_path})
    definitional.append({"quantity": "b1(K3) = 0", "reason": "simply-connectedness, read from D-tda/exports.json (D.b1) and assumed by Track C",
                          "source": D_path, "D_b1_value": b1_D})
    definitional.append({"quantity": "chi_top = 24 used by Track C's Kahler formula b2+ = 2*h^{2,0}+1",
                          "reason": "the formula itself (not just its numeric input chi_top) is a standard Hodge-theory identity assumed, not derived, by Track C",
                          "source": C_path})

    all_consistent = all(l["consistent"] for l in links)

    out = {
        "links": links,
        "all_consistent": all_consistent,
        "definitional_inputs_not_computed_by_any_track": definitional,
        "files_read": [A_path, B_path, B_part1_path, C_path, D_path],
    }

    out_path = os.path.join(HERE, "chain_check_results.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)

    print(json.dumps(out, indent=2))
    print("\n=== SUMMARY ===")
    for l in links:
        status = "CONSISTENT" if l["consistent"] else "INCONSISTENT"
        print(f"[{status}] {l['link']}")
    print(f"\nALL CONSISTENT: {all_consistent}")
    print(f"Results written to: {out_path}")


if __name__ == "__main__":
    main()
