"""Independent reproduction check against the earlier genetics validation.

The earlier run (a different worktree, different code, different author session)
published Rips cross-check numbers for Caulobacter, E. coli and the influenza
segments.  This block recomputed them from the same raw files with independently
written code.  This script compares the two, digit by digit, and records the result
as a control in TopoDB.

It asserts nothing about biology; it checks that two independent implementations of
"Rips on the given metric" agree.  A disagreement would be reported as such.

Usage: python check_reproduction.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, "/mnt/disks/disk-socrateai-local-1/wt-topo-bio/topodb_runs/biology")
from bio_common import RESULTS, Block  # noqa: E402

REPORT = Path("/mnt/disks/disk-socrateai-local-1/dualscale-wt-tdaval/audit/tda_validation/"
              "genetics/report.json")
SCRIPT = "topodb_runs/biology/check_reproduction.py"
COMMAND = ("timeout 300 prlimit --as=8589934592 -- "
           "/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python "
           "topodb_runs/biology/check_reproduction.py")
TOL = 1e-4          # relative agreement required to call a value reproduced


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-12)


def main() -> int:
    ref = json.loads(REPORT.read_text())
    comps = []

    # ---- Hi-C, full resolution (the coarsened path is a different statistic) ----
    hic_map = {
        "caulo_full": ("test_2_Caulobacter_HiC", "full_ring"),
        "caulo_cut_ter": ("test_2_Caulobacter_HiC", "ter_cut_control"),
        "caulo_cut_ori": ("test_2_Caulobacter_HiC", "ori_cut_control"),
        "gm12878_chr1q": ("test_2_Caulobacter_HiC", "external_linear_control_GM12878_chr1q"),
        "ecoli_full": ("test_2b_Ecoli_3Cseq_addendum", "full_ring"),
        "ecoli_cut_ter": ("test_2b_Ecoli_3Cseq_addendum", "ter_cut_control"),
    }
    for case, (test, sub) in hic_map.items():
        f = RESULTS / f"hic_{case}.json"
        if not f.exists():
            comps.append({"quantity": f"hic {case} rips P1/P2", "status": "ABSENT (not computed)"})
            continue
        mine = json.loads(f.read_text())["meta"]["full_resolution_observed"]["h1_dominance_P1_over_P2"]
        theirs = ref["tests"][test][sub]["rips_P1_over_P2"]
        comps.append({"quantity": f"hic {case} Rips P1/P2 (full resolution)",
                      "mine": round(mine, 8), "earlier_validation": theirs,
                      "relative_difference": round(rel(mine, theirs), 10),
                      "agrees_to_1e-4": bool(rel(mine, theirs) < TOL)})

    # ---- influenza ----
    f = RESULTS / "influenza_summary.json"
    if f.exists():
        mine = json.loads(f.read_text())["observed"]
        t3 = ref["tests"]["test_3_influenza_reassortment"]
        for s, nm in ((1, "PB2"), (2, "PB1"), (3, "PA"), (4, "HA"), (5, "NP"), (6, "NA"),
                      (7, "M"), (8, "NS")):
            m = mine[f"seg{s}_{nm}"]
            t = t3["single_segments"][f"{s}_{nm}"]
            for key, theirs_key in (("A_max_h1_persistence", "A_max_h1_pers"),
                                    ("B_n_h1_bars_ge_0.005", "B_n_bars_ge_0.005"),
                                    ("n_h1_bars", "n_h1")):
                comps.append({"quantity": f"influenza seg{s} {nm} {theirs_key}",
                              "mine": round(m[key], 8), "earlier_validation": t[theirs_key],
                              "relative_difference": round(rel(m[key], t[theirs_key]), 10)
                              if t[theirs_key] else (0.0 if m[key] == 0 else 1.0),
                              "agrees_to_1e-4": bool(
                                  abs(m[key] - t[theirs_key]) <= max(TOL * abs(t[theirs_key]),
                                                                     5e-5))})
        mc, tc = mine["concatenated"], t3["concatenated"]
        for key, theirs_key in (("A_max_h1_persistence", "A_max_h1_pers"),
                                ("B_n_h1_bars_ge_0.005", "B_n_bars_ge_0.005"),
                                ("n_h1_bars", "n_h1")):
            comps.append({"quantity": f"influenza concatenated {theirs_key}",
                          "mine": round(mc[key], 8), "earlier_validation": tc[theirs_key],
                          "relative_difference": round(rel(mc[key], tc[theirs_key]), 10),
                          "agrees_to_1e-4": bool(
                              abs(mc[key] - tc[theirs_key]) <= max(TOL * abs(tc[theirs_key]), 5e-5))})
    else:
        comps.append({"quantity": "influenza", "status": "ABSENT (not computed)"})

    checked = [c for c in comps if "agrees_to_1e-4" in c]
    n_ok = sum(c["agrees_to_1e-4"] for c in checked)
    out = {"n_compared": len(checked), "n_agreeing_to_1e-4": n_ok,
           "n_disagreeing": len(checked) - n_ok,
           "reference": str(REPORT), "tolerance_relative": TOL, "comparisons": comps,
           "disagreements": [c for c in checked if not c["agrees_to_1e-4"]]}
    (RESULTS / "reproduction_check.json").write_text(json.dumps(out, indent=1))

    blk = Block("reproduction_check", SCRIPT, COMMAND)
    blk.dataset(id="biology/reproduction_check_vs_earlier_validation", domain="biology",
                title="Cross-implementation reproduction check: this block vs the 2026-09-19 "
                      "genetics validation",
                source=str(REPORT), provenance="simulation", n_objects=len(checked),
                notes="Not a biological dataset: an analysis-level record comparing two "
                      "independently written implementations of 'Rips on the given metric' over "
                      "the same raw files.")
    blk.run(dataset_id="biology/reproduction_check_vs_earlier_validation", method="rips",
            coeff_field=2, max_dim=1,
            params={"check": "numerical agreement with the earlier validation's Rips cross-check",
                    "tolerance_relative": TOL, "quantities_compared": len(checked),
                    "metric": "n/a -- this run compares previously computed numbers"},
            preprocessing="none", tier="X",
            stats=[{"name": "n_quantities_compared", "value": float(len(checked))},
                   {"name": "n_agreeing_to_relative_1e-4", "value": float(n_ok)},
                   {"name": "n_disagreeing", "value": float(len(checked) - n_ok)}],
            controls=[{"kind": "known_answer",
                       "description": "two independent implementations must give the same "
                                      "persistence statistics on the same raw files",
                       "passed": bool(n_ok == len(checked) and checked),
                       "detail": f"{n_ok}/{len(checked)} quantities agree to relative {TOL}; "
                                 f"disagreements: {json.dumps(out['disagreements'])[:400]}"}],
            findings=[{"claim": f"{n_ok} of {len(checked)} persistence statistics computed here "
                                f"agree with the 2026-09-19 genetics validation to a relative "
                                f"{TOL}, using independently written code on the same raw files.",
                       "verdict": "recovered" if (checked and n_ok == len(checked))
                       else "inconclusive",
                       "tier": "X",
                       "caveat": "this checks implementation agreement, not biological truth; both "
                                 "implementations could share a conceptual error, and the earlier "
                                 "run's own verdicts (Caulobacter FAIL, E. coli PASS, influenza "
                                 "PASS) are unchanged by it",
                       "reference": str(REPORT)}])
    print(json.dumps({k: v for k, v in out.items() if k != "comparisons"}, indent=1)[:2000])
    blk.write()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
