"""Assemble the biology block's report from the per-block result files and TopoDB.

Reads only what the compute scripts wrote; computes nothing new except counts.
Opens the shared database briefly for the tallies, then closes it.

Usage: python make_report.py
"""
from __future__ import annotations

import json
import sqlite3
import sys
from pathlib import Path

sys.path.insert(0, "/mnt/disks/disk-socrateai-local-1/wt-topo-bio")
sys.path.insert(0, "/mnt/disks/disk-socrateai-local-1/wt-topo-bio/topodb_runs/biology")
from bio_common import RESULTS, git_head  # noqa: E402
from topodb.api import TopoDB  # noqa: E402

OUT = Path("/mnt/disks/disk-socrateai-local-1/wt-topo-bio/topodb_runs/biology/report.json")


def load(name):
    p = RESULTS / name
    return json.loads(p.read_text()) if p.exists() else None


def main() -> int:
    db = TopoDB()
    try:
        rows = db.con.execute(
            "SELECT d.id, d.domain, d.provenance, count(r.id) FROM dataset d "
            "LEFT JOIN run r ON r.dataset_id = d.id GROUP BY d.id ORDER BY d.id").fetchall()
        datasets = [{"dataset": r[0], "domain": r[1], "provenance": r[2], "runs": r[3]}
                    for r in rows]
        summary = db.summary()
        by_verdict = dict(db.con.execute(
            "SELECT verdict, count(*) FROM finding GROUP BY verdict").fetchall())
        n_stats_with_p = db.con.execute(
            "SELECT count(*) FROM statistic WHERE p_value IS NOT NULL").fetchone()[0]
        n_stats_without_p = db.con.execute(
            "SELECT count(*) FROM statistic WHERE p_value IS NULL").fetchone()[0]
        controls = dict(db.con.execute(
            "SELECT CASE passed WHEN 1 THEN 'passed' WHEN 0 THEN 'failed' ELSE 'not_run' END, "
            "count(*) FROM control GROUP BY passed").fetchall())
    finally:
        db.con.close()

    prot = load("proteins_table.json")
    ext = load("proteins_ext_table.json")
    rep = {
        "title": "Persistent homology over a wide variety of real biological datasets",
        "branch": "topo/bio", "git_commit": git_head(),
        "preregistration": "topodb_runs/biology/expectations.json, committed alone at dde7a07 "
                           "before any biological persistence diagram was computed",
        "binding_method_rules": [
            "(a) A non-Euclidean metric (Hi-C contact, scRNA correlation, sequence Hamming) is "
            "given to Rips DIRECTLY. No MDS-then-alpha anywhere: that route produced a FALSE loop "
            "of dominance 11.1 on the tree-like influenza HA segment in the earlier validation. "
            "bio_common.py contains no MDS function.",
            "(b) Alpha is used only on genuinely Euclidean coordinates (protein/RNA atoms).",
            "gudhi: create_simplex_tree(max_dimension = max_hom_dim + 1), since a d-skeleton "
            "yields homology only to d-1. Step 0 runs through the same helpers, so a regression "
            "would show there as b1 = 0 on the circle.",
        ],
        "step0_known_answer": load("step0_verdicts.json"),
        "database": {"path": "/mnt/disks/disk-socrateai-local-1/topodb/topodb.sqlite",
                     "datasets_with_run_counts": datasets, "summary": summary,
                     "findings_by_verdict": by_verdict, "controls": controls,
                     "statistics_with_p_value": n_stats_with_p,
                     "statistics_without_p_value_no_null_defined": n_stats_without_p},
        "protein_fold_separation": {
            "registered": prot["separation"] if prot else None,
            "post_hoc_extension": ext["analysis"] if ext else None,
        },
        "hic_circular_vs_linear": load("hic_summary.json"),
        "influenza": load("influenza_summary.json"),
        "scrna": load("scrna_summary.json"),
        "ecg": load("ecg_summary.json"),
        "synthetic_control": load("synthetic_control.json"),
        "reproduction_check": (lambda r: {k: v for k, v in r.items() if k != "comparisons"}
                               if r else None)(load("reproduction_check.json")),
    }
    OUT.write_text(json.dumps(rep, indent=1))
    print(f"wrote {OUT}")
    print(json.dumps({"datasets": len(datasets), "summary": summary,
                      "findings_by_verdict": by_verdict, "controls": controls}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
