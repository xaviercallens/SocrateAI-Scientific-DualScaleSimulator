"""Backfill of TopoDB from the committed results of this project.

Nothing in this package computes persistent homology. Each module reads a
committed JSON file from a source worktree (read-only) and writes what that file
records into the database, with its provenance and its own verdict.

    python -m topodb.ingest            # wipe and re-ingest everything
    python -m topodb.ingest --notes topodb/BACKFILL_NOTES.json
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json

from ..api import TopoDB
from . import _common as C

MODULES = [
    ("a1", "simple_suite"),
    ("a2", "tda_fixed"),
    ("b1", "quantum_fluid"),
    ("b2", "genetics"),
    ("b3", "crystallography"),
    ("c", "k3t2_tda"),
    ("d1", "cosmic_vorticity"),
    ("d2", "reverse_zero"),
    ("e", "k3t2_rigidity"),
]


def run_all(db: TopoDB | None = None, wipe: bool = True, only: list[str] | None = None) -> dict:
    import importlib
    db = db or TopoDB()
    removed = None
    if wipe and not only:
        # ownership-scoped: this deletes only rows a previous backfill created and
        # never rows another session wrote into this shared database
        removed = C.wipe_owned(db)
    reports = []
    owned = {"datasets": [], "runs": [], "dataset_level_findings": []}
    for key, mod in MODULES:
        if only and key not in only and mod not in only:
            continue
        m = importlib.import_module(f"{__package__}.{mod}")
        r = m.ingest(db)
        for k, v in (r.pop("_owned", None) or {}).items():
            owned[k].extend(v)
        reports.append(r)
    if not only:
        C.save_manifest({"db_path": db.path,
                         "generated": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
                         **owned})
    foreign_runs = db.con.execute(
        "SELECT count(*) FROM run WHERE params_json NOT LIKE ?", (f'%"{C.MARKER}"%',)).fetchone()[0]
    return {
        "_what_this_is": (
            "A machine-written record of the TopoDB backfill: what each committed source "
            "contributed, what was skipped and why, and every place where a source's own "
            "numbers disagree with its narrative. No number in the database was computed "
            "here; all of them were read from the committed files named below."
        ),
        "generated": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        "db_path": db.path,
        "rules_applied": [
            "expected= is passed to add_betti only where the source itself pre-stated the "
            "expectation, so `matches` is NULL wherever no pre-registration existed",
            "a p_value is stored only with the null model and the null count the same source "
            "states; where a source reports a statistic with no null, the statistic is stored alone",
            "runs whose source records an error, timeout or crash are listed under `skipped`, "
            "not ingested as results",
            "api.add_run has no `notes` column, so a superseded run is marked in params_json "
            "under 'superseded', in the dataset notes, and in the caveat of its finding",
            "api.add_bars truncates to the 50 longest bars per dimension; several sources had "
            "already truncated further before committing, and each module says so",
            "the word 'proved' is asserted nowhere; where it appears it is a verbatim "
            "quotation of a source refusing the word (\"never 'proved'\", \"nothing here is proved\")",
            "the database is shared with other sessions, so the reset before an ingest is "
            "ownership-scoped: it deletes only runs carrying the `_backfill_source` marker and "
            "only datasets and dataset-level findings named in the sidecar manifest that nothing "
            "else refers to",
        ],
        "shared_database": {
            "note": "this database lives outside git and other sessions write to it",
            "rows_removed_by_the_ownership_scoped_reset_before_this_ingest": removed,
            "runs_present_that_this_backfill_does_not_own": foreign_runs,
            "owned_manifest": C.OWNED_MANIFEST,
            "manifest_is_load_bearing": (
                "runs carry the `_backfill_source` marker in params_json and can always be "
                "identified, but dataset-level findings (those with run_id NULL) carry no marker: "
                "they are identified only by the sidecar manifest. If that file is lost, the next "
                "reset leaves them behind as orphans and adds a fresh set."),
        },
        "sources": reports,
        "totals": {k: sum(r["ingested"][k] for r in reports) for k in
                   ("datasets", "runs", "betti", "bars", "statistics", "controls", "findings")},
    }


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="topodb.ingest", description=__doc__)
    p.add_argument("--db", default=None)
    p.add_argument("--notes", default=None, help="write the backfill notes JSON here")
    p.add_argument("--only", nargs="*", help="module keys to run (skips the wipe)")
    p.add_argument("--no-wipe", action="store_true")
    a = p.parse_args(argv)
    db = TopoDB(a.db) if a.db else TopoDB()
    notes = run_all(db, wipe=not a.no_wipe, only=a.only)
    if a.notes:
        with open(a.notes, "w") as fh:
            json.dump(notes, fh, indent=2)
            fh.write("\n")
        print(f"wrote {a.notes}")
    print(json.dumps(notes["totals"], indent=2))
    for r in notes["sources"]:
        print(f"  {r['source']:<24} {r['ingested']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
