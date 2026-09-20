"""Ingest one results json into TopoDB.

Opens the shared SQLite file, writes, closes.  Compute never holds the lock.

Idempotency: TopoDB.add_run/add_bars/add_finding are plain INSERTs, so a re-ingest
would duplicate them.  This script therefore first deletes every run (and its
children) previously written by the SAME script path, then re-inserts.  Re-running
an ingest is safe and the run counts stay truthful.

Usage:  python ingest.py results/<block>.json [--dry-run]
"""
from __future__ import annotations

import json
import sqlite3
import sys
import time
from pathlib import Path

sys.path.insert(0, "/mnt/disks/disk-socrateai-local-1/wt-topo-bio")
from topodb.api import TopoDB  # noqa: E402

WT = "/mnt/disks/disk-socrateai-local-1/wt-topo-bio"


def _retry(fn, *a, **k):
    for attempt in range(8):
        try:
            return fn(*a, **k)
        except sqlite3.OperationalError as exc:
            if "locked" not in str(exc) and "busy" not in str(exc):
                raise
            time.sleep(0.5 * (attempt + 1))
    raise RuntimeError("database stayed locked after 8 attempts")


def purge(db: TopoDB, script: str) -> int:
    ids = [r[0] for r in db.con.execute("SELECT id FROM run WHERE script=?", (script,)).fetchall()]
    if not ids:
        return 0
    qs = ",".join("?" * len(ids))
    for tbl in ("bar", "betti", "statistic", "control", "finding"):
        db.con.execute(f"DELETE FROM {tbl} WHERE run_id IN ({qs})", ids)
    db.con.execute(f"DELETE FROM run WHERE id IN ({qs})", ids)
    db.con.commit()
    return len(ids)


def main() -> int:
    path = Path(sys.argv[1])
    dry = "--dry-run" in sys.argv
    payload = json.loads(path.read_text())
    script = payload["script"]

    if dry:
        print(f"[dry] {path.name}: {len(payload['datasets'])} datasets, {len(payload['runs'])} runs")
        return 0

    db = TopoDB()
    try:
        n_purged = _retry(purge, db, script)
        for ds in payload["datasets"]:
            _retry(db.add_dataset, **ds)
        n_runs = 0
        for r in payload["runs"]:
            rid = _retry(db.add_run, dataset_id=r["dataset_id"], method=r["method"],
                         coeff_field=r["coeff_field"], params=r["params"], script=r["script"],
                         command=r["command"], tier=r["tier"], max_dim=r.get("max_dim"),
                         preprocessing=r.get("preprocessing"), seed=r.get("seed"),
                         wall_sec=r.get("wall_sec"), peak_mb=r.get("peak_mb"), repo=WT)
            n_runs += 1
            if r.get("betti"):
                exp = {int(k): int(v) for k, v in (r.get("expected_betti") or {}).items()} or None
                _retry(db.add_betti, rid, {int(k): int(v) for k, v in r["betti"].items()}, expected=exp)
            for dim, bars in (r.get("bars") or {}).items():
                if bars:
                    _retry(db.add_bars, rid, int(dim), [(b, d) for b, d in bars], top=20)
            for s in r.get("stats") or []:
                _retry(db.add_statistic, rid, s["name"], s["value"], null_model=s.get("null_model"),
                       n_null=s.get("n_null"), p_value=s.get("p_value"),
                       p_method=s.get("p_method"), multiplicity=s.get("multiplicity"))
            for c in r.get("controls") or []:
                _retry(db.add_control, rid, c["kind"], c["description"], c.get("passed"), c.get("detail"))
            for f in r.get("findings") or []:
                _retry(db.add_finding, run_id=rid, dataset_id=r["dataset_id"], claim=f["claim"],
                       verdict=f["verdict"], tier=f["tier"], caveat=f.get("caveat"),
                       reference=f.get("reference"))
        print(f"[ingest] {path.name}: purged {n_purged} old runs, wrote "
              f"{len(payload['datasets'])} datasets / {n_runs} runs")
    finally:
        db.con.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
