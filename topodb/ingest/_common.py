"""Shared machinery for the backfill modules.

Rules every ingest module in this package obeys:

* Every number written to the DB is read out of a committed source JSON at
  ingest time. No module hard-codes a measured value, and none of them recompute
  anything: this package does no persistent homology.
* `expected=` is passed to `add_betti` ONLY where the source file itself carries
  a pre-stated expectation (an expectations.json / registration written before
  the run, or an `expected` field the source declares as pre-registered). Where
  the source labels a result post-hoc, no expectation is passed, so `matches`
  stays NULL.
* A p_value is carried only together with the null model and the number of null
  draws that the SAME source states. Where a source gives a statistic with no
  null, the statistic is stored alone.
* The word "proved" is never written.

Re-runnability, and NOT destroying anyone else's rows. `add_run` and `add_bars`
in api.py are plain INSERTs, so running the backfill twice would duplicate rows.
The database also lives outside git and is SHARED: other sessions write their own
runs into it (protein-structure runs under `topodb_runs/` were found there while
this backfill was being written). So the reset is OWNERSHIP-SCOPED, never global:

  * every run this package creates carries `_backfill_source` in its params_json;
  * `run_all` writes a sidecar manifest of every id it created, next to the
    database at OWNED_MANIFEST;
  * `wipe_owned()` deletes ONLY runs carrying that marker, their dependent rows,
    and the datasets and dataset-level findings named in the manifest that have
    nothing else pointing at them.

Nothing that this package did not create is ever deleted.
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import subprocess
from pathlib import Path

from ..api import TopoDB

DISK = "/mnt/disks/disk-socrateai-local-1"
OWNED_MANIFEST = os.environ.get(
    "TOPODB_BACKFILL_MANIFEST",
    "/mnt/disks/disk-socrateai-local-1/topodb/backfill_owned.json")
MARKER = "_backfill_source"
WORKTREES = {
    "tdasimple": f"{DISK}/dualscale-wt-tdasimple",
    "tdaval": f"{DISK}/dualscale-wt-tdaval",
    "tdak3t2": f"{DISK}/dualscale-wt-tdak3t2",
    "reverse": f"{DISK}/dualscale-wt-reverse",
    "k3t2": f"{DISK}/dualscale-wt-k3t2",
}

def load_manifest(path: str = OWNED_MANIFEST) -> dict:
    try:
        with open(path) as fh:
            return json.load(fh)
    except Exception:
        return {}


def save_manifest(manifest: dict, path: str = OWNED_MANIFEST) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(manifest, fh, indent=2)
        fh.write("\n")


def wipe_owned(db: TopoDB, manifest_path: str = OWNED_MANIFEST) -> dict:
    """Delete ONLY the rows a previous run of this backfill created.

    The database is shared with other sessions, so a global wipe would destroy
    their work. Ownership is established two ways, and both are required to be
    safe: a run must carry the `_backfill_source` marker in its params_json, and
    a dataset or dataset-level finding must be named in the sidecar manifest AND
    have nothing else referring to it.
    """
    owned_runs = [int(r[0]) for r in db.con.execute(
        "SELECT id FROM run WHERE params_json LIKE ?", (f'%"{MARKER}"%',)).fetchall()]
    removed = {"runs": len(owned_runs), "datasets": 0, "findings": 0}
    if owned_runs:
        qs = ",".join("?" * len(owned_runs))
        for t in ("bar", "betti", "statistic", "control"):
            db.con.execute(f"DELETE FROM {t} WHERE run_id IN ({qs})", owned_runs)
        cur = db.con.execute(f"DELETE FROM finding WHERE run_id IN ({qs})", owned_runs)
        removed["findings"] += cur.rowcount if cur.rowcount and cur.rowcount > 0 else 0
        db.con.execute(f"DELETE FROM run WHERE id IN ({qs})", owned_runs)

    man = load_manifest(manifest_path)
    for fid in man.get("dataset_level_findings") or []:
        cur = db.con.execute("DELETE FROM finding WHERE id=?", (fid,))
        removed["findings"] += cur.rowcount if cur.rowcount and cur.rowcount > 0 else 0
        db.con.execute("DELETE FROM search WHERE kind='finding' AND ref_id=?", (str(fid),))
    for did in man.get("datasets") or []:
        still_used = db.con.execute(
            "SELECT (SELECT count(*) FROM run WHERE dataset_id=?) + "
            "(SELECT count(*) FROM finding WHERE dataset_id=?)", (did, did)).fetchone()[0]
        if still_used:
            continue
        cur = db.con.execute("DELETE FROM dataset WHERE id=?", (did,))
        removed["datasets"] += cur.rowcount if cur.rowcount and cur.rowcount > 0 else 0
        db.con.execute("DELETE FROM search WHERE kind='dataset' AND ref_id=?", (did,))
    db.con.commit()
    return removed


def load(path: str):
    with open(path) as fh:
        return json.load(fh)


def sha256(path: str) -> str | None:
    try:
        h = hashlib.sha256()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        return h.hexdigest()
    except Exception:
        return None


def head(worktree: str) -> str | None:
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], cwd=worktree, capture_output=True,
                              text=True, timeout=10).stdout.strip() or None
    except Exception:
        return None


_NUM = re.compile(r"=\s*(-?(?:inf|[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?))")


def rule_value(rule: str | None) -> float | None:
    """Pull the numeric value out of a rule string like 'p_1 / p_2 = 17.6 >= 5'.

    Returns None when the string carries no parseable number; the full rule text
    is always kept verbatim in params_json regardless.
    """
    if not isinstance(rule, str):
        return None
    m = _NUM.search(rule)
    if not m:
        return None
    try:
        return float(m.group(1))
    except ValueError:
        return None


def betti_map(seq) -> dict[int, int]:
    """[1, 2, 1] -> {0: 1, 1: 2, 2: 1}; drops non-integers."""
    out = {}
    for i, v in enumerate(seq or []):
        if isinstance(v, bool):
            continue
        if isinstance(v, int):
            out[i] = v
    return out


class Source:
    """A counted wrapper around TopoDB: one per backfill source, for the report."""

    def __init__(self, db: TopoDB, name: str, worktree_key: str):
        self.db = db
        self.name = name
        self.repo = WORKTREES[worktree_key]
        self.git_head = head(self.repo)
        self.counts = dict(datasets=0, runs=0, betti=0, bars=0, statistics=0,
                           controls=0, findings=0)
        self.skipped: list[dict] = []
        self.discrepancies: list[dict] = []
        # ids created here, so a re-run can delete exactly these and nothing else
        self.owned_datasets: list[str] = []
        self.owned_runs: list[int] = []
        self.owned_dataset_level_findings: list[int] = []

    # -- provenance helpers ------------------------------------------------
    def rel(self, *parts: str) -> str:
        return str(Path(self.repo).joinpath(*parts))

    def skip(self, what: str, why: str, **extra) -> None:
        self.skipped.append({"what": what, "why": why, **extra})

    def note_discrepancy(self, where: str, detail: str) -> None:
        self.discrepancies.append({"where": where, "detail": detail})

    # -- thin, counting proxies -------------------------------------------
    def dataset(self, **kw) -> str:
        self.counts["datasets"] += 1
        did = self.db.add_dataset(**kw)
        self.owned_datasets.append(did)
        return did

    def run(self, **kw) -> int:
        kw.setdefault("repo", self.repo)
        # ownership marker: wipe_owned deletes only runs carrying it
        kw["params"] = dict(kw.get("params") or {}, **{MARKER: self.name})
        rid = self.db.add_run(**kw)
        self.counts["runs"] += 1
        self.owned_runs.append(rid)
        return rid

    def betti(self, run_id: int, betti: dict, expected: dict | None = None) -> None:
        if not betti:
            return
        self.db.add_betti(run_id, betti, expected)
        self.counts["betti"] += len(betti)

    def bars(self, run_id: int, dim: int, bars) -> None:
        bars = [b for b in bars if b is not None]
        if not bars:
            return
        before = self.db.con.execute("SELECT count(*) FROM bar WHERE run_id=?", (run_id,)).fetchone()[0]
        self.db.add_bars(run_id, dim, bars)
        after = self.db.con.execute("SELECT count(*) FROM bar WHERE run_id=?", (run_id,)).fetchone()[0]
        self.counts["bars"] += after - before

    def stat(self, run_id: int, name: str, value, **kw) -> None:
        if value is None:
            return
        try:
            value = float(value)
        except (TypeError, ValueError):
            return
        if value != value:      # NaN
            return
        self.db.add_statistic(run_id, name, value, **kw)
        self.counts["statistics"] += 1

    def control(self, run_id: int, kind: str, description: str, passed, detail=None) -> None:
        self.db.add_control(run_id, kind, description, passed, detail)
        self.counts["controls"] += 1

    def finding(self, **kw) -> int:
        fid = self.db.add_finding(**kw)
        self.counts["findings"] += 1
        if kw.get("run_id") is None:
            # not reachable through an owned run, so it goes in the manifest
            self.owned_dataset_level_findings.append(fid)
        return fid

    def report(self) -> dict:
        return {"source": self.name, "worktree": self.repo, "git_head": self.git_head,
                "ingested": dict(self.counts), "skipped": self.skipped,
                "discrepancies": self.discrepancies, "_owned": self.owned()}

    def owned(self) -> dict:
        return {"datasets": self.owned_datasets, "runs": self.owned_runs,
                "dataset_level_findings": self.owned_dataset_level_findings}
