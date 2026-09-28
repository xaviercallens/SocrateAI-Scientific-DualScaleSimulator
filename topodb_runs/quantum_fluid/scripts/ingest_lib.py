"""Short-transaction TopoDB ingestion helpers.

Three agents share /mnt/disks/disk-socrateai-local-1/topodb/topodb.sqlite in WAL
mode, so: open, write, close. Never hold a connection across a computation.

`add_run` is a plain INSERT with AUTOINCREMENT, so re-running an ingest script
would silently duplicate every run, betti, bar and statistic. `IdGuard` records
the run_id under a stable key in results/ingest_ids.json and refuses to insert
a second time unless the caller passes force=True.
"""
from __future__ import annotations

import json
import os
import sys

sys.path.insert(0, "/mnt/disks/disk-socrateai-local-1/wt-topo-qfluid")
from topodb.api import TopoDB  # noqa: E402

import qf_lib as q  # noqa: E402

IDS = os.path.join(q.RESULTS, "ingest_ids.json")


class IdGuard:
    def __init__(self, path: str = IDS):
        self.path = path
        self.d = json.load(open(path)) if os.path.exists(path) else {}

    def has(self, key: str) -> bool:
        return key in self.d

    def get(self, key: str):
        return self.d.get(key)

    def put(self, key: str, run_id: int):
        self.d[key] = run_id
        with open(self.path, "w") as fh:
            json.dump(self.d, fh, indent=1, sort_keys=True)


def db() -> TopoDB:
    return TopoDB()


def close(d: TopoDB):
    try:
        d.con.close()
    except Exception:
        pass


def add_stats(d: TopoDB, run, items):
    """items: list of dicts {name, value, [null_model, n_null, p_value, p_method, multiplicity]}."""
    for it in items:
        if it.get("value") is None:
            continue
        d.add_statistic(run, it["name"], it["value"],
                        null_model=it.get("null_model"), n_null=it.get("n_null"),
                        p_value=it.get("p_value"), p_method=it.get("p_method"),
                        multiplicity=it.get("multiplicity"))
