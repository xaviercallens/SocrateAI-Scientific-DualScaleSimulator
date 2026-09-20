"""TopoDB ingestion and query API.

One rule governs this module: a record may only assert what a control checked.
- `p_value` is refused unless `null_model` and `n_null` are given.
- `matches` is refused unless `expected` is given.
- `tier` is mandatory and free text is not accepted for it.

Usage (ingestion):
    from topodb.api import TopoDB
    db = TopoDB()                      # default path, or TopoDB(path)
    db.add_dataset(id="astro/desi_dr1_bgs_ngc", domain="astro", title=..., source=..., provenance="observation", ...)
    run = db.add_run(dataset_id=..., method="alpha", coeff_field=2, params={...}, script=..., command=..., tier="X", ...)
    db.add_betti(run, {0: 1, 1: 2, 2: 1}, expected={0: 1, 1: 2, 2: 1})
    db.add_bars(run, dim=1, bars=[(birth, death), ...])       # ranked by persistence automatically
    db.add_statistic(run, "iqr_over_median", 0.454, null_model="500 spectrum-matched Gaussian sims",
                     n_null=500, p_value=0.535, p_method="rank", multiplicity="Bonferroni/6")
    db.add_control(run, "known_answer", "circle recovers b1=1", passed=True, detail="ratio 47.9")
    db.add_finding(run_id=run, claim="...", verdict="null", tier="X", caveat="...")

Query:
    db.search("vortex")                 # full-text over titles and notes
    db.betti_table(domain="biology")
    db.summary()
"""
from __future__ import annotations

import datetime as _dt
import json
import os
import sqlite3
import subprocess
from pathlib import Path

DEFAULT_DB = os.environ.get("TOPODB_PATH", "/mnt/disks/disk-socrateai-local-1/topodb/topodb.sqlite")
SCHEMA = Path(__file__).with_name("schema.sql")

TIERS = {"A", "B", "L", "C", "X"}
METHODS = {"alpha", "rips", "sparse_rips", "witness", "cubical", "lower_star_graph", "chain_complex", "mapper"}


def _git(cwd: str | None, *args: str) -> str | None:
    try:
        return subprocess.run(["git", *args], cwd=cwd, capture_output=True, text=True, timeout=10).stdout.strip() or None
    except Exception:
        return None


class TopoDB:
    def __init__(self, path: str = DEFAULT_DB):
        self.path = path
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        self.con = sqlite3.connect(path, timeout=60)
        self.con.row_factory = sqlite3.Row
        self.con.execute("PRAGMA foreign_keys=ON")
        if SCHEMA.exists():
            self.con.executescript(SCHEMA.read_text())
        self.con.commit()

    # ---------------------------------------------------------------- ingest
    def add_dataset(self, *, id: str, domain: str, title: str, source: str, provenance: str,
                    n_objects: int | None = None, ambient_dim: int | None = None, units: str | None = None,
                    sha256: str | None = None, local_path: str | None = None, notes: str | None = None) -> str:
        self.con.execute(
            "INSERT OR REPLACE INTO dataset(id,domain,title,source,provenance,n_objects,ambient_dim,units,sha256,local_path,notes)"
            " VALUES (?,?,?,?,?,?,?,?,?,?,?)",
            (id, domain, title, source, provenance, n_objects, ambient_dim, units, sha256, local_path, notes))
        self._index("dataset", id, title, " ".join(filter(None, [source, notes, units])))
        self.con.commit()
        return id

    def add_run(self, *, dataset_id: str, method: str, coeff_field: int, params: dict, script: str, command: str,
                tier: str, max_dim: int | None = None, preprocessing: str | None = None, seed: str | None = None,
                wall_sec: float | None = None, peak_mb: float | None = None, repo: str | None = None) -> int:
        if tier not in TIERS:
            raise ValueError(f"tier must be one of {sorted(TIERS)}, got {tier!r}")
        if method not in METHODS:
            raise ValueError(f"method must be one of {sorted(METHODS)}, got {method!r}")
        cur = self.con.execute(
            "INSERT INTO run(dataset_id,method,coeff_field,max_dim,params_json,preprocessing,script,command,seed,"
            "git_commit,branch,wall_sec,peak_mb,tier,created) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
            (dataset_id, method, coeff_field, max_dim, json.dumps(params, sort_keys=True), preprocessing, script,
             command, seed, _git(repo, "rev-parse", "HEAD"), _git(repo, "rev-parse", "--abbrev-ref", "HEAD"),
             wall_sec, peak_mb, tier, _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds")))
        self.con.commit()
        return int(cur.lastrowid)

    def add_betti(self, run_id: int, betti: dict[int, int], expected: dict[int, int] | None = None) -> None:
        for dim, value in betti.items():
            exp = None if expected is None else expected.get(dim)
            self.con.execute("INSERT OR REPLACE INTO betti(run_id,dim,value,expected,matches) VALUES (?,?,?,?,?)",
                             (run_id, dim, int(value), exp, None if exp is None else int(exp == value)))
        self.con.commit()

    def add_bars(self, run_id: int, dim: int, bars, top: int | None = 50) -> None:
        """bars: iterable of (birth, death); death=None or inf means infinite. Ranked by persistence."""
        rows = []
        for b, d in bars:
            d = None if d is None or d == float("inf") else float(d)
            rows.append((float(b), d))
        rows.sort(key=lambda bd: (float("inf") if bd[1] is None else bd[1] - bd[0]), reverse=True)
        if top is not None:
            rows = rows[:top]
        self.con.executemany("INSERT INTO bar(run_id,dim,birth,death,rank) VALUES (?,?,?,?,?)",
                             [(run_id, dim, b, d, i + 1) for i, (b, d) in enumerate(rows)])
        self.con.commit()

    def add_statistic(self, run_id: int, name: str, value: float, *, null_model: str | None = None,
                      n_null: int | None = None, p_value: float | None = None, p_method: str | None = None,
                      multiplicity: str | None = None) -> None:
        if p_value is not None and (null_model is None or n_null is None):
            raise ValueError("a p_value requires null_model and n_null: an unexplained p-value is not storable")
        self.con.execute("INSERT OR REPLACE INTO statistic(run_id,name,value,null_model,n_null,p_value,p_method,multiplicity)"
                         " VALUES (?,?,?,?,?,?,?,?)",
                         (run_id, name, float(value), null_model, n_null, p_value, p_method, multiplicity))
        self.con.commit()

    def add_control(self, run_id: int, kind: str, description: str, passed: bool | None, detail: str | None = None) -> None:
        self.con.execute("INSERT OR REPLACE INTO control(run_id,kind,description,passed,detail) VALUES (?,?,?,?,?)",
                         (run_id, kind, description, None if passed is None else int(passed), detail))
        self.con.commit()

    def add_finding(self, *, claim: str, verdict: str, tier: str, run_id: int | None = None,
                    dataset_id: str | None = None, caveat: str | None = None, reference: str | None = None) -> int:
        if tier not in TIERS:
            raise ValueError(f"tier must be one of {sorted(TIERS)}")
        cur = self.con.execute("INSERT INTO finding(run_id,dataset_id,claim,verdict,tier,caveat,reference)"
                               " VALUES (?,?,?,?,?,?,?)", (run_id, dataset_id, claim, verdict, tier, caveat, reference))
        self._index("finding", str(cur.lastrowid), claim, " ".join(filter(None, [caveat, reference])))
        self.con.commit()
        return int(cur.lastrowid)

    def _index(self, kind: str, ref_id: str, title: str, text: str | None) -> None:
        self.con.execute("DELETE FROM search WHERE kind=? AND ref_id=?", (kind, ref_id))
        self.con.execute("INSERT INTO search(kind,ref_id,title,text) VALUES (?,?,?,?)", (kind, ref_id, title, text or ""))

    # ----------------------------------------------------------------- query
    def search(self, query: str, limit: int = 20):
        return self.con.execute("SELECT kind,ref_id,title,snippet(search,3,'[',']','...',12) AS hit"
                                " FROM search WHERE search MATCH ? LIMIT ?", (query, limit)).fetchall()

    def betti_table(self, domain: str | None = None, limit: int = 100):
        q = ("SELECT d.domain, d.id AS dataset, r.id AS run, r.method, r.coeff_field,"
             " group_concat(b.dim || ':' || b.value, ' ') AS betti, r.tier"
             " FROM run r JOIN dataset d ON d.id=r.dataset_id LEFT JOIN betti b ON b.run_id=r.id")
        args: tuple = ()
        if domain:
            q += " WHERE d.domain=?"
            args = (domain,)
        q += " GROUP BY r.id ORDER BY d.domain, d.id LIMIT ?"
        return self.con.execute(q, (*args, limit)).fetchall()

    def summary(self):
        return {
            "datasets_by_domain": dict(self.con.execute(
                "SELECT domain, count(*) FROM dataset GROUP BY domain").fetchall()),
            "runs_by_method": dict(self.con.execute(
                "SELECT method, count(*) FROM run GROUP BY method").fetchall()),
            "runs_by_tier": dict(self.con.execute(
                "SELECT tier, count(*) FROM run GROUP BY tier").fetchall()),
            "controls_passed": dict(self.con.execute(
                "SELECT CASE passed WHEN 1 THEN 'passed' WHEN 0 THEN 'failed' ELSE 'not_run' END, count(*)"
                " FROM control GROUP BY passed").fetchall()),
            "findings_by_verdict": dict(self.con.execute(
                "SELECT verdict, count(*) FROM finding GROUP BY verdict").fetchall()),
            "statistics_with_p": self.con.execute(
                "SELECT count(*) FROM statistic WHERE p_value IS NOT NULL").fetchone()[0],
            "bars": self.con.execute("SELECT count(*) FROM bar").fetchone()[0],
        }
