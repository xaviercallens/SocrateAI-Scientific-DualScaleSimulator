"""Shared helpers for the astro persistent-homology campaign.

Rules this module enforces mechanically (so no individual script can forget):
  * the TopoDB connection is OPENED, WRITTEN and CLOSED around each ingest
    (three agents share the file in WAL mode; write transactions stay short);
  * every run records repo=WORKTREE, so `git rev-parse HEAD` is taken in the
    astro worktree and not in whatever directory the shell happens to sit in;
  * a rank p-value carries its own floor, 1/(n_null+1), in the p_method string,
    so a p at the floor can never be read as "significant";
  * peak memory and wall time are measured, not guessed.

Every number produced downstream of this module is tier X (numerics).
"""
from __future__ import annotations

import contextlib
import hashlib
import json
import os
import resource
import subprocess
import sys
import time
from pathlib import Path

WORKTREE = "/mnt/disks/disk-socrateai-local-1/wt-topo-astro"
ASTRO = os.path.join(WORKTREE, "topodb_runs", "astro")
RESULTS = os.path.join(ASTRO, "results")
DB_PATH = "/mnt/disks/disk-socrateai-local-1/topodb/topodb.sqlite"

sys.path.insert(0, WORKTREE)


# ----------------------------------------------------------------- provenance
def sha256_file(path: str, _chunk: int = 1 << 22) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(_chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def peak_mb() -> float:
    """Peak RSS of this process, in MB (ru_maxrss is KB on Linux)."""
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


class Timer:
    def __enter__(self):
        self.t0 = time.time()
        return self

    def __exit__(self, *a):
        self.wall = time.time() - self.t0

    @property
    def sec(self) -> float:
        return getattr(self, "wall", time.time() - self.t0)


# --------------------------------------------------------------------- TopoDB
@contextlib.contextmanager
def topodb():
    """Open the shared DB, yield it, always close it."""
    from topodb.api import TopoDB
    db = TopoDB(DB_PATH)
    try:
        yield db
    finally:
        try:
            db.con.close()
        except Exception:
            pass


def add_run(db, **kw) -> int:
    """add_run with repo pinned to the astro worktree (advisor item 3)."""
    kw.setdefault("repo", WORKTREE)
    return db.add_run(**kw)


def rank_p(observed: float, null_values, tail: str = "two"):
    """Distribution-free rank p-value of `observed` in a null sample.

    tail='high'  : P(null >= obs)          -- more features than the null
    tail='low'   : P(null <= obs)
    tail='two'   : 2 * min(high, low), capped at 1.
    The (r+1)/(n+1) form is used, so the smallest attainable value is
    1/(n+1): that FLOOR is returned with the p-value and must be quoted.
    """
    import numpy as np
    x = np.asarray(null_values, dtype=float)
    n = x.size
    hi = (np.sum(x >= observed) + 1) / (n + 1)
    lo = (np.sum(x <= observed) + 1) / (n + 1)
    if tail == "high":
        p = hi
    elif tail == "low":
        p = lo
    else:
        p = min(1.0, 2 * min(hi, lo))
    return float(p), float(1.0 / (n + 1))


def p_method_str(n_null: int, tail: str) -> str:
    return f"rank ({tail}-tail, (r+1)/(n+1); floor 1/{n_null + 1} = {1.0 / (n_null + 1):.4g})"


# --------------------------------------------------- pre-declared Betti rule
def gap_betti(persistences, min_persistence: float = 0.0):
    """Largest-multiplicative-gap rule; PARAMETER-FREE and declared in
    PRE_DECLARED_STATISTICS.md before any dataset was read.

    Sort finite persistences descending; the number of bars above the largest
    multiplicative gap between consecutive entries is the reported Betti
    number.  Returns (betti, gap_ratio).

    KNOWN LIMITATION, measured elsewhere and repeated here so it is never
    forgotten: docs/FUTURE_OBSERVATIONAL_TARGETS.md records that a naive
    persistence-ratio rule FIRED ON 9 OF 12 density-matched nulls.  This rule
    is therefore used ONLY in step0_known_answers.py, where the true topology
    is known and `expected` is pre-declared.  It is NEVER used as a detector
    on real data; there the COUNT statistic against an explicit null is the
    primary record.
    """
    import numpy as np
    p = np.sort(np.asarray([q for q in persistences if np.isfinite(q) and q > min_persistence],
                           dtype=float))[::-1]
    if p.size == 0:
        return 0, 0.0
    if p.size == 1:
        return 1, float("inf")
    ratios = p[:-1] / np.maximum(p[1:], 1e-300)
    k = int(np.argmax(ratios))
    return k + 1, float(ratios[k])


def alpha_bars(simplex_tree, dim: int, sqrt_scale: bool = True):
    """(birth, death) pairs in dimension `dim` from a GUDHI simplex tree.

    GUDHI's AlphaComplex filtration value is the SQUARED circumradius.  With
    sqrt_scale=True the bars are returned in LENGTH units (the units of the
    input coordinates), which is what every statistic below uses.
    """
    import numpy as np
    d = np.asarray(simplex_tree.persistence_intervals_in_dimension(dim), dtype=float).reshape(-1, 2)
    if sqrt_scale and d.size:
        d = np.sqrt(np.maximum(d, 0.0))
    return d


def save_json(name: str, obj) -> str:
    Path(RESULTS).mkdir(parents=True, exist_ok=True)
    path = os.path.join(RESULTS, name)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=1, sort_keys=True, default=float)
    return path


def git_head() -> str:
    return subprocess.run(["git", "rev-parse", "HEAD"], cwd=WORKTREE,
                          capture_output=True, text=True).stdout.strip()
