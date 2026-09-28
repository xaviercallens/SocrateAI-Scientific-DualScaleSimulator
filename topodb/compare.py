"""Compare stored persistence diagrams between TopoDB runs.

WHAT THESE DISTANCES ARE COMPUTED ON, AND WHY THEY ARE APPROXIMATE
------------------------------------------------------------------
These helpers do NOT recompute persistent homology. They read the bars that
were stored in the `bar` table at ingestion time, and those bars are lossy in
two declared ways:

1. TRUNCATION. `topodb.api.TopoDB.add_bars` keeps only the `top` longest bars
   per (run, dimension), and its default is `top=50`. Every diagram in the DB
   that originally had more than 50 bars in a dimension is therefore a
   truncated diagram: its short bars, which sit near the diagonal, are gone.
   Bottleneck and Wasserstein distances are NOT invariant under dropping
   near-diagonal points. Dropping a bar of persistence p can move the
   bottleneck distance by up to p/2 and the q-Wasserstein distance by up to
   (p/2) in the q-norm contribution. So: the numbers here are a LOWER-ish
   approximation of the distance between the true diagrams, with an error
   bounded (per dropped bar) by half that bar's persistence -- which is not
   recorded once the bar is dropped. Treat them as a similarity heuristic for
   *finding* related runs, never as a measured invariant of the underlying
   spaces.

2. INFINITE BARS. The schema stores an infinite bar as `death IS NULL`.
   `gudhi.bottleneck_distance` returns `inf` whenever an essential class in one
   diagram has no partner in the other, so any dim-0 comparison would be `inf`
   almost everywhere. The default policy here is `infinite='drop'`: essential
   bars are removed before the distance is taken, and the count of dropped bars
   is reported alongside the distance so the caller can see what was discarded.
   `infinite='cap'` replaces the death of every essential bar by
   `cap_value` (default: max finite death over both diagrams) instead.

3. INCOMPARABLE FILTRATIONS. Runs in this DB come from alpha complexes, Rips
   complexes, cubical complexes on maps, lower-star filtrations on graphs and
   chain complexes with no filtration at all. Their filtration values are in
   different, unrelated units (squared alpha radius, Euclidean distance, map
   temperature, |contact|...). A small distance between two runs with different
   `method` or different input units is NOT evidence that the spaces are
   similar. `nearest()` reports the method and units of every neighbour for
   exactly this reason and, by default, refuses to rank across methods
   (`same_method=True`).

4. CHAIN-LEVEL RUNS HAVE NO BARS. Exact Mayer-Vietoris / chain-complex runs
   store Betti numbers but no filtration, hence no rows in `bar`. Those runs
   are excluded from `nearest()`, and the number excluded is reported.

Nothing here is a proof of anything: these are numerical comparisons of
truncated diagrams.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

from .api import TopoDB

try:  # gudhi is required for the distances, but not for importing the module
    import gudhi as _gudhi
except Exception:  # pragma: no cover - exercised only where gudhi is absent
    _gudhi = None

try:
    import gudhi.hera as _hera
except Exception:  # pragma: no cover
    _hera = None


CAVEAT = (
    "stored bars are truncated to the top 50 per dimension by api.add_bars, so these "
    "distances are approximate; infinite bars are dropped by default; filtration units "
    "differ between methods and are not commensurable"
)


def _require_gudhi():
    if _gudhi is None:
        raise RuntimeError("gudhi is not importable; the comparison helpers need it")


@dataclass
class Diagram:
    """A diagram as it was stored: finite points plus a count of what was removed."""

    run_id: int
    dim: int
    points: list[tuple[float, float]] = field(default_factory=list)
    n_infinite_dropped: int = 0
    n_stored: int = 0
    method: str | None = None
    dataset_id: str | None = None
    tier: str | None = None

    def as_array(self):
        return [[b, d] for b, d in self.points]


def load_diagram(db: TopoDB, run_id: int, dim: int, *, infinite: str = "drop",
                 cap_value: float | None = None) -> Diagram:
    """Read the stored bars of `run_id` in dimension `dim`.

    infinite='drop'  -- essential bars are discarded (default; see module docstring)
    infinite='cap'   -- essential bars die at `cap_value`, or at the max finite death
    infinite='keep'  -- essential bars are kept with death=inf (distances become inf)
    """
    if infinite not in {"drop", "cap", "keep"}:
        raise ValueError("infinite must be 'drop', 'cap' or 'keep'")
    meta = db.con.execute(
        "SELECT r.method, r.dataset_id, r.tier FROM run r WHERE r.id=?", (run_id,)).fetchone()
    if meta is None:
        raise KeyError(f"no run {run_id} in {db.path}")
    rows = db.con.execute(
        "SELECT birth, death FROM bar WHERE run_id=? AND dim=? ORDER BY rank", (run_id, dim)).fetchall()
    finite = [(float(r["birth"]), float(r["death"])) for r in rows if r["death"] is not None]
    infinites = [float(r["birth"]) for r in rows if r["death"] is None]
    points = list(finite)
    if infinite == "cap" and infinites:
        cap = cap_value
        if cap is None:
            cap = max((d for _, d in finite), default=None)
        if cap is None:
            raise ValueError(
                f"run {run_id} dim {dim}: infinite='cap' needs a cap_value; the run has no finite bar")
        points += [(b, float(cap)) for b in infinites]
    elif infinite == "keep":
        points += [(b, math.inf) for b in infinites]
    return Diagram(run_id=run_id, dim=dim, points=points,
                   n_infinite_dropped=len(infinites) if infinite == "drop" else 0,
                   n_stored=len(rows), method=meta["method"], dataset_id=meta["dataset_id"],
                   tier=meta["tier"])


def bottleneck(run_a: int, run_b: int, dim: int, *, db: TopoDB | None = None,
               infinite: str = "drop", cap_value: float | None = None, e: float | None = None):
    """Bottleneck distance between the STORED diagrams of two runs in one dimension.

    Returns (distance, info) where info records what was dropped. The distance is
    approximate: see the module docstring -- the stored bars are the top 50 per
    dimension only, and near-diagonal bars that were dropped at ingestion can move
    the true bottleneck distance by up to half their persistence.
    """
    _require_gudhi()
    db = db or TopoDB()
    da = load_diagram(db, run_a, dim, infinite=infinite, cap_value=cap_value)
    dbb = load_diagram(db, run_b, dim, infinite=infinite, cap_value=cap_value)
    if e is None:
        dist = _gudhi.bottleneck_distance(da.as_array(), dbb.as_array())
    else:
        dist = _gudhi.bottleneck_distance(da.as_array(), dbb.as_array(), e)
    return float(dist), _info(da, dbb, "bottleneck")


def wasserstein(run_a: int, run_b: int, dim: int, *, db: TopoDB | None = None, order: float = 1.0,
                internal_p: float = math.inf, infinite: str = "drop", cap_value: float | None = None):
    """q-Wasserstein distance (via gudhi.hera) between the STORED diagrams.

    Same approximation caveats as `bottleneck`, and worse: Wasserstein sums over
    every point, so discarding the bars below rank 50 removes terms from the sum
    outright. The result is a similarity heuristic, not a measured invariant.
    """
    if _hera is None:
        raise RuntimeError("gudhi.hera is not importable; wasserstein() needs it")
    import numpy as np

    db = db or TopoDB()
    da = load_diagram(db, run_a, dim, infinite=infinite, cap_value=cap_value)
    dbb = load_diagram(db, run_b, dim, infinite=infinite, cap_value=cap_value)
    A = np.array(da.as_array(), dtype=float).reshape(-1, 2)
    B = np.array(dbb.as_array(), dtype=float).reshape(-1, 2)
    dist = _hera.wasserstein_distance(A, B, order=order, internal_p=internal_p)
    info = _info(da, dbb, f"wasserstein(order={order})")
    return float(dist), info


def _info(da: Diagram, dbb: Diagram, metric: str) -> dict:
    return {
        "metric": metric,
        "dim": da.dim,
        "run_a": {"id": da.run_id, "method": da.method, "dataset": da.dataset_id,
                  "n_points_used": len(da.points), "n_bars_stored": da.n_stored,
                  "n_infinite_dropped": da.n_infinite_dropped},
        "run_b": {"id": dbb.run_id, "method": dbb.method, "dataset": dbb.dataset_id,
                  "n_points_used": len(dbb.points), "n_bars_stored": dbb.n_stored,
                  "n_infinite_dropped": dbb.n_infinite_dropped},
        "methods_comparable": da.method == dbb.method,
        "caveat": CAVEAT,
    }


def runs_with_bars(db: TopoDB, dim: int) -> list[int]:
    return [int(r[0]) for r in db.con.execute(
        "SELECT DISTINCT run_id FROM bar WHERE dim=? ORDER BY run_id", (dim,)).fetchall()]


def nearest(run_id: int, k: int = 5, *, dim: int = 1, db: TopoDB | None = None,
            metric: str = "bottleneck", same_method: bool = True, infinite: str = "drop",
            cap_value: float | None = None) -> dict:
    """The k runs whose STORED diagram in dimension `dim` is closest to `run_id`.

    This is a heuristic over lossy data. Every caveat in the module docstring
    applies, and two more are reported in the result:
      * `excluded_no_bars`: runs with no stored bars in this dimension (every
        exact chain-complex / Mayer-Vietoris run) are not comparable at all and
        are skipped, not ranked last;
      * `excluded_other_method`: with same_method=True (the default) runs built
        from a different filtration are skipped, because their filtration values
        are in different units and a small distance between them means nothing.
    """
    db = db or TopoDB()
    ref = load_diagram(db, run_id, dim, infinite=infinite, cap_value=cap_value)
    total_runs = db.con.execute("SELECT count(*) FROM run").fetchone()[0]
    candidates = [r for r in runs_with_bars(db, dim) if r != run_id]
    excluded_no_bars = total_runs - len(candidates) - 1
    fn = bottleneck if metric == "bottleneck" else wasserstein
    rows, skipped_method = [], 0
    for other in candidates:
        cand = load_diagram(db, other, dim, infinite=infinite, cap_value=cap_value)
        if same_method and cand.method != ref.method:
            skipped_method += 1
            continue
        if not ref.points and not cand.points:
            continue
        d, _ = fn(run_id, other, dim, db=db, infinite=infinite, cap_value=cap_value)
        rows.append({"run_id": other, "distance": d, "dataset": cand.dataset_id,
                     "method": cand.method, "tier": cand.tier,
                     "n_points_used": len(cand.points), "n_bars_stored": cand.n_stored})
    rows.sort(key=lambda r: (math.inf if math.isnan(r["distance"]) else r["distance"]))
    return {
        "run_id": run_id, "dim": dim, "metric": metric,
        "reference": {"dataset": ref.dataset_id, "method": ref.method, "tier": ref.tier,
                      "n_points_used": len(ref.points), "n_bars_stored": ref.n_stored,
                      "n_infinite_dropped": ref.n_infinite_dropped},
        "neighbours": rows[:k],
        "n_compared": len(rows),
        "excluded_no_bars": max(excluded_no_bars, 0),
        "excluded_other_method": skipped_method,
        "caveat": CAVEAT,
    }
