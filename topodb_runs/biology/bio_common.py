"""Shared persistent-homology helpers for the biology block.

Design rules (binding, from audit/tda_validation/genetics/report.json and LeanFlow CLAUDE.md):

  (a) NEVER embed a non-Euclidean distance into 3-D and then run alpha.  The earlier
      genetics validation produced a FALSE dominant loop (alpha on MDS-3 of the HA
      p-distance matrix, dominance 11.1) on a phylogenetically tree-like segment.
      For Hi-C, scRNA correlation/cosine and sequence Hamming metrics this module
      offers ONLY `rips_from_distance`.  There is no MDS function here on purpose.
  (b) Alpha is legitimate only where the point cloud is genuinely Euclidean
      (protein/RNA C-alpha / P atom coordinates in Angstrom).

  gudhi note that governs every result here: `RipsComplex.create_simplex_tree(
  max_dimension=d)` builds the d-skeleton, which yields homology only up to d-1.
  `_diagrams` therefore always passes max_dimension = max_hom_dim + 1.  The Step 0
  known-answer control runs through these same functions, so a regression in this
  convention is caught there rather than silently reported as b1 = 0.
"""
from __future__ import annotations

import hashlib
import json
import os
import resource
import subprocess
import time
from pathlib import Path

# three agents share 8 cores: keep BLAS single-threaded
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np  # noqa: E402
import gudhi  # noqa: E402

BIO_DATA = Path("/mnt/disks/disk-socrateai-local-1/dualscale-data-r3/tda_validation/biology")
GEN_DATA = Path("/mnt/disks/disk-socrateai-local-1/dualscale-data-r3/tda_validation/genetics")
WT = Path("/mnt/disks/disk-socrateai-local-1/wt-topo-bio")
RESULTS = WT / "topodb_runs" / "biology" / "results"


# ------------------------------------------------------------------ persistence
def _diagrams(st, max_hom_dim: int, coeff: int = 2) -> dict[int, list[tuple[float, float]]]:
    st.compute_persistence(homology_coeff_field=coeff, persistence_dim_max=False)
    out = {d: [] for d in range(max_hom_dim + 1)}
    for d in range(max_hom_dim + 1):
        for b, dth in st.persistence_intervals_in_dimension(d):
            out[d].append((float(b), float(dth)))
    return out


def rips_from_points(points, max_hom_dim: int = 1, max_edge: float | None = None,
                     coeff: int = 2) -> dict[int, list[tuple[float, float]]]:
    """Vietoris-Rips on Euclidean coordinates."""
    pts = np.asarray(points, dtype=float)
    kw = {"points": pts}
    if max_edge is not None:
        kw["max_edge_length"] = float(max_edge)
    rc = gudhi.RipsComplex(**kw)
    st = rc.create_simplex_tree(max_dimension=max_hom_dim + 1)  # d-skeleton -> H_{d-1}
    return _diagrams(st, max_hom_dim, coeff)


def rips_from_distance(dmat, max_hom_dim: int = 1, max_edge: float | None = None,
                       coeff: int = 2) -> dict[int, list[tuple[float, float]]]:
    """Vietoris-Rips on a GIVEN metric (lesson (a): no embedding step)."""
    d = np.asarray(dmat, dtype=float)
    kw = {"distance_matrix": d}
    if max_edge is not None:
        kw["max_edge_length"] = float(max_edge)
    rc = gudhi.RipsComplex(**kw)
    st = rc.create_simplex_tree(max_dimension=max_hom_dim + 1)
    return _diagrams(st, max_hom_dim, coeff)


def alpha_from_points(points, max_hom_dim: int = 2, coeff: int = 2):
    """Alpha complex; filtration values are converted from squared to plain radii.

    Legitimate ONLY for genuinely Euclidean coordinates.
    """
    pts = np.asarray(points, dtype=float)
    ac = gudhi.AlphaComplex(points=pts)
    st = ac.create_simplex_tree()
    dg = _diagrams(st, max_hom_dim, coeff)
    return {d: [(float(np.sqrt(b)), float(np.sqrt(x)) if np.isfinite(x) else float("inf"))
                for b, x in bars] for d, bars in dg.items()}


# ------------------------------------------------------------------- statistics
def finite(bars):
    return [(b, d) for b, d in bars if np.isfinite(d)]


def pers(bars):
    return sorted((d - b for b, d in finite(bars)), reverse=True)


def betti_from_bars(bars) -> int:
    """Number of infinite bars (H0: connected components; H1+ normally 0)."""
    return int(sum(1 for _, d in bars if not np.isfinite(d)))


def dominance(bars) -> float:
    """P1/P2 of the finite bars; 0 if none, inf if exactly one."""
    p = pers(bars)
    if not p:
        return 0.0
    if len(p) == 1:
        return float("inf")
    return float(p[0] / p[1]) if p[1] > 0 else float("inf")


def S_stat(bars) -> float:
    """Longest finite bar as a fraction of total finite persistence (0..1)."""
    p = pers(bars)
    tot = float(sum(p))
    return float(p[0] / tot) if tot > 0 else 0.0


def top_bars(bars, k: int = 5):
    return sorted(finite(bars), key=lambda bd: bd[1] - bd[0], reverse=True)[:k]


def max_pers(bars) -> float:
    p = pers(bars)
    return float(p[0]) if p else 0.0


def rank_p(observed: float, null_values, higher_is_extreme: bool = True) -> float:
    """Empirical p with the +1 correction; never returns 0."""
    nv = np.asarray(list(null_values), dtype=float)
    n = nv.size
    if n == 0:
        raise ValueError("empty null")
    ge = int(np.sum(nv >= observed)) if higher_is_extreme else int(np.sum(nv <= observed))
    return float((ge + 1) / (n + 1))


def auc(pos, neg) -> float:
    """Mann-Whitney U / (n_pos n_neg) -- ties count as 1/2."""
    pos = np.asarray(pos, dtype=float)
    neg = np.asarray(neg, dtype=float)
    if pos.size == 0 or neg.size == 0:
        return float("nan")
    gt = (pos[:, None] > neg[None, :]).sum()
    eq = (pos[:, None] == neg[None, :]).sum()
    return float((gt + 0.5 * eq) / (pos.size * neg.size))


# ------------------------------------------------------------------- bookkeeping
def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def peak_mb() -> float:
    """Peak RSS of this process, MB (Linux ru_maxrss is KB)."""
    return round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0, 1)


def git_head(repo=WT) -> str | None:
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], cwd=str(repo),
                              capture_output=True, text=True, timeout=10).stdout.strip() or None
    except Exception:
        return None


class Block:
    """Accumulates dataset/run records for one compute script and writes ONE json.

    Nothing here touches the SQLite file: ingest.py opens the DB briefly, writes and
    closes, so a long computation never holds a lock against the other two agents.
    """

    def __init__(self, name: str, script: str, command: str):
        self.name = name
        self.script = script
        self.command = command
        self.t0 = time.time()
        self.datasets: list[dict] = []
        self.runs: list[dict] = []

    def dataset(self, **kw):
        self.datasets.append(kw)
        return kw["id"]

    def run(self, *, dataset_id, method, coeff_field, params, tier, diagrams=None,
            betti=None, expected_betti=None, bars=None, stats=None, controls=None,
            findings=None, max_dim=None, preprocessing=None, seed=None,
            wall_sec=None, note=None):
        rec = {
            "dataset_id": dataset_id, "method": method, "coeff_field": coeff_field,
            "params": params, "tier": tier, "max_dim": max_dim,
            "preprocessing": preprocessing, "seed": str(seed) if seed is not None else None,
            "script": self.script, "command": self.command,
            "wall_sec": wall_sec, "peak_mb": peak_mb(),
            "betti": betti or {}, "expected_betti": expected_betti,
            "bars": bars or {}, "stats": stats or [], "controls": controls or [],
            "findings": findings or [], "note": note,
        }
        if diagrams is not None:
            rec["bars"] = {str(d): [[float(b), (None if not np.isfinite(x) else float(x))]
                                    for b, x in bs] for d, bs in diagrams.items()}
        self.runs.append(rec)
        return rec

    def write(self):
        RESULTS.mkdir(parents=True, exist_ok=True)
        out = RESULTS / f"{self.name}.json"
        payload = {"block": self.name, "script": self.script, "command": self.command,
                   "git_commit": git_head(), "wall_sec_total": round(time.time() - self.t0, 2),
                   "peak_mb": peak_mb(), "datasets": self.datasets, "runs": self.runs}
        out.write_text(json.dumps(payload, indent=1, default=str))
        print(f"[block] wrote {out}  datasets={len(self.datasets)} runs={len(self.runs)} "
              f"wall={payload['wall_sec_total']}s peak={payload['peak_mb']}MB")
        return out
