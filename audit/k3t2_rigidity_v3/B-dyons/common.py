"""Shared helpers for Track B v3 (paths relative to repo root, exact loaders)."""
import json
from fractions import Fraction as Fr
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]          # <repo>/audit/k3t2_rigidity_v3/B-dyons/common.py -> <repo>


def pf(s):
    if s is None:
        return None
    if isinstance(s, (int,)):
        return Fr(s)
    if "/" in s:
        a, b = s.split("/")
        return Fr(int(a), int(b))
    return Fr(int(s))


def load_inputs():
    with open(HERE / "inputs.json") as f:
        L = json.load(f)
    return {d["name"]: d for d in L}


def load_cache(qmax):
    with open(HERE / f"theta_forms_cache_Q{qmax}.json") as f:
        return json.load(f)


def load2d(cache, key):
    out = {}
    for k, v in cache[key].items():
        n, l = k.split(",")
        out[(int(n), int(l))] = pf(v)
    return out


def load1d(cache, key):
    return {int(k): pf(v) for k, v in cache[key].items()}


def load_cB(cache):
    return {int(k): pf(v) for k, v in cache["cB_table"].items()}


def find_trackD_exports():
    """Prefer a v3 Track D exports.json; fall back to v2 (recorded in the return value)."""
    for tag in ("k3t2_rigidity_v3", "k3t2_rigidity_v2"):
        p = REPO / "audit" / tag / "D-tda" / "exports.json"
        if p.exists():
            with open(p) as f:
                return tag, str(p.relative_to(REPO)), json.load(f)
    return None, None, None
