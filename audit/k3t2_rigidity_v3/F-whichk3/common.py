"""Shared helpers for Track F. Paths are found from the script location (no absolute paths)."""
import json
from pathlib import Path
from fractions import Fraction

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]          # repo root (audit/k3t2_rigidity_v3/F-whichk3 -> 3 levels up)
RELDIR = HERE.relative_to(ROOT).as_posix()


def _load(name, default):
    p = HERE / name
    if p.exists():
        return json.loads(p.read_text())
    return default


def _save(name, obj):
    (HERE / name).write_text(json.dumps(obj, indent=1, sort_keys=False, default=str))


def add_inputs(items):
    d = _load("inputs.json", {})
    for it in items:
        d[it["name"]] = it
    _save("inputs.json", d)


def add_results(items):
    d = _load("results.json", {"results": {}, "rigidity": {}, "could_not_do": {}})
    for it in items:
        d["results"][it["id"]] = it
    _save("results.json", d)


def add_rigidity(items):
    d = _load("results.json", {"results": {}, "rigidity": {}, "could_not_do": {}})
    for it in items:
        d["rigidity"][it["id"]] = it
    _save("results.json", d)


def add_cnd(key, text):
    d = _load("results.json", {"results": {}, "rigidity": {}, "could_not_do": {}})
    d["could_not_do"][key] = text
    _save("results.json", d)


def add_exports(obj):
    d = _load("exports.json", {})
    d.update(obj)
    _save("exports.json", d)


def fr(x):
    return str(x) if isinstance(x, Fraction) else x


def cmd(script, args=""):
    return f"cd {RELDIR} && /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python {script} {args}".strip()


# ---------- integer linear algebra ----------
def int_kernel(A, ncols):
    """Z-basis of {x in Z^ncols : A x = 0}, A a list of rows of ints (saturated)."""
    # rows of [A^T | I]; integer row reduction on the A^T part
    m = len(A)
    rows = [[A[i][j] for i in range(m)] + [1 if k == j else 0 for k in range(ncols)] for j in range(ncols)]
    r = 0
    for c in range(m):
        while True:
            nz = [i for i in range(r, ncols) if rows[i][c] != 0]
            if not nz:
                break
            p = min(nz, key=lambda i: abs(rows[i][c]))
            rows[r], rows[p] = rows[p], rows[r]
            done = True
            for i in range(r + 1, ncols):
                if rows[i][c] != 0:
                    q = rows[i][c] // rows[r][c]
                    rows[i] = [a - q * b for a, b in zip(rows[i], rows[r])]
                    if rows[i][c] != 0:
                        done = False
            if done:
                r += 1
                break
    ker = [row[m:] for row in rows if all(v == 0 for v in row[:m])]
    return ker


def det_int(M):
    import sympy
    return int(sympy.Matrix(M).det())


def lagrange_reduce(p, q, r):
    """Reduce positive-definite binary form [[p,q],[q,r]] (Gram) to canonical: 0<=2|q|<=p<=r, q>=0."""
    assert p * r - q * q > 0 and p > 0
    while True:
        if r < p:
            p, r = r, p
            continue
        # translate: second vector v2 -> v2 - k v1
        k = round(Fraction(q, p))
        if k != 0:
            r = r - 2 * k * q + k * k * p
            q = q - k * p
            continue
        if r < p:
            continue
        break
    if q < 0:
        q = -q
    return (p, q, r)
