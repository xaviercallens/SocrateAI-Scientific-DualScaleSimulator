"""Shared utilities for Track E v3 (flux vacua on K3 x T2/Z2, Tripathy-Trivedi hep-th/0301139).

Everything is relative to the repository root, found from this file's location:
  <root>/audit/k3t2_rigidity_v3/E-flux/scripts/common.py  -> parents[4] == <root>
No absolute paths appear in any script.

Exact arithmetic only: python ints / fractions.Fraction / sympy.
"""
from __future__ import annotations
import json, re, hashlib
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
EDIR = Path(__file__).resolve().parents[1]
RES = EDIR / "results"
RES.mkdir(exist_ok=True)
SRC = ROOT / "audit" / "k3t2_rigidity_v2" / "sources"
TT_TXT = SRC / "hep-th_0301139_TripathyTrivedi.txt"
BL_TXT = SRC / "hep-th_9605184_BLPSSW.txt"

def rel(p: Path) -> str:
    return str(Path(p).resolve().relative_to(ROOT))

def sha256(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

_TT_LINES = None
def tt_lines():
    global _TT_LINES
    if _TT_LINES is None:
        _TT_LINES = TT_TXT.read_text(encoding="utf-8").split("\n")
    return _TT_LINES

def find_lines(snippet: str, lines=None, first=False):
    """1-based line numbers of lines containing snippet (whitespace-normalised)."""
    L = lines if lines is not None else tt_lines()
    norm = lambda s: re.sub(r"\s+", " ", s)
    out = [i + 1 for i, l in enumerate(L) if norm(snippet) in norm(l)]
    if first:
        assert out, f"snippet not found: {snippet!r}"
        return out[0]
    return out

def _joined(lines):
    """whitespace-normalised join of all lines with a char->line map"""
    parts, cmap = [], []
    for i, l in enumerate(lines):
        t = re.sub(r"\s+", " ", l).strip()
        if not t:
            continue
        parts.append(t)
        cmap.extend([i + 1] * (len(t) + 1))
    return " ".join(parts), cmap

def quote(snippet: str, lines=None):
    """Return {line, text} : first line containing snippet; if the snippet wraps across lines, the line where it STARTS
    and the wrapped text (asserts existence)."""
    L = lines if lines is not None else tt_lines()
    out = find_lines(snippet, L)
    if out:
        return {"line": out[0], "text": re.sub(r"\s+", " ", L[out[0] - 1]).strip()}
    joined, cmap = _joined(L)
    sn = re.sub(r"\s+", " ", snippet).strip()
    k = joined.find(sn)
    assert k >= 0, f"snippet not found: {snippet!r}"
    start = cmap[k]; end = cmap[min(k + len(sn) - 1, len(cmap) - 1)]
    txt = " ".join(re.sub(r"\s+", " ", L[i]).strip() for i in range(start - 1, end))
    return {"line": start, "lines": f"{start}-{end}", "text": txt}

def int_rows(start_line: int, end_line: int, width: int):
    """Parse integer rows (with the unicode minus) from TT text lines in [start,end] (1-based)."""
    rows = []
    for l in tt_lines()[start_line - 1:end_line]:
        s = re.sub(r"[^0-9\-−\s]", " ", l).replace("−", "-")
        toks = s.split()
        if len(toks) == width and all(re.fullmatch(r"-?\d+", t) for t in toks):
            rows.append([int(t) for t in toks])
    return rows

def parse_H33():
    a = find_lines("where the matrix H3,3 is defined as", first=True)
    b = find_lines("and E8 is the Catran matrix", first=True)
    rows = int_rows(a, b, 6)
    assert len(rows) == 6, rows
    return rows

def parse_E8():
    a = find_lines("and E8 is the Catran matrix", first=True)
    b = find_lines("The basis vectors (e1", first=True)
    rows = int_rows(a, b, 8)
    assert len(rows) == 8, rows
    return rows

def dotm(u, H, v):
    return sum(u[i] * H[i][j] * v[j] for i in range(len(u)) for j in range(len(v)) if H[i][j])

def write_json(name: str, obj):
    p = RES / name
    with open(p, "w") as f:
        json.dump(obj, f, indent=2, sort_keys=False, default=str)
    return p

def rel_command(script: str, args: str = "") -> str:
    d = rel(EDIR / "scripts")
    return (f"cd {d} && /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python {script} {args}").strip()
