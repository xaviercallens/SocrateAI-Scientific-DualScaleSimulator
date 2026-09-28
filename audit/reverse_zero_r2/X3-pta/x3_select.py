"""Data-blind pulsar selection for X3. Uses only tim MJDs. Writes selection.json.
Run: cd audit/reverse_zero_r2/X3-pta && <venv-pta python> x3_select.py
"""
import json
import numpy as np
from x3_common import *

names = all_names()
rows, dropped = [], []
for n in names:
    par, tim = find_files(n)
    m = np.array(tim_mjds(tim))
    span_yr = (m.max() - m.min()) / 365.25
    rec = dict(name=n, ntoa=int(len(m)), mjd_min=float(m.min()), mjd_max=float(m.max()), span_yr=float(span_yr))
    if is_split_duplicate(n, names):
        dropped.append(dict(rec, reason="per-telescope split of a combined entry"))
    elif span_yr < MIN_SPAN_YR:
        dropped.append(dict(rec, reason=f"span < {MIN_SPAN_YR} yr"))
    else:
        rows.append(rec)
tmin = min(r["mjd_min"] for r in rows)
tmax = max(r["mjd_max"] for r in rows)
out = dict(criterion="drop per-telescope split duplicates (ao/gbt when combined exists), then span >= 3 yr; uses tim MJDs only",
           n_selected=len(rows), tspan_common_s=(tmax - tmin) * 86400.0, mjd_min=tmin, mjd_max=tmax,
           selected=rows, dropped=dropped)
(HERE / "selection.json").write_text(json.dumps(out, indent=1))
print(len(names), "in release;", len(rows), "selected; dropped:", [(d["name"], d["reason"]) for d in dropped])
print("Tspan yr", (tmax - tmin) / 365.25)
