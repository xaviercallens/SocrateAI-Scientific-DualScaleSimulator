"""Content forensics of dual_scale_datasets/domain06_supercon_abrikosov_vortices.npz
(manifest source string: 'NIMS SuperCon Database & Scanning SQUID Microscopy
Benchmark (YBCO)'). Decides from content whether this is a measurement or a
synthetic lattice. NOT used as validation data under any outcome.

Checks:
 1. keys, shapes, dtypes; z column constant?
 2. fit to an ideal triangular lattice with spacing 1 (row pitch sqrt3/2,
    alternate-row offset 0.5): integer site indices, uniqueness, storage order
 3. residual distribution (mean, std per axis, D'Agostino normality p-values)
 4. alpha-path signature (same as a jittered synthetic lattice?) using the
    pipeline under test
 5. sibling files: same raster / exact-lattice fingerprint in other 'domain'
    files with (N,3) coordinate arrays.
Command: prlimit --as=8589934592 -- .venv-tda/bin/python assess_domain06.py
Output: ../results/domain06_assessment.json
"""
import glob
import hashlib
import json
import os
import sys

import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qf_common as qc  # noqa: E402

D = "/mnt/disks/disk-socrateai-local-1/dual_scale_datasets"
F = os.path.join(D, "domain06_supercon_abrikosov_vortices.npz")


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def main():
    out = {"file": F, "sha256": sha(F), "size_bytes": os.path.getsize(F)}
    man = json.load(open(os.path.join(D, "dual_scale_datasets_manifest.json")))
    out["manifest_entry"] = man.get("domain_06_superconductor_vortices")
    d = np.load(F)
    out["keys"] = {k: {"shape": list(d[k].shape), "dtype": str(d[k].dtype)} for k in d.files}
    c = d["coords"]
    out["pinning_strength_value"] = float(d["pinning_strength"])
    out["z_column_unique_values"] = np.unique(c[:, 2]).tolist()
    xy = c[:, :2]
    j = np.round(xy[:, 1] / (np.sqrt(3) / 2)).astype(int)
    i = np.round(xy[:, 0] - 0.5 * (j % 2)).astype(int)
    ideal = np.stack([i + 0.5 * (j % 2), j * np.sqrt(3) / 2], 1)
    res = xy - ideal
    out["lattice_fit"] = {
        "i_range": [int(i.min()), int(i.max())], "j_range": [int(j.min()), int(j.max())],
        "n_unique_sites": len(set(zip(i.tolist(), j.tolist()))), "n_points": int(len(xy)),
        "storage_order_is_i_major_j_minor_raster": bool(np.all(i * 15 + j == np.arange(len(xy)))),
        "residual_mean": res.mean(0).tolist(), "residual_std": res.std(0).tolist(),
        "max_abs_residual": float(np.abs(res).max()),
        "normaltest_p": [float(stats.normaltest(res[:, 0]).pvalue), float(stats.normaltest(res[:, 1]).pvalue)],
        "units": "none: spacing is exactly 1.0 (dimensionless); a real SQUID/STM/decoration image has lengths in nm or um set by B",
    }
    web = qc.web()
    info, by_dim = web.alpha_persistence(xy, 9.0, "domain06", None)
    d0 = qc.h0_finite_deaths(by_dim)
    h1 = np.array(by_dim[1])
    out["alpha_signature"] = {"h0": qc.spread_stats(d0), "h1_death_median": float(np.median(h1[:, 1])),
                              "sqrt3_median_h1_death": float(np.sqrt(3) * np.median(h1[:, 1]))}
    sib = {}
    for p in sorted(glob.glob(os.path.join(D, "domain*.npz"))):
        z = np.load(p)
        sib[os.path.basename(p)] = {k: [list(z[k].shape), str(z[k].dtype)] for k in z.files}
    out["sibling_files_keys"] = sib
    out["verdict"] = None  # filled in report.json from these facts
    with open(os.path.join(qc.RESULTS, "domain06_assessment.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
