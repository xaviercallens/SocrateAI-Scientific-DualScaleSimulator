"""Content-based assessment of /mnt/disks/disk-socrateai-local-1/bio_datasets/*.npy and
/mnt/disks/disk-socrateai-local-1/dual_scale_datasets/*.npz (+ manifest): real measurements or
synthetic arrays under real-source names? NOT used as validation data.
Command: python assess_local_datasets.py   (no randomness)"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import glob
import json
import sys
import urllib.request

import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tda_common import sha256, dump

HERE = os.path.dirname(os.path.abspath(__file__))
BIO = "/mnt/disks/disk-socrateai-local-1/bio_datasets"
DS = "/mnt/disks/disk-socrateai-local-1/dual_scale_datasets"
REF = "/mnt/disks/disk-socrateai-local-1/dualscale-data-r3/tda_validation/genetics/rcsb_reference"


def describe(a):
    d = {"shape": list(a.shape), "dtype": str(a.dtype)}
    if a.dtype.kind in "fiub" and a.size:
        f = a.astype(float).ravel()
        fin = f[np.isfinite(f)]
        d.update({"min": float(fin.min()), "max": float(fin.max()), "mean": float(fin.mean()), "std": float(fin.std()),
                  "frac_integer_valued": float(np.mean(fin == np.round(fin))),
                  "frac_exact_zero": float(np.mean(fin == 0)),
                  "n_unique": int(np.unique(fin).size),
                  "skew": float(stats.skew(fin)), "excess_kurtosis": float(stats.kurtosis(fin))})
        if fin.size > 50:
            z = (fin - fin.mean()) / (fin.std() + 1e-300)
            d["ks_vs_normal_stat"] = float(stats.kstest(z, "norm").statistic)
        if a.ndim == 2 and a.shape[0] == a.shape[1]:
            d["max_abs_asymmetry"] = float(np.nanmax(np.abs(a - a.T)))
            d["numerical_rank"] = int(np.linalg.matrix_rank(np.nan_to_num(a.astype(float))))
            dg = [float(np.nanmean(np.diagonal(a, s))) for s in (1, 2, 5, 10, 50) if s < a.shape[0]]
            d["mean_diagonal_offsets_1_2_5_10_50"] = dg
        if a.ndim == 2 and a.shape[1] in (2, 3, 4) and a.shape[0] > 10:
            c = a.astype(float) - a.astype(float).mean(0)
            r = np.linalg.norm(c[:, :3], axis=1)
            d["radius_from_centroid_cv"] = float(r.std() / r.mean()) if r.mean() > 0 else None
            steps = np.linalg.norm(np.diff(a[:, :3].astype(float), axis=0), axis=1)
            d["consecutive_step_mean_cv"] = [float(steps.mean()), float(steps.std() / steps.mean()) if steps.mean() > 0 else None]
    return d


def ca_from_pdb(path):
    xyz = []
    with open(path) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA" and line[16] in (" ", "A"):
                xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.array(xyz)


def specific_checks():
    from scipy.spatial.distance import pdist, squareform
    c = {}
    C = np.load(os.path.join(BIO, "01_hic_contact_map.npy"))
    X = np.load(os.path.join(BIO, "01_hic_chromatin_coords.npy"))
    D = squareform(pdist(X))
    iu = np.triu_indices(len(C), 1)
    lc, ld = np.log(C[iu]), np.log(D[iu])
    slope, icpt = np.polyfit(ld, lc, 1)
    resid = lc - (icpt + slope * ld)
    c["hic_contact_vs_coords"] = {"loglog_slope": float(slope), "loglog_intercept": float(icpt),
                                  "spearman_contact_vs_distance": float(stats.spearmanr(C[iu], D[iu]).correlation),
                                  "resid_std_loglog": float(resid.std()),
                                  "max_contact": float(C.max()), "log10_max_contact": float(np.log10(C.max())),
                                  "frac_offdiag_zero": float(np.mean(C[iu] == 0)),
                                  "diag_values_unique": np.unique(np.round(np.diagonal(C), 12)).tolist()[:5]}
    E = np.load(os.path.join(BIO, "02_scrna_expression.npy"))
    c["scrna_expression"] = {"frac_zero": float(np.mean(E == 0)), "frac_integer": float(np.mean(E == np.round(E))),
                             "note": "real scRNA-seq count/normalised matrices are dominated by exact zeros (dropout)"}
    B = np.load(os.path.join(BIO, "05_dna_methylation_beta.npy"))
    h, edges = np.histogram(B, bins=10, range=(0, 1))
    c["methylation_beta_histogram_10bins"] = (h / h.sum()).tolist()
    V = np.load(os.path.join(BIO, "04_spatial_visium_coords.npy"))
    c["visium_coords"] = {"n_unique_x": int(np.unique(V[:, 0]).size), "n_unique_y": int(np.unique(V[:, 1]).size),
                          "x_spacing_unique": np.unique(np.round(np.diff(np.unique(V[:, 0])), 9)).tolist()}
    S = np.load(os.path.join(BIO, "10_metabolic_stoichiometry.npy"))
    c["stoichiometry"] = {"nnz_per_column_unique": np.unique((S != 0).sum(0)).tolist(),
                          "column_sums_unique": np.unique(S.sum(0)).tolist()}
    F = np.load(os.path.join(BIO, "06_cell_features.npy"))
    c["cell_features_col_mean_std"] = [F.mean(0).round(4).tolist(), F.std(0).round(4).tolist()]
    z = np.load(os.path.join(DS, "domain08_earthscope_seismic_rupture.npz"))
    m = z["magnitudes"]
    c["seismic_magnitudes"] = {"n_ge_10": int(np.sum(m >= 10)), "max": float(m.max()),
                               "note": "moment magnitudes above ~9.5 have never been recorded"}
    z = np.load(os.path.join(DS, "domain06_supercon_abrikosov_vortices.npz"))
    q = z["coords"]
    c["vortex_coords_unique_per_axis"] = [int(np.unique(np.round(q[:, k], 9)).size) for k in range(3)]
    return c


def main():
    out = {"bio_datasets": {}, "dual_scale_datasets": {}, "notes": []}
    files = sorted(glob.glob(os.path.join(BIO, "*")))
    mt = {os.path.basename(p): os.path.getmtime(p) for p in files}
    out["bio_mtime_span_sec"] = float(max(mt.values()) - min(mt.values()))
    os.makedirs(REF, exist_ok=True)
    for p in files:
        name = os.path.basename(p)
        rec = {"sha256": sha256(p), "size_bytes": os.path.getsize(p)}
        if name.endswith(".npy"):
            a = np.load(p, allow_pickle=False)
            rec.update(describe(a))
        elif name.endswith(".pdb"):
            pid = name.split("_")[1].split(".")[0]
            url = f"https://files.rcsb.org/download/{pid}.pdb"
            ref = os.path.join(REF, f"{pid}.pdb")
            try:
                if not os.path.exists(ref):
                    urllib.request.urlretrieve(url, ref)
                rec["rcsb_url"] = url
                rec["rcsb_sha256"] = sha256(ref)
                rec["identical_to_rcsb_download"] = rec["rcsb_sha256"] == rec["sha256"]
                la = [l for l in open(p) if l.startswith(("ATOM", "HETATM"))]
                lb = [l for l in open(ref) if l.startswith(("ATOM", "HETATM"))]
                rec["atom_records_identical_to_rcsb"] = la == lb
                rec["n_atom_records"] = len(la)
            except Exception as e:  # network failure is recorded, not hidden
                rec["rcsb_error"] = repr(e)
            ca = ca_from_pdb(p)
            npy = os.path.join(BIO, f"03_{pid}_ca_coords.npy")
            if os.path.exists(npy):
                b = np.load(npy)
                rec["ca_npy_shape"] = list(b.shape)
                rec["n_CA_in_pdb_altloc_A_or_blank"] = int(len(ca))
                if b.shape == ca.shape:
                    rec["ca_npy_max_abs_diff_vs_pdb"] = float(np.max(np.abs(b - ca)))
        out["bio_datasets"][name] = rec

    man = json.load(open(os.path.join(DS, "dual_scale_datasets_manifest.json")))
    out["dual_scale_manifest"] = {k: {"source": v["source"], "shape": v["shape"]} for k, v in man.items()}
    for p in sorted(glob.glob(os.path.join(DS, "*.npz"))):
        z = np.load(p, allow_pickle=False)
        rec = {"sha256": sha256(p), "size_bytes": os.path.getsize(p), "arrays": {}}
        for k in z.files:
            rec["arrays"][k] = describe(z[k])
        out["dual_scale_datasets"][os.path.basename(p)] = rec
    out["specific_checks"] = specific_checks()
    print(json.dumps(out["specific_checks"], indent=1))
    dump(out, os.path.join(HERE, "local_datasets_assessment.json"))
    for k, v in out["bio_datasets"].items():
        print(k, {kk: v.get(kk) for kk in ("shape", "dtype", "min", "max", "mean", "std", "frac_integer_valued", "n_unique", "skew",
                                           "excess_kurtosis", "ks_vs_normal_stat", "max_abs_asymmetry", "numerical_rank",
                                           "identical_to_rcsb_download", "atom_records_identical_to_rcsb", "ca_npy_max_abs_diff_vs_pdb",
                                           "radius_from_centroid_cv", "consecutive_step_mean_cv", "mean_diagonal_offsets_1_2_5_10_50") if kk in v})
    for k, v in out["dual_scale_datasets"].items():
        print(k)
        for a, d in v["arrays"].items():
            print("   ", a, {kk: (round(x, 4) if isinstance(x, float) else x) for kk, x in d.items() if kk in ("shape", "min", "max", "mean", "std", "frac_integer_valued", "n_unique", "ks_vs_normal_stat", "radius_from_centroid_cv", "consecutive_step_mean_cv", "max_abs_asymmetry", "numerical_rank")})
    print("mtime span", out["bio_mtime_span_sec"])


if __name__ == "__main__":
    main()
