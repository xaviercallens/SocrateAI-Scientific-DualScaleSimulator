#!/usr/bin/env python3
"""CAMELS IllustrisTNG LH neutral-hydrogen maps: cubical persistence vs cosmology.

DATA.  `params_2/test_LH.npy` = 750 maps of 256x256 float32, with
`params_2/test_labels_LH_2.npy` giving (Omega_m, sigma_8) per map.  The raw
source is the CAMELS Multifield Dataset map `Maps_HI_IllustrisTNG_LH_z=0.00.npy`;
`process_hi_maps.py` (shipped beside the arrays, sha256 recorded) applied
  renorm(log1p(minmax(HI)))
with GLOBAL constants over the whole array.  That is one monotone map applied
identically to every map, so (a) the sublevel topology is exactly the topology of
the raw HI field and (b) persistence VALUES remain comparable between maps.  No
further transform is applied here -- see AMENDMENT A2.

EFFECTIVE N, checked before any correlation was computed (this is the check that
can invalidate the headline result):
  750 maps carry only 555 DISTINCT label rows -- 397 singletons, 126 pairs,
  27 triples, 5 quadruples.  Maps sharing a label row come from the same
  simulation, so a per-map permutation null would be anti-conservative.
  293 label rows also appear in the val split, i.e. the splits share simulations.
THE FIX USED HERE: sample one map per SIMULATION, so every point in the
correlation is an independent simulation and N_eff = N_maps = 200 exactly.

Tier X (numerics).
"""
from __future__ import annotations

import collections
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (Timer, add_run, peak_mb, save_json, sha256_file, topodb)  # noqa: E402

SCRIPT = "topodb_runs/astro/run_camels.py"
PY = "/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python"
COMMAND = f"prlimit --as=8589934592 -- {PY} {SCRIPT}"
SEED = 20260920
B = "/mnt/disks/disk-socrateai-local-1/dualscale-data-r3/camels_hi_hf/"
MAPS = B + "params_2/test_LH.npy"
LAB2 = B + "params_2/test_labels_LH_2.npy"
LAB6 = B + "params_6/test_labels_LH.npy"
DS = "astro/camels_hi_illustristng_lh"

N_SAMPLE = 200
N_PERM = 10000
STATS = ["n_h0_bars", "n_h1_bars", "total_persistence_h0",
         "total_persistence_h1", "max_persistence_h1", "mean_persistence_h1"]
PARAMS = ["Omega_m", "sigma_8"]
N_TESTS = len(STATS) * len(PARAMS)            # 12, declared in advance
BONF = 0.05 / N_TESTS


def map_stats(img):
    """The six pre-declared per-map statistics, sublevel cubical persistence."""
    import gudhi
    cc = gudhi.CubicalComplex(top_dimensional_cells=np.ascontiguousarray(img, dtype=float))
    cc.compute_persistence(homology_coeff_field=2)
    d0 = np.array(cc.persistence_intervals_in_dimension(0), float).reshape(-1, 2)
    d1 = np.array(cc.persistence_intervals_in_dimension(1), float).reshape(-1, 2)
    f0 = d0[np.isfinite(d0[:, 1])]
    f1 = d1[np.isfinite(d1[:, 1])]
    p0 = f0[:, 1] - f0[:, 0] if f0.size else np.zeros(0)
    p1 = f1[:, 1] - f1[:, 0] if f1.size else np.zeros(0)
    return {"n_h0_bars": float(d0.shape[0]), "n_h1_bars": float(d1.shape[0]),
            "total_persistence_h0": float(p0.sum()), "total_persistence_h1": float(p1.sum()),
            "max_persistence_h1": float(p1.max()) if p1.size else 0.0,
            "mean_persistence_h1": float(p1.mean()) if p1.size else 0.0}, (f0, f1)


def spearman_rho(x, y):
    from scipy.stats import rankdata
    rx, ry = rankdata(x), rankdata(y)
    rx = rx - rx.mean(); ry = ry - ry.mean()
    d = np.sqrt((rx ** 2).sum() * (ry ** 2).sum())
    return float((rx * ry).sum() / d) if d else 0.0


def perm_p(x, y, n_perm, rng):
    """Two-sided permutation p for Spearman rho, (r+1)/(n+1)."""
    obs = spearman_rho(x, y)
    y = np.asarray(y, float)
    cnt = sum(abs(spearman_rho(x, rng.permutation(y))) >= abs(obs) for _ in range(n_perm))
    return obs, float((cnt + 1) / (n_perm + 1))


def main():
    rng = np.random.default_rng(SEED)
    lab = np.load(LAB2)
    lab6 = np.load(LAB6)
    maps = np.load(MAPS, mmap_mode="r")
    assert np.allclose(lab6[:, :2], lab), "params_6 labels disagree with params_2 on (Om, s8)"

    # ---- effective-N diagnostic, BEFORE anything is correlated
    groups = collections.defaultdict(list)
    for i, row in enumerate(map(tuple, lab)):
        groups[row].append(i)
    gsizes = collections.Counter(len(v) for v in groups.values())
    diag = {"n_maps": int(lab.shape[0]), "n_distinct_simulations": len(groups),
            "group_size_histogram": {str(k): int(v) for k, v in sorted(gsizes.items())},
            "val_split_shared_label_rows": int(len(set(map(tuple, lab)) &
                                                set(map(tuple, np.load(B + "params_2/val_labels_LH_2.npy")))))}
    print("EFFECTIVE-N DIAGNOSTIC:", diag)

    # ---- stratified sample: ONE map per simulation, strata over Omega_m
    keys = sorted(groups, key=lambda r: r[0])                   # sorted by Omega_m
    edges = np.linspace(0, len(keys), N_SAMPLE + 1).astype(int)
    chosen = []
    for a, b in zip(edges[:-1], edges[1:]):
        if b > a:
            k = keys[a + int(rng.integers(b - a))]
            chosen.append(groups[k][0])                          # first map of that simulation
    chosen = np.array(sorted(set(chosen)))
    n = chosen.size
    print(f"sampled {n} maps, one per simulation, stratified into {N_SAMPLE} Omega_m strata")

    om, s8 = lab[chosen, 0], lab[chosen, 1]
    print(f"Omega_m range {om.min():.4f}-{om.max():.4f}; sigma_8 range {s8.min():.4f}-{s8.max():.4f}")

    # ---- per-map persistence
    rows, first_bars = [], None
    with Timer() as t:
        for j, i in enumerate(chosen):
            s, bars = map_stats(np.asarray(maps[i]))
            if first_bars is None:
                first_bars = (int(i), bars)
            s["map_index"] = int(i)
            s["Omega_m"] = float(lab[i, 0]); s["sigma_8"] = float(lab[i, 1])
            for c, nm in enumerate(["A_SN1", "A_AGN1", "A_SN2", "A_AGN2"]):
                s[nm] = float(lab6[i, 2 + c])
            rows.append(s)
            if (j + 1) % 50 == 0:
                print(f"  {j + 1}/{n}")
    wall, mem = t.sec, peak_mb()
    print(f"persistence on {n} maps: {wall:.1f}s, peak {mem:.0f} MB")

    X = {k: np.array([r[k] for r in rows]) for k in STATS}

    # ---- correlations, with the permutation null and the shuffled-label control
    res, shuffled = {}, {}
    for par, y in (("Omega_m", om), ("sigma_8", s8)):
        for st in STATS:
            rho, p = perm_p(X[st], y, N_PERM, np.random.default_rng(SEED + 7))
            res[f"{st}~{par}"] = {"rho": rho, "p_perm": p, "survives_bonferroni": bool(p < BONF)}
            ysh = np.random.default_rng(SEED + 99).permutation(y)
            rho_s, p_s = perm_p(X[st], ysh, 2000, np.random.default_rng(SEED + 7))
            shuffled[f"{st}~{par}"] = {"rho": rho_s, "p_perm": p_s,
                                       "survives_bonferroni": bool(p_s < BONF)}
            print(f"  {st:24s} ~ {par:8s} rho={rho:+.4f} p={p:.5f} "
                  f"{'PASS-BONF' if p < BONF else ''}   [shuffled rho={rho_s:+.4f} p={p_s:.4f}]")

    n_sig = sum(v["survives_bonferroni"] for v in res.values())
    n_sig_sh = sum(v["survives_bonferroni"] for v in shuffled.values())
    ctrl_ok = (n_sig_sh == 0)
    print(f"\n{n_sig}/{N_TESTS} survive Bonferroni {BONF:.5f}; shuffled-label control: "
          f"{n_sig_sh}/{N_TESTS} survive -> {'PASS' if ctrl_ok else 'FAIL'}")

    # ---------------------------------------------------------------- ingest
    notes = ("750 maps / 555 DISTINCT simulations (397 singletons, 126 pairs, 27 triples, 5 quads); "
             "293 simulations also appear in the val split. Stored values are "
             "renorm(log1p(minmax(HI))) with GLOBAL constants (process_hi_maps.py, sha256 "
             f"{sha256_file(B + 'process_hi_maps.py')[:16]}...), one monotone map applied identically "
             "to every map, so sublevel topology equals that of the raw HI field and persistence "
             "values are comparable across maps. Upstream raw source: CAMELS Multifield Dataset "
             "Maps_HI_IllustrisTNG_LH_z=0.00.npy. The HF repo card itself says 'Set the final license "
             "and citation before making this dataset public' -- the redistribution is a derived copy, "
             "so the CAMELS papers, not this repo, are the citable source.")
    with topodb() as db:
        db.add_dataset(id=DS, domain="astro",
                       title="CAMELS IllustrisTNG LH neutral-hydrogen maps (256x256, z=0)",
                       source=("https://huggingface.co/datasets/collins909/DDPM-HI-CAMELS-LH "
                               "(params_2/test_LH.npy); derived from the CAMELS Multifield Dataset"),
                       provenance="simulation", n_objects=750, ambient_dim=2,
                       units="dimensionless (globally log1p-renormalised HI column density)",
                       sha256=sha256_file(MAPS), local_path=MAPS, notes=notes)

        rid = add_run(db, dataset_id=DS, method="cubical", coeff_field=2, max_dim=1,
                      params={"n_maps_sampled": int(n), "n_strata": N_SAMPLE,
                              "sampling": ("one map per SIMULATION (first map of each distinct label "
                                           "row), simulations sorted by Omega_m and split into 200 "
                                           "equal strata, one simulation drawn per stratum, seed "
                                           f"{SEED}; N_eff = N_maps = {n} with no clustering"),
                              "grid": "256x256", "filtration": "sublevel cubical",
                              "statistics": STATS, "n_permutations": N_PERM,
                              "n_tests": N_TESTS, "bonferroni_threshold": BONF,
                              "effective_n_diagnostic": diag},
                      script=SCRIPT, command=COMMAND, tier="X",
                      preprocessing=("no transform applied here; the stored array is already "
                                     "renorm(log1p(minmax(HI))) with global constants"),
                      seed=str(SEED), wall_sec=wall, peak_mb=mem)
        mi, (f0, f1) = first_bars
        db.add_bars(rid, 0, [(a, b) for a, b in f0], top=25)
        db.add_bars(rid, 1, [(a, b) for a, b in f1], top=25)
        db.add_betti(rid, {0: int(rows[0]["n_h0_bars"]), 1: int(rows[0]["n_h1_bars"])})
        for key, v in res.items():
            db.add_statistic(rid, f"spearman_{key}", v["rho"],
                             null_model=(f"{N_PERM} random permutations of the {n} simulation labels "
                                         "(one map per simulation, so the permutation is at the "
                                         "simulation level and N_eff = 200)"),
                             n_null=N_PERM, p_value=v["p_perm"],
                             p_method=f"permutation, two-sided |rho|, (r+1)/(n+1); floor {1/(N_PERM+1):.2e}",
                             multiplicity=f"Bonferroni/{N_TESTS} -> threshold {BONF:.5f}")
        db.add_control(rid, "shuffle",
                       "Spearman on SHUFFLED labels must show no correlation surviving Bonferroni/12",
                       passed=ctrl_ok,
                       detail=(f"{n_sig_sh}/{N_TESTS} survive on shuffled labels "
                               f"(2000 permutations each); unshuffled {n_sig}/{N_TESTS}. "
                               f"max |rho| shuffled = "
                               f"{max(abs(v['rho']) for v in shuffled.values()):.4f}"))
        db.add_control(rid, "known_answer",
                       "Step 0 K4 recovered b0 = 7 for a grid with 7 placed wells on this same "
                       "GUDHI CubicalComplex pipeline (run 10)",
                       passed=True, detail="topodb_runs/astro/results/step0_known_answers.json")
        db.add_control(rid, "negative",
                       "effective-N check: 750 maps carry only 555 distinct simulations, so a naive "
                       "per-map permutation would be anti-conservative; avoided by sampling one map "
                       "per simulation",
                       passed=True, detail=str(diag))
        best = max(res.items(), key=lambda kv: abs(kv[1]["rho"]))
        db.add_finding(run_id=rid, dataset_id=DS,
                       claim=(f"Sublevel cubical persistence of CAMELS HI maps correlates with cosmology: "
                              f"{n_sig} of {N_TESTS} pre-declared (statistic, parameter) pairs survive "
                              f"Bonferroni/12 at p<{BONF:.5f} over {n} independent simulations. Strongest: "
                              f"{best[0]} Spearman rho = {best[1]['rho']:+.4f}, permutation p = "
                              f"{best[1]['p_perm']:.5f}."),
                       verdict="recovered" if n_sig else "null", tier="X",
                       caveat=("Tier X numerics, not a proof and not a measurement of Omega_m. It is a "
                               "correlation within ONE simulation suite (IllustrisTNG LH) at ONE "
                               "redshift, on maps that a third party had already min-max and log1p "
                               "normalised globally; the astrophysical feedback parameters "
                               "A_SN1/A_AGN1/A_SN2/A_AGN2 vary simultaneously and are NOT controlled "
                               "for, so part of any correlation may be feedback, not cosmology. "
                               "550 of the 750 maps were not used, to keep N_eff = N_maps."),
                       reference="topodb_runs/astro/results/camels.json")

    save_json("camels.json", {"seed": SEED, "run_id": rid, "n_maps": int(n), "diagnostic": diag,
                              "correlations": res, "shuffled_control": shuffled,
                              "n_significant": n_sig, "n_significant_shuffled": n_sig_sh,
                              "bonferroni": BONF, "per_map": rows})
    print(f"run_id {rid}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
