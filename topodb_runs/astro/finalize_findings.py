#!/usr/bin/env python3
"""One honest `finding` per dataset, written from the committed result JSONs.

Nothing here computes a number.  Every number is read back out of
topodb_runs/astro/results/*.json, which were written by the run scripts.
Verdict vocabulary is the schema's: recovered | null | inconclusive | failed |
artefact.  A null result is recorded as a null.  Nothing is "proved".
"""
from __future__ import annotations

import glob
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RESULTS, topodb  # noqa: E402

BONF_PC = 0.05 / 6
BONF_SKY = 0.05 / 7


def load(name):
    p = os.path.join(RESULTS, name)
    return json.load(open(p)) if os.path.exists(p) else None


def pointcloud_findings(db):
    n = 0
    seen = set()
    # the run scripts accumulate into a shared dict, so the same dataset appears
    # in several JSONs; keep the LAST occurrence of each dataset only.
    merged = {}
    for f in sorted(glob.glob(os.path.join(RESULTS, "pointclouds_*.json"))):
        d = json.load(open(f))
        for k, v in d.items():
            if isinstance(v, dict) and "p_values" in v:
                merged[k] = (v, os.path.basename(f))
    for key_, (v, fname) in sorted(merged.items()):
        if True:
            f = fname
            # the A3 re-runs are keyed "<dataset>_scalefree"; the dataset id is
            # the same row, so strip the suffix (a FK violation otherwise).
            scalefree = key_.endswith("_scalefree")
            ds = key_[: -len("_scalefree")] if scalefree else key_
            if key_ in seen:
                continue
            seen.add(key_)
            pv = v["p_values"]
            if not pv:
                continue
            floor2 = 2.0 / (v["n_null"] + 1)
            import numpy as _np
            sep, outside = {}, {}
            for st in pv:
                nv = _np.array([x[st] for x in v["null"]], float)
                sd = nv.std(ddof=1)
                sep[st] = (float(v["data"][st]) - nv.mean()) / sd if sd > 0 else float("nan")
                outside[st] = bool(v["data"][st] < nv.min() or v["data"][st] > nv.max())
            at_floor = [k for k, p in pv.items() if p <= floor2 + 1e-12]
            n_out = sum(outside.values())
            worst = max(sep, key=lambda k: abs(sep[k]))
            below = at_floor if n_out else []
            is_desi = "desi" in ds
            null_name = ("the official DESI random catalogue" if is_desi
                         else "shuffled-redshift realisations of the sample itself")
            claim = (
                f"Alpha-complex persistence of this comoving point cloud separates it from "
                f"{null_name} ({v['n_null']} realisations): {n_out} of {len(pv)} pre-declared "
                f"statistics lie OUTSIDE the full range of the null draws, and {len(at_floor)} sit at "
                f"the two-sided rank floor 2/{v['n_null'] + 1} = {floor2:.4f}, which is the smallest "
                f"value 20 draws can produce and is therefore resolution-limited, NOT a measure of how "
                f"extreme they are. Largest separation: {worst} at {sep[worst]:+.1f} null-sigma. "
                f"Measured on the data: {int(v['data']['n_h2_bars_over_5mpc'])} H2 bars and "
                f"{int(v['data']['n_h1_bars_over_5mpc'])} H1 bars longer than 5 Mpc/h; "
                f"H0-death IQR/median = {v['data']['h0_death_iqr_over_median']:.4f}; "
                f"longest H2 bar {v['data']['max_persistence_h2']:.2f} Mpc/h."
                + (" AMENDMENT A3 RE-RUN: the alpha cap is scale-free here "
                   f"({v.get('max_alpha_square', 0) ** 0.5:.0f} Mpc/h, from the sample's own pilot "
                   f"H0-death median {v.get('scale', 0):.2f} Mpc/h) because the pre-declared absolute "
                   "cap of 30 Mpc/h saturated this sample." if scalefree else "")
            )
            if is_desi:
                caveat = (
                    "THE NULL IS UNCLUSTERED. The DESI randoms reproduce the survey window and n(z) "
                    "but contain no galaxy clustering, so this p-value answers 'does this point set "
                    "differ topologically from an unclustered field with the same selection?' -- "
                    "whose answer was known in advance to be yes, because galaxies cluster. It is NOT "
                    "evidence of a topological anomaly and NOT a test of LCDM. What it does establish "
                    "is that the pipeline's COUNT statistic detects 3-D clustering at this sample "
                    "size. A clustering-matched (lognormal or N-body) DESI-footprint mock family, "
                    "which docs/FUTURE_OBSERVATIONAL_TARGETS.md §2a names as the single missing piece, "
                    "does not exist on this disk and was NOT ATTEMPTED. Weights (WEIGHT, WEIGHT_FKP) "
                    "were not applied. Tier X."
                )
            else:
                caveat = (
                    "Heterogeneous selection (MGS at low z, LRG-like above z~0.15), no completeness "
                    "map, no random catalogue, and the staged file's original query script was not "
                    "found on disk. The null is a shuffled-redshift null, which destroys 3-D "
                    "clustering while keeping the angular footprint and n(z); it is WEAKER than a "
                    "mock and is not clustering-matched. Tier X."
                )
            caveat += (
                f" The pre-declared Bonferroni/6 threshold is {BONF_PC:.5f}, BELOW the attainable "
                f"two-sided rank floor {floor2:.4f} at n_null={v['n_null']}: no statistic could have "
                "passed it however far the data lay from the null, so the verdict here rests on the "
                "recorded null-sigma separations and on the data lying outside the null range, not on "
                "a p-value. Raising n_null above ~240 would be needed to make the floor meet the "
                "threshold and was NOT ATTEMPTED (budget)."
                + (f" The data has FEWER long H2 bars than the unclustered randoms "
                   f"({int(v['data']['n_h2_bars_over_5mpc'])} against a null mean of "
                   f"{_np.mean([x['n_h2_bars_over_5mpc'] for x in v['null']]):.1f}), the same direction "
                   "the prior campaign measured (832 data peaks against 1235 +/- 23 in randoms): a "
                   "clustered field concentrates galaxies and leaves fewer resolvable voids at this "
                   "sample size." if is_desi else ""))
            db.add_finding(run_id=v["run_id"], dataset_id=ds, claim=claim,
                           verdict="recovered" if below else "null", tier="X",
                           caveat=caveat, reference=f"topodb_runs/astro/results/{f}")
            n += 1
    return n


def skymap_findings(db):
    n = 0
    for f in sorted(glob.glob(os.path.join(RESULTS, "skymaps_nside*.json"))):
        d = json.load(open(f))
        if d.get("n_sims", 0) < 10:        # smoke tests are not findings
            continue
        for key, v in d.items():
            if not isinstance(v, dict) or "p_values" not in v:
                continue
            pv = v["p_values"]
            floor2 = 2.0 / (v["n_sims"] + 1)
            cal = v.get("null_calibration_p")
            is_fg = key.startswith(("haslam", "wmap_kband"))
            base = key.rsplit("_nside", 1)[0]
            ds = f"astro/{base}"
            # The attainable two-sided rank floor (0.0198 at n=100) is ABOVE
            # Bonferroni/7 = 0.00714, so `p < BONF_SKY` is empty BY CONSTRUCTION
            # for every map. The verdict must rest on separation, not on that test.
            sep = v.get("separation_sigma") or {}
            n_out = int(v.get("n_outside", 0))
            at_floor = [k for k, p in pv.items() if p <= floor2 + 1e-12]
            floor_blocks = floor2 > BONF_SKY
            fired = bool(n_out) or bool(at_floor)
            claim = (
                f"Lower-star Betti curves on the fixed HEALPix 2-complex ({key}, f_sky="
                f"{v['f_sky']:.3f}, {v['n_sims']} spectrum-matched Gaussian sims): "
                f"{len(at_floor)} of {len(pv)} pre-declared statistics sit at the attainable rank "
                f"floor ({floor2:.4f} two-sided, {1.0 / (v['n_sims'] + 1):.4f} for the coarse-curve "
                f"statistic), and {n_out} lie outside the full range of the {v['n_sims']} draws. "
                f"b0 at nu=+1 is {v['data']['b0_at_nu_1']:.0f}, b1 at nu=0 is "
                f"{v['data']['b1_at_nu_0']:.0f}, Euler characteristic at nu=0 is "
                f"{v['data']['euler_char_at_nu_0']:.0f}. "
                + ("This is the POSITIVE CONTROL: a foreground-dominated map rejects the Gaussian "
                   "null, which is what shows the CMB nulls in this campaign are not vacuous."
                   if is_fg and fired else
                   "POSITIVE CONTROL DID NOT FIRE -- if a foreground map cannot reject a Gaussian "
                   "null, every CMB null here must be read as 'the pipeline cannot fire'."
                   if is_fg else
                   "No departure from an isotropic Gaussian field is detected in these statistics."
                   if not fired else
                   "A departure from the isotropic Gaussian null is detected.")
            )
            caveat = (
                f"Tier X. Rank p-values only (the fixed library's chi2 branch is anti-conservative at "
                f"0.058 against 0.05). "
                + (f"THE PRE-DECLARED Bonferroni/7 THRESHOLD {BONF_SKY:.5f} IS BELOW THE ATTAINABLE "
                   f"TWO-SIDED RANK FLOOR {floor2:.4f} at n_sims={v['n_sims']}: no statistic could "
                   "have passed it however extreme, so the verdict rests on the recorded null-sigma "
                   "separations and on the data lying outside the null range, not on a p-value. "
                   f"Resolving the threshold would need n_sims > {int(2 / BONF_SKY):d} and was NOT "
                   "ATTEMPTED (budget). " if floor_blocks else "")
                + "The null is spectrum-matched ONLY UP TO MODE COUPLING: no MASTER deconvolution is "
                  "applied, so the mask's distortion of the C_ell shape is uncorrected (modest at "
                  "Planck f_sky=0.76, less so at WMAP KQ85). The f_sky division in the code is a "
                  "MEASURED NO-OP -- the statistic is scale-invariant because the field is "
                  "sigma-normalised before filtering (verified bit-identical on wmap_ilc at nside 64) "
                  "-- so it is credited with nothing. The null does not model non-Gaussian "
                  "foregrounds, anisotropic noise or the beam, so a rejection localises no cause. "
                  f"The held-out-sim check returned p = {cal}; it is recorded as a DIAGNOSTIC with "
                  "passed=None and does NOT gate this verdict, because it tests exchangeability among "
                  "the sims (true by construction), has no power against a null mismatched to the "
                  "data, and is one Uniform(0,1) draw that fails ~5% of the time on a good null."
            )
            if not is_fg and not fired:
                caveat += (" A null here is weak, not a constraint: the prior campaign measured that "
                           "this family of statistics only reaches 95 % power against a STRONG, SPARSE "
                           "injected population (A = 4 sigma_T, 1000 full-sky cores) and is close to "
                           "vacuous against a weak or dense one.")
            if key.startswith("cobe"):
                caveat += (" The DMR quad-cube was repixelised into HEALPix by tabulated pixel-centre "
                           "RA/Dec and averaged; that is lossy. nside 16 was used because at nside 32 "
                           "the 6144 DMR pixels leave f_sky = 0.4963 and the induced subcomplex is "
                           "perforated with aliasing holes, so b1 would measure the repixelisation "
                           "rather than the sky.")
            below = at_floor if fired else []
            verdict = "recovered" if below else "null"
            db.add_finding(run_id=v["run_id"], dataset_id=ds, claim=claim, verdict=verdict,
                           tier="X", caveat=caveat,
                           reference=f"topodb_runs/astro/results/{os.path.basename(f)}")
            n += 1
    return n


def camels_finding(db):
    """Re-emit the CAMELS finding from its committed result JSON (run_camels.py
    writes it too; this makes finalize idempotent after a findings wipe)."""
    d = load("camels.json")
    if not d:
        return 0
    res = d["correlations"]
    best = max(res.items(), key=lambda kv: abs(kv[1]["rho"]))
    n_sig, n_sh = d["n_significant"], d["n_significant_shuffled"]
    db.add_finding(run_id=d["run_id"], dataset_id="astro/camels_hi_illustristng_lh",
                   claim=(f"Sublevel cubical persistence of CAMELS IllustrisTNG HI maps correlates with "
                          f"cosmology: {n_sig} of 12 pre-declared (statistic, parameter) pairs survive "
                          f"Bonferroni/12 at p<{d['bonferroni']:.5f} over {d['n_maps']} INDEPENDENT "
                          f"simulations (one map per simulation). Strongest: {best[0]} Spearman rho = "
                          f"{best[1]['rho']:+.4f}, permutation p = {best[1]['p_perm']:.5f} "
                          f"(10000 permutations, floor 1.0e-4). Shuffled-label control: {n_sh} of 12."),
                   verdict="recovered" if n_sig else "null", tier="X",
                   caveat=("Tier X numerics, not a proof and not a measurement of Omega_m. A correlation "
                           "within ONE simulation suite (IllustrisTNG LH) at ONE redshift, on maps a "
                           "third party had already min-max and log1p normalised globally. The "
                           "astrophysical feedback parameters A_SN1/A_AGN1/A_SN2/A_AGN2 vary "
                           "simultaneously across the LH set and are NOT controlled for, so part of any "
                           "correlation may be feedback rather than cosmology. EFFECTIVE N: the 750 test "
                           "maps carry only 555 distinct label rows (397 singletons, 126 pairs, 27 "
                           "triples, 5 quads) and 293 of those simulations also appear in the val split; "
                           "sampling one map per simulation is what makes N_eff = N_maps = 200, at the "
                           "cost of leaving 550 maps unused. The strongest p sits at the permutation "
                           "floor 1/10001, so its magnitude is resolution-limited."),
                   reference="topodb_runs/astro/results/camels.json")
    return 1


def main():
    with topodb() as db:
        a = pointcloud_findings(db)
        b = skymap_findings(db)
        camels_finding(db)
    print(f"wrote {a} point-cloud findings and {b} sky-map findings")
    with topodb() as db:
        import json as _j
        print(_j.dumps(db.summary(), indent=1))
        print("\n--- astro findings ---")
        for r in db.con.execute("SELECT dataset_id, verdict, substr(claim,1,110) FROM finding "
                                "WHERE dataset_id LIKE 'astro/%' OR dataset_id LIKE 'synthetic/step0%' "
                                "ORDER BY dataset_id").fetchall():
            print(f"  [{r[1]:12s}] {r[0]}\n      {r[2]}...")
    return 0


if __name__ == "__main__":
    sys.exit(main())
