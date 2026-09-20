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
    for f in sorted(glob.glob(os.path.join(RESULTS, "pointclouds_*.json"))):
        d = json.load(open(f))
        for ds, v in d.items():
            if not isinstance(v, dict) or "p_values" not in v:
                continue
            pv = v["p_values"]
            if not pv:
                continue
            floor = 1.0 / (v["n_null"] + 1)
            at_floor = [k for k, p in pv.items() if p <= floor + 1e-12]
            below = [k for k, p in pv.items() if p < BONF_PC]
            is_desi = "desi" in ds
            null_name = ("the official DESI random catalogue" if is_desi
                         else "shuffled-redshift realisations of the sample itself")
            claim = (
                f"On {v['n_null']}-realisation comparison, {len(below)} of {len(pv)} pre-declared "
                f"alpha-complex statistics of this comoving point cloud differ from "
                f"{null_name} below Bonferroni/6 = {BONF_PC:.4f}. "
                f"Measured on the data: {int(v['data']['n_h2_bars_over_5mpc'])} H2 bars and "
                f"{int(v['data']['n_h1_bars_over_5mpc'])} H1 bars longer than 5 Mpc/h; "
                f"H0-death IQR/median = {v['data']['h0_death_iqr_over_median']:.4f}; "
                f"longest H2 bar {v['data']['max_persistence_h2']:.2f} Mpc/h."
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
            if at_floor:
                caveat += (f" {len(at_floor)} statistic(s) sit AT the rank floor {floor:.4f}, which is "
                           "resolution-limited, not a measure of how extreme they are.")
            db.add_finding(run_id=v["run_id"], dataset_id=ds, claim=claim,
                           verdict="recovered" if below else "null", tier="X",
                           caveat=caveat, reference=f"topodb_runs/astro/results/{os.path.basename(f)}")
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
            floor = 1.0 / (v["n_sims"] + 1)
            below = {k: p for k, p in pv.items() if p < BONF_SKY}
            cal = v.get("null_calibration_p")
            cal_ok = cal is not None and cal > 0.05
            is_fg = key.startswith(("haslam", "wmap_kband"))
            base = key.rsplit("_nside", 1)[0]
            ds = f"astro/{base}"
            claim = (
                f"Lower-star Betti curves on the fixed HEALPix 2-complex ({key}, f_sky="
                f"{v['f_sky']:.3f}, {v['n_sims']} spectrum-matched Gaussian sims): "
                f"{len(below)} of {len(pv)} pre-declared statistics reject the Gaussian null below "
                f"Bonferroni/7 = {BONF_SKY:.4f}. b0 at nu=+1 is {v['data']['b0_at_nu_1']:.0f}, "
                f"b1 at nu=0 is {v['data']['b1_at_nu_0']:.0f}, Euler characteristic at nu=0 is "
                f"{v['data']['euler_char_at_nu_0']:.0f}. "
                + ("As expected for a foreground-dominated map, and this is what shows the CMB nulls "
                   "are not vacuous." if is_fg else
                   "No departure from an isotropic Gaussian field is detected in these statistics."
                   if not below else "")
            )
            caveat = (
                f"Tier X. Rank p-values only (the fixed library's chi2 branch is anti-conservative at "
                f"0.058 against 0.05); the floor is 1/{v['n_sims'] + 1} = {floor:.4f} and a p there is "
                "resolution-limited. The null is spectrum-matched and isotropic: it does not model "
                "non-Gaussian foregrounds, anisotropic noise or the beam, so a rejection localises no "
                "cause. The pseudo-C_ell was divided by f_sky before synthesis; the null-calibration "
                f"control (held-out sim vs the other n-1) returned p = {cal}, "
                + ("which is not extreme, so the normalisation is not manufacturing the result."
                   if cal_ok else "WHICH IS EXTREME -- this run's p-values should not be trusted.")
            )
            if not is_fg and not below:
                caveat += (" A null here is weak, not a constraint: the prior campaign measured that "
                           "this family of statistics only reaches 95 % power against a STRONG, SPARSE "
                           "injected population (A = 4 sigma_T, 1000 cores) and is close to vacuous "
                           "against a weak or dense one.")
            if key.startswith("cobe"):
                caveat += (" The DMR quad-cube was repixelised into HEALPix by tabulated pixel-centre "
                           "RA/Dec and averaged; that is lossy and the resulting nside 32 complex is "
                           "coarse.")
            verdict = "recovered" if below else "null"
            if not cal_ok:
                verdict = "inconclusive"
            db.add_finding(run_id=v["run_id"], dataset_id=ds, claim=claim, verdict=verdict,
                           tier="X", caveat=caveat,
                           reference=f"topodb_runs/astro/results/{os.path.basename(f)}")
            n += 1
    return n


def main():
    with topodb() as db:
        a = pointcloud_findings(db)
        b = skymap_findings(db)
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
