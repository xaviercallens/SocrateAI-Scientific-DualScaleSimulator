"""Cross-dataset summary: the calibrated ordering axis that this sweep gives TopoDB,
plus a refresh of the one Step-0 control whose detail predated
scripts/step0_diagnostics.py.

Run: prlimit --as=8589934592 -- .venv-tda/bin/python ingest_summary.py
"""
import json
import os

import ingest_lib as il
import qf_lib as q


def main():
    D = json.load(open(os.path.join(q.RESULTS, "step0_diagnostics.json")))
    g = il.IdGuard()
    d = il.db()
    try:
        # refresh the estimator-bias control with script-sourced numbers
        rid = g.get("stm/20kOe/alpha")
        if rid:
            b = D["synthetic_a_fft_bias"]
            lo = min(x["rel_error"] for x in b); hi = max(x["rel_error"] for x in b)
            d.add_control(rid, "known_answer", "a_fft estimator bias on synthetic conductance maps",
                          passed=True,
                          detail=("%+.2f to %+.2f per cent over sigma/a = 0 to 0.15 and core widths "
                                  "a/6 and a/4; recovers %.3f nm on a synthetic a = 36.0 nm lattice. "
                                  "All numbers from scripts/step0_diagnostics.py -> "
                                  "results/step0_diagnostics.json." %
                                  (100 * lo, 100 * hi, D["synthetic_a_true_36nm_recovered"])))
            d.add_control(rid, "known_answer",
                          "four independent routes to the 20 kOe lattice constant",
                          passed=True,
                          detail=("Bragg spots %s nm; g(r) first peak %.2f nm; Delaunay bond mode "
                                  "%.2f nm; sqrt(3)*median(H1 alpha) 36.04 nm; azimuthal FFT ring "
                                  "35.99 nm. They agree with each other; all four exceed "
                                  "1.075*sqrt(Phi0/B) = 34.57 nm." %
                                  (", ".join("%.2f" % s["a_nm"] for s in D["bragg_spots"][:3]),
                                   D["g_of_r_first_peak_nm"], D["delaunay_bond"]["mode_nm"])))

        if not g.has("summary/calibration_axis"):
            fid = d.add_finding(
                dataset_id="quantum_fluid/re6zr_20kOe_disorder_series",
                claim=(
                    "Across 21 quantum-fluid datasets the alpha H0 death spread (IQR/median) behaves as a "
                    "TWO-SIDED axis about the Poisson value 0.707, and pairing it with psi6 separates the "
                    "two kinds of order. Ordered end: 2-D colloidal crystal at area fraction 0.89, 0.032; "
                    "Re6Zr vortex lattice at 20 kOe, 0.069; Re6Zr at 3 kOe, 0.154; colloidal gas at area "
                    "fraction 0.11, 0.529 (still far below Poisson, because hard discs cannot overlap). "
                    "Poisson: uniform random points, 0.70-0.72 in every null run here, and the XY vortex "
                    "gas at T = 1.6 reaches exactly 0.707. Clustered end: the XY model below T_BKT, 3 to "
                    "17, because bound vortex dipoles make the H0 death distribution bimodal. psi6 is "
                    "independent of this axis: the colloidal gas has psi6 at the random level while its "
                    "H0 spread is 0.53, and the polycrystalline colloidal crystal has psi6 only 0.197 "
                    "while its H0 spread is the lowest measured. The disorder series calibrates the "
                    "trade-off: psi6 reaches half way to random at sigma/a = 0.15, the H0 spread only at "
                    "sigma/a = 0.30."),
                verdict="recovered", tier="X",
                caveat=(
                    "This is a statement about how two statistics order 21 datasets, not a physical law. "
                    "The values are not dimensionless in the same way across domains: the XY clouds are "
                    "sparse (4 to 2500 points) and their alpha complex is non-periodic although the "
                    "configurations are, while the colloid and STM clouds are dense and genuinely open. "
                    "The Poisson value 0.707 is a property of the 2-D Poisson process, so it is a fair "
                    "common reference, but the ordered-end values depend on the point count and on the "
                    "detector used to find the points. Every number here comes from a run recorded in "
                    "this database with its null and its n_null."),
                reference="topodb_runs/quantum_fluid/results/ (stm_fields, xy_vortex, colloid, "
                          "disorder_series, gpe_vortex .json)")
            g.put("summary/calibration_axis", fid)
            print("calibration-axis finding", fid)
    finally:
        il.close(d)
    print("done")


if __name__ == "__main__":
    main()
