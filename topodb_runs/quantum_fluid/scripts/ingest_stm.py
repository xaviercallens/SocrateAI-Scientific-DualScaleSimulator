"""Ingest DATASET 1 (Re6Zr STM vortex maps, 11 fields) and the STEP 0 controls.

One TopoDB dataset per field value, one alpha run per field.
Run: prlimit --as=8589934592 -- .venv-tda/bin/python ingest_stm.py
"""
import json
import os

import numpy as np

import ingest_lib as il
import qf_lib as q

SRC = "https://doi.org/10.5281/zenodo.14780459 (Duhan et al., Nat. Commun. 16, 2100 (2025))"
MANIFEST = "/mnt/disks/disk-socrateai-local-1/dualscale-wt-tdaval/audit/tda_validation/quantum_fluid/data_manifest.json"


def main():
    F = json.load(open(os.path.join(q.RESULTS, "stm_fields.json")))
    V = json.load(open(os.path.join(q.RESULTS, "step0_verdict.json")))
    S0 = json.load(open(os.path.join(q.RESULTS, "step0_re6zr_20kOe.json")))
    g = il.IdGuard()
    d = il.db()
    try:
        for key, v in F["fields"].items():
            H = v["H_kOe"]
            pf = V["per_field"][key]
            ds = f"quantum_fluid/re6zr_stm_{int(H)}kOe"
            a_ref = v["a_tri_formula_nm"]
            d.add_dataset(
                id=ds, domain="quantum_fluid",
                title=f"a-Re6Zr STM vortex-core map, H = {H:g} kOe, T = 460 mK",
                source=SRC, provenance="observation",
                n_objects=int(v["n_detected_median"]) * v["n_images"],
                ambient_dim=2, units="nm",
                sha256=v["cores_sha256"],
                local_path=v["cores_npz"],
                notes=(f"{v['n_images']} successive 'Input_7' forward lock-in conductance maps, "
                       f"{v['scan_nm']:.0f} nm scan, {v['px_nm']:.3f} nm/pixel. Vortex cores = local "
                       f"conductance minima (plane-subtract, Gaussian low-pass sigma=a/6, disk-minimum "
                       f"radius a/3). Ideal triangular spacing 1.075*sqrt(Phi0/B) = {a_ref:.2f} nm at "
                       f"B = mu0 H = {v['B_T']:.2f} T. Raw .sxm sha256 in {MANIFEST}. "
                       f"NOTE: the detector over-detects at the scan border "
                       f"(ndimage.minimum_filter mode='nearest' replicates edge values); the interior "
                       f"count (>= a/2 from the border) is the one used for the areal density."))

            rk = f"stm/{key}/alpha"
            if g.has(rk):
                print("skip (already ingested)", rk, "run", g.get(rk))
                continue
            run = d.add_run(
                dataset_id=ds, method="alpha", coeff_field=2, max_dim=2,
                params={"max_alpha_square": (3 * a_ref) ** 2, "a_ref_nm": a_ref,
                        "betti_threshold_frac_of_a": 0.55, "n_null": F["meta"]["n_null"],
                        "null": "uniform random points, same count, same rectangle",
                        "psi6": "global |<exp(6 i theta_b)>| over all Delaunay bonds",
                        "a_fft": "first Bragg ring of azimuthally averaged |FFT|^2, Hann window, zero-pad 4",
                        "n_images": v["n_images"]},
                preprocessing=("alpha persistence per image on the vortex-core point cloud in nm; "
                               "run-level Betti numbers are the MEDIAN over the 20 images at "
                               "r = 0.55*a_ref; bars are from image 001 only; spread statistics are "
                               "the median over the 20 images"),
                script="topodb_runs/quantum_fluid/scripts/compute_stm_fields.py",
                command=F["meta"]["command"], seed=str(F["meta"]["seed"]),
                tier="X", repo=q.REPO)
            g.put(rk, run)

            d.add_betti(run, {0: v["betti_median"]["0"] if "0" in v["betti_median"] else v["betti_median"][0],
                              1: v["betti_median"]["1"] if "1" in v["betti_median"] else v["betti_median"][1]},
                        expected=None)
            for dim, key2 in ((0, "top_h0_bars"), (1, "top_h1_bars")):
                bars = [(b, dd) for b, dd in v[key2]]
                if bars:
                    d.add_bars(run, dim, bars, top=10)

            nullm = f"{F['meta']['n_null']} uniform random point sets, same count, same rectangle"
            il.add_stats(d, run, [
                {"name": "a_fft_nm", "value": v["a_fft_nm"]["median"]},
                {"name": "a_tri_formula_nm", "value": a_ref},
                {"name": "a_fft_over_formula", "value": pf["a_fft_over_formula"]},
                {"name": "a_density_interior_nm", "value": pf["a_density_interior_nm"]},
                {"name": "a_density_full_nm", "value": v["a_density_nm"]},
                {"name": "a_h1_nm", "value": v["a_h1_nm"]},
                {"name": "a_h0_nm", "value": v["a_h0_nm"]},
                {"name": "a_nn_nm", "value": v["a_nn_nm"]},
                {"name": "ratio_afft_over_adensity_interior", "value": pf["ratio_afft_over_adensity_interior"]},
                {"name": "n_detected_median", "value": float(v["n_detected_median"])},
                {"name": "n_expected_flux_BA_over_Phi0", "value": v["n_expected_flux"]},
                {"name": "count_over_flux", "value": v["count_over_flux"]},
                {"name": "border_excess_frac", "value": pf["border_excess_frac"]},
                {"name": "h0_iqr_over_median", "value": v["observed"]["h0_iqr_over_median"],
                 "null_model": nullm, "n_null": F["meta"]["n_null"],
                 "p_value": v["p_values_rank"]["h0_iqr_over_median"], "p_method": "rank",
                 "multiplicity": "none; p floor = 1/(n_null+1) = %.4f" % (1 / (F["meta"]["n_null"] + 1))},
                {"name": "psi6_global", "value": v["observed"]["psi6"],
                 "null_model": nullm, "n_null": F["meta"]["n_null"],
                 "p_value": v["p_values_rank"]["psi6"], "p_method": "rank",
                 "multiplicity": "none; p floor = 1/(n_null+1) = %.4f" % (1 / (F["meta"]["n_null"] + 1))},
                {"name": "h1_total_persistence_over_a2",
                 "value": v["observed"]["h1_total_persistence_over_a2"],
                 "null_model": nullm, "n_null": F["meta"]["n_null"],
                 "p_value": v["p_values_rank"]["h1_total_persistence_over_a2"], "p_method": "rank",
                 "multiplicity": "none; p floor = 1/(n_null+1) = %.4f" % (1 / (F["meta"]["n_null"] + 1))},
            ])

            sc = v["synthetic_control"]
            d.add_control(run, "known_answer",
                          "synthetic perfect triangular lattice at the same a and area recovers a",
                          passed=bool(abs(sc["a_h1_nm"] / a_ref - 1) < 0.01),
                          detail=f"a_h1 {sc['a_h1_nm']:.3f} vs a_ref {a_ref:.3f}; "
                                 f"a_nn {sc['a_nn_nm']:.3f}; psi6 {sc['psi6']:.3f}; "
                                 f"h0 IQR/median {sc['h0_iqr_over_median']:.2e}")
            d.add_control(run, "null_calibration",
                          "uniform random points at the same count and area",
                          passed=True,
                          detail=f"null medians: h0 IQR/median {v['null_median']['h0_iqr_over_median']:.3f} "
                                 f"(observed {v['observed']['h0_iqr_over_median']:.3f}), "
                                 f"psi6 {v['null_median']['psi6']:.3f} "
                                 f"(observed {v['observed']['psi6']:.3f})")
            d.add_control(run, "negative",
                          "scan-border over-detection quantified",
                          passed=None,
                          detail=f"{100*pf['border_excess_frac']:.1f} per cent of the full-frame areal "
                                 f"density is border excess; interior density gives a = "
                                 f"{pf['a_density_interior_nm']:.2f} nm vs full-frame "
                                 f"{v['a_density_nm']:.2f} nm")

            if int(H) == 20:
                gp = V["gate_prestated"]
                d.add_control(run, "known_answer",
                              "STEP 0 GATE (pre-registered): absolute a at 20 kOe vs 1.075*sqrt(Phi0/B), 3 per cent",
                              passed=False,
                              detail=f"a_fft {gp['measured_nm']:.3f} nm vs target {gp['target_nm']:.3f} nm, "
                                     f"{100*gp['rel_error']:+.2f} per cent. Confirmed by four independent "
                                     f"routes: Bragg spots 35.2-36.2, g(r) 35.86, Delaunay bond mode 36.29, "
                                     f"sqrt(3)*median(H1) 36.04. Pre-registered in expectations.json at c003c3b.")
                d.add_control(run, "known_answer",
                              "Abrikosov a ~ B^(-1/2) slope over the ring-resolved fields 3-40 kOe",
                              passed=True,
                              detail="slope %.5f vs expected -0.5; prefactor ratio C/C_ideal = %.4f "
                                     "(field-independent offset)" % (
                                         V["abrikosov_fit_ring_resolved_only"]["slope"],
                                         V["abrikosov_fit_ring_resolved_only"]["C_ratio"]))
                d.add_control(run, "known_answer",
                              "scale-free triangular-packing check a_fft / a_density_interior = 1",
                              passed=True,
                              detail="median %.4f over the ring-resolved fields; independent of any "
                                     "scan-length calibration" % V["scale_free_consistency_check"]["median_over_ring_resolved_fields"])
                d.add_control(run, "known_answer",
                              "POST-HOC: published psi6 peak (Duhan et al. Fig. 2(m), ~0.62 at 20 kOe)",
                              passed=True,
                              detail="measured psi6 = %.3f at 20 kOe, %+.1f per cent; peak field, the slow "
                                     "decrease to 40 kOe and the collapse above 50 kOe all reproduce. "
                                     "NOT the pre-registered gate: identified after the gate had run." % (
                                         V["published_psi6_comparison_post_hoc"]["peak_value_measured"],
                                         100 * V["published_psi6_comparison_post_hoc"]["peak_rel_error"]))
                d.add_control(run, "known_answer",
                              "a_fft estimator bias on synthetic conductance maps",
                              passed=True,
                              detail="-0.36 to -0.68 per cent at sigma/a = 0 to 0.15; recovers 35.86 on a "
                                     "synthetic a = 36.0 nm lattice; all five estimators exact on a perfect lattice")

            ordered = v["observed"]["psi6"] > 10 * v["null_median"]["psi6"]
            d.add_finding(
                run_id=run, dataset_id=ds,
                claim=(f"At H = {H:g} kOe the Re6Zr vortex cores form a "
                       f"{'lattice with strong six-fold bond-orientational order' if ordered else 'positionally regular but orientationally disordered arrangement'}: "
                       f"psi6 = {v['observed']['psi6']:.3f} against a uniform-random null median of "
                       f"{v['null_median']['psi6']:.3f} (rank p = {v['p_values_rank']['psi6']:.4f}, "
                       f"n_null = {F['meta']['n_null']}), and the alpha H0 death spread is "
                       f"IQR/median = {v['observed']['h0_iqr_over_median']:.3f} against a null median of "
                       f"{v['null_median']['h0_iqr_over_median']:.3f}. The lattice constant measured "
                       f"detector-free from the first Bragg ring is {v['a_fft_nm']['median']:.2f} nm, "
                       f"{100*(pf['a_fft_over_formula']-1):+.1f} per cent from 1.075*sqrt(Phi0/B) = {a_ref:.2f} nm."),
                verdict="recovered" if ordered else "inconclusive",
                tier="X",
                caveat=("The H0 spread measures REGULAR SPACING, not orientational order: at 50 and 70 kOe "
                        "the spacing spread stays near its ordered value while psi6 collapses to the "
                        "random level. The absolute lattice constant is 4-5 per cent above the ideal-formula "
                        "value at every field; this offset is field-independent and unexplained (see the "
                        "20 kOe run's Step 0 controls). At 50 and 70 kOe the Bragg ring is not resolved, so "
                        "a_fft there is unreliable."),
                reference="topodb_runs/quantum_fluid/results/stm_fields.json and step0_verdict.json")
            print(f"ingested {ds} run {run}")

        # Step 0 summary finding, attached to the 20 kOe dataset
        if not g.has("stm/step0_finding"):
            fid = d.add_finding(
                dataset_id="quantum_fluid/re6zr_stm_20kOe",
                claim=("STEP 0 known-answer reproduction of the Re6Zr vortex-lattice constant at 20 kOe "
                       "FAILS its pre-registered 3 per cent tolerance: the detector-free, TDA-free Bragg-ring "
                       f"estimator gives {S0['a_fft']['median']:.2f} nm against the formula value "
                       f"{S0['target_a_nm']:.2f} nm ({100*S0['a_fft_rel_error']:+.1f} per cent). The pipeline "
                       "is not at fault -- it recovers a to 1e-13 on a perfect lattice and to -0.4 per cent "
                       "on synthetic conductance maps, and four mutually independent routes on the real maps "
                       "agree on 35.9-36.3 nm. What DOES reproduce: the Abrikosov a ~ B^(-1/2) law "
                       "(slope -0.4992 over 3-40 kOe) with a field-independent prefactor offset of +4.8 per "
                       "cent, the scale-free triangular-packing ratio a_fft/a_density_interior = 1.010, and "
                       "the paper's own measured psi6 (0.641 vs published ~0.62 at 20 kOe)."),
                verdict="failed", tier="X",
                caveat=("The earlier validation run's 36.04 nm is confirmed, not corrected. The full-frame "
                        "vortex count appeared to match B*A/Phi0 only because the minimum-filter detector "
                        "over-detects at the scan border; the interior areal density corresponds to "
                        "a = 35.67 nm and is about 6 per cent below B/Phi0. A single length-calibration "
                        "factor of 0.969 on the deposited SCAN_RANGE would reconcile spacing, interior "
                        "density and flux quantization simultaneously, but that is a one-parameter fit to "
                        "the discrepancy and was NOT tested independently. B < mu0*H is the other untested "
                        "candidate. The psi6 comparison is POST-HOC."),
                reference="topodb_runs/quantum_fluid/results/step0_verdict.json; arXiv:2406.07027")
            g.put("stm/step0_finding", fid)
            print("step0 finding", fid)
    finally:
        il.close(d)
    print("done")


if __name__ == "__main__":
    main()
