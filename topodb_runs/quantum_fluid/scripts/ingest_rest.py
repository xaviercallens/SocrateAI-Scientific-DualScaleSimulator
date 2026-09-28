"""Ingest DATASETS 2-5: XY model, rotating GPE, 2-D colloidal dispersions, disorder series.

Run: prlimit --as=8589934592 -- .venv-tda/bin/python ingest_rest.py
"""
import json
import os

import numpy as np

import ingest_lib as il
import qf_lib as q

R = q.RESULTS
NULL_UNI = "uniform random points, same count, same box"


def J(name):
    return json.load(open(os.path.join(R, name)))


# ---------------------------------------------------------------- XY model
def do_xy(d, g):
    X = J("xy_vortex.json")
    M = X["meta"]
    hc = M["helicity_crosscheck_not_recomputed"]
    for Ls, blk in X["by_L"].items():
        L = int(Ls)
        ds = f"quantum_fluid/xy_model_L{L}"
        temps = blk["temperatures"]
        d.add_dataset(
            id=ds, domain="quantum_fluid",
            title=f"2-D XY model vortex point clouds, L = {L}, {len(temps)} temperatures",
            source=M["source"], provenance="simulation",
            n_objects=int(sum(t.get("n_vortices_median") or 0 for t in temps.values())),
            ambient_dim=2, units="lattice spacings",
            local_path=os.path.join(q.DATA_ROOT, f"xy/L{L}"),
            notes=(f"Vortex cores recomputed here from the stored phase configurations by plaquette "
                   f"phase winding (periodic), {M['n_config']} configurations per temperature; the "
                   f"earlier run's cached TDA files were not used. NON-TDA CROSS-CHECK (measured by the "
                   f"earlier run, not re-derived here): helicity-modulus crossing T_Ups = "
                   f"{hc['T_Ups'][Ls]:.4f} at this L; the three sizes extrapolate as T_inf + c/(ln L)^2 "
                   f"to {hc['T_inf_extrapolated']:.4f} against the published T_BKT = "
                   f"{hc['T_BKT_published']}."))
        for tk, t in sorted(temps.items()):
            if "observed" not in t:
                continue
            rk = f"xy/L{L}/T{tk}/alpha"
            if g.has(rk):
                continue
            run = d.add_run(
                dataset_id=ds, method="alpha", coeff_field=2, max_dim=2,
                params={"T": t["T"], "L": L, "n_config": M["n_config"],
                        "max_alpha_square": (3 * t["a_ref_from_density"]) ** 2,
                        "a_ref_from_density": t["a_ref_from_density"],
                        "betti_threshold_frac_of_a": 0.55, "n_null": t["n_null"],
                        "null": NULL_UNI, "vortex_detector": M["vortex_detector"]},
                preprocessing=("alpha persistence per configuration on the vortex point cloud in lattice "
                               "units; the alpha complex is NOT periodic although the configuration is, "
                               "so points near the box edge see a boundary. Run-level numbers are medians "
                               "over the configurations; bars are from configuration 0."),
                script="topodb_runs/quantum_fluid/scripts/compute_xy.py",
                command=M["command"], seed=str(M["seed"]), tier="X", repo=q.REPO)
            g.put(rk, run)
            d.add_betti(run, {0: t["betti_median"]["0"], 1: t["betti_median"]["1"]})
            for dim, kk in ((0, "top_h0_bars"), (1, "top_h1_bars")):
                if t[kk]:
                    d.add_bars(run, dim, [(b, dd) for b, dd in t[kk]], top=10)
            nm = f"{t['n_null']} {NULL_UNI}"
            il.add_stats(d, run, [
                {"name": "n_vortices_median", "value": t["n_vortices_median"]},
                {"name": "vortex_density_per_site", "value": t["vortex_density_per_site"]},
                {"name": "h0_median", "value": t["h0_median"]},
                {"name": "psi6_global", "value": t["observed"]["psi6"], "null_model": nm,
                 "n_null": t["n_null"], "p_value": t["p_values_rank"]["psi6"], "p_method": "rank",
                 "multiplicity": "none; p floor = %.4f" % (1 / (t["n_null"] + 1))},
                {"name": "h0_iqr_over_median", "value": t["observed"]["h0_iqr_over_median"],
                 "null_model": nm, "n_null": t["n_null"],
                 "p_value": t["p_values_rank"]["h0_iqr_over_median"], "p_method": "rank",
                 "multiplicity": "none; p floor = %.4f" % (1 / (t["n_null"] + 1))},
            ])
            d.add_control(run, "null_calibration", NULL_UNI, passed=True,
                          detail=f"null medians psi6 {t['null_median']['psi6']:.3f}, "
                                 f"h0 IQR/median {t['null_median']['h0_iqr_over_median']:.3f}; "
                                 f"observed {t['observed']['psi6']:.3f} and "
                                 f"{t['observed']['h0_iqr_over_median']:.3f}")

        # lower-star / cubical run carrying the MANDATORY site-shuffle control
        sh = blk["site_shuffle_control"]
        rk = f"xy/L{L}/cubical_shuffle"
        if not g.has(rk) and sh:
            run = d.add_run(
                dataset_id=ds, method="cubical", coeff_field=2, max_dim=1,
                params={"field": "cos(theta) on the L x L site grid", "threshold": 0.0,
                        "statistic": "b0 of the sublevel set at threshold 0, per site",
                        "temperatures": sorted(sh.keys()), "n_config": M["n_config"]},
                preprocessing="sublevel-set cubical persistence on the raw phase field; "
                              "site-shuffle = random permutation of the field values across sites",
                script="topodb_runs/quantum_fluid/scripts/compute_xy.py",
                command=M["command"], seed=str(M["seed"]), tier="X", repo=q.REPO)
            g.put(rk, run)
            for tk, s in sorted(sh.items()):
                il.add_stats(d, run, [
                    {"name": f"lower_star_b0_per_site_T{tk}", "value": s["real_b0_per_site"]},
                    {"name": f"lower_star_b0_per_site_shuffled_T{tk}", "value": s["shuffled_b0_per_site"]},
                    {"name": f"shuffle_diff_over_se_T{tk}", "value": s["diff_over_se"]},
                ])
                d.add_control(run, "shuffle", f"site-shuffle control at T = {tk}",
                              passed=bool(s["separates"]),
                              detail=f"real {s['real_b0_per_site']:.5f} vs shuffled "
                                     f"{s['shuffled_b0_per_site']:.5f} per site, diff/se = "
                                     f"{s['diff_over_se']:+.2f} over {s['n_config']} configurations")
            signs = [np.sign(s["diff_over_se"]) for s in sh.values()]
            d.add_finding(
                run_id=run, dataset_id=ds,
                claim=(f"For the XY model at L = {L}, the lower-star statistic b0 of the sublevel set of "
                       f"cos(theta) at threshold 0 does differ from a site-shuffled field, but the sign of "
                       f"the difference REVERSES between low and high temperature: "
                       + ", ".join(f"T={tk} diff/se {s['diff_over_se']:+.2f}" for tk, s in sorted(sh.items()))
                       + ". It is therefore not a monotone order parameter and no topological "
                         "interpretation of its temperature dependence is attributable."),
                verdict="artefact" if len(set(signs)) > 1 else "inconclusive", tier="X",
                caveat="The binding limit inherited from the earlier validation run is confirmed and "
                       "sharpened: the lower-star path on rough XY fields is sensitive to the field's "
                       "value histogram, not only to its spatial arrangement. The shuffle control "
                       "separates, so the statistic is not pure noise, but the sign reversal near "
                       "T ~ 1 means a crossing in this statistic cannot be attributed to the BKT "
                       "transition.",
                reference="topodb_runs/quantum_fluid/results/xy_vortex.json")

        # one finding per XY dataset
        usable = {k: v for k, v in temps.items() if "observed" in v}
        if usable and not g.has(f"xy/L{L}/finding"):
            lo = min(usable, key=lambda k: usable[k]["T"]); hi = max(usable, key=lambda k: usable[k]["T"])
            fid = d.add_finding(
                dataset_id=ds,
                claim=(f"XY model at L = {L}: the vortex point cloud shows NO six-fold orientational order "
                       f"at any temperature (psi6 is statistically indistinguishable from a uniform-random "
                       f"null at every T, rank p between "
                       f"{min(v['p_values_rank']['psi6'] for v in usable.values()):.3f} and "
                       f"{max(v['p_values_rank']['psi6'] for v in usable.values()):.3f}) -- it is a vortex "
                       f"gas, not a lattice. What the alpha H0 spread DOES track is vortex pairing: "
                       f"IQR/median falls monotonically from {usable[lo]['observed']['h0_iqr_over_median']:.2f} "
                       f"at T = {usable[lo]['T']:.2f} to "
                       f"{usable[hi]['observed']['h0_iqr_over_median']:.3f} at T = {usable[hi]['T']:.2f}, "
                       f"reaching the Poisson value ~0.707 at the highest temperature. Below T_BKT the "
                       f"few vortices sit in tightly bound dipoles, which makes the H0 death distribution "
                       f"strongly bimodal and the spread LARGE; above it the free vortex gas is Poisson."),
                verdict="recovered", tier="X",
                caveat=("This is a one-sided statement about a statistic, not a measurement of T_BKT: no "
                        "crossing temperature was fitted from TDA here, and the alpha complex is not "
                        "periodic although the configurations are. The helicity-modulus cross-check that "
                        "does locate T_BKT was measured by the earlier run and is quoted in the dataset "
                        "notes, not re-derived. At the lowest temperatures there are too few vortices "
                        "(0-14) for the spread to be stable."),
                reference="topodb_runs/quantum_fluid/results/xy_vortex.json")
            g.put(f"xy/L{L}/finding", fid)
        print("XY L", L, "done")


# --------------------------------------------------------------------- GPE
def do_gpe(d, g):
    G = J("gpe_vortex.json")
    M = G["meta"]
    for ok, v in sorted(G["by_omega"].items()):
        Om = v["Omega"]
        ds = f"quantum_fluid/gpe_rotating_Omega{Om:.2f}"
        d.add_dataset(
            id=ds, domain="quantum_fluid",
            title=f"Rotating-BEC Gross-Pitaevskii ground state, Omega = {Om:.2f}",
            source=M["source"], provenance="simulation",
            n_objects=int(np.prod(v["grid"])), ambient_dim=2,
            units="harmonic-oscillator length", sha256=v["sha256"], local_path=v["file"],
            notes=(f"{v['grid'][0]}x{v['grid'][1]} complex field, dx = {v['dx']:.4f}, Lbox = "
                   f"{v['Lbox']:.1f}, mu = {v['mu']:.4f}, converged = {v['converged']}. Condensate mask: "
                   f"smoothed |psi|^2 > {M['mask_frac']} of its maximum, then 3 erosions. "
                   f"Cubical core count uses a persistence threshold of {M['pers_thresh']} in units of "
                   f"max |psi|^2."))
        # cubical run
        rk = f"gpe/{ok}/cubical"
        if not g.has(rk):
            run = d.add_run(
                dataset_id=ds, method="cubical", coeff_field=2, max_dim=1,
                params={"Omega": Om, "mask_frac": M["mask_frac"], "pers_thresh": M["pers_thresh"],
                        "outside_mask_value": 10.0, "field": "|psi|^2 / max|psi|^2"},
                preprocessing=M["cubical"],
                script="topodb_runs/quantum_fluid/scripts/compute_gpe.py",
                command=M["command"], seed=str(M["seed"]), tier="X", repo=q.REPO)
            g.put(rk, run)
            d.add_betti(run, {0: v["n_cores_cubical_H0"]},
                        expected={0: v["n_cores_phase_winding"]})
            il.add_stats(d, run, [
                {"name": "n_cores_cubical_H0", "value": float(v["n_cores_cubical_H0"])},
                {"name": "n_cores_phase_winding", "value": float(v["n_cores_phase_winding"])},
                {"name": "crosscheck_difference", "value": float(v["crosscheck_difference"])},
                {"name": "mask_area", "value": v["mask_area"]},
            ])
            s = v["site_shuffle_control"]
            d.add_control(run, "shuffle", "site-shuffle of |psi|^2 values inside the mask",
                          passed=bool(s["separates"]),
                          detail=f"real {s['real_n_cores_cubical']} cores vs shuffled "
                                 f"{s['shuffled_mean']:.1f} +- {s['shuffled_std']:.1f} over "
                                 f"{s['n_shuffles']} shuffles. NOTE this control is weak: shuffling "
                                 f"destroys smoothness, so the shuffled field has thousands of local "
                                 f"minima and is trivially distinguishable. It rules out only the most "
                                 f"trivial artefact.")
            d.add_control(run, "known_answer",
                          "cubical H0 core count vs an independently recomputed phase-winding count",
                          passed=bool(v["crosscheck_agrees"]),
                          detail=f"cubical {v['n_cores_cubical_H0']} vs winding "
                                 f"{v['n_cores_phase_winding']} (difference "
                                 f"{v['crosscheck_difference']:+d}); the cached vortex_positions_*.npy "
                                 f"were not used")
        # alpha run
        rk = f"gpe/{ok}/alpha"
        if not g.has(rk) and "observed" in v:
            run = d.add_run(
                dataset_id=ds, method="alpha", coeff_field=2, max_dim=2,
                params={"Omega": Om, "n_null": v["n_null"], "null": "uniform random inside the mask",
                        "a_ref": v["a_density_nm_units"], "betti_threshold_frac_of_a": 0.55},
                preprocessing="alpha persistence on the phase-winding vortex cores in physical units",
                script="topodb_runs/quantum_fluid/scripts/compute_gpe.py",
                command=M["command"], seed=str(M["seed"]), tier="X", repo=q.REPO)
            g.put(rk, run)
            d.add_betti(run, {0: v["betti_alpha"]["0"], 1: v["betti_alpha"]["1"]})
            for dim, kk in ((0, "top_h0_bars"), (1, "top_h1_bars")):
                if v[kk]:
                    d.add_bars(run, dim, [(b, dd) for b, dd in v[kk]], top=10)
            nm = f"{v['n_null']} uniform random point sets inside the condensate mask, same count"
            il.add_stats(d, run, [
                {"name": "n_cores", "value": float(v["n_cores_phase_winding"])},
                {"name": "a_density", "value": v["a_density_nm_units"]},
                {"name": "a_h1", "value": v.get("a_h1")},
                {"name": "a_h0", "value": v.get("a_h0")},
                {"name": "psi6_global", "value": v["observed"]["psi6"], "null_model": nm,
                 "n_null": v["n_null"], "p_value": v["p_values_rank"]["psi6"], "p_method": "rank",
                 "multiplicity": "none; p floor = %.4f" % (1 / (v["n_null"] + 1))},
                {"name": "h0_iqr_over_median", "value": v["observed"]["h0_iqr_over_median"],
                 "null_model": nm, "n_null": v["n_null"],
                 "p_value": v["p_values_rank"]["h0_iqr_over_median"], "p_method": "rank",
                 "multiplicity": "none; p floor = %.4f" % (1 / (v["n_null"] + 1))},
            ])
            d.add_control(run, "null_calibration", "uniform random inside the condensate mask",
                          passed=True,
                          detail=f"null medians psi6 {v['null_median']['psi6']:.3f}, h0 IQR/median "
                                 f"{v['null_median']['h0_iqr_over_median']:.3f}")
        if not g.has(f"gpe/{ok}/finding"):
            agree = v["crosscheck_agrees"]
            fid = d.add_finding(
                dataset_id=ds,
                claim=(f"Rotating BEC at Omega = {Om:.2f}: the cubical (sublevel-set) H0 core count and an "
                       f"independently recomputed plaquette phase-winding count DISAGREE -- "
                       f"{v['n_cores_cubical_H0']} against {v['n_cores_phase_winding']}, a difference of "
                       f"{v['crosscheck_difference']:+d} ({100*v['crosscheck_difference']/max(v['n_cores_phase_winding'],1):+.0f} "
                       f"per cent). The excess grows with rotation rate (+3, +5, +9 at Omega = 0.70, 0.80, "
                       f"0.90) and a local-minimum proxy locates it at the condensate edge: unmatched "
                       f"density minima sit at a median distance of 0.06 from the mask boundary while "
                       f"matched vortex cores sit 1.8 to 2.5 inside it."
                       if not agree else
                       f"Rotating BEC at Omega = {Om:.2f}: the cubical H0 core count agrees with the "
                       f"independently recomputed phase-winding count "
                       f"({v['n_cores_cubical_H0']} = {v['n_cores_phase_winding']})."),
                verdict="recovered" if agree else "failed", tier="X",
                caveat=("The disagreement is a property of the sublevel-set filtration on a field with a "
                        "boundary, not of the physics: |psi|^2 falls off at the condensate edge and the "
                        "mask edge therefore creates density minima that are not vortices. A phase-winding "
                        "count is the correct core counter here; the cubical count is an upper bound. "
                        "The edge-localisation evidence comes from a local-minimum proxy (127-229 minima) "
                        "rather than from the cubical bars themselves, so it indicates the mechanism "
                        "without measuring it exactly. The site-shuffle control passes but is weak."
                        if not agree else
                        "At Omega = 0 the cubical count of 1 is the condensate component itself, not a "
                        "vortex; the winding count of 0 is correct."),
                reference="topodb_runs/quantum_fluid/results/gpe_vortex.json")
            g.put(f"gpe/{ok}/finding", fid)
        print("GPE", ok, "done")


# ----------------------------------------------------------------- colloid
def do_colloid(d, g):
    C = J("colloid.json")
    M = C["meta"]
    for pk, v in sorted(C["by_phi"].items()):
        ds = f"quantum_fluid/colloid_2d_phi{pk}"
        ordered = pk == "0.89"
        d.add_dataset(
            id=ds, domain="quantum_fluid",
            title=f"Quasi-2-D colloidal dispersion of 2.8 um particles, area fraction {pk}",
            source=f"{v['url']} (Zenodo 10.5281/zenodo.14518422, CC-BY-4.0)",
            provenance="experiment",
            n_objects=v["n_particles_median"] * v["n_frames"], ambient_dim=2, units="pixels",
            sha256=v["zip_sha256"], local_path=v["clouds_npz"],
            notes=(f"{v['n_frames']} frames sampled evenly from the tracked feature list; "
                   f"{v['n_particles_median']} particles per frame over {v['area_px2']:.0f} px^2. "
                   f"Published particle diameter 2.8 um and published area fraction {pk}. "
                   f"Contact distance (1st percentile of nearest-neighbour distance) "
                   f"{v['d_contact_px']:.2f} px, which implies a pixel size of "
                   f"{1000*v['implied_pixel_size_um']:.0f} nm. Derived point clouds sha256 "
                   f"{v['clouds_sha256']}."))
        rk = f"colloid/{pk}/alpha"
        if g.has(rk):
            continue
        run = d.add_run(
            dataset_id=ds, method="alpha", coeff_field=2, max_dim=2,
            params={"phi_published": v["phi_published"], "n_frames": M["n_frames"],
                    "n_null": v["n_null"], "null": NULL_UNI,
                    "a_ref": v["a_density_px"], "betti_threshold_frac_of_a": 0.55,
                    "particle_diameter_um": M["particle_diameter_um"]},
            preprocessing="alpha persistence per frame on the tracked particle centres in pixels; "
                          "run-level numbers are medians over the frames, bars from frame 0",
            script="topodb_runs/quantum_fluid/scripts/compute_colloid.py",
            command=M["command"], seed=str(M["seed"]), tier="X", repo=q.REPO)
        g.put(rk, run)
        d.add_betti(run, {0: v["betti_median"]["0"], 1: v["betti_median"]["1"]})
        for dim, kk in ((0, "top_h0_bars"), (1, "top_h1_bars")):
            if v[kk]:
                d.add_bars(run, dim, [(b, dd) for b, dd in v[kk]], top=10)
        nm = f"{v['n_null']} {NULL_UNI}"
        il.add_stats(d, run, [
            {"name": "n_particles_median", "value": float(v["n_particles_median"])},
            {"name": "d_contact_px", "value": v["d_contact_px"]},
            {"name": "a_density_px", "value": v["a_density_px"]},
            {"name": "a_nn_px", "value": v["a_nn_px"]},
            {"name": "a_h1_px", "value": v["a_h1_px"]},
            {"name": "a_h0_px", "value": v["a_h0_px"]},
            {"name": "a_density_over_d_contact", "value": v["a_density_over_d_contact"]},
            {"name": "a_over_d_predicted_from_phi", "value": v["a_over_d_predicted_from_phi"]},
            {"name": "phi_implied_from_geometry", "value": v["phi_implied_from_geometry"]},
            {"name": "implied_pixel_size_um", "value": v["implied_pixel_size_um"]},
            {"name": "psi6_global", "value": v["observed"]["psi6"], "null_model": nm,
             "n_null": v["n_null"], "p_value": v["p_values_rank"]["psi6"], "p_method": "rank",
             "multiplicity": "none; p floor = %.4f" % (1 / (v["n_null"] + 1))},
            {"name": "h0_iqr_over_median", "value": v["observed"]["h0_iqr_over_median"],
             "null_model": nm, "n_null": v["n_null"],
             "p_value": v["p_values_rank"]["h0_iqr_over_median"], "p_method": "rank",
             "multiplicity": "none; p floor = %.4f" % (1 / (v["n_null"] + 1))},
        ])
        d.add_control(run, "null_calibration", NULL_UNI, passed=True,
                      detail=f"null medians psi6 {v['null_median']['psi6']:.3f}, h0 IQR/median "
                             f"{v['null_median']['h0_iqr_over_median']:.3f}")
        d.add_control(run, "known_answer",
                      "scale-free hard-disc relation a/d = sqrt(pi/(2 sqrt(3) phi)) at the published phi",
                      passed=False,
                      detail=f"measured a_density/d_contact = {v['a_density_over_d_contact']:.4f} against "
                             f"the predicted {v['a_over_d_predicted_from_phi']:.4f} "
                             f"({100*v['known_answer_rel_error']:+.1f} per cent); the implied area "
                             f"fraction is {v['phi_implied_from_geometry']:.3f} against the published "
                             f"{pk}. The contact distance taken as the 1st percentile of the "
                             f"nearest-neighbour distance is biased LOW by tracking noise, which biases "
                             f"a/d high; the sign and the similar size of the miss at both "
                             f"concentrations (+21.7 and +18.6 per cent) are consistent with that "
                             f"single cause rather than with a structural discrepancy.")
        d.add_finding(
            run_id=run, dataset_id=ds,
            claim=(f"The 2-D colloidal dispersion at area fraction {pk} is "
                   + (f"the most positionally regular point set in this whole sweep: the alpha H0 death "
                      f"spread is IQR/median = {v['observed']['h0_iqr_over_median']:.3f}, below the "
                      f"Re6Zr vortex lattice's best value of 0.069 and far below the uniform-random "
                      f"{v['null_median']['h0_iqr_over_median']:.3f} (rank p = "
                      f"{v['p_values_rank']['h0_iqr_over_median']:.4f}). Its global psi6 is only "
                      f"{v['observed']['psi6']:.3f}, but significantly above the random level "
                      f"{v['null_median']['psi6']:.3f} (p = {v['p_values_rank']['psi6']:.4f}): the sample "
                      f"is polycrystalline, so grains with different orientations cancel in a global "
                      f"bond-orientational average."
                      if ordered else
                      f"a disordered gas: psi6 = {v['observed']['psi6']:.3f} is indistinguishable from "
                      f"the uniform-random level {v['null_median']['psi6']:.3f} (rank p = "
                      f"{v['p_values_rank']['psi6']:.4f}), while the H0 spread "
                      f"{v['observed']['h0_iqr_over_median']:.3f} is still well below the random "
                      f"{v['null_median']['h0_iqr_over_median']:.3f} because hard discs cannot overlap. "
                      f"Excluded volume alone produces spacing regularity with no orientational order.")),
            verdict="recovered", tier="X",
            caveat=("The published scale-free check a/d = sqrt(pi/(2 sqrt(3) phi)) is missed by about "
                    "20 per cent at BOTH concentrations, in the same direction; the contact-distance "
                    "estimator of the particle diameter is biased low by tracking noise, so this is "
                    "most likely an estimator bias rather than a property of the data, but it was not "
                    "resolved. No independent pixel calibration was available in the deposit."),
            reference="topodb_runs/quantum_fluid/results/colloid.json")
        print("colloid", pk, "done")


# ---------------------------------------------------------- disorder series
def do_disorder(d, g):
    D = J("disorder_series.json")
    M = D["meta"]
    ds = "quantum_fluid/re6zr_20kOe_disorder_series"
    d.add_dataset(
        id=ds, domain="quantum_fluid",
        title="Calibrated disorder axis: the real Re6Zr 20 kOe vortex lattice degraded in steps",
        source="derived from quantum_fluid/re6zr_stm_20kOe by "
               "topodb_runs/quantum_fluid/scripts/compute_disorder.py",
        provenance="synthetic_control", n_objects=122 * 20, ambient_dim=2, units="nm",
        notes=(f"Base: the 20 images of Re6Zr vortex cores at 20 kOe, a_ref = {M['a_ref_nm']:.2f} nm, "
               f"scan {M['scan_nm']:.0f} nm. Axis A displaces every core by an isotropic Gaussian of "
               f"width sigma (count unchanged). Axis B adds a fraction f of uniformly placed extra "
               f"points (real cores unchanged). The uniform-random end of the axis: psi6 = "
               f"{D['uniform_random_reference']['psi6']:.3f}, h0 IQR/median = "
               f"{D['uniform_random_reference']['h0_iqr_over_median']:.3f}."))
    for axis, blk, pname in (("A", D["axis_A_positional_noise"], "sigma_over_a"),
                             ("B", D["axis_B_spurious_points"], "spurious_fraction")):
        for lk, r in sorted(blk.items()):
            rk = f"disorder/{axis}/{lk}"
            if g.has(rk):
                continue
            run = d.add_run(
                dataset_id=ds, method="alpha", coeff_field=2, max_dim=2,
                params={"axis": axis, pname: r[pname], "n_null": D["n_null"],
                        "null": NULL_UNI, "a_ref_nm": M["a_ref_nm"],
                        "betti_threshold_frac_of_a": 0.55, "n_images": 20},
                preprocessing=("axis A: isotropic Gaussian displacement of every core, count unchanged"
                               if axis == "A" else
                               "axis B: uniformly placed extra points added, real cores unchanged"),
                script="topodb_runs/quantum_fluid/scripts/compute_disorder.py",
                command=M["command"], seed=str(M["seed"]), tier="X", repo=q.REPO)
            g.put(rk, run)
            d.add_betti(run, {0: int(r["betti0"]), 1: int(r["betti1"])})
            if r["top_h1_bars"]:
                d.add_bars(run, 1, [(b, dd) for b, dd in r["top_h1_bars"]], top=10)
            nm = f"{D['n_null']} {NULL_UNI}"
            il.add_stats(d, run, [
                {"name": "n_points", "value": r["n"]},
                {"name": pname, "value": r[pname]},
                {"name": "a_nn_nm", "value": r["a_nn_nm"]},
                {"name": "a_h1_nm", "value": r["a_h1_nm"]},
                {"name": "a_h0_nm", "value": r["a_h0_nm"]},
                {"name": "a_density_nm", "value": r["a_density_nm"]},
                {"name": "a_h1_over_a_density", "value": r["a_h1_nm"] / r["a_density_nm"]},
                {"name": "psi6_global", "value": r["psi6"], "null_model": nm, "n_null": D["n_null"],
                 "p_value": r["psi6_p_rank"], "p_method": "rank",
                 "multiplicity": "none; p floor = %.4f" % (1 / (D["n_null"] + 1))},
                {"name": "h0_iqr_over_median", "value": r["h0_iqr_over_median"], "null_model": nm,
                 "n_null": D["n_null"], "p_value": r["h0_iqr_p_rank"], "p_method": "rank",
                 "multiplicity": "none; p floor = %.4f" % (1 / (D["n_null"] + 1))},
            ])
            d.add_control(run, "null_calibration", NULL_UNI, passed=True,
                          detail=f"random end of the axis: psi6 "
                                 f"{D['uniform_random_reference']['psi6']:.3f}, h0 IQR/median "
                                 f"{D['uniform_random_reference']['h0_iqr_over_median']:.3f}")
    if not g.has("disorder/finding"):
        o = D["degradation_order"]
        A = D["axis_A_positional_noise"]; Bx = D["axis_B_spurious_points"]
        fid = d.add_finding(
            dataset_id=ds,
            claim=(f"PSI6 DEGRADES FIRST, as pre-stated. Under positional noise psi6 reaches half way "
                   f"from its ordered value {o['ordered_values']['psi6']:.3f} to the uniform-random "
                   f"{o['random_values']['psi6']:.3f} at sigma/a = "
                   f"{o['sigma_at_which_psi6_reaches_half_way_to_random']}, while the alpha H0 spread "
                   f"reaches half way from {o['ordered_values']['h0_iqr_over_median']:.3f} to "
                   f"{o['random_values']['h0_iqr_over_median']:.3f} only at sigma/a = "
                   f"{o['sigma_at_which_h0_iqr_reaches_half_way_to_random']}, and never gets 90 per cent "
                   f"of the way even at sigma/a = 0.30. At sigma/a = 0.30 psi6 is already statistically "
                   f"indistinguishable from random (rank p = {A['0.30']['psi6_p_rank']:.4f}) while the "
                   f"H0 spread is still highly significant (p = {A['0.30']['h0_iqr_p_rank']:.4f}). This "
                   f"is the quantitative form of the binding limit inherited from the earlier run: the "
                   f"alpha H0 spread measures REGULAR SPACING, and regular spacing survives roughly twice "
                   f"as much positional disorder as orientational order does. "
                   f"AXIS B separates the two length estimators: adding spurious points leaves the local "
                   f"spacing almost untouched (a_h1 {Bx['0.00']['a_h1_nm']:.2f} -> "
                   f"{Bx['0.40']['a_h1_nm']:.2f} nm, -3 per cent at f = 0.40) while the density-derived "
                   f"spacing collapses ({Bx['0.00']['a_density_nm']:.2f} -> "
                   f"{Bx['0.40']['a_density_nm']:.2f} nm, -16 per cent), so the ratio a_h1/a_density "
                   f"rises from {Bx['0.00']['a_h1_nm']/Bx['0.00']['a_density_nm']:.3f} to "
                   f"{Bx['0.40']['a_h1_nm']/Bx['0.40']['a_density_nm']:.3f}. Positional noise alone "
                   f"cannot produce that split (a_density is constant along axis A by construction)."),
            verdict="recovered", tier="X",
            caveat=("The base cloud is real data, so the sigma = 0 point already carries the sample's own "
                    "disorder: these are relative degradation rates from an already imperfect lattice, "
                    "not absolute thresholds. Axis B reproduces the spacing-versus-density split seen in "
                    "the Step 0 analysis of the real maps, which supports the border-over-detection "
                    "diagnosis, but it does not by itself prove that the real excess points are border "
                    "artefacts. The two axes were degraded with one seed (%s); no repeat-seed spread is "
                    "recorded." % M["seed"]),
            reference="topodb_runs/quantum_fluid/results/disorder_series.json")
        g.put("disorder/finding", fid)
    print("disorder done")


def main():
    g = il.IdGuard()
    d = il.db()
    try:
        do_xy(d, g)
        do_gpe(d, g)
        do_colloid(d, g)
        do_disorder(d, g)
    finally:
        il.close(d)
    print("all done")


if __name__ == "__main__":
    main()
