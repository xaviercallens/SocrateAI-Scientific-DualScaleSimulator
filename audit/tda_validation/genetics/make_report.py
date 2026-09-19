"""Assemble report.json from the result JSONs (all numbers are read from files written by the test scripts).
Command: python make_report.py"""
import json
import os

H = os.path.dirname(os.path.abspath(__file__))


def L(name):
    return json.load(open(os.path.join(H, name)))


def r4(x):
    return None if x is None else (round(x, 4) if isinstance(x, float) else x)


def bars(b, k=3):
    return [[round(u, 4), round(v, 4)] for u, v in b[:k]]


san = L("sanity_circle.json")
t1u, t1m = L("test1_u2os_results.json"), L("test1_mesc_results.json")
t2 = L("test2_hic_results.json")
t3 = L("test3_influenza_results.json")
ph = L("posthoc_diagnostics.json")
loc = L("local_datasets_assessment.json")
man = L("data_manifest.json")


def t1_block(r, name, data):
    o = r["observed_alpha_pca3"]
    n = r["null_gene_permutation"]
    p = r["positive_check_circular_correlation"]
    ng = r["negative_control_random_noncc_genes"]
    rp = r["rips_crosscheck_full_dim"]
    b = {"test": name, "data": data,
         "n_cells_after_qc": r["info"].get("n_cells_after_qc"), "cell_cycle_genes_used": r["info"]["cc_genes_used"],
         "pipeline_function": "cosmic_web_tda_scaled.alpha_persistence (max_alpha_sq=inf) + top_bars(r_trunc explicit) on PCA-3 of z-scored cell-cycle genes",
         "pca_variance_fraction_top3": [r4(x) for x in o["pca_var_frac_top3"]],
         "observed_top_h1_bars_birth_death": bars(o["top_h1_bars_birth_death"], 5),
         "P1_over_P2": r4(o["P1_over_P2"]), "S": r4(o["S"]), "dominant_birth_over_death": r4(o["dominant_birth_over_death"]),
         "null": {"type": "gene-wise permutation", "n": n["n_null"], "p_S": r4(n["p_S"]),
                  "null_S_q50_95_99_max": [r4(x) for x in n["null_S_quantiles_50_95_99_max"]]},
         "negative_control": {"type": "random expressed non-cell-cycle genes, same count", "draw0_P1_over_P2": r4(ng["primary_draw0"]["P1_over_P2"]),
                              "draw0_p_S_vs_own_null": r4(ng["primary_draw0"]["p_S_vs_own_null"]),
                              "frac_of_20_draws_with_P1_over_P2_ge_2": ng["frac_draws_P1_over_P2_ge_2"],
                              "median_abs_rho_cc_20_draws": r4(ng["median_abs_rho_cc"])},
         "positive_check": {k: r4(v) if isinstance(v, float) else v for k, v in p.items()},
         "rips_crosscheck_full_dim": {"n_points": rp["n_points"], "top_h1_bars": bars(rp["top_h1_bars_birth_death"]),
                                      "P1_over_P2": r4(rp["P1_over_P2"]), "p_S_rips": r4(rp["p_S_rips"]), "n_null": rp["n_null"]},
         "criteria": r["criteria"], "verdict": r["verdict"]}
    if "continuity_diagnostic" in r:
        b["continuity_diagnostic_prestated_nonbinding"] = r["continuity_diagnostic"]
    return b


def t2_case(c):
    s, rp, ln = c["alpha_path_mds3"], c["rips_full_matrix"], c["linear_null"]
    return {"alpha_top_h1_bars": bars(s["top_h1_bars_birth_death"]), "alpha_P1_over_P2": r4(s["P1_over_P2"]), "alpha_S": r4(s["S"]),
            "alpha_p_S_linear_null": r4(ln["p_S_alpha"]), "alpha_null_n": ln["n_null_alpha"],
            "mds_frac_negative_eig_mass": r4(s["mds_frac_negative_eig_mass"]),
            "rips_top_h1_bars": bars(rp["top_h1_bars_birth_death"]), "rips_P1_over_P2": r4(rp["P1_over_P2"]),
            "rips_p_S_linear_null": r4(ln["p_S_rips"]), "rips_null_n": ln["n_null_rips"],
            "exponent1_sensitivity_alpha_P1_over_P2": r4(c["alpha_path_mds3_exponent1_sensitivity"]["P1_over_P2"])}


cases = t2["cases"]
seg = t3["single_segments"]
conc = t3["concatenated_8_segments"]
report = {
    "title": "Validation of the DualScaleSimulator TDA pipeline against published genetics/genomics topology",
    "date": "2026-09-19", "branch": "loop/tda-validation",
    "pipeline_under_test": {"file": "audit/reverse_zero/E5-cosmic-web-tda-scaled/cosmic_web_tda_scaled.py",
                            "sha256": t2["pipeline_file_sha256"],
                            "functions_exercised": ["alpha_persistence", "top_bars (always with explicit r_trunc)"],
                            "functions_not_exercised": ["betti_curve", "euler_curve (not needed for the pre-stated statistics)",
                                                        "cmb_tda.betti_curves_from_topology (no dataset is a scalar field on a graph; stated in expectations.json)"],
                            "import": "loaded by path with importlib; the module only runs under __main__, so main() is not executed"},
    "commands": {
        "python": "/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python (gudhi 3.13), every call wrapped as: timeout 590 prlimit --as=8589934592 -- <python> <script> <args>; OMP/OPENBLAS/MKL threads = 1",
        "python_alignment_and_hic_fetch": "/mnt/disks/disk-socrateai-local-1/venv-tdaval/bin/python (hic-straw, pyfamsa 0.7.0, biopython)",
        "runs": ["sanity_circle.py",
                 "test1_cellcycle.py --dataset u2os --n-null 1000 --n-rips-null 100 --n-neg 20 --n-perm 10000 --budget-sec 540 (re-invoked until complete; cached seeds)",
                 "test1_cellcycle.py --dataset mesc --n-null 1000 --n-rips-null 100 --n-neg 20 --n-perm 10000 --n-tri 200 --budget-sec 530",
                 "fetch_gm12878_chr1q.py (venv-tdaval)",
                 "test2_hic.py --n-null 500 --n-rips-null 50 --budget-sec 530 --cases <case> (one process per case, re-invoked until complete), then test2_hic.py --n-null 500 --n-rips-null 50 --budget-sec 580",
                 "test3_influenza.py prep (venv-tdaval); test3_influenza.py boot --segment s --budget-sec 530 for s = 1..8; test3_influenza.py final",
                 "posthoc_diagnostics.py", "assess_local_datasets.py", "write_data_manifest.py", "make_report.py"],
        "seeds": "see each script docstring: nulls use seeds 0..n-1; subsample seeds 0 (Rips cells), 20260919 (influenza genomes)"},
    "preregistration": {"expectations.json": "commit 440d371 (alone, before any computation)",
                        "expectations_addendum_ecoli.json": "commit 5363e54 (alone, after the Caulobacter result, before any E. coli computation)"},
    "sanity_known_answer": {"circle_n200": {"top_bars": bars(san["circle_n200_sigma0.05_seed0"]["alpha"]["top_h1_bars_birth_death"]),
                                            "P1_over_P2": r4(san["circle_n200_sigma0.05_seed0"]["alpha"]["P1_over_P2"]), "passed": san["passed"]},
                            "gaussian_blob_P1_over_P2": r4(san["gaussian_blob_n200_seed1_alpha"]["P1_over_P2"]),
                            "synthetic_noiseless_ring_contact_matrix_alpha_P1_over_P2": ph["hic"]["synthetic_ring_noiseless"]["alpha_mds3"]["P1_over_P2"],
                            "synthetic_noiseless_chain_contact_matrix_n_h1": len(ph["hic"]["synthetic_chain_noiseless"]["alpha_mds3"]["top_h1_bars_birth_death"])},
    "tests": {
        "test_1a_cell_cycle_U2OS_FUCCI": t1_block(t1u, "cell-cycle loop, U2OS FUCCI scRNA-seq", "GEO GSE146773 (Mahdessian et al. 2021 Nature)"),
        "test_1b_cell_cycle_mESC": t1_block(t1m, "cell-cycle loop, mESC FACS-sorted scRNA-seq", "E-MTAB-2805 (Buettner et al. 2015 Nat Biotechnol)"),
        "test_2_Caulobacter_HiC": {
            "data": "GEO GSM1120445 (Le et al. 2013 Science), 405 x 10 kb bins, iteratively corrected",
            "pipeline_function": "alpha_persistence on classical-MDS-3 of d = C^(-1/3); Rips on full d as cross-check",
            "full_ring": t2_case(cases["caulo_full"]), "ter_cut_control": t2_case(cases["caulo_cut_ter"]),
            "ori_cut_control": t2_case(cases["caulo_cut_ori"]), "external_linear_control_GM12878_chr1q": t2_case(cases["gm12878_chr1q"]),
            "criteria": t2["test2_caulobacter"], "verdict": t2["test2_caulobacter"]["verdict"],
            "deviations": t2["deviations_from_expectations"]},
        "test_2b_Ecoli_3Cseq_addendum": {
            "data": "GEO GSM2870407 (Lioy et al. 2018 Cell), 5 kb raw counts -> 10 kb, ICE",
            "pipeline_function": "same as test 2",
            "full_ring": t2_case(cases["ecoli_full"]), "ter_cut_control": t2_case(cases["ecoli_cut_ter"]),
            "external_linear_control_GM12878_chr1q": t2_case(cases["gm12878_chr1q"]),
            "clean_info": t2["inputs"]["ecoli"]["clean_info"],
            "criteria": t2["test2b_ecoli"], "verdict": t2["test2b_ecoli"]["verdict"]},
        "test_3_influenza_reassortment": {
            "data": "NCBI Influenza Virus Resource genomeset.dat / influenza.fna (avian, complete 8-segment sets), 300 genomes",
            "n_candidates": t3["prep_info"]["n_unique_concatenated"],
            "pipeline_function": "BINDING path = plain gudhi Rips on full p-distance matrix (the published method; not the pipeline's alpha function). Pipeline alpha_persistence on classical-MDS-3 reported.",
            "single_segments": {k: {"A_max_h1_pers": r4(v["A_max_h1_persistence"]), "B_n_bars_ge_0.005": v["B_n_h1_ge_0.005"],
                                    "n_h1": v["n_h1_bars"], "median_pdist": r4(v["median_pdist"]),
                                    "alpha_mds3_P1_over_P2": r4(v["alpha_path_mds3"]["P1_over_P2"])} for k, v in seg.items()},
            "concatenated": {"A_max_h1_pers": r4(conc["A_max_h1_persistence"]), "B_n_bars_ge_0.005": conc["B_n_h1_ge_0.005"],
                             "n_h1": conc["n_h1_bars"], "top_h1_bars": bars(conc["top_h1_bars_birth_death"], 5),
                             "alpha_mds3_P1_over_P2": r4(conc["alpha_path_mds3"]["P1_over_P2"])},
            "null": t3["null_pooled_single_segment_site_bootstrap"],
            "criteria": t3["criteria_rips_binding"], "verdict": t3["verdict"]},
    },
    "post_hoc_diagnostics_NOT_prestated": {
        "u2os": {k: {kk: r4(vv) if isinstance(vv, float) else vv for kk, vv in ph["u2os"][k].items() if kk in ("P1_over_P2", "dominant_birth_over_death", "p_S", "p_P1_over_P2")}
                 for k in ("alpha_pc12_z0", "alpha_pca3_whitened", "alpha_pc12_whitened_z0", "alpha_pc12_whitened_z0_null200")},
        "u2os_rips_pc12_P1_over_P2": r4(ph["u2os"]["rips_pc12"]["P1_over_P2"]),
        "u2os_subsample250_pca3_P1_over_P2_range": [r4(min(x["pca3_P1_over_P2"] for x in ph["u2os"]["subsample_250_seeds0_9"])),
                                                    r4(max(x["pca3_P1_over_P2"] for x in ph["u2os"]["subsample_250_seeds0_9"]))],
        "caulobacter_arm_juxtaposition": {"mean_contact_mirror_pairs": ph["hic"]["caulo_contact_mirror_pairs_i_Nminus1minusi_20_180"],
                                          "mean_contact_offset20": ph["hic"]["caulo_contact_mean_offset20"],
                                          "mean_contact_offset150": ph["hic"]["caulo_contact_mean_offset150"],
                                          "mds_frac_negative_eig_mass": ph["hic"]["caulo_mds_frac_negative_eig_mass"]},
        "caulobacter_top_rips_h1_generators": ph["hic"]["caulo_rips_top5_h1_generators_bins"][:3]},
    "local_datasets_assessment": {
        "verdict": "MIXED: the three PDB files (1CRN, 1UBQ, 4OBE) are byte-identical to RCSB downloads and the *_ca_coords.npy arrays equal their C-alpha coordinates exactly (real). Every other bio_datasets array and every dual_scale_datasets array has content that points to synthetic generation under real-source names; none was used as validation data.",
        "evidence": {
            "all_23_bio_files_written_within_sec": loc["bio_mtime_span_sec"],
            "pdb": {k: {kk: v.get(kk) for kk in ("identical_to_rcsb_download", "ca_npy_max_abs_diff_vs_pdb")} for k, v in loc["bio_datasets"].items() if k.endswith(".pdb")},
            "01_hic_contact_map": loc["specific_checks"]["hic_contact_vs_coords"],
            "02_scrna_expression": {"shape": loc["bio_datasets"]["02_scrna_expression.npy"]["shape"], **loc["specific_checks"]["scrna_expression"]},
            "04_spatial_visium_coords": loc["specific_checks"]["visium_coords"],
            "05_dna_methylation_beta_hist10": loc["specific_checks"]["methylation_beta_histogram_10bins"],
            "06_cell_features_col_mean_std": loc["specific_checks"]["cell_features_col_mean_std"],
            "10_metabolic_stoichiometry": loc["specific_checks"]["stoichiometry"],
            "dual_scale_domain08_seismic_magnitudes": loc["specific_checks"]["seismic_magnitudes"],
            "dual_scale_domain06_vortex_coords_unique_per_axis": loc["specific_checks"]["vortex_coords_unique_per_axis"],
            "dual_scale_domain06_pinning_strength": "scalar 1.5",
            "dual_scale_grids": "domain01 k_grid has 24 unique values on [-pi, pi]; domain05 lats are 128 values on [-60, 60]; domain10 k_path 50 values on [0, 3]"},
        "interpretation": [
            "01_hic_contact_map is an exactly monotone function of the pairwise distances of 01_hic_chromatin_coords (Spearman -1.0; log-log slope and residual std in evidence), has no zero entries and a constant diagonal equal to 10^1.2: a contact map computed from synthetic coordinates, not a measured Hi-C matrix.",
            "02_scrna_expression has no exact zeros and no integer values; measured scRNA-seq matrices are dominated by zeros.",
            "04_spatial_visium_coords is a regular 40 x 40 grid on [0, 10], not Visium spot coordinates.",
            "06_cell_features columns have mean ~0 and std ~1 with a near-normal KS statistic: standard-normal draws.",
            "10_metabolic_stoichiometry: every reaction column has exactly one -1 and one +1 (sum 0), i.e. a random graph incidence matrix, not a metabolic network.",
            "dual_scale domain08 'magnitudes' include values >= 10 (max ~36), which no recorded earthquake has; domain06 vortex coordinates lie in a single z-plane with a scalar 'pinning_strength'.",
            "05_dna_methylation_beta is not the strongly 0/1-bimodal shape typical of array beta values; this alone is weaker evidence and is reported as such."]},
    "absent": {"none_of_the_pre-stated_tests": "all three pre-stated tests (1a, 1b, 2, 3) and the addendum (2b) were run",
               "urls_failed": man["failed_or_unused_urls"]},
    "verdict_summary": {
        "test_1a_U2OS": t1u["verdict"], "test_1b_mESC": t1m["verdict"], "test_2_Caulobacter": t2["test2_caulobacter"]["verdict"],
        "test_2b_Ecoli_addendum": t2["test2b_ecoli"]["verdict"], "test_3_influenza_Rips": t3["verdict"]},
    "narrative": [
        "The pipeline's alpha_persistence/top_bars code is correct on known answers (noisy circle, noiseless synthetic ring vs chain contact matrices).",
        "Cell cycle (1a, 1b): FAIL as pre-stated. The phase angle is recovered (|rho_cc| 0.31 U2OS vs FUCCI angle, 0.44 vs FUCCI time; 0.71 mESC vs FACS labels, permutation p ~1e-4), but the loop is not a DOMINANT H1 bar (P1/P2 about 1.5-1.6) and not significant against a gene-wise permutation null. Post hoc, it stays non-dominant in PC1-PC2, after whitening, and in 250-cell subsamples, so the failure is not an embedding-dimension or density artefact of the pipeline: at this preprocessing, the published ellipse is too thick/filled to register as a persistent H1 class.",
        "Hi-C: FAIL for Caulobacter (no dominant bar; the chromosome arms are juxtaposed, as Le et al. report, and the contact distance is strongly non-Euclidean), PASS for E. coli (pre-stated addendum): P1/P2 3.6 on the alpha path, p = 0.002 (floor of 500 nulls) against a linear-polymer null, ter-cut and human chr1q controls not dominant. The Rips cross-check on E. coli has P1/P2 2.06 but its p = 0.0196 equals the floor 1/51 of the reduced 50-draw null, so it cannot meet p <= 0.01 by construction: the observed value exceeded all 50 null draws. Caveat: the E. coli PASS depends on the pre-stated exponent 1/3; with d = C^(-1) (reported sensitivity) the alpha-path P1/P2 drops to about 1.0, and classical MDS discards about 46% negative-eigenvalue mass, so the result is not robust to the distance transform.",
        "Influenza (Chan et al. 2013): PASS with plain gudhi Rips on the full p-distance matrix: concatenated genomes have max H1 persistence 0.042 and 101 bars >= 0.005 versus <= 0.012 and <= 4 for any single segment; p = 0.0025 (floor of 400 site-bootstrap draws). The pipeline's alpha path on a 3-D MDS embedding does NOT reproduce this and gives a spurious dominant loop in the single HA segment (P1/P2 about 11), a concrete warning that 3-D embeddings of non-Euclidean metrics can create loops.",
        "Overall: on the pipeline's own alpha path the functions recover the pre-stated known topology in 0 of 3 original tests (1a, 1b, 2) and in the 1 addendum test (2b, E. coli); test 3 passes only with plain gudhi Rips on the full metric, not with the pipeline's 3-D alpha path. A clean dominant loop (E. coli ring, synthetic ring, circle) is recovered; noisy biological loops (cell cycle) and a folded ring (Caulobacter) are not. The pipeline is not validated as a general detector of biological loops at these settings, and 3-D embeddings of non-Euclidean data can create spurious loops."],
}
json.dump(report, open(os.path.join(H, "report.json"), "w"), indent=1)
print(json.dumps(report["verdict_summary"]))
