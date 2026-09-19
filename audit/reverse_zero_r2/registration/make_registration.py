"""Writes registration.json for reverse-to-zero round 2 (X1-X4).

Loads NO data: only os.stat (file names and sizes) on data paths. Quote assertions run git/grep on
LeanMaster and PRE_REGISTRATION (text files, not data).
Run: cd audit/reverse_zero_r2/registration && \
  /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python make_registration.py
"""
import json, os, subprocess, pathlib
from scipy.stats import norm, chi2

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parents[2]
LM = "/home/callensxavier_gmail_com/SocrateAI-Scientific-Agora-LeanMaster"
PRE = "/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/audit/PRE_REGISTRATION.md"

def sizes(paths):
    out = {}
    for p in paths:
        fp = REPO / p
        out[p] = os.stat(fp).st_size if fp.exists() else "NOT ON DISK AT REGISTRATION"
    return out

def quote_check():
    q = []
    lm = subprocess.run(["git", "-C", LM, "show", "eb791e7:docs/STREAM8_WHICH_K3.md"], capture_output=True, text=True).stdout
    s1 = "Observables: none"
    q.append({"source": "LeanMaster eb791e7 docs/STREAM8_WHICH_K3.md", "quote": s1, "found_fixed_string": s1 in lm,
              "command": f"git -C {LM} show eb791e7:docs/STREAM8_WHICH_K3.md | grep -n -F -- '{s1}'"})
    pre = open(PRE, encoding="utf-8").read()
    s2 = "DESI DR3 or final BAO combined with CMB and any one major SN compilation"
    q.append({"source": "PRE_REGISTRATION.md (P1, the cosmological constant)", "quote": s2, "found_fixed_string": s2 in pre,
              "command": f"grep -n -F -- '{s2}' {PRE}"})
    assert all(x["found_fixed_string"] for x in q)
    return q

z3, z5 = 3.0, 5.0
def sig_to_dchi2(dof, s):
    return float(chi2.isf(2 * norm.sf(s), dof))

R = {
 "schema": "reverse-to-zero round 2 registration; written BEFORE any agent loads data",
 "date": "2026-09-19",
 "repo_relative_path": "audit/reverse_zero_r2/registration/registration.json",
 "generated_by": "audit/reverse_zero_r2/registration/make_registration.py (loads no data)",
 "no_data_loaded_statement": "This script and registration.json were produced from file names, os.stat sizes, and text-file quote checks only. No FITS, CSV, covariance, catalogue or timing file was opened.",
 "amendment_rule": "registration.json is frozen at its commit. Any later change goes in registration_amendment_<n>.json with its own commit and is reported as a deviation; a threshold typed into an analysis script is not pre-registered.",
 "framing": {
   "M0": "frozen flat LCDM with the imported Planck18 Omega_Lambda = 0.68885 (tier L), GR tensor sector c4_pta_product = 0, no symmetron sector. M0 is a HYPOTHESIS CHANGE, NOT a derivation from K3 x T2. It carries profiled nuisances (r_d*h, the SN offset); H0 is not frozen (PRE_REGISTRATION R1).",
   "rule": "Never write that K3 x T2 predicts or fixes mu_sym, c4, Omega_Lambda or any TDA outcome. Verbatim, LeanMaster Stream 8 records the phrase quoted below.",
   "quote_checks": quote_check(),
   "two_P1s": "LeanMaster Stream 6 P1 is the extra-dimension radius. PRE_REGISTRATION P1 is the cosmological constant. Unqualified 'P1' is not used.",
   "wording": "'not rejected' / 'consistent with' only; never 'confirms'. Sigma from Delta chi2: scipy chi2.sf then norm.isf(p/2).",
 },
 "environment_facts_at_registration": {
   "python": "/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python",
   "importable_checked_by_ModuleNotFound_test": {"healpy": True, "camb": True, "enterprise": False, "enterprise_extensions": False, "libstempo": False, "pint": False, "pymaster(NaMaster)": False, "emcee": False},
   "consequence_X1": "NaMaster is absent: X1 uses the ITERATED-CORRECTION route (primary), not MASTER.",
   "consequence_X3": "data/real2/pulsar_timing/ is EMPTY at registration (size 0 entries) and enterprise is not installed. The NANOGrav 15-yr timing data is ABSENT on disk; MANIFEST.json (round 1) records a Zenodo record id 16051178 with a 638 MB raw tarball (recorded there, not re-verified here). X3 must fetch, or mark ABSENT after 3 attempts with URLs listed. Installation of enterprise into a venv OUTSIDE the repo is NOT ATTEMPTED at registration.",
   "data_hash_rule": "Analysis agents verify sha256 against data/real2/MANIFEST.json after this commit; hashes are not recomputed here.",
 },
 "multiplicity_across_tests": "X1-X4 are separate hypotheses about M0 with separate datasets; no cross-test correction is applied and none is claimed. Each test's correction is stated within it.",
}

# ---------------------------------------------------------------- X1
X1 = {
 "title": "CMB TDA: Betti/Euler curves of WMAP9 ILC (and Planck SMICA if fetched) vs a spectrum-matched Gaussian null",
 "tier": "X",
 "maps": {
   "primary": "data/real2/cmb/wmap_ilc_9yr_v5.fits with mask data/real2/cmb/wmap_temperature_kq85_analysis_mask_r9_9yr_v5.fits (both nside 512, HEALPix)",
   "conditional_replication": "Planck SMICA (PR3) temperature map and its common mask, IF X1 fetches them within 3 attempts per file; candidate URLs are NOT verified at registration and are to be listed by X1. If not fetched: ABSENT, list URLs tried. SMICA is a separate family; results are never pooled with WMAP and no cross-map correction is applied.",
 },
 "pipeline_common_to_data_and_every_simulation": {
   "step_1_data_alm": "unmasked full-sky ILC map at nside 512 -> healpy map2alm(lmax=1536 (=3*512), iter=3) -> truncate at lmax=384 (=3*128)",
   "step_2_work_map": "alm2map at NSIDE_WORK=128 from those alm (band-limited at 384, no additional smoothing). The nside-512 pixel window remains inside the data alm; its size at ell=384 is computed and reported by hp.pixwin(512) at run time. Simulations are built by the identical band-limited construction: synfast(C_in, nside=128, lmax=384, pixwin=False), because C_in is defined on the alm representation that already contains the data's pixel window. Requirement lmax >= 3*nside holds (384 = 3*128).",
   "step_3_mask": "KQ85 mask ud_grade to nside 128 (power mode) then binarised at >= 0.5 (unmasked=1). The mask is applied ONCE, to the working map after step 2; sims are never pre-masked.",
   "step_4_pseudo_cl": "monopole and dipole of the unmasked-pixel field removed by subtracting the unmasked mean (ell=0,1 zeroed); pseudo-C_ell = healpy anafast(masked map, lmax=384), NOT divided by fsky, on the data and on every sim identically.",
 },
 "null_model": {
   "definition": "Gaussian, statistically isotropic random fields on the sphere with full-sky spectrum C_in(ell), built by pipeline steps 2-4 (masked pseudo-C_ell of the null MATCHES the data's).",
   "calibration_of_C_in": "ITERATED CORRECTION (MASTER not used, NaMaster absent). C_in^(0)(ell) = anafast of the UNMASKED data alm spectrum smoothed with a 5-point running mean over ell, floor 0. Iterate up to 6 times: draw N_cal=300 sims from C_in^(k), compute mean pseudo-C_ell in each gate band, multiply C_in by the ratio (data band pseudo-C_ell / sim mean) interpolated linearly in log ell between band centres (constant beyond the end centres). Stop when all bands have |data - sim mean| <= 1 * sd_sim/sqrt(N_cal) or after 6 iterations. C_in is then FROZEN.",
   "gate_bands_ell": [2, 6, 12, 24, 48, 96, 144, 192, 256, 320, 385],
   "spectrum_match_gate": {
     "sims": "N_gate=1000 sims independent of the calibration sims (own seeds), from frozen C_in",
     "statistic": "z_b = (P_data,b - mean_sim,b) / sd_sim,b for each of 10 ell bands, where P_b is the band-mean pseudo-C_ell",
     "threshold": f"PASS iff max_b |z_b| <= z_gate, with z_gate = norm.isf(0.025/10) = {float(norm.isf(0.025/10)):.4f} (Bonferroni over 10 bands, two-sided 5%)",
     "if_fail": "NO p-values are computed or reported for that map; the run is reported as GATE FAILED with the z_b table. No fallback null is substituted after the fact.",
     "order": "the gate is evaluated BEFORE any TDA statistic is computed on the data",
   },
 },
 "topology_statistic": {
   "filtration": "GUDHI SimplexTree on the unmasked-pixel graph at nside 128: vertices = unmasked pixels; edges = HEALPix neighbour pairs both unmasked; triangles = 3-cliques among mutual neighbours; simplex value = max of vertex values. Homology dims 0 and 1 (b0, b1) of a 2-complex.",
   "normalisation": "each map (data or sim) is standardised by its own unmasked mean and sd: nu = (f - mean)/sd",
   "threshold_grid": "41 equally spaced nu in [-4, 4]",
   "curves": "sublevel: b_k(nu) for the filtration by f; superlevel: b_k(nu) for the filtration by -f evaluated at -nu (i.e. excursion sets {f >= nu}); k in {0,1}. euler_chi(nu) = b0(nu) - b1(nu) (a deterministic function of the two curves).",
   "scalar_statistics": {
     "family_of_4": ["sub_b0", "sub_b1", "sup_b0", "sup_b1"],
     "definition": "T = sum over the 41 nu of (b(nu) - mu(nu))^2 / sigma^2(nu), mu and sigma from the N_null null sims (sigma floor 1 to avoid division by zero at nu where the curve is degenerate)",
     "p_value": "p = (1 + #{null sims with T_sim >= T_data}) / (1 + N_null), T_sim computed with the same mu, sigma",
     "informational_not_in_family": "euler_chi curve T computed and reported uncorrected; not counted (linear function of b0 and b1 of the same filtration)",
   },
 },
 "num_simulations": {"N_null": 2000, "N_cal_per_iteration": 300, "N_gate": 1000, "N_string_per_amplitude": 200, "min_reportable_p": "1/2001"},
 "seeds": {
   "rule": "numpy.random.seed(seed) before each healpy synfast; numpy default_rng(seed) for string placement",
   "null": "2000000 + k, k = 0..1999",
   "calibration": "1000000 + 1000*iteration + k",
   "gate": "3000000 + k, k = 0..999",
   "string_injection": "4000000 + 10000*amplitude_index + k, k = 0..199 (Gaussian part drawn with the same seeds as the tag 'string_gauss' = 5000000 + same offsets)",
 },
 "effective_number_of_tests": {
   "N_eff": 2,
   "reason_fixed_now": "Four scalar statistics come from ONE map. By duality on a closed surface, b1 of a sublevel set and b0 of its complement (the superlevel set) carry the same information up to a constant and mask effects, so {sub_b0, sup_b1} and {sub_b1, sup_b0} are two near-duplicate pairs: 2 effective tests. euler_chi = b0 - b1 is excluded from the family.",
   "correction": "Sidak: p_corr = 1 - (1 - p_min)^2 with p_min the smallest of the 4 uncorrected p-values.",
   "diagnostics_not_used_for_decision": "(i) N_eff from the null correlation matrix of the four T's, (sum lambda)^2/sum lambda^2, reported; (ii) Bonferroni with 4 tests reported as the conservative variant. Neither overrides N_eff = 2.",
 },
 "cosmic_string_sensitivity": {
   "purpose": "power of THE SAME 4-statistic family, same N_eff=2 rule, same null, to a string-like signal",
   "signal": "K = 20 great-circle arcs, random pole (uniform on S2) and random start, arc length uniform in [20, 90] degrees; each arc is a unit step (+A on one side, -A on the other over a 2-degree half-width ramp) added to a Gaussian null draw at nside 128 and then smoothed with a 1-degree FWHM Gaussian beam BEFORE the mask and normalisation; amplitude A/sigma_map in {0.05, 0.1, 0.2, 0.4} where sigma_map is the rms of the unmasked Gaussian null map.",
   "physical_conversion": "Kaiser-Stebbins: dT/T = 8*pi*G*mu*v*gamma; with the convention v*gamma = 0.4 (a convention, tier X) X1 reports the G*mu equivalent of each A/sigma_map; the mapping is illustrative, not a limit.",
   "output": "detection power vs amplitude (fraction with p_corr < 0.05). Smallest amplitude reaching power >= 0.80 is reported; if none does, the report states that a non-rejection is UNINFORMATIVE about string-like steps at A/sigma <= 0.4.",
 },
 "decision_rule": {
   "alpha": 0.05,
   "if_gate_fails": "No p-values. M0 outcome: NOT TESTED.",
   "if_p_corr_ge_0.05": "Gaussian null NOT rejected: M0's Gaussian prediction is consistent with the map at the sensitivity established by the string power table. Not evidence for M0 over alternatives; no K3xT2 statement.",
   "if_p_corr_lt_0.05": "The spectrum-matched Gaussian null is rejected by that statistic for that map and pipeline. Causes not separable: non-Gaussianity, foreground/ILC residual, mask leakage, or a null-model defect. Reported as an anomaly, not as a failure of M0 and not as support for any alternative. A rejection in WMAP but not in SMICA (or vice versa) is reported as such.",
 },
 "data_files": sizes(["data/real2/cmb/wmap_ilc_9yr_v5.fits", "data/real2/cmb/wmap_temperature_kq85_analysis_mask_r9_9yr_v5.fits"]),
 "reused_round1_code": "audit/reverse_zero/E5-cmb-tda/cmb_tda.py (topology construction); its null (masked ILC anafast/fsky, lmax 2*nside) is replaced as above.",
}

# ---------------------------------------------------------------- X2
X2 = {
 "title": "Cosmic-web TDA: SDSS DR17 vs redshift-space lognormal mocks, gated by a powered full-covariance xi(r) chi2 on a disjoint r range",
 "tier": "X",
 "data": {
   "catalogue": "data/real2/cosmic_web/sdss_dr17_galaxies_ra140_220_dec0_50.csv (ra, dec, z, zErr; SpecObj class GALAXY, zWarning=0, 0.02<=z<=0.12, ra in [140,220], dec in [0,50]; row count 193536 recorded in FETCH_COMMANDS.md, not re-verified here). ALL rows are used (no subsampling).",
   "randoms": "N_R = 5 x N_D, built with the round-1 recipe (data ra,dec kept, z shuffled from data) via data/real2/cosmic_web/make_random_catalogue.py with seed 20260920 (a new seed; the existing file, seed 20260919, is not used for the primary).",
   "distance": "comoving distance from flat LCDM with Omega_m = 1 - 0.68885 = 0.31115 (M0), h = 1 units (Mpc/h)",
 },
 "footprint_mask": {
   "nside": 256,
   "definition": "pixel is in the footprint if it lies in the ra/dec box and its data galaxy count >= 0.3 x median count over box pixels having >= 1 galaxy (0.3 fixed here). Randoms and mocks use this mask.",
 },
 "mocks": {
   "type": "lognormal, redshift space",
   "linear_power": "CAMB matter P(k) at z_eff = median data z, Planck18 flat LCDM (astropy Planck18 h, Omega_m, Omega_b, Tcmb0, Neff), n_s = 0.9649, sigma8(z=0) fixed at 0.811 (tier L convention; amplitude is degenerate with the fitted bias)",
   "field": "real-space Gaussian field with P_G from b^2 P_m (lognormal mapping xi_G = ln(1 + xi_g)), box 2 x enclosing the volume of the footprint, cell 4 Mpc/h; density n_bar(z) x mask x (exp(G - var/2)); Poisson sampling to the data n(z) (smoothed data dN/dz, 10 Mpc/h shells) and footprint.",
   "kaiser": "line-of-sight displacement s = x + (v.x_hat / (aH)) x_hat with v from the linear velocity field of the SAME Gaussian field (f = Omega_m(z_eff)^0.55, beta = f/b); the observer-centred radial line of sight is used for each galaxy (not plane-parallel).",
   "velocity_dispersion": "additional Gaussian smear along the line of sight, sigma_v = 4.0 Mpc/h (about 400 km/s); stated a priori. Variants sigma_v = 2 and 6 Mpc/h are run as informational robustness and never used to choose the primary.",
   "N_mocks_production": 500,
   "bias_fit": {
     "range_F": "r in [8, 20) Mpc/h, 3 bins of width 4",
     "method": "grid b = 0.8..2.4 step 0.1; 10 mocks per grid point (seeds 6000000 + 100*i + k); linear interpolation of the mock-mean xi_0(r) over the grid; chi2 vs data with diagonal variances from the mocks at each grid point; b_hat = minimiser. Fitted ONCE on range F only; production mocks generated at b_hat.",
   },
 },
 "xi_estimator": "gridded (FFT) Landy-Szalay-type pair estimator on 4 Mpc/h voxels with weights w = 1 (no FKP), spherically averaged over the separation vector in bins of width 4 Mpc/h (monopole of the redshift-space xi, same estimator for data and mocks, randoms shared).",
 "gate": {
   "range_G": "r in [20, 60) Mpc/h, 10 bins of width 4, DISJOINT from F",
   "covariance": "full 10 x 10 covariance from the 500 production mocks, inverse with the Hartlap factor (N - p - 2)/(N - 1) with N = 500, p = 10",
   "statistic": "chi2 = (xi - xi_mock_mean)^T C^-1 (xi - xi_mock_mean) on G",
   "power_shown_first": {
     "alt_xi_zero": "200 Poisson mocks (same n(z), footprint) as the xi = 0 alternative",
     "alt_xi_times_3": "200 mocks at b = sqrt(3) b_hat (same velocity model); the realised amplitude ratio of mock-mean xi on G to the fiducial is reported",
     "criterion": "rejection rate at the empirical 95th percentile of the fiducial mocks' chi2 must be >= 0.95 for BOTH alternatives",
     "if_power_fails": "topology p-values are NOT computed; outcome INCONCLUSIVE (gate has no power)",
   },
   "calibration_criterion": "data chi2 empirical two-sided p (relative to the 500 fiducial mocks) >= 0.05; a data chi2 that itself rejects xi = 0 (data-vs-Poisson) is reported for the power section",
   "seeds": "fiducial mocks 7000000 + k; Poisson 7100000 + k; xi x 3 mocks 7200000 + k",
 },
 "topology_statistic": {
   "restricted_to": "scales the gate validates: Gaussian smoothing radius R_s = 20 Mpc/h (>= lower edge of G). No topology is computed below 20 Mpc/h.",
   "field": "smoothed density contrast delta_s = smooth(D - alpha R) / smooth(alpha R) on 4 Mpc/h voxels, Gaussian smoothing R_s = 20 Mpc/h; analysis region = voxels where smooth(alpha R) >= 0.5 x its 90th percentile; voxels outside the region are set to +large so they never enter a sublevel set (identical treatment in mocks).",
   "subvolumes": "4 non-overlapping tiles: ra in {[140,180),[180,220)} x dec in {[0,25),[25,50)}, full z range; each tile analysed separately (tile boundary treated as region edge).",
   "curves": "GUDHI CubicalComplex Betti b0, b1, b2 of sublevel sets of delta_s and of -delta_s (superlevel), on 31 nu points in [-3, 3] where nu = delta_s / sd(delta_s in the tile region)",
   "T": "sum over tiles and nu of (b - mu)^2 / max(sigma^2, 1), mu/sigma from the 500 mocks; p = (1 + #(T_mock >= T_data)) / 501",
 },
 "effective_number_of_tests": {
   "N_eff": 3,
   "reason_fixed_now": "6 statistics (b0,b1,b2 x sub/super) from one field; by Alexander duality in 3D, b_k(sub) and b_{2-k}(super) carry the same information: 3 distinct. The Euler characteristic is excluded from the family (linear in b0 - b1 + b2).",
   "correction": "Sidak: p_corr = 1 - (1 - p_min)^3, p_min over the 6 uncorrected p-values",
 },
 "decision_rule": {
   "alpha": 0.05,
   "gate_fail": "no topology verdict (INCONCLUSIVE); reason recorded",
   "gate_pass_and_p_corr_ge_0.05": "not rejected: the data topology at R_s = 20 Mpc/h is consistent with the lognormal redshift-space null (which encodes Gaussian initial conditions with linear bias, i.e. M0's frozen-LCDM assumptions).",
   "gate_pass_and_p_corr_lt_0.05": "the null is rejected at R_s = 20 Mpc/h. This is a statement about the mock model (lognormal + linear bias + FoG) at least as much as about M0; no attribution to K3xT2 or to any alternative.",
   "note": "any excess at r < 20 Mpc/h is outside the validated range and cannot enter a verdict",
 },
 "data_files": sizes(["data/real2/cosmic_web/sdss_dr17_galaxies_ra140_220_dec0_50.csv", "data/real2/cosmic_web/make_random_catalogue.py"]),
 "seeds_topology": "same fiducial mocks as the gate",
 "reused_round1_code": "audit/reverse_zero/E5-cosmic-web-tda-scaled/e5b_lognormal_and_poisson_null.py (lognormal machinery only); its point-cloud Rips statistic and undersized gate are replaced.",
}

# ---------------------------------------------------------------- X3
X3 = {
 "title": "PTA: bound on c4_pta_product from binned Gamma(theta) with covariance derived from NANOGrav 15-yr timing data",
 "tier": "X (a bound on a parametrisation; L for the HD curve and the OS method)",
 "model": "Gamma(theta) = A^2 [ HD(theta) + c4 * P4(cos theta) ] with c4 = c4_pta_product, HD the Hellings-Downs overlap reduction function, P4 the l = 4 Legendre polynomial (the form used in audit/reverse_zero/e2_pta_and_tda_extension/e2_pta_sensitivity_scan.py). c4 = 0 is the M0 value.",
 "data": {
   "files_expected": "NANOGrav 15-yr timing data (.par/.tim, and the release's noise dictionary) under data/real2/pulsar_timing/NANOGrav15yr_PulsarTiming_v2.1.0/",
   "status_at_registration": "ABSENT on disk (directory empty). Fetch candidates recorded in MANIFEST.json: Zenodo record id 16051178 (638 MB raw tarball), and github.com/nanograv/15yr_stochastic_analysis (round 1 found no ready table there). X3 must try at most 3 sources, list URLs, and mark ABSENT if none works.",
   "not_a_substitute": "data/real/pulsar_timing/nanograv_kde_freespectrum.zip is a frequency-domain free-spectrum posterior and carries no angular information; it is not used to bound c4.",
 },
 "statistic": {
   "estimator": "noise-weighted optimal statistic pair estimators rho_ab and sigma_ab for every pulsar pair, computed with enterprise/enterprise_extensions (installed in a venv outside the repo; installation NOT ATTEMPTED at registration) with fixed per-pulsar noise from the release, spectrum gamma = 13/3, first 14 frequencies 1/T..14/T (chosen a priori; 10 and 30 frequencies reported as informational variants only).",
   "binning": "7 angular bins with equal numbers of pairs (bin edges from the pair angles only; no dependence on rho).",
   "covariance": "7 x 7 covariance of the binned estimator including pair-pair correlations from shared pulsars (the Allen 2023 'Variance of the Hellings-Downs correlation' formalism; the arXiv id is to be verified by X3 at fetch time, not asserted here).",
   "fit": "GLS on binned rho_bin = a*HD_bin + b*P4_bin, with a = A^2 and b = a*c4 (linear).",
   "bound": "Fieller interval for c4 = b/a: the set of c4 with (b_hat - c4 a_hat)^2 / Var(b_hat - c4 a_hat) <= chi2.isf(0.05, 1) = 3.8415; if the set is unbounded (a not significantly nonzero), the report states 'no bound' rather than a number.",
 },
 "calibration_gates": {
   "covariance_check": "1000 sky scrambles (random permutations of pulsar sky positions, seeds 8000000 + k) run through the same estimator and binning: mean of chi2(rho_bin, C)/7 must lie in [0.8, 1.2]; otherwise the covariance is flagged and no bound is reported.",
   "hd_recovery": "the fitted a with c4 fixed at 0 is reported with its uncertainty; the OS SNR vs HD is reported.",
 },
 "num_tests": 1,
 "effective_number_of_tests": {"N_eff": 1, "reason": "one parameter, one fit; no multiplicity correction. Frequency-count variants are informational and outside the decision."},
 "decision_rule": "Report |c4| < x at 95% (x from the Fieller interval) as a bound, or 'no bound', or ABSENT. M0 (c4 = 0) is 'consistent with' the data iff 0 lies inside the interval. No statement that M0 or K3xT2 predicts c4.",
 "data_files": sizes(["data/real2/pulsar_timing/NANOGrav15yr_PulsarTiming_v2.1.0"]),
}

# ---------------------------------------------------------------- X4
X4 = {
 "title": "DESI DR2 BAO + compressed CMB (distance priors) + Pantheon+ full covariance: CPL vs Lambda (informational for PRE_REGISTRATION P1, the cosmological constant)",
 "tier": "X",
 "status_for_P1": "INFORMATIONAL / pipeline action. PRE_REGISTRATION P1 names DR3 or final BAO + CMB + SN, Euclid or Rubin; DR2 is none of them. No P1 verdict follows from X4, whatever the number.",
 "data": {
   "bao": "data/real2/dark_energy/desi_dr2/desi_gaussian_bao_ALL_GCcomb_mean.txt and ..._cov.txt (13 measurements, 13 x 13 covariance)",
   "sn": "data/real/dark_energy/pantheon_plus_sh0es.dat with data/real2/dark_energy/Pantheon+SH0ES_STAT+SYS.cov; ONE cut everywhere: zHD > 0.01, zHD <= 2.4, IS_CALIBRATOR == 0, giving 1580 SNe, the covariance sliced by the same boolean mask. The 1590 variant (no IS_CALIBRATOR filter) is DISCLOSED as a secondary run, never mixed with the primary. Sample size is asserted equal to 1580 (primary) or 1590 (variant) before the fit; otherwise the run stops.",
   "cmb": "compressed CMB distance priors (R, l_A, omega_b) with their covariance for Planck 2018 (TT,TE,EE+lowE) from Chen, Huang & Wang, 'Distance priors from Planck final release' (JCAP 2019; arXiv id 1808.05724 recalled at registration, to be verified by fetch). X4 must copy the values from the paper's table, cite table and page, and mark ABSENT (no CMB result) if the table cannot be fetched after 3 attempts. No numbers are inserted from memory. The table variant (LCDM vs wCDM) used is stated by X4; priors were derived for restricted models, so applying them to CPL is an approximation (tier X, stated).",
 },
 "models": {
   "M0": "flat LCDM, Omega_Lambda = 0.68885 frozen (Omega_m = 0.31115); free: h, omega_b",
   "Lambda": "flat LCDM; free: Omega_m, h, omega_b",
   "CPL": "flat w0waCDM; free: Omega_m, h, omega_b, w0, wa (Lambda is the nested point w0 = -1, wa = 0)",
   "nuisances": "SN absolute-magnitude offset profiled analytically; BAO uses D/r_d with r_d computed from (omega_b, omega_m) by one stated fitting formula, NOT profiled when the CMB priors are included (X4 also reports the r_d*h-profiled BAO+SN fit as a disclosed variant, since this is the M0 nuisance convention).",
   "cmb_observables": "R = sqrt(Omega_m H0^2) D_M(z*)/c, l_A = pi D_M(z*)/r_s(z*), with z* from the Hu & Sugiyama fitting formula as used in the distance-prior papers; radiation from Tcmb = 2.7255 K, Neff = 3.046; massless neutrinos (approximation stated).",
 },
 "statistics": {
   "delta_chi2_CPL_vs_Lambda": "chi2_min(Lambda) - chi2_min(CPL), 2 dof",
   "delta_chi2_Lambda_vs_M0": "chi2_min(M0) - chi2_min(Lambda), 1 dof",
   "sigma": "p = scipy.stats.chi2.sf(delta, dof); sigma = scipy.stats.norm.isf(p/2)",
   "minimiser": "scipy.optimize with 20 seeded random starts (seed 9000000 + start) per model; best chi2 kept; all starts reported",
 },
 "num_tests": 2,
 "effective_number_of_tests": {"N_eff": 1, "reason": "CPL vs Lambda is the primary comparison; Lambda vs M0 is reported alongside and is not corrected. The result is informational, not a rejection test."},
 "decision_rule": {
   "thresholds_for_context_only": {"2dof_3sigma": round(sig_to_dchi2(2, 3), 2), "2dof_5sigma": round(sig_to_dchi2(2, 5), 2), "1dof_3sigma": round(sig_to_dchi2(1, 3), 2), "1dof_5sigma": round(sig_to_dchi2(1, 5), 2), "source": "PRE_REGISTRATION.md threshold table, recomputed with scipy here; it is a table, not a rule for this comparison"},
   "reporting": "X4 reports Delta chi2, the p-value and sigma for each comparison, primary (1580) and disclosed (1590) cuts, and says 'DR2, informational for PRE_REGISTRATION P1; not a P1 verdict'. Words such as 'falsified' or 'confirmed' are not used. Because the CMB priors are an approximation, X4 also states that the published DR2 combination (PRE_REGISTRATION R4, 3.1 sigma with the full CMB likelihood) is not re-analysed here.",
   "for_M0": "M0 vs Lambda: Delta chi2 with 1 dof; M0 'not rejected' if p >= 0.05, else 'disfavoured relative to fitted Lambda at that DR2 combination', still no P1 verdict.",
 },
 "data_files": sizes(["data/real2/dark_energy/desi_dr2/desi_gaussian_bao_ALL_GCcomb_mean.txt", "data/real2/dark_energy/desi_dr2/desi_gaussian_bao_ALL_GCcomb_cov.txt", "data/real/dark_energy/pantheon_plus_sh0es.dat", "data/real2/dark_energy/Pantheon+SH0ES_STAT+SYS.cov"]),
}

R["tests"] = {"X1": X1, "X2": X2, "X3": X3, "X4": X4}
(HERE / "registration.json").write_text(json.dumps(R, indent=1))
print("written", HERE / "registration.json")
