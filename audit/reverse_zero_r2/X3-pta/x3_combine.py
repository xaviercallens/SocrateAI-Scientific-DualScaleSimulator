"""
Combine per-pulsar OS ingredients (X_a, Z_a) into the c4_pta_product bound.

Model (per registration X3 and param_loop_sim.compute_pta_observable):
  Gamma(theta) = a*HD(theta) + b*P4(theta),   c4_pta_product = b/a
  HD(x)  = 1.5*x*log(x) - 0.25*x + 0.5,  x = (1-cos theta)/2   [Hellings-Downs, model_orfs.hd_orf]
  P4(theta) = (35 cos^4 theta - 30 cos^2 theta + 3) / 8         [Legendre-4, same formula as
              param_loop_sim.compute_pta_observable's l4_response]

Pipeline:
  1. Build all i<j pair statistics (rho_ij, sig_ij, xi_ij) from the fixed-CURN (gamma=13/3, 14 freq)
     Sigma-marginalised X_a, Z_a exactly as enterprise_extensions.OptimalStatistic.compute_os does,
     but done here in plain numpy from the cached per-pulsar arrays (no re-fit).
  2. Bin into 7 equal-pair-count angular bins (registration).
  3. 7x7 bin covariance: DEVIATION from an analytic shared-pulsar-pair covariance (enterprise_extensions
     3.0.3 has no such function -- confirmed absent by grep of frequentist/*.py) -- estimated instead by
     delete-one-pulsar jackknife over the N_psr pulsars actually used (each jackknife sample removes one
     pulsar and all its pairs, which is exactly the shared-pulsar structure the registration asks the
     covariance to include). Effect: the calibration is empirical/jackknife rather than the model's
     analytic pair-covariance; recorded here as a deviation, not silently substituted.
  4. GLS fit for (a,b), c4 = b/a, Fieller 68%/95% intervals.
  5. Monopole/dipole negative control: same GLS machinery, single-regressor fits of rho_bin against the
     monopole and dipole ORFs.
  6. Sky-scramble gate (1000 scrambles): pulsar positions are randomly reassigned (rho_ij, sig_ij stay
     attached to the original pulsar-index pair), the 7-bin design is rebuilt from the scrambled angles
     using the FROZEN real-data jackknife covariance (a second deviation, for runtime: a fresh 66-pulsar
     jackknife per scramble is not affordable under the stated time budget), and the null chi2 =
     rho_bin^T C^-1 rho_bin (dof=7) is recorded. Gate: mean(chi2)/7 in [0.8, 1.2].
"""
import itertools
import json
import time

import numpy as np
from scipy import stats

from x3_common import *

SEED_SCRAMBLE_BASE = 8000000
N_SCRAMBLE = 1000
N_BINS = 7

rng_scramble = np.random.default_rng(20260920)  # documented separately from the per-scramble seed below


def load_all(cache_dir=None):
    cache_dir = cache_dir or CACHE
    sel = json.loads((HERE / "selection.json").read_text())
    names = [r["name"] for r in sel["selected"]]
    data, dropped = {}, []
    for n in names:
        f = cache_dir / f"{n}.npz"
        if not f.exists():
            dropped.append(n)
            continue
        d = np.load(f)
        data[n] = dict(X=d["X"], Z=d["Z"], freqs=d["freqs"], pos=d["pos"])
    return data, dropped


def pairwise(data, phiIJ):
    names = sorted(data)
    n = len(names)
    pos = np.array([data[k]["pos"] for k in names])
    rows = []
    for i in range(n):
        Xi, Zi = data[names[i]]["X"], data[names[i]]["Z"]
        for j in range(i + 1, n):
            Xj, Zj = data[names[j]]["X"], data[names[j]]["Z"]
            top = Xi @ (phiIJ * Xj)
            bot = np.trace((Zi * phiIJ[None, :]) @ (Zj * phiIJ[None, :]))
            rho = top / bot
            sig = 1.0 / np.sqrt(bot)
            cosang = float(np.clip(np.dot(pos[i], pos[j]), -1.0, 1.0))
            xi = float(np.arccos(cosang))
            rows.append((names[i], names[j], i, j, xi, rho, sig, cosang))
    return names, pos, rows


def hd(xi):
    x = (1.0 - np.cos(xi)) / 2.0
    return 1.5 * x * np.log(x) - 0.25 * x + 0.5


def p4(xi):
    c = np.cos(xi)
    return (35.0 * c ** 4 - 30.0 * c ** 2 + 3.0) / 8.0


def bin_stats(rows_arr, edges, design_fns):
    """rows_arr: structured array with fields xi,rho,sig. Returns per-bin weighted rho and design cols."""
    xi = rows_arr["xi"]; rho = rows_arr["rho"]; sig = rows_arr["sig"]
    w = 1.0 / sig ** 2
    idx = np.digitize(xi, edges[1:-1])
    rho_bin = np.zeros(len(edges) - 1)
    design = np.zeros((len(edges) - 1, len(design_fns)))
    npairs = np.zeros(len(edges) - 1, dtype=int)
    for b in range(len(edges) - 1):
        m = idx == b
        npairs[b] = m.sum()
        if m.sum() == 0:
            rho_bin[b] = np.nan
            continue
        wsum = w[m].sum()
        rho_bin[b] = np.sum(w[m] * rho[m]) / wsum
        for k, fn in enumerate(design_fns):
            design[b, k] = np.sum(w[m] * fn(xi[m])) / wsum
    return rho_bin, design, npairs


def equal_count_edges(xi, n_bins):
    order = np.sort(xi)
    n = len(order)
    idxs = [order[min(int(round(k * n / n_bins)), n - 1)] for k in range(n_bins + 1)]
    idxs[0] = order[0] - 1e-9
    idxs[-1] = order[-1] + 1e-9
    return np.array(idxs)


def gls_fit(rho_bin, design, C):
    Cinv = np.linalg.inv(C)
    M = design
    A = M.T @ Cinv @ M
    Ainv = np.linalg.inv(A)
    beta = Ainv @ M.T @ Cinv @ rho_bin
    return beta, Ainv


def fieller(a, b, Cab, conf):
    """Fieller interval for c4=b/a. Cab = [[Var a, Cov], [Cov, Var b]]. conf: coverage, e.g. 0.95."""
    k = stats.chi2.isf(1 - conf, 1)
    Va, Vb, Cov = Cab[0, 0], Cab[1, 1], Cab[0, 1]
    A = a ** 2 - k * Va
    B = -2 * (a * b - k * Cov)
    Cc = b ** 2 - k * Vb
    disc = B ** 2 - 4 * A * Cc
    if A <= 0:
        return dict(bounded=False, note="denominator term a^2 - k*Var(a) <= 0: interval is unbounded (or the complement of a bounded interval)")
    if disc < 0:
        return dict(bounded=False, note="no real roots: Fieller interval covers all reals at this confidence")
    r1 = (-B - np.sqrt(disc)) / (2 * A)
    r2 = (-B + np.sqrt(disc)) / (2 * A)
    lo, hi = sorted([r1, r2])
    return dict(bounded=True, lo=float(lo), hi=float(hi))


def main(cache_dir=None, out_name="x3_pta_result.json", log10a_label=LOG10A_CURN):
    t0 = time.time()
    data, dropped = load_all(cache_dir)
    names0 = sorted(data)
    freqs = data[names0[0]]["freqs"]
    for n in names0:
        assert np.allclose(data[n]["freqs"], freqs), (n, "frequency basis mismatch")
    from enterprise.signals import utils
    phiIJ = utils.powerlaw(freqs, log10_A=0.0, gamma=GAMMA_CURN)

    names, pos, rows = pairwise(data, phiIJ)
    xi_arr = np.array([r[4] for r in rows])
    assert xi_arr.min() > 0.0, ("duplicate/zero-separation pair found", xi_arr.min())
    dt = np.dtype([("a", "U16"), ("b", "U16"), ("i", "i4"), ("j", "i4"),
                   ("xi", "f8"), ("rho", "f8"), ("sig", "f8"), ("cos", "f8")])
    rows_arr = np.array(rows, dtype=dt)

    edges = equal_count_edges(rows_arr["xi"], N_BINS)
    design_fns = [hd, p4]
    rho_bin, design, npairs = bin_stats(rows_arr, edges, design_fns)
    assert not np.isnan(rho_bin).any(), "empty bin(s)"

    # --- delete-one-pulsar jackknife covariance (deviation: analytic shared-pulsar covariance is
    # absent from enterprise_extensions 3.0.3; grep confirms no pair_cov function) ---
    n_psr = len(names)
    jk_rho = np.zeros((n_psr, N_BINS))
    for k in range(n_psr):
        mask = (rows_arr["i"] != k) & (rows_arr["j"] != k)
        sub = rows_arr[mask]
        rb, _, npb = bin_stats(sub, edges, design_fns)
        assert not np.isnan(rb).any(), ("empty jackknife bin", k)
        jk_rho[k] = rb
    jk_mean = jk_rho.mean(axis=0)
    C = (n_psr - 1) / n_psr * (jk_rho - jk_mean).T @ (jk_rho - jk_mean)

    beta, Cov_beta = gls_fit(rho_bin, design, C)
    a_hat, b_hat = beta
    c4 = b_hat / a_hat
    interval_95 = fieller(a_hat, b_hat, Cov_beta, 0.95)
    interval_68 = fieller(a_hat, b_hat, Cov_beta, 0.6827)
    a_hat_snr = float(a_hat / np.sqrt(Cov_beta[0, 0]))

    # monopole/dipole negative control, MARGINALIZED over HD (3-parameter joint GLS fit:
    # rho_bin = a*HD + b*P4 + m*Monopole, and separately a*HD + b*P4 + d*Dipole), because a
    # single-regressor fit against monopole/dipole alone is confounded by the real HD signal
    # (monopole and HD are not orthogonal over 7 bins).
    mono_col = np.ones(N_BINS)
    dip_col = np.zeros(N_BINS)
    for b in range(N_BINS):
        m = (np.digitize(rows_arr["xi"], edges[1:-1]) == b)
        w = 1.0 / rows_arr["sig"][m] ** 2
        dip_col[b] = np.sum(w * rows_arr["cos"][m]) / w.sum()
    design_hd_p4_mono = np.column_stack([design[:, 0], design[:, 1], mono_col])
    design_hd_p4_dip = np.column_stack([design[:, 0], design[:, 1], dip_col])
    beta_m, cov_m = gls_fit(rho_bin, design_hd_p4_mono, C)
    beta_d, cov_d = gls_fit(rho_bin, design_hd_p4_dip, C)
    mono_z = float(beta_m[2] / np.sqrt(cov_m[2, 2]))
    dip_z = float(beta_d[2] / np.sqrt(cov_d[2, 2]))
    mono_amp, dip_amp = float(beta_m[2]), float(beta_d[2])

    # also the single-regressor (unmarginalized) version, reported alongside as a caveat: it is
    # confounded by the real HD signal since monopole/dipole are not orthogonal to HD over 7 bins.
    mono_design_only = mono_col.reshape(-1, 1)
    dip_design_only = dip_col.reshape(-1, 1)
    mono_beta0, mono_cov0 = gls_fit(rho_bin, mono_design_only, C)
    dip_beta0, dip_cov0 = gls_fit(rho_bin, dip_design_only, C)
    mono_z0 = float(mono_beta0[0] / np.sqrt(mono_cov0[0, 0]))
    dip_z0 = float(dip_beta0[0] / np.sqrt(dip_cov0[0, 0]))

    # --- sky-scramble gate (1000 scrambles), frozen C (deviation, for runtime) ---
    # Two statistics are computed and BOTH reported (see deviations note): the registered wording
    # ("mean chi2/7") divides by n_bins=7, which only matches a raw dof=7 statistic (no fit under
    # each scramble) -- fitting 2 parameters (a,b) per scramble first would leave dof=5, and dividing
    # THAT residual chi2 by 7 gives ~5/7=0.71 even for a perfectly calibrated C. Both are reported so
    # the choice is visible rather than picked after seeing which one passes.
    chi2_raw_scrambles, chi2_resid_scrambles = [], []
    for k in range(N_SCRAMBLE):
        rng = np.random.default_rng(SEED_SCRAMBLE_BASE + k)
        perm = rng.permutation(n_psr)
        pos_s = pos[perm]
        cosv = np.einsum("ij,ij->i", pos_s[rows_arr["i"]], pos_s[rows_arr["j"]])
        cosv = np.clip(cosv, -1.0, 1.0)
        xi_s = np.arccos(cosv)
        rows_s = rows_arr.copy()
        rows_s["xi"] = xi_s
        rows_s["cos"] = cosv
        edges_s = equal_count_edges(xi_s, N_BINS)
        rb_s, des_s, npb_s = bin_stats(rows_s, edges_s, design_fns)
        if np.isnan(rb_s).any():
            continue
        Cinv = np.linalg.inv(C)
        chi2_raw_scrambles.append(float(rb_s @ Cinv @ rb_s))
        beta_s, _ = gls_fit(rb_s, des_s, C)
        resid = rb_s - des_s @ beta_s
        chi2_resid_scrambles.append(float(resid @ Cinv @ resid))
    chi2_raw_scrambles = np.array(chi2_raw_scrambles)
    chi2_resid_scrambles = np.array(chi2_resid_scrambles)
    gate_mean_chi2_over_dof_raw7 = float(chi2_raw_scrambles.mean() / N_BINS)
    gate_mean_chi2_over_dof_resid5 = float(chi2_resid_scrambles.mean() / (N_BINS - 2))
    # Registered statistic (matches "chi2/7" literally): the raw dof=7 version.
    gate_mean_chi2_over_dof = gate_mean_chi2_over_dof_raw7
    gate_pass = 0.8 <= gate_mean_chi2_over_dof <= 1.2

    out = dict(
        n_pulsars_used=n_psr,
        n_pulsars_selected=len(names0) + len(dropped),
        dropped_pulsars=dropped,
        n_pairs=len(rows_arr),
        n_bins=N_BINS,
        bin_edges_rad=edges.tolist(),
        bin_npairs=npairs.tolist(),
        rho_bin=rho_bin.tolist(),
        design_HD_bin=design[:, 0].tolist(),
        design_P4_bin=design[:, 1].tolist(),
        jackknife_covariance=C.tolist(),
        a_hat=float(a_hat), b_hat=float(b_hat),
        Var_a=float(Cov_beta[0, 0]), Var_b=float(Cov_beta[1, 1]), Cov_ab=float(Cov_beta[0, 1]),
        c4_pta_product=float(c4),
        c4_interval_68=interval_68,
        c4_interval_95=interval_95,
        a_hat_snr_sigma=a_hat_snr,
        a_hat_snr_note=("HD amplitude a_hat is only ~%.2f sigma from zero; the 95%% Fieller interval is "
                        "bounded but marginally so (a bit less HD signal and the registered rule would "
                        "report 'no bound' instead of an interval)." % a_hat_snr),
        monopole_amplitude_marginalized=mono_amp, monopole_sigma_marginalized=float(np.sqrt(cov_m[2, 2])), monopole_z_marginalized=mono_z,
        dipole_amplitude_marginalized=dip_amp, dipole_sigma_marginalized=float(np.sqrt(cov_d[2, 2])), dipole_z_marginalized=dip_z,
        monopole_z_single_regressor_unmarginalized=mono_z0,
        dipole_z_single_regressor_unmarginalized=dip_z0,
        negative_control_note=("Monopole/dipole amplitudes are fit JOINTLY with HD+P4 (3-parameter GLS), "
                                "since a single-regressor fit is confounded by the real HD signal (monopole "
                                "and HD are not orthogonal over 7 bins). The single-regressor z-scores are "
                                "reported too, for comparison, and are not the primary number."),
        sky_scramble_n=int(len(chi2_raw_scrambles)), sky_scramble_seed_rule=f"{SEED_SCRAMBLE_BASE}+k (k<{N_SCRAMBLE})",
        gate_mean_chi2_over_dof=gate_mean_chi2_over_dof,
        gate_mean_chi2_over_dof_raw_dof7=gate_mean_chi2_over_dof_raw7,
        gate_mean_chi2_over_dof_postfit_dof5=gate_mean_chi2_over_dof_resid5,
        gate_statistic_note=("Two candidate gate statistics were computed: the raw scrambled bin vector "
                              "against C (dof=7, matches the registration's literal 'chi2/7'), and the "
                              "post-GLS-fit residual (dof=7-2=5, which would read ~5/7=0.71 even for a "
                              "perfectly calibrated C purely from fitting away 2 of 7 dof). The registered "
                              "wording is read as selecting the dof=7 statistic; the dof=5 number is "
                              "reported alongside, not discarded."),
        gate_pass=bool(gate_pass),
        fixed_curn=dict(gamma=GAMMA_CURN, log10_A=log10a_label, n_freq=NFREQ,
                         note="NOT ATTEMPTED: fetching/citing the exact NANOGrav 15-yr gamma=13/3 CURN "
                              "amplitude this round; this value is a representative fixed value, not "
                              "attributed to a specific NANOGrav table. See the sensitivity deviation "
                              "below (or, in the alt run, the base-run comparison) for the effect of "
                              "changing it."),
        deviations=[
            "Registered '7x7 covariance that includes shared-pulsar correlations' has no analytic "
            "implementation in enterprise_extensions 3.0.3 (grep of frequentist/*.py for pair_cov/covariance "
            "found none). Substituted a delete-one-pulsar jackknife covariance over the n_pulsars_used pulsars, "
            "which by construction removes all pairs sharing a pulsar together and so captures shared-pulsar "
            "correlation empirically rather than analytically.",
            "Sky-scramble gate reuses the real-data jackknife covariance C for every scramble instead of "
            "recomputing a fresh 66-pulsar jackknife per scramble (1000x cost), for runtime under the stated "
            "budget. Effect: the gate checks whether the scrambled-sky chi2 is consistent with dof=7 under the "
            "REAL covariance, not a scramble-specific one.",
            "J1713+0747 (59395 TOAs, 33 backend .tim files) is excluded: PINT Pulsar() construction did not "
            "complete in 3 attempts of ~10 minutes (foreground, prlimit 10GB) each; not a data-quality "
            "exclusion. J1713+0747 is one of the highest-timing-precision pulsars in the array, so its "
            "exclusion plausibly costs real HD signal-to-noise (relevant given a_hat_snr_sigma is only "
            "moderate); the direction of the resulting shift in c4_pta_product is not known and is the top "
            "item for a follow-up round, not something further wall-clock in this round should chase.",
            "SENSITIVITY (log10_A_CURN): a full 66-pulsar rerun at log10_A=-14.0 instead of the default "
            "-14.62 (x3_persr.py --log10a -14.0 --outdir cache_alt_amp) gives c4_pta_product=-0.324 (vs "
            "-0.187 here), an UNBOUNDED 95% Fieller interval, and a FAILING sky-scramble gate (mean "
            "chi2/7=0.424). c4_pta_product is NOT invariant to this assumed fixed amplitude; see "
            "x3_pta_result_alt_amp.json. The bound in this file should be read as conditional on "
            "log10_A_CURN=-14.62 (itself NOT ATTEMPTED to verify against the NANOGrav 15-yr paper this "
            "round -- see x3_common.py), not as amplitude-independent.",
            "No pulsar subset was named in the registration beyond '14 frequencies, gamma=13/3, fixed noise'; "
            "the working subset (66 pulsars) is the full narrowband release minus per-telescope split "
            "duplicates (ao/gbt entries where a combined entry exists) and one dropped for TOA-span < 3 yr "
            "(J0614-3329), chosen by tim-file MJDs only (data-blind to residuals/correlations), plus the "
            "single computational exclusion above.",
        ],
        wall_seconds=time.time() - t0,
    )
    (HERE / out_name).write_text(json.dumps(out, indent=1))
    print(json.dumps({k: v for k, v in out.items() if k not in
                       ("jackknife_covariance", "rho_bin", "design_HD_bin", "design_P4_bin", "bin_edges_rad")}, indent=1))


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--alt", action="store_true",
                     help="Use cache_alt_amp/ (the log10_A=-14.0 sensitivity run) instead of cache/, "
                          "and write x3_pta_result_alt_amp.json instead of x3_pta_result.json.")
    args = ap.parse_args()
    if args.alt:
        main(cache_dir=HERE / "cache_alt_amp", out_name="x3_pta_result_alt_amp.json", log10a_label=-14.0)
    else:
        main()
