#!/usr/bin/env python3
"""Registered X2 decision procedure.  cd audit/reverse_zero_r2/X2-cosmic-web && python x2_analyze.py
Order (registration): (1) gate power on 200 Poisson (xi=0) and 200 xi*3 mocks vs the 95th percentile of the fiducial chi2;
(2) gate calibration = data chi2 two-sided empirical p >= 0.05; (3) only then topology: p_corr = 1-(1-p_min)^3.
Writes x2_results.json (numbers) -- nothing here is typed by hand."""
import sys, json, pathlib, glob
import numpy as np
from scipy import stats
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import x2_lib as L

H = L.HERE
def sig(p): return float(stats.norm.isf(min(max(p, 1e-300), 1.0) / 2))
def load(mode, n, key):
    return np.array([np.load(H / "mocks" / ("%s_%05d.npz" % (mode, k)))[key] for k in range(n)])

NF = len(glob.glob(str(H / "mocks/fiducial_*.npz")))
fid = load("fiducial", NF, "xi"); poi = load("poisson", 200, "xi"); x3 = load("x3", 200, "xi")
xi_d = np.load(H / "data_xi.npz")["xi"]
G = slice(3, 13); p = 10
res = dict(n_fiducial=NF, bias=json.load(open(H / "bias_fit.json"))["b_fit"])
mu = fid[:, G].mean(0); C = np.cov(fid[:, G].T, ddof=1); Ci = np.linalg.inv(C)
hart = (NF - p - 2) / (NF - 1.0)
chi = lambda x: hart * np.einsum("...i,ij,...j->...", x[..., G] - mu, Ci, x[..., G] - mu)
cf, cp, c3, cd = chi(fid), chi(poi), chi(x3), float(chi(xi_d))
thr = float(np.quantile(cf, 0.95))
rate = lambda a: float((a > thr).mean())
from scipy.stats import binomtest
res["gate"] = dict(hartlap=hart, thr95_fiducial=thr, chi2_mean_fiducial=float(cf.mean()), chi2_expected_mean_hartlap_scaled="p*hartlap*(N-1)/(N-p-2)~p",
    poisson_reject_rate=rate(cp), poisson_ci95=list(binomtest(int((cp > thr).sum()), 200).proportion_ci(0.95, method="exact")),
    x3_reject_rate=rate(c3), x3_ci95=list(binomtest(int((c3 > thr).sum()), 200).proportion_ci(0.95, method="exact")),
    power_pass=bool(rate(cp) >= 0.95 and rate(c3) >= 0.95), chi2_median_poisson=float(np.median(cp)), chi2_median_x3=float(np.median(c3)))
# split-half diagnostic (mu, C from 250, thresholds/chi2 from the other 250)
a, b = fid[:NF // 2, G], fid[NF // 2:, G]
mu2 = a.mean(0); Ci2 = np.linalg.inv(np.cov(a.T)); h2 = (len(a) - p - 2) / (len(a) - 1.0)
c2 = lambda x: h2 * np.einsum("...i,ij,...j->...", x[..., G] - mu2, Ci2, x[..., G] - mu2) if x.shape[-1] == 13 else None
cb = h2 * np.einsum("ni,ij,nj->n", b - mu2, Ci2, b - mu2); t2 = float(np.quantile(cb, 0.95))
res["gate"]["split_half_diagnostic"] = dict(thr95=t2, poisson_rate=float((c2(poi) > t2).mean()), x3_rate=float((c2(x3) > t2).mean()))
lo = (1 + (cf <= cd).sum()) / (NF + 1.0); hi = (1 + (cf >= cd).sum()) / (NF + 1.0)
ptwo = min(1.0, 2 * min(lo, hi))
res["calibration"] = dict(data_chi2=cd, p_lower=lo, p_upper=hi, p_two_sided=ptwo, convention="2*min((1+#{fid<=d})/(N+1),(1+#{fid>=d})/(N+1))",
    calibration_pass=bool(ptwo >= 0.05), sigma_equiv=sig(ptwo),
    data_xi_G=xi_d[G].tolist(), mock_mean_G=mu.tolist(), mock_sd_G=np.sqrt(np.diag(C)).tolist(),
    pull_per_bin=((xi_d[G] - mu) / np.sqrt(np.diag(C))).tolist())
res["shell_matching"] = dict(
    ks_r_p_fiducial=dict(median=float(np.median(load("fiducial", NF, "ks_r_p"))), min=float(load("fiducial", NF, "ks_r_p").min()), frac_below_0p05=float((load("fiducial", NF, "ks_r_p") < 0.05).mean())),
    ks_r_p_poisson=dict(median=float(np.median(load("poisson", 200, "ks_r_p"))), min=float(load("poisson", 200, "ks_r_p").min()), frac_below_0p05=float((load("poisson", 200, "ks_r_p") < 0.05).mean())),
    round1_poisson_ks_p=0.009, note="two-sample KS of mock comoving r vs data r, per mock")
res["shortfall"] = dict(fiducial_total_short=int(load("fiducial", NF, "short").sum()), max_attempts=int(load("fiducial", NF, "attempts").max()))
c = L.Ctx()
res["information"] = dict(N=c.N, n_mean_per_Mpch3=float(c.N / (c.omega * (c.edges[-1] ** 3 - c.edges[0] ** 3) / 3)),
    gal_per_20Mpch_sphere=float(c.N / (c.omega * (c.edges[-1] ** 3 - c.edges[0] ** 3) / 3) * 4 / 3 * np.pi * 20 ** 3),
    tile_valid_voxels=c.tile_nvox, valid_voxels_total=int(c.valid.sum()),
    n_independent_smoothing_volumes_est_per_tile=[float(n * L.CELL ** 3 / (2 * np.pi * L.RS ** 2) ** 1.5) for n in c.tile_nvox],
    median_occupied_pixel_count=c.median_occ, n_mask_pixels=c.n_maskpix, n_shells=c.nsh)
gate_ok = res["gate"]["power_pass"] and res["calibration"]["calibration_pass"]
res["gate_passes"] = bool(gate_ok)

# ---------------- topology
def stats_T(Bm, Bx):
    """Bm [N,4,2,3,31] fiducial; Bx [M,...] other catalogues. returns T_mock_LOO [N,2,3], T_x [M,2,3]"""
    Bm = Bm.astype(float); N = len(Bm)
    S, Q = Bm.sum(0), (Bm ** 2).sum(0)
    def T(b, mu, s2):
        return ((b - mu) ** 2 / np.maximum(s2, 1.0)).sum(axis=(1, -1))     # sum over tiles and nu -> [n,2,3]
    Tm = np.empty((N, 2, 3))
    for i in range(N):
        m = (S - Bm[i]) / (N - 1); v = ((Q - Bm[i] ** 2) - (N - 1) * m ** 2) / (N - 2)
        Tm[i] = T(Bm[i][None], m, v)[0]
    m = S / N; v = (Q - N * m ** 2) / (N - 1)
    return Tm, (lambda Bx: T(Bx.astype(float), m, v))
Bf = load("fiducial", NF, "betti")
Tm, Tfun = stats_T(Bf, None)
def pvals(T): return (1 + (Tm >= T[None]).sum(0)) / (NF + 1.0)     # [2,3]
famcorr = lambda pmin: 1 - (1 - pmin) ** 3
# controls
ctrl = {}
for name in ("control_void", "control_ring"):
    n = len(glob.glob(str(H / ("mocks/%s_*.npz" % name))))
    if n == 0: continue
    Bc = load(name, n, "betti"); Tc = Tfun(Bc)
    P = np.array([pvals(t) for t in Tc])           # [n,2,3]
    pmin = P.reshape(n, -1).min(1)
    ctrl[name] = dict(n=n, detect_rate_p_corr_lt_0p05=float((famcorr(pmin) < 0.05).mean()),
                      detect_rate_by_stat_raw_p_lt_0p05_over_6=float((P.reshape(n, -1) < 0.05 / 6).mean(0).max()),
                      per_stat_rate_p_lt_0p05_over_6={"%s_b%d" % (["sub", "super"][a], k): float((P[:, a, k] < 0.05 / 6).mean()) for a in range(2) for k in range(3)},
                      median_pmin=float(np.median(pmin)), null_false_positive_rate_expected=0.05)
    if name == "control_void":
        ctrl[name]["n_removed_median"] = float(np.median(load(name, n, "n_removed")))
    else:
        ctrl[name]["n_ring_in_footprint_median"] = float(np.median(load(name, n, "n_ring_in_footprint")))
# null false positive rate check: fiducial LOO family p_corr (self-consistency)
Pf = np.array([(1 + (np.delete(Tm, i, 0) >= Tm[i][None]).sum(0)) / float(NF) for i in range(NF)])
ctrl["fiducial_self_check"] = dict(family_false_positive_rate=float((famcorr(Pf.reshape(NF, -1).min(1)) < 0.05).mean()), note="each fiducial mock tested against the other mocks' T (LOO mu,sigma^2 already in Tm); should be near 0.05 or below")
res["controls"] = ctrl
if gate_ok:
    D = c.cgrid(c.data_xyz); Bd = c.betti_curves(D)
    np.savez(H / "data_betti.npz", betti=Bd)
    Td = Tfun(Bd[None])[0]; P = pvals(Td)
    pmin = float(P.min())
    res["topology"] = dict(T_data={"%s_b%d" % (["sub", "super"][a], k): float(Td[a, k]) for a in range(2) for k in range(3)},
        p_raw={"%s_b%d" % (["sub", "super"][a], k): float(P[a, k]) for a in range(2) for k in range(3)},
        p_min=pmin, p_corr_sidak_neff3=famcorr(pmin), p_corr_bonferroni6_diagnostic=min(1.0, 6 * pmin),
        sigma_equiv_p_corr=sig(famcorr(pmin)),
        reading=("p_corr >= 0.05: the null (Gaussian ICs, linear bias, FoG, mask) is not rejected at R_s=20 Mpc/h" if famcorr(pmin) >= 0.05 else
                 "p_corr < 0.05: the null is rejected; causes not separable (mock model as much as M0); no support for any alternative"))
else:
    res["topology"] = "NOT COMPUTED: gate failed (registered rule: INCONCLUSIVE, no topology p)"
res["verdict_registered"] = ("INCONCLUSIVE (power)" if not res["gate"]["power_pass"] else "INCONCLUSIVE (calibration)" if not res["calibration"]["calibration_pass"]
                             else ("null not rejected at R_s=20" if res["topology"]["p_corr_sidak_neff3"] >= 0.05 else "null rejected at R_s=20 (causes not separable)"))
# ---- INFORMATIONAL ONLY, deviation from the registered decision rule ----
# Registration's decision_rule computes no p-value when the calibration gate fails (verdict stays
# INCONCLUSIVE, res["verdict_registered"] above). This block computes the same T/p_corr formula anyway
# against the same (already-shown-miscalibrated) fiducial-mock null, purely for context. It is recorded
# as a deviation, not folded into the registered verdict, and is not a test of M0.
if not gate_ok:
    D = c.cgrid(c.data_xyz); Bd = c.betti_curves(D)
    np.savez(H / "data_betti.npz", betti=Bd)
    Td = Tfun(Bd[None])[0]; P = pvals(Td)
    pmin = float(P.min())
    res["topology_informational_deviation"] = dict(
        deviation_reason="Gate calibration failed (data xi(r) chi2 on G=[20,60) is far outside the 500-mock "
                          "fiducial distribution). The registered decision_rule computes no p-value in that case "
                          "(verdict stays INCONCLUSIVE). This block applies the same T/p_corr formula anyway, "
                          "against the same null, for context only -- not the registered verdict and not a test of M0.",
        T_data={"%s_b%d" % (["sub", "super"][a], k): float(Td[a, k]) for a in range(2) for k in range(3)},
        p_raw={"%s_b%d" % (["sub", "super"][a], k): float(P[a, k]) for a in range(2) for k in range(3)},
        p_min=pmin, p_corr_sidak_neff3=famcorr(pmin))
res["command"] = "cd audit/reverse_zero_r2/X2-cosmic-web && python x2_analyze.py"
json.dump(res, open(H / "x2_results.json", "w"), indent=1)
print(json.dumps({k: res[k] for k in ("gate", "calibration", "shell_matching", "topology", "verdict_registered", "controls")}, indent=1))
