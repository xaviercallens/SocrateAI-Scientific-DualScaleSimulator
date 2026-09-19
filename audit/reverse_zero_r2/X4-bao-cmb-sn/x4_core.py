"""X4 core: unified likelihood (DESI DR2 ALL_GCcomb BAO + CMB distance priors + Pantheon+ full cov).
Merges audit/reverse_zero/E1-desi-dr2/e1_desi_dr2.py and audit/reverse_zero/E3-M0/m0_model.py
with ONE SN cut function (1580 primary: zHD in (0.01,2.4], IS_CALIBRATOR==0; 1590 disclosed:
same without the IS_CALIBRATOR filter). Paths are found from this file's location.

FRAMING: M0 (Omega_Lambda = 0.68885 imported from Planck18, GR tensor sector, no symmetron) is a
hypothesis change, not a derivation from K3xT2. Nothing here fixes mu_sym, c4 or Omega_Lambda from
K3xT2.
"""
import pathlib
import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid, quad

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[2]                      # worktree root (audit/reverse_zero_r2/X4-bao-cmb-sn -> root)
D2 = ROOT / "data" / "real2" / "dark_energy"
D1 = ROOT / "data" / "real" / "dark_energy"

OMEGA_LAMBDA_M0 = 0.68885                   # imported Planck18 value (tier L)
C = 299792.458                              # km/s
TCMB = 2.7255
NEFF = 3.046
OGH2 = 3.0 / (4.0 * 31500.0) * (TCMB / 2.7) ** 4       # Omega_gamma h^2 from the paper's 3/(4 Og h^2)=31500 (T/2.7)^-4
ORH2 = OGH2 * (1.0 + 0.2271 * NEFF)                    # massless-neutrino radiation
RB_COEF = 31500.0 * (TCMB / 2.7) ** -4                 # 3/(4 Omega_gamma h^2)

# Chen, Huang & Wang 2018 (arXiv:1808.05724) base LCDM Planck TT,TE,EE+lowE
CMB_D = np.array([1.750235, 301.4707, 0.02235976])     # distance.ini: r, la, omegabh2
CMB_INVC = np.array([[94392.3971, -1360.4913, 1664517.2916],
                     [-1360.4913, 161.4349, 3671.6180],
                     [1664517.2916, 3671.6180, 79719182.5162]])   # Distance_invcov.txt
# table 68% errors and correlation (for cross-check of the inverse covariance)
CMB_SIG_TAB = np.array([0.0046, 0.0895, 0.00015])
CMB_CORR_TAB = np.array([[1.0, 0.46, -0.66], [0.46, 1.0, -0.33], [-0.66, -0.33, 1.0]])

ZG = np.linspace(0, 2.5, 25001)


# ---------------------------------------------------------------- background
def Or_paper(Om, h):
    """Paper prescription (Chen, Huang & Wang, eq. after E(z)): Omega_r = Omega_m/(1+z_eq),
    z_eq = 2.5e4 Omega_m h^2 (T/2.7K)^-4."""
    zeq = 2.5e4 * Om * h ** 2 * (TCMB / 2.7) ** -4
    return Om / (1.0 + zeq)


def om_m0(h):
    """M0: Omega_Lambda = 0.68885 frozen, flat, Omega_m = 1 - Omega_Lambda - Omega_r(paper prescription)."""
    Om = 1.0 - OMEGA_LAMBDA_M0
    for _ in range(30):
        Om = (1.0 - OMEGA_LAMBDA_M0) / (1.0 + 1.0 / (1.0 + 2.5e4 * Om * h ** 2 * (TCMB / 2.7) ** -4))
    return Om


def E2_fn(Om, Or, w0=-1.0, wa=0.0):
    Ode = 1.0 - Om - Or
    p = 3.0 * (1.0 + w0 + wa)

    def E2(z):
        z = np.asarray(z, dtype=float)
        return (Om * (1 + z) ** 3 + Or * (1 + z) ** 4
                + Ode * (1 + z) ** p * np.exp(-3.0 * wa * z / (1 + z)))
    return E2


def zstar_hs(ombh2, ommh2):
    """Hu & Sugiyama fit as coded in the paper's CosmoMC module (g1 coefficient 0.0783; the paper's
    prose has 0.0738 -- both are recorded in x4_results.json, the code value is primary)."""
    g1 = 0.0783 * ombh2 ** -0.238 / (1 + 39.5 * ombh2 ** 0.763)
    g2 = 0.560 / (1 + 21.1 * ombh2 ** 1.81)
    return 1048.0 * (1 + 0.00124 * ombh2 ** -0.738) * (1 + g1 * ommh2 ** g2)


_GLX, _GLW = np.polynomial.legendre.leggauss(48)


def _gl(f, lo, hi, panels):
    """composite Gauss-Legendre of a vectorised f on [lo, hi] with equal panels."""
    edges = np.linspace(lo, hi, panels + 1)
    tot = 0.0
    for a, b in zip(edges[:-1], edges[1:]):
        x = 0.5 * (b - a) * _GLX + 0.5 * (a + b)
        tot += 0.5 * (b - a) * np.sum(_GLW * f(x))
    return tot


def comoving_over_c_H0(E2, zmax):
    """int_0^zmax dz/E  (dimensionless, i.e. D_M H0/c), composite Gauss-Legendre in x = ln(1+z)
    (10 panels x 48 nodes; checked against scipy quad to < 1e-9 in x4_camb_check.py)."""
    g = lambda x: np.exp(x) / np.sqrt(E2(np.exp(x) - 1.0))
    return _gl(g, 0.0, np.log(1.0 + zmax), 10)


def comoving_over_c_H0_quad(E2, zmax):
    f = lambda z: 1.0 / np.sqrt(E2(z))
    a, _ = quad(f, 0, 10, epsabs=0, epsrel=1e-11, limit=200)
    g = lambda x: np.exp(x) / np.sqrt(E2(np.exp(x) - 1.0))
    b, _ = quad(g, np.log(11.0), np.log(1 + zmax), epsabs=0, epsrel=1e-11, limit=200)
    return a + b


def _rs_integrand(Or, Om, Ode, w0, wa, ombh2):
    p = 3.0 * (1.0 + w0 + wa)

    def integrand(a):
        z = 1.0 / a - 1.0
        de = Ode * a ** 4 * (1 + z) ** p * np.exp(-3.0 * wa * z / (1 + z))
        return 1.0 / (np.sqrt(Om * a + Or + de) * np.sqrt(3.0 * (1.0 + RB_COEF * ombh2 * a)))
    return integrand


def rs_over_c_H0(Or, Om, Ode, w0, wa, ombh2, zeval):
    """sound horizon r_s(z) H0/c = int_0^{a} da / (a^2 E sqrt(3(1+Rb a))), Rb a = 3 ombh2/(4 og h2) a;
    composite Gauss-Legendre (4 panels x 48 nodes)."""
    return _gl(_rs_integrand(Or, Om, Ode, w0, wa, ombh2), 0.0, 1.0 / (1 + zeval), 4)


def rs_over_c_H0_quad(Or, Om, Ode, w0, wa, ombh2, zeval):
    v, _ = quad(_rs_integrand(Or, Om, Ode, w0, wa, ombh2), 0, 1.0 / (1 + zeval), epsabs=0, epsrel=1e-11, limit=200)
    return v


def cmb_paper(Om, h, ombh2, w0=-1.0, wa=0.0):
    """R and l_A by the paper's own prescription (eqs. la, Rz, r_s; Omega_r = Omega_m/(1+z_eq); z* by
    the Hu-Sugiyama fit as coded in the paper). This is the definition under which the published
    table reproduces Planck-like parameters (see x4_camb_check.py)."""
    Or = Or_paper(Om, h)
    Ode = 1.0 - Om - Or
    E2 = E2_fn(Om, Or, w0, wa)
    zs = zstar_hs(ombh2, Om * h ** 2)
    I = comoving_over_c_H0(E2, zs)
    rs = rs_over_c_H0(Or, Om, Ode, w0, wa, ombh2, zs)
    return {"R": float(np.sqrt(Om) * I), "lA": float(np.pi * I / rs), "zstar": float(zs), "Or": Or}


_RD = {}


def rd_table():
    """CAMB r_drag(omega_b, omega_m) table (LCDM early universe; mnu=0.06, Neff=3.046, T=2.7255) built by
    x4_camb_check.build_rd_table(); loaded from data/camb_rd_table.npz."""
    if not _RD:
        from scipy.interpolate import RectBivariateSpline
        t = np.load(HERE / "data" / "camb_rd_table.npz")
        _RD["spl"] = RectBivariateSpline(np.log(t["wb"]), np.log(t["wm"]), np.log(t["rd"]), kx=3, ky=3)
        _RD["range"] = (t["wb"].min(), t["wb"].max(), t["wm"].min(), t["wm"].max())
    return _RD


def rd_camb_interp(ombh2, ommh2):
    t = rd_table()
    lo_b, hi_b, lo_m, hi_m = t["range"]
    assert lo_b <= ombh2 <= hi_b and lo_m <= ommh2 <= hi_m, (ombh2, ommh2)
    return float(np.exp(t["spl"](np.log(ombh2), np.log(ommh2))[0, 0]))


def chi2_cmb_from(R, lA, ombh2, invC=CMB_INVC):
    d = np.array([R, lA, ombh2]) - CMB_D
    return float(d @ invC @ d)


# ---------------------------------------------------------------- data
def load_bao():
    d = D2 / "desi_dr2"
    rows = [l.split() for l in open(d / "desi_gaussian_bao_ALL_GCcomb_mean.txt") if l.strip() and not l.startswith("#")]
    z = np.array([float(r[0]) for r in rows]); v = np.array([float(r[1]) for r in rows]); q = [r[2] for r in rows]
    cov = np.loadtxt(d / "desi_gaussian_bao_ALL_GCcomb_cov.txt")
    assert cov.shape == (len(v), len(v))
    return z, v, q, np.linalg.inv(cov)


def sn_cut_mask(sn, cut):
    m = (sn.zHD > 0.01) & (sn.zHD <= 2.4) & np.isfinite(sn.MU_SH0ES) & (sn.MU_SH0ES_ERR_DIAG > 0)
    if cut == 1580:
        m = m & (sn.IS_CALIBRATOR == 0)
    elif cut != 1590:
        raise ValueError(cut)
    return m.to_numpy()


def load_sn(cut):
    """ONE SN cut function. 1580: IS_CALIBRATOR==0 (primary). 1590: disclosed variant (calibrators with
    zHD>0.01 kept). Data vector and covariance are sliced with the same boolean mask."""
    sn = pd.read_csv(D1 / "pantheon_plus_sh0es.dat", sep=r"\s+")
    with open(D2 / "Pantheon+SH0ES_STAT+SYS.cov") as f:
        n = int(f.readline().strip())
    assert n == len(sn) == 1701
    cov = np.loadtxt(D2 / "Pantheon+SH0ES_STAT+SYS.cov", skiprows=1).reshape(n, n)
    m = sn_cut_mask(sn, cut)
    assert int(m.sum()) == cut, (int(m.sum()), cut)   # registered: sample size asserted before the fit
    sz = sn.zHD.to_numpy()[m]; smu = sn.MU_SH0ES.to_numpy()[m]
    Cinv = np.linalg.inv(cov[np.ix_(m, m)])
    return sz, smu, Cinv, int(m.sum())


class Likelihood:
    def __init__(self, cut):
        self.bz, self.bv, self.bq, self.BCinv = load_bao()
        self.sz, self.smu, self.SCinv, self.nsn = load_sn(cut)
        self.cut = cut
        self.isDM = np.array([q == "DM_over_rs" for q in self.bq])
        self.isDH = np.array([q == "DH_over_rs" for q in self.bq])
        self.ones = np.ones(self.nsn); self.Sd = self.SCinv @ self.ones; self.S11 = float(self.ones @ self.Sd)

    def _bao_sn_geometry(self, E):
        dc = cumulative_trapezoid(1 / E, ZG, initial=0)            # D_C H0/c
        DM = np.interp(self.bz, ZG, dc); DH = np.interp(self.bz, ZG, 1 / E)
        DV = (self.bz * DM ** 2 * DH) ** (1 / 3)
        base = np.where(self.isDM, DM, np.where(self.isDH, DH, DV))
        mu0 = 5 * np.log10((1 + self.sz) * np.interp(self.sz, ZG, dc))
        return base, mu0

    def sn_chi2(self, mu0):
        d = self.smu - mu0
        off = float((d @ self.Sd) / self.S11)
        r = d - off
        return float(r @ self.SCinv @ r), off

    def chi2_nocmb(self, Om, w0=-1.0, wa=0.0):
        """Regression mode (reproduces E1/E3): no radiation, no CMB; BAO amplitude u (r_d h) and SN
        offset are profiled analytically."""
        E = np.sqrt(E2_fn(Om, 0.0, w0, wa)(ZG))
        base, mu0 = self._bao_sn_geometry(E)
        u = (base @ self.BCinv @ self.bv) / (base @ self.BCinv @ base)
        rb = self.bv - u * base
        cb = float(rb @ self.BCinv @ rb)
        cs, off = self.sn_chi2(mu0)
        return {"chi2": cb + cs, "bao": cb, "sn": cs, "sn_offset": off}

    def chi2_full(self, Om, h, ombh2, w0=-1.0, wa=0.0, cmb_mode="paper"):
        """BAO (r_d from CAMB r_drag(omega_b, omega_m); not profiled) + CMB distance priors + SN (offset
        profiled). cmb_mode 'paper': R, l_A by the paper's own formulas (primary). 'camb': R, l_A, r_d
        recomputed directly with CAMB (theta*-based), a sensitivity variant."""
        Or = Or_paper(Om, h)
        E = np.sqrt(E2_fn(Om, Or, w0, wa)(ZG))
        base, mu0 = self._bao_sn_geometry(E)
        if cmb_mode == "paper":
            c = cmb_paper(Om, h, ombh2, w0, wa)
            rd = rd_camb_interp(ombh2, Om * h ** 2)
        else:
            c = camb_exact(Om, h, ombh2, w0, wa)
            rd = c["rd"]
        pred = base * (C / (100.0 * h)) / rd      # base = D_M H0/c, D_H H0/c or D_V H0/c
        rb = self.bv - pred
        cb = float(rb @ self.BCinv @ rb)
        cs, off = self.sn_chi2(mu0)
        cc = chi2_cmb_from(c["R"], c["lA"], ombh2)
        return {"chi2": cb + cs + cc, "bao": cb, "sn": cs, "cmb": cc, "sn_offset": off,
                "R": c["R"], "lA": c["lA"], "rd": rd, "zstar": c["zstar"]}


def camb_exact(Om, h, ombh2, w0=-1.0, wa=0.0):
    """R, l_A (theta*-based: pi*DM(z*)/r_*), r_drag, z* straight from CAMB (2.0.4). Om = total non-DE matter
    (b + cdm + the 0.06 eV neutrino), the CosmoMC convention used for R. 1 massive neutrino, mnu=0.06,
    Neff=3.046, T=2.7255 (the paper's CosmoMC settings, distance.ini)."""
    import camb
    omnu = 0.06 / 93.14
    p = camb.CAMBparams()
    p.set_cosmology(H0=100 * h, ombh2=ombh2, omch2=Om * h * h - ombh2 - omnu, mnu=0.06, nnu=3.046,
                    num_massive_neutrinos=1, TCMB=TCMB, omk=0)
    p.set_dark_energy(w=w0, wa=wa, dark_energy_model="ppf")
    r = camb.get_background(p)
    d = r.get_derived_params()
    DM = d["DAstar"] * 1000.0                      # CAMB 'DAstar' is the comoving angular-diameter distance, Gpc
    Omtot = p.omegam
    return {"R": float(np.sqrt(Omtot) * 100 * h * DM / C), "lA": float(np.pi * DM / d["rstar"]),
            "zstar": float(d["zstar"]), "rd": float(d["rdrag"]), "rstar": float(d["rstar"]), "DMstar": float(DM)}
