#!/usr/bin/env python3
"""X2 shared library (reverse-to-zero round 2).  Registered spec: registration.json id X2.
Specification choices fixed BEFORE looking at data T (not deviations, registration silent):
  S1 grid: 4 Mpc/h cells, shape (128,160,96), origin = data bbox min - 40 Mpc/h (periodic box only for generation).
  S2 xi estimator: gridded (NGP) Landy-Szalay-equivalent  sum C_FF / (alpha^2 sum C_RR), F = D - alpha R, C by FFT
     correlation with zero padding; bins 4 Mpc/h in [8,60): 3 F bins + 10 G bins. (Pair-count cross-check in x2_paircount_check.py.)
  S3 shell-by-shell matching: EVERY catalogue (data, randoms, Poisson, lognormal, xi*3, controls) has the data's per-shell counts
     (4 Mpc/h shells in comoving distance from the observed z) -- exactly N_D in total.
  S4 nu = per-tile (delta_s - mean)/std over the tile's valid voxels; valid voxel = footprint cell with Gaussian-kernel footprint
     fraction >= 0.9; voxels outside are excluded (filtration +1e9, never inside the nu range).
  S5 tiles: 2x2, RA split at 180, Dec split at 25.
  S6 RSD: Psi_k = i k G_k/(b k^2) from the GAUSSIAN field G (before exponentiation), f = Omega_m(z_eff)^0.55, + N(0,4 Mpc/h) along LOS.
  S7 T: LOO mu/sigma^2 from the fiducial mocks for mock T; data uses all mocks.
"""
import json, sys, pathlib, hashlib
import numpy as np
import healpy as hp
from scipy import fft as sfft
from scipy.ndimage import gaussian_filter

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parents[2]
DATA = REPO / "data/real2/cosmic_web/sdss_dr17_galaxies_ra140_220_dec0_50.csv"
CELL = 4.0
SHAPE = (128, 160, 96)
PAD = (144, 176, 112)
NSIDE = 256
SHELL = 4.0
RS_EDGES = np.arange(8.0, 60.01, 4.0)          # 13 bins; 0..2 = F [8,20), 3..12 = G [20,60)
NU = np.linspace(-3, 3, 31)
RS = 20.0
Z_EFF = 0.0742
PLANCK = dict(H0=67.66, ombh2=0.02237, omch2=0.1200, mnu=0.06, omk=0.0, tau=0.0544, As=2.100e-9, ns=0.9649)
SIGMA8 = 0.811
SIGV = 4.0
SEEDS = dict(biasgrid=6000000, fiducial=7000000, poisson=7100000, x3=7200000, randoms=20260920,
             control_void=7300000, control_ring=7400000)


def radecz_to_r(z):
    from astropy.cosmology import Planck18
    return Planck18.comoving_distance(z).value * Planck18.h


def to_xyz(ra, dec, r):
    a, d = np.radians(ra), np.radians(dec)
    return np.column_stack([r * np.cos(d) * np.cos(a), r * np.cos(d) * np.sin(a), r * np.sin(d)])


def xyz_to_rad(xyz):
    r = np.sqrt((xyz ** 2).sum(1))
    dec = np.degrees(np.arcsin(np.clip(xyz[:, 2] / np.maximum(r, 1e-12), -1, 1)))
    ra = np.degrees(np.arctan2(xyz[:, 1], xyz[:, 0])) % 360.0
    return ra, dec, r


def pix_of(ra, dec):
    return hp.ang2pix(NSIDE, np.radians(90.0 - dec), np.radians(ra))


class Ctx:
    """Everything deterministic and shared: data, mask, shells, randoms grid, lag bins, valid voxels."""
    def __init__(self, cache=True):
        import pandas as pd
        df = pd.read_csv(DATA, comment="#")
        self.ra, self.dec, self.z = df.ra.values, df.dec.values, df.z.values
        self.r = radecz_to_r(self.z)
        self.N = len(self.r)
        pix = pix_of(self.ra, self.dec)
        cnt = np.bincount(pix, minlength=hp.nside2npix(NSIDE))
        occ = cnt[cnt > 0]
        self.median_occ = float(np.median(occ))
        self.maskpix = cnt >= 0.3 * self.median_occ
        self.n_maskpix = int(self.maskpix.sum())
        self.pixarea = hp.nside2pixarea(NSIDE)
        self.omega = self.n_maskpix * self.pixarea
        self.data_in_mask = int(self.maskpix[pix].sum())
        # data points: all 193536 are kept (registration: all galaxies). Those in non-mask pixels: by construction none
        # (count>=1 pixels are in the mask when 0.3*median <= 1). Assert.
        assert self.data_in_mask == self.N, (self.data_in_mask, self.N)
        self.r_lo = np.floor(self.r.min() / SHELL) * SHELL
        self.r_hi = np.ceil(self.r.max() / SHELL) * SHELL
        self.edges = np.arange(self.r_lo, self.r_hi + 1e-6, SHELL)
        self.nsh = len(self.edges) - 1
        self.shell_cnt = np.bincount(np.clip(np.searchsorted(self.edges, self.r, "right") - 1, 0, self.nsh - 1),
                                     minlength=self.nsh)
        vol = self.omega * (self.edges[1:] ** 3 - self.edges[:-1] ** 3) / 3.0
        self.prof = self.shell_cnt / vol                       # galaxies per (Mpc/h)^3 in the mask
        xyz = to_xyz(self.ra, self.dec, self.r)
        self.origin = xyz.min(0) - 40.0
        self.data_xyz = xyz
        ext = self.origin + np.array(SHAPE) * CELL
        assert (xyz.max(0) < ext - 40.0 + 1e-9).all(), (xyz.max(0), ext)
        self._grid_setup()
        self._lag_setup()
        self._randoms()
        self._valid()

    def _grid_setup(self):
        ii = [(np.arange(n) + 0.5) * CELL + o for n, o in zip(SHAPE, self.origin)]
        X, Y, Z = np.meshgrid(*ii, indexing="ij")
        r = np.sqrt(X ** 2 + Y ** 2 + Z ** 2)
        ra = np.degrees(np.arctan2(Y, X)) % 360
        dec = np.degrees(np.arcsin(np.clip(Z / np.maximum(r, 1e-9), -1, 1)))
        self.cell_r, self.cell_ra, self.cell_dec = r.astype(np.float32), ra.astype(np.float32), dec.astype(np.float32)
        inm = np.zeros(SHAPE, bool)
        for dx in (0, -2, 2):
            for dy in (0, -2, 2):
                for dz in (0, -2, 2):
                    if (dx, dy, dz) != (0, 0, 0) and not (abs(dx) + abs(dy) + abs(dz) == 6 or True):
                        continue
                    if (dx, dy, dz) != (0, 0, 0) and (dx == 0 or dy == 0 or dz == 0):
                        continue           # centre + 8 corners
                    xx, yy, zz = X + dx, Y + dy, Z + dz
                    rr = np.sqrt(xx ** 2 + yy ** 2 + zz ** 2)
                    a = np.degrees(np.arctan2(yy, xx)) % 360
                    d = np.degrees(np.arcsin(np.clip(zz / np.maximum(rr, 1e-9), -1, 1)))
                    inrect = (a >= 140) & (a <= 220) & (d >= 0) & (d <= 50)
                    p = pix_of(np.where(inrect, a, 180.0), np.where(inrect, d, 25.0))
                    inm |= inrect & self.maskpix[p]
        del X, Y, Z
        sh = np.searchsorted(self.edges, r, "right") - 1
        ok = (sh >= 0) & (sh < self.nsh)
        w = np.where(ok, self.prof[np.clip(sh, 0, self.nsh - 1)], 0.0) * inm
        self.cell_w = w.astype(np.float64)                     # expected density weight (real space) incl. footprint
        self.foot = self.cell_w > 0

    def _lag_setup(self):
        lags = [np.fft.fftfreq(n, 1.0 / n) for n in PAD]
        LX, LY, LZ = np.meshgrid(*lags, indexing="ij")
        rr = np.sqrt(LX ** 2 + LY ** 2 + LZ ** 2) * CELL
        b = np.digitize(rr, RS_EDGES) - 1
        b[(rr < RS_EDGES[0]) | (rr >= RS_EDGES[-1])] = -1
        self.lagbin = b.ravel()
        self.lagsel = self.lagbin >= 0
        self.lagbin_sel = self.lagbin[self.lagsel]

    def cgrid(self, xyz):
        idx = np.floor((xyz - self.origin) / CELL).astype(np.int64)
        ok = ((idx >= 0) & (idx < np.array(SHAPE))).all(1)
        assert ok.all()
        flat = np.ravel_multi_index(idx.T, SHAPE)
        return np.bincount(flat, minlength=int(np.prod(SHAPE))).reshape(SHAPE).astype(np.float64)

    def corr_bins(self, F):
        Fk = sfft.rfftn(F, s=PAD)
        C = sfft.irfftn(np.abs(Fk) ** 2, s=PAD).ravel()
        return np.bincount(self.lagbin_sel, weights=C[self.lagsel], minlength=len(RS_EDGES) - 1)

    def draw_shell_points(self, counts, rng, mult=1):
        """Poisson-in-volume, shell-matched: counts[s]*mult points uniform in volume in shell s, uniform on the mask."""
        rs, ras, decs = [], [], []
        for s in range(self.nsh):
            n = int(counts[s] * mult)
            if n == 0:
                continue
            lo, hi = self.edges[s], self.edges[s + 1]
            keep_ra, keep_dec = [], []
            need = n
            while need > 0:
                m = int(need * 1.3) + 50
                ra = rng.uniform(140, 220, m)
                sd = rng.uniform(0, np.sin(np.radians(50)), m)
                dec = np.degrees(np.arcsin(sd))
                ok = self.maskpix[pix_of(ra, dec)]
                keep_ra.append(ra[ok]); keep_dec.append(dec[ok])
                need -= ok.sum()
            ra = np.concatenate(keep_ra)[:n]; dec = np.concatenate(keep_dec)[:n]
            r = np.cbrt(rng.uniform(lo ** 3, hi ** 3, n))
            rs.append(r); ras.append(ra); decs.append(dec)
        r = np.concatenate(rs); ra = np.concatenate(ras); dec = np.concatenate(decs)
        return to_xyz(ra, dec, r)

    def _randoms(self):
        rng = np.random.default_rng(SEEDS["randoms"])
        rxyz = self.draw_shell_points(self.shell_cnt, rng, mult=5)
        self.NR = len(rxyz)
        self.alpha = self.N / self.NR
        self.Rg = self.cgrid(rxyz)
        del rxyz
        self.RR = self.corr_bins(self.Rg)
        self.norm = self.alpha ** 2 * self.RR

    def xi_from_grid(self, D):
        F = D - self.alpha * self.Rg
        return self.corr_bins(F) / self.norm

    def _valid(self):
        sig = RS / CELL
        gf = gaussian_filter(self.foot.astype(np.float64), sig, mode="constant")
        self.valid = self.foot & (gf >= 0.9)
        self.den = gaussian_filter(self.alpha * self.Rg, sig, mode="constant")
        ra, dec = self.cell_ra, self.cell_dec
        self.tile = (ra >= 180).astype(int) + 2 * (dec >= 25).astype(int)
        self.tile_valid = [self.valid & (self.tile == t) for t in range(4)]
        self.tile_bbox = []
        for t in range(4):
            idx = np.argwhere(self.tile_valid[t])
            self.tile_bbox.append((idx.min(0), idx.max(0) + 1))
        self.tile_nvox = [int(tv.sum()) for tv in self.tile_valid]

    def delta_s(self, D):
        F = D - self.alpha * self.Rg
        num = gaussian_filter(F, RS / CELL, mode="constant")
        return num / np.where(self.den > 0, self.den, 1.0)

    def betti_curves(self, D):
        """returns int16 array [tile(4), type(2: sub,super), dim(3), nu(31)]"""
        import gudhi
        ds = self.delta_s(D)
        out = np.zeros((4, 2, 3, len(NU)), np.int16)
        for t in range(4):
            lo, hi = self.tile_bbox[t]
            sl = tuple(slice(a, b) for a, b in zip(lo, hi))
            v = self.tile_valid[t][sl]
            f = ds[sl]
            vals = f[v]
            nu = (f - vals.mean()) / vals.std()
            for ty, sgn in enumerate((1.0, -1.0)):
                arr = np.where(v, sgn * nu, 1e9)
                cc = gudhi.CubicalComplex(top_dimensional_cells=arr)
                cc.compute_persistence(homology_coeff_field=2, min_persistence=0.0)
                for k in range(3):
                    iv = np.asarray(cc.persistence_intervals_in_dimension(k))
                    if len(iv) == 0:
                        continue
                    bth, dth = iv[:, 0], iv[:, 1]
                    grid = sgn * NU if sgn > 0 else -NU[::-1]
                    # level set {sgn*nu <= t}; for superlevel use t = -nu  => nu >= -t
                    tt = NU if sgn > 0 else -NU
                    cnt = ((bth[None, :] <= tt[:, None]) & (tt[:, None] < dth[None, :])).sum(1)
                    out[t, ty, k] = cnt
        return out


# ------------------------------------------------------------------ CAMB / lognormal
def camb_tables():
    cache = HERE / "camb_tables.npz"
    if cache.exists():
        d = np.load(cache)
        return {k: d[k] for k in d.files}
    import camb
    from astropy.cosmology import Planck18
    p = PLANCK
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=p["H0"], ombh2=p["ombh2"], omch2=p["omch2"], mnu=p["mnu"], omk=p["omk"], tau=p["tau"])
    pars.InitPower.set_params(As=p["As"], ns=p["ns"])
    pars.set_matter_power(redshifts=[0.0], kmax=15.0)
    pars.NonLinear = camb.model.NonLinear_none
    res = camb.get_results(pars)
    kh, _, pk0 = res.get_matter_power_spectrum(minkh=1e-4, maxkh=12.0, npoints=1500)
    s8 = float(res.get_sigma8()[0])
    resc = (SIGMA8 / s8) ** 2
    pars.set_matter_power(redshifts=[Z_EFF], kmax=15.0)
    res2 = camb.get_results(pars)
    kh2, _, pkz = res2.get_matter_power_spectrum(minkh=1e-4, maxkh=12.0, npoints=1500)
    pk = pkz[0] * resc
    out = dict(kh=kh2, pk=pk, sigma8_raw=np.array(s8), rescale=np.array(resc),
               f=np.array(float(Planck18.Om(Z_EFF)) ** 0.55), omz=np.array(float(Planck18.Om(Z_EFF))))
    np.savez(cache, **out)
    return out


_XI = {}
def xi_lin_table():
    if "xi" in _XI:
        return _XI["xi"]
    t = camb_tables()
    kf = np.arange(2e-4, 8.0, 1e-3 if False else 2e-4)
    pk = np.exp(np.interp(np.log(kf), np.log(t["kh"]), np.log(t["pk"])))
    taper = np.where(kf > 3.0, np.exp(-((kf - 3.0) / 1.0) ** 2), 1.0)
    pk = pk * taper
    r = np.arange(0.25, 2000.0, 0.25)
    xi = np.empty_like(r)
    for i in range(0, len(r), 400):
        rr = r[i:i + 400]
        xi[i:i + 400] = (kf[None, :] ** 2 * pk[None, :] * np.sinc(kf[None, :] * rr[:, None] / np.pi)).sum(1) * 2e-4 / (2 * np.pi ** 2)
    _XI["xi"] = (r, xi)
    return _XI["xi"]


KGRID = np.geomspace(4e-3, 1.6, 500)
def pg_table(b):
    r, xi = xi_lin_table()
    xg = np.log1p(np.maximum(b ** 2 * xi, -0.99))
    pg = np.empty_like(KGRID)
    for i, k in enumerate(KGRID):
        pg[i] = 4 * np.pi * np.trapezoid(r ** 2 * xg * np.sinc(r * k / np.pi), r)
    return np.maximum(pg, 0.0)


class Gen:
    """Per-process lognormal generator."""
    def __init__(self, ctx):
        self.c = ctx
        kx = 2 * np.pi * np.fft.fftfreq(SHAPE[0], CELL)
        ky = 2 * np.pi * np.fft.fftfreq(SHAPE[1], CELL)
        kz = 2 * np.pi * np.fft.rfftfreq(SHAPE[2], CELL)
        self.KX, self.KY, self.KZ = np.meshgrid(kx, ky, kz, indexing="ij")
        self.kmag = np.sqrt(self.KX ** 2 + self.KY ** 2 + self.KZ ** 2)
        self.k2 = np.where(self.kmag > 0, self.kmag ** 2, 1.0)
        kzf = 2 * np.pi * np.fft.fftfreq(SHAPE[2], CELL)
        Kf = np.sqrt(np.add.outer(np.add.outer(kx ** 2, ky ** 2), kzf ** 2))
        self.kfull = Kf
        self.f = float(camb_tables()["f"])
        self._pg = {}
        self.Vbox = float(np.prod(SHAPE)) * CELL ** 3
        self.cum_cache = None

    def pgrid(self, b):
        key = round(b, 6)
        if key not in self._pg:
            pg = pg_table(b)
            pkh = np.interp(self.kmag, KGRID, pg, left=pg[0], right=0.0); pkh[self.kmag == 0] = 0
            pkf = np.interp(self.kfull, KGRID, pg, left=pg[0], right=0.0); pkf[self.kfull == 0] = 0
            self._pg[key] = (pkh, float(pkf.sum() / self.Vbox))
        return self._pg[key]

    def field(self, b, seed):
        pkh, sig2 = self.pgrid(b)
        noise = np.random.RandomState(seed).standard_normal(SHAPE)
        Gk = sfft.rfftn(noise) * np.sqrt(pkh / CELL ** 3)
        G = sfft.irfftn(Gk, s=SHAPE)
        psi = [sfft.irfftn(Gk * (1j * K / self.k2) / b, s=SHAPE) for K in (self.KX, self.KY, self.KZ)]
        return G, sig2, psi

    def mock(self, b, seed, mult=2.5, xi_scale_b=None):
        """Lognormal redshift-space mock, shell-matched to the data; returns s-space xyz (N,3)."""
        c = self.c
        rng = np.random.default_rng(seed)
        G, sig2, psi = self.field(b, seed)
        rho = np.exp(G - 0.5 * sig2)
        p = (c.cell_w * rho).ravel()
        cum = np.cumsum(p)
        M = int(mult * c.N)
        xs = []
        for attempt in range(8):
            u = rng.uniform(0, cum[-1], M)
            flat = np.minimum(np.searchsorted(cum, u), len(cum) - 1)
            idx = np.array(np.unravel_index(flat, SHAPE)).T
            x = c.origin + (idx + rng.uniform(0, 1, (M, 3))) * CELL
            rn = np.sqrt((x ** 2).sum(1)); los = x / rn[:, None]
            psi_pt = np.column_stack([psi[i].ravel()[flat] for i in range(3)])
            vlos = self.f * (psi_pt * los).sum(1) + rng.normal(0, SIGV, M)
            xs.append(x + vlos[:, None] * los)
            out = self.finish(np.vstack(xs), rng)
            if self.last_short == 0:
                break
        self.last_attempts = attempt + 1
        return out

    def finish(self, x, rng):
        """apply exact pixel mask + rectangle, then thin to exact per-shell data counts"""
        c = self.c
        ra, dec, r = xyz_to_rad(x)
        inrect = (ra >= 140) & (ra <= 220) & (dec >= 0) & (dec <= 50)
        ok = inrect & c.maskpix[pix_of(np.where(inrect, ra, 180.0), np.where(inrect, dec, 25.0))]
        sh = np.searchsorted(c.edges, r, "right") - 1
        ok &= (sh >= 0) & (sh < c.nsh)
        x, sh = x[ok], sh[ok]
        order = rng.permutation(len(sh))
        x, sh = x[order], sh[order]
        srt = np.argsort(sh, kind="stable")
        x, sh = x[srt], sh[srt]
        starts = np.searchsorted(sh, np.arange(c.nsh))
        ends = np.searchsorted(sh, np.arange(c.nsh), "right")
        take = []
        short = 0
        for s in range(c.nsh):
            avail = ends[s] - starts[s]
            need = c.shell_cnt[s]
            n = min(avail, need)
            take.append(np.arange(starts[s], starts[s] + n))
            short += need - n
        out = x[np.concatenate(take)]
        self.last_short = int(short)
        return out
