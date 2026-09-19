"""Per-pulsar optimal-statistic ingredients X_a (28), Z_a (28x28) for the fixed CURN (gamma=13/3, 14 freq).
Noise: efac/t2equad/ecorr per backend + red noise (30 freqs) = medians of the released noise chains (25% burn).
Timing model: full design matrix (SVD), non-marginalised (required by enterprise_extensions OS class).
Usage: x3_persr.py [--budget SECONDS] [--workers N] [names...]   (resumable: skips cached pulsars)
"""
import argparse, json, os, re, sys, time, pathlib
import multiprocessing as mp
import numpy as np
from x3_common import *


def process(name, log10a=LOG10A_CURN, cache_dir=None):
    cache_dir = cache_dir or CACHE
    out = cache_dir / f"{name}.npz"
    if out.exists():
        return name, "cached", 0.0
    t0 = time.time()
    import scipy.linalg as sl
    from enterprise.pulsar import Pulsar
    from enterprise.signals import parameter, selections, white_signals, gp_signals, signal_base, gp_priors
    from enterprise_extensions.frequentist.optimal_statistic import OptimalStatistic

    sel = json.loads((HERE / "selection.json").read_text())
    tspan = sel["tspan_common_s"]
    par, tim = find_files(name)
    import pickle
    pk = pathlib.Path("/mnt/disks/disk-socrateai-local-1/venv-pta/psr_cache") / f"psr_{name}.pkl"
    if pk.exists():
        psr = pickle.load(open(pk, "rb"))
    else:
        psr = Pulsar(str(par), str(tim), ephem=EPHEM, timing_package="pint")
        pickle.dump(psr, open(pk, "wb"))
    load_s = time.time() - t0
    m = np.array(tim_mjds(tim))
    assert len(psr.toas) == len(m), (name, len(psr.toas), len(m))
    # psr.toas are clock-corrected TDB seconds; tim MJDs are observatory UTC -> offset ~70-140 s expected
    toa_offset_max = float(psr.toas.max() - m.max() * 86400)
    toa_offset_min = float(psr.toas.min() - m.min() * 86400)
    assert abs(toa_offset_max) < 1e5 and abs(toa_offset_min) < 1e5, (toa_offset_max, toa_offset_min)

    nz = chain_median_noise(name)
    backends = sorted(set(psr.flags["f"]))
    for b in backends:  # every backend must have chain values
        for k in ("efac", "log10_ecorr", "log10_equad"):
            assert f"{name}_{b}_{k}" in nz, (name, b, k)

    selection = selections.Selection(selections.by_backend)
    efac = parameter.Constant()
    equad = parameter.Constant()
    ecorr = parameter.Constant()
    ef = white_signals.MeasurementNoise(efac=efac, log10_t2equad=equad, selection=selection)
    ec = gp_signals.EcorrBasisModel(log10_ecorr=ecorr, selection=selection)
    rn_pl = gp_priors.powerlaw(log10_A=parameter.Constant(), gamma=parameter.Constant())
    rn = gp_signals.FourierBasisGP(rn_pl, components=30, Tspan=float(psr.toas.max() - psr.toas.min()), name="red_noise")
    gw_pl = gp_priors.powerlaw(log10_A=parameter.Constant(log10a), gamma=parameter.Constant(GAMMA_CURN))
    gw = gp_signals.FourierBasisGP(gw_pl, components=NFREQ, Tspan=tspan, name="gw")
    tm = gp_signals.TimingModel(use_svd=True)
    model = ef + ec + rn + gw + tm
    pta = signal_base.PTA([model(psr)])

    params = {}
    from enterprise.signals.parameter import ConstantParameter
    cnames = sorted({p.name for sc in pta._signalcollections for sg in sc._signals for p in sg._params.values() if isinstance(p, ConstantParameter)})
    for pn in cnames:
        if pn in ("gw_log10_A", "gw_gamma") or pn.startswith(f"{name}_gw_"):
            continue
        if pn.endswith("_red_noise_log10_A"):
            params[pn] = nz[f"{name}_red_noise_log10_A"]
        elif pn.endswith("_red_noise_gamma"):
            params[pn] = nz[f"{name}_red_noise_gamma"]
        else:
            mm = re.match(rf"^{re.escape(name)}_(?:basis_ecorr_)?(.+?)_(efac|log10_t2equad|log10_ecorr)$", pn)
            assert mm, ("unparsed parameter", pn)
            b, k = mm.groups()
            k = "log10_equad" if k == "log10_t2equad" else k
            params[pn] = nz[f"{name}_{b}_{k}"]
    assert len(params) == len(cnames) - sum(1 for c in cnames if 'gw_' in c and ('log10_A' in c or 'gamma' in c and 'red_noise' not in c)), (len(params), cnames)

    pta.set_default_params(params)
    os_ = OptimalStatistic([psr], pta=pta, orf="hd", gamma_common=GAMMA_CURN)
    freqs = os_.freqs
    assert len(freqs) == 2 * NFREQ, len(freqs)
    TNr = os_.get_TNr(params=params)[0]
    TNT = os_.get_TNT(params=params)[0]
    FNr = os_.get_FNr(params=params)[0]
    FNF = os_.get_FNF(params=params)[0]
    FNT = os_.get_FNT(params=params)[0]
    phiinv = pta.get_phiinv(params, logdet=False)[0]
    Sigma = TNT + (np.diag(phiinv) if phiinv.ndim == 1 else phiinv)
    cf = sl.cho_factor(Sigma)
    X = FNr - FNT @ sl.cho_solve(cf, TNr)
    Z = FNF - FNT @ sl.cho_solve(cf, FNT.T)
    np.savez(out, X=X, Z=Z, freqs=freqs, pos=np.array(psr.pos), ntoa=len(psr.toas), tspan=tspan,
             load_seconds=load_s, toa_offset_max=toa_offset_max, toa_offset_min=toa_offset_min, total_seconds=time.time() - t0, T_cols=TNT.shape[0])
    return name, "done", time.time() - t0


def _w(args):
    n, log10a, cache_dir = args
    try:
        return process(n, log10a=log10a, cache_dir=cache_dir)
    except Exception as e:  # report but keep going
        import traceback
        return n, "ERROR " + repr(e) + traceback.format_exc()[-600:], 0.0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--log10a", type=float, default=LOG10A_CURN,
                     help="Override the fixed CURN log10 amplitude (sensitivity check).")
    ap.add_argument("--outdir", type=str, default=None,
                     help="Override output cache directory (default: ./cache).")
    ap.add_argument("names", nargs="*")
    a = ap.parse_args()
    cache_dir = pathlib.Path(a.outdir) if a.outdir else CACHE
    cache_dir.mkdir(exist_ok=True)
    sel = json.loads((HERE / "selection.json").read_text())
    names = a.names or [r["name"] for r in sorted(sel["selected"], key=lambda r: -r["ntoa"])]
    todo = [n for n in names if not (cache_dir / f"{n}.npz").exists()]
    print(len(todo), "to do of", len(names), "log10a=", a.log10a, "outdir=", cache_dir, flush=True)
    with mp.get_context("spawn").Pool(a.workers) as pool:
        for r in pool.imap_unordered(_w, [(n, a.log10a, cache_dir) for n in todo]):
            print(r, flush=True)
