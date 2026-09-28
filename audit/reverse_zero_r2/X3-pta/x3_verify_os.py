"""Verify the numpy reimplementation of the pairwise OS loop (x3_combine.pairwise) against
enterprise_extensions.OptimalStatistic.compute_os on a small (3-pulsar) PTA built the same way as
x3_persr.py. Also confirms no '<param> not set!' warnings are raised (they are logger.warning, not
python warnings, so captured via logging capture)."""
import json
import logging
import re
import numpy as np
import scipy.linalg as sl

from x3_common import *

names = ["J0030+0451", "J0023+0923", "J0340+4130"]

log_records = []


class Capture(logging.Handler):
    def emit(self, record):
        log_records.append(record.getMessage())


logging.getLogger().addHandler(Capture())
logging.getLogger().setLevel(logging.WARNING)

import pickle
from enterprise.signals import parameter, selections, white_signals, gp_signals, signal_base, gp_priors
from enterprise_extensions.frequentist.optimal_statistic import OptimalStatistic

sel = json.loads((HERE / "selection.json").read_text())
tspan = sel["tspan_common_s"]

from enterprise.signals.parameter import ConstantParameter

psrs = [pickle.load(open(f"/mnt/disks/disk-socrateai-local-1/venv-pta/psr_cache/psr_{n}.pkl", "rb")) for n in names]

signal_collections = []
for name, psr in zip(names, psrs):
    selection = selections.Selection(selections.by_backend)
    ef = white_signals.MeasurementNoise(efac=parameter.Constant(), log10_t2equad=parameter.Constant(), selection=selection)
    ec = gp_signals.EcorrBasisModel(log10_ecorr=parameter.Constant(), selection=selection)
    rn_pl = gp_priors.powerlaw(log10_A=parameter.Constant(), gamma=parameter.Constant())
    rn = gp_signals.FourierBasisGP(rn_pl, components=30, Tspan=float(psr.toas.max() - psr.toas.min()), name="red_noise")
    gw_pl = gp_priors.powerlaw(log10_A=parameter.Constant(LOG10A_CURN), gamma=parameter.Constant(GAMMA_CURN))
    gw = gp_signals.FourierBasisGP(gw_pl, components=NFREQ, Tspan=tspan, name="gw")
    tm = gp_signals.TimingModel(use_svd=True)
    model = ef + ec + rn + gw + tm
    signal_collections.append(model(psr))

pta = signal_base.PTA(signal_collections)

params_all = {}
for name in names:
    nz = chain_median_noise(name)
    cnames = sorted({p.name for sc in pta._signalcollections for sg in sc._signals
                      for p in sg._params.values() if isinstance(p, ConstantParameter) and p.name.startswith(name)})
    for pn in cnames:
        if pn.startswith(f"{name}_gw_"):
            continue
        if pn.endswith("_red_noise_log10_A"):
            params_all[pn] = nz[f"{name}_red_noise_log10_A"]
        elif pn.endswith("_red_noise_gamma"):
            params_all[pn] = nz[f"{name}_red_noise_gamma"]
        else:
            mm = re.match(rf"^{re.escape(name)}_(?:basis_ecorr_)?(.+?)_(efac|log10_t2equad|log10_ecorr)$", pn)
            assert mm, pn
            b, k = mm.groups()
            k = "log10_equad" if k == "log10_t2equad" else k
            params_all[pn] = nz[f"{name}_{b}_{k}"]

pta.set_default_params(params_all)
assert not any("not set" in m for m in log_records), [m for m in log_records if "not set" in m]

os_ = OptimalStatistic(psrs, pta=pta, orf="hd", gamma_common=GAMMA_CURN)
xi_ref, rho_ref, sig_ref, OS, OS_sig = os_.compute_os(params=params_all)

# Now the x3_combine.pairwise reimplementation, using the SAME cached X_a, Z_a route as x3_persr.py
from enterprise.signals import utils as ent_utils

TNrs = os_.get_TNr(params=params_all)
TNTs = os_.get_TNT(params=params_all)
FNrs = os_.get_FNr(params=params_all)
FNFs = os_.get_FNF(params=params_all)
FNTs = os_.get_FNT(params=params_all)
phiinvs = pta.get_phiinv(params_all, logdet=False)
Xs, Zs = [], []
for TNr, TNT, FNr, FNF, FNT, phiinv in zip(TNrs, TNTs, FNrs, FNFs, FNTs, phiinvs):
    Sigma = TNT + (np.diag(phiinv) if phiinv.ndim == 1 else phiinv)
    cf = sl.cho_factor(Sigma)
    Xs.append(FNr - FNT @ sl.cho_solve(cf, TNr))
    Zs.append(FNF - FNT @ sl.cho_solve(cf, FNT.T))

freqs = os_.freqs
phiIJ = ent_utils.powerlaw(freqs, log10_A=0.0, gamma=GAMMA_CURN)
pos = [p.pos for p in psrs]
mine_xi, mine_rho, mine_sig = [], [], []
for i in range(3):
    for j in range(i + 1, 3):
        top = Xs[i] @ (phiIJ * Xs[j])
        bot = np.trace((Zs[i] * phiIJ[None, :]) @ (Zs[j] * phiIJ[None, :]))
        mine_rho.append(top / bot)
        mine_sig.append(1.0 / np.sqrt(bot))
        mine_xi.append(np.arccos(np.clip(np.dot(pos[i], pos[j]), -1, 1)))

mine_xi, mine_rho, mine_sig = map(np.array, (mine_xi, mine_rho, mine_sig))
d_xi = np.max(np.abs(mine_xi - xi_ref))
d_rho = np.max(np.abs((mine_rho - rho_ref) / rho_ref))
d_sig = np.max(np.abs((mine_sig - sig_ref) / sig_ref))
out = dict(pulsars=names, xi_ref=xi_ref.tolist(), xi_mine=mine_xi.tolist(),
           rho_ref=rho_ref.tolist(), rho_mine=mine_rho.tolist(),
           sig_ref=sig_ref.tolist(), sig_mine=mine_sig.tolist(),
           max_abs_diff_xi=float(d_xi), max_rel_diff_rho=float(d_rho), max_rel_diff_sig=float(d_sig),
           no_param_not_set_warnings=True,
           conclusion=("AGREES to <1e-8 relative" if max(d_rho, d_sig) < 1e-8 and d_xi < 1e-10
                       else "DISAGREES -- reimplementation has a bug"))
(HERE / "x3_verify_os_result.json").write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
