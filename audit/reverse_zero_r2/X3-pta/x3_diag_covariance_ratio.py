import json, numpy as np
import sys
sys.path.insert(0, '.')
from x3_combine import load_all, pairwise, equal_count_edges, bin_stats, hd, p4, N_BINS
from x3_common import GAMMA_CURN
from enterprise.signals import utils

def run(cache_dir, label):
    data, dropped = load_all(cache_dir)
    names0 = sorted(data)
    freqs = data[names0[0]]["freqs"]
    phiIJ = utils.powerlaw(freqs, log10_A=0.0, gamma=GAMMA_CURN)
    names, pos, rows = pairwise(data, phiIJ)
    dt = np.dtype([("a","U16"),("b","U16"),("i","i4"),("j","i4"),("xi","f8"),("rho","f8"),("sig","f8"),("cos","f8")])
    rows_arr = np.array(rows, dtype=dt)
    edges = equal_count_edges(rows_arr["xi"], N_BINS)
    design_fns=[hd,p4]
    rho_bin, design, npairs = bin_stats(rows_arr, edges, design_fns)
    n_psr=len(names)
    jk_rho=np.zeros((n_psr,N_BINS))
    for k in range(n_psr):
        mask=(rows_arr["i"]!=k)&(rows_arr["j"]!=k)
        rb,_,_=bin_stats(rows_arr[mask], edges, design_fns)
        jk_rho[k]=rb
    jk_mean=jk_rho.mean(axis=0)
    C=(n_psr-1)/n_psr*(jk_rho-jk_mean).T@(jk_rho-jk_mean)

    pos_arr=pos
    rb_scrambles=[]
    for k in range(1000):
        rng=np.random.default_rng(8000000+k)
        perm=rng.permutation(n_psr)
        pos_s=pos_arr[perm]
        cosv=np.einsum("ij,ij->i", pos_s[rows_arr["i"]], pos_s[rows_arr["j"]])
        cosv=np.clip(cosv,-1,1)
        xi_s=np.arccos(cosv)
        rows_s=rows_arr.copy(); rows_s["xi"]=xi_s; rows_s["cos"]=cosv
        edges_s=equal_count_edges(xi_s, N_BINS)
        rb_s,_,_=bin_stats(rows_s, edges_s, design_fns)
        if np.isnan(rb_s).any(): continue
        rb_scrambles.append(rb_s)
    rb_scrambles=np.array(rb_scrambles)
    var_scramble_diag = rb_scrambles.var(axis=0)
    diagC = np.diag(C)
    ratio = diagC.mean()/var_scramble_diag.mean()
    print(label, "mean diag(C)=", diagC.mean(), "mean var(scramble rho_bin)=", var_scramble_diag.mean(), "ratio=", ratio)

run(None, "base(-14.62)")
import pathlib as _pl
run(_pl.Path("cache_alt_amp"), "alt(-14.0)")

# Result recorded manually in x3_diag_covariance_ratio_result.json:
# base(-14.62):  mean diag(C)=1.2419e-60  mean var(scramble rho_bin)=3.767e-61  ratio=3.297
# alt(-14.0):    mean diag(C)=8.9868e-60  mean var(scramble rho_bin)=2.909e-60  ratio=3.089
# Ratios are close (3.30 vs 3.09) despite the gate flipping pass->fail between the two runs,
# so the gate's amplitude-sensitivity is NOT explained by a simple diagonal-scale/normalization
# artifact of the frozen jackknife covariance C.
