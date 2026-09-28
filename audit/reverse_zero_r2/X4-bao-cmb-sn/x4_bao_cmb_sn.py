"""X4 - DESI DR2 BAO + CMB distance priors (Chen, Huang & Wang, arXiv:1808.05724) + Pantheon+ full cov.

Registered in audit/reverse_zero_r2/registration/registration.json (commit 9e6705e), test X4.
Subcommands (run from this directory, venv python, under prlimit):
  fit      --cut {1580,1590} --model {m0,lcdm,cpl} [--starts 20]   -> parts/fit_<cut>_<model>.json
  polish   --cut C --model M                                       -> parts/polish_<cut>_<model>.json (CAMB-exact CMB, sensitivity)
  regress  --cut C                                                 -> parts/regress_<cut>.json (no CMB, u profiled; reproduces E1/E3)
  assemble                                                         -> x4_results.json
Seeds: start i of every model uses numpy.random.default_rng(9000000 + i) (registered).
FRAMING: M0 is a hypothesis change (frozen Planck18 Omega_Lambda = 0.68885, GR tensor, no symmetron),
not a derivation from K3xT2. Result is INFORMATIONAL for PRE_REGISTRATION P1 (the cosmological
constant): DR2 BAO + CMB priors + Pantheon+ is not DR3/final BAO, Euclid or Rubin (rule at lines 40-44).
"""
import argparse, json, pathlib, subprocess, time, hashlib
import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.stats import chi2 as chi2dist, norm
import x4_core as X

HERE = pathlib.Path(__file__).resolve().parent
PARTS = HERE / "parts"; PARTS.mkdir(exist_ok=True)
SEED0 = 9000000
BOUNDS = {"Om": (0.15, 0.55), "h": (0.55, 0.85), "wb": (0.019, 0.027), "w0": (-2.5, 0.5), "wa": (-4.0, 2.0)}
NAMES = {"m0": ["h", "wb"], "lcdm": ["Om", "h", "wb"], "cpl": ["Om", "h", "wb", "w0", "wa"]}


def sigma_of(d, df):
    """Delta chi2 -> Gaussian-equivalent two-sided sigma with scipy (chi2.sf then norm.isf(p/2))."""
    if d <= 0:
        return 0.0, 1.0
    p = chi2dist.sf(d, df)
    return float(norm.isf(p / 2.0)), float(p)


def thresholds():
    return {f"{df}dof_{k}sigma": float(chi2dist.isf(2 * norm.sf(k), df)) for df in (1, 2) for k in (3, 5)}


def unpack(model, p):
    d = dict(zip(NAMES[model], p))
    if model == "m0":
        d["Om"] = X.om_m0(d["h"])
    d.setdefault("w0", -1.0); d.setdefault("wa", 0.0)
    return d


def make_obj(L, model, mode):
    def obj(p):
        d = unpack(model, p)
        if model == "cpl" and d["w0"] + d["wa"] >= 0:      # DESI-style prior w0 + wa < 0
            return 1e6 + 1e4 * (d["w0"] + d["wa"])
        wm = d["Om"] * d["h"] ** 2
        if mode == "paper" and not (0.0405 <= wm <= 0.419):    # inside the CAMB r_drag table
            return 1e6
        return L.chi2_full(d["Om"], d["h"], d["wb"], d["w0"], d["wa"], cmb_mode=mode)["chi2"]
    return obj


def nm(obj, x0, bounds, maxfev):
    r = minimize(obj, x0, method="Nelder-Mead", bounds=bounds,
                 options=dict(xatol=1e-9, fatol=1e-9, maxfev=maxfev, maxiter=maxfev, adaptive=True))
    return r


def cmd_fit(a):
    L = X.Likelihood(a.cut)
    names = NAMES[a.model]
    bnds = [BOUNDS[n] for n in names]
    lo = np.array([b[0] for b in bnds]); hi = np.array([b[1] for b in bnds])
    obj = make_obj(L, a.model, "paper")
    runs = []
    t0 = time.time()
    for i in range(a.start_lo, a.start_hi):
        rng = np.random.default_rng(SEED0 + i)
        u = rng.uniform(size=5)[:len(names)]
        x0 = lo + u * (hi - lo)
        if a.model == "cpl" and x0[3] + x0[4] >= 0:            # start inside the w0 + wa < 0 region
            x0[4] = -abs(x0[4]) - abs(x0[3]) - 0.01
            x0[4] = max(x0[4], BOUNDS["wa"][0])
        r = nm(obj, x0, bnds, a.maxfev)
        r = nm(obj, r.x, bnds, a.maxfev)                        # restart from the result
        runs.append({"start_index": i, "seed": SEED0 + i, "x0": x0.tolist(), "x": r.x.tolist(), "chi2": float(r.fun),
                     "nfev_last": int(r.nfev), "success": bool(r.success)})
        print(a.cut, a.model, i, round(r.fun, 4), round(time.time() - t0, 1), flush=True)
    best = min(runs, key=lambda r: r["chi2"])
    d = unpack(a.model, best["x"])
    full = L.chi2_full(d["Om"], d["h"], d["wb"], d["w0"], d["wa"], cmb_mode="paper")
    out = {"cut": a.cut, "n_sn": L.nsn, "model": a.model, "params_names": names, "starts": runs,
           "start_range": [a.start_lo, a.start_hi], "best": {**best, "params": d, "components": full,
           "at_bound": [n for n, v in zip(names, best["x"]) if abs(v - BOUNDS[n][0]) < 1e-6 or abs(v - BOUNDS[n][1]) < 1e-6]},
           "n_starts_within_1e-3_of_best": int(sum(r["chi2"] < best["chi2"] + 1e-3 for r in runs)),
           "bounds": {n: BOUNDS[n] for n in names}, "seconds": time.time() - t0}
    (PARTS / f"fit_{a.cut}_{a.model}_{a.start_lo}_{a.start_hi}.json").write_text(json.dumps(out, indent=1))


def load_best(cut, model):
    fs = sorted(PARTS.glob(f"fit_{cut}_{model}_*.json"))
    assert fs, (cut, model)
    parts = [json.loads(f.read_text()) for f in fs]
    runs = sum((p["starts"] for p in parts), [])
    idx = sorted(r["start_index"] for r in runs)
    best = min(runs, key=lambda r: r["chi2"])
    return parts[0], runs, best, idx


def cmd_polish(a):
    L = X.Likelihood(a.cut)
    p0, runs, best, _ = load_best(a.cut, a.model)
    names = NAMES[a.model]
    bnds = [BOUNDS[n] for n in names]
    obj = make_obj(L, a.model, "camb")
    r = nm(obj, np.array(best["x"]), bnds, a.maxfev)
    r = nm(obj, r.x, bnds, a.maxfev)
    d = unpack(a.model, r.x)
    d0 = unpack(a.model, best["x"])
    out = {"cut": a.cut, "model": a.model, "start_from_paper_best": d0, "camb_exact_best": d, "chi2_camb_best": float(r.fun),
           "chi2_camb_at_paper_best": float(obj(np.array(best["x"]))), "nfev": int(r.nfev),
           "components": L.chi2_full(d["Om"], d["h"], d["wb"], d["w0"], d["wa"], cmb_mode="camb"),
           "note": "R, l_A, r_drag straight from CAMB 2.0.4 (theta*-based), 1 massive nu 0.06 eV; sensitivity variant, not the registered primary"}
    (PARTS / f"polish_{a.cut}_{a.model}.json").write_text(json.dumps(out, indent=1))
    print(json.dumps({k: out[k] for k in ("chi2_camb_best", "chi2_camb_at_paper_best", "camb_exact_best")}, indent=1))


def cmd_regress(a):
    """No-CMB regression against E1/E3: Omega_r = 0, BAO amplitude and SN offset profiled analytically."""
    L = X.Likelihood(a.cut)
    Om_m0 = 1.0 - X.OMEGA_LAMBDA_M0
    m0 = L.chi2_nocmb(Om_m0)
    r = minimize_scalar(lambda o: L.chi2_nocmb(o)["chi2"], bounds=(0.02, 0.98), method="bounded", options={"xatol": 1e-10})
    lcdm = {"Om": float(r.x), **L.chi2_nocmb(r.x)}
    # E1-style: unbounded Nelder-Mead from (Om=0.31115, w0=-1, wa=0), exactly as e1_desi_dr2.py
    e1 = minimize(lambda p: L.chi2_nocmb(*p)["chi2"], [Om_m0, -1.0, 0.0], method="Nelder-Mead",
                  options={"xatol": 1e-8, "fatol": 1e-9, "maxiter": 10000, "maxfev": 10000})
    cpl_e1 = {"Om": float(e1.x[0]), "w0": float(e1.x[1]), "wa": float(e1.x[2]), **L.chi2_nocmb(*e1.x)}
    # bounded, 20 seeded starts (Om in [0.02,0.98], w0 in [-3,1], wa in [-5,3])
    bnds = [(0.02, 0.98), (-3.0, 1.0), (-5.0, 3.0)]
    lo = np.array([b[0] for b in bnds]); hi = np.array([b[1] for b in bnds])
    best = None
    for i in range(20):
        x0 = lo + np.random.default_rng(SEED0 + i).uniform(size=3) * (hi - lo)
        rr = minimize(lambda p: L.chi2_nocmb(*p)["chi2"], x0, method="Nelder-Mead", bounds=bnds,
                      options=dict(xatol=1e-8, fatol=1e-9, maxfev=4000, maxiter=4000))
        if best is None or rr.fun < best.fun:
            best = rr
    cpl = {"Om": float(best.x[0]), "w0": float(best.x[1]), "wa": float(best.x[2]), **L.chi2_nocmb(*best.x)}
    out = {"cut": a.cut, "n_sn": L.nsn, "M0": m0, "LCDM": lcdm, "CPL_Om_free_E1_style_unbounded": cpl_e1, "CPL_Om_free_bounded_20_seeded_starts": cpl,
           "dchi2_M0_minus_LCDM_1dof": m0["chi2"] - lcdm["chi2"],
           "dchi2_LCDM_minus_CPL_2dof_E1_style": lcdm["chi2"] - cpl_e1["chi2"],
           "dchi2_LCDM_minus_CPL_2dof_bounded_seeded": lcdm["chi2"] - cpl["chi2"],
           "note": "no radiation, no CMB, u and SN offset profiled, no w0+wa<0 constraint (E1 convention); an earlier unbounded seeded run drifted to Omega_m < 0 and was discarded, hence the bounds"}
    (PARTS / f"regress_{a.cut}.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


def gitq(*args, cwd):
    return subprocess.run(["git", *args], cwd=cwd, capture_output=True, text=True)


def grepF(path_or_blob, quote, git_repo=None, commit=None, file=None):
    """fixed-string assertion; records the command."""
    if git_repo:
        cmd = f"git -C {git_repo} show {commit}:{file} | grep -F -- '{quote}'"
        src = subprocess.run(["git", "-C", git_repo, "show", f"{commit}:{file}"], capture_output=True, text=True).stdout
    else:
        cmd = f"grep -F -- '{quote}' {path_or_blob}"
        src = pathlib.Path(path_or_blob).read_text()
    ok = quote in src
    return {"quote": quote, "command": cmd, "matched": bool(ok)}


def cmd_assemble(a):
    R = {"id": "X4", "script": "audit/reverse_zero_r2/X4-bao-cmb-sn/x4_bao_cmb_sn.py (+ x4_core.py, x4_camb_check.py, x4_fetch_prior.py)",
         "commands": [
             "cd audit/reverse_zero_r2/X4-bao-cmb-sn && <venv python> x4_fetch_prior.py",
             "cd audit/reverse_zero_r2/X4-bao-cmb-sn && prlimit --as=10737418240 -- <venv python> x4_camb_check.py",
             "cd audit/reverse_zero_r2/X4-bao-cmb-sn && prlimit --as=10737418240 -- <venv python> x4_bao_cmb_sn.py regress --cut {1580,1590}",
             "cd audit/reverse_zero_r2/X4-bao-cmb-sn && prlimit --as=10737418240 -- <venv python> x4_bao_cmb_sn.py fit --cut {1580,1590} --model {m0,lcdm,cpl} --start-lo I --start-hi J   (starts 0..19)",
             "cd audit/reverse_zero_r2/X4-bao-cmb-sn && prlimit --as=10737418240 -- <venv python> x4_bao_cmb_sn.py polish --cut {1580,1590} --model {m0,lcdm,cpl}",
             "cd audit/reverse_zero_r2/X4-bao-cmb-sn && prlimit --as=10737418240 -- <venv python> x4_bao_cmb_sn.py assemble"],
         "seeds": "minimiser starts default_rng(9000000 + start_index), start_index 0..19 (registered); r_d table validation default_rng(9100000); no other randomness",
         "framing": "M0 = frozen flat LCDM with the imported Planck18 Omega_Lambda = 0.68885, GR tensor sector, no symmetron sector: a hypothesis change, NOT a derivation from K3xT2. It carries profiled/free nuisances (h, omega_b free in the CMB-included fit; SN offset profiled). LeanMaster Stream 8 records 'Observables: none'. Nothing here says K3xT2 predicts or fixes mu_sym, c4, Omega_Lambda or any TDA outcome.",
         "tier": "X (data fit)"}
    # ---- provenance and registration checks
    repo = str(X.ROOT)
    reg_rel = "audit/reverse_zero_r2/registration/registration.json"
    show = gitq("show", f"9e6705e940f35aef8b377ec3631f62ba3d52aea4:{reg_rel}", cwd=repo)
    reg_sha_commit = hashlib.sha256(show.stdout.encode()).hexdigest()
    reg_sha_file = hashlib.sha256((X.ROOT / reg_rel).read_bytes()).hexdigest()
    anc = gitq("merge-base", "--is-ancestor", "9e6705e940f35aef8b377ec3631f62ba3d52aea4", "HEAD", cwd=repo)
    R["registration"] = {"commit": "9e6705e940f35aef8b377ec3631f62ba3d52aea4", "file": reg_rel,
                         "file_sha256_equals_commit_blob": reg_sha_commit == reg_sha_file, "commit_is_ancestor_of_HEAD": anc.returncode == 0,
                         "registered_test_X4_full": json.loads(show.stdout)["tests"]["X4"],
                         "deviations": [
                             {"what": "CPL minimiser rejects w0+wa >= 0 with a smooth penalty (a DESI-style prior), which the registered spec for CPL ('free: Omega_m, h, omega_b, w0, wa') does not impose.",
                              "reason": "carried over from the E1-style convention (round 2) to keep the minimiser away from a numerically pathological corner; not requested by X4's own registration.",
                              "effect": "inactive at every reported CPL best fit: w0+wa = -1.368 (1580) and -1.363 (1590), both well inside the unconstrained region, so the penalty changes none of the reported chi2/sigma numbers."},
                             {"what": "r_drag(omega_b, omega_m) is supplied by a CAMB-2.0.4-built bicubic-spline interpolation table (data/camb_rd_table.npz), not a hand-coded closed-form fitting formula (e.g. Eisenstein & Hu 1998).",
                              "reason": "the registered spec asks for r_d 'computed from (omega_b, omega_m) by one stated fitting formula'; CAMB is available in this environment and is more accurate than any closed-form fit.",
                              "effect": "validated to a maximum 1.4e-5 relative error against direct CAMB calls at 40 seeded points (seed 9100000, x4_camb_check.json.rd_spline_vs_camb); a strengthening, not a weakening, of the registered method."},
                             {"what": "the sample-size assert in x4_core.py (assert int(m.sum()) == cut) was added to the source file after this run's fits had already completed.",
                              "reason": "the registered spec requires the sample size be asserted before the fit; it is a no-op guard (both cuts already gave exactly 1580/1590), added for future reruns of this same script.",
                              "effect": "none on the numbers reported here, since assemble does not reload the SN data; disclosed for transparency about the order of operations in this run."}]}
    LM = "/home/callensxavier_gmail_com/SocrateAI-Scientific-Agora-LeanMaster"
    PR = "/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/audit/PRE_REGISTRATION.md"
    R["quotes_verbatim"] = [
        grepF(PR, "- 2 dof: Δχ² = 11.83 is 3σ, and 28.74 is 5σ."),
        grepF(PR, "- 1 dof: Δχ² = 9.00 is 3σ, and 25.00 is 5σ;"),
        grepF(PR, "- **Falsified if** any of these, analysed with full published covariances, prefers the w₀ > −1, wₐ < 0 quadrant over w = −1 at **≥ 5σ (Δχ² ≥ 28.74, 2 dof)**:"),
        grepF(PR, "  - DESI DR3 or final BAO combined with CMB and any one major SN compilation;"),
        grepF(PR, "| R4 | w = −1 exactly | DESI DR2 BAO + CMB + SNe (arXiv:2503.14738, 2503.14743), as recorded in LeanMaster `DualScaleValidation/Observables.lean` module docstring | **Disfavored at 3.1σ**"),
        grepF(None, "preferred over $\\Lambda$CDM at $3.1\\sigma$ for DESI BAO combined with CMB and supernovae", LM, "eb791e7", "DualScaleValidation/Observables.lean"),
        grepF(None, "Observables: none", LM, "eb791e7", "docs/STREAM8_WHICH_K3.md"),
    ]
    tex = X.HERE / "data" / "distance-priors-2018.tex"
    T = "audit/reverse_zero_r2/X4-bao-cmb-sn/data/distance-priors-2018.tex"
    R["paper"] = {"citation": "Chen, Huang & Wang, 'Distance Priors from Planck Final Release', arXiv:1808.05724, JCAP 02 (2019) 028 (doi 10.1088/1475-7516/2019/02/028), from arXiv metadata in data/fetch_log.json",
                  "fetch_log": json.loads((X.HERE / "data" / "fetch_log.json").read_text()),
                  "tex_source": T, "tex_sha256": hashlib.sha256(tex.read_bytes()).hexdigest(),
                  "verbatim_lines": {}, "quotes_verbatim": []}
    lines = tex.read_text().splitlines()
    want = {"table_LCDM_R_row": "$R$                                           & $1.7502\\pm0.0046$       & $1.0$&    $0.46$&   $-0.66$&   $-0.74$ \\\\",
            "table_LCDM_lA_row": "$l_\\textrm{A}$                                & $301.471^{+0.089}_{-0.090}$        & $0.46$&    $1.0$&   $-0.33$&   $-0.35$ \\\\",
            "table_LCDM_wb_row": "$\\Omega_b h^2$                              & $0.02236\\pm0.00015$      & $-0.66$&   $-0.33$&  $1.0$&    $0.46$\\\\",
            "distance_ini_r": "r=1.750235", "distance_ini_la": "la=301.4707", "distance_ini_wb": "omegabh2=0.02235976",
            "invcov_row1": "     94392.3971   -1360.4913   1664517.2916   ",
            "invcov_row2": "    -1360.4913   161.4349   3671.6180   ",
            "invcov_row3": "     1664517.2916   3671.6180   79719182.5162  ",
            "chi2_definition": "\\chi^2_{\\rm distance\\ priors}=\\sum(x_i-d_i)(C^{-1})_{ij}(x_j-d_j),",
            "z_star_use": "Here we use the approximate formula of $z_*$ to calculate $x_i$ \\cite{Hu:1995en}",
            "g1_prose": "g_1 &=& \\frac{0.0738(\\Omega_b h^2)^{-0.238}}{1+39.5(\\Omega_b h^2)^{0.763}}\\ ,\\\\",
            "g1_code": "g_1 = 0.0783D0*(CMB%ombh2)**(-0.238D0)/(1.0D0+39.5D0*(CMB%ombh2)**0.763D0)",
            "Omega_r_prescription": "\\Omega_{r}=\\frac{\\Omega_{m}}{1+z_{\\textrm{eq}}}, ~~z_{\\textrm{eq}}=2.5\\times10^{4}\\Omega_{m}h^2\\left(T_{CMB}/2.7\\textrm{K}\\right)^{-4}\\ .",
            "likelihood_lA_in_code": "l_a=const_pi/this%Calculator%CMBToTheta(CMB)"}
    for k, s in want.items():
        hits = [i + 1 for i, l in enumerate(lines) if s in l]
        R["paper"]["verbatim_lines"][k] = {"tex_line": hits[0] if hits else None, "text": lines[hits[0] - 1] if hits else None, "matched": bool(hits)}
        R["paper"]["quotes_verbatim"].append({"quote": s, "command": f"grep -nF -- '{s}' {T}", "matched": bool(hits)})
    R["paper"]["table_caption_dataset"] = "Planck TT,TE,EE+lowE (base LCDM); the published table also lists wCDM, LCDM+Omega_k and LCDM+A_L, not used"
    R["paper"]["values_used"] = {"means": X.CMB_D.tolist(), "inverse_covariance_before_normalisation": X.CMB_INVC.tolist(),
                                 "z_star": "Hu-Sugiyama fit with the code coefficient 0.0783 (the prose prints 0.0738; a 0.0738 variant is not the primary)"}
    R["camb_check"] = json.loads((HERE / "x4_camb_check.json").read_text())
    R["thresholds_recomputed_scipy"] = thresholds()
    # ---- regression
    E1 = json.loads((X.ROOT / "audit/reverse_zero/E1-desi-dr2/e1_desi_dr2_report.json").read_text())
    E3 = json.loads((X.ROOT / "audit/reverse_zero/E3-M0/m0_result.json").read_text())
    R["regression_no_cmb"] = {}
    for cut in (1580, 1590):
        f = PARTS / f"regress_{cut}.json"
        if f.exists():
            R["regression_no_cmb"][str(cut)] = json.loads(f.read_text())
    R["regression_no_cmb"]["E1_json_DR2_dchi2_CPL_vs_Lambda_2dof_1590"] = E1["DR2"]["delta_chi2_CPL_vs_Lambda_2dof"]
    R["regression_no_cmb"]["E3_json_dchi2_M0_minus_LCDM_1dof_1580"] = E3["delta_chi2_M0_minus_LCDM_fitted"]
    R["regression_no_cmb"]["ERRATA_quoted_CPL_dchi2_1580_and_1590"] = "ERRATA.md: 'DR2 Delta chi2 is 4.547 and 4.735' (1590 and 1580); not re-derived from a JSON"
    # ---- fits
    R["fits"] = {}
    for cut in (1580, 1590):
        F = {"n_sn": None}
        best = {}
        Lc = X.Likelihood(cut)
        for m in ("m0", "lcdm", "cpl"):
            try:
                p0, runs, b, idx = load_best(cut, m)
            except AssertionError:
                continue
            d = unpack(m, b["x"])
            # components at the TRUE global best (min over all chunk files), not a chunk-local "best"
            comps = Lc.chi2_full(d["Om"], d["h"], d["wb"], d["w0"], d["wa"], cmb_mode="paper")
            assert abs(comps["chi2"] - b["chi2"]) < 1e-4, (comps["chi2"], b["chi2"])
            n_conv = int(sum(r["chi2"] < b["chi2"] + 1e-3 for r in runs))
            best[m] = {"chi2": b["chi2"], "params": d, "n_starts": len(runs), "start_indices": idx,
                       "n_starts_within_1e-3_of_best": n_conv,
                       "non_converged_start_indices": [r["start_index"] for r in runs if r["chi2"] >= b["chi2"] + 1e-3],
                       "chi2_spread_all_starts": [min(r["chi2"] for r in runs), float(np.median([r["chi2"] for r in runs])), max(r["chi2"] for r in runs)],
                       "at_bound": [n for n, v in zip(NAMES[m], b["x"]) if abs(v - BOUNDS[n][0]) < 1e-6 or abs(v - BOUNDS[n][1]) < 1e-6],
                       "components_at_best": {k: comps[k] for k in ("bao", "sn", "cmb")}}
            F["n_sn"] = p0["n_sn"]
            pf = PARTS / f"polish_{cut}_{m}.json"
            if pf.exists():
                best[m]["camb_exact_polish"] = json.loads(pf.read_text())
        F["models"] = best
        if all(m in best for m in ("m0", "lcdm", "cpl")):
            c0, cl, cc = best["m0"]["chi2"], best["lcdm"]["chi2"], best["cpl"]["chi2"]
            def comp(d, df, label):
                s, p = sigma_of(d, df)
                return {"delta_chi2": d, "dof": df, "p_wilks": p, "sigma": s, "label": label}
            cp = best["cpl"]["params"]
            F["CPL_vs_Lambda_2dof_primary"] = {**comp(cl - cc, 2, "CPL vs flat Lambda, Omega_m free in both (registered primary)"),
                                               "w0": cp["w0"], "wa": cp["wa"], "in_w0_gt_-1_wa_lt_0_quadrant": bool(cp["w0"] > -1 and cp["wa"] < 0),
                                               "ge_3sigma_threshold_2dof": bool(cl - cc >= thresholds()["2dof_3sigma"]),
                                               "ge_5sigma_threshold_2dof": bool(cl - cc >= thresholds()["2dof_5sigma"])}
            F["M0_vs_fitted_Lambda_1dof"] = {**comp(c0 - cl, 1, "M0 minus fitted flat Lambda (reported alongside, uncorrected; not a rejection test)"),
                                             "not_rejected_if_p_ge_0.05": bool(chi2dist.sf(max(c0 - cl, 0), 1) >= 0.05),
                                             "M0_h": best["m0"]["params"]["h"], "Lambda_Om": best["lcdm"]["params"]["Om"], "Lambda_h": best["lcdm"]["params"]["h"]}
            F["M0_vs_CPL_3dof_informational"] = comp(c0 - cc, 3, "M0 minus CPL (3 extra parameters)")
            F["published_3.1sigma_comparison"] = {"published_sigma_recorded_in_PRE_REG_R4": 3.1,
                                                  "this_fit_sigma": F["CPL_vs_Lambda_2dof_primary"]["sigma"],
                                                  "shortfall_sigma": 3.1 - F["CPL_vs_Lambda_2dof_primary"]["sigma"],
                                                  "note": "3.1 is the number recorded in PRE_REGISTRATION R4 / LeanMaster Observables.lean for DESI BAO + CMB + SNe (their own analysis, different CMB likelihood and SN sample handling). This fit uses compressed CMB priors and one SN sample; the comparison is informational."}
        F["sensitivity_camb_exact"] = {}
        for m in ("m0", "lcdm", "cpl"):
            if m in best and "camb_exact_polish" in best[m]:
                F["sensitivity_camb_exact"][m] = best[m]["camb_exact_polish"]["chi2_camb_best"]
        if len(F["sensitivity_camb_exact"]) == 3:
            s = F["sensitivity_camb_exact"]
            F["sensitivity_camb_exact"]["dchi2_CPL_vs_Lambda_2dof"] = s["lcdm"] - s["cpl"]
            F["sensitivity_camb_exact"]["sigma_CPL_vs_Lambda"] = sigma_of(s["lcdm"] - s["cpl"], 2)[0]
            F["sensitivity_camb_exact"]["dchi2_M0_minus_Lambda_1dof"] = s["m0"] - s["lcdm"]
            F["sensitivity_camb_exact"]["sigma_M0_vs_Lambda"] = sigma_of(s["m0"] - s["lcdm"], 1)[0]
        R["fits"][str(cut)] = F
    R["N_data"] = {"BAO": 13, "CMB_priors": 3, "SN_1580": 1580, "SN_1590": 1590}
    R["p1_status"] = ("INFORMATIONAL. PRE_REGISTRATION P1 (the cosmological constant) names DESI DR3 or final BAO + CMB + one SN compilation, Euclid, or Rubin/LSST. "
                      "DESI DR2 BAO + compressed Planck-2018 distance priors + Pantheon+ is none of these, so no P1 verdict follows. The words 'falsified' and 'confirmed' are not used.")
    (HERE / "x4_results.json").write_text(json.dumps(R, indent=1, default=float))
    print(json.dumps({k: R["fits"][k].get(f) for k in R["fits"] for f in ("CPL_vs_Lambda_2dof_primary", "M0_vs_fitted_Lambda_1dof")}, indent=1, default=float))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("fit", "polish", "regress"):
        s = sub.add_parser(name)
        s.add_argument("--cut", type=int, required=True, choices=[1580, 1590])
        if name != "regress":
            s.add_argument("--model", required=True, choices=list(NAMES))
        if name == "fit":
            s.add_argument("--start-lo", type=int, default=0); s.add_argument("--start-hi", type=int, default=20)
            s.add_argument("--maxfev", type=int, default=2500)
        if name == "polish":
            s.add_argument("--maxfev", type=int, default=600)
    sub.add_parser("assemble")
    a = ap.parse_args()
    {"fit": cmd_fit, "polish": cmd_polish, "regress": cmd_regress, "assemble": cmd_assemble}[a.cmd](a)
