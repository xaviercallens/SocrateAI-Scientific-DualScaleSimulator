# Round-2 skeptic ledger (independent re-runs; every number below came from a command listed here)

Python: /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python (gudhi 3.13.0, numpy 2.2.6). Seed 42.

## 1. Decisive experiment re-run (copy of decisive_experiment.py, REPO_ROOT patched, run in scratchpad)
chi2_total: LCDM Om fitted 702.6753004830802 | LCDM Om frozen 0.31115 702.9066083231721 | Om=0.5 856.1487994018894 |
CPL(Om frozen) 696.6192624361269 | CPL(all free) 696.5836572924258. Diff vs committed report: 0.0 on all five.
DESI rows hard-coded in the script were compared line-by-line with data/real/dark_energy/desi_2024_bao_all.txt (sha256 dd2873a0...): identical.
Tighter negative controls (same chi2 code): Om=0.28 -> 712.9134 (+10.24); Om=0.30 -> 704.7109; Om=0.33 -> 703.7895; Om=0.35 -> 709.4778 (+6.80); Om=0.40 -> 741.0223.
1-sigma (dchi2=1) Om interval from this pipeline: [0.30499, 0.32930]; 3-sigma (dchi2=9): [0.28220, 0.35520]. Frozen 0.31115 is inside the 1-sigma interval.
Caveat: Pantheon+ uses diagonal errors only (chi2_SN/N = 688/1590 = 0.43), so absolute chi2 thresholds are not calibrated; the test is weaker, not stronger, than nominal.

## 2. Harness x10 / x0.1 probes at nominal (skeptic/probe.py -> probe_report.json), chi2 = same BAO+SN machinery
nominal:        w0=-1.00103 wa=0.00624 H(z=1)/H0=1.00056 chi2_total=6429.02 (matter fully diluted at a0=a(t_max): effectively de Sitter)
a_pot x10:      w0=-1.00022 chi2=6452.09 ; a_pot x0.1: w0=-1.0000 chi2=6458.80
b_pot x10:      w0=-1.0000  chi2=6458.80 ; b_pot x0.1: w0=-0.94090 wa=1.07746 H(z=1)/H0=1.42902 chi2=1072.96
=> a_pot, b_pot ARE observable (chi2 moves 6429 -> 1073 -> 703.89 at the sweep best point). Mechanism label "delete_unobservable" is WRONG for them.
lambda_sym x10: ssf 1.1811735180570304e-04 ; x0.1: 1.1811734233485311e-04 ; nominal 1.1811720750749558e-04 (rel 1.2e-6, non-monotonic) ; DE and PTA blocks byte-identical.
mu_sym x0.1 (positive control): ssf 3.2785215660758164e-05 (factor 3.60) ; mu_sym x10: NaN (numerically_stable False).
Exact algebra checked by reading workshopcosmo.py:316-372: phi = phi_0*psi with phi_0^2 = mu^2/lambda turns lambda*phi^3 into phi_0*mu^2*psi^3; BCs (psi'(eps)=0, psi(r_max)=1) and the initial guess also scale with phi_0. lambda_sym cancels exactly.

## 3. Sweep (audit/zero_param_loop/round2/sweep_chi2.csv)
513 rows, 378 finite; min chi2_total 703.8946371287691; 0 rows below frozen-LCDM 702.9066; 0 rows below 703.5. Top-5 rows all have b_pot ~ 4e-4..6e-4, w0 ~ -0.70, wa ~ 0.95..1.02.

## 4. TDA re-run (copy of tda_gudhi.py in scratchpad, same sweep_chi2.csv)
import gudhi confirmed (gudhi.RipsComplex / bottleneck_distance used). Betti model [2,0], real [1,0]: identical to committed. All six bottleneck numbers: diff 0.0.
Circle known-answer: 25 H1 bars, top 0.69086, second 0.01962, ratio 35.216 -> PASS.
DEGENERACY (new finding): every X_vs_real H1 bottleneck equals exactly half of X's own largest H1 bar:
 model 0.011094914869650904 = 0.02218982973930181/2 ; null_box 0.03871767485515332 = 0.07743534971030665/2 ;
 null_vm 0.027568786206406695 = 0.05513757241281339/2 ; null_poisson 0.03563970809302948 = 0.07127941618605896/2.
 The real diagram (max H1 persistence 0.013562) never enters any of them. "model_vs_real smaller than all three nulls" is a statement about the model cloud's largest bar, not about the real data.
Null Poisson (real side) Betti [1,0] == real Betti [1,0]: the null control is NOT distinguishable from the real data by any statistic computed here.

## 5. LeanMaster (READ-ONLY, HEAD 3e60cd0, tag v3.16.1)
git show v3.16.1:DualScaleCosmology/DarkEnergyScale.lean | sed -n 84p -> `def omegaLambda : ℝ := 0.68885`
git show v3.16.1:DualScaleCosmology/SelfDualCutoff.lean | sed -n 100p -> `def hubbleRadius_m : ℝ := 1.3672e26`
CKNInstance.lean:63 `theorem planckEnergy_gt_ckn_horizon_cutoff : planckEnergy_eV ≥ (10 : ℝ) ^ 30 * cknLambdaHorizon_eV := by unfold ...; norm_num` (present)
CharactersAll.lean:215 `theorem ratio_fails_at_every_class : ∀ j : Fin 26, j.val ≠ 0 → ¬ ratioTwines ... := by decide` (present)
DualScaleDyons/Immortal.lean:86 immortal_m1, ImmortalHigher.lean:88/95 immortal_m2/m3 (present).
"587 theorems": NOT independently re-derived with LeanMaster's own gate; a line-start `git grep -c "^theorem\|^lemma"` over tracked .lean gives 667 declarations (different counting method). Reported as ABSENT-verification, not as a contradiction.
H0 = c/hubbleRadius_m = 299792458/1.3672e26 * 3.0856775814913673e22/1000 = 67.661122494938 km/s/Mpc (recomputed).
