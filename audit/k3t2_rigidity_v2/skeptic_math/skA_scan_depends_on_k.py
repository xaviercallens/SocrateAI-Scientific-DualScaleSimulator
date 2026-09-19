"""Skeptic check of Track A 'rigidity scan (a)': does the slice-agreement selector pin N=24
independently of the overall normalisation k of Z = k*phi01?  Reuses Track A's own
build_blocks/extract_H (imported read-only) but rescales Z*eta^3 by k/2.
Command: cd audit/k3t2_rigidity_v2/skeptic_math && /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python skA_scan_depends_on_k.py
"""
import sys, os, json
from fractions import Fraction as Fr
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", "A-genus"))
os.chdir(os.path.join(here, "..", "A-genus"))
import run_track_a_v2 as A
from series2d import scal
CUT = 6
imax, ZK3, phi01, T1sq, e3, ZK3_eta3, y12Psi = A.build_blocks(CUT, CUT + 10)
res = {}
for k in (Fr(1), Fr(2), Fr(3)):
    Zk_eta3 = scal(ZK3_eta3, k / 2)
    sel = [N for N in range(6, 40) if A.extract_H(N, Zk_eta3, y12Psi, T1sq, imax, 16)[1]]
    res[str(k)] = sel
out = {"cutoff_q": CUT, "selected_N_by_k": res,
       "reading": "slice agreement selects N = 12k = Z(tau,0): it recovers the normalisation of Z, it does not force 24 independently of k"}
json.dump(out, open(os.path.join(here, "skA_scan_depends_on_k_results.json"), "w"), indent=1)
print(out)
