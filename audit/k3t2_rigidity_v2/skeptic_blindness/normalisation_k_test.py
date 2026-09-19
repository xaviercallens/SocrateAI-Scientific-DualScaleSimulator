"""Skeptic (blindness/circularity): is Track A's 'RIGID N=24' scan selecting 24
structurally, or selecting 12*k where k=2 is the typed overall factor in
ZK3 = scal(phi01, 2) (A-genus/run_track_a_v2.py:81,209)?

Rebuilds Track A's blocks via its own functions, replaces ZK3 = k*phi01 for
k in {1,2,3}, and reruns Track A's own extract_H slice-agreement selector over
N in [12k-4, 12k+4]. No target value is used: we only report which N pass.
Run: cd audit/k3t2_rigidity_v2/skeptic_blindness &&
  /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python normalisation_k_test.py
"""
import json, os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "A-genus"))
import run_track_a_v2 as A
from series2d import mul, scal, y_to_1

H_CUTOFF, MARGIN = 14, 16
imax, ZK3_2, phi01, T1sq, e3, _, y12Psi = A.build_blocks(H_CUTOFF, H_CUTOFF + 10)
phi0_const = y_to_1(phi01).get(0)
out = {"phi01_tau_z0_constant_computed": str(phi0_const), "rows": {}}
for k in (1, 2, 3):
    ZK3_k = scal(phi01, k)
    ZK3_eta3_k = mul(ZK3_k, e3, imax)
    centre = int(k * phi0_const)   # 12k from computed phi01(tau,0); a scan window, not a selector
    passing = []
    for N in range(centre - 4, centre + 5):
        _, agree, _, _, _ = A.extract_H(N, ZK3_eta3_k, y12Psi, T1sq, imax, MARGIN)
        if agree:
            passing.append(N)
    out["rows"][str(k)] = {"typed_factor_k": k, "N_window": [centre - 4, centre + 4],
                           "N_passing_slice_agreement": passing,
                           "Z_k(tau,0)": str(k * phi0_const)}
out["conclusion"] = ("slice agreement selects N = Z(tau,0) = k*phi01(tau,0) for every typed k; "
                     "it pins N to the normalisation, so N=24 is conditional on the typed k=2")
out["conclusion_holds"] = all(r["N_passing_slice_agreement"] == [int(Fr(r["Z_k(tau,0)"]))] for r in out["rows"].values())
json.dump(out, open(os.path.join(HERE, "normalisation_k_test_results.json"), "w"), indent=2)
print(json.dumps(out, indent=2))
