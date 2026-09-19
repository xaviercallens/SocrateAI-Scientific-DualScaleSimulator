"""Round-3 skeptic: independent check of commit 8b8911c (lambda_sym deletion).
(1) psi-form (current workshopcosmo) vs verbatim legacy phi-form at lambda/mu values
    NOT used by test_lambda_deletion.py.
(2) Mutations of the psi-form, injected by source-patching the CURRENT function, must be
    detected by the same comparison (power of the test).
Writes lambda_mutation_probe.json next to this file."""
import inspect, json, os, sys, textwrap
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "audit", "zero_param_loop", "lambda_deletion"))
import workshopcosmo as wc
from legacy_phi_form import legacy_phi_form
KEYS = ["screening_suppression_factor", "phi_center_ratio", "phi_surface_ratio"]
def rel(a, b): return abs(a - b) / max(abs(a), abs(b), 1e-300)
src = textwrap.dedent(inspect.getsource(wc.run_symmetron_screening_simulation))
def patched(old, new):
    assert src.count(old) >= 1, old
    ns = dict(wc.__dict__); exec(src.replace(old, new), ns); return ns["run_symmetron_screening_simulation"]
MUTANTS = {
    "cubic_x1.5": patched("+ (mu_sym ** 2) * (psi ** 3)\n", "+ 1.5 * (mu_sym ** 2) * (psi ** 3)\n"),
    "bc_psi_eq_mu": patched("yb[0] - 1.0", "yb[0] - mu_sym"),
    "cubic_no_mu2": patched("+ (mu_sym ** 2) * (psi ** 3)\n", "+ (psi ** 3)\n"),
}
out = {"equivalence": [], "mutants": {}}
for mu in (0.5, 1.5):
    new = wc.run_symmetron_screening_simulation(mu_sym=mu)
    for lam in (0.003, 0.37, 7.3, 250.0):
        old = legacy_phi_form(mu_sym=mu, lambda_sym=lam)
        out["equivalence"].append({"mu": mu, "lam": lam, **{k: [new[k], old[k], rel(new[k], old[k])] for k in KEYS}})
    for name, f in MUTANTS.items():
        m = f(mu_sym=mu); old = legacy_phi_form(mu_sym=mu, lambda_sym=0.37)
        out["mutants"].setdefault(name, []).append({"mu": mu, "ssf_rel_diff_vs_legacy": rel(m["screening_suppression_factor"], old["screening_suppression_factor"]),
                                                   "pcr_rel_diff_vs_legacy": rel(m["phi_center_ratio"], old["phi_center_ratio"])})
out["max_equiv_rel_diff"] = max(max(e[k][2] for k in KEYS) for e in out["equivalence"])
out["min_mutant_ssf_rel_diff"] = {n: min(x["ssf_rel_diff_vs_legacy"] for x in v) for n, v in out["mutants"].items()}
json.dump(out, open(os.path.join(HERE, "lambda_mutation_probe.json"), "w"), indent=1)
print("max_equiv_rel_diff", out["max_equiv_rel_diff"]); print("min_mutant_ssf_rel_diff", out["min_mutant_ssf_rel_diff"])
