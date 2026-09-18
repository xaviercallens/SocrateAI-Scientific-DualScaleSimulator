#!/usr/bin/env python3
"""
REVERSE LOOP step 1: build the REDUCED model.

Applies ONLY the round-1 tier-B accepted reduction (see
../reduce_proposals.json accepted_this_round_tier_B): the pair
(pta_suppression, c4_c0_ratio) enters every observable this harness
computes ONLY through the product c4_c0_ratio*pta_suppression
(workshopcosmo.py:766; round-0 Check C negative control, max diff
1.110e-16; round-1 jacobian_report.json: the DIFFERENCE direction
ln(pta_suppression)-ln(c4_c0_ratio) is an EXACT machine-precision null
direction -- singular value 0.0, loadings exactly -0.7071/+0.7071 -- at
ALL 6 probed points). That pair is therefore collapsed to ONE parameter,
c4_pta_product.

lambda_sym is NOT removed here. Round-1 found it only "0.5-evidence /
fixed by convention" (insensitivity, not a derived replacement value) --
per the loop's ground rule (c), insensitivity alone is not derivation, so
it STAYS in the free-parameter count.

REDUCED FREE PARAMETERS (5): a_pot, b_pot, mu_sym, lambda_sym, c4_pta_product.

The product is split back into (pta_suppression, c4_c0_ratio) = (sqrt(P),
sqrt(P)) before calling the untouched param_loop_sim.evaluate_point. This
split is ARBITRARY -- any positive (x, P/x) pair gives byte-identical
results, because no observable this harness computes depends on x itself,
only on the product P (that is exactly what "exact null direction, SV=0.0"
means). Verified below in a small self-test rather than asserted.
"""
import os
import sys
from typing import Any, Dict

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
sys.path.insert(0, REPO_ROOT)
sys.path.insert(0, os.path.join(REPO_ROOT, "scripts"))
import param_loop_sim as pls  # noqa: E402

REDUCED_PARAM_NAMES = ["a_pot", "b_pot", "mu_sym", "lambda_sym", "c4_pta_product"]

# Nominal reduced defaults: identical physical point to pls.DEFAULT_PARAMS
# (c4_pta_product = 16.07 * 0.005 = 0.08035 exactly).
REDUCED_DEFAULT_PARAMS: Dict[str, float] = {
    "a_pot": pls.DEFAULT_PARAMS["a_pot"],
    "b_pot": pls.DEFAULT_PARAMS["b_pot"],
    "mu_sym": pls.DEFAULT_PARAMS["mu_sym"],
    "lambda_sym": pls.DEFAULT_PARAMS["lambda_sym"],
    "c4_pta_product": pls.DEFAULT_PARAMS["c4_c0_ratio"] * pls.DEFAULT_PARAMS["pta_suppression"],
}


def evaluate_reduced_point(params: Dict[str, float]) -> Dict[str, Any]:
    merged = dict(REDUCED_DEFAULT_PARAMS)
    for k in params:
        if k not in REDUCED_DEFAULT_PARAMS:
            raise SystemExit(f"Unknown reduced parameter '{k}'. Valid keys: {sorted(REDUCED_DEFAULT_PARAMS)}")
    merged.update(params)
    product = float(merged["c4_pta_product"])
    if product <= 0:
        split = float("nan")
    else:
        split = product ** 0.5  # arbitrary symmetric split; see module docstring
    full_params = {
        "a_pot": merged["a_pot"],
        "b_pot": merged["b_pot"],
        "mu_sym": merged["mu_sym"],
        "lambda_sym": merged["lambda_sym"],
        "pta_suppression": split,
        "c4_c0_ratio": split,
    }
    result = pls.evaluate_point(full_params)
    result["reduced_params"] = merged
    result["split_pta_suppression"] = split
    result["split_c4_c0_ratio"] = split
    return result


def _selftest() -> bool:
    """Verify the split is inconsequential: two different splits of the SAME
    product must give byte-identical pta.gamma_theta, and the reduced model
    at its nominal point must reproduce the full model's nominal PTA output
    exactly (not just approximately)."""
    base_full = pls.evaluate_point({})  # full-model nominal
    base_reduced = evaluate_reduced_point({})  # reduced-model nominal
    import numpy as np
    diff_gamma = float(np.max(np.abs(
        np.array(base_full["pta"]["gamma_theta"]) - np.array(base_reduced["pta"]["gamma_theta"])
    )))
    diff_prod = abs(base_full["pta"]["c4_pta_product"] - base_reduced["pta"]["c4_pta_product"])

    # Now try a DIFFERENT arbitrary split of the same product (0.001, product/0.001)
    # by calling pls.compute_pta_observable directly, and confirm identical gamma_theta.
    product = REDUCED_DEFAULT_PARAMS["c4_pta_product"]
    alt = pls.compute_pta_observable(0.001, product / 0.001)
    diff_alt = float(np.max(np.abs(np.array(alt["gamma_theta"]) - np.array(base_reduced["pta"]["gamma_theta"]))))

    ok = (diff_gamma < 1e-12) and (diff_prod < 1e-12) and (diff_alt < 1e-12)
    print(f"[reduced_model selftest] diff_gamma(full vs reduced, nominal)={diff_gamma:.3e}")
    print(f"[reduced_model selftest] diff_c4_pta_product={diff_prod:.3e}")
    print(f"[reduced_model selftest] diff_gamma(reduced-split vs alt-split 0.001/{product/0.001:.4f})={diff_alt:.3e}")
    print(f"[reduced_model selftest] ALL_PASS={ok}")
    return ok


if __name__ == "__main__":
    ok = _selftest()
    raise SystemExit(0 if ok else 1)
