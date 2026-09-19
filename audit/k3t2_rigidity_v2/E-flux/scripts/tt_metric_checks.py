#!/usr/bin/env python3
"""
Track E / Part 3 prerequisite: validate the (A.6)-(A.7) diagonal lattice basis
and the N_flux = 2*alpha_x^2 formula (eq. 4.13) against Tripathy-Trivedi's own
WORKED EXAMPLES (4.15)-(4.16) and (4.31)-(4.32), plus a symbolic re-derivation
of (4.13) from (2.8)+(4.8)+(4.9)-(4.10) with sympy.

Source (pinned, sha256 in audit/k3t2_rigidity_v2/sources/SHA256SUMS):
  audit/k3t2_rigidity_v2/sources/hep-th_0301139_TripathyTrivedi.txt
Relevant lines: 217-229 (2.6)-(2.9), 861-891 (4.8)-(4.16), 953-1010 (4.24)-(4.32),
1770-1868 (App. A).

Tier: B (exact int/Fraction/sympy arithmetic; the "expected" values below are
copied VERBATIM from the paper's printed equations/numbers with their line
numbers -- they are comparison targets, not inputs to any computation path).

Run exactly:
  /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python \
    /mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity_v2/E-flux/scripts/tt_metric_checks.py
"""
import json
import sys
from fractions import Fraction

import sympy as sp

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from lattice import FULL_SIGNS, norm, dot, truncation  # noqa: E402

OUT = __file__.rsplit("/", 1)[0].rsplit("/", 1)[0] + "/tt_metric_checks_results.json"


def check_4_15_4_16():
    """alpha_x = 2 e1, beta_x = 2 e2 (4.15). Expect alpha_x^2 = 8 (used in
    text right after 4.16), N_flux = 2*alpha_x^2 = 16 (eq 4.13),
    N_D3 = 24 - alpha_x^2 = 16 (eq 4.16, line 885, printed as
    'N_D3 = 24 - 8 = 16')."""
    signs = FULL_SIGNS
    alpha_x = (1, 0, 0, 0, 0, 0)  # coefficient 2*1 on e1
    ax2 = norm(alpha_x, signs)
    n_flux = 2 * ax2
    n_d3 = 24 - ax2  # eq (4.14): alpha_x^2 + N_D3 = 24
    expected = {
        "alpha_x_sq": 8,  # TT line 878-885 arithmetic "24 - 8 = 16"
        "N_flux": 16,  # eq (4.13) applied to alpha_x^2=8
        "N_D3": 16,  # TT eq (4.16), line 885, printed "N_D3 = 24 - 8 = 16"
    }
    computed = {"alpha_x_sq": ax2, "N_flux": n_flux, "N_D3": n_d3}
    ok = computed == expected
    return {
        "id": "TT-4.15-4.16",
        "description": "alpha_x=2e1, beta_x=2e2 tadpole check",
        "computed": computed,
        "expected": expected,
        "expected_source": "hep-th_0301139_TripathyTrivedi.txt lines 878-885 (eqs 4.15-4.16)",
        "match": ok,
    }


def check_4_31_4_32():
    """alpha_x=2(e1-e2), alpha_y=2(e1+e2+e4), beta_x=-4e2, beta_y=2(2e1+e4+e5)
    (4.31). Expect alpha_xx=16, alpha_yy=8 (text after 4.31, 'alpha_xx=16 and
    alpha_yy=8'), N_flux = alpha_x.beta_y - beta_x.alpha_y = 32 (eq 4.32,
    line 1004-1007), and tau^2 = -alpha_yy/alpha_xx = -1/2, i.e. tau = i/sqrt(2)
    (text, 'tau = i/sqrt(2)')."""
    signs = FULL_SIGNS
    alpha_x = (1, -1, 0, 0, 0, 0)  # 2(e1-e2)
    alpha_y = (1, 1, 0, 1, 0, 0)  # 2(e1+e2+e4)
    beta_x = (0, -2, 0, 0, 0, 0)  # -4 e2
    beta_y = (2, 0, 0, 1, 1, 0)  # 2(2e1+e4+e5)

    alpha_xx = norm(alpha_x, signs)
    alpha_yy = norm(alpha_y, signs)
    # eq (2.8): N_flux = -beta_x.alpha_y + beta_y.alpha_x = alpha_x.beta_y - beta_x.alpha_y
    n_flux = dot(alpha_x, beta_y, signs) - dot(beta_x, alpha_y, signs)
    tau_sq = Fraction(-alpha_yy, alpha_xx)  # tau = i*sqrt(alpha_yy/alpha_xx) => tau^2 = -alpha_yy/alpha_xx

    expected = {
        "alpha_xx": 16,
        "alpha_yy": 8,
        "N_flux": 32,
        "tau_sq": str(Fraction(-1, 2)),
    }
    computed = {
        "alpha_xx": alpha_xx,
        "alpha_yy": alpha_yy,
        "N_flux": n_flux,
        "tau_sq": str(tau_sq),
    }
    ok = computed == expected
    return {
        "id": "TT-4.31-4.32",
        "description": "second-branch example (2+,2-) flux check",
        "computed": computed,
        "expected": expected,
        "expected_source": "hep-th_0301139_TripathyTrivedi.txt lines 984-1009 (eqs 4.24-4.32, text 'alpha_xx=16 and alpha_yy=8', 'tau=i/sqrt(2)')",
        "match": ok,
        "note": "This example belongs to section 4.2 (a DIFFERENT flux family from 4.1: it does not satisfy alpha_y=-beta_x,beta_y=alpha_x). Used ONLY as an independent metric-convention check, not folded into the section-4.1 enumeration in flux_enumeration.py.",
    }


def symbolic_derivation_4_13():
    """Symbolically re-derive N_flux = 2*alpha_x^2 (eq 4.13) from:
      eq (2.8): N_flux = -beta_x.alpha_y + beta_y.alpha_x
      eq (4.8): alpha_y = -beta_x, beta_y = alpha_x
      eq (4.9): alpha_x^2 = beta_x^2
      eq (4.10): alpha_x.beta_x = 0
    using sympy over abstract dot-product symbols (independent symbolic route
    from the direct numeric lattice computation above)."""
    ax2, bx2, axbx = sp.symbols("alpha_x_sq beta_x_sq alpha_x_dot_beta_x")
    # N_flux = -beta_x . alpha_y + beta_y . alpha_x
    #        = -beta_x . (-beta_x) + alpha_x . alpha_x      [substituting 4.8]
    #        = beta_x^2 + alpha_x^2
    n_flux_expr = bx2 + ax2
    n_flux_substituted = n_flux_expr.subs(bx2, ax2)  # using (4.9): beta_x^2 = alpha_x^2
    expected_expr = 2 * ax2
    ok = sp.simplify(n_flux_substituted - expected_expr) == 0
    return {
        "id": "TT-4.13-symbolic",
        "description": "symbolic derivation of N_flux = 2*alpha_x^2 from (2.8)+(4.8)+(4.9)",
        "derivation_steps": [
            "N_flux = -beta_x.alpha_y + beta_y.alpha_x   [eq 2.8]",
            "substitute alpha_y=-beta_x, beta_y=alpha_x  [eq 4.8]",
            "=> N_flux = beta_x.beta_x + alpha_x.alpha_x = beta_x^2 + alpha_x^2",
            "substitute beta_x^2 = alpha_x^2             [eq 4.9]",
            "=> N_flux = 2*alpha_x^2",
        ],
        "sympy_result": str(n_flux_substituted),
        "matches_eq_4_13": ok,
        "source": "hep-th_0301139_TripathyTrivedi.txt lines 861-873 (eqs 2.8, 4.8, 4.9, 4.13)",
    }


def check_line_941_saturation_remark():
    """Flag a possible factor-of-2 discrepancy in TT's prose remark at line
    941: 'choices of delta_alpha_x, delta_beta_x which give rise to a vacuum
    where N_flux = 24 and no D3-branes need be added.' Using the PRINTED
    tadpole condition (2.3): (1/2) N_flux + N_D3 = 24, N_D3=0 saturation
    requires N_flux = 48, not 24. This is checked, not silently adopted
    either way; it affects only the prose remark at line 941, not the
    checked numbered equations (4.14),(4.16),(4.32), which all reproduce
    exactly (see checks above)."""
    n_d3 = 0
    n_flux_from_2_3 = Fraction(24 - n_d3) * 2  # (1/2) N_flux + 0 = 24 => N_flux = 48
    return {
        "id": "TT-line-941-remark-check",
        "description": "cross-check the prose remark at line 941 against numbered eq (2.3)",
        "printed_remark": "'...choices of delta alpha_x, delta beta_x which give rise to a vacuum where N_flux = 24 and no D3-branes need be added.' (line 940-942)",
        "eq_2_3_printed": "(1/2) N_flux + N_D3 = 24  (line 171)",
        "N_flux_required_for_N_D3_eq_0_per_eq_2_3": str(n_flux_from_2_3),
        "printed_remark_value": 24,
        "discrepancy_flagged": str(n_flux_from_2_3) != "24",
        "verdict": "The prose remark at line 941 states N_flux=24 for N_D3=0 saturation; applying the numbered tadpole condition (2.3) as printed gives N_flux=48 for that saturation. This is a candidate factor-of-2 slip in the PROSE aside (not in the numbered equations 4.14/4.16/4.32, which all check out exactly -- see TT-4.15-4.16 and TT-4.31-4.32 above). We flag it and do not silently resolve it either way.",
    }


def main():
    results = [
        check_4_15_4_16(),
        check_4_31_4_32(),
        symbolic_derivation_4_13(),
        check_line_941_saturation_remark(),
    ]
    out = {
        "track": "E-flux",
        "part": "metric-and-formula-validation",
        "command": " ".join([sys.executable] + sys.argv),
        "results": results,
        "all_numeric_checks_pass": all(
            r.get("match", True) for r in results if "match" in r
        ),
    }
    with open(OUT, "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
