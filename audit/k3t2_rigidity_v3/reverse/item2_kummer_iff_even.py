"""
Item 2 -- reverse-tests DualScaleDyons.WhichK3.kummer_iff_even beyond its Lean range.

Lean (WhichK3.lean, `decide +kernel`) checks, for every D in 0..200 and every reduced form
(a,b,c) of discriminant -D: isKummerForm(a,b,c) == (a,b,c all even). This script re-derives
reducedForms and the even/odd characterisation independently (common.py) and checks the same
biconditional for every D up to D_max_kummer_iff_even, plus a negative control that perturbs the
selecting condition itself (require only a,c even, not b) and shows it is NOT equivalent.

Run:
    cd audit/k3t2_rigidity_v3/reverse && \
    /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python item2_kummer_iff_even.py
"""
import json
from pathlib import Path

from common import reduced_forms, is_kummer_form

HERE = Path(__file__).resolve().parent
LEAN_RANGE_MAX = 200


def all_even(f):
    a, b, c = f
    return a % 2 == 0 and b % 2 == 0 and c % 2 == 0


def a_c_even_only(f):
    a, b, c = f
    return a % 2 == 0 and c % 2 == 0


def main() -> None:
    inputs = json.loads((HERE / "inputs.json").read_text())
    D_max = next(i["value"] for i in inputs if i["name"] == "D_max_kummer_iff_even")

    mismatches = []
    control_mismatches = []
    for D in range(0, D_max + 1):
        for f in reduced_forms(D):
            if is_kummer_form(f) != all_even(f):
                mismatches.append((D, f))
            if D > LEAN_RANGE_MAX and is_kummer_form(f) != a_c_even_only(f):
                control_mismatches.append((D, f))

    result = {
        "theorem": "DualScaleDyons.WhichK3.kummer_iff_even",
        "lean_range": [0, LEAN_RANGE_MAX],
        "scanned_range": [0, D_max],
        "holds_for_all_scanned": len(mismatches) == 0,
        "num_mismatches": len(mismatches),
        "sample_mismatches": mismatches[:5],
        "negative_control": {
            "description": "Perturbs the SAME selecting condition (which of a,b,c must be even) "
                            "by dropping the requirement on b, checked only beyond the Lean range "
                            "(D>200); this must NOT be equivalent to isKummerForm, i.e. some "
                            "reduced form must have a,c even, b odd.",
            "num_forms_where_control_disagrees": len(control_mismatches),
            "control_disagrees_as_expected": len(control_mismatches) > 0,
            "sample": control_mismatches[:5],
        },
    }
    (HERE / "item2_results.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2)[:2000])


if __name__ == "__main__":
    main()
