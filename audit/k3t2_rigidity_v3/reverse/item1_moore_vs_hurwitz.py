"""
Item 1 -- reverse-tests DualScaleDyons.WhichK3.moore_vs_hurwitz beyond its Lean range.

Lean (WhichK3.lean, `decide +kernel`) checks, for every D in 0..400:
    12 * nForms(D) == h12(D) + 6*[D = 4f^2] + 8*[D = 3f^2]      (for D >= 3, D % 4 in {0,3})
i.e. Moore's count of attractive-K3 / dyon-charge classes of discriminant -D agrees with the
Hurwitz class number 12H(D) except at the two self-dual-torus discriminants, where it disagrees
by exactly the extra-automorphism weight.

This script re-derives nForms and h12 from scratch (common.py, not imported from Lean) and checks
the SAME identity for every D up to D_max_moore_hurwitz (from inputs.json), then runs a negative
control: dropping the correction terms should make the identity fail on some D in the new range.

Run:
    cd audit/k3t2_rigidity_v3/reverse && \
    /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python item1_moore_vs_hurwitz.py
"""
import json
from pathlib import Path

from common import n_forms, h12, is_k_square

HERE = Path(__file__).resolve().parent
LEAN_RANGE_MAX = 400  # WhichK3.moore_vs_hurwitz's own range, for reporting only


def identity_holds(D: int) -> bool:
    if D < 3 or D % 4 in (1, 2):
        return True  # vacuously true guard, exactly as in the Lean statement
    lhs = 12 * n_forms(D)
    rhs = h12(D) + (6 if is_k_square(4, D) else 0) + (8 if is_k_square(3, D) else 0)
    return lhs == rhs


def identity_holds_no_correction(D: int) -> bool:
    if D < 3 or D % 4 in (1, 2):
        return True
    return 12 * n_forms(D) == h12(D)


def main() -> None:
    inputs = json.loads((HERE / "inputs.json").read_text())
    D_max = next(i["value"] for i in inputs if i["name"] == "D_max_moore_hurwitz")

    failures = []
    for D in range(0, D_max + 1):
        if not identity_holds(D):
            failures.append(D)

    beyond_lean_checked = D_max - LEAN_RANGE_MAX
    beyond_lean_failures = [D for D in failures if D > LEAN_RANGE_MAX]

    # Negative control: perturb the SAME parameter the identity is about (drop the two
    # correction terms) and confirm it now fails somewhere in the newly scanned range.
    control_failures = [D for D in range(LEAN_RANGE_MAX + 1, D_max + 1)
                         if not identity_holds_no_correction(D)]

    result = {
        "theorem": "DualScaleDyons.WhichK3.moore_vs_hurwitz",
        "lean_range": [0, LEAN_RANGE_MAX],
        "scanned_range": [0, D_max],
        "newly_scanned_count": beyond_lean_checked,
        "identity_holds_for_all_scanned": len(failures) == 0,
        "failures_full_range": failures,
        "failures_beyond_lean_range": beyond_lean_failures,
        "negative_control": {
            "description": "Same identity with the two correction terms (6*[D=4f^2], "
                            "8*[D=3f^2]) removed, scanned only beyond the Lean range.",
            "num_failures": len(control_failures),
            "control_fails_as_expected": len(control_failures) > 0,
            "sample_failures": control_failures[:10],
        },
    }
    (HERE / "item1_results.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2)[:2000])


if __name__ == "__main__":
    main()
