"""
Item 3 -- reverse-tests DualScaleDyons.WhichK3.most_attractive beyond its Lean range.

Lean (WhichK3.lean, `decide +kernel`) checks only D=3 and D=4: each carries exactly one reduced
form (a,b,c), i.e. exactly one attractive K3 of that discriminant. This script computes
nForms(D) for every valid discriminant D (D % 4 in {0,3}, D >= 3) up to D_max_class_number_one and
finds ALL D with nForms(D) == 1, using the SAME reducedForms code as item1/item2 (common.py) --
nothing about which D should come out is hard-coded into the scan.

Only AFTER that scan is a Tier-L "expected" list consulted (inputs.json's
heegner_discriminants_from_memory: the nine discriminants with primitive class number 1), to see
whether the computed list matches the literature prediction. A negative control checks that
neighbouring valid discriminants (nearest D % 4 in {0,3} on each side of every hit) do NOT have
nForms == 1, i.e. the property is not automatically true nearby.

Run:
    cd audit/k3t2_rigidity_v3/reverse && \
    /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python item3_class_number_one.py
"""
import json
from pathlib import Path

from common import n_forms

HERE = Path(__file__).resolve().parent
LEAN_RANGE = [3, 4]


def valid_discriminant(D: int) -> bool:
    return D >= 3 and D % 4 in (0, 3)


def neighbours(D: int, all_valid: list[int]) -> list[int]:
    idx = all_valid.index(D)
    out = []
    if idx > 0:
        out.append(all_valid[idx - 1])
    if idx + 1 < len(all_valid):
        out.append(all_valid[idx + 1])
    return out


def main() -> None:
    inputs = json.loads((HERE / "inputs.json").read_text())
    D_max = next(i["value"] for i in inputs if i["name"] == "D_max_class_number_one")
    heegner_expected = next(
        i["value"] for i in inputs if i["name"] == "heegner_discriminants_from_memory"
    )

    all_valid = [D for D in range(3, D_max + 1) if valid_discriminant(D)]
    computed_hits = [D for D in all_valid if n_forms(D) == 1]

    matches_literature = sorted(computed_hits) == sorted(
        [d for d in heegner_expected if d <= D_max]
    )

    control_failures = []
    for D in computed_hits:
        for nb in neighbours(D, all_valid):
            if n_forms(nb) == 1:
                control_failures.append({"hit": D, "neighbour": nb})

    result = {
        "theorem": "DualScaleDyons.WhichK3.most_attractive",
        "lean_range_D": LEAN_RANGE,
        "scanned_range_D": [3, D_max],
        "computed_nForms_eq_1": computed_hits,
        "expected_from_literature_FROM_MEMORY": sorted(
            d for d in heegner_expected if d <= D_max
        ),
        "computed_matches_literature_list": matches_literature,
        "negative_control": {
            "description": "For every D with nForms(D)==1, its two nearest valid "
                            "discriminants (D % 4 in {0,3}) must NOT also have nForms==1 -- "
                            "otherwise 'nForms==1' would not be a rare, selecting property.",
            "num_neighbour_violations": len(control_failures),
            "control_passes_as_expected": len(control_failures) == 0,
            "violations": control_failures,
        },
    }
    (HERE / "item3_results.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2)[:2500])


if __name__ == "__main__":
    main()
