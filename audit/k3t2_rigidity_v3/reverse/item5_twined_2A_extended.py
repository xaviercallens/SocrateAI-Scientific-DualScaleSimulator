"""
Item 5 -- reverse-tests DualScaleMoonshine.Twining.twined_div24 / twined_2A beyond Lean's q^9.

Lean (Twining.lean, `decide`) computes, from CDH's frame-shape formula (eq. 4.18, Table 3, Table
14) H_g = (chi_g/24) H + F_g/eta^3 with F_2A = -16*Lambda_2, chi_2A = 8, and checks through q^9
that (a) every coefficient of 24*H_2A is divisible by 24 (twined_div24) and (b) the quotient
matches CDH's printed Table 20 column (twined_2A).

This script reimplements the whole pipeline independently (mulTrunc/divTrunc/eta3/sigma/f2Coeff/
numer from QSeries.lean's definitions, lambda24/twined24 from Twining.lean's definitions -- read
from the Lean source, not executed) in plain Python exact integers, and:
  1. regression-checks the reimplementation against table2A through q^9 (sanity against the value
     Lean itself proves, read verbatim from Twining.lean's `table2A` -- see inputs.json);
  2. extends the STRUCTURAL claim (divisibility by 24) from q^9 to q^N_max_twined_2A, where no
     literature table is available to compare against, so the extension is a structural
     prediction, not a further table match;
  3. characterises exactly which frame-shape coefficients c (in F_2A = c*Lambda_2) preserve the
     divisibility-by-24 property. Working through div_trunc by hand shows the q^0 coefficient is
     combined[0] = chi*(-2) + c*Nlev*(Nlev-1) = -16 + 2c, eta3's constant term is 1, and every
     later quotient coefficient is a linear combination of earlier ones plus (chi*numer(n) +
     c*24*Nlev*(...)), which is automatically a multiple of 24 for every integer c once n>=1 --
     so divisibility-by-24-for-all-n reduces to the single condition combined[0] % 24 == 0, i.e.
     c == 8 (mod 12). This is checked directly (not assumed) below over a spread of c values,
     replacing the originally planned single-c "negative control", which a first attempt (c=-15)
     could not distinguish from a genuine selection principle: -15 fails, but so would infinitely
     many wrong values, and infinitely many OTHER values (c=-4, c=8, ...) pass exactly as -16
     does. The check is windowed from q^0, not q^10, because that is where the c-dependence
     actually lives.

Run:
    cd audit/k3t2_rigidity_v3/reverse && \
    /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python item5_twined_2A_extended.py
"""
import json
from pathlib import Path

from common import sigma

HERE = Path(__file__).resolve().parent
LEAN_RANGE_MAX = 9
# Read verbatim from DualScaleMoonshine/Twining.lean's own `table2A` definition (its own
# transcription of CDH Table 20, column 2A) -- not retyped from memory of the paper.
TABLE2A_FROM_LEAN_SOURCE = [-2, -6, 14, -28, 42, -56, 86, -138, 188, -238]


def mul_trunc(N: int, a: list[int], b: list[int]) -> list[int]:
    def get(lst, i):
        return lst[i] if 0 <= i < len(lst) else 0
    return [sum(get(a, k) * get(b, n - k) for k in range(n + 1)) for n in range(N + 1)]


def div_trunc(N: int, r: list[int], p: list[int]) -> list[int]:
    def get(lst, i):
        return lst[i] if 0 <= i < len(lst) else 0
    acc: list[int] = []
    for n in range(N + 1):
        val = get(r, n) - sum(get(acc, k) * get(p, n - k) for k in range(n))
        acc.append(val)
    return acc


def one_minus_q(N: int, n: int) -> list[int]:
    out = [0] * (N + 1)
    out[0] = 1
    if 0 < n <= N:
        out[n] = -1
    return out


def eta3(N: int) -> list[int]:
    acc = one_minus_q(N, 0)
    for i in range(N):
        f = one_minus_q(N, i + 1)
        acc = mul_trunc(N, mul_trunc(N, mul_trunc(N, acc, f), f), f)
    return acc


def f2_coeff(n: int) -> int:
    total = 0
    for s in range(1, 2 * n + 1):
        if (2 * n) % s != 0:
            continue
        r = 2 * n // s
        if not (s < r):
            continue
        if (r - s) % 2 != 1:
            continue
        total += (-s if r % 2 == 1 else s)
    return total


def numer(N: int) -> list[int]:
    out = [0] * (N + 1)
    out[0] = -2
    for n in range(1, N + 1):
        out[n] = 48 * (sigma(n) + f2_coeff(n))
    return out


def lambda24(N: int, Nlev: int, c: int) -> list[int]:
    out = []
    for n in range(N + 1):
        if n == 0:
            out.append(c * Nlev * (Nlev - 1))
        else:
            term = sigma(n) - (Nlev * sigma(n // Nlev) if n % Nlev == 0 else 0)
            out.append(c * 24 * Nlev * term)
    return out


def twined24(N: int, chi: int, Nlev: int, c: int) -> list[int]:
    num = numer(N)
    lam = lambda24(N, Nlev, c)
    combined = [chi * num[i] + lam[i] for i in range(N + 1)]
    return div_trunc(N, combined, eta3(N))


def main() -> None:
    inputs = json.loads((HERE / "inputs.json").read_text())
    N_max = next(i["value"] for i in inputs if i["name"] == "N_max_twined_2A")
    params = next(i["value"] for i in inputs if i["name"] == "twined_2A_chi_Nlev_c_from_CDH")
    chi, Nlev, c = params["chi"], params["Nlev"], params["c"]

    series_full = twined24(N_max, chi, Nlev, c)
    series9 = series_full[: LEAN_RANGE_MAX + 1]

    regression_ok = [x // 24 for x in series9] == TABLE2A_FROM_LEAN_SOURCE and all(
        x % 24 == 0 for x in series9
    )

    div_beyond = [i for i in range(LEAN_RANGE_MAX + 1, N_max + 1) if series_full[i] % 24 != 0]

    # Characterise divisibility-by-24 as a function of c, checked from q^0 (not windowed to
    # q>9): does it hold for EVERY coefficient in [0, N_max], for a spread of c values including
    # some predicted (by the -16+2c == 0 mod 24 hand-derivation above) to pass and some to fail?
    c_candidates = [-16, -15, -4, 8, -3, 20, 32]
    c_scan = []
    for c_val in c_candidates:
        s = twined24(N_max, chi, Nlev, c_val)
        bad = [i for i in range(N_max + 1) if s[i] % 24 != 0]
        c_scan.append({
            "c": c_val,
            "predicted_by_congruence_c_eq_8_mod_12": (c_val - 8) % 12 == 0,
            "divisible_by_24_for_all_n_in_range": len(bad) == 0,
            "first_failing_index": (bad[0] if bad else None),
        })
    congruence_confirmed = all(
        row["predicted_by_congruence_c_eq_8_mod_12"] == row["divisible_by_24_for_all_n_in_range"]
        for row in c_scan
    )

    result = {
        "theorem": "DualScaleMoonshine.Twining.twined_div24 / twined_2A",
        "lean_range": [0, LEAN_RANGE_MAX],
        "scanned_range": [0, N_max],
        "regression_matches_table2A_through_q9": regression_ok,
        "quotient_series_full": [x // 24 for x in series_full],
        "divisibility_by_24_holds_beyond_q9": len(div_beyond) == 0,
        "failing_indices_beyond_q9": div_beyond,
        "c_dependence_scan": {
            "description": "Divisibility-by-24 (checked from q^0, the coefficient that "
                            "actually carries the c-dependence) as a function of c in "
                            "F_2A = c*Lambda_2, for c = -16 (CDH's value) plus six others.",
            "rows": c_scan,
            "matches_hand_derived_congruence_c_eq_8_mod_12": congruence_confirmed,
        },
    }
    (HERE / "item5_results.json").write_text(json.dumps(result, indent=2) + "\n")
    printable = {k: v for k, v in result.items() if k != "quotient_series_full"}
    print(json.dumps(printable, indent=2))
    print("quotient series (q^0 .. q^N_max):", result["quotient_series_full"])


if __name__ == "__main__":
    main()
