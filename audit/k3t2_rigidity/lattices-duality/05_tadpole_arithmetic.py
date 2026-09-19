"""
Track C, item (5): flux tadpole arithmetic.

chi(K3) from Betti numbers b0,b1,b2,b3,b4 = 1,0,22,0,1 (given, standard for K3):
    chi(K3) = sum (-1)^i b_i = b0 - b1 + b2 - b3 + b4.

chi(K3 x K3) via Kunneth: chi(X x Y) = chi(X) * chi(Y) (multiplicativity of
Euler characteristic under products -- derived here from Kunneth's formula:
b_k(X x Y) = sum_{i+j=k} b_i(X) b_j(Y), summed with signs, which factors as
a product of the two Poincare-polynomial-at-(-1) evaluations). Then
chi(K3 x K3) / 24 is reported as an exact Fraction.

chi(K3 x T2) via Kunneth: T2 (a 2-torus) has Betti numbers b0,b1,b2 = 1,2,1,
so chi(T2) = 1 - 2 + 1 = 0. Kunneth gives chi(K3 x T2) = chi(K3) * chi(T2).
"""
import json
import sys
from fractions import Fraction

sys.path.insert(0, "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity/lattices-duality")


def euler_char(betti):
    return sum(((-1) ** i) * b for i, b in enumerate(betti))


def poincare_at_minus1_product_check(betti_x, betti_y):
    """Verify chi(X x Y) = chi(X) * chi(Y) directly from the Kunneth Betti
    numbers b_k(XxY) = sum_{i+j=k} b_i(X) b_j(Y), computed explicitly (not
    assumed), then Euler-characteristic-summed."""
    max_k = (len(betti_x) - 1) + (len(betti_y) - 1)
    betti_prod = []
    for k in range(max_k + 1):
        s = 0
        for i in range(len(betti_x)):
            j = k - i
            if 0 <= j < len(betti_y):
                s += betti_x[i] * betti_y[j]
        betti_prod.append(s)
    chi_direct = euler_char(betti_prod)
    chi_factored = euler_char(betti_x) * euler_char(betti_y)
    return betti_prod, chi_direct, chi_factored, (chi_direct == chi_factored)


# K3 Betti numbers (given, standard: b0=1,b1=0,b2=22,b3=0,b4=1; rank-22 H^2
# with signature (3,19), matching the K3 lattice signature computed in item 2)
BETTI_K3 = [1, 0, 22, 0, 1]
chi_K3 = euler_char(BETTI_K3)

# T2 = S1 x S1 Betti numbers: b0=1, b1=2 (two independent 1-cycles), b2=1
BETTI_T2 = [1, 2, 1]
chi_T2 = euler_char(BETTI_T2)

betti_K3xK3, chi_K3xK3_direct, chi_K3xK3_factored, kunneth_check_K3K3 = \
    poincare_at_minus1_product_check(BETTI_K3, BETTI_K3)

betti_K3xT2, chi_K3xT2_direct, chi_K3xT2_factored, kunneth_check_K3T2 = \
    poincare_at_minus1_product_check(BETTI_K3, BETTI_T2)

chi_K3xK3_over_24 = Fraction(chi_K3xK3_direct, 24)
chi_K3xT2_over_24 = Fraction(chi_K3xT2_direct, 24)

# negative control: perturb one Betti number of K3 by +1 (breaking the true
# K3 Hodge numbers) and confirm chi changes (i.e. the computation is
# sensitive to the input, not a hardcoded constant)
BETTI_K3_PERTURBED = [1, 0, 23, 0, 1]  # b2 = 23 instead of 22 (not K3)
chi_K3_perturbed = euler_char(BETTI_K3_PERTURBED)
negative_control_passes = (chi_K3_perturbed != chi_K3)

results = {
    "track": "C",
    "item": "5_tadpole_arithmetic",
    "method": "Euler characteristic from given Betti numbers via alternating sum; "
               "Kunneth formula for products verified two ways (direct Betti-number "
               "convolution vs. the chi(X)*chi(Y) factorization) and shown to agree, "
               "rather than assuming multiplicativity.",
    "K3": {
        "betti_numbers_b0_b4": BETTI_K3,
        "euler_characteristic": chi_K3,
        "expected": {"value": 24, "source": "standard K3 fact, from memory, unverified prior to computation"},
    },
    "T2": {
        "betti_numbers_b0_b2": BETTI_T2,
        "euler_characteristic": chi_T2,
        "expected": {"value": 0, "source": "standard T^2 fact (torus Euler char = 0), from memory, unverified prior to computation"},
    },
    "K3_x_K3": {
        "betti_numbers_via_kunneth": betti_K3xK3,
        "chi_direct_from_kunneth_betti": chi_K3xK3_direct,
        "chi_via_multiplicativity_chiX_times_chiY": chi_K3xK3_factored,
        "two_methods_agree": kunneth_check_K3K3,
        "chi_over_24_exact": str(chi_K3xK3_over_24),
    },
    "K3_x_T2": {
        "betti_numbers_via_kunneth": betti_K3xT2,
        "chi_direct_from_kunneth_betti": chi_K3xT2_direct,
        "chi_via_multiplicativity_chiX_times_chiY": chi_K3xT2_factored,
        "two_methods_agree": kunneth_check_K3T2,
        "chi_over_24_exact": str(chi_K3xT2_over_24),
        "note": "chi(K3 x T2) = chi(K3) * chi(T2) = chi(K3) * 0 = 0, since T2 has "
                "vanishing Euler characteristic (it is a group manifold / has a "
                "nowhere-vanishing vector field); this is exact, not a numerical "
                "coincidence, and follows directly from Kunneth applied to the "
                "given Betti numbers.",
    },
    "negative_control": {
        "description": "perturbing K3's b2 from 22 to 23 (not the real K3 Hodge number) "
                        "must change chi",
        "chi_with_perturbed_betti": chi_K3_perturbed,
        "chi_with_true_betti": chi_K3,
        "control_passes": negative_control_passes,
    },
}

out_path = "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity/lattices-duality/05_tadpole_arithmetic_result.json"
with open(out_path, "w") as f:
    json.dump(results, f, indent=2)

print(json.dumps(results, indent=2))
