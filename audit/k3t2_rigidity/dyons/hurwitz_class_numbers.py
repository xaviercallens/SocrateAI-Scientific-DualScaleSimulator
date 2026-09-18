"""
Hurwitz class numbers H(D), D>0, D = 0 or 3 mod 4, computed by counting reduced
positive-definite integral binary quadratic forms (a,b,c), disc = b^2-4ac = -D
(imprimitive forms included -- this automatically sums over all f^2|D, which is
the definition of the Hurwitz class number), reduced by:
    -a < b <= a <= c,  and b >= 0 if a==c or a==b.
weighted 1, except weight 1/2 for the (unique, when it occurs) form equivalent
to a*(x^2+y^2) [i.e. b=0,a=c], and weight 1/3 for the form equivalent to
a*(x^2+xy+y^2) [i.e. a=b=c].
H(0) := -1/12 by convention (not from this counting method; stated as given).
"""
import json
from fractions import Fraction as Fr

def H(D):
    if D == 0:
        return Fr(-1, 12)
    if D % 4 not in (0, 3):
        return Fr(0)
    total = Fr(0)
    forms = []
    amax = int((D / 3) ** 0.5) + 2
    for a in range(1, amax + 1):
        for b in range(-a + 1, a + 1):
            num = b * b + D
            if num % (4 * a) != 0:
                continue
            c = num // (4 * a)
            if c < a:
                continue
            if c == a and b < 0:
                continue
            if a == b and c > a:
                # a==b but not a==c: still a legitimately reduced form, generic weight
                pass
            # weight
            if b == 0 and a == c:
                w = Fr(1, 2)
            elif a == b and a == c:
                w = Fr(1, 3)
            else:
                w = Fr(1)
            total += w
            forms.append((a, b, c, str(w)))
    return total, forms

targets = [3, 4, 7, 8, 11, 12, 15]
result = {}
for D in targets:
    val, forms = H(D)
    result[D] = {"H": str(val), "reduced_forms": forms}

# sanity: H(0) by convention
result[0] = {"H": str(Fr(-1, 12)), "note": "by convention, not from the counting method"}

out = {str(k): v for k, v in sorted(result.items())}
with open("hurwitz_class_numbers_results.json", "w") as f:
    json.dump(out, f, indent=1)
print(json.dumps(out, indent=1))
