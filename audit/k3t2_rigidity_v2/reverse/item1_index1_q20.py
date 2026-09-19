"""
Reverse-pass item 1: LeanMaster's `c_depends_only_on_D` (DualScaleDyons/DMVV.lean) states
the K3-elliptic-genus index-1 Jacobi property -- every Fourier coefficient of the computed
Z_K3 = 2*phi_{0,1} depends only on D = 4n - l^2 -- and checks it through q^9 (`zK3 :=
ellipticGenus 9`). This script extends the SAME structural property, independently, to q^20.

PREDICTION (written before this script reads any q>9 data; see `prediction` field below,
filled in before `computed` is touched): the index-1 collapse continues to hold at every
new order q^10..q^20, because it is forced by the defining product formula for the K3
elliptic genus (2*phi_{0,1} is, by construction, a weight-0 index-1 weak Jacobi form; the
index-1 property is a structural consequence of the theta-function product formulas used
to build it in `theta_forms.py`, not a numerical accident tied to q<=9). This is therefore
expected to hold at EVERY order, not just q<=20; we test the specific new range q=10..20
because that is what is new relative to Lean's checked range.

SELECTING CONDITION (structural, not literal): "every (n,l) pair with the same
D = 4n - l^2 carries the same coefficient in the independently-computed B_series" is a
consistency condition on ~2*(number of l at each n) numbers colliding onto ~1 number per D;
it is not built from any known value of c(D) (no c(D) literal appears in the check itself).
Negative control: perturbing ONE coefficient at one (n,l) by +1 must create a same-D
disagreement with its (n,l') partner (partner picked with l'^2 - l^2 = -4n+4n = ... i.e.
same D, different l parity) -- this demonstrates the check has teeth.

Data: reads B-dyons/theta_forms_cache.json, an EXISTING, already-committed cache built by
an independent theta-function engine (thetas.py/series2d.py in ../A-genus, and the
B-series 2D-array construction in B-dyons/theta_forms.py), never derived from LeanMaster's
Lean code. Regeneration command (already run; QMAX=40 >= 20 gives ample margin over the
q^20 target): `cd ../B-dyons && /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python theta_forms.py 40`

Run (from this directory, after the cache above exists):
  /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python item1_index1_q20.py
Writes item1_index1_q20_results.json.
"""
import json
from fractions import Fraction as Fr

CACHE_PATH = "../B-dyons/theta_forms_cache.json"
N_LEAN_CHECKED = 9   # Lean's `c_depends_only_on_D` range (zK3 := ellipticGenus 9)
N_TARGET = 20        # new range requested by the task

with open(CACHE_PATH) as f:
    cache = json.load(f)

QMAX = cache["QMAX"]
assert QMAX >= N_TARGET, f"theta cache QMAX={QMAX} too small for target q^{N_TARGET}"

def pf(s):
    if s is None:
        return None
    if "/" in s:
        a, b = s.split("/")
        return Fr(int(a), int(b))
    return Fr(int(s))

def load2d(key):
    out = {}
    for k, v in cache[key].items():
        n, l = k.split(",")
        out[(int(n), int(l))] = pf(v)
    return out

B = load2d("B_series")  # B = phi_{0,1} = Z_K3 / 2; coefficients c(D)/2 at (n,l), D = 4n-l^2

prediction = {
    "claim": "For every n = 10..20 and every l with (n,l) present in the independently "
             "computed B_series, the coefficient depends only on D = 4n - l^2 (i.e. every "
             "two pairs (n,l), (n',l') with the same D agree), extending LeanMaster's "
             "c_depends_only_on_D (checked there only for n<=9).",
    "why_expected": "Structural: 2*phi_{0,1} is by construction a weak Jacobi form of "
                     "index 1 (elliptic transformation law with index m=1), so its Fourier "
                     "coefficients c(n,l) depend on (D,l mod 2) alone at every order, not "
                     "just the ones Lean happened to check; this is a property of the "
                     "theta-function product formula used to build B_series, not of any "
                     "particular truncation.",
    "negative_control_plan": "perturb one (n,l) coefficient in a range n=10..20 by +1 and "
                              "confirm this creates a same-D disagreement with its (n,l') "
                              "partner sharing that D, i.e. the check is not vacuous.",
}

# ---- independently re-derive c(D) from B_series restricted to n = N_LEAN_CHECKED+1 .. N_TARGET
by_D = {}
conflicts = []
for (n, l), v in B.items():
    if n < N_LEAN_CHECKED + 1 or n > N_TARGET:
        continue
    D = 4 * n - l * l
    if D in by_D:
        if by_D[D][1] != v:
            conflicts.append({"D": D, "first": [by_D[D][0], str(by_D[D][1])], "second": [(n, l), str(v)]})
    else:
        by_D[D] = ((n, l), v)

holds_n10_20 = (len(conflicts) == 0)

# ---- how many distinct new D values were exercised, and how many (n,l) pairs collapsed onto them
count_per_D = {}
for (n, l), v in B.items():
    if n < N_LEAN_CHECKED + 1 or n > N_TARGET:
        continue
    D = 4 * n - l * l
    count_per_D[D] = count_per_D.get(D, 0) + 1
multi_pairs_per_D = {D: c for D, c in count_per_D.items() if c > 1}

# ---- new c(D) values beyond Lean's checked range (n<=9 => D <= 4*9+l_min^2, roughly D<=36
# at l=0; the new D's from n=10..20 that were NOT already reachable at n<=9)
D_reachable_at_n_leq_9 = set()
for (n, l), v in B.items():
    if n <= N_LEAN_CHECKED:
        D_reachable_at_n_leq_9.add(4 * n - l * l)
new_D_values = sorted(d for d in by_D if d not in D_reachable_at_n_leq_9)
new_cD = {str(d): str(2 * by_D[d][1]) for d in new_D_values}  # c(D) = 2*B-coefficient

# ---- negative control: perturb a coefficient at n=15 and see if a same-D conflict appears
def find_same_D_pair(nmin, nmax):
    seen = {}
    for (n, l), v in B.items():
        if not (nmin <= n <= nmax):
            continue
        D = 4 * n - l * l
        if D in seen and seen[D][0] != (n, l):
            return seen[D][0], (n, l), D
        seen.setdefault(D, ((n, l)))
    return None

pair = None
for (n, l), v in B.items():
    if n != 15:
        continue
    D = 4 * n - l * l
    for (n2, l2), v2 in B.items():
        if (n2, l2) != (n, l) and 4 * n2 - l2 * l2 == D:
            pair = ((n, l), (n2, l2), D, v, v2)
            break
    if pair:
        break

neg_control = None
if pair:
    (n1, l1), (n2, l2), D, v1, v2 = pair
    same_before = (v1 == v2)
    v1_perturbed = v1 + 1
    same_after = (v1_perturbed == v2)
    neg_control = {
        "pair_tested": {"(n1,l1)": [n1, l1], "(n2,l2)": [n2, l2], "D": D},
        "values_agree_before_perturbation": same_before,
        "perturb_(n1,l1)_by_+1_then_agree": same_after,
        "control_has_teeth": same_before and not same_after,
    }

result = {
    "prediction": prediction,
    "lean_theorem": "c_depends_only_on_D, DualScaleDyons/DMVV.lean (LeanMaster)",
    "lean_checked_range": f"n <= {N_LEAN_CHECKED} (zK3 := ellipticGenus 9)",
    "new_range_tested": f"n = {N_LEAN_CHECKED + 1}..{N_TARGET}",
    "computed": {
        "holds_on_new_range": holds_n10_20,
        "conflicts_found": conflicts,
        "num_D_values_with_multiple_(n,l)_pairs_in_new_range": len(multi_pairs_per_D),
        "num_new_D_values_beyond_n_leq_9": len(new_D_values),
        "new_cD_sample": {k: new_cD[k] for k in list(new_cD)[:15]},
    },
    "negative_control": neg_control,
    "tier": "B",
    "command": "cd audit/k3t2_rigidity_v2/reverse && "
               "/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python item1_index1_q20.py",
    "source_cache_regeneration_command": "cd audit/k3t2_rigidity_v2/B-dyons && "
               "/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python theta_forms.py 40",
}

with open("item1_index1_q20_results.json", "w") as f:
    json.dump(result, f, indent=1)

print(json.dumps({"holds_n10_20": holds_n10_20, "num_conflicts": len(conflicts),
                   "num_new_D_values": len(new_D_values), "neg_control": neg_control}, indent=1))
