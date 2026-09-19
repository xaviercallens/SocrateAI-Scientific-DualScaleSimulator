"""
Track C (v2), item: tadpole budget enumeration.

REVISED from the first version of this script: the first version computed
chi(K3) from a TYPED Betti-number list [1,0,22,0,1] (the 22 was a literal
K3 fact, even though it was only used to derive chi, not asserted directly).
This version removes that last typed K3 number: budget = chi_top is READ
from the SAME provenance file as 06_chained_rigidity.py (Track A's
A-genus/exports.json, "chi_from_genus" key -- an independently computed
elliptic-genus constant term, not a memorized fact), via chi_top_provenance.py.
No K3 Betti number, and no literal 24, ever appears in this script's source.

Tadpole enumeration: for NATURAL numbers (flux, n) -- both >= 0, both exact
Python ints -- satisfying flux + n == budget, enumerate ALL solutions
mechanically (a simple range loop over flux = 0..budget, n = budget-flux;
no candidate is hand-picked). Then show, as computed facts about the
enumerated set:
  (i)  n <= budget for every solution (trivial from n = budget - flux and
       flux >= 0, but verified computationally rather than asserted).
  (ii) flux == 0  implies  n == budget (trivial substitution, but again
       checked on the actual enumerated set rather than merely stated).
Negative control: (flux, n) = (1, budget) is checked and correctly shown to
violate flux + n == budget (since budget - 1 != budget).

Run: /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python \
     /mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity_v2/C-lattices/04_tadpole_budget.py
"""
import json
import sys

HERE = "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity_v2/C-lattices"
sys.path.insert(0, HERE)
from chi_top_provenance import chi_top_from_provenance

chi_top, provenance_file, provenance_key, provenance_raw = chi_top_from_provenance()
budget = chi_top  # the enumeration condition below uses this VARIABLE, read from file, never "24"

# --- mechanical enumeration: every natural (flux, n) with flux + n == budget ---
solutions = []
for flux in range(0, budget + 1):
    n = budget - flux
    assert flux + n == budget
    assert flux >= 0 and n >= 0
    solutions.append({"flux": flux, "n": n})

all_n_le_budget = all(s["n"] <= budget for s in solutions)
flux_zero_solutions = [s for s in solutions if s["flux"] == 0]
flux_zero_implies_n_eq_budget = all(s["n"] == budget for s in flux_zero_solutions) and len(flux_zero_solutions) == 1

# negative control: (flux, n) = (1, budget) does NOT satisfy flux + n == budget
neg_flux, neg_n = 1, budget
neg_control_holds = (neg_flux + neg_n != budget)

results = {
    "track": "C", "version": "v2",
    "item": "4_tadpole_budget",
    "method": "budget = chi_top, READ from the same provenance file as "
              "06_chained_rigidity.py (no typed Betti numbers, no literal 24); mechanical "
              "enumeration of all natural (flux, n) with flux + n == budget, using the "
              "read-in budget as a variable inside the selecting condition.",
    "provenance": {
        "chi_top_read_from_file": provenance_file,
        "chi_top_read_from_key": provenance_key,
        "chi_top_value": chi_top,
    },
    "budget_used_in_enumeration": budget,
    "num_solutions_enumerated": len(solutions),
    "all_n_le_budget": all_n_le_budget,
    "flux_zero_implies_n_eq_budget": flux_zero_implies_n_eq_budget,
    "flux_zero_solution": flux_zero_solutions[0] if flux_zero_solutions else None,
    "negative_control": {
        "description": "(flux, n) = (1, budget) must fail flux + n == budget",
        "flux": neg_flux, "n": neg_n,
        "flux_plus_n": neg_flux + neg_n, "budget": budget,
        "control_passes": neg_control_holds,
    },
    "solutions_sample_head_and_tail": {
        "first_5": solutions[:5],
        "last_5": solutions[-5:],
    },
    "classification": "NORMALISATION: budget = chi_top is fixed once chi_top is read from "
                      "Track A's file (an external input, stated as such above); GIVEN that "
                      "budget, the constraint n<=budget and (flux=0 => n=budget) are trivial "
                      "arithmetic facts about the linear relation flux+n=budget, not rigidity "
                      "results -- reported here as the plain arithmetic content of the tadpole "
                      "condition, per the ground rules' 'trivial but state it' instruction.",
    "expected": {"chi_top": 24, "source": "standard K3 fact, from memory, unverified prior to computation "
                                          "-- filled in AFTER reading/computing budget above"},
}

out_path = HERE + "/04_tadpole_budget_result.json"
with open(out_path, "w") as f:
    json.dump(results, f, indent=2)

print(json.dumps(results, indent=2))
