#!/usr/bin/env python3
"""
ROUND 3, step 4 (REDUCE, exploratory candidate for mu_sym): the harness's
NEW CONTEXT states that what survives in DarkEnergyScale.lean after the
CKN dual-scale hypothesis was falsified at cosmological scales is "the
self-dual length sqrt(l_P * L_H) and the MVV micron-scale dark
dimension" (tier L/C). This script computes that self-dual length EXACTLY
from the two already-verified Lean constants (planckLength_m,
SelfDualCutoff.lean:96; hubbleRadius_m, SelfDualCutoff.lean:100; both
pinned v3.9.0 per round1's independent verification) and checks whether
it lands near the "47 micron" short-range-gravity scale this loop's own
ABSENT list already names (fifth_force_screening: no dataset).

This is NOT a derivation of mu_sym (mu_sym is a dimensionless BVP
parameter internal to the symmetron ODE, with no length calibration
anywhere in this model -- the same missing-calibration gap that already
forces every dark-energy distance in this harness to be reported
dimensionless, H0*D). It IS a candidate physical length scale that a
FUTURE calibration bridge (mu_sym <-> 1/length, in whatever units the
model would need to acquire) could be tested against, IF a real
short-range-gravity/fifth-force dataset is ever fetched (still ABSENT).
Reported here as tier C, explicitly not counted as a reduction.

Writes: mu_sym_length_bridge_report.json
"""
import json
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))

PLANCK_LENGTH_M = 1.616255e-35   # SelfDualCutoff.lean:96 (Lean def, quoted verbatim)
HUBBLE_RADIUS_M = 1.3672e26      # SelfDualCutoff.lean:100 (Lean def, quoted verbatim, = c/H0)

self_dual_length_m = math.sqrt(PLANCK_LENGTH_M * HUBBLE_RADIUS_M)
self_dual_length_micron = self_dual_length_m * 1e6

# Reference short-range-gravity torsion-balance scale this loop's own
# ABSENT list already names ("short-range gravity bound (47 micron)");
# NOT re-derived or fetched here -- reported as the named comparison
# point, tier L provenance of that number itself is NOT independently
# verified in this session (no dataset fetched), stated honestly.
NAMED_REFERENCE_MICRON = 47.0

report = {
    "generated": "2026-09-18",
    "inputs_quoted_verbatim": {
        "planckLength_m": {"value": PLANCK_LENGTH_M, "file_line": "SelfDualCutoff.lean:96"},
        "hubbleRadius_m": {"value": HUBBLE_RADIUS_M, "file_line": "SelfDualCutoff.lean:100"},
    },
    "self_dual_length_m": self_dual_length_m,
    "self_dual_length_micron": self_dual_length_micron,
    "named_reference_short_range_gravity_micron": NAMED_REFERENCE_MICRON,
    "ratio_computed_to_reference": self_dual_length_micron / NAMED_REFERENCE_MICRON,
    "agreement_within_10pct": bool(abs(self_dual_length_micron / NAMED_REFERENCE_MICRON - 1.0) < 0.10),
    "identification_tier": "C",
    "status": "NOT a reduction. mu_sym has no length calibration anywhere in workshopcosmo.py's symmetron ODE or param_loop_sim.py (it is a bare dimensionless coefficient of the BVP); this number is a candidate target for a future calibration, not a value FOR mu_sym.",
    "missing_bridge": "A formula relating mu_sym (dimensionless, current model) to 1/self_dual_length_m (or any other physical inverse-length) does not exist in this codebase; would need to be constructed and then tested against a real torsion-balance/short-range-gravity dataset, which is itself ABSENT from data/real/MANIFEST.json (fifth_force_screening entry).",
    "absent_dataset_needed_to_test": "fifth_force_screening (data/real/MANIFEST.json: status ABSENT)",
}

self_dual_length_micron_rounded = round(self_dual_length_micron, 3)
print(json.dumps(report, indent=2))
with open(os.path.join(HERE, "mu_sym_length_bridge_report.json"), "w") as f:
    json.dump(report, f, indent=2)
