#!/usr/bin/env python3
"""
Chain-check for the K3 x T2 rigidity v3 tracks.

Reads every value it compares from the tracks' own exports.json / results.json
/ inputs.json files under audit/k3t2_rigidity_v3/, using paths built from this
script's own location with pathlib (works from any clone location, never an
absolute /mnt/... path). No physics number is typed in this script: every
number reported is read out of a JSON file that a track's own script wrote.

For each of the six CHAIN CHECK links named in the v3 task:
  1. A (elliptic genus, Z(tau,0) at the declared k) vs D (GUDHI chi)
  2. B (DMVV product at z=0, Goettsche comparison) vs D (chi)
  3. C (Noether + Hodge index: chi_top, b2) vs D (chi, b2)
  4. E (n_singular consumed from D) vs D (n_singular)
  5. F (M24 cycle shapes + order) vs A (its own, separately re-derived, M24
     cycle shapes + order)
  6. F internal: T(A) / T(Km A) discriminants vs F's own binary-form
     enumeration
this script reports:
  - CONSISTENT / INCONSISTENT (values agree or not)
  - CROSS_METHOD (the two ends were computed by genuinely different methods
    with no file-level import between them) or POINTER (one end's script
    literally reads the other track's exports file, detected by grepping the
    committed .py source for the other track's directory name -- not asserted
    by hand)
  - the declared inputs (by name, from each track's own inputs.json) that the
    two ends share, if any, plus any note a track's own results.json already
    attaches admitting a value was tuned/solved to match the other track
    (e.g. B's own "k_solved_exactly_from_chi_D" field), which downgrades an
    apparent CROSS_METHOD agreement to a non-independent one even without a
    file-level POINTER.

A link whose two ends share a declared normalisation, or whose "independent"
value was in fact solved from the other end, is flagged as NOT independent
corroboration even when the numbers agree.

Run from a CLEAN CLONE:
  git clone /mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2 /tmp/k3t2v3_clone
  cd /tmp/k3t2v3_clone
  /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python \
      audit/k3t2_rigidity_v3/chain/chain_check.py
"""
import json
import re
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent            # .../audit/k3t2_rigidity_v3/chain
V3 = HERE.parent                                   # .../audit/k3t2_rigidity_v3
AUDIT = V3.parent                                  # .../audit
REPO = AUDIT.parent                                # repo root (works from any clone)

A_DIR = "audit/k3t2_rigidity_v3/A-genus"
B_DIR = "audit/k3t2_rigidity_v3/B-dyons"
C_DIR = "audit/k3t2_rigidity_v3/C-lattices"
D_DIR = "audit/k3t2_rigidity_v3/D-tda"
E_DIR = "audit/k3t2_rigidity_v3/E-flux"
F_DIR = "audit/k3t2_rigidity_v3/F-whichk3"


def load(relpath):
    p = REPO / relpath
    with open(p) as f:
        return json.load(f), str(Path(relpath))


def rat(x):
    """Exact rational from an int, a numeral string, or a numeral string
    embedded in text like 'chi(O)=2, chi_top=24 (...)' is NOT handled here;
    callers pass the numeral field itself."""
    return Fraction(str(x))


def input_names(inputs_json):
    if isinstance(inputs_json, list):
        return set(x.get("name") for x in inputs_json if isinstance(x, dict))
    if isinstance(inputs_json, dict):
        return set(inputs_json.keys())
    return set()


def py_files(reldir):
    return sorted((REPO / reldir).glob("*.py"))


def find_file_pointer(reldir, needle):
    """Grep this track's own committed .py scripts for a literal reference to
    another track's directory name (e.g. 'D-tda'). Returns the first script
    path that references it, or None. This is how POINTER is detected: by
    inspecting the actual source, not by asserting it."""
    hits = []
    for p in py_files(reldir):
        try:
            text = p.read_text()
        except OSError:
            continue
        if needle in text:
            hits.append(str(p.relative_to(REPO)))
    return hits


def extract_ints(s):
    return [int(x) for x in re.findall(r"-?\d+", str(s))]


def rigidity_entry(results_json, needle):
    """Find the rigidity-table entry (top-level 'rigidity' list) whose
    'parameter_inserted' field contains `needle`."""
    for entry in results_json.get("rigidity", []) or []:
        if needle in entry.get("parameter_inserted", ""):
            return entry
    return None


def input_by_name(inputs_json, name):
    """Look an input up BY NAME (never positionally) from a track's own
    inputs.json (a list of {name, value, tier, why} dicts)."""
    if isinstance(inputs_json, list):
        for x in inputs_json:
            if isinstance(x, dict) and x.get("name") == name:
                return x
    return None


def main():
    links = []

    # ---- raw loads (every value below is read, not typed) ----
    A_exports, A_exports_p = load(f"{A_DIR}/exports.json")
    A_results, A_results_p = load(f"{A_DIR}/results.json")
    A_inputs, A_inputs_p = load(f"{A_DIR}/inputs.json")
    A_m24, A_m24_p = load(f"{A_DIR}/m24_shapes.json")

    B_exports, B_exports_p = load(f"{B_DIR}/exports.json")
    B_part1, B_part1_p = load(f"{B_DIR}/part1_euler_results.json")
    B_inputs, B_inputs_p = load(f"{B_DIR}/inputs.json")

    C_exports, C_exports_p = load(f"{C_DIR}/exports.json")
    C_inputs, C_inputs_p = load(f"{C_DIR}/inputs.json")

    D_exports, D_exports_p = load(f"{D_DIR}/exports.json")
    D_inputs, D_inputs_p = load(f"{D_DIR}/inputs.json")

    E_poll, E_poll_p = load(f"{E_DIR}/results/00_track_d_poll.json")
    E_inputs, E_inputs_p = load(f"{E_DIR}/inputs.json")

    F_exports, F_exports_p = load(f"{F_DIR}/exports.json")
    F_results, F_results_p = load(f"{F_DIR}/results.json")
    F_inputs, F_inputs_p = load(f"{F_DIR}/inputs.json")

    d_chi = rat(D_exports["chi"])
    d_b2 = rat(D_exports["b2"])
    d_n_singular = rat(D_exports["n_singular"])

    # ================================================================
    # LINK 1: A.chi_from_genus (declared k) vs D.chi
    # ================================================================
    a_chi = rat(A_exports["chi_from_genus"])
    consistent1 = (a_chi == d_chi)
    pointer1 = find_file_pointer(A_DIR, "D-tda")
    shared1 = sorted(input_names(A_inputs) & input_names(D_inputs))
    k_entry = rigidity_entry(A_results, "k (overall factor")
    k_solution_ints = extract_ints(k_entry["solution_set"]) if k_entry else []
    k_non_unique = len(set(k_solution_ints)) > 1
    k_input = input_by_name(A_inputs, "k")
    links.append({
        "link": "A.chi_from_genus (Z(tau,0) at declared k) == D.chi",
        "A_value": str(a_chi), "A_source": A_exports_p,
        "D_value": str(d_chi), "D_source": D_exports_p,
        "consistent": consistent1,
        "file_pointer_A_to_D": pointer1,
        "method": "POINTER" if pointer1 else "CROSS_METHOD",
        "shared_declared_inputs": shared1,
        "independent_corroboration": consistent1 and not pointer1 and not k_non_unique,
        "caveat": (
            "A's own rigidity table classifies k as NORMALISATION with a "
            "non-unique admissible solution set (%s); k=%s was one admissible "
            "choice, informed by k=2 being the literature K3 elliptic-genus "
            "factor (see A/inputs.json 'k'.why: %r), not uniquely forced by "
            "Track A's internal tests. So agreement with D, though not a "
            "file-level POINTER, is not full independent corroboration of "
            "the numeral 24." % (
                k_entry["solution_set"] if k_entry else "?",
                k_input.get("value") if k_input else "?",
                k_input.get("why") if k_input else "?",
            )
            if k_non_unique else None
        ),
    })

    # ================================================================
    # LINK 2: B (DMVV product at z=0, y=1, q^0) vs D.chi, via Goettsche
    # ================================================================
    b_chi_own = rat(B_part1["chi_from_own_series = k*cB_sum"])
    b_goettsche_ok = str(B_part1["product_equals_goettsche"]).strip() == "True"
    trackD_source = B_part1.get("trackD_exports_source", {})
    b_chi_D_seen = rat(trackD_source.get("chi"))
    k_solved_equals_declared = str(B_part1.get("declared_k_equals_solved_k")).strip() == "True"
    consistent2 = (b_chi_own == d_chi) and b_goettsche_ok and (b_chi_D_seen == d_chi)
    pointer2 = find_file_pointer(B_DIR, "D-tda") or find_file_pointer(B_DIR, "trackD")
    shared2 = sorted(input_names(B_inputs) & input_names(D_inputs))
    links.append({
        "link": "B (DMVV product at z=0 == Goettsche with chi read from D) == D.chi",
        "B_own_series_chi": str(b_chi_own),
        "B_product_equals_goettsche": b_goettsche_ok,
        "B_chi_of_D_as_seen_by_B": str(b_chi_D_seen),
        "B_source": B_part1_p,
        "D_value": str(d_chi), "D_source": D_exports_p,
        "consistent": consistent2,
        "file_pointer_B_to_D": pointer2,
        "method": "POINTER",
        "shared_declared_inputs": shared2,
        "independent_corroboration": False,
        "caveat": (
            "B's declared factor k was reverse-solved from D's chi "
            "(k_solved_exactly_from_chi_D == declared k: %s, see "
            "part1_euler_results.json 'k_solved_exactly_from_chi_D' and "
            "'declared_k_equals_solved_k'), and B's own script "
            "(common.py:find_trackD_exports) opens D-tda/exports.json "
            "directly. So this link is a POINTER, not an independent "
            "cross-check of the numeral 24, even though the DMVV-vs-Goettsche "
            "IDENTITY at fixed k is itself a real, target-free computation "
            "(see B's own rigidity table, 'coefficients of G_2,...' entries)."
            % k_solved_equals_declared
        ),
    })

    # ================================================================
    # LINK 3: C (chi_top, b2 from Noether + Hodge index) vs D (chi, b2)
    # ================================================================
    c_chi_top = rat(C_exports["chi_top"]["value"])
    c_b2 = rat(C_exports["signature"]["b2"])
    consistent3a = (c_chi_top == d_chi)
    consistent3b = (c_b2 == d_b2)
    pointer3 = find_file_pointer(C_DIR, "D-tda")
    shared3 = sorted(input_names(C_inputs) & input_names(D_inputs))
    # Read each end's OWN admission of internal (non-)independence, rather
    # than asserting it: C's C1_chi_top result and D's exports.status.
    C_results, C_results_p = load(f"{C_DIR}/results.json")
    c1_note = next((r.get("note") for r in C_results.get("results", [])
                     if r.get("id") == "C1_chi_top"), None)
    d_status_note = D_exports.get("status")
    links.append({
        "link": "C.chi_top (Noether) == D.chi, and C.b2 (Hodge index) == D.b2",
        "C_chi_top": str(c_chi_top), "C_b2": str(c_b2), "C_source": C_exports_p,
        "D_chi": str(d_chi), "D_b2": str(d_b2), "D_source": D_exports_p,
        "consistent": consistent3a and consistent3b,
        "file_pointer_C_to_D": pointer3,
        "method": "POINTER" if pointer3 else "CROSS_METHOD",
        "shared_declared_inputs": shared3,
        "independent_corroboration": consistent3a and consistent3b and not pointer3 and not shared3,
        "caveat": (
            ("C's scripts reference D-tda directly: %s. " % pointer3 if pointer3 else "")
            + "Neither end is unconditional even though no shared declared "
              "input or file pointer was found between them: C's own note on "
              "C1_chi_top reads %r (a typed-constant x computed-number "
              "pattern, not a free derivation of 24), and D's own exports.json "
              "'status' field reads %r (a declared local-model input, tier L). "
              "C.b2=22 also follows arithmetically from C's own chi_top=24 "
              "given b0=b4=1, b1=b3=0 (not a second independent number)."
            % (c1_note, d_status_note)
        ),
    })

    # ================================================================
    # LINK 3b: C's LATTICE SELECTION (m,n) in the mU+n(-E8) family --
    # its signature must equal C's own derived Hodge-index signature, and
    # its rank must equal D.b2. Read from C-lattices/01_lattices.json's
    # structured fields (not from prose), per the task's own chain spec:
    # "C gives chi_top, then the signature, THEN THE LATTICE."
    # ================================================================
    C_lat, C_lat_p = load(f"{C_DIR}/01_lattices.json")
    selected = C_lat["selected"][0]
    c_sig = (rat(selected["signature_ldl"][0]), rat(selected["signature_ldl"][1]))
    c_derived_sig = tuple(rat(x) for x in C_exports["signature"]["value"])
    consistent3c = (c_sig == c_derived_sig)
    consistent3d = (rat(selected["rank"]) == d_b2)
    pointer3b = find_file_pointer(C_DIR, "D-tda")
    typed_gram_note = C_lat.get("typed_gram_3U_2mE8", {}).get("note")
    links.append({
        "link": "C's SELECTED lattice m*U+n(-E8) signature == C's own Hodge-index "
                "signature, and its rank == D.b2",
        "selected_m_n": [selected["m"], selected["n"]],
        "selected_lattice_signature": [str(c_sig[0]), str(c_sig[1])],
        "C_derived_hodge_signature": [str(c_derived_sig[0]), str(c_derived_sig[1])],
        "selected_lattice_rank": str(selected["rank"]), "D_b2": str(d_b2),
        "C_source": C_lat_p, "D_source": D_exports_p,
        "consistent": consistent3c and consistent3d,
        "file_pointer_C_to_D": pointer3b,
        "method": "POINTER",
        "shared_declared_inputs": ["k3_lattice_identification"],
        "independent_corroboration": False,
        "caveat": "the lattice FAMILY mU+n(-E8) is itself a declared input "
                  "(k3_lattice_identification, tier L, FROM MEMORY); given that "
                  "family, (m,n)=(3,2) is picked because it is the only one in "
                  "0<=m,n<=30 whose signature matches C's own Hodge-index result "
                  "-- a within-Track-C selection, not a second independent route "
                  "to chi=24. C's typed-Gram cross-check note reads %r."
                  % typed_gram_note,
    })

    # ================================================================
    # LINK 4a: E's imported n_singular == D.n_singular (POINTER by
    # construction: E polls D's exports.json file directly).
    # LINK 4b: E's own, separately-computed count (exact 2-torsion point
    # enumeration on (R/Z)^4) vs D.n_singular (CROSS_METHOD: a different,
    # elementary combinatorial computation, not read from D).
    # ================================================================
    e_imported = rat(E_poll["T4Z2_fixed_points_from_track_D"])
    e_own = rat(E_poll["T4Z2_fixed_points_computed_here_tierB"])
    consistent4a = (e_imported == d_n_singular)
    consistent4b = (e_own == d_n_singular)
    pointer4 = find_file_pointer(f"{E_DIR}/scripts", "D-tda")
    shared4 = sorted(input_names(E_inputs) & input_names(D_inputs))
    links.append({
        "link": "E.n_singular (imported from D) == D.n_singular",
        "E_value": str(e_imported), "E_source": E_poll_p,
        "D_value": str(d_n_singular), "D_source": D_exports_p,
        "consistent": consistent4a,
        "file_pointer_E_to_D": pointer4,
        "method": "POINTER",
        "shared_declared_inputs": shared4,
        "independent_corroboration": False,
        "caveat": "E reads this number directly from D-tda/exports.json (poll_track_d.py); "
                  "sha256 of the file it read is recorded in E's own results "
                  "(00_track_d_poll.json 'sha256_of_export_used').",
    })
    links.append({
        "link": "E's own combinatorial count of order-2 fixed points on (R/Z)^4 == D.n_singular",
        "E_value": str(e_own), "E_source": E_poll_p,
        "D_value": str(d_n_singular), "D_source": D_exports_p,
        "consistent": consistent4b,
        "file_pointer_E_to_D": None,
        "method": "CROSS_METHOD",
        "shared_declared_inputs": shared4,
        "independent_corroboration": consistent4b and not shared4,
        "caveat": "genuinely different method (exact enumeration of 2-torsion points of "
                  "(R/Z)^4, not GUDHI/Mayer-Vietoris), but the underlying fact "
                  "(2^4=16 fixed points of an order-2 involution on a 4-torus) is an "
                  "elementary counting identity, not a deep independent check.",
    })

    # ================================================================
    # LINK 5: F's M24 cycle shapes + group order vs A's own, separately
    # written, M24 cycle shapes + group order (two independent
    # constructions of the Golay code from QR mod 23 and its generators;
    # neither script imports the other -- checked below).
    # ================================================================
    f_shapes = sorted({s for lst in F_exports["M24_cycle_shapes_all_orders"].values() for s in lst})
    a_shapes = sorted({c["shape_str"] for c in A_m24["classes"]})
    f_order = rat(F_exports["M24_order"])
    a_order = rat(A_m24["group_order_computed"])
    consistent5 = (set(f_shapes) == set(a_shapes)) and (f_order == a_order)
    pointer5a = find_file_pointer(A_DIR, "F-whichk3")
    pointer5b = find_file_pointer(F_DIR, "A-genus")
    shared5 = sorted(input_names(A_inputs) & input_names(F_inputs))
    a_sample_size = A_m24.get("n_random_samples")
    f_sample_size = F_results["results"]["F4c_M24"]["computed"].get("sample_size")
    links.append({
        "link": "F.M24_cycle_shapes_all_orders (+order) == A.m24_shapes classes (+order)",
        "F_n_shapes": len(f_shapes), "F_order": str(f_order), "F_source": F_exports_p,
        "F_sample_size": f_sample_size,
        "A_n_shapes": len(a_shapes), "A_order": str(a_order), "A_source": A_m24_p,
        "A_sample_size": a_sample_size,
        "shapes_match": set(f_shapes) == set(a_shapes),
        "consistent": consistent5,
        "file_pointer_A_to_F": pointer5a, "file_pointer_F_to_A": pointer5b,
        "method": "POINTER" if (pointer5a or pointer5b) else "CROSS_METHOD",
        "shared_declared_inputs": shared5,
        "shared_unnamed_premise": "Golay code from QR mod 23, Aut(G24)=M24 (each FROM MEMORY, "
                                   "under differently named inputs: A='golay_and_generators', "
                                   "F='QR_golay_construction'/'Aut_G24_is_M24')",
        "independent_corroboration": consistent5 and not pointer5a and not pointer5b,
        "caveat": (
            "the two implementations were written independently (A-genus/m24_shapes.py, "
            "F-whichk3/f4_codes.py) and neither script's source references the other "
            "track's directory, so this IS independent corroboration of the group's "
            "structure -- but the cycle-SHAPE SETS were each built from a finite random "
            "sample (A: %s elements, F: %s elements), not an exhaustive class enumeration, "
            "so 'all shapes agree' means 'every shape either sampler happened to hit agrees', "
            "not 'the full conjugacy-class list was independently enumerated twice'. Both "
            "ends also share the unnamed literature premise above (see "
            "'shared_unnamed_premise')." % (a_sample_size, f_sample_size)
        ),
    })

    # ================================================================
    # LINK 6: F internal -- T(A)/T(Km A) discriminants consistent with
    # F's own binary-form enumeration, for the rank-2-T cases (the only
    # ones the enumeration covers; rank-4-T cases are out of its scope
    # and are reported separately, not counted as failures).
    # ================================================================
    f_ns_t = F_results["results"]["F2_NS_T"]["computed"]
    rank2_checks = [e for e in f_ns_t if e.get("rank_T") == 2]
    rank4_cases = [e for e in f_ns_t if e.get("rank_T") != 2]
    rank2_ok = all(e.get("T_in_binary_list") and e.get("TKm_in_binary_list") for e in rank2_checks)
    links.append({
        "link": "F: T(A)/T(Km A) discriminants (rank-2 T cases) found in F's own reduced "
                "binary-form enumeration (F2_binary_forms)",
        "n_rank2_T_cases_checked": len(rank2_checks),
        "all_rank2_T_and_TKm_in_binary_list": rank2_ok,
        "n_rank4_T_cases_out_of_scope": len(rank4_cases),
        "binary_form_enumeration": F_results["results"]["F2_binary_forms"]["computed"],
        "source": F_results_p,
        "consistent": rank2_ok,
        "method": "POINTER",  # both computed by the same script, f2_abelian.py
        "shared_declared_inputs": [],
        "independent_corroboration": False,
        "caveat": "internal self-consistency check within Track F (both quantities computed "
                  "by f2_abelian.py); rank-4-T cases (mixed CM points, e.g. tau=i, tau'=e^{i pi/3}) "
                  "are outside the scope of the rank-2 binary-form list by construction and "
                  "correctly report T_in_binary_list=False -- not a failure.",
    })

    all_consistent = all(l["consistent"] for l in links)
    independent_links = [l for l in links if l.get("independent_corroboration")]

    out = {
        "links": links,
        "all_consistent": all_consistent,
        "n_links": len(links),
        "n_independent_corroborations": len(independent_links),
        "independent_link_names": [l["link"] for l in independent_links],
        "command": "cd audit/k3t2_rigidity_v3/chain && "
                   "/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python "
                   "chain_check.py",
        "files_read": [
            A_exports_p, A_results_p, A_inputs_p, A_m24_p,
            B_exports_p, B_part1_p, B_inputs_p,
            C_exports_p, C_inputs_p, C_results_p, C_lat_p,
            D_exports_p, D_inputs_p,
            E_poll_p, E_inputs_p,
            F_exports_p, F_results_p, F_inputs_p,
        ],
    }

    out_path = HERE / "chain_check_results.json"
    with open(out_path, "w") as f:
        json.dump(out, f, indent=1)

    print(json.dumps(out, indent=1))
    print("\n=== SUMMARY ===")
    for l in links:
        status = "CONSISTENT" if l["consistent"] else "INCONSISTENT"
        indep = "INDEPENDENT" if l.get("independent_corroboration") else "NOT independent (" + l["method"] + ")"
        print(f"[{status}] [{indep}] {l['link']}")
    print(f"\nALL CONSISTENT: {all_consistent}")
    print(f"INDEPENDENT CROSS-METHOD CORROBORATIONS: {len(independent_links)} / {len(links)}")
    print(f"Results written to: {out_path.relative_to(REPO)}")


if __name__ == "__main__":
    main()
