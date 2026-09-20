"""POST-HOC replication and extension of the protein fold-separation result.

This is explicitly NOT the pre-registered test.  The pre-registered test
(run_proteins.py, commit c060fe9, expectations.json dde7a07) gave AUC = 1.000 with
permutation p = 0.0009, but on only 4 all-beta chains, which is its main weakness.
Everything here was decided AFTER seeing that result and is labelled post_hoc in
TopoDB so it can never be mistaken for the registered test:

  1. 19 further beta-rich entries, chosen from the literature fold literature by
     name only (barrels, immunoglobulin sandwiches, OB folds, propellers), fetched
     and classified by the SAME unchanged HELIX/SHEET rule before any persistence
     was computed on them.  AUC recomputed on the enlarged panel.
  2. A continuous test that does not depend on the 0.30/0.05 thresholds at all:
     Spearman rho between h1_small_bar_density and helix_frac across EVERY protein
     chain in the panel, with a permutation null, and the same for n_CA as the size
     confound control.

Seed 20260920.  Usage: python run_proteins_ext.py
"""
from __future__ import annotations

import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

sys.path.insert(0, "/mnt/disks/disk-socrateai-local-1/wt-topo-bio/topodb_runs/biology")
from bio_common import (BIO_DATA, RESULTS, Block, alpha_from_points, auc, betti_from_bars,  # noqa: E402
                        dominance, finite, rank_p, sha256, top_bars)
from run_proteins import (LARGE_BAND_A, SMALL_BAND_A, fold_class, h1_stats,  # noqa: E402
                          parse_pdb, ss_fractions)

SEED = 20260920
N_PERM = 10000
SCRIPT = "topodb_runs/biology/run_proteins_ext.py"
COMMAND = ("timeout 590 prlimit --as=8589934592 -- "
           "/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python "
           "topodb_runs/biology/run_proteins_ext.py")
URL = "https://files.rcsb.org/download/{}.pdb"

# beta-rich entries chosen by literature fold name only, before any computation here
EXTRA_BETA = {
    "1TIT": "titin I27, immunoglobulin beta sandwich",
    "1CSP": "cold-shock protein B, OB fold, all beta",
    "1MJC": "major cold-shock protein CspA, OB fold",
    "1WIT": "twitchin Ig domain, beta sandwich",
    "1FNA": "fibronectin type III 10th domain, beta sandwich",
    "1TTG": "fibronectin type III, beta sandwich",
    "1BEB": "beta-lactoglobulin, lipocalin beta barrel",
    "2SOD": "Cu,Zn superoxide dismutase, Greek-key beta barrel",
    "1REI": "immunoglobulin VL domain (Bence-Jones), beta sandwich",
    "1HOE": "tendamistat alpha-amylase inhibitor, all beta",
    "7AHL": "staphylococcal alpha-hemolysin, beta-barrel pore",
    "2OMF": "OmpF porin, 16-stranded beta barrel",
    "1BXW": "OmpA transmembrane domain, 8-stranded beta barrel",
    "1PRN": "porin from Rhodobacter, beta barrel",
    "1TBG": "G protein beta subunit, 7-bladed beta propeller",
    "1PIN": "Pin1 WW domain, three-stranded beta sheet",
    "1E0L": "FBP28 WW domain, beta sheet",
    "1CD8": "CD8 alpha, immunoglobulin beta sandwich",
    "1AF6": "maltoporin LamB, 18-stranded beta barrel",
}


def fetch(pid):
    out = BIO_DATA / "pdb" / f"{pid}.pdb"
    url = URL.format(pid)
    if out.exists() and out.stat().st_size > 1000:
        return {"id": pid, "kind": "pdb", "url": url, "path": str(out),
                "bytes": out.stat().st_size, "sha256": sha256(out), "cached": True}, None
    tried = []
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=60) as fh:
                data = fh.read()
            if len(data) < 1000:
                raise ValueError(f"short: {len(data)} bytes")
            out.write_bytes(data)
            return {"id": pid, "kind": "pdb", "url": url, "path": str(out), "bytes": len(data),
                    "sha256": sha256(out), "cached": False,
                    "fetched_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}, None
        except (urllib.error.URLError, urllib.error.HTTPError, ValueError, OSError) as exc:
            tried.append(f"attempt {attempt + 1}: {type(exc).__name__}: {exc}")
            time.sleep(2)
    return None, {"id": pid, "urls_tried": [url] * 3, "errors": tried, "status": "ABSENT"}


def perm_auc_p(vals, labels, n_perm=N_PERM, seed=SEED):
    obs = auc(vals[labels == 1], vals[labels == 0])
    rng = np.random.default_rng(seed)
    null = np.empty(n_perm)
    for i in range(n_perm):
        lab = rng.permutation(labels)
        null[i] = auc(vals[lab == 1], vals[lab == 0])
    return float(obs), float(rank_p(abs(obs - 0.5), np.abs(null - 0.5))), null


def main() -> int:
    blk = Block("proteins_ext", SCRIPT, COMMAND)
    base = json.loads((RESULTS / "proteins_table.json").read_text())
    table = [dict(r, panel="registered") for r in base["table"]]

    absent = []
    new_rows = []
    for pid in EXTRA_BETA:
        ent, miss = fetch(pid)
        if miss:
            absent.append(miss)
            print(f"  ABSENT {pid}")
            continue
        parsed = parse_pdb(Path(ent["path"]), " CA ")
        if parsed is None or len(parsed[0]) < 12:
            print(f"  skip {pid}: too few CA")
            continue
        pts, ch, helix, sheet, resseqs, n_chains = parsed
        n = len(pts)
        hf, sf, annotated = ss_fractions(resseqs, helix, sheet)
        fc = fold_class(hf, sf, annotated)
        ds_id = f"biology/rcsb_{pid.lower()}"
        blk.dataset(id=ds_id, domain="biology",
                    title=f"RCSB {pid} chain {ch}: {n} C-alpha atoms (post-hoc beta extension)",
                    source=ent["url"], provenance="experiment", n_objects=n, ambient_dim=3,
                    units="angstrom", sha256=ent["sha256"], local_path=ent["path"],
                    notes=(f"POST-HOC extension panel, chosen after the registered test to add "
                           f"beta-rich chains ({EXTRA_BETA[pid]}); chain {ch} of {n_chains}; "
                           f"helix_frac={hf:.3f} sheet_frac={sf:.3f} fold_class={fc} by the same "
                           f"unchanged HELIX/SHEET rule"))
        t = time.time()
        dga = alpha_from_points(pts, max_hom_dim=2)
        wa = round(time.time() - t, 2)
        sa = h1_stats(dga[1], n)
        blk.run(dataset_id=ds_id, method="alpha", coeff_field=2, max_dim=2,
                params={"max_hom_dim": 2, "metric": "euclidean_coordinates_angstrom",
                        "filtration": "plain radius = sqrt(alpha^2)", "atom": "CA", "chain": ch,
                        "small_band_death_lt_A": SMALL_BAND_A, "large_band_death_ge_A": LARGE_BAND_A,
                        "n_points": n, "panel": "post_hoc_beta_extension"},
                preprocessing=f"first model, altLoc blank/A, one atom per residue, chain {ch}",
                seed=SEED, tier="X", wall_sec=wa, diagrams=dga,
                betti={0: betti_from_bars(dga[0]), 1: betti_from_bars(dga[1]),
                       2: betti_from_bars(dga[2])}, expected_betti={0: 1, 1: 0, 2: 0},
                stats=([{"name": k, "value": float(v)} for k, v in sa.items()]
                       + [{"name": "n_ca_atoms", "value": float(n)},
                          {"name": "helix_frac", "value": float(hf)},
                          {"name": "sheet_frac", "value": float(sf)}]),
                controls=[{"kind": "known_answer",
                           "description": "alpha complex of a connected chain must end with b0 = 1",
                           "passed": bool(betti_from_bars(dga[0]) == 1),
                           "detail": f"b0={betti_from_bars(dga[0])}"}],
                findings=[{"claim": f"{pid} ({fc}, post-hoc extension): h1_small_bar_density "
                                    f"{sa['h1_small_bar_density']:.3f}, {sa['h1_bar_count']} finite "
                                    f"H1 bars over {n} C-alpha atoms.",
                           "verdict": "recovered", "tier": "X",
                           "caveat": "post-hoc extension panel, not part of the pre-registered test"}])
        row = {"id": pid, "kind": "pdb", "chain": ch, "n": n, "helix_frac": hf, "sheet_frac": sf,
               "fold_class": fc, "annotated": annotated, "sha256": ent["sha256"], "url": ent["url"],
               "panel": "post_hoc_extension",
               "alpha": {k: float(v) for k, v in sa.items()},
               "alpha_top_h1_bars": [[round(b, 4), round(d, 4)] for b, d in top_bars(dga[1], 5)]}
        new_rows.append(row)
        print(f"  {pid:5s} {fc:12s} n={n:5d} h={hf:.2f} s={sf:.2f} "
              f"dens={sa['h1_small_bar_density']:.3f}")
    table += new_rows

    prot = [r for r in table if r["fold_class"] != "rna"]
    A = [r for r in prot if r["fold_class"] == "all_alpha"]
    B = [r for r in prot if r["fold_class"] == "all_beta"]
    labels = np.array([1] * len(A) + [0] * len(B))
    out = {"n_all_alpha": len(A), "n_all_beta": len(B),
           "all_alpha_ids": [r["id"] for r in A], "all_beta_ids": [r["id"] for r in B],
           "n_beta_added_by_extension": sum(1 for r in new_rows if r["fold_class"] == "all_beta"),
           "absent": absent, "n_perm": N_PERM, "seed": SEED}

    prim_vals = np.array([r["alpha"]["h1_small_bar_density"] for r in A + B])
    a_prim, p_prim, _ = perm_auc_p(prim_vals, labels)
    # Does the pre-registered 4.5 A band do any work?  The same test on the UNBANDED
    # count of H1 bars per residue answers that directly.  On the registered panel the
    # two gave an identical AUC of 1.000, which says the band was decorative; this
    # records the comparison on the enlarged panel too, so the honest reading of the
    # primary statistic is in the database rather than only in a commit message.
    unb_vals = np.array([r["alpha"]["h1_bar_count_per_residue"] for r in A + B])
    a_unb, p_unb, _ = perm_auc_p(unb_vals, labels)
    size_vals = np.array([float(r["n"]) for r in A + B])
    a_size, p_size, _ = perm_auc_p(size_vals, labels)
    out["enlarged_panel"] = {
        "auc_h1_small_bar_density": round(a_prim, 4), "p_perm": round(p_prim, 5),
        "auc_size_only_n_CA": round(a_size, 4), "p_perm_size": round(p_size, 5),
        "beats_size_baseline": bool(abs(a_prim - 0.5) > abs(a_size - 0.5)),
        "auc_unbanded_h1_bar_count_per_residue": round(a_unb, 4), "p_perm_unbanded": round(p_unb, 5),
        "band_adds_nothing": bool(abs(a_unb - a_prim) < 0.02),
        "alpha_values": sorted(round(r["alpha"]["h1_small_bar_density"], 3) for r in A),
        "beta_values": sorted(round(r["alpha"]["h1_small_bar_density"], 3) for r in B)}

    # ---- continuous test over ALL protein chains, threshold-free ----
    hv = np.array([r["helix_frac"] for r in prot])
    sv = np.array([r["sheet_frac"] for r in prot])
    dv = np.array([r["alpha"]["h1_small_bar_density"] for r in prot])
    nv = np.array([float(r["n"]) for r in prot])
    cont = {}
    rng = np.random.default_rng(SEED)
    for name, x in (("helix_frac", hv), ("sheet_frac", sv), ("n_CA_size_control", nv)):
        rho = float(spearmanr(x, dv).statistic)
        null = np.array([float(spearmanr(rng.permutation(x), dv).statistic) for _ in range(2000)])
        cont[name] = {"spearman_rho_vs_h1_small_bar_density": round(rho, 4),
                      "p_perm": round(float(rank_p(abs(rho), np.abs(null))), 5),
                      "n_perm": 2000, "n_chains": int(len(prot))}
    out["continuous_test_all_protein_chains"] = cont

    blk.dataset(id="biology/rcsb_fold_panel_extended", domain="biology",
                title=f"Extended RCSB fold panel: {len(prot)} protein chains "
                      f"({len(A)} all-alpha, {len(B)} all-beta)",
                source="https://files.rcsb.org/download/ (per-entry sha256 in pdb_manifest.json "
                       "and proteins_ext_manifest.json)",
                provenance="experiment", n_objects=len(prot), ambient_dim=3, units="angstrom",
                notes="POST-HOC panel: the registered 39-chain panel plus 19 beta-rich entries "
                      "added after the registered result was known. Every analysis on this row "
                      "is post hoc and is labelled so.")
    blk.run(dataset_id="biology/rcsb_fold_panel_extended", method="alpha", coeff_field=2, max_dim=2,
            params={"test": "POST-HOC replication of the fold separation on an enlarged panel, "
                            "plus a threshold-free continuous test",
                    "primary_statistic": "h1_small_bar_density", "n_perm": N_PERM,
                    "small_band_death_lt_A": SMALL_BAND_A,
                    "metric": "euclidean_coordinates_angstrom",
                    "prespecified": False,
                    "registered_result_being_replicated": "AUC 1.000, p 0.0009, 14 vs 4 chains "
                                                          "(run_proteins.py, commit c060fe9)"},
            preprocessing="alpha diagrams of every protein chain in the extended panel",
            seed=SEED, tier="X",
            stats=[{"name": "auc_h1_small_bar_density_extended", "value": round(a_prim, 4),
                    "null_model": f"{N_PERM} permutations of the all-alpha/all-beta labels",
                    "n_null": N_PERM, "p_value": round(p_prim, 5), "p_method": "rank",
                    "multiplicity": "post hoc; not corrected, and not pre-registered"},
                   {"name": "auc_unbanded_h1_bar_count_per_residue_extended", "value": round(a_unb, 4),
                    "null_model": f"{N_PERM} permutations of the same labels", "n_null": N_PERM,
                    "p_value": round(p_unb, 5), "p_method": "rank",
                    "multiplicity": "post hoc; reported to show whether the pre-registered 4.5 A "
                                    "band contributes anything over an unbanded bar count"},
                   {"name": "auc_size_only_n_CA_extended", "value": round(a_size, 4),
                    "null_model": f"{N_PERM} permutations of the same labels", "n_null": N_PERM,
                    "p_value": round(p_size, 5), "p_method": "rank"},
                   {"name": "spearman_helix_frac_vs_h1_small_bar_density",
                    "value": cont["helix_frac"]["spearman_rho_vs_h1_small_bar_density"],
                    "null_model": "2000 permutations of helix_frac against the same densities",
                    "n_null": 2000, "p_value": cont["helix_frac"]["p_perm"], "p_method": "rank"},
                   {"name": "spearman_sheet_frac_vs_h1_small_bar_density",
                    "value": cont["sheet_frac"]["spearman_rho_vs_h1_small_bar_density"],
                    "null_model": "2000 permutations of sheet_frac against the same densities",
                    "n_null": 2000, "p_value": cont["sheet_frac"]["p_perm"], "p_method": "rank"},
                   {"name": "spearman_n_CA_vs_h1_small_bar_density",
                    "value": cont["n_CA_size_control"]["spearman_rho_vs_h1_small_bar_density"],
                    "null_model": "2000 permutations of n_CA against the same densities",
                    "n_null": 2000, "p_value": cont["n_CA_size_control"]["p_perm"],
                    "p_method": "rank"}],
            controls=[{"kind": "negative", "description":
                       "size-only confound control on the enlarged panel: AUC from n_CA alone",
                       "passed": out["enlarged_panel"]["beats_size_baseline"],
                       "detail": f"topology AUC {a_prim:.3f} vs size AUC {a_size:.3f}"},
                      {"kind": "injection", "description":
                       "band-necessity control: the same test on the UNBANDED count of H1 bars per "
                       "residue. If it matches the banded AUC, the pre-registered 4.5 A cut-off is "
                       "doing no work and the statistic is simply 'H1 bars per residue'.",
                       "passed": bool(abs(a_unb - a_prim) < 0.02),
                       "detail": f"banded AUC {a_prim:.3f} vs unbanded {a_unb:.3f} "
                                 f"(p {p_unb:.4f}); on the registered panel both were 1.000"},
                      {"kind": "shuffle", "description":
                       "threshold-free control: the same association measured as a Spearman "
                       "correlation across all protein chains, so it does not depend on the "
                       "0.30/0.05 fold-class thresholds",
                       "passed": bool(cont["helix_frac"]["p_perm"] <= 0.01),
                       "detail": f"rho(helix_frac) = "
                                 f"{cont['helix_frac']['spearman_rho_vs_h1_small_bar_density']}, "
                                 f"p {cont['helix_frac']['p_perm']}; rho(n_CA) = "
                                 f"{cont['n_CA_size_control']['spearman_rho_vs_h1_small_bar_density']}, "
                                 f"p {cont['n_CA_size_control']['p_perm']}"}],
            findings=[{"claim": f"POST-HOC: on an enlarged panel of {len(A)} all-alpha vs {len(B)} "
                                f"all-beta chains, h1_small_bar_density gives AUC {a_prim:.3f} "
                                f"(permutation p {p_prim:.4f}) against a size-only AUC of "
                                f"{a_size:.3f}. Threshold-free, across all {len(prot)} protein "
                                f"chains, h1_small_bar_density correlates with helix_frac at "
                                f"Spearman rho "
                                f"{cont['helix_frac']['spearman_rho_vs_h1_small_bar_density']} "
                                f"(p {cont['helix_frac']['p_perm']}) and with chain length at rho "
                                f"{cont['n_CA_size_control']['spearman_rho_vs_h1_small_bar_density']} "
                                f"(p {cont['n_CA_size_control']['p_perm']}).",
                       "verdict": "recovered" if (abs(a_prim - 0.5) > abs(a_size - 0.5)
                                                  and p_prim <= 0.01) else "inconclusive",
                       "tier": "X",
                       "caveat": "POST HOC. The panel was enlarged after the registered result was "
                                 "known. This replicates, it does not independently confirm. The "
                                 "fold classes are this pipeline's HELIX/SHEET counts, not CATH or "
                                 "SCOP.",
                       "reference": "topodb_runs/biology/run_proteins_ext.py"}])

    (BIO_DATA / "proteins_ext_manifest.json").write_text(json.dumps(
        {"extra_beta_literature_labels": EXTRA_BETA, "absent": absent,
         "entries": [{"id": r["id"], "url": r["url"], "sha256": r["sha256"]} for r in new_rows]},
        indent=1))
    (RESULTS / "proteins_ext_table.json").write_text(json.dumps(
        {"seed": SEED, "new_rows": new_rows, "analysis": out}, indent=1))
    print(json.dumps(out["enlarged_panel"], indent=1))
    print(json.dumps(cont, indent=1))
    blk.write()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
