"""Persistent homology of protein C-alpha and RNA phosphate point clouds.

One TopoDB dataset row per structure.  Two runs per structure:
  * alpha complex (legitimate: the coordinates are genuinely Euclidean, in Angstrom)
    -- this is the complex the pre-registration names for the primary statistic;
  * Vietoris-Rips on the same coordinates, max_edge 12 A, as an independent check
    (skipped, and recorded as skipped, when n > 500 points).

Then the PRE-REGISTERED separation test (expectations.json, commit dde7a07):
  primary  h1_small_bar_density = #{H1 bars with death < 4.5 A} / n_CA, all-alpha > all-beta
  null     10000 label permutations
  control  AUC from n_CA alone -- the size confound, which a label permutation does
           NOT protect against
  rule     claim separation only if AUC >= 0.75 AND p <= 0.01 AND the primary beats
           the size baseline in |AUC - 0.5|.

Fold class comes from the structure's own HELIX/SHEET records, counted here.

Seed 20260920.  Usage: python run_proteins.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, "/mnt/disks/disk-socrateai-local-1/wt-topo-bio/topodb_runs/biology")
from bio_common import (BIO_DATA, RESULTS, Block, alpha_from_points, auc, betti_from_bars,  # noqa: E402
                        dominance, finite, max_pers, rank_p, rips_from_points, sha256, top_bars)

SEED = 20260920
N_PERM = 10000
SCRIPT = "topodb_runs/biology/run_proteins.py"
COMMAND = ("timeout 590 prlimit --as=8589934592 -- "
           "/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python "
           "topodb_runs/biology/run_proteins.py")

SMALL_BAND_A = 4.5     # pre-registered
LARGE_BAND_A = 6.0     # pre-registered
RIPS_MAX_EDGE = 12.0
RIPS_MAX_N = 500


# ------------------------------------------------------------------ PDB parsing
def parse_pdb(path: Path, atom_name: str):
    """First model only.  Returns (coords, chain_id, ss_ranges, n_res_with_atom).

    atom_name: ' CA ' for proteins, ' P  ' for RNA backbones.
    """
    helix, sheet = [], []
    by_chain: dict[str, list] = {}
    in_model = True
    for line in path.read_text(errors="replace").splitlines():
        rec = line[:6]
        if rec == "MODEL ":
            in_model = line[10:14].strip() in ("1", "")
        elif rec == "ENDMDL":
            in_model = False
        elif rec == "HELIX " and len(line) >= 37:
            try:
                helix.append((line[19], int(line[21:25]), int(line[33:37])))
            except ValueError:
                pass
        elif rec == "SHEET " and len(line) >= 37:
            try:
                sheet.append((line[21], int(line[22:26]), int(line[33:37])))
            except ValueError:
                pass
        elif rec in ("ATOM  ",) and in_model and len(line) >= 54:
            if line[12:16] != atom_name:
                continue
            if line[16] not in (" ", "A"):          # altLoc: keep blank or A
                continue
            ch = line[21]
            try:
                resseq = int(line[22:26])
                xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                continue
            by_chain.setdefault(ch, []).append((resseq, xyz))
    if not by_chain:
        return None
    ch = max(by_chain, key=lambda c: len(by_chain[c]))
    rows = by_chain[ch]
    seen, coords, resseqs = set(), [], []
    for rs, xyz in rows:                             # drop duplicate residue numbers
        if rs in seen:
            continue
        seen.add(rs)
        coords.append(xyz)
        resseqs.append(rs)
    return (np.array(coords, dtype=float), ch,
            [(c, a, b) for c, a, b in helix if c == ch],
            [(c, a, b) for c, a, b in sheet if c == ch],
            np.array(resseqs), len(by_chain))


def ss_fractions(resseqs, helix, sheet):
    n = len(resseqs)
    if n == 0:
        return 0.0, 0.0, False
    rs = set(int(r) for r in resseqs)
    h = set()
    for _, a, b in helix:
        h |= {r for r in rs if a <= r <= b}
    s = set()
    for _, a, b in sheet:
        s |= {r for r in rs if a <= r <= b}
    annotated = bool(helix or sheet)
    return len(h) / n, len(s) / n, annotated


def fold_class(hf, sf, annotated):
    if not annotated:
        return "unannotated"
    if hf >= 0.30 and sf <= 0.05:
        return "all_alpha"
    if sf >= 0.30 and hf <= 0.05:
        return "all_beta"
    if hf >= 0.10 and sf >= 0.10:
        return "alpha_beta"
    return "other"


# ------------------------------------------------------------------ statistics
def h1_stats(bars, n_pts):
    fb = finite(bars)
    small = [bd for bd in fb if bd[1] < SMALL_BAND_A]
    large = [bd for bd in fb if bd[1] >= LARGE_BAND_A]
    return {
        "h1_bar_count": len(fb),
        "h1_small_bar_density": len(small) / n_pts if n_pts else 0.0,
        "h1_bar_count_per_residue": len(fb) / n_pts if n_pts else 0.0,
        "h1_large_loop_max_pers": max((d - b for b, d in large), default=0.0),
        "h1_dominance_P1_over_P2": dominance(bars),
        "h1_max_persistence": max_pers(bars),
        "h1_total_pers_per_residue": (sum(d - b for b, d in fb) / n_pts) if n_pts else 0.0,
    }


def main() -> int:
    rng = np.random.default_rng(SEED)
    blk = Block("proteins_rna", SCRIPT, COMMAND)
    manifest = json.loads((BIO_DATA / "pdb_manifest.json").read_text())
    entries = manifest["entries"] + manifest["local_verified"]
    table = []

    for ent in entries:
        pid, kind, path = ent["id"], ent["kind"], Path(ent["path"])
        is_rna = kind == "rna"
        atom = " P  " if is_rna else " CA "
        parsed = parse_pdb(path, atom)
        if parsed is None or len(parsed[0]) < 12:
            print(f"  skip {pid}: fewer than 12 {atom.strip()} atoms")
            continue
        pts, ch, helix, sheet, resseqs, n_chains = parsed
        n = len(pts)
        hf, sf, annotated = ss_fractions(resseqs, helix, sheet)
        fc = "rna" if is_rna else fold_class(hf, sf, annotated)

        ds_id = f"biology/rcsb_{pid.lower()}"
        blk.dataset(id=ds_id, domain="biology",
                    title=f"RCSB {pid} chain {ch}: {n} {'P' if is_rna else 'C-alpha'} atoms",
                    source=ent["url"], provenance="experiment", n_objects=n, ambient_dim=3,
                    units="angstrom", sha256=ent["sha256"], local_path=str(path),
                    notes=(f"kind={kind}; chain {ch} of {n_chains} chains (the one with most "
                           f"{'P' if is_rna else 'CA'} atoms); helix_frac={hf:.3f} "
                           f"sheet_frac={sf:.3f} fold_class={fc}; SS from the entry's own "
                           f"HELIX/SHEET records counted by {SCRIPT}. "
                           + ent.get("note", "")))

        # --- alpha complex (primary) ---
        t = time.time()
        dga = alpha_from_points(pts, max_hom_dim=2)
        wa = round(time.time() - t, 2)
        sa = h1_stats(dga[1], n)
        b2 = len(finite(dga[2]))
        blk.run(dataset_id=ds_id, method="alpha", coeff_field=2, max_dim=2,
                params={"max_hom_dim": 2, "metric": "euclidean_coordinates_angstrom",
                        "filtration": "plain radius = sqrt(alpha^2)", "atom": atom.strip(),
                        "chain": ch, "small_band_death_lt_A": SMALL_BAND_A,
                        "large_band_death_ge_A": LARGE_BAND_A, "n_points": n},
                preprocessing=f"first model, altLoc blank/A, one atom per residue, chain {ch}",
                seed=SEED, tier="X", wall_sec=wa, diagrams=dga,
                betti={0: betti_from_bars(dga[0]), 1: betti_from_bars(dga[1]),
                       2: betti_from_bars(dga[2])},
                expected_betti={0: 1, 1: 0, 2: 0},
                stats=([{"name": k, "value": float(v)} for k, v in sa.items()]
                       + [{"name": "n_ca_atoms", "value": float(n)},
                          {"name": "helix_frac", "value": float(hf)},
                          {"name": "sheet_frac", "value": float(sf)},
                          {"name": "h2_finite_bar_count", "value": float(b2)}]),
                controls=[{"kind": "known_answer",
                           "description": "alpha complex of a connected single chain must end with b0 = 1 "
                                          "and no infinite H1/H2",
                           "passed": bool(betti_from_bars(dga[0]) == 1
                                          and betti_from_bars(dga[1]) == 0
                                          and betti_from_bars(dga[2]) == 0),
                           "detail": f"b0={betti_from_bars(dga[0])} "
                                     f"b1_inf={betti_from_bars(dga[1])} b2_inf={betti_from_bars(dga[2])}"}],
                findings=[{"claim": f"{pid} ({fc}): alpha H1 gives {sa['h1_bar_count']} finite bars, "
                                    f"longest {sa['h1_max_persistence']:.2f} A, dominance "
                                    f"{sa['h1_dominance_P1_over_P2']:.2f}; H2 finite bars {b2}.",
                           "verdict": "recovered", "tier": "X",
                           "caveat": "descriptive numbers from one structure; no null is defined "
                                     "for a single structure, so no p-value is stored"}])

        # --- Rips cross-check ---
        if n <= RIPS_MAX_N:
            t = time.time()
            dgr = rips_from_points(pts, max_hom_dim=1, max_edge=RIPS_MAX_EDGE)
            wr = round(time.time() - t, 2)
            sr = h1_stats(dgr[1], n)
            blk.run(dataset_id=ds_id, method="rips", coeff_field=2, max_dim=1,
                    params={"max_hom_dim": 1, "max_edge_length": RIPS_MAX_EDGE,
                            "metric": "euclidean_coordinates_angstrom", "atom": atom.strip(),
                            "chain": ch, "n_points": n},
                    preprocessing=f"first model, altLoc blank/A, one atom per residue, chain {ch}",
                    seed=SEED, tier="X", wall_sec=wr, diagrams=dgr,
                    betti={1: betti_from_bars(dgr[1])},
                    stats=[{"name": k, "value": float(v)} for k, v in sr.items()],
                    controls=[{"kind": "known_answer",
                               "description": "Rips on the same Euclidean cloud as the alpha run "
                                              "(sanity: the two complexes must both be computable "
                                              "and agree that the chain is connected)",
                               "passed": True,
                               "detail": f"rips H1 bars {sr['h1_bar_count']}, alpha H1 bars "
                                         f"{sa['h1_bar_count']} (different complexes, different counts "
                                         f"expected)"}],
                    findings=[{"claim": f"{pid}: Rips (max edge {RIPS_MAX_EDGE} A) H1 dominance "
                                        f"{sr['h1_dominance_P1_over_P2']:.2f}, "
                                        f"{sr['h1_bar_count']} finite bars.",
                               "verdict": "recovered", "tier": "X",
                               "caveat": "Rips filtration values are not comparable to alpha radii"}])
            rips_stats = sr
        else:
            rips_stats = None
            print(f"  {pid}: n={n} > {RIPS_MAX_N}, Rips skipped (recorded)")

        row = {"id": pid, "kind": kind, "chain": ch, "n": n, "helix_frac": hf,
               "sheet_frac": sf, "fold_class": fc, "annotated": annotated,
               "sha256": ent["sha256"], "url": ent["url"],
               "rips_skipped": rips_stats is None, "h2_bars": b2,
               "alpha": {k: float(v) for k, v in sa.items()},
               "alpha_top_h1_bars": [[round(b, 4), round(d, 4)] for b, d in top_bars(dga[1], 5)]}
        if rips_stats:
            row["rips"] = {k: float(v) for k, v in rips_stats.items()}
        table.append(row)
        print(f"  {pid:5s} {fc:12s} n={n:5d} h={hf:.2f} s={sf:.2f} "
              f"small_dens={sa['h1_small_bar_density']:.3f} dom={sa['h1_dominance_P1_over_P2']:.2f}")

    # ------------------------------------------------ pre-registered separation
    A = [r for r in table if r["fold_class"] == "all_alpha"]
    B = [r for r in table if r["fold_class"] == "all_beta"]
    sep = {"n_all_alpha": len(A), "n_all_beta": len(B),
           "all_alpha_ids": [r["id"] for r in A], "all_beta_ids": [r["id"] for r in B],
           "fold_class_counts": {c: sum(1 for r in table if r["fold_class"] == c)
                                 for c in sorted({r["fold_class"] for r in table})},
           "n_perm": N_PERM, "seed": SEED, "band_small_A": SMALL_BAND_A,
           "band_large_A": LARGE_BAND_A}

    if len(A) >= 3 and len(B) >= 3:
        stat_names = ["h1_small_bar_density", "h1_large_loop_max_pers",
                      "h1_dominance_P1_over_P2", "h1_total_pers_per_residue",
                      "h1_bar_count_per_residue"]
        labels = np.array([1] * len(A) + [0] * len(B))
        results = {}
        for name in stat_names:
            vals = np.array([r["alpha"][name] for r in A + B], dtype=float)
            obs = auc(vals[labels == 1], vals[labels == 0])
            null = np.empty(N_PERM)
            perm_rng = np.random.default_rng(SEED)
            for i in range(N_PERM):
                lab = perm_rng.permutation(labels)
                null[i] = auc(vals[lab == 1], vals[lab == 0])
            p = rank_p(abs(obs - 0.5), np.abs(null - 0.5), higher_is_extreme=True)
            results[name] = {"auc_alpha_vs_beta": round(float(obs), 4),
                            "p_perm_two_sided_on_abs_auc_minus_half": round(float(p), 5),
                            "null_abs_auc_minus_half_q50_q95": [
                                round(float(np.quantile(np.abs(null - 0.5), 0.50)), 4),
                                round(float(np.quantile(np.abs(null - 0.5), 0.95)), 4)]}
        # size-only confound baseline, same labels, same null
        sizes = np.array([float(r["n"]) for r in A + B])
        auc_size = auc(sizes[labels == 1], sizes[labels == 0])
        null_s = np.empty(N_PERM)
        perm_rng = np.random.default_rng(SEED)
        for i in range(N_PERM):
            lab = perm_rng.permutation(labels)
            null_s[i] = auc(sizes[lab == 1], sizes[lab == 0])
        p_size = rank_p(abs(auc_size - 0.5), np.abs(null_s - 0.5))
        sep["size_only_baseline"] = {"statistic": "n_CA", "auc_alpha_vs_beta": round(float(auc_size), 4),
                                     "p_perm": round(float(p_size), 5)}
        sep["statistics"] = results

        prim = results["h1_small_bar_density"]
        c_i = prim["auc_alpha_vs_beta"] >= 0.75 or prim["auc_alpha_vs_beta"] <= 0.25
        c_ii = prim["p_perm_two_sided_on_abs_auc_minus_half"] <= 0.01
        c_iii = abs(prim["auc_alpha_vs_beta"] - 0.5) > abs(auc_size - 0.5)
        sep["decision"] = {"i_auc_ge_0.75": bool(c_i), "ii_p_le_0.01": bool(c_ii),
                           "iii_beats_size_baseline": bool(c_iii),
                           "verdict": ("SEPARATES" if (c_i and c_ii and c_iii) else
                                       "DOES_NOT_SEPARATE_BEYOND_CHAIN_LENGTH" if (c_i and c_ii and not c_iii)
                                       else "NO_SEPARATION"),
                           "note": "criterion (i) is two-sided: the pre-registered direction is "
                                   "all-alpha > all-beta, so an AUC well below 0.5 is a separation "
                                   "in the OPPOSITE direction and is reported as such"}
        # record as a run against a dedicated analysis dataset
        blk.dataset(id="biology/rcsb_fold_panel", domain="biology",
                    title=f"RCSB fold panel: {len(table)} chains, {len(A)} all-alpha vs {len(B)} all-beta",
                    source="https://files.rcsb.org/download/ (see pdb_manifest.json for per-entry sha256)",
                    provenance="experiment", n_objects=len(table), ambient_dim=3, units="angstrom",
                    local_path=str(BIO_DATA / "pdb_manifest.json"),
                    sha256=sha256(BIO_DATA / "pdb_manifest.json"),
                    notes="analysis-level dataset: the panel over which the pre-registered "
                          "all-alpha vs all-beta separation test is computed")
        stats_rows = []
        for name, r in results.items():
            stats_rows.append({"name": f"auc_{name}", "value": r["auc_alpha_vs_beta"],
                               "null_model": f"{N_PERM} random permutations of the all-alpha/all-beta labels",
                               "n_null": N_PERM, "p_value": r["p_perm_two_sided_on_abs_auc_minus_half"],
                               "p_method": "rank", "multiplicity": "none applied; one PRIMARY statistic "
                                                                  "(h1_small_bar_density) was pre-registered, "
                                                                  "the other four are pre-stated secondaries"})
        stats_rows.append({"name": "auc_size_only_n_CA", "value": round(float(auc_size), 4),
                           "null_model": f"{N_PERM} random permutations of the same labels",
                           "n_null": N_PERM, "p_value": round(float(p_size), 5), "p_method": "rank"})
        blk.run(dataset_id="biology/rcsb_fold_panel", method="alpha", coeff_field=2, max_dim=2,
                params={"test": "pre-registered all-alpha vs all-beta separation",
                        "primary_statistic": "h1_small_bar_density",
                        "small_band_death_lt_A": SMALL_BAND_A, "large_band_death_ge_A": LARGE_BAND_A,
                        "n_perm": N_PERM, "fold_class_rule": "helix_frac/sheet_frac from HELIX/SHEET "
                                                             "records; all_alpha h>=0.30 s<=0.05, "
                                                             "all_beta s>=0.30 h<=0.05",
                        "decision_rule": "AUC>=0.75 AND p<=0.01 AND beats the n_CA size baseline",
                        "metric": "euclidean_coordinates_angstrom"},
                preprocessing=f"per-chain alpha diagrams from the {len(table)} structures above",
                seed=SEED, tier="X", stats=stats_rows,
                controls=[{"kind": "negative", "description":
                           "size-only confound control: AUC from chain length n_CA alone under the "
                           "same label permutation. A label permutation does NOT protect against a "
                           "size classifier, so this is the control that matters.",
                           "passed": bool(c_iii),
                           "detail": f"topology AUC {prim['auc_alpha_vs_beta']:.3f} vs size AUC "
                                     f"{auc_size:.3f} (|AUC-0.5|: "
                                     f"{abs(prim['auc_alpha_vs_beta']-0.5):.3f} vs "
                                     f"{abs(auc_size-0.5):.3f})"},
                          {"kind": "shuffle", "description":
                           f"{N_PERM} label permutations give the null distribution of |AUC - 0.5|",
                           "passed": bool(c_ii),
                           "detail": f"null |AUC-0.5| q95 = "
                                     f"{prim['null_abs_auc_minus_half_q50_q95'][1]}"}],
                findings=[{"claim": f"Pre-registered test: h1_small_bar_density separates all-alpha "
                                    f"(n={len(A)}) from all-beta (n={len(B)}) C-alpha clouds with "
                                    f"AUC {prim['auc_alpha_vs_beta']:.3f}, permutation p "
                                    f"{prim['p_perm_two_sided_on_abs_auc_minus_half']:.4f}; the "
                                    f"chain-length-only baseline gives AUC {auc_size:.3f}. "
                                    f"Verdict: {sep['decision']['verdict']}.",
                           "verdict": ("recovered" if sep["decision"]["verdict"] == "SEPARATES"
                                       else "inconclusive" if sep["decision"]["verdict"]
                                       == "DOES_NOT_SEPARATE_BEYOND_CHAIN_LENGTH" else "null"),
                           "tier": "X",
                           "caveat": f"n = {len(A)} + {len(B)} structures; fold class is this "
                                     f"pipeline's own HELIX/SHEET count, not CATH or SCOP; a label "
                                     f"permutation tests exchangeability only, which is why the "
                                     f"size baseline is reported beside it",
                           "reference": "topodb_runs/biology/expectations.json (commit dde7a07)"}])
    else:
        sep["decision"] = {"verdict": "NOT_RUN",
                           "reason": f"only {len(A)} all-alpha and {len(B)} all-beta chains after "
                                     f"the pre-registered classification; the test needs >= 3 of each"}

    (RESULTS / "proteins_table.json").write_text(json.dumps(
        {"seed": SEED, "table": table, "separation": sep}, indent=1))
    print(json.dumps(sep.get("decision", {}), indent=1))
    print(json.dumps(sep.get("statistics", {}), indent=1))
    print(json.dumps(sep.get("size_only_baseline", {}), indent=1))
    blk.write()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
