"""Test 3: topology of viral evolution (Chan, Carlsson & Rabadan 2013), pre-stated in expectations.json.

prep runs with /mnt/disks/disk-socrateai-local-1/venv-tdaval/bin/python (pyfamsa 0.7.0);
boot/final with the .venv-tda python (gudhi 3.13).
Usage (resumable; site-bootstrap caches per segment; run 'prep' once, then 'boot' per segment
(can be parallel), then 'final'):
  python test3_influenza.py prep
  python test3_influenza.py boot --segment 1 ... --segment 8   (one process per segment)
  python test3_influenza.py final
Seeds: genome subsample numpy RandomState(20260919); site bootstrap for segment s, replicate r uses
RandomState(50*(s-1) + r), r = 0..49 (global seeds 0..399).
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import argparse
import gzip
import json
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    from tda_common import (DATA_ROOT, alpha_h1_summary, rips_h1_summary, classical_mds, pval_ge, dump,
                            pipeline_sha256, resumable_seed_loop, Budget)
except ImportError:
    # 'prep' (alignment) runs in the separate venv /mnt/disks/disk-socrateai-local-1/venv-tdaval that has
    # pyfamsa but not gudhi; boot/final run in .venv-tda (gudhi 3.13) where this import succeeds.
    DATA_ROOT = "/mnt/disks/disk-socrateai-local-1/dualscale-data-r3/tda_validation/genetics"

    class Budget(Exception):
        pass

HERE = os.path.dirname(os.path.abspath(__file__))
FLU = os.path.join(DATA_ROOT, "influenza")
WORK = os.path.join(DATA_ROOT, "influenza", "work")
N_SUB = 300
SEED_SUB = 20260919
N_BOOT = 50
DELTA_B = 0.005
SEG_NAMES = {1: "PB2", 2: "PB1", 3: "PA", 4: "HA", 5: "NP", 6: "NA", 7: "M", 8: "NS"}


def parse_genomeset():
    groups, cur = [], []
    with gzip.open(os.path.join(FLU, "genomeset.dat.gz"), "rt") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                if cur:
                    groups.append(cur); cur = []
                continue
            p = line.split("\t")
            cur.append({"acc": p[0], "host": p[1], "seg": p[2], "subtype": p[3], "len": int(p[6]) if p[6].isdigit() else -1,
                        "name": p[7], "gid": p[-1]})
    if cur:
        groups.append(cur)
    return groups


def prep():
    groups = parse_genomeset()
    cand = []
    for g in groups:
        if not all(r["host"] == "Avian" and r["name"].startswith("Influenza A virus") for r in g):
            continue
        segs = {}
        for r in g:
            if r["seg"].isdigit():
                segs.setdefault(int(r["seg"]), []).append(r)
        if sorted(segs) == list(range(1, 9)) and all(len(v) == 1 for v in segs.values()):
            cand.append({s: segs[s][0] for s in range(1, 9)})
    med = {s: float(np.median([c[s]["len"] for c in cand])) for s in range(1, 9)}
    elig = [c for c in cand if all(c[s]["len"] >= 0.9 * med[s] for s in range(1, 9))]
    need = {c[s]["acc"] for c in elig for s in range(1, 9)}
    seqs, cur_acc, buf = {}, None, []
    with gzip.open(os.path.join(FLU, "influenza.fna.gz"), "rt") as f:
        for line in f:
            if line.startswith(">"):
                if cur_acc in need:
                    seqs[cur_acc] = "".join(buf).upper()
                parts = line[1:].split("|")
                cur_acc = parts[3].split(".")[0] if len(parts) > 3 else None
                buf = []
            else:
                buf.append(line.strip())
        if cur_acc in need:
            seqs[cur_acc] = "".join(buf).upper()
    elig = [c for c in elig if all(c[s]["acc"] in seqs for s in range(1, 9))]
    seen, uniq = set(), []
    for c in elig:
        key = "|".join(seqs[c[s]["acc"]] for s in range(1, 9))
        if key not in seen:
            seen.add(key); uniq.append(c)
    rng = np.random.RandomState(SEED_SUB)
    pick = sorted(rng.choice(len(uniq), size=N_SUB, replace=False).tolist())
    sub = [uniq[i] for i in pick]
    from pyfamsa import Aligner, Sequence
    import pyfamsa
    os.makedirs(WORK, exist_ok=True)
    info = {"n_genome_sets_total": len(groups), "n_avian_8segment_sets": len(cand), "median_len": med,
            "n_after_length_filter_and_seq_found": len(elig), "n_unique_concatenated": len(uniq), "n_sub": N_SUB,
            "seed_sub": SEED_SUB, "pyfamsa_version": pyfamsa.__version__,
            "strains": [c[1]["name"] for c in sub], "accessions": [[c[s]["acc"] for s in range(1, 9)] for c in sub],
            "subtypes": [c[1]["subtype"] for c in sub]}
    for s in range(1, 9):
        aligner = Aligner()
        msa = aligner.align([Sequence(str(i).encode(), seqs[c[s]["acc"]].encode()) for i, c in enumerate(sub)])
        rows = {int(m.id.decode()): m.sequence.decode() for m in msa}
        A = np.array([list(rows[i]) for i in range(N_SUB)])
        np.save(os.path.join(WORK, f"aln_seg{s}.npy"), A.astype("S1"))
        info[f"aln_len_seg{s}"] = int(A.shape[1])
    json.dump(info, open(os.path.join(WORK, "prep_info.json"), "w"), indent=1)
    print({k: v for k, v in info.items() if k not in ("strains", "accessions", "subtypes")})


def counts(A, w=None):
    """mismatch and shared-site counts between all pairs; only A/C/G/T count as sites."""
    L = A.shape[1]
    w = np.ones(L, dtype=np.float64) if w is None else w.astype(np.float64)
    V = np.zeros(A.shape, dtype=np.float64)
    M = np.zeros((A.shape[0], A.shape[0]))
    for b in (b"A", b"C", b"G", b"T"):
        Bb = (A == b).astype(np.float64)
        V += Bb
        M += (Bb * w) @ Bb.T
    shared = (V * w) @ V.T
    return shared - M, shared


def pdist_from(mm, sh):
    D = np.where(sh > 0, mm / np.maximum(sh, 1e-12), 1.0)
    np.fill_diagonal(D, 0.0)
    return D


def stats_from_D(D):
    r = rips_h1_summary(distance_matrix=D, k=10)
    pers = np.array(r["all_h1_persistence_sorted"])
    return {"A_max_h1_persistence": float(pers[0]) if pers.size else 0.0,
            "B_n_h1_ge_0.005": int(np.sum(pers >= DELTA_B)),
            "n_h1_bars": int(pers.size), "top_h1_bars_birth_death": r["top_h1_bars_birth_death"][:5],
            "n_h1_infinite_death": r["n_h1_infinite_death"]}


def load_aln(s):
    return np.load(os.path.join(WORK, f"aln_seg{s}.npy"))


def boot(seg, budget):
    A = load_aln(seg)
    L = A.shape[1]

    def f(r):
        rng = np.random.RandomState(N_BOOT * (seg - 1) + r)
        w = np.bincount(rng.randint(0, L, size=L), minlength=L)
        mm, sh = counts(A, w)
        st = stats_from_D(pdist_from(mm, sh))
        return [st["A_max_h1_persistence"], st["B_n_h1_ge_0.005"]]
    os.makedirs(os.path.join(DATA_ROOT, "cache"), exist_ok=True)
    v = resumable_seed_loop(os.path.join(DATA_ROOT, "cache", f"test3_boot_seg{seg}_n{N_BOOT}.json"), N_BOOT, f, budget)
    print("seg", seg, "done", len(v))


def final():
    info = json.load(open(os.path.join(WORK, "prep_info.json")))
    MM, SH = 0.0, 0.0
    seg_res = {}
    for s in range(1, 9):
        mm, sh = counts(load_aln(s))
        MM = MM + mm; SH = SH + sh
        D = pdist_from(mm, sh)
        st = stats_from_D(D)
        X3, _ = classical_mds(D, 3)
        a = alpha_h1_summary(X3)
        st["alpha_path_mds3"] = {k: a[k] for k in ("P1", "P2", "P1_over_P2", "S", "top_h1_bars_birth_death")}
        st["median_pdist"] = float(np.median(D[np.triu_indices(len(D), 1)]))
        seg_res[f"{s}_{SEG_NAMES[s]}"] = st
    Dc = pdist_from(MM, SH)
    conc = stats_from_D(Dc)
    X3, _ = classical_mds(Dc, 3)
    a = alpha_h1_summary(X3)
    conc["alpha_path_mds3"] = {k: a[k] for k in ("P1", "P2", "P1_over_P2", "S", "top_h1_bars_birth_death")}
    conc["median_pdist"] = float(np.median(Dc[np.triu_indices(len(Dc), 1)]))
    nullA, nullB = [], []
    for s in range(1, 9):
        v = json.load(open(os.path.join(DATA_ROOT, "cache", f"test3_boot_seg{s}_n{N_BOOT}.json")))
        assert len(v) == N_BOOT
        nullA += [x[0] for x in v]; nullB += [x[1] for x in v]
    A_c = conc["A_max_h1_persistence"]
    crit = {"i_A_concat_p_le_0.01_vs_pooled_single_segment_bootstrap": pval_ge(nullA, A_c) <= 0.01,
            "ii_A_HA_and_A_NA_le_half_A_concat": bool(seg_res["4_HA"]["A_max_h1_persistence"] <= 0.5 * A_c and
                                                      seg_res["6_NA"]["A_max_h1_persistence"] <= 0.5 * A_c),
            "iii_B_concat_gt_max_single_B": bool(conc["B_n_h1_ge_0.005"] > max(v["B_n_h1_ge_0.005"] for v in seg_res.values()))}
    out = {"pipeline_file_sha256": pipeline_sha256(), "prep_info": info, "single_segments": seg_res,
           "concatenated_8_segments": conc,
           "null_pooled_single_segment_site_bootstrap": {"n": len(nullA), "p_A_concat": pval_ge(nullA, A_c),
                                                         "A_q50_95_99_max": np.quantile(nullA, [0.5, 0.95, 0.99, 1]).tolist(),
                                                         "B_q50_95_max": np.quantile(nullB, [0.5, 0.95, 1]).tolist()},
           "criteria_rips_binding": crit, "verdict": "PASS" if all(crit.values()) else "FAIL"}
    dump(out, os.path.join(HERE, "test3_influenza_results.json"))
    print(out["verdict"], crit)
    for k, v in list(seg_res.items()) + [("CONCAT", conc)]:
        print(k, "A %.4f B %d nH1 %d medp %.3f" % (v["A_max_h1_persistence"], v["B_n_h1_ge_0.005"], v["n_h1_bars"], v["median_pdist"]),
              "alpha P1/P2 %.2f" % v["alpha_path_mds3"]["P1_over_P2"])
    print(out["null_pooled_single_segment_site_bootstrap"])


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["prep", "boot", "final"])
    ap.add_argument("--segment", type=int)
    ap.add_argument("--budget-sec", type=float, default=530.0)
    a = ap.parse_args()
    try:
        if a.stage == "prep":
            prep()
        elif a.stage == "boot":
            boot(a.segment, a.budget_sec)
        else:
            final()
    except Budget as e:
        print("INCOMPLETE (budget reached, re-invoke with the same arguments to resume):", e)
        sys.exit(3)
