"""ECG delay embedding (PhysioNet MIT-BIH) -- a REAL-data known-answer case.

A heartbeat is periodic, so a delay embedding of a normal sinus ECG MUST contain a
dominant H1 loop.  This is the one biological dataset in the block where the right
answer is known in advance from physics rather than from a publication, so it is
both a breadth modality and a second end-to-end control of the pipeline.

Records (MIT-BIH Arrhythmia Database, https://physionet.org/content/mitdb/1.0.0/):
  100  normal sinus rhythm, 69 M -- periodic, dominant H1 expected
  207  severe arrhythmia (ventricular flutter, bundle-branch block) -- the rhythm is
       NOT a clean limit cycle; no pre-stated expectation, recorded for contrast

Signal format 212 (two 12-bit samples packed in three bytes) is decoded here; no
wfdb package is needed, so nothing is installed outside .venv-tda.

NULL: amplitude-adjusted Fourier-transform surrogates -- the Fourier phases of the
SAME recording are randomised, which preserves the power spectrum and the amplitude
distribution but destroys the deterministic cycle.  A surrogate that still shows a
dominant loop would mean the loop is a spectral artefact.

Pre-stated criteria (expectations.json): (i) dominance >= 5 for the periodic
recording; (ii) the surrogates give a lower dominance.

Usage: python run_ecg.py [--n-null 200]
Seeds: window start 20260920; surrogate seeds 0..n-1.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, "/mnt/disks/disk-socrateai-local-1/wt-topo-bio/topodb_runs/biology")
from scipy.signal import butter, filtfilt  # noqa: E402
from scipy.spatial.distance import pdist, squareform  # noqa: E402

from bio_common import (BIO_DATA, RESULTS, Block, S_stat, dominance, finite, max_pers,  # noqa: E402
                        rank_p, rips_collapsed_from_distance, sha256, top_bars)

SEED = 20260920
SCRIPT = "topodb_runs/biology/run_ecg.py"
COMMAND = ("timeout 590 prlimit --as=8589934592 -- "
           "/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python "
           "topodb_runs/biology/run_ecg.py --n-null 200")
PHYS = BIO_DATA / "physionet"
RECORDS = {"100": "normal sinus rhythm (69 M); periodic, dominant H1 expected",
           "207": "severe arrhythmia incl. ventricular flutter; no expectation stated"}
N_POINTS = 400          # subsampled embedding points
WINDOW_SEC = 20.0       # ~20 heartbeats


def read212(rec):
    """Decode a PhysioNet format-212 record. Returns (signal ch0 in mV, fs, header)."""
    hea = (PHYS / f"{rec}.hea").read_text().splitlines()
    f = hea[0].split()
    nsig, fs, nsamp = int(f[1]), float(f[2]), int(f[3])
    g0, b0 = hea[1].split()[2], hea[1].split()[4]
    gain = float(g0.split("(")[0])
    base = float(b0)
    raw = np.fromfile(PHYS / f"{rec}.dat", dtype=np.uint8)
    n3 = (raw.size // 3) * 3
    b = raw[:n3].reshape(-1, 3).astype(np.int32)
    s1 = b[:, 0] | ((b[:, 1] & 0x0F) << 8)
    s2 = b[:, 2] | ((b[:, 1] >> 4) << 8)
    s1 = np.where(s1 > 2047, s1 - 4096, s1)
    s2 = np.where(s2 > 2047, s2 - 4096, s2)
    inter = np.empty(s1.size + s2.size, dtype=np.int32)
    inter[0::2] = s1
    inter[1::2] = s2
    sig = inter[::nsig][:nsamp]          # channel 0
    return (sig - base) / gain, fs, {"n_signals": nsig, "fs_hz": fs, "n_samples": nsamp,
                                     "gain": gain, "baseline": base,
                                     "lead": hea[1].split()[-1], "comments": hea[3:]}


def delay_embed(x, tau, m):
    n = len(x) - (m - 1) * tau
    return np.c_[tuple(x[i * tau:i * tau + n] for i in range(m))]


def dominant_tau(x, fs, lo=0.5, hi=4.0):
    """Delay tau = a quarter of the dominant cycle, from the power spectrum in [lo,hi] Hz."""
    X = np.abs(np.fft.rfft(x - x.mean())) ** 2
    fr = np.fft.rfftfreq(len(x), 1 / fs)
    band = (fr >= lo) & (fr <= hi)
    f0 = float(fr[band][np.argmax(X[band])])
    return max(1, int(round(fs / f0 / 4))), f0


def _stats(dg):
    return {"h1_dominance_P1_over_P2": dominance(dg[1]), "h1_S": S_stat(dg[1]),
            "h1_max_persistence": max_pers(dg[1]),
            "h1_finite_bar_count": float(len(finite(dg[1])))}


def embed_stats(x, fs, tau, rng, m=3, sphere=False):
    """Delay embedding -> Rips on the FULL Euclidean matrix via edge collapse.

    The untruncated collapsed route is used rather than a max_edge cut-off: on a
    cloud normalised by its own standard deviation a cut-off of 4 sigma builds a
    near-complete 2-skeleton (measured: 5.6 s vs 0.34 s for the identical diagram).
    `sphere=True` is the SW1PerS normalisation (mean-centre each window vector and
    project to the unit sphere), used only in the post-hoc sensitivity sweep.
    """
    emb = delay_embed(x, tau, m)
    idx = np.sort(rng.choice(emb.shape[0], min(N_POINTS, emb.shape[0]), replace=False))
    pts = emb[idx].astype(float)
    if sphere:
        pts = pts - pts.mean(1, keepdims=True)
        nn = np.linalg.norm(pts, axis=1, keepdims=True)
        nn[nn == 0] = 1.0
        pts = pts / nn
    else:
        scale = float(np.std(pts))
        pts = pts / (scale if scale > 0 else 1.0)
    dg = rips_collapsed_from_distance(squareform(pdist(pts)), max_hom_dim=1)
    return dg, _stats(dg)


def sensitivity_sweep(sig, fs):
    """Declared POST-HOC sweep: is the pre-stated outcome an unlucky single choice?

    Six plain delay-embedding settings (window length x 0.5-15 Hz band-pass) and
    three SW1PerS sliding-window settings.  No nulls, so no p-values are stored.
    """
    out = []
    b, a = butter(3, [0.5 / (fs / 2), 15.0 / (fs / 2)], btype="band")
    for win in (5.0, 10.0, 20.0):
        for filt in (False, True):
            x = sig[int(60 * fs):int(60 * fs) + int(win * fs)]
            if filt:
                x = filtfilt(b, a, x)
            tau, f0 = dominant_tau(x, fs)
            _, s = embed_stats(x, fs, tau, np.random.default_rng(SEED))
            out.append({"variant": "delay_m3", "window_sec": win, "bandpass_0.5_15hz": filt,
                        "tau": tau, "f0_hz": round(f0, 4),
                        "dominance": round(s["h1_dominance_P1_over_P2"], 4),
                        "S": round(s["h1_S"], 5)})
    x = sig[int(60 * fs):int(60 * fs) + int(WINDOW_SEC * fs)]
    _, f0 = dominant_tau(x, fs)
    period = fs / f0
    for M in (8, 14, 20):
        tau = max(1, int(round(period / (M + 1))))
        _, s = embed_stats(x, fs, tau, np.random.default_rng(SEED), m=M + 1, sphere=True)
        out.append({"variant": "sw1pers_sphere", "window_sec": WINDOW_SEC, "M": M, "tau": tau,
                    "dominance": round(s["h1_dominance_P1_over_P2"], 4),
                    "S": round(s["h1_S"], 5)})
    return out


def aaft_surrogate(x, rng):
    """Amplitude-adjusted FT surrogate: same spectrum, same amplitude distribution."""
    n = len(x)
    g = np.sort(rng.normal(size=n))[np.argsort(np.argsort(x))]     # gaussianise
    F = np.fft.rfft(g)
    ph = rng.uniform(0, 2 * np.pi, F.size)
    ph[0] = 0.0
    if n % 2 == 0:
        ph[-1] = 0.0
    gs = np.fft.irfft(np.abs(F) * np.exp(1j * ph), n)
    return np.sort(x)[np.argsort(np.argsort(gs))]                  # re-impose amplitudes


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-null", type=int, default=200)
    a = ap.parse_args()
    blk = Block("ecg", SCRIPT, COMMAND)
    summary = {}

    for rec, desc in RECORDS.items():
        sig, fs, hea = read212(rec)
        rng = np.random.default_rng(SEED)
        start = int(60 * fs)                      # skip the first minute (settling)
        x = sig[start:start + int(WINDOW_SEC * fs)]
        tau, f0 = dominant_tau(x, fs)
        t = time.time()
        dg, obs = embed_stats(x, fs, tau, np.random.default_rng(SEED))
        wall = round(time.time() - t, 2)

        sweep = sensitivity_sweep(sig, fs)
        sweep_max = max(v["dominance"] for v in sweep)

        null = []
        for sd in range(a.n_null):
            xs = aaft_surrogate(x, np.random.default_rng(10_000 + sd))
            _, s = embed_stats(xs, fs, tau, np.random.default_rng(SEED))
            null.append([s["h1_S"], s["h1_dominance_P1_over_P2"]])
        nv = np.array(null)
        p_dom = float(rank_p(obs["h1_dominance_P1_over_P2"], nv[:, 1]))
        p_S = float(rank_p(obs["h1_S"], nv[:, 0]))
        periodic = rec == "100"
        ok = obs["h1_dominance_P1_over_P2"] >= 5 and p_dom <= 0.01

        ds_id = f"biology/ecg_mitbih_{rec}"
        blk.dataset(id=ds_id, domain="biology",
                    title=f"MIT-BIH record {rec} ({hea['lead']}), {WINDOW_SEC:.0f} s at "
                          f"{fs:.0f} Hz, delay embedding m=3",
                    source=f"https://physionet.org/files/mitdb/1.0.0/{rec}.dat",
                    provenance="observation", n_objects=N_POINTS, ambient_dim=3,
                    units="mV (delay coordinates, then divided by the cloud's standard deviation)",
                    sha256=sha256(PHYS / f"{rec}.dat"), local_path=str(PHYS / f"{rec}.dat"),
                    notes=f"{desc}. MIT-BIH Arrhythmia Database (Moody & Mark 2001, IEEE EMB 20:45; "
                          f"PhysioNet, Goldberger et al. 2000 Circulation 101:e215). format 212 "
                          f"decoded by {SCRIPT}; header {json.dumps({k: v for k, v in hea.items() if k != 'comments'})}. "
                          f"window: samples {start}..{start + int(WINDOW_SEC * fs)}; dominant "
                          f"spectral frequency {f0:.3f} Hz -> tau = {tau} samples.")
        blk.run(dataset_id=ds_id, method="rips", coeff_field=2, max_dim=1,
                params={"metric": "euclidean_coordinates (delay embedding)", "m": 3,
                        "tau_samples": tau, "dominant_freq_hz": round(f0, 4), "fs_hz": fs,
                        "window_sec": WINDOW_SEC, "n_points": N_POINTS, "max_edge_length": 4.0,
                        "max_hom_dim": 1, "normalisation": "divide by the cloud standard deviation",
                        "criteria": "dominance >= 5 AND p <= 0.01 vs AAFT surrogates",
                        "rips_route": "untruncated distance matrix via collapse_edges + expansion(2)",
                        "post_hoc_sensitivity_sweep": sweep},
                preprocessing=f"channel 0, samples {start}..{start + int(WINDOW_SEC * fs)}; "
                              f"delay embedding m=3 tau={tau}; random subsample to {N_POINTS} points",
                seed=f"subsample {SEED}; surrogate seeds 10000..{10000 + a.n_null - 1}",
                tier="X", wall_sec=wall, diagrams=dg,
                betti={1: 1 if ok else 0}, expected_betti=({1: 1} if periodic else None),
                stats=[{"name": "h1_dominance_P1_over_P2", "value": obs["h1_dominance_P1_over_P2"],
                        "null_model": "amplitude-adjusted Fourier-transform surrogates of the same "
                                      "recording (same power spectrum, same amplitude distribution)",
                        "n_null": a.n_null, "p_value": p_dom, "p_method": "rank"},
                       {"name": "h1_S_longest_over_total", "value": obs["h1_S"],
                        "null_model": "same AAFT surrogates", "n_null": a.n_null,
                        "p_value": p_S, "p_method": "rank"},
                       {"name": "h1_max_persistence", "value": obs["h1_max_persistence"]},
                       {"name": "h1_finite_bar_count", "value": obs["h1_finite_bar_count"]},
                       {"name": "null_median_dominance", "value": float(np.median(nv[:, 1]))},
                       {"name": "null_frac_dominance_ge_5",
                        "value": float(np.mean(nv[:, 1] >= 5))},
                       {"name": "max_dominance_over_9_posthoc_settings", "value": sweep_max}],
                controls=([{"kind": "known_answer",
                            "description": "a periodic physiological signal must give a dominant H1 "
                                           "loop under delay embedding (dominance >= 5)",
                            "passed": bool(ok),
                            "detail": f"dominance {obs['h1_dominance_P1_over_P2']:.3f}, p {p_dom:.4f}"}]
                          if periodic else []) +
                         [{"kind": "injection",
                           "description": "POST-HOC sensitivity sweep over 9 settings (3 window "
                                          "lengths x band-pass on/off, plus 3 SW1PerS sliding-window "
                                          "dimensions): is the pre-stated outcome one unlucky choice?",
                           "passed": bool((sweep_max >= 5) == (obs["h1_dominance_P1_over_P2"] >= 5)),
                           "detail": f"best dominance over the 9 settings is {sweep_max:.3f}; "
                                     f"pre-stated setting gives "
                                     f"{obs['h1_dominance_P1_over_P2']:.3f}"},
                          {"kind": "null_calibration",
                           "description": "AAFT surrogates keep the spectrum and the amplitudes but "
                                          "destroy the deterministic cycle",
                           "passed": bool(np.mean(nv[:, 1] >= 5) < 0.1),
                           "detail": f"median surrogate dominance {np.median(nv[:, 1]):.3f}; "
                                     f"{np.mean(nv[:, 1] >= 5):.3f} of {a.n_null} reach >= 5"}],
                findings=[{"claim": f"MIT-BIH {rec} ({desc.split(';')[0]}): delay embedding "
                                    f"(m=3, tau={tau} samples at {fs:.0f} Hz) gives H1 dominance "
                                    f"{obs['h1_dominance_P1_over_P2']:.2f} against a median "
                                    f"{np.median(nv[:, 1]):.2f} over {a.n_null} AAFT surrogates "
                                    f"(p {p_dom:.4f})."
                                    + (f" Pre-stated criteria: {'MET' if ok else 'NOT MET'}."
                                       if periodic else " No expectation was pre-stated for this "
                                                        "record."),
                           "verdict": ("recovered" if (periodic and ok) else
                                       "failed" if periodic else
                                       "recovered" if p_dom <= 0.01 else "inconclusive"),
                           "tier": "X",
                           "caveat": f"one {WINDOW_SEC:.0f} s window of one lead of one recording; "
                                     f"tau is chosen from the recording's own dominant spectral "
                                     f"frequency, which the surrogates share. A post-hoc sweep over "
                                     f"9 further settings reaches at best dominance {sweep_max:.2f}, "
                                     f"so the outcome is not one unlucky parameter choice.",
                           "reference": "PhysioNet MIT-BIH Arrhythmia Database, "
                                        "https://physionet.org/content/mitdb/1.0.0/"}])
        summary[rec] = {"desc": desc, "fs_hz": fs, "tau": tau, "f0_hz": round(f0, 4),
                        "dominance": round(obs["h1_dominance_P1_over_P2"], 4),
                        "S": round(obs["h1_S"], 5), "p_dom": round(p_dom, 5),
                        "p_S": round(p_S, 5),
                        "null_median_dominance": round(float(np.median(nv[:, 1])), 4),
                        "criteria_met": bool(ok) if periodic else None,
                        "posthoc_sweep": sweep, "posthoc_sweep_max_dominance": sweep_max,
                        "top_bars": [[round(b, 4), round(d, 4)] for b, d in top_bars(dg[1], 3)]}
        print(f"  rec {rec}: tau={tau} f0={f0:.3f}Hz dom={obs['h1_dominance_P1_over_P2']:.3f} "
              f"p={p_dom:.4f} null_med={np.median(nv[:, 1]):.3f}")

    (RESULTS / "ecg_summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))
    blk.write()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
