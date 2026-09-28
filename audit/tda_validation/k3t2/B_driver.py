#!/usr/bin/env python
"""
Part B driver (tier X). For a named study, scan N upward; at each N run seed 0, and seeds 1, 2
only if seed 0 passes (N_min = smallest grid N at which seeds 0,1,2 all pass, pre-registered rule).
Primary field for the criterion: F3 (F2 also computed and reported).
Each run is a separate process: prlimit --as=8 GiB, timeout PER_RUN_TIMEOUT s.
Stops at N_min, or at the first run that fails by time/memory (largest feasible N recorded).
Run: B_TIMEOUT=560 B_CALL_BUDGET=300 timeout 595 prlimit --as=8589934592 -- <venv-python> B_driver.py <study>
  (resumable: re-run until stopped_because is not a pause); writes B_<study>.json
"""
import json, os, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
PY = "/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python"
MEM = str(8 * 1024 ** 3)
PER_RUN_TIMEOUT = int(os.environ.get("B_TIMEOUT", "900"))

STUDIES = {
    "T2_alpha": dict(space="T2", method="alpha", tau=1.8, expected=[1, 2, 1], grid=[25, 50, 100, 200, 400, 800, 1600, 3200]),
    "T2_rips": dict(space="T2", method="rips", tau=1.2, expected=[1, 2, 1], grid=[50, 100, 200, 400, 800, 1600, 3200, 6400]),
    "T3_rips": dict(space="T3", method="rips", tau=1.2, expected=[1, 3, 3, 1], grid=[400, 800, 1600, 3200, 6400, 12800, 25600, 51200]),
    "T3_rips_fps": dict(space="T3", method="rips", tau=1.2, expected=[1, 3, 3, 1], fps=True, grid=[400, 800, 1600, 3200, 6400, 12800]),
    "T4_rips": dict(space="T4", method="rips", tau=1.2, expected=[1, 4, 6, 4, 1], grid=[800, 1600, 3200, 6400, 12800, 25600]),
    "T4_rips_fps": dict(space="T4", method="rips", tau=1.2, expected=[1, 4, 6, 4, 1], fps=True, grid=[800, 1600, 3200, 6400, 12800]),
    "orbifold_rips": dict(space="orbifold", method="rips", tau=1.2, expected=[1, 0, 6], grid=[400, 800, 1600, 3200, 6400, 12800, 25600]),
    "orbifold_rips_fps": dict(space="orbifold", method="rips", tau=1.2, expected=[1, 0, 6], fps=True, grid=[400, 800, 1600, 3200, 6400, 12800]),
    "K3_rips": dict(space="K3", method="rips", tau=0.8, expected=[1, 0, 22], grid=[500, 1000, 2000, 4000, 8000, 16000, 32000],
                    save_diagram_dims=[2]),
    "K3_rips_fps": dict(space="K3", method="rips", tau=0.8, expected=[1, 0, 22], fps=True, grid=[500, 1000, 2000, 4000, 8000, 16000],
                        save_diagram_dims=[2]),
    # added after the advisor review (post Part-B commit 69baa5d): same-pipeline positive control for K3
    "quadric_rips": dict(space="quadric", method="rips", tau=0.8, expected=[1, 0, 2], grid=[500, 1000, 2000, 4000],
                         save_diagram_dims=[2], all_seeds=True, no_stop=True),
    # added after the advisor review: strong witness complexes (N = witnesses, landmarks = N/20 by farthest-point)
    "T2_witness": dict(space="T2", method="witness", witness_ratio=20, tau=0.9, expected=[1, 2, 1], grid=[1000, 2000, 4000, 8000]),
    "T3_witness": dict(space="T3", method="witness", witness_ratio=20, tau=0.9, expected=[1, 3, 3, 1], grid=[4000, 8000, 16000, 32000, 64000]),
    "T4_witness": dict(space="T4", method="witness", witness_ratio=20, tau=0.9, expected=[1, 4, 6, 4, 1],
                       grid=[8000, 16000, 32000, 64000, 128000, 256000]),
    "K3_witness": dict(space="K3", method="witness", witness_ratio=20, tau=0.8, expected=[1, 0, 22],
                       grid=[10000, 20000, 40000, 80000, 160000], save_diagram_dims=[2]),
    "null4m_rips": dict(space="null4m", method="rips", tau=0.8, expected=[1, 0, 0], grid=[500, 1000, 2000, 4000],
                        save_diagram_dims=[2], all_seeds=True, no_stop=True),
    "null4_rips": dict(space="null4", method="rips", tau=0.8, expected=[1, 0, 0], null_D=16, null_diam=1.414,
                       grid=[500, 1000, 2000, 4000, 8000, 16000, 32000], save_diagram_dims=[2], all_seeds=True, no_stop=True),
}


def run(spec):
    cmd = ["timeout", str(PER_RUN_TIMEOUT), "prlimit", f"--as={MEM}", "--", PY, os.path.join(HERE, "B_ph_worker.py"), json.dumps(spec)]
    t = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True)
    wall = round(time.time() - t, 1)
    if p.returncode != 0:
        return {"spec": spec, "failed": True, "returncode": p.returncode, "wall": wall,
                "stderr_tail": p.stderr[-600:], "reason": "timeout" if p.returncode == 124 else "error/memory"}
    r = json.loads(p.stdout.strip().splitlines()[-1])
    r["wall"] = wall
    return r


def main(name):
    st = STUDIES[name]
    base = {k: v for k, v in st.items() if k not in ("grid", "all_seeds", "no_stop")}
    base.setdefault("fields", [3, 2])
    out = {"study": name, "tier": "X", "config": st, "per_run_timeout_s": PER_RUN_TIMEOUT, "mem_cap_bytes": MEM,
           "runs": [], "N_min": None, "largest_N_completed": None, "stopped_because": None}
    path = os.path.join(HERE, f"B_{name}.json")
    done = {}
    if os.path.exists(path):  # resume: reuse completed (N, seed) runs (calls are chunked to fit tool time limits)
        old = json.load(open(path))
        for r in old["runs"]:  # failures are reused too (a timeout/memory failure is a result, not retried)
            done[(r["spec"]["N"], r["spec"]["seed"])] = r
    t_call = time.time()
    budget = float(os.environ.get("B_CALL_BUDGET", "1e9"))
    for N in st["grid"]:
        passes = []
        for seed in (0, 1, 2):
            spec = dict(base, N=N, seed=seed)
            if (N, seed) in done:
                r = done[(N, seed)]
            elif time.time() - t_call > budget:
                out["stopped_because"] = "call budget exhausted (resume by re-running)"
                json.dump(out, open(path, "w"), indent=1)
                print(name, "paused", flush=True)
                return out
            else:
                r = run(spec)
            out["runs"].append(r)
            json.dump(out, open(path, "w"), indent=1)
            if r.get("failed"):
                out["stopped_because"] = f"N={N} seed={seed}: {r['reason']} (rc={r['returncode']}, {r['wall']} s)"
                json.dump(out, open(path, "w"), indent=1)
                print(name, out["stopped_because"], flush=True)
                return out
            ok = r["by_field"]["F3"]["recovered"]
            passes.append(ok)
            out["largest_N_completed"] = N
            f3 = r["by_field"]["F3"]
            print(name, N, seed, "F3 recovered" if ok else "F3 not", round(f3["best_window_ratio"], 3), f3["max_beta_k"],
                  "F2", r["by_field"]["F2"]["recovered"], r["info"].get("n_simplices"), r["wall"], "s", r["maxrss_MB"], "MB", flush=True)
            if not ok and not st.get("all_seeds"):
                break
        if len(passes) == 3 and all(passes) and not st.get("no_stop"):
            out["N_min"] = N
            out["stopped_because"] = "N_min found"
            break
    if out["stopped_because"] is None:
        out["stopped_because"] = "grid exhausted"
    json.dump(out, open(path, "w"), indent=1)
    print(name, "N_min", out["N_min"], out["stopped_because"])
    return out


if __name__ == "__main__":
    main(sys.argv[1])
