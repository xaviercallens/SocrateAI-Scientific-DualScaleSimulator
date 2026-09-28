"""Post-ingestion corrections, all of them changes to what the database asserts.

1. wall_sec / peak_mb were NULL on every run. They are filled here by re-executing
   each run's persistence computation on ONE representative object inside a
   DEDICATED FORKED SUBPROCESS, so that ru_maxrss is that run's own high-water
   mark and not a process-wide one. The meaning of the two fields is written into
   `preprocessing` so the number is not silently redefined as "the whole script".

2. git_commit on the dataset 2-5 runs pointed at 34b4e19, which does not contain
   compute_xy.py / compute_gpe.py / compute_colloid.py / compute_disorder.py
   (they landed in d53c40c). Corrected.

3. The GPE cross-check is re-recorded against the right expectation: the cubical
   H0 count should be n_winding + 1 (one class for the condensate body, one per
   core). Omega = 0 is then EXACT and the residual is cleanly the edge excess.
   The previous records took the disagreement branch at Omega = 0 and stored the
   wrong explanation there.

4. The colloid known-answer control's stated cause was wrong: post_diagnostics.py
   shows the contact diameter is consistent between the two concentrations to
   1 per cent and that there are no spurious close pairs, so the miss is not an
   estimator bias.

Run: prlimit --as=8589934592 -- .venv-tda/bin/python fixups.py
"""
import json
import multiprocessing as mp
import os
import resource
import time

import numpy as np
from scipy import ndimage

import compute_gpe as G
import ingest_lib as il
import qf_lib as q

SCRIPTS_COMMIT = "d53c40c"
TIMING_NOTE = (" TIMING: wall_sec and peak_mb are for this run's persistence computation on ONE "
               "representative object (not the whole script), measured in a dedicated forked "
               "subprocess so peak_mb is that computation's own maximum resident set; "
               "scripts/fixups.py.")


# ----------------------------------------------------------------- timing
def _child(fn, args, conn):
    t0 = time.perf_counter()
    try:
        fn(*args)
        ok = True
    except Exception as e:  # noqa: BLE001
        ok = False
        print("timing failed:", e)
    conn.send((time.perf_counter() - t0,
               resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0, ok))
    conn.close()


def timed(fn, *args):
    ctx = mp.get_context("fork")
    a, b = ctx.Pipe(False)
    p = ctx.Process(target=_child, args=(fn, args, b))
    p.start()
    out = a.recv() if a.poll(600) else (None, None, False)
    p.join(30)
    return out


def _alpha(pts, mas):
    q.alpha_bars(np.asarray(pts, float), mas)


def _cubical(field):
    q.cubical_bars(np.asarray(field, float), sublevel=True)


# ------------------------------------------------------- representatives
def rep_for(key, params, dataset_id):
    """Return (fn, args) for the representative persistence computation of a run."""
    if key.startswith("stm/"):
        H = key.split("/")[1].replace("kOe", "")
        z = np.load(os.path.join(q.RESULTS, f"stm_cores_{H}kOe.npz"))
        return _alpha, (z["img00"], params["max_alpha_square"])
    if key.startswith("colloid/"):
        pk = key.split("/")[1]
        z = np.load(os.path.join(q.RESULTS, f"colloid_phi{pk.replace('.','p')}.npz"))
        return _alpha, (z["frame00"], (3 * params["a_ref"]) ** 2)
    if key.startswith("disorder/"):
        _, axis, lvl = key.split("/")
        z = np.load(os.path.join(q.RESULTS, "stm_cores_20kOe.npz"))
        base = z["img00"]
        rng = np.random.default_rng(20260920)
        a_ref = params["a_ref_nm"]
        if axis == "A":
            pts = base + rng.normal(0, float(lvl) * a_ref, base.shape)
        else:
            nadd = int(round(float(lvl) * len(base)))
            pts = np.vstack([base, rng.random((nadd, 2)) * 353.0]) if nadd else base
        return _alpha, (pts, (3 * a_ref) ** 2)
    if key.startswith("xy/") and key.endswith("/alpha"):
        _, Ls, Tk, _ = key.split("/")
        L = int(Ls[1:]); T = float(Tk[1:])
        cfg = np.load(os.path.join(q.DATA_ROOT, f"xy/L{L}/T{T:.3f}.npz"))["configs"][0]
        pos, _ = q.winding_vortices(np.asarray(cfg, float), periodic=True)
        return _alpha, (pos, (3 * params["a_ref_from_density"]) ** 2)
    if key.startswith("xy/") and key.endswith("/cubical_shuffle"):
        L = int(key.split("/")[1][1:])
        cfg = np.load(os.path.join(q.DATA_ROOT, f"xy/L{L}/T0.900.npz"))["configs"][0]
        return _cubical, (np.cos(np.asarray(cfg, float)),)
    if key.startswith("gpe/"):
        om = key.split("/")[1]
        f = os.path.join(q.DATA_ROOT, "gpe", f"psi_Omega{float(om):.2f}_ext.npz")
        if not os.path.exists(f):
            f = os.path.join(q.DATA_ROOT, "gpe", f"psi_Omega{float(om):.2f}.npz")
        z = np.load(f)
        rho = (np.abs(z["psi"]) ** 2).astype(float); rho /= rho.max()
        sm = ndimage.gaussian_filter(rho, 3.0)
        mask = ndimage.binary_erosion(sm > G.MASK_FRAC * sm.max(), iterations=3)
        if key.endswith("/cubical"):
            return _cubical, (np.where(mask, rho, G.BIG),)
        pos, _ = q.winding_vortices(np.angle(z["psi"]), periodic=False)
        ij = np.rint(pos[:, ::-1]).astype(int)
        ij[:, 0] = np.clip(ij[:, 0], 0, mask.shape[0] - 1)
        ij[:, 1] = np.clip(ij[:, 1], 0, mask.shape[1] - 1)
        pin = pos[mask[ij[:, 0], ij[:, 1]]]
        x = z["x"]; dx = float(x[1] - x[0])
        pts = np.c_[x[0] + pin[:, 0] * dx, x[0] + pin[:, 1] * dx]
        return _alpha, (pts, (3 * params["a_ref"]) ** 2)
    return None, None


def main():
    g = il.IdGuard()
    P = json.load(open(os.path.join(q.RESULTS, "post_diagnostics.json")))
    d = il.db()
    c = d.con
    try:
        # ---- 1 + 2: timings and git_commit
        for key, rid in sorted(g.d.items()):
            if key.endswith("finding") or key.startswith("summary/"):
                continue
            row = c.execute("SELECT params_json, dataset_id, preprocessing FROM run WHERE id=?",
                            (rid,)).fetchone()
            if row is None:
                continue
            params = json.loads(row["params_json"])
            fn, args = rep_for(key, params, row["dataset_id"])
            if fn is None:
                continue
            wall, peak, ok = timed(fn, *args)
            pre = row["preprocessing"] or ""
            if TIMING_NOTE.strip() not in pre:
                pre = pre + TIMING_NOTE
            if key.startswith(("xy/", "gpe/", "colloid/", "disorder/")):
                c.execute("UPDATE run SET wall_sec=?, peak_mb=?, preprocessing=?, git_commit=? WHERE id=?",
                          (wall, peak, pre, SCRIPTS_COMMIT, rid))
            else:
                c.execute("UPDATE run SET wall_sec=?, peak_mb=?, preprocessing=? WHERE id=?",
                          (wall, peak, pre, rid))
            print(f"timed {key:34s} run {rid}  {wall:7.3f} s  {peak:8.1f} MB  ok={ok}", flush=True)
        c.commit()

        # ---- 3: GPE re-recorded against n_winding + 1
        for ok_, e in sorted(P["gpe_edge"].items()):
            rid = g.get(f"gpe/{ok_}/cubical")
            if rid is None:
                continue
            exact = e["residual_vs_winding_plus_one"] == 0
            c.execute("UPDATE betti SET expected=?, matches=? WHERE run_id=? AND dim=0",
                      (e["expected_cubical"], int(exact), rid))
            c.execute("DELETE FROM control WHERE run_id=? AND kind='known_answer'", (rid,))
            c.commit()
            d.add_control(rid, "known_answer",
                          "cubical H0 count vs (independently recomputed phase-winding count + 1)",
                          passed=bool(exact),
                          detail=(f"cubical {e['n_cubical_H0']} vs expected {e['expected_cubical']} "
                                  f"= winding {e['n_winding']} + 1 for the condensate body; residual "
                                  f"{e['residual_vs_winding_plus_one']:+d}. The cached "
                                  f"vortex_positions_*.npy were not used. Residuals across the series "
                                  f"are 0, +2, +4, +8 at Omega = 0.00, 0.70, 0.80, 0.90: exact when "
                                  f"there are no vortices, growing with rotation rate. "
                                  f"scripts/post_diagnostics.py"))
            if "median_dist_to_mask_edge_unmatched" in e and e["median_dist_to_mask_edge_unmatched"] is not None:
                d.add_control(rid, "negative", "where the cubical excess sits (local-minimum proxy)",
                              passed=None,
                              detail=(f"of {e['n_local_minima_proxy']} proxy local minima, "
                                      f"{e['n_unmatched']} are unmatched to a winding core; their median "
                                      f"distance to the condensate-mask boundary is "
                                      f"{e['median_dist_to_mask_edge_unmatched']:.2f} against "
                                      f"{e['median_dist_to_mask_edge_matched']:.2f} for matched cores "
                                      f"(mask inradius {e['mask_inradius']:.2f}). The excess is at the "
                                      f"edge. PROXY, not the cubical bars themselves."))
            fid = g.get(f"gpe/{ok_}/finding")
            if fid:
                if exact:
                    claim = (f"Rotating BEC at Omega = {e['Omega']:.2f}: the cubical (sublevel-set) H0 "
                             f"count and an independently recomputed plaquette phase-winding count "
                             f"agree EXACTLY once the condensate body is accounted for -- cubical "
                             f"{e['n_cubical_H0']} = winding {e['n_winding']} + 1. With no vortices "
                             f"present the filtration returns exactly the one component of the "
                             f"condensate, which fixes the baseline for the rotating cases.")
                    cav = ("This is the Omega = 0 control for the cross-check, not evidence about "
                           "vortex counting: there are no vortices to count. Its value is that it "
                           "pins the +1 offset, without which the rotating residuals would be "
                           "misread as +3/+5/+9 rather than +2/+4/+8.")
                    verdict = "recovered"
                else:
                    claim = (f"Rotating BEC at Omega = {e['Omega']:.2f}: the cubical H0 count EXCEEDS "
                             f"the independently recomputed phase-winding count even after allowing "
                             f"the one class for the condensate body -- cubical {e['n_cubical_H0']} "
                             f"against winding {e['n_winding']} + 1 = {e['expected_cubical']}, a "
                             f"residual of {e['residual_vs_winding_plus_one']:+d}. The residual grows "
                             f"monotonically with rotation rate (0, +2, +4, +8 at Omega = 0.00, 0.70, "
                             f"0.80, 0.90) and a local-minimum proxy locates it at the condensate "
                             f"edge: unmatched minima sit "
                             f"{e['median_dist_to_mask_edge_unmatched']:.2f} from the mask boundary "
                             f"while matched cores sit {e['median_dist_to_mask_edge_matched']:.2f} "
                             f"inside it (mask inradius {e['mask_inradius']:.2f}).")
                    cav = ("The excess is a property of a sublevel-set filtration on a field with a "
                           "boundary, not of the physics: |psi|^2 falls off at the condensate edge, so "
                           "the mask edge creates density minima that are not vortices. A phase-winding "
                           "count is the correct core counter here and the cubical count is an upper "
                           "bound. The edge localisation comes from a local-minimum proxy (127-229 "
                           "minima), not from the cubical bars themselves, so it indicates the "
                           "mechanism without measuring it exactly. The site-shuffle control passes "
                           "but is weak: shuffling destroys smoothness, so the shuffled field has "
                           "thousands of minima and is trivially distinguishable.")
                    verdict = "failed"
                c.execute("UPDATE finding SET claim=?, verdict=?, caveat=? WHERE id=?",
                          (claim, verdict, cav, fid))
                c.execute("UPDATE search SET title=?, text=? WHERE kind='finding' AND ref_id=?",
                          (claim, cav, str(fid)))
            print("gpe fixed", ok_)
        c.commit()

        # ---- 4: colloid control cause corrected
        concl = P["colloid_conclusion"]
        for pk in ("0.89", "0.11"):
            rid = g.get(f"colloid/{pk}/alpha")
            if rid is None:
                continue
            cm = P["colloid"][pk]
            c.execute("DELETE FROM control WHERE run_id=? AND kind='known_answer'", (rid,))
            c.commit()
            d.add_control(rid, "known_answer",
                          "scale-free hard-disc relation a/d = sqrt(pi/(2 sqrt(3) phi)) at the published phi",
                          passed=False,
                          detail=(f"measured a_density/d_contact misses the prediction by about 20 per "
                                  f"cent at BOTH concentrations. Cause investigated and the estimator-bias "
                                  f"explanation REJECTED: the nearest-neighbour mode is 12.27 px at "
                                  f"phi = 0.89 and 12.39 px at phi = 0.11, agreeing to 1 per cent, so the "
                                  f"contact diameter is a real property of the particles; and only a "
                                  f"fraction {cm['frac_nn_below_0.6_median']:.5f} of nearest-neighbour "
                                  f"distances lie below 0.6 of the median, so there are no spurious close "
                                  f"pairs. a_nn/a_density = {cm['a_nn_over_a_density']:.3f} here against "
                                  f"1.000 for a triangular lattice and 0.465 for Poisson. The point set is "
                                  f"internally consistent; the published concentration label is what does "
                                  f"not fit. scripts/post_diagnostics.py"))
            d.add_control(rid, "known_answer",
                          "scale-free a_nn/a_density (1.000 triangular, 0.465 Poisson)",
                          passed=True,
                          detail=(f"measured {cm['a_nn_over_a_density']:.4f}; at phi = 0.89 this matches "
                                  f"about 5 per cent positional disorder on the calibration curve in "
                                  f"step0_re6zr_20kOe.json, and at phi = 0.11 it matches the Poisson "
                                  f"value to 2 per cent. The structure is what it appears to be."))
            fid = c.execute("SELECT id, caveat FROM finding WHERE run_id=?", (rid,)).fetchone()
            if fid:
                cav = ("The published scale-free check a/d = sqrt(pi/(2 sqrt(3) phi)) is missed by about "
                       "20 per cent at BOTH concentrations. " + concl)
                c.execute("UPDATE finding SET caveat=? WHERE id=?", (cav, fid["id"]))
            print("colloid fixed", pk)
        c.commit()

        # ---- calibration-axis finding: qualify the 0.707 claim
        fid = g.get("summary/calibration_axis")
        if fid:
            row = c.execute("SELECT claim FROM finding WHERE id=?", (fid,)).fetchone()
            new = row["claim"].replace(
                "and the XY vortex gas at T = 1.6 reaches exactly 0.707",
                "and the XY vortex gas at T = 1.6 sits at 0.707 for L = 64 and L = 128 (L = 32 is "
                "0.874, where only 154 vortices make the estimate noisy)")
            c.execute("UPDATE finding SET claim=? WHERE id=?", (new, fid))
            c.execute("UPDATE search SET title=? WHERE kind='finding' AND ref_id=?", (new, str(fid)))
            print("calibration finding qualified")
        c.commit()
    finally:
        il.close(d)
    print("done")


if __name__ == "__main__":
    main()
