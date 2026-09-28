"""Calibrated coarse-bin statistics for Betti/Euler curves (fixed replacement
for cmb_tda.coarse_stats).

DEFECT FIXED (D1, found 2026-09-19 by the known-answer suite on branch
loop/tda-simple, audit/tda_validation/simple_suite/results.json case F3):
the original ``audit/reverse_zero/E5-cmb-tda/cmb_tda.py::coarse_stats``
returned chi2 p-values that are not uniform under the null.  Three measured
causes:

  (a) ``df`` was always ``n_bins`` (8), even when a coarse bin was identically
      constant across the whole simulation ensemble.  Measured: b0 bin 7 took
      the value 1.0 in 100/100 sims (std 0.0); b1 bin 0 took 0.0 in 100/100.
  (b) ``cov = np.cov(...) + 1e-8 * np.eye(n_bins)`` followed by
      ``np.linalg.pinv`` *inverted* the ridge instead of dropping the dead
      direction, giving that direction a weight of about 1e8.
  (c) bins that are atomic but not dead were treated as Gaussian.  Measured:
      b1 bin 1 took the value 0.0 in 74/100 sims (std 0.1202), so a single
      extra count in that bin gives |z| ~ 6 and a chi2 survival p ~ 1e-6.

The fix is declared in ``fix_expectations.json`` (committed 1bc7ebc, before
this file was written).  In one sentence: drop dead and atomic bins using the
ENSEMBLE only, condition the kept block by eigenvalue truncation instead of a
ridge, set df to the retained rank, and return a distribution-free
pooled-exchangeable rank p-value alongside the chi2 p-value.

Tier of every number this module produces: X (numerics).  Nothing here is
"proved"; the rank p-value is uniform only under exchangeability of the
pooled candidates, and the negative control in
``run_validation.py`` shows that it does reject when that fails.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import chi2 as _chi2dist

__all__ = [
    "coarsen",
    "select_bins",
    "conditioned_inverse",
    "mahalanobis",
    "coarse_stats_fixed",
]

# ---------------------------------------------------------------- constants
# Fixed in fix_expectations.json (commit 1bc7ebc) BEFORE any acceptance number
# was computed.  They are not tuned afterwards.
ATOMIC_MODAL_MASS_MAX = 0.50
"""Drop a coarse bin whose most frequent single ensemble value occurs in more
than this fraction of the sims.  0.50 = 'a majority of the ensemble sits on
one value', so the Gaussian/chi2 reference is not applicable to that bin.
On the F3 ensemble the measured modal masses are
{1.00, 0.74, 0.17, 0.16, 0.14, 0.12, 0.08, ...}: any threshold in (0.18, 0.74]
selects the same bins.  run_validation.py reports that sensitivity band."""

EIGEN_TRUNCATION_TAU = 1e-10
"""Retain eigen-directions of the kept-block covariance with
lambda_i > tau * lambda_max.  This is a numerical-rank safeguard, not a
statistical one; it REPLACES the 1e-8 additive ridge, which inverted
degenerate directions instead of dropping them."""


# ------------------------------------------------------------------ helpers
def coarsen(curve, n_bins=8):
    """Average a curve into ``n_bins`` contiguous bins.

    Identical to ``cmb_tda.coarsen`` (no defect was found here); reimplemented
    so the fixed library has no import dependency on the original file.
    """
    c = np.asarray(curve, dtype=float)
    n = len(c)
    edges = np.linspace(0, n, n_bins + 1).astype(int)
    return np.array([c[edges[i]:edges[i + 1]].mean() for i in range(n_bins)])


def select_bins(ensemble, atomic_modal_mass_max=ATOMIC_MODAL_MASS_MAX):
    """Choose which coarse bins carry a usable Gaussian statistic.

    ``ensemble`` is ``(n_candidates, n_bins)``.  Selection uses the ensemble
    ONLY -- no test vector is passed in, so the kept set cannot be chosen to
    suit a particular answer.

    Returns ``(kept, dropped, reasons, values)``:
      kept     sorted list of kept bin indices
      dropped  sorted list of dropped bin indices
      reasons  {bin index: "zero_variance" | "atomic_modal_mass"}
      values   {bin index: the constant / modal ensemble value}

    Rule 1 (dead bin): ensemble std (ddof=1) exactly 0.0.
    Rule 2 (atomic bin): the most frequent single value occurs in a fraction
    of the ensemble strictly greater than ``atomic_modal_mass_max``.
    """
    E = np.asarray(ensemble, dtype=float)
    n, p = E.shape
    dropped, reasons, values = [], {}, {}
    for j in range(p):
        col = E[:, j]
        if col.std(ddof=1) == 0.0:
            dropped.append(j)
            reasons[j] = "zero_variance"
            values[j] = float(col[0])
            continue
        vals, counts = np.unique(col, return_counts=True)
        k = int(np.argmax(counts))
        modal_mass = counts[k] / float(n)
        if modal_mass > atomic_modal_mass_max:
            dropped.append(j)
            reasons[j] = "atomic_modal_mass"
            values[j] = float(vals[k])
    kept = [j for j in range(p) if j not in reasons]
    return kept, dropped, reasons, values


def conditioned_inverse(cov, tau=EIGEN_TRUNCATION_TAU):
    """Pseudo-inverse of a symmetric covariance by eigenvalue truncation.

    No ridge is added.  Eigen-directions with ``lambda_i <= tau * lambda_max``
    are DROPPED (weight 0), not inverted.  Returns
    ``(cov_inv, retained_rank, condition_number_of_retained_block)``.
    """
    C = np.asarray(cov, dtype=float)
    C = 0.5 * (C + C.T)
    w, V = np.linalg.eigh(C)
    if w.size == 0 or w.max() <= 0.0:
        return np.zeros_like(C), 0, float("inf")
    keep = w > tau * w.max()
    k = int(keep.sum())
    if k == 0:
        return np.zeros_like(C), 0, float("inf")
    Vk = V[:, keep]
    wk = w[keep]
    cov_inv = (Vk / wk) @ Vk.T
    cond = float(wk.max() / wk.min())
    return cov_inv, k, cond


def mahalanobis(vec, mean, cov_inv):
    d = np.asarray(vec, float) - np.asarray(mean, float)
    return float(d @ cov_inv @ d)


def _score(vec, ref, kept, tau):
    """Chi2 score of ``vec`` against reference set ``ref`` on the kept bins.

    Both are already restricted to the kept columns by the caller.  Returns
    ``(chi2, retained_rank, hartlap, cond)``.  The Hartlap factor uses the
    RETAINED RANK, not the nominal bin count.
    """
    n = ref.shape[0]
    mean = ref.mean(axis=0)
    cov = np.cov(ref, rowvar=False)
    cov = np.atleast_2d(cov)
    cov_inv, k, cond = conditioned_inverse(cov, tau)
    hartlap = (n - k - 2) / (n - 1) if n > k + 2 else float("nan")
    f = hartlap if np.isfinite(hartlap) else 1.0
    return mahalanobis(vec, mean, cov_inv * f), k, hartlap, cond


def coarse_stats_fixed(sim_curves, data_curve, n_bins=8,
                       atomic_modal_mass_max=ATOMIC_MODAL_MASS_MAX,
                       tau=EIGEN_TRUNCATION_TAU,
                       drop_atomic=True,
                       with_rank_p=True,
                       full_curve_diagnostics=True):
    """Calibrated replacement for ``cmb_tda.coarse_stats``.

    Parameters
    ----------
    sim_curves : (n, L) the null ensemble's full curves
    data_curve : (L,)   the curve under test
    drop_atomic : if False, apply ONLY the minimal enumerated fix (dead-bin
        drop + df = kept bins + eigenvalue conditioning).  This is column B of
        the validation table; the pre-registration predicts it does not meet
        the acceptance for b1.  Default True = the full declared fix.

    Returns a dict with BOTH ``p_value_chi2_survival`` and
    ``empirical_rank_p``, the kept-bin count, the dropped-bin list with
    reasons, and the flag ``test_differs_in_dropped_bin`` -- a dropped bin is
    reported, never silently discarded.

    The rank p-value is the pooled-exchangeable construction: the n sims and
    the data are pooled into n+1 candidates; the kept-bin set is chosen once
    from the POOLED set (symmetric in the data and the sims); each candidate
    is then scored against the other n members, with mean, covariance,
    conditioning and Hartlap all recomputed from those n.  Under
    exchangeability of the pool this is uniform on {1/(n+1), ..., 1}.
    """
    S = np.asarray(sim_curves, dtype=float)
    d = np.asarray(data_curve, dtype=float)
    sim_coarse = np.array([coarsen(c, n_bins) for c in S])
    data_coarse = coarsen(d, n_bins)
    n = sim_coarse.shape[0]

    # --- bin selection.  For the chi2 branch the ensemble alone defines the
    # kept set (the test vector is never consulted).  For the rank branch the
    # POOLED set defines it, so the data and every sim are treated identically.
    amm = atomic_modal_mass_max if drop_atomic else 1.0  # 1.0 => rule 2 never fires
    kept, dropped, reasons, dropped_values = select_bins(sim_coarse, amm)

    test_differs = {}
    for j in dropped:
        if data_coarse[j] != dropped_values[j]:
            test_differs[str(j)] = {"ensemble_value": dropped_values[j],
                                    "test_value": float(data_coarse[j]),
                                    "reason_bin_was_dropped": reasons[j]}

    out = {
        "fix": "tda_fixed.stats.coarse_stats_fixed",
        "fix_declared_in": "audit/tda_validation/tda_fixed/fix_expectations.json (commit 1bc7ebc)",
        "tier": "X",
        "n_sims": int(n),
        "n_bins_nominal": int(n_bins),
        "n_bins_kept": len(kept),
        "kept_bin_indices": list(map(int, kept)),
        "dropped_bin_indices": list(map(int, dropped)),
        "dropped_bin_reasons": {str(j): reasons[j] for j in dropped},
        "dropped_bin_ensemble_values": {str(j): dropped_values[j] for j in dropped},
        "test_differs_in_dropped_bin": test_differs,
        "atomic_modal_mass_max_used": float(amm),
        "eigen_truncation_tau": float(tau),
        "drop_atomic": bool(drop_atomic),
    }

    if not kept:
        out.update({"data_chi2_hartlap": float("nan"), "p_value_chi2_survival": None,
                    "df": 0, "retained_rank": 0, "hartlap_factor": float("nan"),
                    "empirical_rank_p": None,
                    "status": "no bin survived selection; no chi2 is reported (fail closed)"})
        return out

    ref = sim_coarse[:, kept]
    vec = data_coarse[kept]
    chi2v, k, hartlap, cond = _score(vec, ref, kept, tau)
    pval = float(_chi2dist.sf(chi2v, df=k)) if k > 0 else None
    out.update({
        "data_chi2_hartlap": float(chi2v),
        "p_value_chi2_survival": pval,
        "df": int(k),
        "retained_rank": int(k),
        "hartlap_factor": float(hartlap),
        "kept_block_condition_number": float(cond),
        "df_note": "df = retained rank of the kept block, NOT the nominal bin count "
                   "(the original used df = n_bins = %d)" % n_bins,
    })

    if with_rank_p:
        # --- pooled-exchangeable rank p.  Kept-bin set from the POOLED n+1.
        pool = np.vstack([sim_coarse, data_coarse[None, :]])
        kept_pool, dropped_pool, reasons_pool, _ = select_bins(pool, amm)
        out["kept_bin_indices_pooled"] = list(map(int, kept_pool))
        out["dropped_bin_indices_pooled"] = list(map(int, dropped_pool))
        if not kept_pool:
            out["empirical_rank_p"] = None
            out["rank_p_note"] = "no bin survived pooled selection (fail closed)"
        else:
            P = pool[:, kept_pool]
            m = P.shape[0]  # n + 1
            scores = np.empty(m)
            for i in range(m):
                others = np.delete(P, i, axis=0)
                scores[i] = _score(P[i], others, kept_pool, tau)[0]
            s_data = scores[-1]
            s_sims = scores[:-1]
            out["empirical_rank_p"] = float((np.sum(s_sims >= s_data) + 1) / (n + 1))
            out["pooled_exchangeable_score_data"] = float(s_data)
            out["pooled_exchangeable_score_sims_mean"] = float(s_sims.mean())
            out["pooled_exchangeable_score_sims_std"] = float(s_sims.std())
            out["rank_p_note"] = ("pooled-exchangeable: each of the n+1 candidates is scored "
                                  "against the other n, with mean/cov/conditioning/Hartlap "
                                  "recomputed each time; the kept-bin set comes from the pooled set")

    if full_curve_diagnostics:
        mu = S.mean(axis=0)
        sd = S.std(axis=0)
        out["sim_ensemble_full_curve_mean"] = mu.tolist()
        out["sim_ensemble_full_curve_std"] = sd.tolist()
        out["standardized_residual_full_curve"] = (
            (d - mu) / np.where(sd == 0, np.nan, sd)).tolist()
    return out
