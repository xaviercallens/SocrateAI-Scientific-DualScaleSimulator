"""Regression guard for defect D1 (miscalibrated chi2 p-values).

Every test here is written so that it FAILS against
``audit/reverse_zero/E5-cmb-tda/cmb_tda.py::coarse_stats`` and PASSES against
``tda_fixed.stats.coarse_stats_fixed``.  The tests that assert the old
behaviour are marked ``test_OLD_*``: they pin down the defect, so if someone
re-introduces it the new-behaviour tests below fail loudly.

Run:
  OMP_NUM_THREADS=1 timeout 590 prlimit --as=8589934592 -- \
    <venv>/bin/python -m pytest audit/tda_validation/tda_fixed/tests -q
"""
import importlib.util
import json
import os
import sys

import numpy as np
import pytest
from scipy.stats import kstest

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
WT = os.path.abspath(os.path.join(PKG, "..", "..", ".."))
sys.path.insert(0, os.path.dirname(PKG))

from tda_fixed import stats as fs  # noqa: E402

N_BINS = 8


@pytest.fixture(scope="module")
def original():
    spec = importlib.util.spec_from_file_location(
        "cmb_tda_original_for_tests",
        os.path.join(WT, "audit/reverse_zero/E5-cmb-tda/cmb_tda.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def make_ensemble(n=100, seed=0, dead_bin=7, atomic_bin=1, atomic_rate=0.74):
    """An ensemble with the structure the suite measured in F3: one bin that
    is identically constant, one bin that is a point mass in 74% of the sims,
    and six ordinary correlated Gaussian bins.

    Curves are returned at the raw length used by coarse_stats (8 * 5 points,
    constant inside each bin), so coarsen() reproduces the bin values exactly.
    """
    rng = np.random.default_rng(seed)
    A = rng.normal(size=(N_BINS, N_BINS))
    Sig = A @ A.T + 3.0 * np.eye(N_BINS)
    L = np.linalg.cholesky(Sig)

    def draw(m, rg):
        X = rg.normal(size=(m, N_BINS)) @ L.T
        X[:, dead_bin] = 1.0                                   # dead bin
        a = rg.random(m) < atomic_rate
        X[:, atomic_bin] = np.where(a, 0.0, np.abs(X[:, atomic_bin]) * 0.2)
        return X

    ens = draw(n, rng)
    test = draw(200, rng)
    expand = lambda X: np.repeat(X, 5, axis=1)  # noqa: E731
    return expand(ens), expand(test)


# ------------------------------------------------------------------ D1.1
def test_OLD_gives_a_dead_bin_a_weight_of_about_1e8(original):
    """The defect itself: a test vector that differs ONLY in a bin that is
    identically constant across the ensemble gets a p-value of ~0 because the
    1e-8 ridge is inverted instead of the direction being dropped."""
    ens, _ = make_ensemble(seed=1)
    v = ens.mean(axis=0).copy()
    v[7 * 5:8 * 5] += 1.0                     # differ only in the dead bin
    old = original.coarse_stats(ens, v, n_bins=N_BINS)
    assert old["data_chi2_hartlap"] > 1e6, old["data_chi2_hartlap"]
    assert old["p_value_chi2_survival"] < 1e-12
    assert old["n_bins"] == N_BINS            # df never adapts


def test_FIXED_reports_the_dead_bin_instead_of_inverting_it():
    ens, _ = make_ensemble(seed=1)
    v = ens.mean(axis=0).copy()
    v[7 * 5:8 * 5] += 1.0
    new = fs.coarse_stats_fixed(ens, v, n_bins=N_BINS, full_curve_diagnostics=False)
    assert new["data_chi2_hartlap"] < 1e3, new["data_chi2_hartlap"]
    assert new["p_value_chi2_survival"] > 1e-6
    # reported, not silently discarded (CLAUDE.md rule 3)
    assert "7" in new["dropped_bin_indices"] or 7 in new["dropped_bin_indices"]
    assert "7" in new["test_differs_in_dropped_bin"]
    assert new["test_differs_in_dropped_bin"]["7"]["test_value"] != \
        new["test_differs_in_dropped_bin"]["7"]["ensemble_value"]


# ------------------------------------------------------------------ D1.2
def test_df_is_the_retained_rank_not_the_nominal_bin_count(original):
    ens, test = make_ensemble(seed=2)
    old = original.coarse_stats(ens, test[0], n_bins=N_BINS)
    new = fs.coarse_stats_fixed(ens, test[0], n_bins=N_BINS, full_curve_diagnostics=False)
    assert old["n_bins"] == 8, "the original always used df = n_bins"
    assert new["df"] == new["n_bins_kept"] < 8
    assert new["df"] == new["retained_rank"]


def test_dead_and_atomic_bins_are_both_dropped_with_a_reason():
    ens, _ = make_ensemble(seed=3)
    new = fs.coarse_stats_fixed(ens, ens[0], n_bins=N_BINS, full_curve_diagnostics=False)
    assert new["dropped_bin_reasons"]["7"] == "zero_variance"
    assert new["dropped_bin_reasons"]["1"] == "atomic_modal_mass"
    assert new["n_bins_kept"] == 6


def test_no_ridge_is_added_to_the_covariance():
    """The 1e-8 ridge made the smallest eigenvalue exactly 1e-8 for a dead
    direction.  The fixed conditioning drops the direction instead."""
    C = np.diag([4.0, 1.0, 0.0])
    cinv, k, cond = fs.conditioned_inverse(C)
    assert k == 2
    assert cinv[2, 2] == 0.0, "a degenerate direction must get weight 0, not 1e8"
    assert cond == pytest.approx(4.0)


# ------------------------------------------------------------------ D1.3
def test_pvalues_are_uniform_after_the_fix_and_were_not_before(original):
    """The acceptance in miniature, POOLED over independent ensembles.

    Pooling matters: within one ensemble the test p-values are correlated
    (they share its estimated mean and covariance), so a KS test on a single
    ensemble is anti-conservative and would make this guard flaky.

    Measured over 12 independent ensembles of 40 test vectors
    (seeds 100..111), 480 p-values:
        OLD  KS p = 3.3e-23, KS stat 0.232, median p 0.718, rate 0.0583
        NEW  KS p = 1.7e-2,  KS stat 0.070, median p 0.557, rate 0.0396
        per-ensemble KS p < 0.05: old 11/12, new 0/12

    Note what this asserts and what it does NOT.  The fixed chi2 branch is
    approximately, not exactly, uniform: with an estimated covariance the
    exact small-sample reference is Hotelling's T^2 (an F distribution), not
    chi2, and the Hartlap factor corrects the mean of the inverse covariance
    rather than the whole distribution.  ``--step calibcheck`` measures the
    residual on exactly Gaussian data: 0.0577 instead of 0.05 at the nominal
    5% level.  The test therefore asserts a large, robust IMPROVEMENT and a
    bounded residual, not exact uniformity.  The distribution-free rank
    p-value has no such residual (calibcheck: 0.0511 vs an exact 0.0495).
    """
    p_old, p_new, ks_old_each, ks_new_each = [], [], [], []
    for s in range(12):
        ens, test = make_ensemble(n=100, seed=100 + s)
        test = test[:40]
        po = np.array([original.coarse_stats(ens, t, n_bins=N_BINS)["p_value_chi2_survival"]
                       for t in test])
        pn = np.array([fs.coarse_stats_fixed(ens, t, n_bins=N_BINS, with_rank_p=False,
                                             full_curve_diagnostics=False)["p_value_chi2_survival"]
                       for t in test])
        p_old.append(po); p_new.append(pn)
        ks_old_each.append(kstest(po, "uniform").pvalue)
        ks_new_each.append(kstest(pn, "uniform").pvalue)
    po = np.concatenate(p_old); pn = np.concatenate(p_new)
    ks_old = kstest(po, "uniform").pvalue
    ks_new = kstest(pn, "uniform").pvalue
    assert ks_old < 1e-6, "the old statistic is expected to be grossly miscalibrated (%g)" % ks_old
    assert ks_new > 1e-3, "the fixed statistic is far from uniform (%g)" % ks_new
    assert ks_new / ks_old > 1e6, "the fix must improve calibration by orders of magnitude"
    # shape, not just the tail count: the old median sits at 0.72, not 0.5
    assert abs(np.median(po) - 0.5) > 0.15, np.median(po)
    assert abs(np.median(pn) - 0.5) < 0.10, np.median(pn)
    assert kstest(pn, "uniform").statistic < 0.5 * kstest(po, "uniform").statistic
    n_bad_old = sum(k < 0.05 for k in ks_old_each)
    n_bad_new = sum(k < 0.05 for k in ks_new_each)
    assert n_bad_old >= 8, n_bad_old
    assert n_bad_new <= 3, n_bad_new


# ------------------------------------------------------------------ D1.4
def test_rank_p_direction_and_normalisation():
    ens, _ = make_ensemble(n=100, seed=5)
    at_centre = fs.coarse_stats_fixed(ens, ens.mean(axis=0), n_bins=N_BINS,
                                      full_curve_diagnostics=False)
    far = fs.coarse_stats_fixed(ens, ens.mean(axis=0) + 50.0, n_bins=N_BINS,
                                full_curve_diagnostics=False)
    assert at_centre["empirical_rank_p"] == pytest.approx(1.0)
    assert far["empirical_rank_p"] == pytest.approx(1.0 / 101)
    assert 0.0 < far["empirical_rank_p"] <= at_centre["empirical_rank_p"] <= 1.0


def test_rank_p_kept_bin_set_is_symmetric_in_the_data_and_the_sims():
    """Exchangeability requires the kept-bin set to come from the POOLED set.
    If it came from the sims alone, the data would be handled differently."""
    ens, test = make_ensemble(n=100, seed=6)
    r = fs.coarse_stats_fixed(ens, test[0], n_bins=N_BINS, full_curve_diagnostics=False)
    assert "kept_bin_indices_pooled" in r
    pool = np.vstack([np.array([fs.coarsen(c, N_BINS) for c in ens]),
                      fs.coarsen(test[0], N_BINS)[None, :]])
    assert r["kept_bin_indices_pooled"] == fs.select_bins(pool)[0]


# ------------------------------------------------------------------ D1.5
def test_negative_control_the_statistic_does_reject():
    """A check that cannot fail is not a check (CLAUDE.md rule 2)."""
    ens, test = make_ensemble(n=100, seed=7)
    # a shift of 5 is INSIDE the ensemble spread (measured: max chi2 p = 0.64),
    # so the statistic correctly does not reject; 10 is outside it.
    inside = [fs.coarse_stats_fixed(ens, t + 5.0, n_bins=N_BINS, full_curve_diagnostics=False)
              for t in test[:20]]
    assert max(p["p_value_chi2_survival"] for p in inside) > 0.05, \
        "a shift inside the ensemble spread must NOT be rejected"
    outside = [fs.coarse_stats_fixed(ens, t + 10.0, n_bins=N_BINS, full_curve_diagnostics=False)
               for t in test[:20]]
    assert all(p["p_value_chi2_survival"] < 1e-3 for p in outside)
    assert all(p["empirical_rank_p"] <= 1.0 / 101 + 1e-12 for p in outside)


def test_fail_closed_when_no_bin_survives():
    """All bins dead -> no chi2 is reported, rather than a confident number."""
    ens = np.ones((50, 40))
    r = fs.coarse_stats_fixed(ens, np.ones(40), n_bins=N_BINS, full_curve_diagnostics=False)
    assert r["n_bins_kept"] == 0
    assert r["p_value_chi2_survival"] is None
    assert "fail closed" in r["status"]


# ------------------------------------------------- the real F3 ensemble
def _committed_f3():
    import glob
    maps = []
    for fn in sorted(glob.glob(os.path.join(PKG, "..", "simple_suite", "cases", "F3chunk_*.json"))):
        maps += json.load(open(fn))["maps"]
    return sorted(maps, key=lambda m: m["i"])


@pytest.mark.parametrize("key", ["b0", "b1"])
def test_F3_acceptance_on_the_committed_curves(original, key):
    """The suite's own failing case: b0 and b1 chi2 p-values on the committed
    F3 curves.  Old: KS against uniform fails.  Fixed: it passes."""
    maps = _committed_f3()
    if len(maps) != 200:
        pytest.skip("committed F3 curves not present")
    S = np.array([m[key] for m in maps if m["i"] < 100], float)
    T = np.array([m[key] for m in maps if m["i"] >= 100], float)
    p_old = np.array([original.coarse_stats(S, t, n_bins=N_BINS)["p_value_chi2_survival"] for t in T])
    p_new = np.array([fs.coarse_stats_fixed(S, t, n_bins=N_BINS, with_rank_p=False,
                                            full_curve_diagnostics=False)["p_value_chi2_survival"]
                      for t in T])
    assert kstest(p_old, "uniform").pvalue <= 0.05, "the old statistic failed the suite's gate here"
    assert kstest(p_new, "uniform").pvalue > 0.05
    assert int(np.sum(p_new < 0.05)) <= 9
