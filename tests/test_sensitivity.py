"""Tests for Rambachan and Roth smoothness-restricted inference."""

from __future__ import annotations

import numpy as np
import pytest

from rxinc.diagnostics import detrend_event_study
from rxinc.estimators import EventStudyResult
from rxinc.sensitivity import (
    _folded_normal_cv,
    _second_differences,
    breakdown_m,
    flci,
    observed_pre_curvature,
)

REL = np.array([-4, -3, -2, 0, 1, 2, 3])  # reference -1 omitted
TIMELINE = np.arange(-4, 4)


def _result(coefs, sd=0.01, rel=REL, corr=0.0):
    k = len(rel)
    vcov = sd**2 * ((1 - corr) * np.eye(k) + corr * np.ones((k, k)))
    return EventStudyResult(rel_periods=rel, coefs=np.asarray(coefs, float),
                            ses=np.sqrt(np.diag(vcov)), vcov=vcov, n_obs=1, n_clusters=1)


def test_folded_normal_cv_matches_normal_when_unbiased():
    assert _folded_normal_cv(0.0, 0.05) == pytest.approx(1.959964, abs=1e-5)
    assert _folded_normal_cv(3.0, 0.05) > 3.0


def test_second_difference_matrix_pins_reference_at_zero():
    d = _second_differences(REL, -1)
    # Rows centred on -3, -2, -1, 0, 1, 2. The row centred on -1 is d(-2) + d(0).
    assert d.shape == (6, len(REL))
    row = d[2]
    assert row[list(REL).index(-2)] == 1 and row[list(REL).index(0)] == 1
    assert np.count_nonzero(row) == 2


def test_rejects_non_consecutive_periods():
    with pytest.raises(ValueError, match="consecutive"):
        flci(_result(np.zeros(4), rel=np.array([-4, -2, 0, 1])), 0.0)


def test_rejects_wrong_target_length():
    with pytest.raises(ValueError, match="post-period weights"):
        flci(_result(np.zeros(len(REL))), 0.0, target=[1.0, 0.0])


def test_m_zero_matches_straight_line_detrending():
    """With independent coefficients, both are the same GLS trend removal."""
    rng = np.random.default_rng(1)
    coefs = 0.004 * (REL + 1) + np.where(REL >= 0, 0.03, 0.0) + rng.normal(0, 0.005, len(REL))
    res = _result(coefs, sd=0.005)
    center = flci(res, 0.0).estimate
    assert center == pytest.approx(detrend_event_study(res).post_average()[0], abs=1e-10)


def test_m_zero_removes_any_straight_line_exactly():
    effect = np.where(REL >= 0, 0.05, 0.0)
    for slope in (-0.02, 0.0, 0.013):
        est = flci(_result(effect + slope * (REL + 1)), 0.0).estimate
        assert est == pytest.approx(0.05, abs=1e-10)


def test_interval_widens_with_m():
    res = _result(np.where(REL >= 0, 0.05, 0.0) + 0.003 * (REL + 1), corr=0.3)
    widths = [flci(res, m).upper - flci(res, m).lower for m in (0.0, 0.002, 0.005, 0.01)]
    assert all(b >= a - 1e-12 for a, b in zip(widths, widths[1:], strict=False))


def _random_violation(m, rng, timeline=TIMELINE, ref=-1):
    """A differential trend through zero at the reference whose slope changes by at most m."""
    n = len(timeline)
    r = int(np.flatnonzero(timeline == ref)[0])
    d = rng.uniform(-m, m, n)               # second difference at each interior point
    g = np.zeros(n - 1)                     # g[i] = path[i + 1] - path[i]
    g[r - 1] = rng.normal(0, 0.01)
    for i in range(r, n - 1):
        g[i] = g[i - 1] + d[i]
    for i in range(r - 2, -1, -1):
        g[i] = g[i + 1] - d[i + 1]
    path = np.zeros(n)
    for i in range(r, n - 1):
        path[i + 1] = path[i] + g[i]
    for i in range(r - 1, -1, -1):
        path[i] = path[i + 1] - g[i]
    lookup = dict(zip(timeline.tolist(), path.tolist(), strict=True))
    return np.array([lookup[int(t)] for t in REL])


@pytest.mark.parametrize("m", [0.0, 0.002, 0.004, 0.01])
def test_exact_coverage_over_every_allowed_violation(m):
    """The estimator is linear, so coverage under a given violation is exact, not simulated.

    Coverage must be at least 95% for every violation the restriction allows, and
    the computed worst-case bias must actually be reached, or the interval is loose.
    An earlier Monte Carlo version of this test produced 92.7% at M = 0 from noise;
    3,000 draws gave 95.1%.
    """
    from scipy.stats import norm
    tau = np.where(REL >= 0, np.array([0.0, 0.0, 0.0, 0.02, 0.04, 0.05, 0.05]), 0.0)
    truth = float(tau[REL >= 0].mean())
    ci = flci(_result(np.zeros(len(REL)), sd=0.01, corr=0.2), m)
    half = (ci.upper - ci.lower) / 2
    s_ = REL + 1.0
    rng = np.random.default_rng(0)
    violations = [0.006 * s_ + 0.5 * m * s_**2, 0.006 * s_ - 0.5 * m * s_**2]
    violations += [_random_violation(m, rng) for _ in range(200)]
    biases = [float(ci.weights @ (tau + v)) - truth for v in violations]
    coverage = [norm.cdf((half - b) / ci.sd) - norm.cdf((-half - b) / ci.sd) for b in biases]
    assert min(coverage) >= 0.95 - 1e-9
    assert max(abs(b) for b in biases) <= ci.max_bias + 1e-12
    assert max(abs(b) for b in biases) == pytest.approx(ci.max_bias, abs=1e-9)


def test_breakdown_is_zero_when_not_significant_even_under_linearity():
    assert breakdown_m(_result(np.zeros(len(REL)), sd=0.05)) == 0.0


def test_breakdown_is_positive_and_brackets_the_boundary():
    res = _result(np.where(REL >= 0, 0.06, 0.0), sd=0.005)
    m_star = breakdown_m(res)
    assert m_star > 0
    assert flci(res, m_star * 0.98).excludes_zero
    assert not flci(res, m_star * 1.05).excludes_zero


def test_observed_pre_curvature():
    res = _result([-0.08, -0.03, -0.01, 0, 0, 0, 0])
    # leads -4,-3,-2 and reference 0: second differences at -3 and -2
    # Centred on -3: -0.01 - 2(-0.03) + (-0.08) = -0.03. Centred on -2: 0 - 2(-0.01) + (-0.03) = -0.01.
    assert np.allclose(observed_pre_curvature(res), [-0.03, -0.01])


def _realistic():
    """Correlated covariance and a concave pre-trend, like the adoption event study."""
    rel = np.array([-5, -4, -3, -2, 0, 1, 2, 3, 4])
    coefs = np.array([-0.049, -0.022, -0.013, -0.004, -0.001, 0.067, 0.079, 0.070, 0.067])
    ses = np.array([0.0063, 0.0038, 0.0027, 0.0017, 0.0014, 0.0017, 0.0021, 0.0024, 0.0030])
    corr = 0.3 * np.ones((9, 9)) + 0.7 * np.eye(9)
    return EventStudyResult(rel, coefs, ses, corr * np.outer(ses, ses), 1, 1)


def test_half_length_never_shrinks_as_m_grows():
    """An earlier local-search optimiser produced a visible kink; the convex one must not."""
    res = _realistic()
    target = [0.5, 0.5, 0.0, 0.0, 0.0]
    hl = [(lambda c: c.upper - c.lower)(flci(res, m, target=target)) for m in np.linspace(0, 0.03, 31)]
    assert all(b >= a - 1e-9 for a, b in zip(hl, hl[1:], strict=False))


def test_no_random_trend_removing_weights_beat_the_optimum():
    from rxinc.sensitivity import _folded_normal_cv, _max_bias, _second_differences
    res = _realistic()
    target = np.array([0.5, 0.5, 0.0, 0.0, 0.0])
    m = 0.008
    best = flci(res, m, target=target)
    best_hl = (best.upper - best.lower) / 2
    rel = res.rel_periods
    pre, post = rel < 0, rel >= 0
    dist = (rel + 1).astype(float)
    d = _second_differences(rel, -1)
    a = dist[pre]
    w0 = a * (-(target @ dist[post])) / (a @ a)
    null = np.linalg.svd(a[None, :])[2][1:].T
    rng = np.random.default_rng(5)
    for _ in range(400):
        c = np.zeros(len(rel))
        c[post] = target
        c[pre] = w0 + null @ rng.normal(0, 0.5, null.shape[1])
        sd = float(np.sqrt(c @ res.vcov @ c))
        hl = sd * _folded_normal_cv(_max_bias(c, d, m) / sd, 0.05)
        assert hl >= best_hl * (1 - 1e-6)
