"""Robust inference when parallel trends may fail (Rambachan and Roth, 2023).

A pre-trend test cannot prove parallel trends, and when it fails it says
nothing about how wrong the post-period estimates are. Rambachan and Roth
replace the assumption that parallel trends hold exactly with a bound on how
far they can be violated, and report what can still be concluded.

This module implements their **smoothness restriction**, Delta^SD(M): the
difference in trends between treated and control may follow any straight line
through the reference period, but its slope may change by at most ``M`` per
period. ``M = 0`` allows only straight-line trends, which is what
:func:`rxinc.diagnostics.detrend_event_study` assumes. Larger ``M`` allows the
trend to curve.

Inference uses the **fixed-length confidence interval** (FLCI) the paper
recommends for this restriction. The estimator is an affine combination of the
event-study coefficients that removes any straight-line trend exactly, with
weights on the pre-period coefficients chosen to make the interval as short as
possible. Its half-length combines sampling error with the worst-case bias over
all trends the restriction allows. Worst-case bias is a linear program; written
through its dual, the weight choice for a given bias budget becomes a convex
quadratic program, so the optimum is global rather than local.

The headline output is the **breakdown value**: the largest ``M`` at which the
interval still excludes zero. It answers "how much would the pre-existing trend
have had to bend for this effect to disappear?"

Reference: Rambachan, A. and Roth, J. (2023), A More Credible Approach to
Parallel Trends, Review of Economic Studies 90(5), 2555-2591. The authors'
implementation is the R package HonestDiD.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from scipy import optimize, stats

from rxinc.estimators import EventStudyResult


@dataclass
class FLCI:
    """A fixed-length confidence interval under Delta^SD(M).

    Attributes:
        m: The smoothness bound used.
        estimate: Point estimate of the target, with straight-line trends removed.
        lower: Lower confidence limit.
        upper: Upper confidence limit.
        sd: Standard deviation of the estimator.
        max_bias: Worst-case bias over all trends the restriction allows.
        weights: Weights on every event-study coefficient, in ``rel_periods`` order.
    """

    m: float
    estimate: float
    lower: float
    upper: float
    sd: float
    max_bias: float
    weights: np.ndarray

    @property
    def excludes_zero(self) -> bool:
        return self.lower > 0 or self.upper < 0


def _timeline(rel_periods: np.ndarray, reference: int) -> np.ndarray:
    periods = np.sort(np.append(rel_periods.astype(int), reference))
    if not np.array_equal(periods, np.arange(periods[0], periods[-1] + 1)):
        raise ValueError(
            "Smoothness restrictions need consecutive relative periods; "
            f"got {periods.tolist()}"
        )
    return periods


def _second_differences(rel_periods: np.ndarray, reference: int) -> np.ndarray:
    """Matrix D with (D delta)_k = delta_{t+1} - 2 delta_t + delta_{t-1}.

    Columns follow ``rel_periods``; the reference period is pinned at zero and
    has no column.
    """
    periods = _timeline(rel_periods, reference)
    col = {int(p): i for i, p in enumerate(rel_periods.astype(int))}
    rows = []
    for t in periods[1:-1]:
        row = np.zeros(len(rel_periods))
        for p, coef in ((t - 1, 1.0), (t, -2.0), (t + 1, 1.0)):
            if p != reference:
                row[col[int(p)]] += coef
        rows.append(row)
    return np.array(rows)


def _max_bias(c: np.ndarray, d: np.ndarray, m: float) -> float:
    """Largest c'delta over {delta : |D delta| <= m}, with delta_ref = 0."""
    if m == 0:
        return 0.0
    res = optimize.linprog(
        -c,
        A_ub=np.vstack([d, -d]),
        b_ub=np.full(2 * len(d), m),
        bounds=[(None, None)] * len(c),
        method="highs",
    )
    if res.status == 3:
        return np.inf
    if not res.success:
        raise RuntimeError(f"bias linear program failed: {res.message}")
    return float(-res.fun)


def _folded_normal_cv(bias_ratio: float, alpha: float) -> float:
    """The 1 - alpha quantile of |N(bias_ratio, 1)|."""
    b = abs(bias_ratio)

    def coverage(x: float) -> float:
        return stats.norm.cdf(x - b) - stats.norm.cdf(-x - b) - (1 - alpha)

    return float(optimize.brentq(coverage, 0.0, b + 10.0))


def _target(rel_periods: np.ndarray, target: Sequence[float] | None) -> np.ndarray:
    post = rel_periods >= 0
    if target is None:
        ell = np.zeros(post.sum())
        ell[:] = 1.0 / post.sum()
    else:
        ell = np.asarray(target, dtype=float)
        if len(ell) != post.sum():
            raise ValueError(f"target needs {post.sum()} post-period weights, got {len(ell)}")
    return ell


def _min_variance_given_bias(
    sigma: np.ndarray, d: np.ndarray, pre: np.ndarray, ell_full: np.ndarray, m: float, budget: float,
    start: np.ndarray | None = None,
) -> np.ndarray | None:
    """Minimum-variance weights whose worst-case bias under Delta^SD(m) is at most ``budget``.

    By linear-programming duality, max {c'delta : |D delta| <= m} <= budget holds
    exactly when c = D'(lp - lm) for some lp, lm >= 0 with m * sum(lp + lm) <= budget.
    That turns the bias constraint into linear constraints, so this is a convex
    quadratic program and its solution is the global optimum.
    """
    n, r = len(ell_full), d.shape[0]
    n_pre = int(pre.sum())

    def unpack(x):
        c = ell_full.copy()
        c[pre] = x[:n_pre]
        return c, x[n_pre:n_pre + r], x[n_pre + r:]

    def objective(x):
        c, _, _ = unpack(x)
        return float(c @ sigma @ c)

    def objective_grad(x):
        c, _, _ = unpack(x)
        g = np.zeros_like(x)
        g[:n_pre] = 2.0 * (sigma @ c)[pre]
        return g

    def eq(x):
        c, lp, lm = unpack(x)
        return d.T @ (lp - lm) - c

    sel = np.zeros((n, n_pre))
    sel[np.flatnonzero(pre), np.arange(n_pre)] = 1.0
    eq_jac = np.hstack([-sel, d.T, -d.T])

    def ineq(x):
        _, lp, lm = unpack(x)
        return np.array([budget - m * float(lp.sum() + lm.sum())])

    ineq_jac = np.concatenate([np.zeros(n_pre), -m * np.ones(2 * r)])[None, :]
    x0 = start if start is not None else np.concatenate([np.zeros(n_pre), np.full(2 * r, 1e-3)])
    res = optimize.minimize(
        objective, x0, jac=objective_grad, method="SLSQP",
        constraints=[{"type": "eq", "fun": eq, "jac": lambda x: eq_jac},
                     {"type": "ineq", "fun": ineq, "jac": lambda x: ineq_jac}],
        bounds=[(None, None)] * n_pre + [(0.0, None)] * (2 * r),
        options={"maxiter": 500, "ftol": 1e-16},
    )
    if not res.success:
        return None
    return res.x


def _min_achievable_bias(d: np.ndarray, pre: np.ndarray, ell_full: np.ndarray, m: float) -> float:
    """Smallest worst-case bias any trend-removing weights can achieve, by linear programming."""
    n, r = len(ell_full), d.shape[0]
    n_pre = int(pre.sum())
    sel = np.zeros((n, n_pre))
    sel[np.flatnonzero(pre), np.arange(n_pre)] = 1.0
    # variables [c_pre, lp, lm]; D'(lp - lm) - S c_pre = ell restricted to post (pre rows are 0)
    a_eq = np.hstack([-sel, d.T, -d.T])
    b_eq = ell_full * (~pre)
    cost = np.concatenate([np.zeros(n_pre), m * np.ones(2 * r)])
    res = optimize.linprog(cost, A_eq=a_eq, b_eq=b_eq,
                           bounds=[(None, None)] * n_pre + [(0.0, None)] * (2 * r), method="highs")
    if not res.success:
        raise RuntimeError(f"minimum-bias linear program failed: {res.message}")
    return float(res.fun)


def flci(
    result: EventStudyResult,
    m: float,
    target: Sequence[float] | None = None,
    alpha: float = 0.05,
    reference: int = -1,
) -> FLCI:
    """Fixed-length confidence interval for a post-period effect under Delta^SD(M).

    The weights minimise the interval's half-length, sd * cv(max bias / sd). For
    each bias budget B, the minimum-variance weights with worst-case bias at most
    B solve a convex quadratic program (:func:`_min_variance_given_bias`); a
    one-dimensional search over B then picks the shortest interval. This is
    globally optimal for each B, so the half-length is non-decreasing in M, which
    the breakdown search relies on.

    Args:
        result: Event study with its full covariance.
        m: Bound on the change in slope of the differential trend per period, in
            the outcome's units per period squared.
        target: Weights on the post-period coefficients defining the target.
            Defaults to their simple average.
        alpha: One minus the coverage level.
        reference: The omitted, normalised relative period.

    Returns:
        The interval, its centre, and the weights that produced it.
    """
    rel = result.rel_periods.astype(int)
    pre, post = rel < 0, rel >= 0
    if not pre.any():
        raise ValueError("Need pre-period coefficients to remove trends.")
    ell = _target(rel, target)
    d = _second_differences(rel, reference)
    sigma = result.vcov
    dist = (rel - reference).astype(float)
    ell_full = np.zeros(len(rel))
    ell_full[post] = ell

    # Minimum-variance weights that remove any straight line through the
    # reference period exactly: optimal when m = 0, and the end of the search
    # range when m > 0 (beyond its bias, extra budget cannot reduce variance).
    a = dist[pre]
    rhs = -float(ell @ dist[post])
    w0 = a * rhs / float(a @ a)
    null = np.linalg.svd(a[None, :])[2][1:].T
    s_pp, s_pq = sigma[np.ix_(pre, pre)], sigma[np.ix_(pre, post)]
    if null.shape[1]:
        z0 = -np.linalg.solve(null.T @ s_pp @ null, null.T @ (s_pp @ w0 + s_pq @ ell))
        w_star = w0 + null @ z0
    else:
        w_star = w0
    c_star = ell_full.copy()
    c_star[pre] = w_star

    def summarise(c: np.ndarray) -> FLCI:
        sd = float(np.sqrt(c @ sigma @ c))
        bias = _max_bias(c, d, m)
        est = float(c @ result.coefs)
        hl = sd * _folded_normal_cv(bias / sd, alpha)
        return FLCI(m=m, estimate=est, lower=est - hl, upper=est + hl, sd=sd, max_bias=bias, weights=c)

    if m == 0:
        return summarise(c_star)

    b_hi = _max_bias(c_star, d, m)
    b_lo = _min_achievable_bias(d, pre, ell_full, m)
    if b_hi <= b_lo * (1 + 1e-9):
        return summarise(c_star)

    cache: dict[float, np.ndarray] = {}

    def length(budget: float) -> float:
        x = _min_variance_given_bias(sigma, d, pre, ell_full, m, budget)
        if x is None:
            return np.inf
        c = ell_full.copy()
        c[pre] = x[: int(pre.sum())]
        cache[budget] = c
        sd = float(np.sqrt(c @ sigma @ c))
        return sd * _folded_normal_cv(budget / sd, alpha)

    grid = np.linspace(b_lo * (1 + 1e-6), b_hi, 24)
    vals = [length(b) for b in grid]
    k = int(np.argmin(vals))
    lo_b, hi_b = grid[max(k - 1, 0)], grid[min(k + 1, len(grid) - 1)]
    if hi_b > lo_b:
        res = optimize.minimize_scalar(length, bounds=(lo_b, hi_b), method="bounded",
                                       options={"xatol": (hi_b - lo_b) * 1e-4})
        best_b = res.x if res.fun <= vals[k] else grid[k]
    else:
        best_b = grid[k]
    if best_b not in cache:
        length(best_b)
    candidates = [cache[best_b], c_star]
    return min((summarise(c) for c in candidates), key=lambda f: f.upper - f.lower)


def breakdown_m(
    result: EventStudyResult,
    target: Sequence[float] | None = None,
    alpha: float = 0.05,
    reference: int = -1,
    tol: float = 1e-6,
) -> float:
    """Largest M at which the FLCI still excludes zero.

    Returns 0.0 when the interval includes zero even for straight-line trends.
    """
    if not flci(result, 0.0, target, alpha, reference).excludes_zero:
        return 0.0
    lo, hi = 0.0, 1e-3
    while flci(result, hi, target, alpha, reference).excludes_zero:
        lo, hi = hi, hi * 2
        if hi > 1e3:
            return np.inf
    while hi - lo > tol * max(1.0, hi):
        mid = (lo + hi) / 2
        if flci(result, mid, target, alpha, reference).excludes_zero:
            lo = mid
        else:
            hi = mid
    return lo


def observed_pre_curvature(result: EventStudyResult, reference: int = -1) -> np.ndarray:
    """Second differences of the estimated leads, for calibrating M.

    Rambachan and Roth suggest judging M against how much the pre-period trend
    itself bends. These are noisy estimates, not bounds.
    """
    rel = result.rel_periods.astype(int)
    pre = rel < 0
    leads = dict(zip(rel[pre].tolist(), result.coefs[pre].tolist(), strict=True))
    leads[reference] = 0.0
    ts = sorted(leads)
    return np.array([leads[t + 1] - 2 * leads[t] + leads[t - 1]
                     for t in ts[1:-1] if t - 1 in leads and t + 1 in leads])
