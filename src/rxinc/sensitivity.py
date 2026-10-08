"""Robust inference when parallel trends may fail (Rambachan and Roth, 2023).

A pre-trend test cannot prove parallel trends, and when it fails it says
nothing about how wrong the post-period estimates are. Rambachan and Roth
replace the assumption that parallel trends hold exactly with a bound on how
far they can be violated, and report what can still be concluded.

The core restriction is **smoothness**, Delta^SD(M): the difference in trends
between treated and control may follow any straight line through the reference
period, but its slope may change by at most ``M`` per period. ``M = 0`` allows
only straight-line trends, which is what
:func:`rxinc.diagnostics.detrend_event_study` assumes. Shape and sign
restrictions (:data:`RESTRICTIONS`) narrow this further, for example to trends
that can only flatten.

Inference uses the **fixed-length confidence interval** the paper recommends for
smoothness, generalised to one-sided restrictions. The estimator is an affine
combination of the event-study coefficients that removes any straight-line
trend exactly; its pre-period weights are chosen to make the interval as short
as possible. Worst-case bias over the allowed trends is a linear program, and the
weights are found by cutting planes over the few free pre-period weights.

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


def _lp_extreme(objective: np.ndarray, g: np.ndarray, h: np.ndarray, m: float) -> np.ndarray:
    """Vertex maximising objective'delta over {G delta <= h}.

    Solved with free variables, which is exact in the ordinary case. Some
    restrictions (M = 0, or sign and monotonicity restrictions) leave an
    unbounded ray along the straight-line direction; trend-removing weights are
    orthogonal to it only to rounding error, and the solver can then fail on an
    almost flat objective over an infinite ray. Only in that case it retries
    inside a box. Under these restrictions some optimal trend always lies within
    about M * T**2 of zero, so a box of that order cannot move the optimum, and
    keeping it small preserves numerical precision: an earlier fixed box of
    +/-10,000 cost about 1e-12 of accuracy on every problem, enough for a
    worst-case bias to fall short of an attainable one.
    """
    res = optimize.linprog(-objective, A_ub=g, b_ub=h, bounds=[(None, None)] * len(objective), method="highs")
    if res.success:
        return res.x
    box = 1.0 + 4.0 * m * (len(objective) + 1) ** 2
    res = optimize.linprog(-objective, A_ub=g, b_ub=h, bounds=[(-box, box)] * len(objective), method="highs")
    if not res.success:
        raise RuntimeError(f"bias linear program failed: {res.message}")
    return res.x

RESTRICTIONS = (
    "smooth", "concave", "convex", "increasing", "decreasing", "positive_bias", "negative_bias",
)
"""Restrictions on the differential trend, each combined with smoothness bound M.

``smooth``          slope changes by at most M per period, either direction (Delta^SD).
``concave``         slope can only fall, by at most M per period: the trend can flatten,
                    never steepen. Appropriate when pre-period bending is consistently
                    negative, as with decelerating adoption.
``convex``          slope can only rise, by at most M per period.
``increasing``      smoothness, and the trend never falls (Rambachan and Roth's Delta^SDI).
``decreasing``      smoothness, and the trend never rises (Delta^SDD).
``positive_bias``   smoothness, and post-period bias is non-negative (Delta^SDPB).
``negative_bias``   smoothness, and post-period bias is non-positive.
"""


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


@dataclass
class RestrictedCI:
    """Confidence interval under a shape or sign restriction.

    Attributes:
        restriction: Which restriction.
        m: The smoothness bound.
        estimate: The trend-removing point estimate.
        lower: Lower confidence limit.
        upper: Upper confidence limit.
        sd: Standard deviation of the estimator.
        bias_low: Smallest bias the restriction allows for these weights.
        bias_high: Largest bias the restriction allows for these weights.
        weights: Weights on every event-study coefficient.
    """

    restriction: str
    m: float
    estimate: float
    lower: float
    upper: float
    sd: float
    bias_low: float
    bias_high: float
    weights: np.ndarray

    @property
    def excludes_zero(self) -> bool:
        return self.lower > 0 or self.upper < 0


# ---------------------------------------------------------------------------
# Building blocks
# ---------------------------------------------------------------------------


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


def _first_differences(rel_periods: np.ndarray, reference: int) -> np.ndarray:
    """Matrix F with (F delta)_k = delta_t - delta_{t-1}, reference pinned at zero."""
    periods = _timeline(rel_periods, reference)
    col = {int(p): i for i, p in enumerate(rel_periods.astype(int))}
    rows = []
    for t in periods[1:]:
        row = np.zeros(len(rel_periods))
        for p, coef in ((t, 1.0), (t - 1, -1.0)):
            if p != reference:
                row[col[int(p)]] += coef
        rows.append(row)
    return np.array(rows)


def _polyhedron(rel: np.ndarray, reference: int, m: float, restriction: str) -> tuple[np.ndarray, np.ndarray]:
    """The restriction as {delta : G delta <= h}."""
    if restriction not in RESTRICTIONS:
        raise ValueError(f"restriction must be one of {RESTRICTIONS}, got {restriction!r}")
    d = _second_differences(rel, reference)
    r = len(d)
    g = np.vstack([d, -d])
    if restriction == "concave":
        return g, np.concatenate([np.zeros(r), np.full(r, m)])
    if restriction == "convex":
        return g, np.concatenate([np.full(r, m), np.zeros(r)])
    h = np.full(2 * r, m)
    if restriction in ("increasing", "decreasing"):
        f = _first_differences(rel, reference)
        sign = -1.0 if restriction == "increasing" else 1.0
        return np.vstack([g, sign * f]), np.concatenate([h, np.zeros(len(f))])
    if restriction in ("positive_bias", "negative_bias"):
        post = np.zeros((int((rel >= 0).sum()), len(rel)))
        post[np.arange(len(post)), np.flatnonzero(rel >= 0)] = 1.0
        sign = -1.0 if restriction == "positive_bias" else 1.0
        return np.vstack([g, sign * post]), np.concatenate([h, np.zeros(len(post))])
    return g, h


def _max_bias(c: np.ndarray, d: np.ndarray, m: float) -> float:
    """Largest c'delta over {delta : |D delta| <= m}, with delta_ref = 0."""
    if m == 0:
        return 0.0
    return float(c @ _lp_extreme(c, np.vstack([d, -d]), np.full(2 * len(d), m), m))


def _bias_extremes(
    c: np.ndarray, g: np.ndarray, h: np.ndarray, m: float = 0.0
) -> tuple[float, float, np.ndarray, np.ndarray]:
    """Smallest and largest c'delta over {G delta <= h}, with the vertices attaining them."""
    d_hi = _lp_extreme(c, g, h, m)
    d_lo = _lp_extreme(-c, g, h, m)
    return float(c @ d_lo), float(c @ d_hi), d_lo, d_hi


def _folded_normal_cv(bias_ratio: float, alpha: float) -> float:
    """The 1 - alpha quantile of |N(bias_ratio, 1)|."""
    b = abs(bias_ratio)

    def coverage(x: float) -> float:
        return stats.norm.cdf(x - b) - stats.norm.cdf(-x - b) - (1 - alpha)

    return float(optimize.brentq(coverage, 0.0, b + 10.0))


def _tight_cv(width_ratio: float, alpha: float) -> float:
    """x solving Phi(x + w) - Phi(-x) = 1 - alpha.

    With bias known to lie in [b_lo, b_hi], the interval
    [est - b_hi - x*sd, est - b_lo + x*sd] has coverage at least 1 - alpha, with
    equality at either end of the bias range. When the range is symmetric this is
    exactly the fixed-length interval, so the two constructions agree.
    """
    w = max(width_ratio, 0.0)

    def coverage(x: float) -> float:
        return stats.norm.cdf(x + w) - stats.norm.cdf(-x) - (1 - alpha)

    return float(optimize.brentq(coverage, -w, 10.0))


def _target(rel_periods: np.ndarray, target: Sequence[float] | None) -> np.ndarray:
    post = rel_periods >= 0
    if target is None:
        return np.full(int(post.sum()), 1.0 / post.sum())
    ell = np.asarray(target, dtype=float)
    if len(ell) != post.sum():
        raise ValueError(f"target needs {post.sum()} post-period weights, got {len(ell)}")
    return ell


def _trend_removing_basis(
    rel: np.ndarray, pre: np.ndarray, post: np.ndarray, ell: np.ndarray, sigma: np.ndarray, reference: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Parameterise pre-period weights that remove straight lines through the reference.

    Returns ``(w0, null, w_star)``: every admissible weight vector is
    ``w0 + null @ z``, and ``w_star`` is the minimum-variance one.
    """
    dist = (rel - reference).astype(float)
    a = dist[pre]
    w0 = a * (-float(ell @ dist[post])) / float(a @ a)
    null = np.linalg.svd(a[None, :])[2][1:].T
    s_pp, s_pq = sigma[np.ix_(pre, pre)], sigma[np.ix_(pre, post)]
    if null.shape[1]:
        z = -np.linalg.solve(null.T @ s_pp @ null, null.T @ (s_pp @ w0 + s_pq @ ell))
        return w0, null, w0 + null @ z
    return w0, null, w0


# ---------------------------------------------------------------------------
# The optimiser
# ---------------------------------------------------------------------------


def _min_width(g: np.ndarray, h: np.ndarray, pre: np.ndarray, ell_full: np.ndarray) -> tuple[float, np.ndarray]:
    """Smallest achievable bias-range width over admissible weights, and weights attaining it.

    By LP duality, max c'delta = min {h'lam : G'lam = c, lam >= 0} and similarly
    for the minimum, so the width is a linear program in the weights and both
    dual vectors. The interval length usually keeps falling as the width
    budget tightens, so this end of the range is often the optimum and must be
    located exactly rather than approached on a grid: an earlier version
    searched from zero, wasted most of its grid on infeasible widths, and ended
    up to 0.6% longer than necessary.
    """
    n, k = len(ell_full), len(g)
    n_pre = int(pre.sum())
    sel = np.zeros((n, n_pre))
    sel[np.flatnonzero(pre), np.arange(n_pre)] = 1.0
    a_eq = np.vstack([np.hstack([-sel, g.T, np.zeros((n, k))]),
                      np.hstack([sel, np.zeros((n, k)), g.T])])
    b_eq = np.concatenate([ell_full, -ell_full])
    cost = np.concatenate([np.zeros(n_pre), h, h])
    res = optimize.linprog(cost, A_eq=a_eq, b_eq=b_eq,
                           bounds=[(None, None)] * n_pre + [(0.0, None)] * (2 * k), method="highs")
    if not res.success:
        raise RuntimeError(f"minimum-width linear program failed: {res.message}")
    c = ell_full.copy()
    c[pre] = res.x[:n_pre]
    return float(res.fun), c


def _min_variance_within_width(
    sigma: np.ndarray, pre: np.ndarray, ell_full: np.ndarray, w0: np.ndarray, null: np.ndarray,
    g: np.ndarray, h: np.ndarray, width: float, cuts: list[np.ndarray], start: np.ndarray,
) -> np.ndarray | None:
    """Minimum-variance admissible weights whose bias range is at most ``width`` wide.

    Cutting planes on the free coordinates ``z``: whenever the current weights
    violate the width budget, the two extreme violations ``d_hi`` and ``d_lo``
    give a linear constraint ``c'(d_hi - d_lo) <= width`` that every feasible
    weight vector satisfies. ``cuts`` is shared across calls because a cut depends
    only on the polytope, not on the budget. The quadratic subproblem has only a
    handful of variables and no dual variables, so it stays well conditioned
    whether or not the budget binds.

    Returns ``None`` if no admissible weights achieve the width.
    """
    n_free = null.shape[1]
    pre_idx = np.flatnonzero(pre)
    base = ell_full.copy()
    base[pre] = w0

    def full(z: np.ndarray) -> np.ndarray:
        c = base.copy()
        c[pre] = w0 + null @ z
        return c

    if n_free == 0:
        b_lo, b_hi, _, _ = _bias_extremes(base, g, h, h.max(initial=0.0))
        return base if b_hi - b_lo <= width * (1 + 1e-9) + 1e-14 else None

    q = null.T @ sigma[np.ix_(pre_idx, pre_idx)] @ null
    lin = null.T @ (sigma @ base)[pre_idx]
    z = start.copy()
    for _ in range(200):
        c = full(z)
        b_lo, b_hi, d_lo, d_hi = _bias_extremes(c, g, h, h.max(initial=0.0))
        if b_hi - b_lo <= width * (1 + 1e-9) + 1e-14:
            return c
        cuts.append(d_hi - d_lo)
        rows = np.array([null.T @ v[pre_idx] for v in cuts])
        rhs = np.array([width - float(base @ v) for v in cuts])
        res = optimize.minimize(
            lambda zz: float(zz @ q @ zz + 2 * lin @ zz), z,
            jac=lambda zz: 2 * (q @ zz + lin), method="SLSQP",
            constraints=[{"type": "ineq",
                          "fun": lambda zz, rows=rows, rhs=rhs: rhs - rows @ zz,
                          "jac": lambda zz, rows=rows: -rows}],
            options={"maxiter": 300, "ftol": 1e-15},
        )
        if not res.success and res.status != 9:
            return None
        z = res.x
    return None


def _shortest_interval(
    result: EventStudyResult, m: float, restriction: str, target: Sequence[float] | None,
    alpha: float, reference: int,
) -> tuple[np.ndarray, float, float, float]:
    """Weights minimising interval length under a restriction, with their bias range.

    The interval [est - b_hi - x*sd, est - b_lo + x*sd], with x from
    :func:`_tight_cv`, has coverage at least 1 - alpha. Its length depends on the
    bias-range width W and the standard deviation. For each W on a grid from
    slack to tight, the minimum-variance weights come from
    :func:`_min_variance_within_width`, warm started from the previous solution;
    a bounded search then refines around the best, and the best candidate seen
    anywhere is kept.

    Returns ``(weights, sd, bias_low, bias_high)``.
    """
    rel = result.rel_periods.astype(int)
    pre, post = rel < 0, rel >= 0
    if not pre.any():
        raise ValueError("Need pre-period coefficients to remove trends.")
    ell = _target(rel, target)
    sigma = result.vcov
    g, h = _polyhedron(rel, reference, m, restriction)
    ell_full = np.zeros(len(rel))
    ell_full[post] = ell
    w0, null, w_star = _trend_removing_basis(rel, pre, post, ell, sigma, reference)
    c_star = ell_full.copy()
    c_star[pre] = w_star

    def length_of(c: np.ndarray) -> tuple[float, float, float, float]:
        sd = float(np.sqrt(c @ sigma @ c))
        b_lo, b_hi, _, _ = _bias_extremes(c, g, h, h.max(initial=0.0))
        return (b_hi - b_lo) + 2 * sd * _tight_cv((b_hi - b_lo) / sd, alpha), sd, b_lo, b_hi

    best_c, best_len = c_star, length_of(c_star)[0]
    b_lo, b_hi, _, _ = _bias_extremes(c_star, g, h, h.max(initial=0.0))
    w_hi = b_hi - b_lo
    if m > 0 and w_hi > 1e-14 and null.shape[1]:
        w_min, c_min = _min_width(g, h, pre, ell_full)
        w_floor = w_min * (1 + 1e-7) + 1e-15   # strictly feasible, so the cutting planes can land
        cuts: list[np.ndarray] = []

        def to_z(c: np.ndarray) -> np.ndarray:
            return np.linalg.lstsq(null, c[pre] - w0, rcond=None)[0]

        def solve(wdt: float, start: np.ndarray) -> np.ndarray | None:
            return _min_variance_within_width(sigma, pre, ell_full, w0, null, g, h, wdt, cuts, start)

        widths = []
        if w_hi > w_floor:
            # Denser near the floor, where the optimum usually sits.
            widths = list(w_floor + (w_hi - w_floor) * np.linspace(0.0, 1.0, 21) ** 2)
        z_prev = to_z(c_star)
        evaluated: list[tuple[float, float, np.ndarray]] = []
        for wdt in sorted(widths, reverse=True):
            c = solve(wdt, z_prev)
            if c is None:
                c = solve(wdt, to_z(c_min))
            if c is None:
                continue
            z_prev = to_z(c)
            length = length_of(c)[0]
            evaluated.append((wdt, length, c))
            if length < best_len:
                best_c, best_len = c, length

        if len(evaluated) >= 2:
            evaluated.sort(key=lambda e: e[0])
            j = int(np.argmin([e[1] for e in evaluated]))
            lo = evaluated[max(j - 1, 0)][0]
            hi = evaluated[min(j + 1, len(evaluated) - 1)][0]
            start = to_z(evaluated[j][2])

            def objective(wdt: float) -> float:
                c = solve(wdt, start)
                return np.inf if c is None else length_of(c)[0]

            if hi > lo:
                res = optimize.minimize_scalar(objective, bounds=(lo, hi), method="bounded",
                                               options={"xatol": (hi - lo) * 1e-4})
                if np.isfinite(res.fun) and res.fun < best_len:
                    c = solve(res.x, start)
                    if c is not None and length_of(c)[0] < best_len:
                        best_c, best_len = c, length_of(c)[0]
        # Polish: the length is a convex function of the free weights (a perspective of
        # a convex critical-value function, composed with convex bias and standard
        # deviation), so a short local search from the best candidate can only move
        # toward the global optimum. It closes the residual gap left by the
        # one-dimensional search over widths.
        z0 = to_z(best_c)

        def length_z(z: np.ndarray) -> float:
            c = ell_full.copy()
            c[pre] = w0 + null @ z
            return length_of(c)[0]

        scale = max(float(np.abs(z0).max()), 1e-3) * 0.02
        simplex = np.vstack([z0] + [z0 + scale * e for e in np.eye(len(z0))])
        res = optimize.minimize(length_z, z0, method="Nelder-Mead",
                                options={"initial_simplex": simplex, "xatol": 1e-10, "fatol": 1e-13,
                                         "maxiter": 400})
        if res.fun < best_len:
            best_c = ell_full.copy()
            best_c[pre] = w0 + null @ res.x
            best_len = res.fun
    _, sd, b_lo, b_hi = length_of(best_c)
    return best_c, sd, b_lo, b_hi


# ---------------------------------------------------------------------------
# Public interface
# ---------------------------------------------------------------------------


def flci(
    result: EventStudyResult,
    m: float,
    target: Sequence[float] | None = None,
    alpha: float = 0.05,
    reference: int = -1,
) -> FLCI:
    """Fixed-length confidence interval for a post-period effect under Delta^SD(M).

    Two earlier optimisers failed here, and the tests now guard against both.
    Nelder-Mead stalled at local points, producing an interval that narrowed as M
    grew. An SLSQP solve over dual variables broke down whenever the bias budget
    stopped binding, because the duals then have no unique value: on one real
    event study it failed at 17 of 24 budgets and missed the optimum by 1.6%.

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
    c, sd, b_lo, b_hi = _shortest_interval(result, m, "smooth", target, alpha, reference)
    bias = max(b_hi, -b_lo, 0.0)
    est = float(c @ result.coefs)
    hl = sd * _folded_normal_cv(bias / sd, alpha)
    return FLCI(m=m, estimate=est, lower=est - hl, upper=est + hl, sd=sd, max_bias=bias, weights=c)


def breakdown_m(
    result: EventStudyResult,
    target: Sequence[float] | None = None,
    alpha: float = 0.05,
    reference: int = -1,
    tol: float = 1e-6,
) -> float:
    """Largest M at which the smoothness FLCI still excludes zero.

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


def restricted_ci(
    result: EventStudyResult,
    m: float,
    restriction: str = "smooth",
    target: Sequence[float] | None = None,
    alpha: float = 0.05,
    reference: int = -1,
) -> RestrictedCI:
    """Confidence interval for a post-period effect under a shape or sign restriction.

    One-sided restrictions shift the interval rather than only widening it: under
    ``concave``, for instance, the counterfactual trend can only fall below the
    straight line, so the lower limit stops depending on M. For the symmetric
    ``smooth`` restriction the interval coincides with :func:`flci`.

    Args:
        result: Event study with its full covariance.
        m: Smoothness bound, per period squared.
        restriction: One of :data:`RESTRICTIONS`.
        target: Weights on the post-period coefficients; defaults to their average.
        alpha: One minus the coverage level.
        reference: The omitted, normalised relative period.

    Returns:
        The interval, its centre, the allowed bias range and the weights.
    """
    c, sd, b_lo, b_hi = _shortest_interval(result, m, restriction, target, alpha, reference)
    est = float(c @ result.coefs)
    x = _tight_cv((b_hi - b_lo) / sd, alpha)
    return RestrictedCI(restriction=restriction, m=m, estimate=est, lower=est - b_hi - x * sd,
                        upper=est - b_lo + x * sd, sd=sd, bias_low=b_lo, bias_high=b_hi, weights=c)


def restricted_breakdown(
    result: EventStudyResult,
    restriction: str = "smooth",
    target: Sequence[float] | None = None,
    alpha: float = 0.05,
    reference: int = -1,
    m_max: float = 1.0,
    tol: float = 1e-6,
) -> float:
    """Largest M at which the restricted interval still excludes zero.

    Returns ``inf`` when it excludes zero for every M up to ``m_max``: under that
    restriction, no amount of bending allowed by it overturns the conclusion.
    """
    if not restricted_ci(result, 0.0, restriction, target, alpha, reference).excludes_zero:
        return 0.0
    if restricted_ci(result, m_max, restriction, target, alpha, reference).excludes_zero:
        return np.inf
    lo, hi = 0.0, m_max
    while hi - lo > tol * max(1.0, hi):
        mid = (lo + hi) / 2
        if restricted_ci(result, mid, restriction, target, alpha, reference).excludes_zero:
            lo = mid
        else:
            hi = mid
    return lo


def pre_period_support(result: EventStudyResult, restriction: str, reference: int = -1) -> dict[str, object]:
    """Whether the estimated leads support a shape restriction.

    A restriction the pre-period data contradict should not be imposed on the
    post period, and one they only weakly favour should not be leaned on. Each
    relevant pre-period difference is tested one-sided at 5%. Verdicts:

    ``contradicted``               some difference is significantly of the wrong sign.
    ``supported``                  some difference is significantly of the right sign,
                                   and at least half the point estimates agree.
    ``mixed``                      some difference is significantly of the right sign,
                                   but most point estimates have the wrong sign.
    ``consistent, not supported``  nothing is significant either way.

    ``mixed`` exists because a single significant value from the noisiest lead,
    with the precisely estimated ones pointing the other way, is not support.
    """
    rel = result.rel_periods.astype(int)
    pre = rel < 0
    if restriction in ("concave", "convex"):
        mat, want = _second_differences(rel, reference), (-1.0 if restriction == "concave" else 1.0)
    elif restriction in ("increasing", "decreasing"):
        mat, want = _first_differences(rel, reference), (1.0 if restriction == "increasing" else -1.0)
    else:
        return {"restriction": restriction, "verdict": "not testable from leads", "values": []}
    rows = [r for r in mat if not np.any(r[~pre])]  # built only from leads and the reference
    vals = []
    for r in rows:
        v, s = float(r @ result.coefs), float(np.sqrt(r @ result.vcov @ r))
        vals.append({"value": v, "se": s, "z": v / s if s > 0 else 0.0})
    crit = stats.norm.ppf(0.95)
    wrong = any(want * x["z"] < -crit for x in vals)
    right = any(want * x["z"] > crit for x in vals)
    agree = sum(want * x["value"] > 0 for x in vals)
    if wrong:
        verdict = "contradicted"
    elif right and agree * 2 >= len(vals):
        verdict = "supported"
    elif right:
        verdict = "mixed"
    else:
        verdict = "consistent, not supported"
    return {"restriction": restriction, "verdict": verdict, "values": vals}


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
