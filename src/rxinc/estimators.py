"""Estimators for the effect of payment relationships on prescribing.

Five designs, ordered by how much confounding they can survive:

:func:`naive_ols`
    Pooled cross-section.  Compares prescribing of paid vs unpaid physicians.
    Absorbs every persistent difference between those groups into the
    coefficient.  This is the estimator most press coverage of Open Payments
    implicitly reports.

:func:`twoway_fe`
    Physician and period fixed effects.  Differences out time-invariant
    confounding, so it survives static targeting.  In a staggered-adoption
    design it still uses already-treated physicians as comparisons, which is
    the Goodman-Bacon problem.

:func:`did_never_treated`
    Cohort-by-period ATT against never-treated physicians only, aggregated with
    cohort weights, in the spirit of Callaway and Sant'Anna.  Avoids
    already-treated comparisons.

:func:`event_study`
    Relative-time coefficients around payment onset.  Its value is not the
    post-period estimate but the *pre*-period ones: they are the only thing in
    this file that can detect dynamic targeting.

:func:`interrupted_time_series`
    Within-physician level and slope change at onset, the design the 2021
    systematic review explicitly called for.

All standard errors are clustered on physician, because the AR(1) shock makes
observations within a physician dependent and unclustered errors badly
overstate precision.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from rxinc.simulate import NEVER_TREATED


@dataclass
class Estimate:
    """A single scalar treatment-effect estimate.

    Attributes:
        name: Human-readable estimator name.
        coef: Point estimate, in log points.
        se: Cluster-robust standard error.
        n_obs: Observations used.
        n_clusters: Physicians used.
        true_effect: Ground truth when known (simulation); ``None`` on real
            data, where bias is not computable.
    """

    name: str
    coef: float
    se: float
    n_obs: int
    n_clusters: int
    true_effect: float | None = None

    @property
    def t_stat(self) -> float:
        return self.coef / self.se if self.se > 0 else np.nan

    @property
    def ci95(self) -> tuple[float, float]:
        return (self.coef - 1.96 * self.se, self.coef + 1.96 * self.se)

    @property
    def bias(self) -> float | None:
        """Point estimate minus truth, or ``None`` if truth is unknown."""
        if self.true_effect is None:
            return None
        return self.coef - self.true_effect

    @property
    def bias_pct(self) -> float | None:
        """Bias as a percentage of the true effect."""
        if self.true_effect is None or self.true_effect == 0:
            return None
        return 100.0 * (self.coef - self.true_effect) / self.true_effect

    @property
    def covers_truth(self) -> bool | None:
        """Whether the 95% interval contains the true effect.

        Returns ``None`` when truth is unknown or the standard error was not
        computed, so that averaging coverage over replications yields ``NaN``
        rather than silently reporting 0%.
        """
        if self.true_effect is None or not np.isfinite(self.se):
            return None
        low, high = self.ci95
        return bool(low <= self.true_effect <= high)

    def __str__(self) -> str:
        low, high = self.ci95
        core = (
            f"{self.name:<28s} {self.coef: .4f}  (se {self.se:.4f})  "
            f"95% CI [{low: .4f}, {high: .4f}]"
        )
        if self.true_effect is None:
            return core
        if self.bias_pct is None:
            # True effect is exactly zero (e.g. a placebo), so relative bias
            # is undefined and only the absolute deviation is meaningful.
            return f"{core}  bias {self.bias:+.4f}"
        return f"{core}  bias {self.bias:+.4f} ({self.bias_pct:+.0f}%)"


@dataclass
class EventStudyResult:
    """Relative-time coefficients around payment onset.

    Attributes:
        rel_periods: Relative period for each coefficient.  The reference
            period (-1) is omitted and carries an implicit zero.
        coefs: Coefficient at each relative period.
        ses: Cluster-robust standard errors.
        vcov: Full covariance matrix, retained so joint pre-trend tests can be
            run without re-fitting.
        n_obs: Observations used.
        n_clusters: Physicians used.
        true_effect: Ground truth when known.
    """

    rel_periods: np.ndarray
    coefs: np.ndarray
    ses: np.ndarray
    vcov: np.ndarray
    n_obs: int
    n_clusters: int
    true_effect: float | None = None

    def to_frame(self) -> pd.DataFrame:
        """Return coefficients as a tidy DataFrame with 95% bounds."""
        return pd.DataFrame(
            {
                "rel_period": self.rel_periods,
                "coef": self.coefs,
                "se": self.ses,
                "ci_low": self.coefs - 1.96 * self.ses,
                "ci_high": self.coefs + 1.96 * self.ses,
                "is_pre": self.rel_periods < 0,
            }
        )

    @property
    def pre_mask(self) -> np.ndarray:
        """Boolean mask selecting pre-treatment (lead) coefficients."""
        return self.rel_periods < 0

    def post_average(self) -> float:
        """Simple mean of post-onset coefficients."""
        post = self.rel_periods >= 0
        return float(self.coefs[post].mean()) if post.any() else np.nan


def _ols(
    y: np.ndarray,
    X: np.ndarray,
    cluster: np.ndarray | None = None,
    absorbed: int = 0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Least squares with optional cluster-robust covariance.

    Args:
        y: Outcome vector, shape ``(n,)``.
        X: Design matrix, shape ``(n, k)``.
        cluster: Cluster id per observation.  ``None`` gives classical errors.
        absorbed: Number of parameters removed by a prior within transform.
            Counted against the residual degrees of freedom so the finite-sample
            correction is right.

    Returns:
        ``(beta, vcov, resid)``.
    """
    n, k = X.shape
    xtx_inv = np.linalg.pinv(X.T @ X)
    beta = xtx_inv @ (X.T @ y)
    resid = y - X @ beta
    k_eff = k + absorbed

    if cluster is None:
        dof = max(n - k_eff, 1)
        vcov = (resid @ resid / dof) * xtx_inv
        return beta, vcov, resid

    codes, _ = pd.factorize(cluster)
    n_groups = codes.max() + 1
    scores = np.zeros((n_groups, k))
    np.add.at(scores, codes, X * resid[:, None])
    meat = scores.T @ scores

    correction = (n_groups / max(n_groups - 1, 1)) * (
        (n - 1) / max(n - k_eff, 1)
    )
    vcov = correction * (xtx_inv @ meat @ xtx_inv)
    return beta, vcov, resid


def _within_transform(
    frame: pd.DataFrame,
    columns: list[str],
    unit: str,
    time: str,
    tol: float = 1e-10,
    max_iter: int = 500,
) -> np.ndarray:
    """Two-way demean ``columns`` by ``unit`` and ``time``.

    Implemented by alternating projections: demean by unit, demean by period,
    repeat until both sets of group means are within ``tol`` of zero.  A single
    pass is exact only on a balanced panel; iterating makes the transform exact
    on unbalanced panels too, which matters because placebo and subgroup
    analyses routinely drop observations. The result is equivalent to including
    full unit and period dummies (Frisch-Waugh-Lovell) at far less memory.

    Args:
        frame: Panel containing ``columns``, ``unit`` and ``time``.
        columns: Columns to transform.
        unit: Unit identifier column.
        time: Period identifier column.
        tol: Convergence tolerance on the largest remaining group mean.
        max_iter: Maximum alternating-projection sweeps.

    Returns:
        Array of shape ``(len(frame), len(columns))`` holding demeaned values.
    """
    values = frame[columns].to_numpy(dtype=float).copy()
    unit_codes, _ = pd.factorize(frame[unit])
    time_codes, _ = pd.factorize(frame[time])
    n_units = unit_codes.max() + 1
    n_times = time_codes.max() + 1

    unit_counts = np.bincount(unit_codes, minlength=n_units).astype(float)
    time_counts = np.bincount(time_codes, minlength=n_times).astype(float)

    values -= values.mean(axis=0, keepdims=True)

    for _ in range(max_iter):
        worst = 0.0
        for j in range(values.shape[1]):
            col = values[:, j]

            unit_mean = np.bincount(unit_codes, weights=col, minlength=n_units)
            unit_mean /= unit_counts
            col = col - unit_mean[unit_codes]

            time_mean = np.bincount(time_codes, weights=col, minlength=n_times)
            time_mean /= time_counts
            col = col - time_mean[time_codes]

            values[:, j] = col
            worst = max(worst, float(np.abs(time_mean).max()))

        if worst < tol:
            break

    return values


def _true_effect(frame: pd.DataFrame) -> float | None:
    return frame.attrs.get("true_effect")


def naive_ols(
    frame: pd.DataFrame, outcome: str = "log_rx", treatment: str = "treated"
) -> Estimate:
    """Pooled cross-sectional OLS of outcome on treatment status.

    No physician controls, so any persistent difference between paid and unpaid
    physicians -- specialty, panel size, practice setting, baseline prescribing
    volume -- is loaded onto the coefficient.

    Args:
        frame: Physician-period panel.
        outcome: Outcome column, in logs.
        treatment: Binary treatment column.

    Returns:
        The estimate, with cluster-robust standard errors.
    """
    y = frame[outcome].to_numpy(dtype=float)
    d = frame[treatment].to_numpy(dtype=float)
    X = np.column_stack([np.ones_like(d), d])
    beta, vcov, _ = _ols(y, X, cluster=frame["physician_id"].to_numpy())
    return Estimate(
        name="Naive pooled OLS",
        coef=float(beta[1]),
        se=float(np.sqrt(vcov[1, 1])),
        n_obs=len(frame),
        n_clusters=frame["physician_id"].nunique(),
        true_effect=_true_effect(frame),
    )


def twoway_fe(
    frame: pd.DataFrame, outcome: str = "log_rx", treatment: str = "treated"
) -> Estimate:
    """Two-way fixed effects: physician and period.

    Removes time-invariant confounding.  Two caveats worth stating plainly:
    with staggered adoption it draws comparisons from already-treated
    physicians, and it cannot touch confounders that move over time -- which is
    exactly what dynamic targeting is.

    Args:
        frame: Physician-period panel.
        outcome: Outcome column, in logs.
        treatment: Binary treatment column.

    Returns:
        The estimate, with cluster-robust standard errors.
    """
    demeaned = _within_transform(
        frame, [outcome, treatment], unit="physician_id", time="period"
    )
    y, d = demeaned[:, 0], demeaned[:, 1:]
    n_units = frame["physician_id"].nunique()
    n_times = frame["period"].nunique()
    beta, vcov, _ = _ols(
        y,
        d,
        cluster=frame["physician_id"].to_numpy(),
        absorbed=n_units + n_times - 1,
    )
    return Estimate(
        name="Two-way fixed effects",
        coef=float(beta[0]),
        se=float(np.sqrt(vcov[0, 0])),
        n_obs=len(frame),
        n_clusters=n_units,
        true_effect=_true_effect(frame),
    )


def event_study(
    frame: pd.DataFrame,
    leads: int = 4,
    lags: int = 6,
    outcome: str = "log_rx",
) -> EventStudyResult:
    """Fit relative-time coefficients around payment onset.

    Relative period -1 is the omitted reference.  Relative periods beyond the
    requested window are binned into the endpoints so that no observation is
    silently dropped.  Never-treated physicians enter with all relative-time
    dummies at zero and identify the period effects.

    The pre-onset coefficients are the payload.  Under static targeting they
    should sit flat at zero; under dynamic targeting they trend upward, which
    is the visible signature that the treated group was already rising.

    Args:
        frame: Physician-period panel with an ``event_time`` column.
        leads: Pre-onset periods to estimate (excluding the reference).
        lags: Post-onset periods to estimate.
        outcome: Outcome column, in logs.

    Returns:
        Coefficients, standard errors and covariance by relative period.
    """
    work = frame.copy()
    rel = work["event_time"].to_numpy(dtype=float)
    is_treated_unit = ~np.isnan(rel)
    rel_clipped = np.clip(rel, -leads, lags)

    rel_grid = [r for r in range(-leads, lags + 1) if r != -1]
    dummies = np.zeros((len(work), len(rel_grid)))
    for j, r in enumerate(rel_grid):
        dummies[:, j] = (is_treated_unit & (rel_clipped == r)).astype(float)

    columns = [f"_rel_{r}" for r in rel_grid]
    for j, name in enumerate(columns):
        work[name] = dummies[:, j]

    demeaned = _within_transform(
        work, [outcome, *columns], unit="physician_id", time="period"
    )
    y, X = demeaned[:, 0], demeaned[:, 1:]
    n_units = work["physician_id"].nunique()
    n_times = work["period"].nunique()
    beta, vcov, _ = _ols(
        y,
        X,
        cluster=work["physician_id"].to_numpy(),
        absorbed=n_units + n_times - 1,
    )
    return EventStudyResult(
        rel_periods=np.array(rel_grid),
        coefs=beta,
        ses=np.sqrt(np.diag(vcov)),
        vcov=vcov,
        n_obs=len(work),
        n_clusters=n_units,
        true_effect=_true_effect(frame),
    )


def did_never_treated(
    frame: pd.DataFrame,
    outcome: str = "log_rx",
    n_boot: int = 200,
    seed: int = 20260928,
) -> Estimate:
    """Cohort-by-period ATT using never-treated physicians as the comparison.

    For each onset cohort ``g`` and post period ``t >= g``, the estimate is the
    change in outcome from ``g - 1`` to ``t`` among cohort ``g``, minus the same
    change among never-treated physicians.  Cohort-period effects are averaged
    with weights proportional to cohort size, giving an overall ATT that never
    uses an already-treated physician as a control.

    Standard errors come from a physician-level bootstrap, which respects the
    within-physician dependence induced by the AR(1) shock.

    Args:
        frame: Physician-period panel.
        outcome: Outcome column, in logs.
        n_boot: Bootstrap replications.  Set to 0 to skip and return ``se=nan``.
        seed: Bootstrap seed.

    Returns:
        The aggregated ATT.

    Raises:
        ValueError: If the panel contains no never-treated physicians, in which
            case this design is not available.
    """
    wide = frame.pivot(index="physician_id", columns="period", values=outcome)
    cohort = frame.groupby("physician_id")["first_treat_period"].first()
    cohort = cohort.reindex(wide.index)

    never = (cohort == NEVER_TREATED).to_numpy()
    if not never.any():
        raise ValueError(
            "did_never_treated requires never-treated physicians; none present. "
            "Lower onset_intercept or use event_study instead."
        )

    values = wide.to_numpy(dtype=float)
    periods = wide.columns.to_numpy()
    cohort_values = cohort.to_numpy()
    treat_cohorts = sorted({int(g) for g in cohort_values if g != NEVER_TREATED})

    def _att(rows: np.ndarray) -> float:
        vals = values[rows]
        coh = cohort_values[rows]
        is_never = coh == NEVER_TREATED
        if not is_never.any():
            return np.nan
        total, weight = 0.0, 0.0
        for g in treat_cohorts:
            in_cohort = coh == g
            n_g = int(in_cohort.sum())
            if n_g == 0 or g - 1 < periods.min():
                continue
            base_idx = int(np.where(periods == g - 1)[0][0])
            for t in periods[periods >= g]:
                t_idx = int(np.where(periods == t)[0][0])
                treated_change = float(
                    np.mean(vals[in_cohort, t_idx] - vals[in_cohort, base_idx])
                )
                control_change = float(
                    np.mean(vals[is_never, t_idx] - vals[is_never, base_idx])
                )
                total += n_g * (treated_change - control_change)
                weight += n_g
        return total / weight if weight else np.nan

    all_rows = np.arange(values.shape[0])
    point = _att(all_rows)

    se = np.nan
    if n_boot > 0:
        rng = np.random.default_rng(seed)
        draws = np.empty(n_boot)
        for b in range(n_boot):
            sample = rng.integers(0, values.shape[0], size=values.shape[0])
            draws[b] = _att(sample)
        se = float(np.nanstd(draws, ddof=1))

    return Estimate(
        name="DiD vs never-treated",
        coef=float(point),
        se=se,
        n_obs=len(frame),
        n_clusters=int(values.shape[0]),
        true_effect=_true_effect(frame),
    )


def interrupted_time_series(
    frame: pd.DataFrame, outcome: str = "log_rx"
) -> tuple[Estimate, Estimate]:
    """Within-physician level and slope change at payment onset.

    Restricted to eventually-treated physicians, this fits, net of physician
    fixed effects, a pre-onset trend, a level shift at onset, and a change in
    slope afterwards.  It is the interrupted time-series design the 2021
    systematic review named as a way past pure association.

    Its weakness is the mirror of its strength: the pre-onset trend is estimated
    from the treated units themselves, so under dynamic targeting the level
    shift absorbs part of a pre-existing rise.

    Args:
        frame: Physician-period panel.
        outcome: Outcome column, in logs.

    Returns:
        ``(level_change, slope_change)`` estimates.

    Raises:
        ValueError: If no eventually-treated physicians are present.
    """
    treated_units = frame[frame["ever_treated"]].copy()
    if treated_units.empty:
        raise ValueError("interrupted_time_series requires treated physicians.")

    rel = treated_units["event_time"].to_numpy(dtype=float)
    post = (rel >= 0).astype(float)
    treated_units["_time"] = rel
    treated_units["_post"] = post
    treated_units["_post_time"] = post * rel

    columns = ["_time", "_post", "_post_time"]
    demeaned = _within_transform(
        treated_units, [outcome, *columns], unit="physician_id", time="period"
    )
    y, X = demeaned[:, 0], demeaned[:, 1:]
    n_units = treated_units["physician_id"].nunique()
    n_times = treated_units["period"].nunique()
    beta, vcov, _ = _ols(
        y,
        X,
        cluster=treated_units["physician_id"].to_numpy(),
        absorbed=n_units + n_times - 1,
    )

    truth = _true_effect(frame)
    level = Estimate(
        name="ITS level change",
        coef=float(beta[1]),
        se=float(np.sqrt(vcov[1, 1])),
        n_obs=len(treated_units),
        n_clusters=n_units,
        true_effect=truth,
    )
    slope = Estimate(
        name="ITS slope change",
        coef=float(beta[2]),
        se=float(np.sqrt(vcov[2, 2])),
        n_obs=len(treated_units),
        n_clusters=n_units,
        true_effect=None,
    )
    return level, slope

def triple_diff(
    frame: pd.DataFrame,
    outcome: str = "log_rx",
    treatment: str = "treated",
    unit: str = "physician_year",
    product: str = "drug_year",
    cluster_on: str = "physician_id",
) -> Estimate:
    """Physician-by-year and drug-by-year fixed effects on a drug-level panel.

    Fits

        log(claims_ijt) = tau * Paid_ijt + alpha_it + delta_jt + e_ijt

    where ``alpha_it`` absorbs everything about a physician in a year, including
    their whole prescribing trajectory, and ``delta_jt`` absorbs national
    drug-level shocks such as launches, guideline changes and competitor entry.

    ``tau`` is therefore identified only from variation within a physician-year
    across drugs: does this physician prescribe more of the drug they were paid
    about than of competing drugs they were not paid about, in the same year.

    That makes the estimator robust to selection on physician momentum, which
    defeats :func:`twoway_fe` and every other physician-level design here. It is
    *not* robust to selection on drug-specific momentum, because that varies
    within the physician-year. Use :func:`rxinc.diagnostics.pretrend_test` on a
    drug-level event study to tell the two apart.

    Args:
        frame: Physician-drug-year panel.
        outcome: Outcome column, in logs.
        treatment: Binary or continuous exposure column.
        unit: Composite physician-by-period key to absorb.
        product: Composite drug-by-period key to absorb.
        cluster_on: Column to cluster standard errors on. Physician is the right
            level, since a physician's drugs share their unobserved shocks.

    Returns:
        The estimate, with cluster-robust standard errors.

    Raises:
        ValueError: If required columns are absent, or the treatment has no
            variation left after absorbing both fixed effects, which means the
            design is not identified on this panel.
    """
    required = {outcome, treatment, unit, product, cluster_on}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"triple_diff needs columns {missing}")

    demeaned = _within_transform(frame, [outcome, treatment], unit=unit, time=product)
    y, d = demeaned[:, 0], demeaned[:, 1:]

    residual_variation = float(np.abs(d).max())
    if residual_variation < 1e-10:
        raise ValueError(
            "Treatment has no variation within physician-year across drugs; "
            "tau is not identified. Every physician was paid about all or none "
            "of their drugs."
        )

    n_unit = frame[unit].nunique()
    n_product = frame[product].nunique()
    beta, vcov, _ = _ols(
        y,
        d,
        cluster=frame[cluster_on].to_numpy(),
        absorbed=n_unit + n_product - 1,
    )
    return Estimate(
        name="Triple diff (physician-yr + drug-yr)",
        coef=float(beta[0]),
        se=float(np.sqrt(vcov[0, 0])),
        n_obs=len(frame),
        n_clusters=int(frame[cluster_on].nunique()),
        true_effect=_true_effect(frame),
    )


def drug_event_study(
    frame: pd.DataFrame,
    leads: int = 3,
    lags: int = 4,
    outcome: str = "log_rx",
    unit: str = "physician_year",
    product: str = "drug_year",
    cluster_on: str = "physician_id",
) -> EventStudyResult:
    """Event study around payment onset, inside the triple-difference design.

    Relative period -1 is omitted and endpoints are binned, as in
    :func:`event_study`, but the absorbed effects are physician-by-year and
    drug-by-year rather than physician and period.

    The leads are the diagnostic that matters. Under selection on physician
    momentum they should be flat, because that momentum is absorbed. Under
    selection on drug-specific momentum they will trend, revealing that the
    within-physician comparison is contaminated too.

    Args:
        frame: Physician-drug-year panel with ``event_time``.
        leads: Pre-onset periods to estimate, excluding the reference.
        lags: Post-onset periods to estimate.
        outcome: Outcome column, in logs.
        unit: Composite physician-by-period key to absorb.
        product: Composite drug-by-period key to absorb.
        cluster_on: Clustering column.

    Returns:
        Coefficients, standard errors and covariance by relative period.
    """
    work = frame.copy()
    rel = work["event_time"].to_numpy(dtype=float)
    is_treated_pair = ~np.isnan(rel)
    rel_clipped = np.clip(rel, -leads, lags)

    rel_grid = [r for r in range(-leads, lags + 1) if r != -1]
    columns = []
    for r in rel_grid:
        name = f"_rel_{r}"
        work[name] = (is_treated_pair & (rel_clipped == r)).astype(float)
        columns.append(name)

    demeaned = _within_transform(work, [outcome, *columns], unit=unit, time=product)
    y, X = demeaned[:, 0], demeaned[:, 1:]
    beta, vcov, _ = _ols(
        y,
        X,
        cluster=work[cluster_on].to_numpy(),
        absorbed=work[unit].nunique() + work[product].nunique() - 1,
    )
    return EventStudyResult(
        rel_periods=np.array(rel_grid),
        coefs=beta,
        ses=np.sqrt(np.diag(vcov)),
        vcov=vcov,
        n_obs=len(work),
        n_clusters=int(work[cluster_on].nunique()),
        true_effect=_true_effect(frame),
    )
