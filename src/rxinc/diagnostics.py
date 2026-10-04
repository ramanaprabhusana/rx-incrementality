"""Tests that tell you whether an estimate can be believed.

The estimators in :mod:`rxinc.estimators` all return a number.  Under dynamic
targeting most of those numbers are wrong, and nothing about the number itself
says so.  This module holds the checks that do:

:func:`pretrend_test`
    Joint Wald test that pre-onset event-study coefficients are zero.  This is
    the workhorse.  A rejection means the treated physicians were already
    diverging before any payment arrived, so the design is not identified.

:func:`placebo_shift`
    Re-estimate with payment onset moved earlier and all genuinely
    post-treatment periods discarded.  Any measured "effect" is pure
    confounding, since no treatment has occurred yet.

:func:`balance_table`
    Pre-period levels and growth rates by eventual treatment status --
    the descriptive version of the same question.

:func:`monte_carlo`
    Repeat an estimator across simulated panels to measure its bias and
    confidence-interval coverage, rather than inferring either from one draw.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass, replace

import numpy as np
import pandas as pd
from scipy import stats

from rxinc.estimators import Estimate, EventStudyResult, event_study, twoway_fe
from rxinc.simulate import NEVER_TREATED, PanelConfig, simulate_panel


@dataclass
class PreTrendTest:
    """Result of a joint test on pre-onset event-study coefficients.

    Attributes:
        statistic: Wald statistic.
        df: Degrees of freedom (number of lead coefficients tested).
        p_value: Upper-tail probability under the chi-squared null.
        max_abs_coef: Largest absolute lead coefficient, as an effect-size
            companion to the p-value.
        slope: OLS slope of lead coefficients on relative period.  A positive
            slope means the treated group was already rising into onset.
    """

    statistic: float
    df: int
    p_value: float
    max_abs_coef: float
    slope: float

    @property
    def passes(self) -> bool:
        """Whether parallel pre-trends survives at the 5% level."""
        return self.p_value >= 0.05

    def __str__(self) -> str:
        verdict = "PASS" if self.passes else "FAIL"
        return (
            f"Pre-trend test: {verdict}  chi2({self.df})={self.statistic:.1f}  "
            f"p={self.p_value:.2e}  max|lead|={self.max_abs_coef:.3f}  "
            f"slope={self.slope:+.4f}/period"
        )


def pretrend_test(result: EventStudyResult) -> PreTrendTest:
    """Jointly test that all pre-onset event-study coefficients are zero.

    Uses the full cluster-robust covariance retained on the event study, so the
    test accounts for correlation between lead coefficients rather than testing
    them one at a time.

    Args:
        result: A fitted event study.

    Returns:
        The test result.

    Raises:
        ValueError: If the event study has no pre-onset coefficients.
    """
    mask = result.pre_mask
    if not mask.any():
        raise ValueError("Event study has no lead coefficients to test.")

    b = result.coefs[mask]
    v = result.vcov[np.ix_(mask, mask)]
    statistic = float(b @ np.linalg.pinv(v) @ b)
    df = int(mask.sum())
    p_value = float(stats.chi2.sf(statistic, df))

    rel = result.rel_periods[mask].astype(float)
    slope = float(np.polyfit(rel, b, 1)[0]) if len(rel) > 1 else np.nan

    return PreTrendTest(
        statistic=statistic,
        df=df,
        p_value=p_value,
        max_abs_coef=float(np.abs(b).max()),
        slope=slope,
    )


@dataclass
class DetrendedEventStudy:
    """Event-study coefficients net of an extrapolated linear pre-trend.

    Attributes:
        rel_periods: Relative periods, reference omitted.
        coefs: Coefficients after subtracting the fitted trend.
        ses: Delta-method standard errors, accounting for estimation of the slope.
        slope: Fitted pre-trend per period, constrained through zero at the
            reference period.
        slope_se: Standard error of the slope.
        vcov: Covariance of the adjusted coefficients.
    """

    rel_periods: np.ndarray
    coefs: np.ndarray
    ses: np.ndarray
    slope: float
    slope_se: float
    vcov: np.ndarray

    def post_average(self) -> tuple[float, float]:
        """Mean adjusted post-onset coefficient and its standard error."""
        post = (self.rel_periods >= 0).astype(float)
        a = post / max(post.sum(), 1.0)
        return float(a @ self.coefs), float(np.sqrt(max(a @ self.vcov @ a, 0.0)))


def detrend_event_study(result: EventStudyResult, reference: int = -1) -> DetrendedEventStudy:
    """Subtract a linear pre-trend extrapolated through the post period.

    When pre-trends fail, post-period coefficients mix any effect with the
    continuation of whatever was already happening. This fits a line through the
    lead coefficients, constrained to pass through zero at the reference period
    (where the normalisation puts it), by generalised least squares using their
    full covariance, so a noisy early lead gets little weight. The fitted line is
    then subtracted from every coefficient.

    It is a sensitivity analysis, not identification: it assumes the
    pre-existing trend would have continued linearly. Rambachan and Roth (2023)
    formalise how far such assumptions can be relaxed.

    Args:
        result: A fitted event study with its covariance.
        reference: The omitted relative period.

    Returns:
        Detrended coefficients with delta-method standard errors.

    Raises:
        ValueError: With no lead coefficients, the trend cannot be fitted.
    """
    rel = result.rel_periods.astype(float)
    pre = result.rel_periods < 0
    if not pre.any():
        raise ValueError("No lead coefficients to fit a pre-trend on.")
    x = rel[pre] - reference
    v_pre = result.vcov[np.ix_(pre, pre)]
    w = np.linalg.pinv(v_pre)
    denom = float(x @ w @ x)
    # slope = g' b_pre with g = W x / (x' W x); write it as a linear map on all coefficients.
    g = np.zeros(len(rel))
    g[pre] = (w @ x) / denom
    slope = float(g @ result.coefs)
    # adjusted = (I - d g') b, where d_k = rel_k - reference
    d = rel - reference
    A = np.eye(len(rel)) - np.outer(d, g)
    adjusted = A @ result.coefs
    vcov = A @ result.vcov @ A.T
    return DetrendedEventStudy(
        rel_periods=result.rel_periods,
        coefs=adjusted,
        ses=np.sqrt(np.clip(np.diag(vcov), 0, None)),
        slope=slope,
        slope_se=float(np.sqrt(max(g @ result.vcov @ g, 0.0))),
        vcov=vcov,
    )


def placebo_shift(
    frame: pd.DataFrame, shift: int = 2, outcome: str = "log_rx"
) -> Estimate:
    """Estimate a treatment effect that cannot exist.

    Payment onset is moved ``shift`` periods earlier and every period from the
    real onset onward is dropped, so no observation in the resulting panel has
    actually been treated.  A two-way fixed effects estimate on this panel is
    therefore an estimate of zero plus whatever confounding remains.

    Args:
        frame: Physician-period panel.
        shift: Periods to move onset earlier.  Must be at least 1.
        outcome: Outcome column, in logs.

    Returns:
        The placebo estimate.  Interpret a significant coefficient as evidence
        that the design is picking up selection rather than treatment.

    Raises:
        ValueError: If ``shift`` is not positive, or if the shift leaves no
            usable pre-period observations.
    """
    if shift < 1:
        raise ValueError(f"shift must be >= 1, got {shift}")

    work = frame.copy()
    real_onset = work["first_treat_period"].to_numpy()
    treated_unit = real_onset != NEVER_TREATED

    # Discard anything from the real onset onward: those periods are genuinely
    # treated and would contaminate a placebo.
    keep = ~(treated_unit & (work["period"].to_numpy() >= real_onset))
    work = work.loc[keep].copy()
    if work.empty:
        raise ValueError("Placebo shift removed every observation.")

    pseudo_onset = np.where(
        work["first_treat_period"].to_numpy() == NEVER_TREATED,
        NEVER_TREATED,
        work["first_treat_period"].to_numpy() - shift,
    )
    work["first_treat_period"] = pseudo_onset
    work["treated"] = (pseudo_onset != NEVER_TREATED) & (
        work["period"].to_numpy() >= pseudo_onset
    )
    work["ever_treated"] = pseudo_onset != NEVER_TREATED
    work["event_time"] = np.where(
        pseudo_onset == NEVER_TREATED,
        np.nan,
        work["period"].to_numpy() - pseudo_onset,
    )

    if not work["treated"].any():
        raise ValueError(
            f"shift={shift} leaves no pseudo-treated observations; "
            "use a smaller shift or a longer burn-in."
        )

    work.attrs["true_effect"] = 0.0
    est = twoway_fe(work, outcome=outcome)
    est.name = f"Placebo (onset -{shift})"
    return est


def balance_table(frame: pd.DataFrame, outcome: str = "log_rx") -> pd.DataFrame:
    """Compare pre-onset levels and growth by eventual treatment status.

    Uses only periods strictly before each physician's own onset, so the
    comparison is untainted by treatment.  Never-treated physicians contribute
    their first periods, matched to the median onset among treated physicians.

    Args:
        frame: Physician-period panel.
        outcome: Outcome column, in logs.

    Returns:
        A two-row DataFrame indexed by ``ever_treated`` with the mean
        pre-period level, the mean within-physician pre-period growth, and the
        number of physicians in each group.
    """
    onset = frame["first_treat_period"].to_numpy()
    treated_unit = onset != NEVER_TREATED
    median_onset = (
        int(np.median(onset[treated_unit])) if treated_unit.any() else frame["period"].max()
    )
    cutoff = np.where(treated_unit, onset, median_onset)
    pre = frame.loc[frame["period"].to_numpy() < cutoff].copy()

    by_unit = pre.sort_values("period").groupby("physician_id")
    summary = pd.DataFrame(
        {
            "ever_treated": by_unit["ever_treated"].first(),
            "pre_level": by_unit[outcome].mean(),
            "pre_growth": by_unit[outcome].apply(
                lambda s: (s.iloc[-1] - s.iloc[0]) / max(len(s) - 1, 1)
            ),
        }
    )
    out = summary.groupby("ever_treated").agg(
        n_physicians=("pre_level", "size"),
        mean_pre_level=("pre_level", "mean"),
        mean_pre_growth=("pre_growth", "mean"),
    )
    return out


def monte_carlo(
    config: PanelConfig,
    estimators: Mapping[str, Callable[[pd.DataFrame], Estimate]] | None = None,
    n_reps: int = 100,
) -> pd.DataFrame:
    """Measure estimator bias and coverage across repeated simulated panels.

    A single simulated panel confuses bias with sampling noise.  This repeats
    the whole exercise on fresh seeds and reports the sampling distribution, so
    a claim like "two-way fixed effects is unbiased under static targeting" is
    backed by a mean over replications and a coverage rate rather than one
    lucky draw.

    Args:
        config: Base data generating process.  The seed is advanced by one per
            replication.
        estimators: Mapping of display name to a function taking the panel and
            returning an :class:`~rxinc.estimators.Estimate`.  Defaults to
            two-way fixed effects alone, which is the cheap one.
        n_reps: Number of replications.

    Returns:
        One row per estimator with ``mean_coef``, ``mean_bias``, ``sd_coef``,
        ``rmse``, ``coverage95`` and ``n_reps``.
    """
    funcs = dict(estimators) if estimators else {"Two-way fixed effects": twoway_fe}
    records: list[dict[str, object]] = []

    for rep in range(n_reps):
        panel = simulate_panel(replace(config, seed=config.seed + rep))
        truth = panel.attrs["true_effect"]
        for name, fn in funcs.items():
            est = fn(panel)
            records.append(
                {
                    "estimator": name,
                    "coef": est.coef,
                    "bias": est.coef - truth,
                    "covers": est.covers_truth,
                }
            )

    raw = pd.DataFrame.from_records(records)
    out = raw.groupby("estimator").agg(
        mean_coef=("coef", "mean"),
        mean_bias=("bias", "mean"),
        sd_coef=("coef", "std"),
        coverage95=("covers", "mean"),
        n_reps=("coef", "size"),
    )
    out["rmse"] = np.sqrt(
        raw.groupby("estimator")["bias"].apply(lambda s: float((s**2).mean()))
    )
    return out[
        ["mean_coef", "mean_bias", "sd_coef", "rmse", "coverage95", "n_reps"]
    ]


def diagnose(frame: pd.DataFrame, leads: int = 4, lags: int = 6) -> dict[str, object]:
    """Run the full diagnostic battery on a panel.

    Args:
        frame: Physician-period panel.
        leads: Pre-onset periods for the event study.
        lags: Post-onset periods for the event study.

    Returns:
        Mapping with keys ``event_study``, ``pretrend``, ``placebo`` and
        ``balance``.
    """
    es = event_study(frame, leads=leads, lags=lags)
    return {
        "event_study": es,
        "pretrend": pretrend_test(es),
        "placebo": placebo_shift(frame),
        "balance": balance_table(frame),
    }
