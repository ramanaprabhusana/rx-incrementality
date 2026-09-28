"""Simulated physician-quarter panels with a known causal effect.

The point of the simulation is to reproduce the confounding structure that the
observational literature flags but cannot rule out.  A 2021 systematic review in
*Annals of Internal Medicine* found 21 of 36 studies at serious risk of bias,
and named the mechanism: dose-response patterns "may also reflect residual
confounding if industry targets clinicians who already have higher baseline
prescribing volumes."

So the generator supports three targeting regimes, in increasing order of how
much trouble they cause:

``none``
    Payments land at random.  Every estimator here should recover the truth.

``static``
    Industry targets physicians with high *time-invariant* prescribing volume.
    Cross-sectional comparisons break; physician fixed effects survive, because
    the confounder is absorbed by the fixed effect.

``dynamic``
    Industry targets physicians whose prescribing has *recently been rising*.
    This is the hard case.  The confounder is time-varying and correlated with
    the outcome's own path, so physician fixed effects no longer rescue you --
    the treated units were on an upward trajectory before any payment arrived.
    Only a design that inspects pre-treatment trends can detect it.

Outcomes are in logs, so the treatment effect ``tau`` reads as an approximate
proportional change in prescribing volume.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

NEVER_TREATED = -1

VALID_TARGETING: tuple[str, ...] = ("none", "static", "dynamic")


@dataclass(frozen=True)
class PanelConfig:
    """Parameters of the physician-quarter data generating process.

    Attributes:
        n_physicians: Number of physicians in the panel.
        n_periods: Number of quarters observed per physician (balanced panel).
        tau: True causal effect of being in a payment relationship, in log
            points.  ``0.05`` is roughly a 5% lift in prescribing volume.
        targeting: One of ``"none"``, ``"static"`` or ``"dynamic"``.  See the
            module docstring for what each regime does to identification.
        sigma_alpha: SD of the physician fixed effect (persistent differences in
            prescribing volume across physicians).
        rho: AR(1) persistence of the idiosyncratic shock.  Persistence is what
            turns dynamic targeting into a pre-trend rather than noise.
        sigma_nu: SD of the AR(1) innovation.
        trend: Per-period drift common to all physicians (market growth).
        base_log_rx: Intercept in log claims, setting the overall volume scale.
        onset_intercept: Baseline log-odds of a first payment in a given period.
            More negative means fewer physicians ever receive a payment.
        onset_on_alpha: Weight on the physician fixed effect in the payment
            hazard.  Drives ``"static"`` and ``"dynamic"`` targeting.
        onset_on_growth: Weight on recent prescribing growth in the payment
            hazard.  Drives ``"dynamic"`` targeting only.
        burn_in: Periods at the start during which no payments occur, so that
            every eventually-treated physician has observable pre-periods.
        seed: Seed for the random number generator.
    """

    n_physicians: int = 2000
    n_periods: int = 16
    tau: float = 0.05
    targeting: str = "dynamic"
    sigma_alpha: float = 0.80
    rho: float = 0.60
    sigma_nu: float = 0.25
    trend: float = 0.01
    base_log_rx: float = 4.0
    onset_intercept: float = -2.80
    onset_on_alpha: float = 0.90
    onset_on_growth: float = 2.50
    burn_in: int = 3
    seed: int = 20260928

    def __post_init__(self) -> None:
        if self.targeting not in VALID_TARGETING:
            raise ValueError(
                f"targeting must be one of {VALID_TARGETING}, "
                f"got {self.targeting!r}"
            )
        if self.n_periods <= self.burn_in + 2:
            raise ValueError(
                "n_periods must exceed burn_in by at least 3 so that "
                "pre-treatment periods are observable"
            )
        if not 0.0 <= self.rho < 1.0:
            raise ValueError(f"rho must be in [0, 1), got {self.rho}")


def _logistic(x: np.ndarray) -> np.ndarray:
    """Numerically stable logistic function."""
    out = np.empty_like(x, dtype=float)
    positive = x >= 0
    out[positive] = 1.0 / (1.0 + np.exp(-x[positive]))
    exp_x = np.exp(x[~positive])
    out[~positive] = exp_x / (1.0 + exp_x)
    return out


def simulate_panel(config: PanelConfig | None = None) -> pd.DataFrame:
    """Generate a balanced physician-quarter panel with a known treatment effect.

    Payment onset is absorbing: once a physician enters a payment relationship
    they stay in it, which matches how these relationships behave in Open
    Payments and makes the panel a staggered-adoption design.

    Crucially, onset in period ``t`` is decided using information through
    ``t - 1`` only.  Outcomes are then drawn given that treatment status, so the
    simulation never leaks contemporaneous outcomes into assignment.

    Args:
        config: Data generating process parameters.  Defaults to
            :class:`PanelConfig` with dynamic targeting.

    Returns:
        A tidy DataFrame with one row per physician-period and columns:

        ``physician_id``, ``period``, ``log_rx``, ``rx_claims``, ``treated``
        (post-onset indicator), ``first_treat_period`` (``-1`` if never
        treated), ``event_time`` (periods since onset, ``NaN`` if never
        treated), ``ever_treated``, and ``alpha`` -- the physician fixed effect,
        which is observable here and *not* observable in real data.  It is
        carried so tests and diagnostics can check whether an estimator is
        merely rediscovering it.

        ``df.attrs["true_effect"]`` holds the ``tau`` used, so downstream code
        can compute bias without the caller passing it around.
    """
    cfg = config or PanelConfig()
    rng = np.random.default_rng(cfg.seed)

    n, t_max = cfg.n_physicians, cfg.n_periods

    alpha = rng.normal(0.0, cfg.sigma_alpha, size=n)
    gamma = cfg.trend * np.arange(t_max, dtype=float)
    innovations = rng.normal(0.0, cfg.sigma_nu, size=(n, t_max))

    eps = np.zeros((n, t_max))
    log_rx = np.zeros((n, t_max))
    treated = np.zeros((n, t_max), dtype=bool)
    first_treat = np.full(n, NEVER_TREATED, dtype=int)

    for t in range(t_max):
        eps[:, t] = (
            innovations[:, t]
            if t == 0
            else cfg.rho * eps[:, t - 1] + innovations[:, t]
        )

        # Decide payment onset for this period using only information the
        # manufacturer could have had, i.e. outcomes through t - 1.
        if t >= cfg.burn_in and cfg.targeting != "none":
            hazard = np.full(n, cfg.onset_intercept)
            hazard += cfg.onset_on_alpha * alpha
            if cfg.targeting == "dynamic" and t >= 2:
                recent_growth = log_rx[:, t - 1] - log_rx[:, t - 2]
                hazard += cfg.onset_on_growth * recent_growth
            onset = (first_treat == NEVER_TREATED) & (
                rng.random(n) < _logistic(hazard)
            )
            first_treat[onset] = t
        elif t >= cfg.burn_in:
            onset = (first_treat == NEVER_TREATED) & (
                rng.random(n) < _logistic(np.full(n, cfg.onset_intercept))
            )
            first_treat[onset] = t

        treated[:, t] = (first_treat != NEVER_TREATED) & (first_treat <= t)
        log_rx[:, t] = (
            cfg.base_log_rx + alpha + gamma[t] + cfg.tau * treated[:, t] + eps[:, t]
        )

    physician_id = np.repeat(np.arange(n), t_max)
    period = np.tile(np.arange(t_max), n)
    first_treat_long = np.repeat(first_treat, t_max)

    event_time = np.where(
        first_treat_long == NEVER_TREATED,
        np.nan,
        period - first_treat_long,
    )

    frame = pd.DataFrame(
        {
            "physician_id": physician_id,
            "period": period,
            "log_rx": log_rx.ravel(),
            "rx_claims": np.rint(np.exp(log_rx.ravel())).astype(int),
            "treated": treated.ravel(),
            "first_treat_period": first_treat_long,
            "event_time": event_time,
            "ever_treated": first_treat_long != NEVER_TREATED,
            "alpha": np.repeat(alpha, t_max),
        }
    )
    frame.attrs["true_effect"] = cfg.tau
    frame.attrs["targeting"] = cfg.targeting
    return frame
