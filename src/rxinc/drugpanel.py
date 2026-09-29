"""Physician-drug-year panels, for validating the triple-difference design.

The single-outcome simulator in :mod:`rxinc.simulate` shows that
trajectory-based selection defeats every physician-level estimator. This module
adds the drug dimension that makes a better design possible, and lets that
design be tested against known truth before it is pointed at CMS data.

Three selection regimes, distinguished by *what* the manufacturer keys on:

``none``
    Payments land at random across physician-drug pairs.

``physician_trajectory``
    The manufacturer targets physicians whose overall prescribing is rising.
    Every drug of a rising physician is equally likely to attract a payment, so
    which drug is paid about is random within the physician-year. This is the
    regime that destroys physician-level estimators and that the triple
    difference should survive, because a physician-by-year effect absorbs the
    physician's whole trajectory.

``drug_trajectory``
    The manufacturer targets physicians whose prescribing *of that specific
    drug* is rising. The confounder now varies within physician-year, which is
    exactly the dimension the triple difference relies on, so it should fail
    here too. Its saving grace is that a drug-level pre-trend test can see it.

The point of separating the last two is honesty. The triple difference is not
a universal solvent, and the simulation is built to show precisely where it
stops working.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

NEVER_TREATED = -1

VALID_DRUG_TARGETING: tuple[str, ...] = (
    "none",
    "physician_trajectory",
    "drug_trajectory",
)


@dataclass(frozen=True)
class DrugPanelConfig:
    """Parameters for the physician-drug-year data generating process.

    Attributes:
        n_physicians: Number of physicians.
        n_drugs: Number of competing drugs in the class.
        n_periods: Number of annual periods.
        tau: True causal effect of a payment relationship on prescribing of the
            drug it concerns, in log points.
        targeting: One of :data:`VALID_DRUG_TARGETING`.
        sigma_physician: SD of the physician level effect.
        sigma_affinity: SD of the physician-drug affinity, capturing persistent
            preference of a physician for particular drugs.
        sigma_growth: SD of physician-specific trend, the trajectory that
            defeats physician-level designs.
        rho: AR(1) persistence of the physician-year shock.
        sigma_shock: SD of the AR(1) innovation at physician-year level.
        sigma_drug_year: SD of the drug-year effect, standing in for launches,
            guideline changes and competitor entry.
        sigma_noise: SD of idiosyncratic physician-drug-year noise.
        rho_pair: AR(1) persistence of the physician-drug shock. This is what
            gives drug-specific momentum real teeth: with ``rho_pair = 0`` a
            drug's recent growth is transitory noise that mean-reverts, so
            selecting on it produces visible pre-trends but no post-period bias.
            Persistence is what turns it into genuine confounding.
        sigma_pair_shock: SD of the physician-drug AR(1) innovation.
        base_log_rx: Intercept in log claims.
        onset_intercept: Baseline log-odds of a first payment.
        onset_on_trajectory: Weight on the relevant trajectory in the payment
            hazard. Drives both trajectory regimes.
        burn_in: Periods before any payment can occur.
        seed: Random seed.
    """

    n_physicians: int = 1500
    n_drugs: int = 6
    n_periods: int = 8
    tau: float = 0.05
    targeting: str = "physician_trajectory"
    sigma_physician: float = 0.70
    sigma_affinity: float = 0.60
    sigma_growth: float = 0.06
    rho: float = 0.55
    sigma_shock: float = 0.20
    sigma_drug_year: float = 0.15
    sigma_noise: float = 0.18
    rho_pair: float = 0.65
    sigma_pair_shock: float = 0.22
    base_log_rx: float = 3.2
    onset_intercept: float = -2.60
    onset_on_trajectory: float = 3.00
    burn_in: int = 2
    seed: int = 20260929

    def __post_init__(self) -> None:
        if self.targeting not in VALID_DRUG_TARGETING:
            raise ValueError(
                f"targeting must be one of {VALID_DRUG_TARGETING}, "
                f"got {self.targeting!r}"
            )
        if self.n_periods <= self.burn_in + 2:
            raise ValueError(
                "n_periods must exceed burn_in by at least 3 so pre-periods exist"
            )
        if not 0.0 <= self.rho_pair < 1.0:
            raise ValueError(f"rho_pair must be in [0, 1), got {self.rho_pair}")
        if self.n_drugs < 2:
            raise ValueError(
                "n_drugs must be at least 2; the design compares drugs within a "
                "physician-year"
            )


def _logistic(x: np.ndarray) -> np.ndarray:
    out = np.empty_like(x, dtype=float)
    pos = x >= 0
    out[pos] = 1.0 / (1.0 + np.exp(-x[pos]))
    ex = np.exp(x[~pos])
    out[~pos] = ex / (1.0 + ex)
    return out


def simulate_drug_panel(config: DrugPanelConfig | None = None) -> pd.DataFrame:
    """Generate a physician-drug-year panel with a known treatment effect.

    Payment onset is absorbing per physician-drug pair, and is decided using
    only information through the previous period, so contemporaneous outcomes
    never leak into assignment.

    Args:
        config: Data generating process parameters.

    Returns:
        A tidy DataFrame with one row per physician-drug-year and columns:
        ``physician_id``, ``drug_id``, ``period``, ``log_rx``, ``rx_claims``,
        ``treated``, ``first_treat_period``, ``event_time``, ``ever_treated``,
        plus ``physician_year`` and ``drug_year`` composite keys ready for the
        triple-difference transform.

        ``df.attrs["true_effect"]`` carries ``tau``.
    """
    cfg = config or DrugPanelConfig()
    rng = np.random.default_rng(cfg.seed)
    n, j, t_max = cfg.n_physicians, cfg.n_drugs, cfg.n_periods

    physician = rng.normal(0.0, cfg.sigma_physician, size=n)
    growth = rng.normal(0.0, cfg.sigma_growth, size=n)
    affinity = rng.normal(0.0, cfg.sigma_affinity, size=(n, j))
    drug_year = rng.normal(0.0, cfg.sigma_drug_year, size=(j, t_max)).cumsum(axis=1)
    shock_innov = rng.normal(0.0, cfg.sigma_shock, size=(n, t_max))
    noise = rng.normal(0.0, cfg.sigma_noise, size=(n, j, t_max))
    pair_innov = rng.normal(0.0, cfg.sigma_pair_shock, size=(n, j, t_max))

    shock = np.zeros((n, t_max))
    pair_shock = np.zeros((n, j, t_max))
    alpha = np.zeros((n, t_max))
    y = np.zeros((n, j, t_max))
    treated = np.zeros((n, j, t_max), dtype=bool)
    first_treat = np.full((n, j), NEVER_TREATED, dtype=int)

    for period in range(t_max):
        shock[:, period] = (
            shock_innov[:, period]
            if period == 0
            else cfg.rho * shock[:, period - 1] + shock_innov[:, period]
        )
        alpha[:, period] = physician + growth[:, None].ravel() * period + shock[:, period]
        pair_shock[:, :, period] = (
            pair_innov[:, :, period]
            if period == 0
            else cfg.rho_pair * pair_shock[:, :, period - 1] + pair_innov[:, :, period]
        )

        if period >= cfg.burn_in:
            hazard = np.full((n, j), cfg.onset_intercept)
            if cfg.targeting == "physician_trajectory" and period >= 2:
                # Physician-level momentum: identical across that physician's
                # drugs, so which drug is paid about stays random within the
                # physician-year.
                recent = alpha[:, period - 1] - alpha[:, period - 2]
                hazard += cfg.onset_on_trajectory * recent[:, None]
            elif cfg.targeting == "drug_trajectory" and period >= 2:
                # Drug-specific momentum: varies within physician-year, which
                # is the dimension the triple difference depends on.
                recent = y[:, :, period - 1] - y[:, :, period - 2]
                hazard += cfg.onset_on_trajectory * recent

            onset = (first_treat == NEVER_TREATED) & (rng.random((n, j)) < _logistic(hazard))
            first_treat[onset] = period

        treated[:, :, period] = (first_treat != NEVER_TREATED) & (first_treat <= period)
        y[:, :, period] = (
            cfg.base_log_rx
            + affinity
            + alpha[:, period][:, None]
            + drug_year[:, period][None, :]
            + cfg.tau * treated[:, :, period]
            + pair_shock[:, :, period]
            + noise[:, :, period]
        )

    physician_id = np.repeat(np.arange(n), j * t_max)
    drug_id = np.tile(np.repeat(np.arange(j), t_max), n)
    period_id = np.tile(np.arange(t_max), n * j)
    onset_long = np.repeat(first_treat.ravel(), t_max)

    frame = pd.DataFrame(
        {
            "physician_id": physician_id,
            "drug_id": drug_id,
            "period": period_id,
            "log_rx": y.reshape(-1),
            "rx_claims": np.rint(np.exp(y.reshape(-1))).astype(int),
            "treated": treated.reshape(-1),
            "first_treat_period": onset_long,
            "ever_treated": onset_long != NEVER_TREATED,
        }
    )
    frame["event_time"] = np.where(
        frame["ever_treated"],
        frame["period"] - frame["first_treat_period"],
        np.nan,
    )
    frame["physician_year"] = (
        frame["physician_id"].astype(str) + "_" + frame["period"].astype(str)
    )
    frame["drug_year"] = (
        frame["drug_id"].astype(str) + "_" + frame["period"].astype(str)
    )
    frame.attrs["true_effect"] = cfg.tau
    frame.attrs["targeting"] = cfg.targeting
    return frame
