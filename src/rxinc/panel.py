"""Assemble a physician-period panel from linked CMS data.

The output schema deliberately matches :func:`rxinc.simulate.simulate_panel`,
so every estimator and diagnostic runs unchanged on real and simulated data.
That symmetry is the point: you calibrate your expectations about an estimator
on data where truth is known, then apply the identical code path to data where
it is not.

Two modelling choices are exposed rather than buried, because both move the
answer:

``payment_threshold``
    Treatment is "in a payment relationship", not "received any transfer of
    value". A single $14 sandwich is a different object from a speaking
    contract, and a threshold of zero makes almost everyone treated.

``min_claims``
    Prescribers with very few claims produce wildly noisy log outcomes and
    Part D suppresses small counts, so a floor keeps the panel estimable.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from rxinc.simulate import NEVER_TREATED


@dataclass
class PanelReport:
    """What survived panel construction.

    Attributes:
        n_physicians: Physicians in the final panel.
        n_periods: Periods in the final panel.
        n_treated: Physicians with an observed payment onset.
        n_never_treated: Physicians never crossing the payment threshold.
        dropped_unbalanced: Physicians dropped for missing periods.
        dropped_low_volume: Physicians dropped by the claims floor.
        dropped_always_treated: Physicians treated in the first period, who
            have no usable pre-period and so identify nothing.
    """

    n_physicians: int
    n_periods: int
    n_treated: int
    n_never_treated: int
    dropped_unbalanced: int
    dropped_low_volume: int
    dropped_always_treated: int

    def __str__(self) -> str:
        return (
            f"Panel: {self.n_physicians:,} physicians x {self.n_periods} periods "
            f"({self.n_treated:,} treated, {self.n_never_treated:,} never treated). "
            f"Dropped {self.dropped_low_volume:,} low-volume, "
            f"{self.dropped_unbalanced:,} unbalanced, "
            f"{self.dropped_always_treated:,} always-treated."
        )


def build_panel(
    prescribing: pd.DataFrame,
    payments: pd.DataFrame,
    npi_col: str = "npi",
    period_col: str = "period",
    claims_col: str = "claims",
    amount_col: str = "amount",
    payment_threshold: float = 100.0,
    min_claims: int = 11,
    require_balanced: bool = True,
) -> tuple[pd.DataFrame, PanelReport]:
    """Build an estimator-ready physician-period panel.

    Args:
        prescribing: One row per physician-period with claim counts. Needs
            ``npi_col``, ``period_col`` and ``claims_col``.
        payments: Payment records with ``npi_col``, ``period_col`` and
            ``amount_col``. Multiple rows per physician-period are summed.
        npi_col: Physician identifier column, present in both frames.
        period_col: Period column, present in both frames. Periods must sort
            in chronological order.
        claims_col: Claim count column in ``prescribing``.
        amount_col: Payment amount column in ``payments``.
        payment_threshold: Total payment value within a period at or above
            which a physician counts as being in a payment relationship.
        min_claims: Drop physicians whose claims fall below this in any
            period. Part D suppresses counts below 11, so the default matches
            the suppression floor.
        require_balanced: Drop physicians not observed in every period.

    Returns:
        ``(panel, report)``. The panel carries the same columns as a simulated
        one: ``physician_id``, ``period``, ``log_rx``, ``rx_claims``,
        ``treated``, ``first_treat_period``, ``event_time``, ``ever_treated``
        -- plus ``payment_amount``. There is no ``alpha`` column, because on
        real data the physician effect is exactly what you cannot observe.

    Raises:
        ValueError: If required columns are missing or the panel is empty
            after filtering.
    """
    for frame, cols, label in (
        (prescribing, (npi_col, period_col, claims_col), "prescribing"),
        (payments, (npi_col, period_col, amount_col), "payments"),
    ):
        missing = [c for c in cols if c not in frame.columns]
        if missing:
            raise ValueError(f"{label} is missing required columns: {missing}")

    rx = prescribing[[npi_col, period_col, claims_col]].copy()
    rx = rx.groupby([npi_col, period_col], as_index=False)[claims_col].sum()

    periods = sorted(rx[period_col].unique())
    period_index = {p: i for i, p in enumerate(periods)}

    n_start = rx[npi_col].nunique()

    low_volume = set(rx.loc[rx[claims_col] < min_claims, npi_col].unique())
    rx = rx[~rx[npi_col].isin(low_volume)]

    dropped_unbalanced = 0
    if require_balanced:
        counts = rx.groupby(npi_col)[period_col].nunique()
        complete = set(counts[counts == len(periods)].index)
        dropped_unbalanced = int(rx[npi_col].nunique() - len(complete))
        rx = rx[rx[npi_col].isin(complete)]

    if rx.empty:
        raise ValueError(
            "No physicians survived filtering; lower min_claims or set "
            "require_balanced=False."
        )

    pay = payments[[npi_col, period_col, amount_col]].copy()
    pay = pay[pay[npi_col].isin(set(rx[npi_col]))]
    pay = pay.groupby([npi_col, period_col], as_index=False)[amount_col].sum()

    panel = rx.merge(pay, on=[npi_col, period_col], how="left")
    panel[amount_col] = panel[amount_col].fillna(0.0)
    panel["_period_ix"] = panel[period_col].map(period_index)
    panel = panel.sort_values([npi_col, "_period_ix"]).reset_index(drop=True)

    in_relationship = panel[amount_col] >= payment_threshold
    onset_ix = (
        panel.loc[in_relationship]
        .groupby(npi_col)["_period_ix"]
        .min()
        .rename("first_treat_period")
    )
    panel = panel.merge(onset_ix, on=npi_col, how="left")
    panel["first_treat_period"] = (
        panel["first_treat_period"].fillna(NEVER_TREATED).astype(int)
    )

    # A physician already in a relationship in period 0 has no pre-period, so
    # contributes no within-physician identifying variation.
    always_treated = set(
        panel.loc[panel["first_treat_period"] == 0, npi_col].unique()
    )
    panel = panel[~panel[npi_col].isin(always_treated)]
    if panel.empty:
        raise ValueError("Every physician was treated in the first period.")

    ever = panel["first_treat_period"] != NEVER_TREATED
    panel["treated"] = ever & (panel["_period_ix"] >= panel["first_treat_period"])
    panel["ever_treated"] = ever
    panel["event_time"] = np.where(
        ever, panel["_period_ix"] - panel["first_treat_period"], np.nan
    )
    panel["log_rx"] = np.log(panel[claims_col].to_numpy(dtype=float))

    out = panel.rename(
        columns={
            npi_col: "physician_id",
            claims_col: "rx_claims",
            amount_col: "payment_amount",
        }
    )[
        [
            "physician_id",
            "_period_ix",
            "log_rx",
            "rx_claims",
            "treated",
            "first_treat_period",
            "event_time",
            "ever_treated",
            "payment_amount",
        ]
    ].rename(columns={"_period_ix": "period"})

    out.attrs["true_effect"] = None  # unknowable on real data, by construction
    out.attrs["period_labels"] = periods

    by_unit = out.groupby("physician_id")["ever_treated"].first()
    report = PanelReport(
        n_physicians=int(out["physician_id"].nunique()),
        n_periods=len(periods),
        n_treated=int(by_unit.sum()),
        n_never_treated=int((~by_unit).sum()),
        dropped_unbalanced=dropped_unbalanced,
        dropped_low_volume=len(low_volume),
        dropped_always_treated=len(always_treated),
    )
    _ = n_start
    return out.reset_index(drop=True), report
