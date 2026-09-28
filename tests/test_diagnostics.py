"""Tests for the diagnostic battery."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from rxinc.diagnostics import (
    balance_table,
    diagnose,
    monte_carlo,
    placebo_shift,
    pretrend_test,
)
from rxinc.estimators import event_study, naive_ols, twoway_fe
from rxinc.simulate import PanelConfig, simulate_panel


def test_pretrend_passes_under_static_targeting():
    df = simulate_panel(PanelConfig(n_physicians=1500, targeting="static"))
    result = pretrend_test(event_study(df))
    assert result.passes
    assert result.p_value >= 0.05


def test_pretrend_fails_loudly_under_dynamic_targeting():
    df = simulate_panel(PanelConfig(n_physicians=1500, targeting="dynamic"))
    result = pretrend_test(event_study(df))
    assert not result.passes
    assert result.p_value < 1e-6
    assert result.max_abs_coef > 0.05


def test_pretrend_requires_lead_coefficients():
    df = simulate_panel(PanelConfig(n_physicians=200))
    study = event_study(df, leads=0, lags=3)
    with pytest.raises(ValueError, match="no lead coefficients"):
        pretrend_test(study)


def test_placebo_is_unbiased_under_random_assignment():
    """Averaged across replications: a single draw cannot separate bias from noise."""
    base = PanelConfig(n_physicians=900, targeting="none")
    coefs = [
        placebo_shift(simulate_panel(replace(base, seed=base.seed + r))).coef
        for r in range(30)
    ]
    mean = float(np.mean(coefs))
    mc_se = float(np.std(coefs, ddof=1) / np.sqrt(len(coefs)))
    assert abs(mean) < 3 * mc_se


def test_placebo_detects_dynamic_confounding():
    base = PanelConfig(n_physicians=900, targeting="dynamic")
    coefs = [
        placebo_shift(simulate_panel(replace(base, seed=base.seed + r))).coef
        for r in range(30)
    ]
    assert float(np.mean(coefs)) > 0.01


def test_placebo_rejects_nonpositive_shift():
    df = simulate_panel(PanelConfig(n_physicians=100))
    with pytest.raises(ValueError, match="shift must be >= 1"):
        placebo_shift(df, shift=0)


def test_placebo_contains_no_genuinely_treated_observations():
    """The whole point of a placebo is that nothing in it has been treated."""
    df = simulate_panel(PanelConfig(n_physicians=300))
    est = placebo_shift(df, shift=2)
    assert est.true_effect == 0.0
    assert est.n_obs < len(df)


def test_balance_table_shape_and_direction():
    df = simulate_panel(PanelConfig(n_physicians=1200, targeting="static"))
    table = balance_table(df)
    assert set(table.index) == {False, True}
    assert list(table.columns) == [
        "n_physicians",
        "mean_pre_level",
        "mean_pre_growth",
    ]
    # Static targeting picks high-volume prescribers, visible before onset.
    assert table.loc[True, "mean_pre_level"] > table.loc[False, "mean_pre_level"]


def test_balance_table_shows_growth_gap_only_under_dynamic():
    static = balance_table(
        simulate_panel(PanelConfig(n_physicians=1200, targeting="static"))
    )
    dynamic = balance_table(
        simulate_panel(PanelConfig(n_physicians=1200, targeting="dynamic"))
    )
    static_gap = (
        static.loc[True, "mean_pre_growth"] - static.loc[False, "mean_pre_growth"]
    )
    dynamic_gap = (
        dynamic.loc[True, "mean_pre_growth"] - dynamic.loc[False, "mean_pre_growth"]
    )
    assert dynamic_gap > static_gap


def test_monte_carlo_reports_expected_columns():
    table = monte_carlo(
        PanelConfig(n_physicians=200, targeting="static"),
        {"Two-way fixed effects": twoway_fe, "Naive pooled OLS": naive_ols},
        n_reps=5,
    )
    assert list(table.columns) == [
        "mean_coef",
        "mean_bias",
        "sd_coef",
        "rmse",
        "coverage95",
        "n_reps",
    ]
    assert (table["n_reps"] == 5).all()
    assert set(table.index) == {"Two-way fixed effects", "Naive pooled OLS"}


def test_monte_carlo_ranks_fe_above_naive_on_rmse():
    table = monte_carlo(
        PanelConfig(n_physicians=400, targeting="static"),
        {"Two-way fixed effects": twoway_fe, "Naive pooled OLS": naive_ols},
        n_reps=10,
    )
    assert (
        table.loc["Two-way fixed effects", "rmse"]
        < table.loc["Naive pooled OLS", "rmse"]
    )


def test_diagnose_returns_full_battery():
    df = simulate_panel(PanelConfig(n_physicians=400))
    out = diagnose(df)
    assert set(out) == {"event_study", "pretrend", "placebo", "balance"}
    assert out["pretrend"].df > 0
