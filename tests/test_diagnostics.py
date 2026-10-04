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


def _synthetic_es(effect_by_e, slope, noise_sd=0.0, seed=0):
    from rxinc.estimators import EventStudyResult
    rel = np.array([r for r in range(-5, 5) if r != -1])
    rng = np.random.default_rng(seed)
    coefs = slope * (rel + 1) + np.array([effect_by_e.get(int(r), 0.0) for r in rel])
    coefs = coefs + rng.normal(0, noise_sd, len(rel))
    ses = np.full(len(rel), max(noise_sd, 1e-3))
    return EventStudyResult(rel_periods=rel, coefs=coefs, ses=ses, vcov=np.diag(ses**2),
                            n_obs=1, n_clusters=1)


def test_detrend_removes_a_pure_linear_pretrend():
    from rxinc.diagnostics import detrend_event_study
    effects = {0: 0.02, 1: 0.05, 2: 0.06, 3: 0.04, 4: 0.03}
    res = detrend_event_study(_synthetic_es(effects, slope=0.012))
    assert res.slope == pytest.approx(0.012, abs=1e-9)
    post = res.rel_periods >= 0
    assert np.allclose(res.coefs[post], [effects[int(e)] for e in res.rel_periods[post]], atol=1e-9)
    assert np.allclose(res.coefs[~post], 0.0, atol=1e-9)


def test_detrend_leaves_flat_pretrends_alone():
    from rxinc.diagnostics import detrend_event_study
    es = _synthetic_es({1: 0.05}, slope=0.0)
    res = detrend_event_study(es)
    assert res.slope == pytest.approx(0.0, abs=1e-12)
    assert np.allclose(res.coefs, es.coefs)


def test_detrend_inflates_post_uncertainty_for_slope_estimation():
    from rxinc.diagnostics import detrend_event_study
    es = _synthetic_es({1: 0.05}, slope=0.01, noise_sd=0.01, seed=3)
    res = detrend_event_study(es)
    far_post = res.rel_periods == 4
    assert res.ses[far_post][0] > es.ses[far_post][0]
    mean, se = res.post_average()
    assert np.isfinite(mean) and se > 0
