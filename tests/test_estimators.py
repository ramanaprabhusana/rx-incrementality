"""Tests for the estimators, including correctness of the within transform."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from rxinc.estimators import (
    Estimate,
    _ols,
    _within_transform,
    did_never_treated,
    event_study,
    interrupted_time_series,
    naive_ols,
    twoway_fe,
)
from rxinc.simulate import PanelConfig, simulate_panel


def _explicit_dummy_fe(frame, outcome="log_rx", treatment="treated"):
    """Reference implementation: full unit and period dummies via lstsq.

    Deliberately naive and memory-hungry, which is exactly what makes it a
    trustworthy check on the alternating-projection transform.
    """
    y = frame[outcome].to_numpy(float)
    d = frame[treatment].to_numpy(float).reshape(-1, 1)
    unit = pd.get_dummies(frame["physician_id"], drop_first=True).to_numpy(float)
    period = pd.get_dummies(frame["period"], drop_first=True).to_numpy(float)
    X = np.column_stack([np.ones(len(frame)), d, unit, period])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return float(beta[1])


def test_within_transform_matches_explicit_dummies_balanced():
    df = simulate_panel(PanelConfig(n_physicians=60, n_periods=8))
    assert twoway_fe(df).coef == pytest.approx(_explicit_dummy_fe(df), abs=1e-8)


def test_within_transform_matches_explicit_dummies_unbalanced():
    """The case that a single demeaning pass gets wrong."""
    df = simulate_panel(PanelConfig(n_physicians=60, n_periods=8, seed=3))
    rng = np.random.default_rng(0)
    keep = rng.random(len(df)) < 0.75
    unbalanced = df.loc[keep].copy()
    # Keep only physicians who still have variation in both dimensions.
    counts = unbalanced.groupby("physician_id").size()
    unbalanced = unbalanced[
        unbalanced["physician_id"].isin(counts[counts >= 4].index)
    ].copy()
    unbalanced.attrs["true_effect"] = df.attrs["true_effect"]

    assert twoway_fe(unbalanced).coef == pytest.approx(
        _explicit_dummy_fe(unbalanced), abs=1e-6
    )


def test_within_transform_removes_both_margins():
    df = simulate_panel(PanelConfig(n_physicians=40, n_periods=6))
    out = _within_transform(df, ["log_rx"], unit="physician_id", time="period")
    tmp = df.assign(z=out[:, 0])
    assert tmp.groupby("physician_id")["z"].mean().abs().max() < 1e-8
    assert tmp.groupby("period")["z"].mean().abs().max() < 1e-8


def test_clustered_se_exceeds_classical_under_within_unit_correlation():
    df = simulate_panel(PanelConfig(n_physicians=400, rho=0.8))
    y = df["log_rx"].to_numpy(float)
    d = df["treated"].to_numpy(float)
    X = np.column_stack([np.ones_like(d), d])
    _, v_classical, _ = _ols(y, X)
    _, v_clustered, _ = _ols(y, X, cluster=df["physician_id"].to_numpy())
    assert np.sqrt(v_clustered[1, 1]) > np.sqrt(v_classical[1, 1])


def test_naive_ols_is_severely_biased_under_static_targeting():
    df = simulate_panel(PanelConfig(n_physicians=1500, targeting="static"))
    est = naive_ols(df)
    assert est.coef > 5 * df.attrs["true_effect"]
    assert est.covers_truth is False


def test_twoway_fe_is_unbiased_under_static_targeting():
    """Averaged over replications, not asserted from a single lucky draw."""
    base = PanelConfig(n_physicians=800, targeting="static")
    coefs = [
        twoway_fe(simulate_panel(replace(base, seed=base.seed + r))).coef
        for r in range(25)
    ]
    assert np.mean(coefs) == pytest.approx(base.tau, abs=0.008)


def test_twoway_fe_is_biased_upward_under_dynamic_targeting():
    base = PanelConfig(n_physicians=800, targeting="dynamic")
    coefs = [
        twoway_fe(simulate_panel(replace(base, seed=base.seed + r))).coef
        for r in range(25)
    ]
    assert np.mean(coefs) > base.tau * 1.5


def test_event_study_omits_reference_period():
    df = simulate_panel(PanelConfig(n_physicians=300))
    result = event_study(df, leads=3, lags=4)
    assert -1 not in set(result.rel_periods.tolist())
    assert len(result.coefs) == len(result.rel_periods) == 7
    assert result.vcov.shape == (7, 7)


def test_event_study_frame_flags_pre_periods():
    df = simulate_panel(PanelConfig(n_physicians=200))
    frame = event_study(df, leads=2, lags=3).to_frame()
    assert frame.loc[frame["rel_period"] < 0, "is_pre"].all()
    assert not frame.loc[frame["rel_period"] >= 0, "is_pre"].any()


def test_did_requires_never_treated_units():
    df = simulate_panel(PanelConfig(n_physicians=100))
    df = df[df["ever_treated"]].copy()
    df.attrs["true_effect"] = 0.05
    with pytest.raises(ValueError, match="never-treated"):
        did_never_treated(df, n_boot=0)


def test_its_requires_treated_units():
    df = simulate_panel(PanelConfig(n_physicians=100))
    df = df[~df["ever_treated"]].copy()
    with pytest.raises(ValueError, match="requires treated"):
        interrupted_time_series(df)


def test_its_returns_level_and_slope():
    df = simulate_panel(PanelConfig(n_physicians=300))
    level, slope = interrupted_time_series(df)
    assert level.name == "ITS level change"
    assert slope.name == "ITS slope change"
    assert slope.true_effect is None  # slope has no counterpart in the DGP


def test_estimate_reports_bias_only_when_truth_known():
    known = Estimate("x", coef=0.08, se=0.01, n_obs=10, n_clusters=5, true_effect=0.05)
    assert known.bias == pytest.approx(0.03)
    assert known.bias_pct == pytest.approx(60.0)

    unknown = Estimate("y", coef=0.08, se=0.01, n_obs=10, n_clusters=5)
    assert unknown.bias is None
    assert unknown.covers_truth is None
    assert "bias" not in str(unknown)


def test_coverage_is_none_when_se_missing():
    est = Estimate("z", coef=0.05, se=np.nan, n_obs=10, n_clusters=5, true_effect=0.05)
    assert est.covers_truth is None


def test_zero_true_effect_formats_without_percentage():
    est = Estimate("p", coef=0.02, se=0.01, n_obs=10, n_clusters=5, true_effect=0.0)
    assert est.bias_pct is None
    rendered = str(est)
    assert rendered.endswith("bias +0.0200")
    assert "%)" not in rendered  # the percentage-bias suffix, not the "95% CI" label
