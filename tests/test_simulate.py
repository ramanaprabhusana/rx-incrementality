"""Tests for the data generating process."""

from __future__ import annotations

import numpy as np
import pytest

from rxinc.simulate import NEVER_TREATED, PanelConfig, simulate_panel


def test_panel_is_balanced():
    cfg = PanelConfig(n_physicians=50, n_periods=10)
    df = simulate_panel(cfg)
    assert len(df) == 500
    assert (df.groupby("physician_id").size() == 10).all()


def test_same_seed_reproduces():
    a = simulate_panel(PanelConfig(n_physicians=40, seed=7))
    b = simulate_panel(PanelConfig(n_physicians=40, seed=7))
    np.testing.assert_allclose(a["log_rx"], b["log_rx"])


def test_different_seed_differs():
    a = simulate_panel(PanelConfig(n_physicians=40, seed=7))
    b = simulate_panel(PanelConfig(n_physicians=40, seed=8))
    assert not np.allclose(a["log_rx"], b["log_rx"])


def test_rejects_unknown_targeting():
    with pytest.raises(ValueError, match="targeting must be one of"):
        PanelConfig(targeting="sideways")


def test_rejects_panel_too_short_for_pre_periods():
    with pytest.raises(ValueError, match="n_periods must exceed burn_in"):
        PanelConfig(n_periods=4, burn_in=3)


def test_rejects_nonstationary_rho():
    with pytest.raises(ValueError, match="rho must be in"):
        PanelConfig(rho=1.0)


def test_treatment_is_absorbing():
    df = simulate_panel(PanelConfig(n_physicians=200))
    for _, unit in df.groupby("physician_id"):
        flags = unit.sort_values("period")["treated"].to_numpy()
        # Once true, never false again.
        assert not np.any(flags[:-1] & ~flags[1:])


def test_no_treatment_during_burn_in():
    cfg = PanelConfig(n_physicians=300, burn_in=4)
    df = simulate_panel(cfg)
    early = df[df["period"] < cfg.burn_in]
    assert not early["treated"].any()


def test_some_physicians_never_treated():
    df = simulate_panel(PanelConfig(n_physicians=400))
    assert (df["first_treat_period"] == NEVER_TREATED).any()


def test_true_effect_recorded_in_attrs():
    cfg = PanelConfig(tau=0.11)
    df = simulate_panel(cfg)
    assert df.attrs["true_effect"] == pytest.approx(0.11)
    assert df.attrs["targeting"] == cfg.targeting


def test_random_targeting_is_uncorrelated_with_physician_effect():
    df = simulate_panel(PanelConfig(n_physicians=3000, targeting="none"))
    by_unit = df.groupby("physician_id").first()
    corr = np.corrcoef(by_unit["alpha"], by_unit["ever_treated"].astype(float))[0, 1]
    assert abs(corr) < 0.06


def test_static_targeting_selects_high_prescribers():
    df = simulate_panel(PanelConfig(n_physicians=3000, targeting="static"))
    by_unit = df.groupby("physician_id").first()
    corr = np.corrcoef(by_unit["alpha"], by_unit["ever_treated"].astype(float))[0, 1]
    assert corr > 0.25


def test_dynamic_targeting_creates_pre_onset_growth_gap():
    df = simulate_panel(PanelConfig(n_physicians=3000, targeting="dynamic"))
    onset = df["first_treat_period"].to_numpy()
    cutoff = np.where(onset == NEVER_TREATED, df["period"].max(), onset)
    pre = df[df["period"].to_numpy() < cutoff]
    growth = pre.sort_values("period").groupby("physician_id").agg(
        ever=("ever_treated", "first"),
        slope=("log_rx", lambda s: (s.iloc[-1] - s.iloc[0]) / max(len(s) - 1, 1)),
    )
    treated = growth.loc[growth["ever"], "slope"].mean()
    control = growth.loc[~growth["ever"], "slope"].mean()
    assert treated > control


def test_event_time_aligns_with_onset():
    df = simulate_panel(PanelConfig(n_physicians=100))
    treated = df[df["ever_treated"]]
    expected = treated["period"] - treated["first_treat_period"]
    np.testing.assert_allclose(treated["event_time"], expected)
    assert df.loc[~df["ever_treated"], "event_time"].isna().all()


def test_rx_claims_are_positive_integers():
    df = simulate_panel(PanelConfig(n_physicians=100))
    assert (df["rx_claims"] > 0).all()
    assert df["rx_claims"].dtype.kind in "iu"
