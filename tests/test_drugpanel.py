"""Tests for the physician-drug panel and the triple-difference estimator."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from rxinc.diagnostics import pretrend_test
from rxinc.drugpanel import NEVER_TREATED, DrugPanelConfig, simulate_drug_panel
from rxinc.estimators import drug_event_study, triple_diff, twoway_fe


def test_panel_shape_and_keys():
    cfg = DrugPanelConfig(n_physicians=40, n_drugs=4, n_periods=6)
    df = simulate_drug_panel(cfg)
    assert len(df) == 40 * 4 * 6
    assert df["physician_year"].nunique() == 40 * 6
    assert df["drug_year"].nunique() == 4 * 6


def test_rejects_single_drug():
    with pytest.raises(ValueError, match="n_drugs must be at least 2"):
        DrugPanelConfig(n_drugs=1)


def test_rejects_unknown_targeting():
    with pytest.raises(ValueError, match="targeting must be one of"):
        DrugPanelConfig(targeting="vibes")


def test_rejects_nonstationary_pair_persistence():
    with pytest.raises(ValueError, match="rho_pair must be in"):
        DrugPanelConfig(rho_pair=1.0)


def test_reproducible_by_seed():
    a = simulate_drug_panel(DrugPanelConfig(n_physicians=30, seed=11))
    b = simulate_drug_panel(DrugPanelConfig(n_physicians=30, seed=11))
    np.testing.assert_allclose(a["log_rx"], b["log_rx"])


def test_treatment_absorbing_within_pair():
    df = simulate_drug_panel(DrugPanelConfig(n_physicians=60, n_drugs=3))
    for _, pair in df.groupby(["physician_id", "drug_id"]):
        flags = pair.sort_values("period")["treated"].to_numpy()
        assert not np.any(flags[:-1] & ~flags[1:])


def test_no_treatment_during_burn_in():
    cfg = DrugPanelConfig(n_physicians=80, burn_in=3, n_periods=8)
    df = simulate_drug_panel(cfg)
    assert not df[df["period"] < cfg.burn_in]["treated"].any()


def test_within_physician_year_variation_exists():
    """Without this, the triple difference has nothing to identify from."""
    df = simulate_drug_panel(DrugPanelConfig(n_physicians=300))
    mixed = df.groupby(["physician_id", "period"])["treated"].apply(
        lambda s: bool(s.any() and not s.all())
    )
    assert mixed.mean() > 0.2


def test_event_time_nan_for_untreated_pairs():
    df = simulate_drug_panel(DrugPanelConfig(n_physicians=50))
    untreated = df[df["first_treat_period"] == NEVER_TREATED]
    assert untreated["event_time"].isna().all()


def test_triple_diff_requires_its_columns():
    df = simulate_drug_panel(DrugPanelConfig(n_physicians=20))
    with pytest.raises(ValueError, match="triple_diff needs columns"):
        triple_diff(df.drop(columns=["drug_year"]))


def test_triple_diff_rejects_unidentified_panel():
    """Every physician paid about all or none of their drugs identifies nothing."""
    df = simulate_drug_panel(DrugPanelConfig(n_physicians=40, n_drugs=3))
    df = df.copy()
    df["treated"] = df["period"] >= 4  # identical across drugs within physician-year
    df.attrs["true_effect"] = 0.05
    with pytest.raises(ValueError, match="no variation within physician-year"):
        triple_diff(df)


def _mean_coef(regime: str, fn, n_reps: int = 12, n_physicians: int = 500) -> float:
    base = DrugPanelConfig(targeting=regime, n_physicians=n_physicians)
    return float(
        np.mean(
            [
                fn(simulate_drug_panel(replace(base, seed=base.seed + r))).coef
                for r in range(n_reps)
            ]
        )
    )


def _pair_fe(df: pd.DataFrame):
    """Physician-drug + period FE: same estimand, no physician-year absorption."""
    work = df.copy()
    work["pair"] = (
        work["physician_id"].astype(str) + "_" + work["drug_id"].astype(str)
    )
    work = work.rename(columns={"physician_id": "_orig", "pair": "physician_id"})
    work.attrs["true_effect"] = df.attrs["true_effect"]
    return twoway_fe(work)


def test_triple_diff_recovers_truth_under_random_assignment():
    truth = DrugPanelConfig().tau
    assert _mean_coef("none", triple_diff) == pytest.approx(truth, abs=0.015)


def test_triple_diff_survives_physician_momentum():
    """The headline claim: it withstands the selection that breaks panel designs."""
    truth = DrugPanelConfig().tau
    assert _mean_coef("physician_trajectory", triple_diff) == pytest.approx(
        truth, abs=0.015
    )


def test_pair_fe_fails_under_physician_momentum():
    """The design it replaces does not survive the same selection."""
    truth = DrugPanelConfig().tau
    assert _mean_coef("physician_trajectory", _pair_fe) > truth * 1.5


def test_triple_diff_is_broken_by_drug_specific_momentum():
    """Stated plainly, because the design is not a universal solvent."""
    truth = DrugPanelConfig().tau
    assert _mean_coef("drug_trajectory", triple_diff) > truth * 1.3


def test_drug_pretrends_flag_drug_momentum_and_clear_physician_momentum():
    """The diagnostic must separate the case it survives from the one it does not."""
    ok = simulate_drug_panel(
        DrugPanelConfig(targeting="physician_trajectory", n_physicians=900)
    )
    bad = simulate_drug_panel(
        DrugPanelConfig(targeting="drug_trajectory", n_physicians=900)
    )
    assert pretrend_test(drug_event_study(ok)).passes
    assert not pretrend_test(drug_event_study(bad)).passes


def test_drug_event_study_omits_reference_period():
    df = simulate_drug_panel(DrugPanelConfig(n_physicians=150))
    result = drug_event_study(df, leads=2, lags=3)
    assert -1 not in set(result.rel_periods.tolist())
    assert len(result.coefs) == 5
