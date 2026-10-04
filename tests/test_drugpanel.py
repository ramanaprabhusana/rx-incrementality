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


def test_three_factor_absorption_matches_explicit_dummies():
    """The k-factor engine against brute-force dummies, unbalanced on purpose."""
    from rxinc.estimators import _absorb
    df = simulate_drug_panel(DrugPanelConfig(n_physicians=25, n_drugs=3, n_periods=5))
    df = df.sample(frac=0.8, random_state=1).reset_index(drop=True)
    factors = ["physician_year", "drug_year", "physician_drug"]
    d = _absorb(df, ["log_rx", "treated"], factors)
    coef = float((d[:, 1] @ d[:, 0]) / (d[:, 1] @ d[:, 1]))
    X = [np.ones(len(df)), df["treated"].astype(float).to_numpy()]
    for f in factors:
        X.append(pd.get_dummies(df[f], drop_first=True).to_numpy(float))
    X = np.column_stack(X)
    beta, *_ = np.linalg.lstsq(X, df["log_rx"].to_numpy(), rcond=None)
    assert coef == pytest.approx(beta[1], abs=1e-6)


def test_nested_factors_excluded_from_dof():
    from rxinc.estimators import _absorbed_dof
    df = simulate_drug_panel(DrugPanelConfig(n_physicians=30, n_drugs=4, n_periods=5))
    # physician_year and physician_drug sit inside physician clusters; drug_year does not.
    dof = _absorbed_dof(df, ["physician_year", "drug_year", "physician_drug"], "physician_id")
    assert dof == df["drug_year"].nunique() - 2


def _three_way(df):
    return triple_diff(df, absorb=("physician_year", "drug_year", "physician_drug"))


def test_affinity_targeting_breaks_two_way_and_three_way_survives():
    """Level selection at the physician-drug grain: the gap the real data exposed."""
    truth = DrugPanelConfig().tau
    two_way = _mean_coef("affinity", triple_diff)
    three_way = _mean_coef("affinity", _three_way)
    assert two_way > truth * 2
    assert three_way == pytest.approx(truth, abs=0.015)


def test_three_way_still_survives_physician_momentum():
    truth = DrugPanelConfig().tau
    assert _mean_coef("physician_trajectory", _three_way) == pytest.approx(truth, abs=0.015)


def test_default_effects_are_homogeneous():
    df = simulate_drug_panel(DrugPanelConfig(n_physicians=80))
    assert np.allclose(df.loc[df.treated, "cell_effect"], DrugPanelConfig().tau)
    assert (df.loc[~df.treated, "cell_effect"] == 0).all()


def test_dynamic_effects_grow_with_time_since_onset():
    df = simulate_drug_panel(DrugPanelConfig(n_physicians=200, tau=0.02, tau_dynamic=0.01))
    by_e = df[df.treated].groupby("event_time")["cell_effect"].mean()
    assert by_e.loc[0] == pytest.approx(0.02) and by_e.loc[2] == pytest.approx(0.04)


def _heterogeneous_cfg(**kw):
    return DrugPanelConfig(targeting="affinity", n_physicians=700, tau=0.03,
                           tau_dynamic=0.02, tau_cohort_slope=-0.008, **kw)


def test_cohort_event_study_recovers_heterogeneous_dynamic_path():
    """Mean error within 3.5 Monte Carlo standard errors of zero at every period.

    A fixed tolerance with few replications mistakes noise for bias: an earlier
    version of this test failed at six replications, and 40 fresh replications
    then showed every error within one standard error of zero.
    """
    from rxinc.estimators import cohort_event_study
    base = _heterogeneous_cfg()
    errs = {}
    for r in range(12):
        df = simulate_drug_panel(replace(base, seed=base.seed + r))
        truth = df[df.treated].groupby("event_time")["cell_effect"].mean()
        res = cohort_event_study(df)
        for e, b in zip(res.rel_periods, res.coefs, strict=True):
            if -3 <= e <= 3:
                errs.setdefault(int(e), []).append(b - (truth.get(e, 0.0) if e >= 0 else 0.0))
    for e, v in errs.items():
        v = np.asarray(v)
        mc_se = v.std(ddof=1) / np.sqrt(len(v))
        assert abs(v.mean()) < 3.5 * mc_se, f"event time {e}: {v.mean():+.4f} vs mc_se {mc_se:.4f}"


def test_cohort_weights_sum_to_one_per_period():
    from rxinc.estimators import cohort_event_study
    res = cohort_event_study(simulate_drug_panel(_heterogeneous_cfg(n_periods=7)))
    w = res.detail[res.detail.identified].groupby("e")["weight"].sum()
    assert np.allclose(w, 1.0)
    assert -1 not in set(res.rel_periods.tolist())


def test_cohort_event_study_requires_never_treated():
    from rxinc.estimators import cohort_event_study
    df = simulate_drug_panel(DrugPanelConfig(n_physicians=60))
    with pytest.raises(ValueError, match="never-treated"):
        cohort_event_study(df[df.ever_treated])


def test_event_study_drops_unreached_relative_periods():
    """Asking for more lags than any pair reaches must not invent a zero estimate."""
    df = simulate_drug_panel(DrugPanelConfig(n_physicians=200, n_periods=6, burn_in=2))
    max_e = int(df["event_time"].max())
    res = drug_event_study(df, leads=2, lags=max_e + 3)
    assert res.rel_periods.max() <= max_e
    assert (res.ses > 0).all()
