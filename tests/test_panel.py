"""Tests for real-data panel construction and linkage."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from rxinc.estimators import twoway_fe
from rxinc.linkage import is_valid_npi, link_payments_to_prescribers, normalize_name
from rxinc.panel import build_panel

# NPIs observed live in the CMS Part D API, so the check digits are genuine.
REAL_NPI_A = "1053753350"
REAL_NPI_B = "1053757047"
REAL_NPI_C = "1053758334"


def test_npi_checksum_accepts_real_and_rejects_corrupted():
    assert is_valid_npi(REAL_NPI_A)
    assert is_valid_npi(REAL_NPI_B)
    assert not is_valid_npi(REAL_NPI_A[:-1] + "1")
    assert not is_valid_npi("0000000000")
    assert not is_valid_npi("123")
    assert not is_valid_npi(None)


def test_normalize_name_strips_credentials_and_punctuation():
    assert normalize_name("O'Brien-Smith, MD") == "OBRIENSMITH"
    assert normalize_name("Sana Jr.") == "SANA"
    assert normalize_name(np.nan) == ""


def _prescribers():
    return pd.DataFrame(
        {
            "Prscrbr_NPI": [REAL_NPI_A, REAL_NPI_B, REAL_NPI_C],
            "Prscrbr_Last_Org_Name": ["Tadang", "Davitto", "Feipel"],
            "Prscrbr_First_Name": ["Brian", "Anna", "Mark"],
            "Prscrbr_State_Abrvtn": ["CA", "IN", "IN"],
        }
    )


def test_linkage_prefers_npi():
    payments = pd.DataFrame(
        {
            "Covered_Recipient_NPI": [REAL_NPI_A, REAL_NPI_B],
            "Covered_Recipient_Last_Name": ["Tadang", "Davitto"],
            "Covered_Recipient_First_Name": ["Brian", "Anna"],
            "Recipient_State": ["CA", "IN"],
        }
    )
    linked, report = link_payments_to_prescribers(payments, _prescribers())
    assert list(linked["match_method"]) == ["npi", "npi"]
    assert report.matched_by_npi == 2
    assert report.match_rate == pytest.approx(1.0)
    assert report.fallback_share == pytest.approx(0.0)


def test_linkage_falls_back_to_name_when_npi_missing():
    payments = pd.DataFrame(
        {
            "Covered_Recipient_NPI": ["", ""],
            "Covered_Recipient_Last_Name": ["Davitto", "Feipel"],
            "Covered_Recipient_First_Name": ["Anna", "Mark"],
            "Recipient_State": ["IN", "IN"],
        }
    )
    linked, report = link_payments_to_prescribers(payments, _prescribers())
    assert list(linked["match_method"]) == ["name", "name"]
    assert report.matched_by_name == 2
    assert report.fallback_share == pytest.approx(1.0)


def test_linkage_rejects_invalid_npi_and_counts_it():
    payments = pd.DataFrame(
        {
            "Covered_Recipient_NPI": ["0000000000"],
            "Covered_Recipient_Last_Name": ["Nobody"],
            "Covered_Recipient_First_Name": ["Nemo"],
            "Recipient_State": ["TX"],
        }
    )
    _, report = link_payments_to_prescribers(payments, _prescribers())
    assert report.invalid_npi == 1
    assert report.unmatched == 1


def test_linkage_leaves_ambiguous_blocks_unmatched():
    """A wrong link fabricates a treatment assignment, so refuse to guess."""
    prescribers = pd.DataFrame(
        {
            "Prscrbr_NPI": [REAL_NPI_A, REAL_NPI_B],
            "Prscrbr_Last_Org_Name": ["Smith", "Smith"],
            "Prscrbr_First_Name": ["John", "James"],
            "Prscrbr_State_Abrvtn": ["CA", "CA"],
        }
    )
    payments = pd.DataFrame(
        {
            "Covered_Recipient_NPI": [""],
            "Covered_Recipient_Last_Name": ["Smith"],
            "Covered_Recipient_First_Name": ["J"],
            "Recipient_State": ["CA"],
        }
    )
    linked, report = link_payments_to_prescribers(payments, prescribers)
    assert linked.loc[0, "match_method"] == "unmatched"
    assert report.ambiguous_name == 1


def _toy_inputs(n=60, t=6, seed=0):
    rng = np.random.default_rng(seed)
    npis = [f"P{i:04d}" for i in range(n)]
    rows = []
    for i, npi in enumerate(npis):
        base = rng.uniform(4.0, 6.0)
        for p in range(t):
            rows.append(
                {"npi": npi, "period": 2018 + p, "claims": int(np.exp(base + 0.02 * p))}
            )
    prescribing = pd.DataFrame(rows)
    pays = [
        {"npi": npis[i], "period": 2018 + (i % (t - 2)) + 2, "amount": 500.0}
        for i in range(0, n, 2)
    ]
    return prescribing, pd.DataFrame(pays)


def test_build_panel_matches_simulation_schema():
    prescribing, payments = _toy_inputs()
    panel, report = build_panel(prescribing, payments)
    for col in (
        "physician_id",
        "period",
        "log_rx",
        "rx_claims",
        "treated",
        "first_treat_period",
        "event_time",
        "ever_treated",
    ):
        assert col in panel.columns
    assert report.n_physicians > 0
    assert report.n_never_treated > 0
    # The whole point: estimators run unchanged on a real-shaped panel.
    est = twoway_fe(panel)
    assert np.isfinite(est.coef)
    assert est.true_effect is None


def test_build_panel_respects_payment_threshold():
    prescribing, payments = _toy_inputs()
    payments["amount"] = 50.0
    panel, report = build_panel(prescribing, payments, payment_threshold=100.0)
    assert report.n_treated == 0
    assert not panel["treated"].any()


def test_build_panel_drops_low_volume_prescribers():
    prescribing, payments = _toy_inputs()
    prescribing.loc[prescribing["npi"] == "P0000", "claims"] = 5
    _, report = build_panel(prescribing, payments, min_claims=11)
    assert report.dropped_low_volume >= 1


def test_build_panel_requires_expected_columns():
    prescribing, payments = _toy_inputs()
    with pytest.raises(ValueError, match="missing required columns"):
        build_panel(prescribing.drop(columns=["claims"]), payments)


def test_build_panel_treatment_is_absorbing():
    prescribing, payments = _toy_inputs()
    panel, _ = build_panel(prescribing, payments)
    for _, unit in panel.groupby("physician_id"):
        flags = unit.sort_values("period")["treated"].to_numpy()
        assert not np.any(flags[:-1] & ~flags[1:])
