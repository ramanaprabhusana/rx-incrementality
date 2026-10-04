"""Tests for the real-data panel builder, on hand-built inputs with known answers."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from rxinc.drugdata import build_drug_panel, load_part_d, load_payments
from rxinc.drugpanel import DrugPanelConfig, simulate_drug_panel
from rxinc.estimators import triple_diff


def _part_d_csv(tmp_path):
    rows = []
    # Two physicians, three years. P1 prescribes Victoza (both pack sizes) and
    # Jardiance every year; P2 prescribes Jardiance every year and Victoza only
    # in 2020, so P2-Victoza is suppressed (absent) in 2019 and 2021.
    for y in (2019, 2020, 2021):
        rows += [
            ("1111111111", "Victoza 2-Pak", 20, y), ("1111111111", "Victoza 3-Pak", 30, y),
            ("1111111111", "Jardiance", 100, y), ("2222222222", "Jardiance", 40, y),
            ("1111111111", "Soliqua 100-33", 50, y),  # excluded family
        ]
    rows.append(("2222222222", "Victoza 3-Pak", 15, 2020))
    df = pd.DataFrame(rows, columns=["Prscrbr_NPI", "Brnd_Name", "Tot_Clms", "data_year"])
    df["Prscrbr_Type"] = "Internal Medicine"
    df["Prscrbr_City"] = "LAFAYETTE"
    df["Prscrbr_State_Abrvtn"] = "IN"
    path = tmp_path / "pd.csv"
    df.to_csv(path, index=False)
    return path


def _payments_csvs(tmp_path):
    cols = ["Covered_Recipient_NPI", "Total_Amount_of_Payment_USDollars"] + [
        f"Name_of_Drug_or_Biological_or_Device_or_Medical_Supply_{i}" for i in range(1, 6)
    ]
    def frame(rows):
        return pd.DataFrame(rows, columns=cols)
    p19 = frame([["2222222222", "30", "VICTOZA", None, None, None, None]])
    # 2020: one $40 lunch naming two products -> $20 each; P1 paid about Jardiance.
    p20 = frame([["1111111111", "40", "JARDIANCE", "SOMEOTHERDRUG", None, None, None]])
    p21 = frame([["1111111111", "10", "JARDIANCE", None, None, None, None],
                 ["1111111111", "10", "JARDIANCE", "JARDIANCE", None, None, None]])
    out = {}
    for y, f in ((2019, p19), (2020, p20), (2021, p21)):
        path = tmp_path / f"op{y}.csv"
        f.to_csv(path, index=False)
        out[y] = path
    return out


@pytest.fixture
def built(tmp_path):
    claims, attrs = load_part_d([_part_d_csv(tmp_path)])
    pay = load_payments(_payments_csvs(tmp_path))
    return build_drug_panel(claims, attrs, pay, min_family_peak_claims=0)


def test_pack_sizes_are_summed_into_one_family(tmp_path):
    claims, _ = load_part_d([_part_d_csv(tmp_path)])
    v = claims[(claims.npi == "1111111111") & (claims.family == "VICTOZA")]
    assert (v["claims"] == 50).all()


def test_excluded_families_never_enter(built):
    panel, report = built
    assert "SOLIQUA" not in set(panel["family"])
    assert report.families == ["JARDIANCE", "VICTOZA"]


def test_grid_is_complete_and_suppressed_cells_are_flagged(built):
    panel, _ = built
    assert len(panel) == 2 * 2 * 3
    cell = panel[(panel.npi == "2222222222") & (panel.family == "VICTOZA")].set_index("year")
    assert cell.loc[2019, "any_rx"] == 0 and np.isnan(cell.loc[2019, "log_rx"])
    assert cell.loc[2020, "any_rx"] == 1 and cell.loc[2020, "log_rx"] == pytest.approx(np.log(15))


def test_amount_split_across_listed_products(built):
    panel, _ = built
    cell = panel[(panel.npi == "1111111111") & (panel.family == "JARDIANCE") & (panel.year == 2020)]
    assert cell["pay_amount"].iloc[0] == pytest.approx(20.0)


def test_same_family_twice_on_one_record_counts_once(built):
    panel, _ = built
    cell = panel[(panel.npi == "1111111111") & (panel.family == "JARDIANCE") & (panel.year == 2021)]
    assert cell["pay_records"].iloc[0] == 2
    assert cell["pay_amount"].iloc[0] == pytest.approx(10.0 + 5.0)


def test_left_censoring_and_clean_onset(built):
    panel, report = built
    p2v = panel[(panel.npi == "2222222222") & (panel.family == "VICTOZA")]
    assert p2v["left_censored"].all() and p2v["event_time"].isna().all()
    p1j = panel[(panel.npi == "1111111111") & (panel.family == "JARDIANCE")].set_index("year")
    assert p1j.loc[2019, "event_time"] == -1 and p1j.loc[2020, "treated"] == 1
    assert report.n_left_censored_pairs == 1 and report.n_clean_onset_pairs == 1


def test_lag_and_lead_align_within_pair(built):
    panel, _ = built
    p1j = panel[(panel.npi == "1111111111") & (panel.family == "JARDIANCE")].set_index("year")
    assert np.isnan(p1j.loc[2019, "pay_any_lag"])
    assert p1j.loc[2021, "pay_any_lag"] == 1
    assert p1j.loc[2019, "pay_any_lead"] == 1
    assert np.isnan(p1j.loc[2021, "pay_any_lead"])


def test_peer_share_is_leave_one_out(built):
    panel, _ = built
    # Same city. In 2019 P2 is paid about Victoza, so P1's peer share is 1 and
    # P2's own payment does not count toward P2's peer share.
    v19 = panel[(panel.family == "VICTOZA") & (panel.year == 2019)].set_index("npi")
    assert v19.loc["1111111111", "peer_pay_share"] == 1.0
    assert v19.loc["2222222222", "peer_pay_share"] == 0.0


def test_controls_reported_and_do_not_disturb_treatment_on_simulation():
    df = simulate_drug_panel(DrugPanelConfig(n_physicians=400))
    rng = np.random.default_rng(0)
    df["noise_ctrl"] = rng.normal(size=len(df))
    base = triple_diff(df)
    with_ctrl = triple_diff(df, controls=["noise_ctrl"])
    assert "noise_ctrl" in with_ctrl.extra and "noise_ctrl_se" in with_ctrl.extra
    assert with_ctrl.coef == pytest.approx(base.coef, abs=0.01)
