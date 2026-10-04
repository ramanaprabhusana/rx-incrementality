"""Tests for the real-data panel builder, on hand-built inputs with known answers."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from rxinc.drugdata import (
    build_drug_panel,
    coverage_class,
    load_part_d,
    load_payments,
    specialty_group,
)
from rxinc.drugpanel import DrugPanelConfig, simulate_drug_panel
from rxinc.estimators import triple_diff

P1, P2, NP, RPH = "1111111111", "2222222222", "3333333333", "4444444444"
SLOTS = [f"Name_of_Drug_or_Biological_or_Device_or_Medical_Supply_{i}" for i in range(1, 6)]
FOOD, SPEAKER = "Food and Beverage", "Compensation for services other than consulting"


def _part_d_csv(tmp_path):
    rows = []
    for y in (2019, 2020, 2021, 2022):
        rows += [
            (P1, "Victoza 2-Pak", 20, y, "Internal Medicine"),
            (P1, "Victoza 3-Pak", 30, y, "Internal Medicine"),
            (P1, "Jardiance", 100, y, "Internal Medicine"),
            (P2, "Jardiance", 40, y, "Family Practice"),
            (P1, "Soliqua 100-33", 50, y, "Internal Medicine"),
            (NP, "Jardiance", 60, y, "Nurse Practitioner"),
            (NP, "Ozempic", 30, y, "Nurse Practitioner"),
            (RPH, "Jardiance", 25, y, "Pharmacist"),
        ]
    rows.append((P2, "Victoza 3-Pak", 15, 2020, "Family Practice"))
    df = pd.DataFrame(rows, columns=["Prscrbr_NPI", "Brnd_Name", "Tot_Clms", "data_year", "Prscrbr_Type"])
    df["Prscrbr_City"] = "LAFAYETTE"
    df["Prscrbr_State_Abrvtn"] = "IN"
    path = tmp_path / "pd.csv"
    df.to_csv(path, index=False)
    return path


def _payments_csvs(tmp_path):
    cols = ["Covered_Recipient_NPI", "Total_Amount_of_Payment_USDollars",
            "Nature_of_Payment_or_Transfer_of_Value", *SLOTS]
    def frame(rows):
        return pd.DataFrame([r + [None] * (len(cols) - len(r)) for r in rows], columns=cols)
    by_year = {
        2019: frame([[P2, "30", FOOD, "VICTOZA"]]),
        # $40 lunch naming two products: $20 each. P1 paid about Jardiance.
        2020: frame([[P1, "40", FOOD, "JARDIANCE", "SOMEOTHERDRUG"]]),
        # Two records for P1-Jardiance, one naming it twice; a speaker fee for P1-Victoza.
        2021: frame([[P1, "10", FOOD, "JARDIANCE"], [P1, "10", FOOD, "JARDIANCE", "JARDIANCE"],
                     [P1, "900", SPEAKER, "VICTOZA"], [NP, "15", FOOD, "OZEMPIC"]]),
        2022: frame([[NP, "15", FOOD, "JARDIANCE"]]),
    }
    out = {}
    for y, f in by_year.items():
        path = tmp_path / f"op{y}.csv"
        f.to_csv(path, index=False)
        out[y] = path
    return out


@pytest.fixture
def inputs(tmp_path):
    claims, attrs = load_part_d([_part_d_csv(tmp_path)])
    return claims, attrs, load_payments(_payments_csvs(tmp_path))


@pytest.fixture
def physicians(inputs):
    return build_drug_panel(*inputs, min_family_peak_claims=0, population="physician")


def _cell(panel, npi, fam, year):
    return panel[(panel.npi == npi) & (panel.family == fam) & (panel.year == year)].iloc[0]


@pytest.mark.parametrize("ptype, cls", [
    ("Internal Medicine", "physician"), ("Endocrinology", "physician"),
    ("Podiatry", "physician"), ("Nurse Practitioner", "npp"),
    ("Physician Assistant", "npp"), ("Certified Clinical Nurse Specialist", "npp"),
    ("Pharmacist", "not_covered"),
    ("Student in an Organized Health Care Education/Training Program", "not_covered"),
    ("Registered Nurse", "not_covered"), (None, "physician"),
])
def test_coverage_class(ptype, cls):
    assert coverage_class(ptype) == cls


def test_specialty_groups():
    assert specialty_group("Family Practice") == "Family or general practice"
    assert specialty_group("Interventional Cardiology") == "Cardiology"
    assert specialty_group("Rheumatology") == "Other"


def test_loader_refuses_missing_payment_nature(tmp_path):
    bad = tmp_path / "bad.csv"
    pd.DataFrame({"Covered_Recipient_NPI": [P1], "Total_Amount_of_Payment_USDollars": ["5"],
                  SLOTS[0]: ["JARDIANCE"]}).to_csv(bad, index=False)
    with pytest.raises(ValueError, match="Nature_of_Payment"):
        load_payments({2021: bad})


def test_pack_sizes_are_summed_into_one_family(inputs):
    claims = inputs[0]
    v = claims[(claims.npi == P1) & (claims.family == "VICTOZA")]
    assert (v["claims"] == 50).all()


def test_physician_population_excludes_npp_and_uncovered(physicians):
    panel, report = physicians
    assert set(panel.npi) == {P1, P2}
    assert report.excluded_by_coverage == {"npp": 1, "not_covered": 1}
    assert "SOLIQUA" not in set(panel.family)


def test_suppressed_cells_flagged(physicians):
    panel, _ = physicians
    c19, c20 = _cell(panel, P2, "VICTOZA", 2019), _cell(panel, P2, "VICTOZA", 2020)
    assert c19.any_rx == 0 and np.isnan(c19.log_rx)
    assert c20.any_rx == 1 and c20.log_rx == pytest.approx(np.log(15))


def test_amount_split_and_duplicate_naming(physicians):
    panel, _ = physicians
    assert _cell(panel, P1, "JARDIANCE", 2020).pay_amount == pytest.approx(20.0)
    c = _cell(panel, P1, "JARDIANCE", 2021)
    assert c.pay_records == 2 and c.pay_amount == pytest.approx(15.0)


def test_food_versus_speaker_split(physicians):
    panel, _ = physicians
    v = _cell(panel, P1, "VICTOZA", 2021)
    assert v.pay_nonfood_any == 1 and v.pay_food_only == 0 and v.pay_amount_nonfood == pytest.approx(900)
    j = _cell(panel, P1, "JARDIANCE", 2021)
    assert j.pay_nonfood_any == 0 and j.pay_food_only == 1


def test_dollar_thresholds(physicians):
    panel, _ = physicians
    j = _cell(panel, P1, "JARDIANCE", 2021)
    assert j.pay_ge25 == 0 and j.pay_ge100 == 0
    v = _cell(panel, P1, "VICTOZA", 2021)
    assert v.pay_ge25 == 1 and v.pay_ge100 == 1


def test_same_manufacturer_other_family():
    claims = pd.DataFrame({"npi": [P1] * 2, "family": ["JARDIANCE", "TRADJENTA"],
                           "year": [2019] * 2, "claims": [50.0, 40.0]})
    attrs = pd.DataFrame({"npi": [P1], "year": [2019], "specialty": ["Internal Medicine"],
                          "city": ["X"], "state": ["IN"]})
    pay = pd.DataFrame({"npi": [P1], "family": ["TRADJENTA"], "year": [2019], "pay_amount": [20.0],
                        "pay_records": [1], "pay_amount_nonfood": [0.0], "pay_records_nonfood": [0]})
    panel, _ = build_drug_panel(claims, attrs, pay, min_family_peak_claims=0)
    # Jardiance and Tradjenta are both Boehringer/Lilly.
    assert _cell(panel, P1, "JARDIANCE", 2019).pay_same_mfr_other == 1
    assert _cell(panel, P1, "TRADJENTA", 2019).pay_same_mfr_other == 0


def test_physician_left_censoring_and_clean_onset(physicians):
    panel, report = physicians
    p2v = panel[(panel.npi == P2) & (panel.family == "VICTOZA")]
    assert p2v.left_censored.all() and p2v.event_time.isna().all()
    p1j = panel[(panel.npi == P1) & (panel.family == "JARDIANCE")].set_index("year")
    assert p1j.loc[2019, "event_time"] == -1 and p1j.loc[2020, "treated"] == 1
    assert report.n_left_censored_pairs == 1


def test_npp_exposure_missing_before_2021_and_censoring_moves(inputs):
    panel, report = build_drug_panel(*inputs, min_family_peak_claims=0, population="npp")
    assert set(panel.npi) == {NP}
    oz = panel[panel.family == "OZEMPIC"].set_index("year")
    assert oz.loc[[2019, 2020], "pay_any"].isna().all()
    assert not oz.loc[2020, "exposure_known"] and oz.loc[2021, "exposure_known"]
    # Paid in 2021, the first observable year, so onset is unknown.
    assert oz["left_censored"].all() and oz["event_time"].isna().all()
    # Paid about Jardiance in 2022 after an observed unpaid 2021: a clean onset.
    ja = panel[panel.family == "JARDIANCE"].set_index("year")
    assert ja.loc[2022, "event_time"] == 0 and ja.loc[2021, "event_time"] == -1
    assert report.n_clean_onset_pairs == 1


def test_lag_and_lead_align_within_pair(physicians):
    panel, _ = physicians
    p1j = panel[(panel.npi == P1) & (panel.family == "JARDIANCE")].set_index("year")
    assert np.isnan(p1j.loc[2019, "pay_any_lag"]) and p1j.loc[2021, "pay_any_lag"] == 1
    assert p1j.loc[2019, "pay_any_lead"] == 1 and np.isnan(p1j.loc[2022, "pay_any_lead"])


def test_peer_share_is_leave_one_out(physicians):
    panel, _ = physicians
    v19 = panel[(panel.family == "VICTOZA") & (panel.year == 2019)].set_index("npi")
    assert v19.loc[P1, "peer_pay_share"] == 1.0
    assert v19.loc[P2, "peer_pay_share"] == 0.0


def test_observed_only_grid_is_smaller_and_consistent(inputs):
    full, _ = build_drug_panel(*inputs, min_family_peak_claims=0, full_grid=True)
    obs, _ = build_drug_panel(*inputs, min_family_peak_claims=0, full_grid=False)
    assert len(obs) == int(full.any_rx.sum()) and obs.any_rx.all()


def test_rejects_unknown_population(inputs):
    with pytest.raises(ValueError, match="population must be one of"):
        build_drug_panel(*inputs, population="dentists")


def test_controls_reported_and_do_not_disturb_treatment_on_simulation():
    df = simulate_drug_panel(DrugPanelConfig(n_physicians=400))
    df["noise_ctrl"] = np.random.default_rng(0).normal(size=len(df))
    base, with_ctrl = triple_diff(df), triple_diff(df, controls=["noise_ctrl"])
    assert "noise_ctrl" in with_ctrl.extra and "noise_ctrl_se" in with_ctrl.extra
    assert with_ctrl.coef == pytest.approx(base.coef, abs=0.01)


def test_pair_always_observed(physicians):
    panel, _ = physicians
    # P1-Jardiance clears the floor every year; P2-Victoza only in 2020.
    assert panel[(panel.npi == P1) & (panel.family == "JARDIANCE")].pair_always_observed.all()
    assert not panel[(panel.npi == P2) & (panel.family == "VICTOZA")].pair_always_observed.any()
