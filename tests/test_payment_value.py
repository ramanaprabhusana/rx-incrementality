"""The per-dollar translation, on cells small enough to check by hand."""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "estimate_real.py"


@pytest.fixture(scope="module")
def payment_value():
    spec = importlib.util.spec_from_file_location("estimate_real", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.payment_value


def _cells() -> pd.DataFrame:
    return pd.DataFrame({
        # paid and observed; paid below the floor; unpaid; exposure not yet observed
        "year": [2021, 2021, 2021, 2020],
        "exposure_known": [True, True, True, False],
        "any_rx": [1, 0, 1, 1],
        "pay_any": [1.0, 1.0, 0.0, np.nan],
        "pay_amount": [10.0, 5.0, 0.0, np.nan],
        "cost": [1000.0, 0.0, 500.0, 700.0],
    })


def test_caused_spending_counts_only_paid_observed_cells(payment_value):
    tau = math.log(1.05)
    v = payment_value(_cells(), {"headline": tau})
    assert v["payment_dollars"] == pytest.approx(15.0)       # includes the cell below the floor
    assert v["payment_dollars_in_observed_cells"] == pytest.approx(10.0)
    assert v["gross_spending_paid_observed"] == pytest.approx(1000.0)
    assert v["naive_per_payment_dollar"] == pytest.approx(1000.0 / 15.0)
    caused = 1000.0 * (1 - 1 / 1.05)                          # spending that would not have happened
    assert v["causal"]["headline"]["caused_gross_spending"] == pytest.approx(caused)
    assert v["causal"]["headline"]["per_payment_dollar"] == pytest.approx(caused / 15.0)
    assert v["naive_by_year"] == {2021: pytest.approx(1000.0 / 15.0)}


def test_no_effect_means_no_caused_spending(payment_value):
    v = payment_value(_cells(), {"none": 0.0})
    assert v["causal"]["none"]["per_payment_dollar"] == 0.0
