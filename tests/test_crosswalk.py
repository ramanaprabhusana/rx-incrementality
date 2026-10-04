"""Tests for brand-family normalisation, the stage most likely to fail silently."""

from __future__ import annotations

import pytest

from rxinc.crosswalk import (
    FAMILY_CLASS,
    FAMILY_MANUFACTURER,
    brand_family,
    crosswalk_report,
    in_class,
)


@pytest.mark.parametrize(
    "raw, family",
    [
        ("Victoza 2-Pak", "VICTOZA"),
        ("Victoza 3-Pak", "VICTOZA"),
        ("VICTOZA", "VICTOZA"),
        ("Bydureon Pen", "BYDUREON"),
        ("Bydureon Bcise", "BYDUREON"),
        ("BYDUREON BCISE", "BYDUREON"),
        ("Janumet Xr", "JANUMET"),
        ("JANUMET", "JANUMET"),
        ("Xigduo Xr", "XIGDUO"),
        ("Kombiglyze Xr", "KOMBIGLYZE"),
        ("Synjardy Xr", "SYNJARDY"),
        ("Soliqua 100-33", "SOLIQUA"),
        ("Xultophy 100-3.6", "XULTOPHY"),
        ("SOLIQUA 100/33", "SOLIQUA"),
        ("XULTOPHY 100/3.6", "XULTOPHY"),
        ("TRIJARDY XR", "TRIJARDY"),
        ("Trijardy Xr", "TRIJARDY"),
        ("XIGDUO", "XIGDUO"),
        ("  jardiance ", "JARDIANCE"),
    ],
)
def test_brand_family_collapses_variants(raw, family):
    assert brand_family(raw) == family


@pytest.mark.parametrize("raw", [None, float("nan"), "", "   "])
def test_brand_family_missing_is_none(raw):
    assert brand_family(raw) is None


def test_victoza_regression():
    """The bug that motivated this module: Part D never says plain 'Victoza'."""
    assert brand_family("Victoza 3-Pak") == brand_family("VICTOZA")


def test_every_class_family_has_a_manufacturer():
    assert set(FAMILY_CLASS) == set(FAMILY_MANUFACTURER)


def test_exclusions_are_not_in_class():
    for fam in ("SOLIQUA", "XULTOPHY", "WEGOVY"):
        assert not in_class(fam)


def test_report_flags_one_sided_families():
    report = crosswalk_report(
        part_d_names=["Jardiance", "Victoza 3-Pak", "Tradjenta", "Soliqua 100-33"],
        payment_names=["JARDIANCE", "VICTOZA", "OZEMPIC", "WIDGETMAB"],
    )
    status = dict(zip(report["raw_name"], report["status"]))
    assert status["Jardiance"] == "matched"
    assert status["Victoza 3-Pak"] == "matched"
    assert status["Tradjenta"] == "part_d_only"
    assert status["OZEMPIC"] == "payments_only"
    assert status["Soliqua 100-33"] == "excluded"
    assert status["WIDGETMAB"] == "out_of_class"
