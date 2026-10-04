"""Tests for CMS catalog parsing, which has changed layout without notice."""

from __future__ import annotations

import pytest

from rxinc import datasets

BASE = datasets.PART_D_BY_PROVIDER_DRUG_TITLE


def _combined_layout():
    """Original layout: one record, one distribution per year."""
    return [
        {
            "title": BASE,
            "temporal": "2013-01-01/2024-12-31",
            "distribution": [
                {"title": f"{BASE} : 2024-12-01",
                 "accessURL": "https://x/api/aaa/data"},
                {"title": f"{BASE} : 2024-12-01",
                 "downloadURL": "https://x/DY24.csv"},
                {"title": f"{BASE} : 2019-12-31",
                 "accessURL": "https://x/api/bbb/data"},
            ],
        }
    ]


def _per_year_layout():
    """Current layout: one record per data year."""
    return [
        {
            "title": f"{BASE} : 2019-12-31",
            "temporal": [{"startDate": "2019-01-01", "endDate": "2019-12-31"}],
            "distribution": [
                {"title": f"{BASE} : 2019-12-31",
                 "downloadURL": "https://x/DY19.csv"},
                {"title": f"{BASE} : 2019-12-31",
                 "accessURL": "https://x/api/bbb/data"},
            ],
        },
        {
            "title": f"{BASE} : 2024-12-01",
            "temporal": [{"startDate": "2024-01-01", "endDate": "2024-12-31"}],
            "distribution": [
                {"title": f"{BASE} : 2024-12-01",
                 "accessURL": "https://x/api/aaa/data"},
            ],
        },
        {"title": "Some Unrelated Dataset", "distribution": []},
    ]


@pytest.mark.parametrize("layout", [_combined_layout, _per_year_layout])
def test_api_by_year_handles_both_catalog_layouts(monkeypatch, layout):
    monkeypatch.setattr(datasets, "cms_catalog", lambda timeout=0: layout())
    by_year = datasets.part_d_api_by_year(BASE)
    assert by_year == {2024: "https://x/api/aaa/data", 2019: "https://x/api/bbb/data"}


def test_distributions_sorted_newest_first(monkeypatch):
    monkeypatch.setattr(datasets, "cms_catalog", lambda timeout=0: _per_year_layout())
    years = [d.year for d in datasets.part_d_distributions(BASE)]
    assert years == sorted(years, reverse=True)


def test_prefix_match_does_not_capture_other_datasets(monkeypatch):
    catalog = _per_year_layout() + [
        {"title": f"{BASE} Supplement : 2020-12-31",
         "distribution": [{"accessURL": "https://x/api/zzz/data"}]}
    ]
    monkeypatch.setattr(datasets, "cms_catalog", lambda timeout=0: catalog)
    assert "https://x/api/zzz/data" not in datasets.part_d_api_by_year(BASE).values()


def test_missing_title_raises_with_guidance(monkeypatch):
    monkeypatch.setattr(datasets, "cms_catalog", lambda timeout=0: [])
    with pytest.raises(LookupError, match="either the combined or per-year layout"):
        datasets.part_d_distributions(BASE)
