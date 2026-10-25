"""Tests for CMS catalog parsing, which has changed layout without notice."""

from __future__ import annotations

import json

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


class _FakeResponse:
    def __init__(self, payload):
        self._payload = payload

    def read(self):
        return json.dumps(self._payload).encode()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def _http_error(code, retry_after=None):
    import email.message
    import urllib.error
    headers = email.message.Message()
    if retry_after is not None:
        headers["Retry-After"] = str(retry_after)
    return urllib.error.HTTPError("https://x", code, "err", headers, None)


def test_get_json_retries_rate_limits_then_succeeds(monkeypatch):
    calls, waits = [], []
    outcomes = [_http_error(429, retry_after=3), _http_error(503), [{"ok": 1}]]

    def fake_urlopen(request, timeout=None):
        calls.append(request)
        result = outcomes.pop(0)
        if isinstance(result, Exception):
            raise result
        return _FakeResponse(result)

    monkeypatch.setattr(datasets.urllib.request, "urlopen", fake_urlopen)
    monkeypatch.setattr(datasets.time, "sleep", waits.append)
    assert datasets._get_json("https://x", backoff=2.0) == [{"ok": 1}]
    assert len(calls) == 3
    assert waits == [3.0, 4.0]   # the server's Retry-After, then backoff * 2**1


def test_get_json_does_not_retry_client_errors(monkeypatch):
    import urllib.error
    calls = []

    def fake_urlopen(request, timeout=None):
        calls.append(request)
        raise _http_error(403)

    monkeypatch.setattr(datasets.urllib.request, "urlopen", fake_urlopen)
    monkeypatch.setattr(datasets.time, "sleep", lambda s: None)
    with pytest.raises(urllib.error.HTTPError):
        datasets._get_json("https://x")
    assert len(calls) == 1


def test_get_json_gives_up_after_the_last_retry(monkeypatch):
    import urllib.error
    calls = []

    def fake_urlopen(request, timeout=None):
        calls.append(request)
        raise _http_error(429)

    monkeypatch.setattr(datasets.urllib.request, "urlopen", fake_urlopen)
    monkeypatch.setattr(datasets.time, "sleep", lambda s: None)
    with pytest.raises(urllib.error.HTTPError):
        datasets._get_json("https://x", retries=2)
    assert len(calls) == 3
