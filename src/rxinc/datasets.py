"""Clients for the two public CMS datasets this analysis runs on.

Both are free, require no API key, and contain no protected health
information: they are provider-level aggregates, not patient records.

Medicare Part D Prescribers (data.cms.gov)
    Actual dispensed-claim counts by prescriber, by year.  Unlike the
    commercial prescriber panels used in industry, these are counts rather
    than projections from a pharmacy sample, which is why validation studies
    use Part D as the benchmark.  Limited to Medicare Part D beneficiaries, so
    it under-represents younger patients.

Open Payments (openpaymentsdata.cms.gov)
    Every reportable payment or transfer of value from drug and device
    manufacturers to clinicians, under the Sunshine Act.  Program Year 2025
    published 17.07 million records totalling $14.67 billion.

Endpoints are discovered from each portal's own catalog rather than
hard-coded, because CMS reissues dataset UUIDs when it republishes a year.
"""

from __future__ import annotations

import json
import urllib.parse
import urllib.request
from collections.abc import Iterator
from dataclasses import dataclass
from typing import Any

import pandas as pd

CMS_CATALOG_URL = "https://data.cms.gov/data.json"
OPEN_PAYMENTS_METASTORE = (
    "https://openpaymentsdata.cms.gov/api/1/metastore/schemas/dataset/items"
)
PART_D_BY_PROVIDER_TITLE = "Medicare Part D Prescribers - by Provider"
PART_D_BY_PROVIDER_DRUG_TITLE = "Medicare Part D Prescribers - by Provider and Drug"

DEFAULT_TIMEOUT = 120
_PAGE_SIZE = 5000

USER_AGENT: str | None = None
"""Optional User-Agent override for CMS requests.

Left as ``None``, requests go out with the standard library default. That is
deliberate: the Open Payments portal sits behind a filter that rejects some
custom agent strings with HTTP 403, and the default is accepted. Set this only
if you have a reason to identify your client differently.
"""


def _headers(accept_json: bool = True) -> dict[str, str]:
    """Build request headers, overriding User-Agent only if configured."""
    headers = {"Accept": "application/json"} if accept_json else {}
    if USER_AGENT:
        headers["User-Agent"] = USER_AGENT
    return headers


@dataclass(frozen=True)
class Distribution:
    """One published year of a CMS dataset.

    Attributes:
        title: Distribution title as CMS publishes it.
        api_url: JSON API endpoint supporting ``size``/``offset``/``filter``.
        download_url: Bulk CSV URL, when the portal exposes one.
        year: Data year, when it can be determined.
    """

    title: str
    api_url: str | None
    download_url: str | None
    year: int | None = None


def _get_json(url: str, timeout: int = DEFAULT_TIMEOUT) -> Any:
    """Fetch and parse JSON from ``url``.

    Args:
        url: Absolute URL.
        timeout: Socket timeout in seconds.

    Returns:
        Parsed JSON.

    Raises:
        urllib.error.URLError: On network failure.
        urllib.error.HTTPError: On a non-2xx response. The Open Payments
            portal answers 403 to some custom User-Agent strings; see
            :data:`USER_AGENT`.
    """
    request = urllib.request.Request(url, headers=_headers())
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return json.loads(response.read().decode("utf-8"))


def cms_catalog(timeout: int = DEFAULT_TIMEOUT) -> list[dict[str, Any]]:
    """Return the data.cms.gov DCAT catalog entries.

    Args:
        timeout: Socket timeout in seconds.

    Returns:
        List of dataset records.
    """
    return _get_json(CMS_CATALOG_URL, timeout=timeout).get("dataset", [])


def _record_year(record: dict[str, Any], title: str) -> int | None:
    """Data year of a catalog record, from ``temporal`` or the title suffix."""
    temporal = record.get("temporal")
    if isinstance(temporal, list) and temporal:
        temporal = temporal[0]
    if isinstance(temporal, dict) and temporal.get("endDate"):
        return int(str(temporal["endDate"])[:4])
    if isinstance(temporal, str) and "/" in temporal:
        # Old layout spans every year in one record, so it says nothing
        # about an individual distribution.
        pass
    tail = title.rsplit(":", 1)[-1].strip()
    return int(tail[:4]) if tail[:4].isdigit() else None


def part_d_distributions(
    title: str = PART_D_BY_PROVIDER_TITLE, timeout: int = DEFAULT_TIMEOUT
) -> list[Distribution]:
    """List published years of a Part D Prescribers dataset.

    CMS has published this catalog in two layouts, and switched between them
    without notice. Originally one record held every year as separate
    distributions. Later each year became its own record titled
    ``"<title> : YYYY-MM-DD"``. Both are handled, so a catalog restructure
    does not silently break acquisition.

    Args:
        title: Base catalog title, without any year suffix. Use
            :data:`PART_D_BY_PROVIDER_TITLE` for prescriber totals or
            :data:`PART_D_BY_PROVIDER_DRUG_TITLE` for per-drug detail.
        timeout: Socket timeout in seconds.

    Returns:
        One :class:`Distribution` per published distribution, newest year
        first.

    Raises:
        LookupError: If no catalog entry matches ``title`` in either layout.
    """
    out: list[Distribution] = []
    for record in cms_catalog(timeout=timeout):
        record_title = record.get("title", "")
        if record_title != title and not record_title.startswith(f"{title} :"):
            continue
        per_year_layout = record_title != title
        for dist in record.get("distribution") or []:
            dist_title = dist.get("title") or record_title
            access = dist.get("accessURL")
            year = (
                _record_year(record, record_title)
                if per_year_layout
                else _record_year({}, dist_title)
            )
            out.append(
                Distribution(
                    title=dist_title,
                    api_url=access if access and access.endswith("/data") else None,
                    download_url=dist.get("downloadURL"),
                    year=year,
                )
            )
    if not out:
        raise LookupError(
            f"No CMS catalog entry titled {title!r}, in either the combined or "
            "per-year layout. List cms_catalog() to find the current title."
        )
    out.sort(key=lambda d: d.year or 0, reverse=True)
    return out


def part_d_api_by_year(
    title: str = PART_D_BY_PROVIDER_DRUG_TITLE, timeout: int = DEFAULT_TIMEOUT
) -> dict[int, str]:
    """Map data year to its queryable API URL.

    Args:
        title: Base catalog title.
        timeout: Socket timeout in seconds.

    Returns:
        ``{year: api_url}`` for every year with a JSON API distribution.
    """
    out: dict[int, str] = {}
    for dist in part_d_distributions(title, timeout=timeout):
        if dist.api_url and dist.year is not None:
            out.setdefault(dist.year, dist.api_url)
    return out


def part_d_row_count(api_url: str, filters: dict[str, str] | None = None) -> int:
    """Return how many rows a filtered Part D query will yield.

    Worth calling before a download: an unfiltered year is ~1.4 million rows.

    Args:
        api_url: Distribution API URL ending in ``/data``.
        filters: Column-name to exact-value filters.

    Returns:
        Number of matching rows.
    """
    stats_url = api_url.rsplit("/data", 1)[0] + "/data-viewer/stats"
    if filters:
        query = urllib.parse.urlencode(
            {f"filter[{k}]": v for k, v in filters.items()}
        )
        stats_url = f"{stats_url}?{query}"
    payload = _get_json(stats_url)
    return int(payload["data"]["found_rows"])


def iter_part_d(
    api_url: str,
    filters: dict[str, str] | None = None,
    max_rows: int | None = None,
    page_size: int = _PAGE_SIZE,
) -> Iterator[dict[str, Any]]:
    """Yield Part D prescriber rows, paging through the API.

    Args:
        api_url: Distribution API URL ending in ``/data``.
        filters: Column-name to exact-value filters, e.g.
            ``{"Prscrbr_State_Abrvtn": "IN"}``.
        max_rows: Stop after this many rows.  ``None`` fetches everything
            matching, which for an unfiltered year is ~1.4 million rows.
        page_size: Rows per request.

    Yields:
        One dict per prescriber-year row.
    """
    base = {f"filter[{k}]": v for k, v in (filters or {}).items()}
    fetched = 0
    offset = 0
    while True:
        remaining = None if max_rows is None else max_rows - fetched
        if remaining is not None and remaining <= 0:
            return
        size = page_size if remaining is None else min(page_size, remaining)
        query = urllib.parse.urlencode({**base, "size": size, "offset": offset})
        rows = _get_json(f"{api_url}?{query}")
        if not rows:
            return
        yield from rows
        fetched += len(rows)
        offset += len(rows)
        if len(rows) < size:
            return


def fetch_part_d(
    api_url: str,
    filters: dict[str, str] | None = None,
    max_rows: int | None = 50_000,
) -> pd.DataFrame:
    """Fetch Part D prescriber rows into a DataFrame.

    Args:
        api_url: Distribution API URL ending in ``/data``.
        filters: Column-name to exact-value filters.
        max_rows: Row cap.  Defaults to 50,000 so an accidental unfiltered
            call does not pull a million rows.

    Returns:
        DataFrame of the requested rows.  Key columns are ``Prscrbr_NPI``,
        ``Prscrbr_Last_Org_Name``, ``Prscrbr_First_Name``,
        ``Prscrbr_State_Abrvtn``, ``Prscrbr_Type`` and ``Tot_Clms``.
    """
    rows = list(iter_part_d(api_url, filters=filters, max_rows=max_rows))
    return pd.DataFrame(rows)


def open_payments_catalog(timeout: int = DEFAULT_TIMEOUT) -> list[dict[str, Any]]:
    """Return Open Payments metastore dataset records.

    Args:
        timeout: Socket timeout in seconds.

    Returns:
        List of dataset records.
    """
    return _get_json(OPEN_PAYMENTS_METASTORE, timeout=timeout)


def open_payments_general_distribution(
    year: int, timeout: int = DEFAULT_TIMEOUT
) -> Distribution:
    """Locate the general-payments CSV for a program year.

    General payments are the non-research, non-ownership transfers: meals,
    travel, speaking and consulting fees.  These are the promotional
    touchpoints relevant to prescribing.

    Args:
        year: Program year, e.g. ``2023``.
        timeout: Socket timeout in seconds.

    Returns:
        The distribution, whose ``download_url`` is a bulk CSV.

    Raises:
        LookupError: If that program year is not published.
    """
    wanted = f"{year} General Payment Data"
    for record in open_payments_catalog(timeout=timeout):
        if record.get("title") == wanted:
            for dist in record.get("distribution") or []:
                data = dist.get("data", dist)
                url = data.get("downloadURL")
                if url:
                    return Distribution(
                        title=data.get("title", wanted),
                        api_url=None,
                        download_url=url,
                    )
    raise LookupError(f"No Open Payments general payment dataset for {year}.")
