"""Stream a year of Open Payments and keep only rows naming a target drug.

The general payments file for a single program year runs 6 to 8.6 GB, and the
query API caps at 500 rows per request at roughly 18 rows/sec, which makes it
unusable at this scale. So this streams the CSV over HTTP, filters in flight,
and never stores the raw file. Peak disk is the filtered output only.

Matched chunks are appended as they are parsed, so an interrupted run leaves
usable partial output rather than nothing.

    python3 scripts/stream_open_payments.py --year 2023 --class diabetes
"""

from __future__ import annotations

import argparse
import sys
import time
import urllib.request
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from rxinc.crosswalk import brand_family, in_class  # noqa: E402
from rxinc.datasets import open_payments_general_distribution  # noqa: E402

PRODUCT_SLOTS = range(1, 6)
PRODUCT_COL = "Name_of_Drug_or_Biological_or_Device_or_Medical_Supply_{}"
NDC_COL = "Associated_Drug_or_Biological_NDC_{}"
CATEGORY_COL = "Product_Category_or_Therapeutic_Area_{}"

KEEP_COLUMNS = [
    "Covered_Recipient_Type",
    "Covered_Recipient_NPI",
    "Covered_Recipient_First_Name",
    "Covered_Recipient_Last_Name",
    "Recipient_State",
    "Applicable_Manufacturer_or_Applicable_GPO_Making_Payment_Name",
    "Total_Amount_of_Payment_USDollars",
    "Date_of_Payment",
    "Number_of_Payments_Included_in_Total_Amount",
    "Nature_of_Payment_or_Transfer_of_Value",
    *[PRODUCT_COL.format(i) for i in PRODUCT_SLOTS],
    *[NDC_COL.format(i) for i in PRODUCT_SLOTS],
    *[CATEGORY_COL.format(i) for i in PRODUCT_SLOTS],
]


class _CountingReader:
    """Wrap a stream and count bytes read, so progress can be reported."""

    def __init__(self, stream, total: int | None):
        self._stream = stream
        self.total = total
        self.bytes_read = 0

    def read(self, size: int = -1) -> bytes:
        chunk = self._stream.read(size)
        self.bytes_read += len(chunk)
        return chunk

    @property
    def pct(self) -> float:
        if not self.total:
            return float("nan")
        return 100.0 * self.bytes_read / self.total


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--year", type=int, required=True)
    parser.add_argument("--class", dest="klass", default="diabetes", choices=["diabetes"])
    parser.add_argument("--chunksize", type=int, default=400_000)
    parser.add_argument("--out-dir", default="data/raw")
    args = parser.parse_args()

    dist = open_payments_general_distribution(args.year)
    url = dist.download_url
    assert url

    out = Path(args.out_dir) / f"openpayments_{args.klass}_{args.year}.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists():
        out.unlink()

    print(f"year={args.year} class={args.klass} matching=brand-family", flush=True)
    print(f"streaming {url.split('/')[-1]}", flush=True)

    request = urllib.request.Request(url)
    started = time.time()
    scanned = matched_total = 0
    wrote_header = False
    family_cache: dict[str, bool] = {}

    with urllib.request.urlopen(request, timeout=300) as response:
        total = response.headers.get("Content-Length")
        reader = _CountingReader(response, int(total) if total else None)

        reader_iter = pd.read_csv(
            reader,
            chunksize=args.chunksize,
            usecols=lambda c: c in set(KEEP_COLUMNS),
            dtype=str,
            low_memory=False,
            on_bad_lines="warn",
        )

        for chunk_no, chunk in enumerate(reader_iter, start=1):
            scanned += len(chunk)

            # Match on brand family, not exact names. Exact matching missed
            # XIGDUO (listed as XIGDUO XR), TRIJARDY XR, STEGLUJAN, SEGLUROMET
            # and variant spellings. Family lookups are memoised per distinct
            # string, so the cost is per unique name, not per row.
            mask = pd.Series(False, index=chunk.index)
            for slot in PRODUCT_SLOTS:
                col = PRODUCT_COL.format(slot)
                if col in chunk.columns:
                    values = chunk[col]
                    for name in values.dropna().unique():
                        if name not in family_cache:
                            family_cache[name] = in_class(brand_family(name))
                    keep = {n for n in values.dropna().unique() if family_cache[n]}
                    mask |= values.isin(keep)

            hits = chunk[mask]
            if len(hits):
                hits.to_csv(out, mode="a", index=False, header=not wrote_header)
                wrote_header = True
                matched_total += len(hits)

            if chunk_no % 5 == 0 or reader.pct >= 99:
                elapsed = time.time() - started
                print(
                    f"  chunk {chunk_no:>3}  {reader.pct:5.1f}%  "
                    f"scanned {scanned:>10,}  matched {matched_total:>8,}  "
                    f"{reader.bytes_read/2**20/max(elapsed,1):5.1f} MB/s  "
                    f"{elapsed/60:5.1f} min",
                    flush=True,
                )

    elapsed = time.time() - started
    print(
        f"\nDONE year={args.year}: scanned {scanned:,} rows, "
        f"matched {matched_total:,} ({matched_total/max(scanned,1):.2%}) "
        f"in {elapsed/60:.1f} min -> {out.name}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
