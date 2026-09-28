# Data sources

Both datasets are public, free, require no API key, and contain **no protected
health information**. They are provider-level aggregates, not patient records.

## Medicare Part D Prescribers

CMS publishes prescriber-level Part D data at `data.cms.gov`. The package
discovers endpoints from the portal's own DCAT catalog rather than hard-coding
UUIDs, because CMS reissues identifiers when it republishes a year.

- **By Provider** — one row per prescriber-year with totals. ~1.42M rows per
  year.
- **By Provider and Drug** — one row per prescriber-drug-year.

Query with `filter[Column]=value`, `size` and `offset`. A stats endpoint
returns matching row counts before you commit to a download.

```bash
rxinc fetch-partd --state IN --max-rows 5000 --out data/raw/part_d_in.csv
```

Key columns: `Prscrbr_NPI`, `Prscrbr_Last_Org_Name`, `Prscrbr_First_Name`,
`Prscrbr_State_Abrvtn`, `Prscrbr_Type`, `Tot_Clms`.

**What it is not.** Part D covers Medicare beneficiaries, so it
under-represents younger patients and commercially insured populations. Counts
below 11 are suppressed, which is why `build_panel` defaults `min_claims=11`.

**Why counts matter.** Commercial prescriber panels project national estimates
from a pharmacy sample — roughly 93% of retail and up to 77% of mail and
long-term-care — using a geospatial method that assumes unsampled pharmacies
resemble nearby sampled ones. Validation work against Part D found projection
errors around −8.4% and −5.6% per physician for female and male patients in
rural areas. Part D claim counts are actual adjudicated claims, which is
precisely why the literature uses them as the benchmark.

## Open Payments

The Sunshine Act registry of manufacturer payments to clinicians, at
`openpaymentsdata.cms.gov`. Program Year 2025 published 17.07 million records
totalling $14.67 billion — the highest annual total since collection began in
2013.

General payments (non-research, non-ownership) are the promotional
touchpoints: meals, travel, speaking and consulting fees. Distributed as bulk
CSV, one file per program year, several gigabytes each.

```python
from rxinc.datasets import download_open_payments, load_open_payments_csv

path = download_open_payments(2023)
payments = load_open_payments_csv(path)  # keeps 10 of 91 columns
```

### NPI is present — a correction worth stating

Secondary sources frequently claim Open Payments contains no NPI and that
linkage to Part D therefore requires name-and-address matching. **That is
outdated.** The current detailed general-payments schema includes
`Covered_Recipient_NPI` alongside `Covered_Recipient_Profile_ID`, verified
against the published data dictionary for the 2023 program year (91 fields).

The claim was true of early program years, and much published work on this
linkage was done under that constraint. `rxinc.linkage` therefore keys on NPI
where present and valid, and keeps a surname-plus-state fallback for older
years and blank fields — reporting how many records took each route, because a
linkage leaning heavily on the fallback is not interchangeable with one that
keyed cleanly.

NPIs are validated against their Luhn check digit (computed over the constant
prefix `80840` plus the first nine digits), which rejects transcription errors
and placeholders like `0000000000` that would otherwise produce confident
wrong joins.

## A note on the network client

The Open Payments portal answers HTTP 403 to some custom `User-Agent` strings.
Requests go out with the standard library default; override via
`rxinc.datasets.USER_AGENT` only if you have reason to.

## Rate and volume

Part D queries page at 5,000 rows. `fetch_part_d` caps at 50,000 rows by
default so an accidental unfiltered call does not pull 1.4 million. Open
Payments downloads stream in 1 MB chunks and skip work if the file is present.
