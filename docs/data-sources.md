# Data sources

Both datasets are public, free, require no API key, and contain **no protected
health information**. They are provider-level aggregates, not patient records.

## Medicare Part D Prescribers

CMS publishes prescriber-level Part D data at `data.cms.gov`. The package
discovers endpoints from the portal's own DCAT catalog rather than hard-coding
UUIDs, because CMS reissues identifiers when it republishes a year.

- **By Provider**: one row per prescriber-year with totals. ~1.42M rows per
  year.
- **By Provider and Drug**: one row per prescriber-drug-year.

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
from a pharmacy sample (roughly 93% of retail and up to 77% of mail and
long-term-care) using a geospatial method that assumes unsampled pharmacies
resemble nearby sampled ones. Validation work against Part D found projection
errors around −8.4% and −5.6% per physician for female and male patients in
rural areas. Part D claim counts are actual adjudicated claims, which is
precisely why the literature uses them as the benchmark.

## Open Payments

The Sunshine Act registry of manufacturer payments to clinicians, at
`openpaymentsdata.cms.gov`. Program Year 2025 published 17.07 million records
totalling $14.67 billion, the highest annual total since collection began in
2013.

General payments (non-research, non-ownership) are the promotional
touchpoints: meals, travel, speaking and consulting fees. Distributed as bulk
CSV, one file per program year, several gigabytes each.

```python
from rxinc.datasets import download_open_payments, load_open_payments_csv

path = download_open_payments(2023)
payments = load_open_payments_csv(path)  # keeps 10 of 91 columns
```

### NPI is present: a correction worth stating

Secondary sources frequently claim Open Payments contains no NPI and that
linkage to Part D therefore requires name-and-address matching. **That is
outdated.** The current detailed general-payments schema includes
`Covered_Recipient_NPI` alongside `Covered_Recipient_Profile_ID`, verified
against the published data dictionary for the 2023 program year (91 fields).

The claim was true of early program years, and much published work on this
linkage was done under that constraint. `rxinc.linkage` therefore keys on NPI
where present and valid, and keeps a surname-plus-state fallback for older
years and blank fields, reporting how many records took each route, because a
linkage leaning heavily on the fallback is not interchangeable with one that
keyed cleanly.

NPIs are validated against their Luhn check digit (computed over the constant
prefix `80840` plus the first nine digits), which rejects transcription errors
and placeholders like `0000000000` that would otherwise produce confident
wrong joins.

## Joining the two sources: brand families

Open Payments reports what a payment was about at the brand level (`VICTOZA`,
`JANUMET`). Part D splits brands across pack sizes, devices and formulations
(`Victoza 2-Pak`, `Victoza 3-Pak`, `Janumet Xr`, `Bydureon Pen`). Joining raw
strings drops prescribing silently. The first version of this pipeline did
exactly that and lost Victoza entirely: 1.58M claims in 2019, more than Ozempic
that year.

`rxinc.crosswalk.brand_family` normalises both sides to a family, and
`crosswalk_report` lists every raw name against its family. Discovery was done
systematically, not by guessing names:

- Part D brands found by generic molecule in the national Geography and Drug
  file, using short stems, because CMS abbreviates triple-combination generics
  (`Empaglifloz/Linaglip/Metformin` for Trijardy Xr).
- Open Payments names found by streaming one full unfiltered year.

Result: 28 of 29 families match on both sides. Kombiglyze has prescribing but no
payments in 2019 to 2024, which is genuine; it was not promoted before
discontinuation. Insulin combinations, obesity indications and unbranded
generics are excluded by design and reported as such.

Open Payments extracts are filtered by family during streaming. Exact-name
filtering had missed XIGDUO (listed as XIGDUO XR), TRIJARDY XR, STEGLUJAN and
SEGLUROMET, about 1% of records.

## Who Open Payments covers, and when

| Prescriber type | Covered recipient | Records in the 2019 and 2020 extracts |
| --- | --- | ---: |
| Physicians, including dentists and podiatrists | from 2013 | 948,651 and 605,178 |
| Nurse practitioners, physician assistants, clinical nurse specialists, nurse midwives, CRNAs | from program year 2021 | 0 and 0 |
| Pharmacists, physicians in training, registered nurses | never | |

From 2021 non-physician practitioners receive about 360,000 to 390,000 records a
year, roughly 37% of the total. Before 2021 their payments were simply not
collected. A pipeline that treats a missing record as no payment will therefore
see spurious 2021 onsets for every nurse practitioner who was already being
visited. `rxinc.drugdata.coverage_class` classifies each prescriber from their
Part D type and the panel marks pre-coverage exposure as missing.

## Provenance

`make provenance` writes `results/provenance.json`: for every input file, its
source endpoint, row count, size, modification time and SHA-256, plus library
versions. Re-fetched data can be checked against it.

## The CMS catalog changes shape

Between two runs on the same day, data.cms.gov moved from one catalog record per
dataset holding every year, to one record per year titled
`<dataset> : YYYY-MM-DD`. The API endpoints did not change. `part_d_distributions`
handles both layouts and `tests/test_datasets.py` covers each.

## A note on the network client

The Open Payments portal answers HTTP 403 to some custom `User-Agent` strings.
Requests go out with the standard library default; override via
`rxinc.datasets.USER_AGENT` only if you have reason to.

## Rate and volume

Part D queries page at 5,000 rows. `fetch_part_d` caps at 50,000 rows by
default so an accidental unfiltered call does not pull 1.4 million. Open
Payments downloads stream in 1 MB chunks and skip work if the file is present.
