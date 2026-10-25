"""Record exactly which CMS data produced the results files.

The raw extracts are not committed (about 1.5 GB, and CMS republishes). This
writes results/provenance.json with, for every input file, the endpoint it came
from, its row count, size, modification time and SHA-256, plus library versions.
Someone re-fetching can compare hashes to know whether they have the same data.
"""

from __future__ import annotations

import hashlib
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rxinc.datasets import (  # noqa: E402
    PART_D_BY_PROVIDER_DRUG_TITLE,
    open_payments_general_distribution,
    part_d_api_by_year,
)

YEARS = range(2019, 2025)
CLASSES = ("diabetes", "anticoagulant")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def describe(path: Path, first_col: str, source: str | None) -> dict:
    rows = len(pd.read_csv(path, usecols=[first_col], dtype=str))
    st = path.stat()
    return {"file": str(path.relative_to(ROOT)), "rows": rows, "bytes": st.st_size,
            "modified": time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(st.st_mtime)),
            "sha256": sha256(path), "source": source}


def main() -> None:
    try:
        part_d_api = part_d_api_by_year(PART_D_BY_PROVIDER_DRUG_TITLE)
    except Exception as exc:  # noqa: BLE001
        print(f"warning: could not reach CMS catalog ({exc}); sources left blank")
        part_d_api = {}
    op_urls = {}
    for y in YEARS:
        try:
            op_urls[y] = open_payments_general_distribution(y).download_url
        except Exception:  # noqa: BLE001
            op_urls[y] = None
    files = []
    for klass in CLASSES:
        part_d = [ROOT / f"data/raw/{klass}_dy{y}.csv" for y in YEARS]
        payments = [ROOT / f"data/raw/openpayments_{klass}_{y}.csv" for y in YEARS]
        if not all(p.exists() for p in part_d + payments):
            print(f"  {klass}: inputs not all present, skipped")
            continue
        for y, pd_path, op_path in zip(YEARS, part_d, payments, strict=True):
            files.append({"dataset": "Medicare Part D Prescribers by Provider and Drug", "class": klass, "year": y,
                          **describe(pd_path, "Prscrbr_NPI", part_d_api.get(y))})
            files.append({"dataset": f"Open Payments general payments ({klass} brand families)", "class": klass,
                          "year": y, **describe(op_path, "Covered_Recipient_NPI", op_urls[y])})
            print(f"  {klass} {y}: hashed")
    est = ROOT / "results" / "estimates.json"
    sens = ROOT / "results" / "sensitivity.json"
    replications = {p.stem.removeprefix("replication_"): sha256(p)
                    for p in sorted((ROOT / "results").glob("replication_*.json"))}
    out = {
        "generated": time.strftime("%Y-%m-%d"),
        "estimates_sha256": sha256(est) if est.exists() else None,
        "sensitivity_sha256": sha256(sens) if sens.exists() else None,
        "replication_sha256": replications,
        "environment": {"python": platform.python_version(), "numpy": np.__version__,
                        "pandas": pd.__version__, "scipy": scipy.__version__,
                        "platform": platform.platform()},
        "files": files,
        "notes": [
            "Part D extracts were fetched through the data.cms.gov JSON API, filtered by Brnd_Name.",
            "Open Payments extracts were streamed from the bulk CSV and filtered by brand family in flight.",
            "Row counts exclude the header.",
        ],
    }
    path = ROOT / "results" / "provenance.json"
    path.write_text(json.dumps(out, indent=2))
    print(f"wrote {path.relative_to(ROOT)}: {len(files)} files, "
          f"{sum(f['rows'] for f in files):,} rows")


if __name__ == "__main__":
    main()
