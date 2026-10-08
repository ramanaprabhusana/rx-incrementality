"""Record exactly which CMS data produced results/estimates.json.

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
    files = []
    for y in YEARS:
        files.append({"dataset": "Medicare Part D Prescribers by Provider and Drug", "year": y,
                      **describe(ROOT / f"data/raw/diabetes_dy{y}.csv", "Prscrbr_NPI", part_d_api.get(y))})
        try:
            op_url = open_payments_general_distribution(y).download_url
        except Exception:  # noqa: BLE001
            op_url = None
        files.append({"dataset": "Open Payments general payments (diabetes brand families)", "year": y,
                      **describe(ROOT / f"data/raw/openpayments_diabetes_{y}.csv",
                                 "Covered_Recipient_NPI", op_url)})
        print(f"  {y}: hashed")
    est = ROOT / "results" / "estimates.json"
    sens = ROOT / "results" / "sensitivity.json"
    out = {
        "generated": time.strftime("%Y-%m-%d"),
        "estimates_sha256": sha256(est) if est.exists() else None,
        "sensitivity_sha256": sha256(sens) if sens.exists() else None,
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
