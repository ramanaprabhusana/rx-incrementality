"""Download a therapeutic class panel from Part D by Provider and Drug.

Writes one CSV per data year to ``data/raw/``. Resumable: a year whose file
already exists and is non-empty is skipped, so an interrupted run can simply be
restarted.

    python3 scripts/fetch_class_panel.py --class diabetes --start 2016 --end 2024
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from rxinc.datasets import (  # noqa: E402
    PART_D_BY_PROVIDER_DRUG_TITLE,
    fetch_part_d,
    part_d_distributions,
)

CLASSES: dict[str, dict[str, list[str]]] = {
    "diabetes": {
        "SGLT2": [
            "Jardiance", "Farxiga", "Invokana", "Steglatro",
            "Synjardy", "Xigduo XR", "Glyxambi", "Invokamet",
        ],
        "DPP-4": [
            "Januvia", "Janumet", "Janumet XR", "Tradjenta",
            "Onglyza", "Jentadueto", "Kombiglyze XR",
        ],
        "GLP-1": [
            "Ozempic", "Trulicity", "Rybelsus", "Victoza",
            "Mounjaro", "Byetta", "Bydureon Bcise", "Adlyxin",
        ],
    },
    "anticoagulant": {
        "DOAC": ["Eliquis", "Xarelto", "Pradaxa", "Savaysa"],
        "VKA": ["Warfarin Sodium", "Jantoven"],
    },
}


def year_of(distribution_title: str) -> int:
    """Extract the data year from a distribution title.

    CMS titles these with the data year end date, e.g. ``... : 2023-12-31``
    means data year 2023.
    """
    return int(distribution_title.split(":")[-1].strip()[:4])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--class", dest="klass", default="diabetes", choices=sorted(CLASSES))
    parser.add_argument("--start", type=int, default=2016)
    parser.add_argument("--end", type=int, default=2024)
    parser.add_argument("--out-dir", default="data/raw")
    args = parser.parse_args()

    classes = CLASSES[args.klass]
    brand_to_class = {b: c for c, bs in classes.items() for b in bs}
    brands = sorted(brand_to_class)

    dists = [d for d in part_d_distributions(PART_D_BY_PROVIDER_DRUG_TITLE) if d.api_url]
    by_year: dict[int, str] = {}
    for dist in dists:
        y = year_of(dist.title)
        by_year.setdefault(y, dist.api_url)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    years = [y for y in sorted(by_year) if args.start <= y <= args.end]
    print(f"class={args.klass} brands={len(brands)} years={years}", flush=True)

    grand_total = 0
    for year in years:
        target = out_dir / f"{args.klass}_dy{year}.csv"
        if target.exists() and target.stat().st_size > 0:
            existing = sum(1 for _ in open(target)) - 1
            print(f"[{year}] skip, exists ({existing:,} rows)", flush=True)
            grand_total += existing
            continue

        start = time.time()
        frames = []
        for brand in brands:
            try:
                frame = fetch_part_d(
                    by_year[year], {"Brnd_Name": brand}, max_rows=None
                )
            except Exception as exc:  # noqa: BLE001
                print(f"[{year}] {brand}: FAILED {type(exc).__name__}: {exc}", flush=True)
                continue
            if len(frame):
                frame["therapeutic_class"] = brand_to_class[brand]
                frame["data_year"] = year
                frames.append(frame)

        if not frames:
            print(f"[{year}] no rows for any brand", flush=True)
            continue

        panel = pd.concat(frames, ignore_index=True)
        panel.to_csv(target, index=False)
        grand_total += len(panel)
        elapsed = time.time() - start
        print(
            f"[{year}] {len(panel):>8,} rows  "
            f"{panel['Prscrbr_NPI'].nunique():>7,} prescribers  "
            f"{panel['Brnd_Name'].nunique():>3} brands  "
            f"{elapsed/60:5.1f} min  -> {target.name}",
            flush=True,
        )

    print(f"\nDONE. {grand_total:,} rows across {len(years)} years.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
