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
    part_d_api_by_year,
)

CLASSES: dict[str, dict[str, list[str]]] = {
    # Exact Part D Brnd_Name strings, discovered by generic molecule from the
    # national Geography and Drug file rather than guessed. Part D splits brand
    # families across pack sizes, devices and formulations (Victoza 2-Pak and
    # 3-Pak, Bydureon Pen and Bcise, XR variants), and guessing missed 1.58M
    # Victoza claims in 2019 alone. Generic-name search alone also misses
    # products, because CMS abbreviates triple-combination generics
    # ("Empaglifloz/Linaglip/Metformin" for Trijardy Xr). Excluded by design: insulin combinations
    # (Soliqua, Xultophy), obesity indications (Wegovy, Saxenda, Zepbound) and
    # unbranded generics, which are never promoted.
    "diabetes": {
        "SGLT2": [
            "Jardiance", "Farxiga", "Invokana", "Steglatro",
            "Synjardy", "Synjardy Xr", "Xigduo Xr", "Glyxambi",
            "Invokamet", "Invokamet Xr", "Qtern", "Segluromet", "Steglujan",
            "Trijardy Xr",
        ],
        "DPP-4": [
            "Januvia", "Janumet", "Janumet Xr", "Tradjenta", "Onglyza",
            "Jentadueto", "Jentadueto Xr", "Kombiglyze Xr",
            "Nesina", "Kazano", "Oseni",
        ],
        "GLP-1": [
            "Ozempic", "Trulicity", "Rybelsus", "Victoza 2-Pak", "Victoza 3-Pak",
            "Mounjaro", "Byetta", "Bydureon Bcise", "Bydureon Pen", "Bydureon",
            "Adlyxin",
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

    by_year = part_d_api_by_year(PART_D_BY_PROVIDER_DRUG_TITLE)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    years = [y for y in sorted(by_year) if args.start <= y <= args.end]
    print(f"class={args.klass} brands={len(brands)} years={years}", flush=True)

    grand_total = 0
    for year in years:
        target = out_dir / f"{args.klass}_dy{year}.csv"
        existing = None
        have: set[str] = set()
        if target.exists() and target.stat().st_size > 0:
            existing = pd.read_csv(target, dtype={"Prscrbr_NPI": str}, low_memory=False)
            have = set(existing["Brnd_Name"].unique())
        # Re-running after the brand list grows fills only the gaps. A brand
        # with zero rows in a year is re-queried, which is cheap.
        todo = [b for b in brands if b not in have]
        if not todo:
            print(f"[{year}] complete, {len(existing):,} rows", flush=True)
            grand_total += len(existing)
            continue

        start = time.time()
        frames = [] if existing is None else [existing]
        added = 0
        for brand in todo:
            try:
                frame = fetch_part_d(by_year[year], {"Brnd_Name": brand}, max_rows=None)
            except Exception as exc:  # noqa: BLE001
                print(f"[{year}] {brand}: FAILED {type(exc).__name__}: {exc}", flush=True)
                continue
            if len(frame):
                frame["therapeutic_class"] = brand_to_class[brand]
                frame["data_year"] = year
                frames.append(frame)
                added += len(frame)
                print(f"[{year}]   + {brand:16s} {len(frame):>7,} rows", flush=True)

        panel = pd.concat(frames, ignore_index=True)
        panel.to_csv(target, index=False)
        grand_total += len(panel)
        print(
            f"[{year}] {len(panel):>8,} rows (+{added:,})  "
            f"{panel['Prscrbr_NPI'].nunique():>7,} prescribers  "
            f"{panel['Brnd_Name'].nunique():>3} brands  "
            f"{(time.time()-start)/60:5.1f} min",
            flush=True,
        )

    print(f"\nDONE. {grand_total:,} rows across {len(years)} years.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
