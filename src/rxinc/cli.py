"""Command line entry point.

    rxinc demo                       reproduce the headline bias result
    rxinc simulate --targeting ...   write a simulated panel to CSV
    rxinc estimate --panel p.csv     run every estimator and diagnostic
    rxinc catalog                    list available CMS data years
    rxinc fetch-partd --state IN     pull Part D prescriber rows to CSV
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from rxinc.diagnostics import balance_table, placebo_shift, pretrend_test
from rxinc.estimators import (
    did_never_treated,
    event_study,
    interrupted_time_series,
    naive_ols,
    twoway_fe,
)
from rxinc.simulate import VALID_TARGETING, PanelConfig, simulate_panel


def _cmd_demo(_: argparse.Namespace) -> int:
    from rxinc.demo import main as run_demo  # noqa: PLC0415

    run_demo()
    return 0


def _cmd_simulate(args: argparse.Namespace) -> int:
    config = PanelConfig(
        n_physicians=args.n_physicians,
        n_periods=args.n_periods,
        tau=args.tau,
        targeting=args.targeting,
        seed=args.seed,
    )
    panel = simulate_panel(config)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    panel.to_csv(out, index=False)
    print(
        f"Wrote {len(panel):,} rows "
        f"({panel['physician_id'].nunique():,} physicians x {config.n_periods} periods) "
        f"to {out}"
    )
    print(f"True effect: {config.tau:+.4f}  targeting: {config.targeting}")
    return 0


def _cmd_estimate(args: argparse.Namespace) -> int:
    panel = pd.read_csv(args.panel)
    if args.true_effect is not None:
        panel.attrs["true_effect"] = args.true_effect

    print(f"Panel: {len(panel):,} rows, "
          f"{panel['physician_id'].nunique():,} physicians\n")

    print("Estimates")
    print("-" * 74)
    for estimate in (
        naive_ols(panel),
        twoway_fe(panel),
        did_never_treated(panel, n_boot=args.n_boot),
        interrupted_time_series(panel)[0],
    ):
        print(f"  {estimate}")

    print("\nDiagnostics")
    print("-" * 74)
    study = event_study(panel, leads=args.leads, lags=args.lags)
    trend = pretrend_test(study)
    print(f"  {trend}")
    print(f"  {placebo_shift(panel)}")
    print()
    print(balance_table(panel).to_string())

    if not trend.passes:
        print(
            "\n  WARNING: pre-trends reject. Treated physicians were already "
            "diverging\n  before any payment. Treat every estimate above as "
            "an upper bound on\n  association, not a causal effect."
        )
    return 0


def _cmd_catalog(_: argparse.Namespace) -> int:
    from rxinc.datasets import open_payments_catalog, part_d_distributions

    dists = part_d_distributions()
    print(f"Medicare Part D Prescribers - by Provider: {len(dists)} distributions")
    for dist in dists[:6]:
        print(f"  {dist.title}")
    if len(dists) > 6:
        print(f"  ... and {len(dists) - 6} more")

    years = sorted(
        record["title"].split()[0]
        for record in open_payments_catalog()
        if record.get("title", "").endswith("General Payment Data")
    )
    print(f"\nOpen Payments general payment years: {', '.join(years)}")
    return 0


def _cmd_fetch_partd(args: argparse.Namespace) -> int:
    from rxinc.datasets import fetch_part_d, part_d_distributions, part_d_row_count

    dists = [d for d in part_d_distributions() if d.api_url]
    if not dists:
        print("No queryable Part D distribution found.", file=sys.stderr)
        return 1
    dist = dists[args.year_index]
    filters = {"Prscrbr_State_Abrvtn": args.state} if args.state else {}

    total = part_d_row_count(dist.api_url, filters)
    print(f"Using: {dist.title}")
    print(f"Matching rows: {total:,} (fetching up to {args.max_rows:,})")

    frame = fetch_part_d(dist.api_url, filters, max_rows=args.max_rows)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(out, index=False)
    print(f"Wrote {len(frame):,} rows x {frame.shape[1]} columns to {out}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    """Construct the argument parser."""
    parser = argparse.ArgumentParser(
        prog="rxinc",
        description=(
            "Causal estimators for industry payments and physician prescribing."
        ),
    )
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("demo", help="reproduce the headline bias result").set_defaults(
        func=_cmd_demo
    )

    sim = sub.add_parser("simulate", help="write a simulated panel to CSV")
    sim.add_argument("--targeting", choices=VALID_TARGETING, default="dynamic")
    sim.add_argument("--n-physicians", type=int, default=2000)
    sim.add_argument("--n-periods", type=int, default=16)
    sim.add_argument("--tau", type=float, default=0.05)
    sim.add_argument("--seed", type=int, default=20260928)
    sim.add_argument("--out", default="data/processed/panel.csv")
    sim.set_defaults(func=_cmd_simulate)

    est = sub.add_parser("estimate", help="run estimators and diagnostics")
    est.add_argument("--panel", required=True, help="panel CSV")
    est.add_argument(
        "--true-effect",
        type=float,
        default=None,
        help="known effect, if the panel is simulated",
    )
    est.add_argument("--leads", type=int, default=4)
    est.add_argument("--lags", type=int, default=6)
    est.add_argument("--n-boot", type=int, default=200)
    est.set_defaults(func=_cmd_estimate)

    sub.add_parser("catalog", help="list available CMS data years").set_defaults(
        func=_cmd_catalog
    )

    fetch = sub.add_parser("fetch-partd", help="download Part D prescriber rows")
    fetch.add_argument("--state", default=None, help="two-letter state, e.g. IN")
    fetch.add_argument("--max-rows", type=int, default=50_000)
    fetch.add_argument(
        "--year-index", type=int, default=0, help="0 is the newest distribution"
    )
    fetch.add_argument("--out", default="data/raw/part_d.csv")
    fetch.set_defaults(func=_cmd_fetch_partd)

    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the CLI.

    Args:
        argv: Argument list, defaulting to ``sys.argv[1:]``.

    Returns:
        Process exit code.
    """
    args = build_parser().parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
