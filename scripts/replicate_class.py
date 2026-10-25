"""Replicate the core design on another drug class.

    python3 scripts/replicate_class.py --class anticoagulant

Runs the specifications that carry the diabetes conclusions (the three-way
fixed-effects estimate, the next-year falsification test across fixed-effect
sets, the always-observed pairs, payment type, the cohort event study and the
dollar translation) on another class, with the same code and thresholds, and
writes ``results/replication_<class>.json``. The diabetes results in
``results/estimates.json`` are untouched.

Every specification carries round ``R``: chosen before the replication was
run, as the subset of the diabetes analysis that its conclusions rest on.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rxinc.drugdata import build_drug_panel, load_part_d, load_payments  # noqa: E402

_spec = importlib.util.spec_from_file_location("estimate_real", ROOT / "scripts" / "estimate_real.py")
er = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(er)

YEARS = er.YEARS


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--class", dest="klass", default="anticoagulant")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    klass = args.klass
    out_path = ROOT / (args.out or f"results/replication_{klass}.json")

    t0 = time.time()
    claims, attrs = load_part_d([ROOT / f"data/raw/{klass}_dy{y}.csv" for y in YEARS], klass=klass)
    pay = load_payments({y: ROOT / f"data/raw/openpayments_{klass}_{y}.csv" for y in YEARS}, klass=klass)
    quality = {"part_d": dict(claims.attrs["quality"]), "open_payments": dict(pay.attrs["quality"])}
    peak = claims.groupby(["family", "year"])["claims"].sum().groupby("family").max()
    families = sorted(peak[peak >= er.MIN_FAMILY_PEAK_CLAIMS].index)
    print(f"{klass}: families {families}; loaded in {time.time()-t0:.0f}s", flush=True)

    phys, report = build_drug_panel(claims, attrs, pay, families=families, population="physician")
    print(f"PANEL: {report}", flush=True)
    obs = phys[phys["any_rx"] == 1]
    er.RESULTS.clear()
    R = "R"

    print("\n== Primary ==", flush=True)
    er.run("intensive: log claims ~ paid", obs, "log_rx", "pay_any", rnd=R)
    er.run("extensive: prescribes at >= 11 claims ~ paid", phys, "any_rx", "pay_any", rnd=R)
    er.run("intensive, pairs above the floor in every year", obs[obs.pair_always_observed], "log_rx", "pay_any", rnd=R)

    print("\n== Which fixed effects matter (next year's payment) ==", flush=True)
    for name, ab in (("2-way: physician-year + drug-year", er.TWO_WAY),
                     ("published: physician-drug + drug-year", er.CAREY),
                     ("3-way: all three", er.THREE_WAY)):
        er.run(name, obs, "log_rx", "pay_any", absorb=ab, rnd=R)
        er.run(f"{name} + next year's payment", obs, "log_rx", "pay_any", controls=("pay_any_lead",), absorb=ab, rnd=R)

    print("\n== Robustness ==", flush=True)
    er.run("food only vs any non-food (speaker, consulting, travel)", obs, "log_rx", "pay_food_only",
           controls=("pay_nonfood_any",), rnd=R)
    if "WARFARIN" in families:
        er.run("promoted drugs only (excluding warfarin)", obs[obs.family != "WARFARIN"], "log_rx", "pay_any", rnd=R)

    print("\n== Event study ==", flush=True)
    events = {"physician": er.event_studies(obs, "physicians", "physician", R, pooled=False)}

    print("\n== Dollars ==", flush=True)
    spend = er.run("intensive: log gross spending ~ paid", obs, "log_cost", "pay_any", rnd=R)
    spend_lead = er.run("intensive: log gross spending ~ paid + next year's payment", obs, "log_cost", "pay_any",
                        controls=("pay_any_lead",), rnd=R)
    value = {"all_years": er.payment_value(phys, {"headline": spend.coef, "conservative": spend_lead.coef})}
    v = value["all_years"]
    print(f"  payments ${v['payment_dollars'] / 1e6:,.1f}M; naive ${v['naive_per_payment_dollar']:.0f}; caused "
          + ", ".join(f"{k} ${c['per_payment_dollar']:.2f}" for k, c in v["causal"].items()), flush=True)

    out = {
        "generated": time.strftime("%Y-%m-%d"),
        "class": klass,
        "families": families,
        "panel": {k: (vv.item() if isinstance(vv, (np.integer, np.floating)) else vv)
                  for k, vv in report.__dict__.items()},
        "data_quality": quality,
        "estimates": er.RESULTS,
        "event_studies": events,
        "payment_value": value,
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(out, indent=2, default=float))
    print(f"\nwrote {out_path.relative_to(ROOT)} in {time.time()-t0:.0f}s, peak {er.mem_gb():.1f} GB", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
