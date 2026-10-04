"""Run the triple difference on the CMS diabetes panel.

    python3 scripts/estimate_real.py                 # full run
    python3 scripts/estimate_real.py --sample 0.1    # 10% of physicians, for a dry run

Writes ``results/estimates.json`` and prints every table. The specification set
is fixed here in code, before results are seen, so that what gets reported is
not chosen after the fact.
"""

from __future__ import annotations

import argparse
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rxinc.crosswalk import crosswalk_report  # noqa: E402
from rxinc.diagnostics import pretrend_test  # noqa: E402
from rxinc.drugdata import build_drug_panel, load_part_d, load_payments  # noqa: E402
from rxinc.estimators import drug_event_study, triple_diff  # noqa: E402

YEARS = list(range(2019, 2025))
TWO_WAY = ("physician_year", "drug_year")
CAREY = ("physician_drug", "drug_year")
THREE_WAY = ("physician_year", "drug_year", "physician_drug")
MIN_FAMILY_PEAK_CLAIMS = 100_000


def mem_gb() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30  # bytes on macOS


def run(label: str, frame: pd.DataFrame, outcome: str, treatment: str,
        controls: tuple[str, ...] = (), note: str = "",
        absorb: tuple[str, ...] = THREE_WAY) -> dict:
    t0 = time.time()
    sub = frame.dropna(subset=[outcome, treatment, *controls])
    try:
        est = triple_diff(sub, outcome=outcome, treatment=treatment,
                          controls=controls, absorb=absorb)
    except ValueError as exc:
        print(f"  {label:44s} SKIPPED: {exc}", flush=True)
        return {"spec": label, "skipped": str(exc), "absorb": list(absorb)}
    row = {
        "spec": label, "absorb": list(absorb), "outcome": outcome, "treatment": treatment,
        "controls": list(controls), "coef": est.coef, "se": est.se,
        "ci_low": est.ci95[0], "ci_high": est.ci95[1], "t": est.t_stat,
        "n_obs": est.n_obs, "n_physicians": est.n_clusters,
        "extra": est.extra, "note": note, "seconds": round(time.time() - t0, 1),
    }
    print(f"  {label:44s} {est.coef:+.4f} (se {est.se:.4f})  t={est.t_stat:6.1f}  "
          f"n={est.n_obs:>10,}  [{row['seconds']}s]", flush=True)
    for k, v in est.extra.items():
        if not k.endswith("_se"):
            print(f"  {'':44s}   {k}: {v:+.4f} (se {est.extra[k + '_se']:.4f})", flush=True)
    return row


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", type=float, default=1.0)
    ap.add_argument("--op-dir", default="data/raw")
    ap.add_argument("--out", default="results/estimates.json")
    args = ap.parse_args()

    t0 = time.time()
    claims, attrs = load_part_d([ROOT / f"data/raw/diabetes_dy{y}.csv" for y in YEARS])
    pay = load_payments({y: ROOT / args.op_dir / f"openpayments_diabetes_{y}.csv" for y in YEARS})
    print(f"loaded: {len(claims):,} claim cells, {len(pay):,} payment cells "
          f"({time.time()-t0:.0f}s, peak {mem_gb():.1f} GB)", flush=True)

    peak = claims.groupby(["family", "year"])["claims"].sum().groupby("family").max()
    families = sorted(peak[peak >= MIN_FAMILY_PEAK_CLAIMS].index)

    if args.sample < 1.0:
        rng = np.random.default_rng(20261004)
        npis = claims["npi"].unique()
        keep = set(rng.choice(npis, size=int(len(npis) * args.sample), replace=False))
        claims, attrs = claims[claims.npi.isin(keep)], attrs[attrs.npi.isin(keep)]
        pay = pay[pay.npi.isin(keep)]
        print(f"sampled {len(keep):,} physicians ({args.sample:.0%})", flush=True)

    panel, report = build_drug_panel(claims, attrs, pay, families=families)
    print(f"\nPANEL: {report}", flush=True)
    print(f"families kept: {', '.join(report.families)}", flush=True)
    print(f"families dropped (low volume): {', '.join(report.dropped_families)}", flush=True)
    print(f"built in {time.time()-t0:.0f}s, peak {mem_gb():.1f} GB\n", flush=True)

    results: list[dict] = []
    obs = panel[panel["any_rx"] == 1]

    print("== Primary (3-way: physician-year + drug-year + physician-drug) ==", flush=True)
    results.append(run("intensive: log claims ~ paid (same year)", obs, "log_rx", "pay_any"))
    results.append(run("extensive: prescribes at >=11 ~ paid", panel, "any_rx", "pay_any"))
    results.append(run("intensive: log claims ~ log(1+$)", obs, "log_rx", "pay_log"))

    print("\n== Same, 2-way (no physician-drug effects), for comparison ==", flush=True)
    results.append(run("intensive ~ paid [2-way]", obs, "log_rx", "pay_any", absorb=TWO_WAY))
    results.append(run("extensive ~ paid [2-way]", panel, "any_rx", "pay_any", absorb=TWO_WAY))
    results.append(run("intensive ~ paid + paid next year [2-way]", obs, "log_rx", "pay_any",
                       controls=("pay_any_lead",), absorb=TWO_WAY))

    print("\n== Published design (Carey, Lieber and Miller): physician-drug + drug-year ==", flush=True)
    print("   [added after the first full run, to answer what physician-year effects add]", flush=True)
    posthoc = "added after first full run; comparison to the published specification"
    results.append(run("intensive ~ paid [Carey-style]", obs, "log_rx", "pay_any",
                       absorb=CAREY, note=posthoc))
    results.append(run("extensive ~ paid [Carey-style]", panel, "any_rx", "pay_any",
                       absorb=CAREY, note=posthoc))
    results.append(run("intensive ~ paid + paid next year [Carey-style]", obs, "log_rx", "pay_any",
                       controls=("pay_any_lead",), absorb=CAREY, note=posthoc))

    print("\n== Timing and anticipation ==", flush=True)
    results.append(run("intensive ~ paid last year", obs, "log_rx", "pay_any_lag",
                       note="2020-2024 only; 2018 payments unobserved"))
    results.append(run("intensive ~ paid this year + paid next year", obs, "log_rx", "pay_any",
                       controls=("pay_any_lead",),
                       note="a large lead coefficient means targeting anticipates prescribing"))

    print("\n== Peer spillover (same city, leave one out) ==", flush=True)
    results.append(run("intensive ~ paid + peer share paid", obs, "log_rx", "pay_any",
                       controls=("peer_pay_share",)))

    print("\n== By drug class ==", flush=True)
    for cls in ("GLP-1", "SGLT2", "DPP-4"):
        sub = obs[obs.drug_class == cls]
        if sub["family"].nunique() < 2:
            print(f"  intensive, {cls} only{'':27s} SKIPPED: fewer than 2 families", flush=True)
            results.append({"spec": f"intensive, {cls} only", "skipped": "fewer than 2 families"})
            continue
        results.append(run(f"intensive, {cls} only", sub, "log_rx", "pay_any"))

    print("\n== GLP-1 era split ==", flush=True)
    g = obs[obs.drug_class == "GLP-1"]
    results.append(run("intensive, GLP-1, 2019-2021", g[g.year <= 2021], "log_rx", "pay_any"))
    results.append(run("intensive, GLP-1, 2022-2024 (Mounjaro, shortages)", g[g.year >= 2022], "log_rx", "pay_any"))

    print("\n== Event study, clean onsets only, 3-way (intensive margin) ==", flush=True)
    es_frame = obs[~obs["left_censored"]].copy()
    es = drug_event_study(es_frame, leads=3, lags=3, outcome="log_rx", absorb=THREE_WAY)
    pt = pretrend_test(es)
    for r, b, s in zip(es.rel_periods, es.coefs, es.ses):
        print(f"  t{int(r):+d}: {b:+.4f} (se {s:.4f})", flush=True)
    print(f"  {pt}", flush=True)

    out = {
        "generated": time.strftime("%Y-%m-%d"),
        "sample": args.sample,
        "panel": {k: (v if not isinstance(v, (np.integer, np.floating)) else v.item())
                  for k, v in report.__dict__.items()},
        "estimates": results,
        "event_study": {"rel_periods": es.rel_periods.tolist(), "coefs": es.coefs.tolist(),
                        "ses": es.ses.tolist(), "pretrend_chi2": pt.statistic,
                        "pretrend_df": pt.df, "pretrend_p": pt.p_value},
    }
    path = ROOT / args.out
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=float))
    print(f"\nwrote {path.relative_to(ROOT)}  total {time.time()-t0:.0f}s, peak {mem_gb():.1f} GB", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
