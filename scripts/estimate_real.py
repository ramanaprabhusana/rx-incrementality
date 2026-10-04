"""Run the triple difference on the CMS diabetes panel.

    python3 scripts/estimate_real.py                 # full run
    python3 scripts/estimate_real.py --sample 0.05   # dry run on 5% of prescribers

Writes ``results/estimates.json``. Every specification carries a ``round``:

1. Fixed before the first full run.
2. Added after it: the 2-way and published-style comparisons.
3. Added after discovering that Open Payments only covers non-physician
   practitioners from 2021, which moved the primary population to physicians,
   plus robustness checks the documentation had promised but not delivered.
4. Added to address selection through Part D's 11-claim reporting floor.

Labelling rounds is the honest substitute for a pre-registration this project
never had.
"""

from __future__ import annotations

import argparse
import json
import math
import resource
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rxinc.diagnostics import detrend_event_study, pretrend_test  # noqa: E402
from rxinc.drugdata import build_drug_panel, load_part_d, load_payments  # noqa: E402
from rxinc.estimators import cohort_event_study, drug_event_study, triple_diff  # noqa: E402

YEARS = list(range(2019, 2025))
TWO_WAY = ("physician_year", "drug_year")
CAREY = ("physician_drug", "drug_year")
THREE_WAY = ("physician_year", "drug_year", "physician_drug")
MIN_FAMILY_PEAK_CLAIMS = 100_000
RESULTS: list[dict] = []


def mem_gb() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30


def pct(b: float) -> float:
    return 100.0 * (math.exp(b) - 1.0)


def run(label, frame, outcome, treatment, controls=(), absorb=THREE_WAY,
        rnd=1, population="physician", note=""):
    t0 = time.time()
    sub = frame.dropna(subset=[outcome, treatment, *controls])
    try:
        est = triple_diff(sub, outcome=outcome, treatment=treatment,
                          controls=controls, absorb=absorb)
    except ValueError as exc:
        print(f"  {label:50s} SKIPPED: {exc}", flush=True)
        RESULTS.append({"spec": label, "round": rnd, "population": population, "skipped": str(exc)})
        return None
    row = {
        "spec": label, "round": rnd, "population": population, "absorb": list(absorb),
        "outcome": outcome, "treatment": treatment, "controls": list(controls),
        "coef": est.coef, "se": est.se, "ci_low": est.ci95[0], "ci_high": est.ci95[1],
        "t": est.t_stat, "n_obs": est.n_obs, "n_clusters": est.n_clusters,
        "pct": pct(est.coef) if outcome == "log_rx" else None,
        "extra": est.extra, "note": note,
    }
    RESULTS.append(row)
    shown = f"{row['pct']:+5.1f}%" if row["pct"] is not None else f"{100*est.coef:+5.2f}pt"
    print(f"  {label:50s} {est.coef:+.4f} (se {est.se:.4f}) {shown}  n={est.n_obs:>9,}  "
          f"[{time.time()-t0:.0f}s]", flush=True)
    for k, v in est.extra.items():
        if not k.endswith("_se"):
            print(f"  {'':50s}   {k}: {v:+.4f} (se {est.extra[k + '_se']:.4f})", flush=True)
    return est


def event_studies(frame, label, population, rnd, outcome="log_rx", pooled=True):
    es_frame = frame[frame["exposure_known"] & ~frame["left_censored"]]
    out = {}
    fns = [("cohort (Sun-Abraham)", lambda f: cohort_event_study(
        f, outcome=outcome, cohort_col="first_treat_period", time_col="year", absorb=THREE_WAY))]
    if pooled:
        fns.append(("pooled dummies", lambda f: drug_event_study(
            f, leads=3, lags=3, outcome=outcome, absorb=THREE_WAY)))
    for name, fn in fns:
        t0 = time.time()
        es = fn(es_frame)
        pt = pretrend_test(es)
        print(f"  [{label}] {name}  ({time.time()-t0:.0f}s)", flush=True)
        for r, b, s in zip(es.rel_periods, es.coefs, es.ses, strict=True):
            shown = f"{pct(b):+5.1f}%" if outcome == "log_rx" else f"{100*b:+5.2f}pt"
            print(f"     t{int(r):+d}: {b:+.4f} (se {s:.4f})  {shown}", flush=True)
        print(f"     {pt}", flush=True)
        dt = detrend_event_study(es)
        post_mean, post_se = dt.post_average()
        print(f"     trend-adjusted: slope {dt.slope:+.4f}/yr (se {dt.slope_se:.4f}); "
              f"post average {post_mean:+.4f} (se {post_se:.4f})", flush=True)
        for r, b, s in zip(dt.rel_periods, dt.coefs, dt.ses, strict=True):
            if r >= 0:
                print(f"       adj t{int(r):+d}: {b:+.4f} (se {s:.4f})", flush=True)
        out[name] = {"rel_periods": es.rel_periods.tolist(), "coefs": es.coefs.tolist(),
                     "ses": es.ses.tolist(), "vcov": es.vcov.tolist(),
                     "detrended": {"slope": dt.slope, "slope_se": dt.slope_se,
                                   "coefs": dt.coefs.tolist(), "ses": dt.ses.tolist(),
                                   "post_average": post_mean, "post_average_se": post_se},
                     "pretrend_chi2": pt.statistic, "pretrend_df": pt.df,
                     "pretrend_p": pt.p_value, "round": rnd, "population": population,
                     "outcome": outcome}
        if es.detail is not None:
            out[name]["cohorts"] = es.detail.to_dict(orient="records")
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", type=float, default=1.0)
    ap.add_argument("--op-dir", default="data/raw")
    ap.add_argument("--out", default="results/estimates.json")
    args = ap.parse_args()

    t0 = time.time()
    claims, attrs = load_part_d([ROOT / f"data/raw/diabetes_dy{y}.csv" for y in YEARS])
    pay = load_payments({y: ROOT / args.op_dir / f"openpayments_diabetes_{y}.csv" for y in YEARS})
    peak = claims.groupby(["family", "year"])["claims"].sum().groupby("family").max()
    families = sorted(peak[peak >= MIN_FAMILY_PEAK_CLAIMS].index)
    print(f"loaded {len(claims):,} claim cells, {len(pay):,} payment cells ({time.time()-t0:.0f}s)", flush=True)

    if args.sample < 1.0:
        rng = np.random.default_rng(20261004)
        npis = claims["npi"].unique()
        keep = set(rng.choice(npis, size=int(len(npis) * args.sample), replace=False))
        claims, attrs, pay = claims[claims.npi.isin(keep)], attrs[attrs.npi.isin(keep)], pay[pay.npi.isin(keep)]
        print(f"sampled {len(keep):,} prescribers", flush=True)

    panels, reports = {}, {}
    for key, kw in (
        ("physician", dict(population="physician")),
        ("npp", dict(population="npp")),
        ("physician_unbalanced", dict(population="physician", require_all_years=False, full_grid=False)),
    ):
        panels[key], reports[key] = build_drug_panel(claims, attrs, pay, families=families, **kw)
        print(f"\nPANEL {key}: {reports[key]}", flush=True)
    print(f"families: {', '.join(reports['physician'].families)}", flush=True)
    print(f"built in {time.time()-t0:.0f}s, peak {mem_gb():.1f} GB", flush=True)

    phys = panels["physician"]
    obs = phys[phys["any_rx"] == 1]

    print("\n== 1. Primary: physicians, 3-way fixed effects ==", flush=True)
    run("intensive: log claims ~ paid", obs, "log_rx", "pay_any", rnd=3)
    run("extensive: prescribes at >= 11 claims ~ paid", phys, "any_rx", "pay_any", rnd=3)
    run("intensive: log claims ~ log(1 + $)", obs, "log_rx", "pay_log", rnd=3)

    print("\n== 1b. Selection through the 11-claim floor ==", flush=True)
    run("intensive, pairs above the floor in every year", obs[obs.pair_always_observed],
        "log_rx", "pay_any", rnd=4,
        note="treatment cannot change whether these pairs are observed")

    print("\n== 2. Which fixed effects matter (falsification: next year's payment) ==", flush=True)
    for name, ab in (("2-way: physician-year + drug-year", TWO_WAY),
                     ("published: physician-drug + drug-year", CAREY),
                     ("3-way: all three", THREE_WAY)):
        run(f"{name}", obs, "log_rx", "pay_any", absorb=ab, rnd=2)
        run(f"{name} + next year's payment", obs, "log_rx", "pay_any",
            controls=("pay_any_lead",), absorb=ab, rnd=2)

    print("\n== 3. Timing ==", flush=True)
    run("paid last year", obs, "log_rx", "pay_any_lag", rnd=1)

    print("\n== 4. Payment type and size ==", flush=True)
    run("food only vs any non-food (speaker, consulting, travel)", obs, "log_rx", "pay_food_only",
        controls=("pay_nonfood_any",), rnd=3, note="coef is food-only; non-food in extra")
    run("payment >= $25", obs, "log_rx", "pay_ge25", rnd=3)
    run("payment >= $100", obs, "log_rx", "pay_ge100", rnd=3)

    print("\n== 5. Spillovers ==", flush=True)
    run("paid + same-city peer share paid", obs, "log_rx", "pay_any",
        controls=("peer_pay_share",), rnd=1)
    run("paid + paid about another drug by same manufacturer", obs, "log_rx", "pay_any",
        controls=("pay_same_mfr_other",), rnd=3)

    print("\n== 6. Heterogeneity ==", flush=True)
    for cls in ("GLP-1", "SGLT2", "DPP-4"):
        sub = obs[obs.drug_class == cls]
        run(f"class: {cls}", sub, "log_rx", "pay_any", rnd=1)
    for grp in ("Family or general practice", "Internal medicine", "Endocrinology",
                "Cardiology", "Nephrology", "Other"):
        run(f"specialty: {grp}", obs[obs.specialty_group == grp], "log_rx", "pay_any", rnd=3)
    g = obs[obs.drug_class == "GLP-1"]
    run("GLP-1, 2019-2021", g[g.year <= 2021], "log_rx", "pay_any", rnd=1)
    run("GLP-1, 2022-2024", g[g.year >= 2022], "log_rx", "pay_any", rnd=1)
    run("GLP-1, 2022-2024, excluding Mounjaro", g[(g.year >= 2022) & (g.family != "MOUNJARO")],
        "log_rx", "pay_any", rnd=3)
    run("all classes excluding Mounjaro (launch)", obs[obs.family != "MOUNJARO"], "log_rx", "pay_any", rnd=3)

    print("\n== 7. Other populations ==", flush=True)
    npp = panels["npp"]
    npp_obs = npp[npp["any_rx"] == 1]
    run("non-physician practitioners, 2021-2024: intensive", npp_obs, "log_rx", "pay_any",
        rnd=3, population="npp")
    run("non-physician practitioners, 2021-2024: extensive", npp, "any_rx", "pay_any",
        rnd=3, population="npp")
    for grp in ("Nurse practitioner", "Physician assistant"):
        run(f"specialty: {grp}", npp_obs[npp_obs.specialty_group == grp], "log_rx", "pay_any",
            rnd=3, population="npp")
    unb = panels["physician_unbalanced"]
    run("physicians, unbalanced panel (all years present not required)", unb, "log_rx", "pay_any",
        rnd=3, population="physician_unbalanced")

    print("\n== 8. Event studies ==", flush=True)
    events = {"physician": event_studies(phys[phys.any_rx == 1], "physicians", "physician", 3),
              "physician_extensive": event_studies(phys, "physicians, extensive margin", "physician", 4,
                                                   outcome="any_rx", pooled=False),
              "npp": event_studies(npp_obs, "non-physician practitioners", "npp", 3)}

    out = {
        "generated": time.strftime("%Y-%m-%d"),
        "sample": args.sample,
        "rounds": {"1": "fixed before first full run",
                   "2": "added after first full run: fixed-effect comparisons",
                   "3": "added after discovering 2021 coverage of non-physician practitioners",
                   "4": "added to address selection through the 11-claim reporting floor"},
        "panels": {k: {kk: (vv.item() if isinstance(vv, (np.integer, np.floating)) else vv)
                       for kk, vv in r.__dict__.items()} for k, r in reports.items()},
        "estimates": RESULTS,
        "event_studies": events,
    }
    path = ROOT / args.out
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=float))
    print(f"\nwrote {path.relative_to(ROOT)}  total {time.time()-t0:.0f}s, peak {mem_gb():.1f} GB", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
