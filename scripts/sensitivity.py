"""Rambachan and Roth smoothness sensitivity for the saved event studies.

    python3 scripts/sensitivity.py

Reads event-study coefficients and covariances from results/estimates.json, so
no re-estimation is needed, and writes results/sensitivity.json.

For each event study it reports breakdown values for every post-payment year,
for the first two years together (the headline target, fixed in advance to
avoid choosing whichever single year looks best), and for all years averaged.
It also records fixed-length confidence intervals over a grid of M for the
headline target, and the observed bending of the pre-period trend with standard
errors, which is what M should be judged against.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rxinc.estimators import EventStudyResult  # noqa: E402
from rxinc.sensitivity import breakdown_m, flci  # noqa: E402

STUDIES = {
    "physician": "Claims, physicians",
    "physician_extensive": "Prescribes at all, physicians",
    "npp": "Claims, non-physician practitioners",
}


def load(es: dict) -> EventStudyResult:
    return EventStudyResult(
        rel_periods=np.array(es["rel_periods"]), coefs=np.array(es["coefs"]),
        ses=np.array(es["ses"]), vcov=np.array(es["vcov"]), n_obs=0, n_clusters=0,
    )


def bending(res: EventStudyResult, reference: int = -1) -> list[dict]:
    """Second differences of the leads, with standard errors from the covariance."""
    rel = res.rel_periods.astype(int)
    idx = {int(t): i for i, t in enumerate(rel) if t < 0}
    times = sorted([*idx, reference])
    out = []
    for t in times[1:-1]:
        a = np.zeros(len(rel))
        for p, c in ((t - 1, 1.0), (t, -2.0), (t + 1, 1.0)):
            if p != reference:
                a[idx[p]] += c
        out.append({"centred_on": t, "value": float(a @ res.coefs),
                    "se": float(np.sqrt(a @ res.vcov @ a))})
    return out


def main() -> None:
    data = json.loads((ROOT / "results" / "estimates.json").read_text())
    out = {"generated": time.strftime("%Y-%m-%d"),
           "restriction": "Delta^SD(M): slope of the differential trend changes by at most M per year",
           "method": "fixed-length confidence interval, Rambachan and Roth (2023)",
           "headline_target": "average of the first two post-payment years (t0, t+1)",
           "units": "outcome units per year squared; x100 gives log points or percentage points",
           "studies": {}}
    for key, label in STUDIES.items():
        res = load(data["event_studies"][key]["cohort (Sun-Abraham)"])
        rel = res.rel_periods.astype(int)
        n_post = int((rel >= 0).sum())
        first_two = np.zeros(n_post)
        first_two[: min(2, n_post)] = 1.0 / min(2, n_post)

        targets = {}
        for j, t in enumerate(rel[rel >= 0]):
            ell = np.zeros(n_post)
            ell[j] = 1.0
            c0 = flci(res, 0.0, target=ell)
            targets[f"t{t:+d}"] = {"breakdown_m": breakdown_m(res, target=ell),
                                   "estimate_m0": c0.estimate, "lower_m0": c0.lower, "upper_m0": c0.upper}
        for name, ell in (("first_two_years", first_two), ("all_years", None)):
            c0 = flci(res, 0.0, target=ell)
            targets[name] = {"breakdown_m": breakdown_m(res, target=ell),
                             "estimate_m0": c0.estimate, "lower_m0": c0.lower, "upper_m0": c0.upper}

        bends = bending(res)
        # Zero to twice the breakdown value: the region where the conclusion changes.
        top = max(targets["first_two_years"]["breakdown_m"] * 2.0, 1e-3)
        grid = []
        for m in np.linspace(0.0, top, 25):
            ci = flci(res, float(m), target=first_two)
            grid.append({"m": float(m), "estimate": ci.estimate, "lower": ci.lower,
                         "upper": ci.upper, "max_bias": ci.max_bias})

        print(f"\n[{label}]")
        print("  pre-period bending (x100): " + ", ".join(
            f"{b['centred_on']:+d}: {100 * b['value']:+.2f} (se {100 * b['se']:.2f})" for b in bends))
        for name, v in targets.items():
            print(f"  {name:16s} breakdown M x100 = {100 * v['breakdown_m']:5.2f}   "
                  f"M=0: {100 * v['estimate_m0']:+.2f} [{100 * v['lower_m0']:+.2f}, {100 * v['upper_m0']:+.2f}]")
        out["studies"][key] = {"label": label, "pre_period_bending": bends,
                               "targets": targets, "grid_first_two_years": grid}

    path = ROOT / "results" / "sensitivity.json"
    path.write_text(json.dumps(out, indent=2))
    print(f"\nwrote {path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
