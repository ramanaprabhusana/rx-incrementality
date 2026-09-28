"""Reproduce the headline result: which estimators survive which targeting.

Run with ``make demo`` or ``python3 examples/bias_demo.py``.

The script simulates three physician panels that differ only in how industry
chooses whom to pay, applies every estimator in the package to each, and
reports the bias against the known true effect.  It then runs the diagnostic
battery, to show that the failure under dynamic targeting is detectable from
the data alone -- which is the practically useful part, since on real data the
true effect is exactly what you do not have.
"""

from __future__ import annotations

from rxinc.diagnostics import monte_carlo, placebo_shift, pretrend_test
from rxinc.estimators import (
    did_never_treated,
    event_study,
    interrupted_time_series,
    naive_ols,
    twoway_fe,
)
from rxinc.simulate import PanelConfig, simulate_panel

REGIMES = {
    "none": "Payments assigned at random",
    "static": "Industry targets persistently high prescribers",
    "dynamic": "Industry targets physicians whose Rx is already rising",
}

RULE = "=" * 76


def _its_level(frame):
    """Adapter so interrupted_time_series fits the single-Estimate interface."""
    return interrupted_time_series(frame)[0]


ESTIMATORS = {
    "Naive pooled OLS": naive_ols,
    "Two-way fixed effects": twoway_fe,
    "DiD vs never-treated": lambda d: did_never_treated(d, n_boot=60),
    "ITS level change": _its_level,
}


def main() -> None:
    print(RULE)
    print("SINGLE PANEL: estimates vs known truth")
    print(RULE)

    for regime, description in REGIMES.items():
        config = PanelConfig(targeting=regime)
        panel = simulate_panel(config)
        print(f"\n[{regime.upper()}] {description}")
        print(f"  true effect = {panel.attrs['true_effect']:+.4f} log points\n")

        for estimate in (
            naive_ols(panel),
            twoway_fe(panel),
            did_never_treated(panel, n_boot=150),
            interrupted_time_series(panel)[0],
        ):
            print(f"    {estimate}")

        study = event_study(panel)
        print(f"\n    {pretrend_test(study)}")
        print(f"    {placebo_shift(panel)}")

    print(f"\n{RULE}")
    print("MONTE CARLO: mean bias over 40 replications")
    print(RULE)

    for regime in REGIMES:
        config = PanelConfig(targeting=regime, n_physicians=1200)
        table = monte_carlo(config, ESTIMATORS, n_reps=40)
        table = table.loc[list(ESTIMATORS)]
        print(f"\n[{regime.upper()}]  true effect = {config.tau:+.4f}")
        print(table.round(4).to_string())


if __name__ == "__main__":
    main()
