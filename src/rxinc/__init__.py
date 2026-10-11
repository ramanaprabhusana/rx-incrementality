"""Causal estimators for industry payments and physician prescribing.

The package exists to answer a narrow question honestly: how much does a given
estimator overstate the effect of pharmaceutical industry payments on
prescribing, when industry targets physicians who were already prescribing
heavily or whose prescribing was already rising?

Because that question cannot be settled on observational data alone, every
estimator is first run against simulated panels whose true effect is known by
construction (:mod:`rxinc.simulate`, :mod:`rxinc.drugpanel`), and its bias and
interval coverage reported. The real-data pipeline (:mod:`rxinc.datasets`,
:mod:`rxinc.crosswalk`, :mod:`rxinc.drugdata`) then builds a
physician-by-drug-by-year panel from CMS Open Payments and Medicare Part D, and
the same estimator code is pointed at it. :mod:`rxinc.diagnostics` and
:mod:`rxinc.sensitivity` test the parallel-trends assumption and bound what
happens when it fails.
"""

from rxinc.drugpanel import DrugPanelConfig, simulate_drug_panel
from rxinc.estimators import (
    Estimate,
    EventStudyResult,
    did_never_treated,
    drug_event_study,
    event_study,
    interrupted_time_series,
    naive_ols,
    triple_diff,
    twoway_fe,
)
from rxinc.simulate import PanelConfig, simulate_panel

__version__ = "0.1.0"

__all__ = [
    "DrugPanelConfig",
    "Estimate",
    "EventStudyResult",
    "PanelConfig",
    "did_never_treated",
    "drug_event_study",
    "event_study",
    "interrupted_time_series",
    "naive_ols",
    "simulate_drug_panel",
    "simulate_panel",
    "triple_diff",
    "twoway_fe",
]
