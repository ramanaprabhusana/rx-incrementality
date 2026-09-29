# Findings

All numbers below come from `make demo`. Panels are 1,200 physicians over 16
quarters, true effect $\tau = +0.05$ log points (roughly a 5% lift in
prescribing volume), averaged over 40 replications.

## Random assignment: the sanity check

| Estimator | Mean estimate | Mean bias | RMSE | 95% coverage |
| --- | ---: | ---: | ---: | ---: |
| Naive pooled OLS | 0.0921 | +0.0421 | 0.0562 | 78% |
| Two-way fixed effects | 0.0484 | −0.0016 | 0.0125 | 93% |
| DiD vs never-treated | 0.0481 | −0.0019 | 0.0170 | 88% |
| ITS level change | 0.0476 | −0.0024 | 0.0141 | 93% |

Everything except naive OLS recovers the truth. Naive OLS is biased *even
under random assignment*, because payment relationships accumulate over time
while the market grows, with no period controls it attributes secular growth
to treatment. Worth noting on its own: the uncontrolled comparison is wrong
before confounding enters the picture.

## Static targeting: fixed effects do their job

| Estimator | Mean estimate | Mean bias | RMSE | 95% coverage |
| --- | ---: | ---: | ---: | ---: |
| Naive pooled OLS | 0.6060 | +0.5560 | 0.5575 | 0% |
| Two-way fixed effects | 0.0495 | −0.0005 | 0.0099 | 100% |
| DiD vs never-treated | 0.0482 | −0.0018 | 0.0148 | 95% |
| ITS level change | 0.0510 | +0.0010 | 0.0128 | 100% |

When industry targets persistently high-volume prescribers, the naive estimate
is **12× the true effect**, it reports +0.61 when the truth is +0.05, and its
confidence interval never contains the truth. The panel designs are all fine,
because the confounder is time-invariant and gets differenced away.

This is the reassuring case, and it is the one most commonly assumed.

## Dynamic targeting: everything fails

| Estimator | Mean estimate | Mean bias | RMSE | 95% coverage |
| --- | ---: | ---: | ---: | ---: |
| Naive pooled OLS | 0.6117 | +0.5617 | 0.5628 | 0% |
| Two-way fixed effects | 0.0984 | +0.0484 | 0.0499 | 0% |
| DiD vs never-treated | −0.0065 | −0.0565 | 0.0585 | 0% |
| ITS level change | 0.1104 | +0.0604 | 0.0626 | 0% |

When industry targets physicians whose prescribing is *already rising*, every
estimator fails, and they fail in different directions:

- Two-way fixed effects **roughly doubles** the true effect (+0.098 vs +0.050).
- ITS is similar (+0.110), for the same reason: the pre-onset trend is
  estimated from treated units, so the level shift absorbs a pre-existing rise.
- DiD against never-treated **erases the effect entirely**, returning −0.007.
  It anchors on period $g-1$, which under dynamic targeting is a transitory
  peak; subsequent mean reversion cancels the real effect. This is an
  Ashenfelter dip inverted, selection on a temporary high rather than a
  temporary low.

**Coverage is 0% for all four.** Every design reports a tight confidence
interval that excludes the truth. Nothing in the point estimate or its
standard error signals a problem.

## The diagnostics do work

What fails silently in the estimates shows up loudly in the diagnostics:

| Regime | Pre-trend test | Placebo (onset −2) | Pre-period growth gap |
| --- | --- | ---: | ---: |
| `none` | PASS (p = 0.13) | +0.018, CI includes 0 | ~0.000 |
| `static` | PASS (p = 0.17) | +0.020, CI includes 0 | +0.004 |
| `dynamic` | **FAIL (p = 4×10⁻¹⁰⁷)** | +0.029, CI excludes 0 | **+0.022** |

The joint pre-trend test separates the identified cases from the unidentified
one decisively. The balance table shows the mechanism directly: under dynamic
targeting, physicians who will later receive payments were already growing at
0.032 log points per period before any payment, against 0.010 for those who
never receive one.

The placebo is unbiased under random assignment, verified across 30
replications, mean +0.0006 against a Monte Carlo standard error of 0.0013.

## What this implies for the literature

The systematic review in *Annals of Internal Medicine* found 21 of 36 studies
at serious risk of bias, and identified the mechanism precisely: dose-response
patterns "may also reflect residual confounding if industry targets clinicians
who already have higher baseline prescribing volumes."

These results sharpen that in a specific way. If targeting keys on *levels*,
panel methods handle it and the published within-physician designs are
defensible. If targeting keys on *momentum*, they do not, and the direction
of the error depends on the design, so a set of studies using different methods
will not converge on the truth. They will disagree, and each will look
internally precise.

Which regime holds is an empirical question about how manufacturers actually
target, and commercial analytics practice, dynamic cohorts updating on
incoming prescription data, points toward momentum rather than levels.

The practical recommendation is narrow and testable: **report the pre-trend
test and the placebo alongside any estimate.** They cost nothing, they are
computable on the same data, and they distinguish the regime where the estimate
means something from the regime where it does not.

## Limitations

- The simulation's dynamic-targeting parameters are chosen to make the problem
  visible, not calibrated to observed manufacturer behaviour. The direction of
  each failure is a property of the design; the magnitudes are not estimates of
  real-world bias.
- The real-data pipeline is verified against live CMS endpoints and unit-tested
  on synthetic frames, but a full Open Payments × Part D panel has **not** been
  built and analysed here. Doing so is the obvious next step.
- No estimator here solves the dynamic case. Instrumental variables and
  policy natural experiments, both named by the review, are not implemented.
