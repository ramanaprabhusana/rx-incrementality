# Methodology

## The estimand

Let $y_{it}$ be log Medicare Part D claims for physician $i$ in period $t$, and
$D_{it}=1$ once physician $i$ has entered a payment relationship with a
manufacturer. The target is the average effect of that relationship on
prescribing volume:

$$\tau = E\left[y_{it}(1) - y_{it}(0)\right]$$

Because payments are not randomly assigned, $\tau$ is identified only under
assumptions about *how* manufacturers choose whom to pay.

## Three targeting regimes

The simulator implements the three cases that matter, because they differ in
which assumption they break.

| Regime | Assignment rule | What breaks |
| --- | --- | --- |
| `none` | Payments land at random | Nothing. A sanity check. |
| `static` | Target persistently high-volume prescribers | Cross-sectional comparison |
| `dynamic` | Target physicians whose Rx is already rising | Fixed effects as well |

Under `static`, the confounder is the physician effect $\alpha_i$, which is
time-invariant and therefore absorbed by a physician fixed effect. Under
`dynamic`, assignment depends on the recent path of the outcome itself. The
confounder now moves over time and correlates with the error, so differencing
out $\alpha_i$ is no longer enough.

The distinction is not academic. Industry targeting explicitly incorporates
prescribing momentum: commercial analytics platforms score prescribers on
propensity and switching likelihood using dynamic cohorts that update as new
prescription data arrives. That is `dynamic`, not `static`.

## Estimators

**Naive pooled OLS.** Regress $y_{it}$ on $D_{it}$ with an intercept. Every
persistent difference between paid and unpaid physicians loads onto the
coefficient. Included because it is what an uncontrolled comparison reports.

**Two-way fixed effects.** Absorb physician and period effects. Implemented by
alternating projections rather than explicit dummies; the transform iterates
demeaning by unit and by period until both margins are within tolerance of
zero, which makes it exact on unbalanced panels too. `tests/test_estimators.py`
checks it against an explicit dummy regression in both the balanced and
unbalanced case.

**DiD against never-treated.** For each onset cohort $g$ and period $t \ge g$:

$$\widehat{ATT}(g,t) = \left[\bar y^{g}_{t} - \bar y^{g}_{g-1}\right] - \left[\bar y^{never}_{t} - \bar y^{never}_{g-1}\right]$$

aggregated with cohort-size weights. This avoids using already-treated
physicians as controls, which is the Goodman-Bacon problem in staggered
two-way fixed effects. Standard errors come from a physician-level bootstrap.

**Event study.** Relative-time dummies around onset, reference period $-1$,
endpoints binned. The post-onset coefficients are not the interesting part.
The *pre*-onset coefficients are the only quantity in the package that can
detect dynamic targeting from the data alone.

**Interrupted time series.** Within treated physicians, a pre-onset trend, a
level shift at onset, and a post-onset slope change. Named by the 2021
systematic review as a design that moves past pure association.

## Inference

All standard errors cluster on physician. The idiosyncratic shock is AR(1), so
observations within a physician are dependent and unclustered errors
substantially overstate precision. `test_clustered_se_exceeds_classical_under_within_unit_correlation`
pins this down.

The cluster-robust covariance uses the standard sandwich with a finite-sample
correction $\frac{G}{G-1}\cdot\frac{N-1}{N-K}$. For the single-outcome
estimators, $K$ counts parameters removed by the within transform as well as
those explicitly estimated. For the drug-level estimators, fixed effects nested
within clusters are excluded from $K$, the convention of high-dimensional
fixed-effect estimators such as reghdfe: physician-by-year and physician-by-drug
effects sit inside physician clusters, and counting them would make $K$ exceed
$N$ on the real panel.

## The drug-level design

The design applied to CMS data compares drugs within a physician-year:

$$\log(\text{claims}_{ijt}) = \tau\, \text{Paid}_{ijt} + \alpha_{it} + \delta_{jt} + \mu_{ij} + \varepsilon_{ijt}$$

with physician-by-year $\alpha_{it}$, drug-by-year $\delta_{jt}$ and
physician-by-drug $\mu_{ij}$ effects, absorbed by alternating projections over
any number of factors (`rxinc.estimators._absorb`, checked against explicit
dummy regression with three overlapping factors on an unbalanced panel).
$\alpha_{it}$ absorbs physician momentum, $\mu_{ij}$ absorbs a physician's
persistent preference for a drug, which is what representatives target, and
$\delta_{jt}$ absorbs national drug shocks. The remaining threat is targeting on
a physician's rising use of one specific drug, which varies within the
physician-year.

## Event studies with staggered onsets

`drug_event_study` fits pooled relative-time dummies. With staggered onsets and
effects that build over time or differ across cohorts, pooled dummies average
cohort-specific effects with weights that can be contaminated, so even
pre-period coefficients can pick up other cohorts' post-period effects.
`cohort_event_study` implements the interaction-weighted estimator of Sun and
Abraham (2021): one coefficient per (onset cohort, relative period) against
never-treated pairs, aggregated per relative period with weights equal to each
cohort's share of that period's treated observations, and delta-method standard
errors with weights held fixed. Relative periods no pair reaches are dropped
rather than reported as zero. Pairs whose onset is unknown must be removed first,
or they would enter as never treated.

## When a pre-trend test fails

`detrend_event_study` fits a line through the lead coefficients, constrained
through zero at the reference period and weighted by their full covariance,
extrapolates it, and subtracts it from every coefficient, with delta-method
standard errors that include the uncertainty in the slope. It is a sensitivity
analysis under the assumption that the pre-existing trend would have continued
linearly, not identification. Rambachan and Roth (2023) formalise how far that
assumption can be relaxed.

## Exposure that is unobserved, not zero

Open Payments covers non-physician practitioners only from 2021, and never covers
pharmacists or physicians in training. `rxinc.drugdata` classifies every
prescriber, marks payment exposure missing before coverage begins, and moves the
left-censoring year to the first covered year. Estimation drops cells with
missing exposure rather than treating them as unpaid.

## Why simulate at all

On real data $\tau$ is unknown, so "how biased is this estimator" is
unanswerable. The simulation supplies a known $\tau$ and a targeting rule you
control, which converts that unanswerable question into an arithmetic one. The
panel schema is identical to the one :func:`rxinc.panel.build_panel` produces
from CMS data, so the same estimator code runs on both.

Claims about bias are made from `monte_carlo`, which repeats the exercise over
fresh seeds. A single panel cannot distinguish bias from sampling noise, an
early version of this repository reported a spurious placebo effect that
disappeared under replication.
