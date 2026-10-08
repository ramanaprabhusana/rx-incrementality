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

## Bounding parallel-trends violations (Rambachan and Roth, 2023)

`rxinc.sensitivity` implements the smoothness restriction $\Delta^{SD}(M)$: the
differential trend $\delta$ satisfies $\delta_{-1} = 0$ and
$|\delta_{t+1} - 2\delta_t + \delta_{t-1}| \le M$ for every period. Inference
uses the fixed-length confidence interval the paper recommends for this class.
The estimator $c'\hat\beta$ fixes the post-period weights to the target and
chooses pre-period weights; any straight line through the reference period
must leave it unbiased. Its half-length is
$\text{sd}(c)\cdot cv_\alpha(\bar b(c)/\text{sd}(c))$, where $\bar b(c)$ is
the worst-case bias, a linear program, and $cv_\alpha(b)$ the $1-\alpha$
quantile of $|N(b,1)|$.

### Choosing the weights

The half-length is convex in the free pre-period weights: the bias and standard
deviation are convex, and $s \cdot cv_\alpha(B/s)$ is jointly convex and
non-decreasing in both arguments, which `tests/test_sensitivity.py` checks
numerically. The optimiser exploits this in three steps:

1. **Minimum achievable bias width**, by a linear program in the weights and the
   dual variables of both bias programs.
2. **Minimum variance for each width budget**, by cutting planes on the 4 or 5
   free weights: whenever the current weights exceed the budget, the extreme
   violations form a linear constraint, and a small quadratic program is
   re-solved. Cuts depend only on the restriction, so they are shared across
   budgets. A one-dimensional search over the budget, anchored at the exact
   minimum width, picks the shortest interval.
3. **A local polish** from the best candidate, safe because the objective is
   convex.

Two earlier optimisers failed, and the tests now guard against both. Nelder-Mead
on the weights stalled and produced an interval that narrowed as M grew, which
is impossible at the optimum and would corrupt the breakdown search. A quadratic
program over the dual variables failed whenever the budget stopped binding,
because the duals then have no unique value: on the claims event study it failed
at 17 of 24 budgets and missed the optimum by 1.6%. Its general version also had
a redundant equality, one per bias program, which SLSQP cannot handle, and took
576 seconds per interval. The current optimiser takes 0.2 to 0.5 seconds, is never
longer than the dual version by more than five parts in a hundred million, and is
1.6% shorter where the dual version failed. Only one published number moved: the
claims five-year-average breakdown, from 0.15 to 0.17.

The bias programs are solved with free variables. Restrictions that leave the
slope free admit an unbounded ray along the straight-line direction; trend-removing
weights are orthogonal to it only to rounding, which occasionally makes the solver
fail. Only then is the program re-solved inside a box sized to the problem. A
fixed box of plus or minus 10,000 was tried first and cost about $10^{-12}$ of
accuracy on every problem, enough for a reported worst-case bias to fall short of
an attainable one.

### One-sided restrictions

`restricted_ci` adds shape and sign restrictions: the trend can only flatten
(concave), only steepen (convex), only rise or only fall (Rambachan and Roth's
monotonicity classes), or the bias can only be positive or negative. These make
the bias range $[b_{lo}, b_{hi}]$ asymmetric, and the interval becomes
$[c'\hat\beta - b_{hi} - x\,\text{sd},\; c'\hat\beta - b_{lo} + x\,\text{sd}]$
with $x$ solving $\Phi(x + w) - \Phi(-x) = 1 - \alpha$, $w = (b_{hi} - b_{lo})/\text{sd}$.
Its coverage is at least $1-\alpha$ everywhere and exactly $1-\alpha$ at either
end of the bias range. For a symmetric range it reduces algebraically to the
fixed-length interval, which the tests confirm to machine precision.

`pre_period_support` decides whether the leads permit a restriction: a
significant bend of the wrong sign contradicts it; support requires a
significant bend of the right sign and at least half the point estimates
agreeing. The middle case, one significant bend while most estimates disagree, is
reported as mixed rather than as support.

### Validation

In `tests/test_sensitivity.py`:

- **Exact coverage.** For any given violation the estimator is normal with known
  bias, so coverage is computed exactly, not simulated, for the two-sided bound
  and for every one-sided restriction. Coverage is at least 95% for every
  violation tested and exactly 95% at the worst case, so the worst-case bias is
  correct and attained. An earlier Monte Carlo version returned 92.7% at M = 0
  from noise; 3,000 draws gave 95.1%.
- **Monotonicity.** Interval length is non-decreasing in M, two-sided and
  one-sided.
- **Optimality.** No random trend-removing weights produce a shorter interval.
- **Consistency.** At M = 0 the estimate equals straight-line detrending, and the
  one-sided construction equals the fixed-length one when symmetric.
- **Convexity.** The critical-value function is convex and its perspective
  non-decreasing in the standard deviation.

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
