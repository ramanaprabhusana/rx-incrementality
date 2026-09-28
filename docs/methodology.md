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
correction $\frac{G}{G-1}\cdot\frac{N-1}{N-K}$, where $K$ counts parameters
removed by the within transform as well as those explicitly estimated.

## Why simulate at all

On real data $\tau$ is unknown, so "how biased is this estimator" is
unanswerable. The simulation supplies a known $\tau$ and a targeting rule you
control, which converts that unanswerable question into an arithmetic one. The
panel schema is identical to the one :func:`rxinc.panel.build_panel` produces
from CMS data, so the same estimator code runs on both.

Claims about bias are made from `monte_carlo`, which repeats the exercise over
fresh seeds. A single panel cannot distinguish bias from sampling noise — an
early version of this repository reported a spurious placebo effect that
disappeared under replication.
