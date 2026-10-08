# Results on CMS data, 2019 to 2024

Every number here comes from `scripts/estimate_real.py`, written to
`results/estimates.json`. The exact input files, with row counts and SHA-256
hashes, are in `results/provenance.json`. Diabetes agents (SGLT2, DPP-4, GLP-1),
Medicare Part D by Provider and Drug linked to Open Payments general payments by
NPI and brand family.

## Populations

Open Payments did not cover nurse practitioners, physician assistants and other
non-physician practitioners until program year 2021: the extracts contain zero
such records in 2019 and 2020, then about 360,000 a year. Their earlier exposure
is unknown, not zero. Pharmacists and physicians in training are never covered,
so their exposure carries no information at all. The primary population is
therefore **physicians**, whose exposure is measured in every year.

| | Physicians (primary) | Non-physician practitioners |
| --- | ---: | ---: |
| Prescribers, observed every year | 71,869 | 16,988 |
| Years with payment exposure observed | 2019 to 2024 | 2021 to 2024 |
| Prescriber-family-year cells | 6,180,734 | 1,002,292 with exposure observed |
| Cells at the 11-claim reporting floor | 1,813,546 (29.3%) | 383,027 |
| Cells with a payment | 13.1% | 13.4% |
| Pairs with a clean payment onset | 139,169 | 22,696 |
| Pairs paid in the first observed year, onset unknown | 166,620 | 36,544 |
| Pairs never paid | 772,246 | 195,580 |

Fifteen brand families clear the 100,000 peak-annual-claims threshold:
Bydureon, Farxiga, Invokana, Janumet, Januvia, Jardiance, Mounjaro, Onglyza,
Ozempic, Rybelsus, Synjardy, Tradjenta, Trulicity, Victoza, Xigduo.

## Headline

Physicians. Fixed effects for physician-by-year, drug-by-year and
physician-by-drug. Standard errors clustered by physician.

| Outcome | Estimate | 95% CI | Pre-trend test |
| --- | ---: | ---: | --- |
| **Claims for the promoted drug** | **+4.0%** | 3.6 to 4.3% | passes, p = 0.33 |
| Same, pairs above the reporting floor every year | +4.1% | 3.6 to 4.7% | |
| Same, controlling for next year's payment | +2.5% | 2.2 to 2.9% | |
| Claims per log(1 + payment dollars) | +1.6% | 1.5 to 1.7% | |
| Prescribes the drug at all (11-claim floor) | +2.5 points | 2.4 to 2.7 | **fails, p < 0.001** |

**The credible result is the intensive margin: about 2.5 to 4.0% more claims**
among physicians already prescribing the drug. It passes its pre-trend test, and
it is unchanged (+4.1%) on pairs above the reporting floor in every year, where
a payment cannot change who is observed, so selection through the floor is not
driving it.

**The extensive margin is not credible as a plain causal estimate.** In the
years before a physician is first paid about a drug, the probability that they
prescribe it was already climbing: -4.9, -2.2, -1.3 and -0.4 points at four,
three, two and one year before the reference year. Representatives target
physicians who are starting to adopt a drug, which is drug-specific targeting,
the one pattern this design cannot absorb. Allowing for that trend with
Rambachan and Roth bounds, the jump in the first two years after payment
survives trend bending as large as anything seen in the two years before onset,
but not the larger bending four years before. See the sensitivity section below.

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="figures/event-study-dark.png">
  <img alt="Left: claims flat before the first payment, rising to about 7.5 percent two years after. Right: the probability of prescribing at all was already rising before the first payment; net of that trend it still rises about 5 points" src="figures/event-study-light.png" width="760">
</picture>

## The specification is most of the answer

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="figures/specification-dark.png">
  <img alt="Effect falls from 23.6 percent to 7.5 to 4.0 percent as fixed effects are added; the impossible next-year effect falls from 13.0 to 2.9 to 0.6 percent" src="figures/specification-light.png" width="760">
</picture>

| Fixed effects | Effect | Next year's payment (log points, should be 0) |
| --- | ---: | ---: |
| Physician-year + drug-year | +23.6% | +0.122 (SE 0.002) |
| Physician-drug + drug-year (Carey, Lieber and Miller) | +7.5% | +0.028 (SE 0.002) |
| **All three** | **+4.0%** | **+0.006 (SE 0.002)** |

A payment in year t+1 cannot cause prescribing in year t, so its coefficient
measures targeting. Omitting physician-drug effects inflates the estimate about
sixfold, because representatives target physicians who already favour their
drug. On annual data the published specification also fails the test; adding
physician-year effects removes most of the failure and lowers the estimate by
about 45%. This does not contradict Carey et al., whose monthly data show flat
pre-trends; at monthly granularity a lead is weeks rather than a year.

## Event studies

Clean onsets only, 3-way fixed effects, physicians, reference year t-1.

| Years from first payment | Claims: cohort estimator | Claims: pooled dummies | Prescribes at all: cohort estimator |
| ---: | ---: | ---: | ---: |
| -5 | -4.8% | | -4.88 pt |
| -4 | -0.4% | | -2.19 pt |
| -3 | -0.8% | -2.1% | -1.26 pt |
| -2 | -0.5% | -0.9% | -0.39 pt |
| 0 | +2.0% | +2.1% | -0.07 pt |
| +1 | +6.9% | +6.8% | +6.71 pt |
| +2 | **+7.5%** | +7.9% | +7.93 pt |
| +3 | +5.3% | +6.1% | +7.01 pt |
| +4 | +3.0% | | +6.69 pt |
| Joint pre-trend test | p = 0.33 | p = 0.06 | **p < 0.001** |

The **cohort estimator** (Sun and Abraham, 2021) fits each onset cohort
separately against never-paid pairs and aggregates with explicit cohort weights.
Pooled relative-time dummies average across cohorts with weights that can be
contaminated when effects build over time, which they do here. On simulated data
with growing, cohort-specific effects the pooled version manufactures a lead
that is not there; on the real data it shows a -2.1% lead at t-3 (t = -2.3) that
the cohort estimator does not (-0.8%, SE 0.9%). The t-5 estimates come from the
2024 cohort alone.

### Net of the pre-trend

`rxinc.diagnostics.detrend_event_study` fits a line through the leads,
constrained through zero at the reference year and weighted by their full
covariance, extrapolates it, and subtracts it. It assumes the pre-existing trend
would have continued linearly; Rambachan and Roth (2023) formalise how far that
can be relaxed.

| Years from first payment | Claims, net of pre-trend | Prescribes at all, net of pre-trend |
| ---: | ---: | ---: |
| Pre-trend slope per year | +0.5% (SE 0.4), not significant | +0.75 pt (SE 0.11) |
| 0 | +1.5% | -0.83 pt |
| +1 | +5.9% | +5.20 pt |
| +2 | +6.0% | +5.67 pt |
| +3 | +3.4% | +4.00 pt |
| +4 | +0.6% | +2.93 pt |
| Post-payment average | **+3.5%** (SE 1.3) | **+3.40 pt** (SE 0.36) |

For claims the adjustment changes little, as it should when pre-trends are flat.
For adoption the leads flatten toward onset, which suggests adoption was
saturating; if so, a straight-line projection overstates the counterfactual trend
and the adjusted figures lean conservative.

## How far can parallel trends fail before the results do?

A pre-trend test cannot prove parallel trends. Following Rambachan and Roth
(2023), `rxinc.sensitivity` instead allows the differential trend to follow any
straight line through the reference year and lets its slope change by up to
**M** per year, then reports a 95% fixed-length confidence interval that holds
under every such violation. The **breakdown value** is the largest M at which the
interval still excludes zero. `make sensitivity` reproduces this from the saved
event-study covariances in `results/sensitivity.json`.

M is judged against how much the pre-payment trend itself bent:

| Pre-period bending (points per year squared) | Claims | Adoption |
| --- | ---: | ---: |
| Centred two years before onset | +0.14 (SE 1.08) | -0.48 (SE 0.30) |
| Centred three years before onset | +0.68 (SE 1.66) | -0.06 (SE 0.42) |
| Centred four years before onset | -4.92 (SE 2.87) | **-1.75 (SE 0.65)** |

Breakdown values, physicians:

| Target | Claims: estimate at M = 0 | Claims: breakdown M | Adoption: estimate at M = 0 | Adoption: breakdown M |
| --- | ---: | ---: | ---: | ---: |
| Onset year (t0) | +1.5% | 0.44 | -0.88 pt (a dip) | 0.16 |
| Year after (t+1) | +5.7% | 1.22 | +5.11 pt | 1.76 |
| Two years after (t+2) | +5.8% | 0.51 | +5.59 pt | 0.97 |
| Three years after (t+3) | +3.3% | 0.03 | +3.96 pt | 0.43 |
| Four years after (t+4) | +0.6%, not significant | 0 | +2.94 pt | 0.20 |
| **First two years** (headline) | **+3.6%** | **1.02** | **+2.12 pt** | **1.11** |
| All five years averaged | +3.4% | 0.17 | +3.35 pt | 0.51 |

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="figures/sensitivity-dark.png">
  <img alt="Confidence bands widening as the allowed bend M grows; the claims band reaches zero at M = 1.02 and the adoption band at M = 1.11" src="figures/sensitivity-light.png" width="760">
</picture>

The headline target, the first two years after first payment, was chosen before
looking at per-year breakdowns so as not to pick whichever year looked best.

- **Short-run effects survive bending like that seen just before payments.** For
  both outcomes the first-two-years breakdown (1.02 and 1.11) exceeds the bending
  in the two years before onset (at most 0.68).
- **Adoption does not survive the largest bending observed, if bends may go
  either way.** Four years before onset the adoption trend bent by -1.75, the
  only bend significant on a two-sided test, which exceeds its 1.11 breakdown.
  Only the year-after jump alone (1.76) survives it. But every adoption bend was
  a deceleration, and allowing for that changes the answer: see below.
- **Long-run effects rest on straight-line extrapolation.** Worst-case bias grows
  with the square of the distance from the reference year, so three and four
  years out the breakdown values fall to 0.03 and 0 for claims. The five-year
  average is fragile for the same reason.
- Claims' four-years-before bending (-4.92) comes from the 2024 cohort alone and
  is not significant on a two-sided test (SE 2.87, z = -1.7).

### If the trend can only flatten

Symmetric bounds let the trend bend up or down. The pre-period data say which
way it actually bent, so `rxinc.sensitivity` also supports one-sided
restrictions, imposed only where the leads do not contradict them:

| Pre-period bends point to | Claims | Adoption | Non-physician claims |
| --- | --- | --- | --- |
| Flattening (all bends negative) | **mixed**: one marginal bend from the noisiest lead, the two precise ones the other way | **supported**: all three negative, -1.75 significant | consistent, not supported |
| Steepening | contradicted | contradicted | consistent, not supported |
| Rising trend | supported | supported | consistent, not supported |
| Falling trend | contradicted | contradicted | consistent, not supported |

| First two years after payment | Claims | Adoption | Non-physician claims |
| --- | ---: | ---: | ---: |
| Breakdown, bends either way | 1.02 | 1.11 | 0.75 |
| Breakdown, trend can only flatten | never | **never** | never |
| Lowest lower limit at any M, can only flatten | +2.05% | **+1.67 pt** | +1.34% |
| Breakdown, trend can only rise | 1.02 | 1.11 | 0.75 |
| Breakdown, bias can only be positive | 1.02 | 1.11 | 0.75 |

*Breakdown values searched up to M = 50 points per year squared.*

- **The supported restriction settles adoption.** If the adoption trend could
  only keep flattening, as it did before payments, the first-two-years effect is
  at least +1.67 points at every M: no amount of flattening overturns it. A
  flattening trend lies below its straight-line projection, so the
  counterfactual can only fall further below the observed path. The assumption
  is that adoption did not re-accelerate on its own just when payments started.
- **Claims cannot lean on it.** The case for flattening there rests on a single
  marginal bend from the 2024 cohort alone, while the two precisely estimated
  bends point the other way, so the claims result stays at its two-sided
  breakdown of 1.02.
- **Sign restrictions change nothing.** "The trend only rises" and "the bias is
  only positive" are supported or plausible, but they constrain the upper end of
  the interval, not the lower, so every breakdown value is unchanged. This was
  expected: a rising trend means the raw estimate overstates the effect, so ruling
  out a falling one cannot help show the effect is positive.

The intervals are validated against exact coverage rather than simulation, for
the two-sided bound and for every one-sided restriction: for any allowed
violation the estimator is normal with known bias, and coverage is at least 95%
for every violation tested and exactly 95% at the worst case. See
[methodology.md](methodology.md).

## Robustness

| Check | Effect | 95% CI |
| --- | ---: | ---: |
| Paid last year rather than this year | +4.4% | 4.1 to 4.7% |
| Food and beverage only | +4.0% | 3.6 to 4.3% |
| Any speaker, consulting or travel payment | +4.9% | 3.4 to 6.4% |
| Payment of at least $25 | +5.4% | 5.0 to 5.8% |
| Payment of at least $100 | +4.9% | 4.4 to 5.5% |
| Controlling for same-city peer exposure | +3.9% | 3.6 to 4.3% |
| Controlling for payments about the same manufacturer's other drugs | +4.1% | 3.8 to 4.5% |
| Excluding Mounjaro | +4.0% | 3.6 to 4.3% |
| Physicians not required to appear every year (217,388) | +4.0% | 3.7 to 4.3% |
| Pairs above the reporting floor in every year | +4.1% | 3.6 to 4.7% |

**Money matters less than contact.** Speaker, consulting and travel payments are
under 2% of records but about half of all dollars, with a median near $940
against about $17 for a meal. Their effect, +4.9%, is barely above a meal's
+4.0%. That matches Carey et al.'s finding that results hold for small payments,
and suggests the effect runs through the relationship more than the transfer.

**Promotion cannibalises the manufacturer's own portfolio.** Being paid about
another of the same manufacturer's drugs lowers prescribing of this one by 0.8%
(SE 0.2%). It also means same-manufacturer comparison drugs are slightly
depressed when the focal drug is promoted, but controlling for it moves the
estimate from +4.0% to +4.1%.

**Peers, by city, show nothing.** Full same-city peer exposure adds +0.4% (SE
0.4%). Agha and Zeltzer found peer spillovers worth about a quarter of the total
effect using shared-patient networks; a city is far too coarse a proxy to
contradict that.

## Heterogeneity

| Group | Effect | 95% CI |
| --- | ---: | ---: |
| GLP-1 agents | +4.7% | 4.0 to 5.4% |
| SGLT2 inhibitors | +3.7% | 2.9 to 4.4% |
| DPP-4 inhibitors | -0.3% | -1.1 to +0.5% |
| GLP-1, 2019 to 2021 | +0.5% | -0.7 to +1.6% |
| GLP-1, 2022 to 2024 | +5.4% | 4.5 to 6.3% |
| GLP-1, 2022 to 2024, excluding Mounjaro | +5.9% | 4.9 to 6.9% |
| Family or general practice | +3.7% | 3.3 to 4.2% |
| Internal medicine | +4.1% | 3.6 to 4.6% |
| Endocrinology | +4.0% | 3.1 to 5.0% |
| Cardiology | +3.3% | 0.3 to 6.3% |
| Nephrology | +3.8% | -0.5 to +8.3% |
| Nurse practitioners, 2021 to 2024 | +3.1% | 2.2 to 4.1% |
| Physician assistants, 2021 to 2024 | +4.9% | 3.4 to 6.4% |

**The GLP-1 effect is concentrated in 2022 to 2024 and is not a launch
artifact.** Excluding Mounjaro raises it to +5.9%, so it is carried by Ozempic,
Rybelsus and Trulicity during the semaglutide surge. **DPP-4 inhibitors show a
precise zero.** The class was mature and losing share throughout; why promotion
does nothing there is not identified here.

**Specialties are strikingly similar**, between +3.3% and +5.2%, and
non-physician practitioners respond like physicians: +3.7% overall (95% CI 2.9
to 4.5%), with a pre-trend test that passes (p = 0.74).

## Against the published literature

| Study | Period | Setting | Estimate |
| --- | --- | --- | --- |
| Carey, Lieber and Miller (2021) | 2013 to 2015 | All promoted Part D drugs, monthly | +2.2% patients, +1.6% days supply, +5.2% spending, first six months |
| Mizik and Jacobson (2004) | | Detailing, three drugs | +3.6 to 11.8% new prescriptions |
| Agha and Zeltzer (2022) | 2014 to 2016 | Anticoagulants | about 10% (working paper); 23% including peer spillovers (published) |
| Grennan et al. (2018) | | Statins, cardiologists | about 70% for a meal |
| **This repository** | **2019 to 2024** | **Diabetes, annual public data** | **+4.0% claims (2.5% conservative); +2.0% in the onset year rising to +7.5%** |

The estimate sits at the low end of the range, beside the most comparable design.
Published magnitudes roughly hold into the GLP-1 era for GLP-1 and SGLT2 agents,
with no detectable effect for DPP-4 agents.

## What changed from the previous version

The previous version of this document reported **+3.9%** on 89,011 prescribers of
every type. That population mixed in nurse practitioners and physician
assistants, whose 2019 and 2020 payments Open Payments did not collect and the
pipeline recorded as zero, and pharmacists and physicians in training, who are
never covered. Restricting to physicians, whose exposure is measured throughout,
gives +4.0%. The estimate barely moved; the population it describes is now
correct. Non-physician practitioners are estimated separately on 2021 to 2024.

Specifications are labelled in `results/estimates.json` by when they were added:
round 1 before the first full run, round 2 after it (the fixed-effect
comparisons), round 3 after the coverage problem was found (the population
change and the robustness checks the documentation had promised but not run),
round 4 to address selection through the reporting floor (the always-observed
pairs and the extensive-margin event study, which is what exposed the adoption
pre-trend).

## What these numbers do not support

- **A clean causal effect on adoption.** The extensive margin fails its
  pre-trend test. Its short-run jump survives bending like that just before
  payment, but not the largest bending observed in the pre-period.
- **Effects three or more years out.** Under smoothness bounds they survive
  almost no departure from a straight-line trend.
- **Causation beyond doubt, even for claims.** Pre-trends pass and the
  falsification test nearly does, but the residual next-year coefficient
  (+0.006, SE 0.002) shows a little anticipatory targeting remains.
- **The whole market.** Medicare Part D only, which matters for GLP-1 agents given
  substantial commercial and cash-pay use.
- **Small prescribers.** Cells under 11 claims are suppressed; the intensive
  margin is conditional on clearing that floor.
- **Precise timing.** Annual data cannot order a payment and prescribing within a
  year.
- **Peer effects.** A city is too crude a network to say anything about them.
  Shared-patient network data for recent years exist (CareSet's DocGraph Hop
  Teaming dataset, editions through 2022) but are available only on request, not
  as a public download, so they are not used here.

Standard errors are clustered by physician, with nested fixed effects excluded
from the small-sample correction. With 15 families, clustering by drug or
two-way would rest on too few clusters to be reliable.
