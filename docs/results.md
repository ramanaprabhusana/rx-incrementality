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

| Outcome | Estimate | 95% CI |
| --- | ---: | ---: |
| Claims for the promoted drug | **+4.0%** | 3.6 to 4.3% |
| Prescribes the drug at the reporting floor | **+2.5 points** | 2.4 to 2.7 (about 9% of a 29.3% base) |
| Claims per log(1 + payment dollars) | +1.6% | 1.5 to 1.7% |
| Claims, controlling for next year's payment | **+2.5%** | 2.2 to 2.9% (conservative) |

The honest range for the effect is about **2.5 to 4.0%**.

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="figures/event-study-dark.png">
  <img alt="Event study: estimates near zero in the years before the first payment, rising to about 7.5 percent two years after" src="figures/event-study-light.png" width="760">
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

## Event study

Clean onsets only, 3-way fixed effects, physicians, reference year t-1.

| Years from first payment | Cohort estimator | Pooled dummies |
| ---: | ---: | ---: |
| -5 | -4.8% | |
| -4 | -0.4% | |
| -3 | -0.8% | -2.1% |
| -2 | -0.5% | -0.9% |
| 0 | +2.0% | +2.1% |
| +1 | +6.9% | +6.8% |
| +2 | **+7.5%** | +7.9% |
| +3 | +5.3% | +6.1% |
| +4 | +3.0% | |
| Joint pre-trend test | p = 0.33 | p = 0.06 |

The **cohort estimator** (Sun and Abraham, 2021) fits each onset cohort
separately against never-paid pairs and aggregates with explicit cohort weights.
Pooled relative-time dummies average across cohorts with weights that can be
contaminated when effects build over time, which they do here. On simulated data
with growing, cohort-specific effects the pooled version manufactures a lead
that is not there; on the real data it shows a -2.1% lead at t-3 (t = -2.3) that
the cohort estimator does not (-0.8%, SE 0.9%). The t-5 estimate comes from the
2024 cohort alone and is correspondingly noisy.

The effect builds over two years and then fades: +2.0% in the onset year, +7.5%
two years later, +3.0% by year four.

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
change and the robustness checks the documentation had promised but not run).

## What these numbers do not support

- **Causation beyond doubt.** Pre-trends pass and the falsification test nearly
  does, but the residual next-year coefficient (+0.006, SE 0.002) shows a little
  anticipatory targeting remains. Drug-specific targeting cannot be excluded.
- **The whole market.** Medicare Part D only, which matters for GLP-1 agents given
  substantial commercial and cash-pay use.
- **Small prescribers.** Cells under 11 claims are suppressed; the intensive
  margin is conditional on clearing that floor.
- **Precise timing.** Annual data cannot order a payment and prescribing within a
  year.
- **Peer effects.** A city is too crude a network to say anything about them.

Standard errors are clustered by physician, with nested fixed effects excluded
from the small-sample correction. With 15 families, clustering by drug or
two-way would rest on too few clusters to be reliable.
