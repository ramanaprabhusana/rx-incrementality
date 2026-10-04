# Results on CMS data, 2019 to 2024

All numbers come from `scripts/estimate_real.py`, written to
`results/estimates.json`. Diabetes agents (SGLT2, DPP-4, GLP-1), Medicare Part D
by Provider and Drug linked to Open Payments general payments by NPI.

## Panel

| | |
| --- | ---: |
| Physicians (observed every year, 2019 to 2024) | 89,011 |
| Brand families | 15 |
| Physician-family-year cells | 7,654,946 |
| Cells observed at the 11-claim reporting floor | 2,199,852 (28.7%) |
| Cells with a payment | 12.3% |
| Physician-family pairs with a clean payment onset (2020 or later) | 198,466 |
| Pairs paid in 2019, onset unknown (excluded from event study) | 166,642 |
| Pairs never paid | 970,057 |
| In-class payment dollars landing on panel physicians | 68.9% |

Families: Bydureon, Farxiga, Invokana, Janumet, Januvia, Jardiance, Mounjaro,
Onglyza, Ozempic, Rybelsus, Synjardy, Tradjenta, Trulicity, Victoza, Xigduo.
Fourteen low-volume families fall below the 100,000 peak-annual-claims
threshold.

## Headline

Fixed effects for physician-by-year, drug-by-year and physician-by-drug.
Standard errors clustered by physician.

| Outcome | Estimate | 95% CI | Reading |
| --- | ---: | ---: | --- |
| Log claims, observed cells | +0.0386 | 0.0357 to 0.0415 | **+3.9%** more claims for the promoted drug |
| Prescribes at the reporting floor | +0.0268 | 0.0254 to 0.0282 | +2.7 points, about 9% relative to a 28.7% base |
| Log claims per log(1 + payment $) | +0.0154 | 0.0144 to 0.0164 | doubling payment value adds about 1.1% |

## The specification is most of the answer

| Fixed effects | Effect | Future payment coefficient (should be zero) |
| --- | ---: | ---: |
| Physician-year + drug-year | +23.3% | +0.122 |
| Physician-drug + drug-year (Carey, Lieber and Miller) | +6.8% | +0.026 |
| **Physician-year + drug-year + physician-drug** | **+3.9%** | **+0.005** |

The future-payment coefficient is a falsification test: a payment in year t+1
cannot cause prescribing in year t, so a non-zero coefficient measures targeting
rather than effect.

- Omitting physician-drug effects inflates the estimate sixfold. Representatives
  target physicians who already favour their drug, the pattern the simulation's
  `affinity` regime reproduces, where the same omission overstates fourteenfold.
- The published specification fails the falsification test on these data
  (t = 16). Adding physician-year effects removes most of the failure and lowers
  the estimate by about 40%.

That second point **does not contradict Carey et al.** Their monthly data show
flat pre-trends, and at monthly granularity a lead is weeks, not a year. It
indicates that with annual public data, physician-level momentum has to be
absorbed explicitly. The 2-way and published-style rows were added after the
first full run, to answer what each set of fixed effects contributes, and are
labelled as such in the results file.

## Event study

Clean onsets only, 3-way fixed effects, intensive margin, reference year t-1.

| Years from first payment | Coefficient | SE | Percent |
| ---: | ---: | ---: | ---: |
| -3 | -0.0019 | 0.0075 | -0.2% |
| -2 | -0.0073 | 0.0041 | -0.7% |
| 0 | +0.0226 | 0.0032 | +2.3% |
| +1 | +0.0586 | 0.0043 | +6.0% |
| +2 | +0.0591 | 0.0052 | +6.1% |
| +3 | +0.0409 | 0.0063 | +4.2% |

Joint pre-trend test: chi-squared(2) = 3.6, p = 0.17. Only two leads are
testable, because onsets begin in 2020 and data in 2019.

## Robustness and heterogeneity

| Specification | Effect | 95% CI |
| --- | ---: | ---: |
| Paid last year instead of this year | +4.1% | |
| Paid this year, controlling for next year's payment | **+2.6%** | 2.2 to 2.9% |
| Controlling for same-city peer exposure | +3.8% | 3.5 to 4.1% |
| GLP-1 only | +4.5% | 3.8 to 5.1% |
| SGLT2 only | +3.7% | 3.0 to 4.4% |
| DPP-4 only | +0.0% | -0.7 to +0.8% |
| GLP-1, 2019 to 2021 | +0.7% | -0.4 to +1.8% |
| GLP-1, 2022 to 2024 | +5.0% | 4.2 to 5.8% |

**Residual anticipation.** The future-payment coefficient under 3-way effects is
+0.005 (SE 0.0018). Small, but not zero. Controlling for it gives +2.6%, which
is the conservative reading. The honest range for the effect is about 2.6 to
3.9%.

**Peers.** Full same-city peer exposure adds about 1.2% (SE 0.4%), and the own
effect barely moves. Agha and Zeltzer measured peer spillovers at about a
quarter of the total effect using shared-patient networks; a city is a much
cruder proxy, so this should not be read as contradicting them.

**DPP-4.** A precise zero. The class was mature and losing share throughout, and
promotion has no detectable effect on prescribing. Whether that reflects
saturation, the shift of promotion toward SGLT2 and GLP-1, or something else is
not identified here.

**The GLP-1 era split is suggestive, not clean.** The 2022 to 2024 effect
includes the Mounjaro launch. Launches are where drug-specific targeting, the one
selection pattern this design cannot absorb, is most plausible: representatives
seek out early adopters. The 2019 to 2021 subsample also has only three years of
within-pair variation, so less power.

## Against the published literature

| Study | Period | Setting | Estimate |
| --- | --- | --- | --- |
| Carey, Lieber and Miller (2021) | 2013 to 2015 | All promoted Part D drugs, monthly | +2.2% patients, +1.6% days supply, +5.2% spending, first six months |
| Mizik and Jacobson (2004) | | Detailing, three drugs | +3.6 to 11.8% new prescriptions |
| Agha and Zeltzer (2022) | 2014 to 2016 | Anticoagulants | about 10% (working paper); 23% including peer spillovers (published) |
| Grennan et al. (2018) | | Statins, cardiologists | about 70% for a meal |
| **This repository** | **2019 to 2024** | **Diabetes, annual public data** | **+3.9% claims (2.6% conservative); +2.3% rising to +6.0% by event time** |

The estimate sits at the low end of the range, in line with the most comparable
design. On the question this repository set out to answer, whether published
effect sizes hold into the GLP-1 era, the answer is approximately yes for GLP-1
and SGLT2 agents, and no detectable effect for DPP-4 agents.

## What these numbers do not support

- **Causation beyond reasonable doubt.** Pre-trends pass and the falsification
  test nearly passes, but drug-specific targeting cannot be excluded, and the
  residual lead suggests a little of it is present.
- **The whole market.** Medicare Part D only, which matters for GLP-1 agents given
  substantial commercial and cash-pay use.
- **All prescribers.** The panel is physicians observed in all six years. About
  31% of in-class payment dollars go to physicians outside it.
- **Small prescribers.** Cells under 11 claims are suppressed; the intensive
  margin is conditional on clearing that floor.
- **Precise timing.** Annual data cannot order a payment and prescribing within
  the same year.

Standard errors are clustered by physician. With 15 families, clustering by drug
or two-way would rest on too few clusters to be reliable.
