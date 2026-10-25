# Replication: oral anticoagulants, 2019 to 2024

Every number here comes from `scripts/replicate_class.py --class anticoagulant`,
written to `results/replication_anticoagulant.json`. The script runs the
specifications the diabetes conclusions rest on, chosen before it was run, with
the same code and thresholds. The diabetes results are unchanged.

## Why this class

Oral anticoagulants are the class Agha and Zeltzer (2022) studied for 2014 to
2016, when the direct oral anticoagulants were new. By 2019 to 2024 the market
had matured: Eliquis dominates, Xarelto is second, Pradaxa is small, and generic
warfarin is the unpromoted alternative. Promotion also tapered, from 490,509
Open Payments records naming these brands in 2019 to 208,353 in 2024.

Prior work on this class: Fleischman et al. (BMJ 2016) related anticoagulant
payments to regional prescribing for 2013 to 2014. Annapureddy et al. (American
Heart Journal 2025), using 2017 cardiology registry data (NCDR PINNACLE), found
no overall association between DOAC-related payments and DOAC prescribing (odds
ratio 0.99, 95% CI 0.95 to 1.03), only among physicians receiving higher-value
payments. As far as a literature search found, no study applies a
within-physician design to this class for 2019 to 2024, so this is an
extension, not a first.

## Data

| | Anticoagulants | Diabetes (main analysis) |
| --- | ---: | ---: |
| Physicians, observed every year | 106,014 | 71,869 |
| Brand families above 100,000 peak annual claims | 4 | 15 |
| Physician-family-year cells | 2,544,336 | 6,180,734 |
| Cells at the 11-claim reporting floor | 1,531,826 | 1,813,546 |
| Cells with a payment | 10.6% | 13.1% |

The families are Eliquis, Xarelto, Pradaxa and warfarin (with Jantoven).
Savaysa falls below the volume threshold. Open Payments records without an NPI
are a larger share here: 1,989 of 1,988,835 records (0.10%), but 5.57% of
dollars.

## Results

Physicians, fixed effects for physician-by-year, drug-by-year and
physician-by-drug, standard errors clustered by physician.

| Outcome | Anticoagulants | 95% CI | Diabetes |
| --- | ---: | ---: | ---: |
| **Claims for the promoted drug** | **+1.9%** | 1.6 to 2.3% | +4.0% |
| Same, controlling for next year's payment | +1.9% | 1.5 to 2.2% | +2.5% |
| Same, pairs above the reporting floor every year | +1.7% | 1.3 to 2.1% | +4.1% |
| Prescribes the drug at all (11-claim floor) | +1.23 points | 1.04 to 1.42 | +2.5 points |
| Medicare gross spending on the drug | +1.7% | 1.4 to 2.0% | +3.7% |
| Promoted drugs only, warfarin excluded | +1.3% | 1.0 to 1.7% | |

### Targeting bias replicates, and is larger

| Fixed effects | Effect | Next year's payment (log points, should be 0) |
| --- | ---: | ---: |
| Physician-year + drug-year | +43.4% | +0.249 (SE 0.003) |
| Physician-drug + drug-year (Carey, Lieber and Miller) | +5.5% | +0.028 (SE 0.002) |
| **All three** | **+1.9%** | **+0.005 (SE 0.002)** |

Without physician-by-drug effects the estimate is 22 times too large, against six
times for diabetes: representatives in this market select physicians who already
favour their drug even more strongly. The published two-way design again fails
the falsification test on annual data, and all three sets of fixed effects again
bring the next-year coefficient close to zero.

### The event study does not show an effect building after first payment

Clean onsets only, cohort estimator (Sun and Abraham, 2021), reference year t-1.

| Years from first payment | Claims |
| ---: | ---: |
| -5 | +2.0% |
| -4 | -0.1% |
| -3 | +0.1% |
| -2 | -0.1% |
| 0 | -0.1% |
| +1 | -0.7% |
| +2 | -1.0% |
| +3 | -1.4% |
| +4 | -2.9% |
| Joint pre-trend test | p = 0.76 |

Pre-trends pass, but unlike diabetes (+2.0% in the onset year rising to +7.5%),
prescribing does not rise after a physician's first payment about a blood
thinner; if anything it drifts down. The static estimate compares paid and
unpaid years within physician-drug pairs, so it also draws on pairs whose payment
status switches on and off within an established relationship. The event study
uses only first payments. A positive static estimate with a flat event study
fits representatives visiting in the years a physician prescribes more, more
than it fits a lasting effect of a new relationship.

**So the anticoagulant evidence supports a small association, about 1.9%, not
a demonstrated causal effect.** The diabetes result is the stronger of the two:
there the event study shows the effect building after first payment.

### Payment type and substitution

- A speaker, consulting or travel payment goes with +5.5% (SE 0.7%), against
  +1.9% for meals only, a wider gap than in diabetes (+4.9% against +4.0%).
- Excluding warfarin lowers the estimate to +1.3%: part of the measured effect
  is patients moved from generic warfarin to the promoted drug, which a
  within-physician comparison counts on both sides.

## What a payment dollar buys

| Per payment dollar, Medicare gross spending, 2019 to 2024 | Anticoagulants | Diabetes |
| --- | ---: | ---: |
| Payments to these physicians about these drugs | $59.1M | $185.0M |
| Medicare gross spending on the paid drug, in paid years | $27.1B | $29.7B |
| **Naive: that spending per payment dollar** | **$459** | $161 |
| **Spending caused, per payment dollar** | **$7.70** ($7.53 conservative) | $5.76 ($3.54 conservative) |

The naive ratio overstates the return 60 times, and it rises from $169 in 2019 to
$1,300 in 2024 as promotion tapered while spending on these drugs grew. The
causal figure rests on the static estimate, which the event study does not
support for this class, so treat it as an upper end. The caveats of the
diabetes figure apply: gross of rebates, Medicare only, excludes the cost of the
sales visit.

## Against the published literature

Agha and Zeltzer (2022) estimate about 10% for anticoagulants in 2014 to 2016
(working paper, before peer spillovers). In the mature market of 2019 to 2024
the within-physician estimate is about a fifth of that and lacks event-study
support. That is consistent with promotion mattering most while a class is
diffusing, though two periods and two designs cannot establish it.

The registry study above points the same way from a different design: no overall
association for DOACs, an association only with higher-value payments. Here,
too, speaker, consulting and travel payments go with a larger effect than meals.

## Sources

- [Agha and Zeltzer (2022), American Economic Journal: Economic Policy 14(2)](https://doi.org/10.1257/pol.20200044)
- [Annapureddy et al. (online 2025), American Heart Journal 291](https://doi.org/10.1016/j.ahj.2025.07.015)
- [Fleischman et al. (2016), BMJ 354](https://doi.org/10.1136/bmj.i4189)
- [Sun and Abraham (2021), Journal of Econometrics 225(2)](https://doi.org/10.1016/j.jeconom.2020.09.006)

## Limitations

- Four families, two of them dominant, give the within-physician comparison
  less variation than the fifteen diabetes families.
- Warfarin as the unpromoted comparator lets substitution enter the estimate
  (see above).
- Medicare only, annual timing, suppression below 11 claims, and peer effects
  unmeasured, as for diabetes.
