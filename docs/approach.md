# Approach

## Why the current estimators are not enough

Every estimator in `rxinc.estimators` compares **treated physicians to
untreated physicians**. That is the source of their vulnerability. If a
manufacturer engages Dr. X because Dr. X's prescribing is already climbing, no
amount of physician fixed effects repairs it, because the confounder moves with
time. The simulation results in [findings.md](findings.md) quantify the damage:
under trajectory-based selection, two-way fixed effects roughly doubles the true
effect, interrupted time series does the same, differences-in-differences
against never-treated erases it, and confidence interval coverage falls to zero
for all four.

## The core idea: compare drugs, not physicians

Open Payments attributes each payment to a specific product, and Part D reports
prescribing by drug. That makes a different comparison available.

Instead of asking whether paid physicians prescribe more than unpaid
physicians, ask whether **a physician prescribes more of the drug they were paid
about than of competing drugs in the same class that they were not paid about,
in the same year**.

## Specification

For physician `i`, drug `j`, year `t`:

```
log(claims_ijt) = tau * Paid_ijt + alpha_it + delta_jt + mu_ij + e_ijt
```

**Revised after the first real-data run.** The original specification omitted
`mu_ij`, the physician-by-drug effect. On CMS data that version reported a
+23% effect, and next year's payment predicted this year's prescribing almost as
strongly as this year's payment did (+0.122), which is impossible as an effect
and diagnostic of targeting. Representatives target physicians who already
favour their drug. A simulated `affinity` regime reproduces this: without
`mu_ij` the estimate is fourteen times the truth, with it the truth is recovered.
All three effects are now absorbed; see [results.md](results.md).

- `mu_ij` is a **physician-by-drug** fixed effect, absorbing a physician's
  persistent preference for a particular drug, which is exactly what
  representatives select on.

- `alpha_it` is a **physician-by-year** fixed effect, absorbing everything about
  that physician in that year: overall prescribing trajectory, panel growth,
  risk mix shift, retirement, practice change. Physician-level momentum, the
  thing that defeats every current estimator, is differenced out by
  construction.
- `delta_jt` is a **drug-by-year** fixed effect, absorbing national drug-level
  shocks: launches, patent expiry, guideline changes, safety warnings,
  competitor entry.
- `tau` is identified only from variation within a physician-year, across drugs.

This is a triple difference, the within-firm across-product identification
standard in industrial organisation.

**This design is not new, and this repository did not originate it.** Carey,
Lieber and Miller (2021, Journal of Public Economics) use fixed effects for each
physician-drug combination on monthly Open Payments and Part D enrollee data for
2013 to 2015, motivated by the same targeting problem, and already report the
pre-trend result. Their monthly granularity is strictly better than the annual
public files used here. See [related-work.md](related-work.md) for what is and
is not contributed here.

## Verified feasibility

The design requires two kinds of within-physician variation. Both were measured
on real CMS data before committing to the approach.

**Prescribing variation.** Indiana direct oral anticoagulants, DY2024, 8,410
prescriber-drug rows across 5,695 physicians:

| Competing drugs prescribed | Prescribers | Share |
| --- | ---: | ---: |
| 1 | 2,999 | 52.7% |
| 2 | 2,678 | 47.0% |
| 3 | 17 | 0.3% |
| 4 | 1 | 0.0% |

47.3% prescribe two or more competitors.

**Exposure variation.** National diabetes payments, 2023, 948,267 payment
records across 132,756 physicians:

| Distinct brands paid about | Physicians | Share |
| --- | ---: | ---: |
| 1 | 53,675 | 40.4% |
| 2 | 27,693 | 20.9% |
| 3 | 17,530 | 13.2% |
| 4 | 14,035 | 10.6% |
| 5 or more | 19,823 | 14.9% |

59.6% are paid about two or more brands.

Both margins have ample variation, so the design is identified in practice and
not merely in principle.

## Defining treatment

The payment size distribution rules out a naive binary indicator. Physician-
brand-year totals for 2023:

| Statistic | Value |
| --- | ---: |
| Median | $33.32 |
| 75th percentile | $69.41 |
| 90th percentile | $124.99 |
| 99th percentile | $332.01 |
| Maximum | $138,520.53 |

97% of records are "Food and Beverage." The typical exposure is two or three
sponsored meals a year. A binary "any payment" variable would label nearly every
engaged physician as treated and discard the intensity that actually varies.

Treatment is therefore specified as `log(1 + amount)` for the primary estimate,
with threshold-based binary variants reported as robustness, and the small tail
of speaking and consulting relationships analysed separately since those differ
in kind, not only degree.

## Remaining threats

1. **Drug-specific targeting.** A manufacturer may pay Dr. X about drug A
   precisely because Dr. X's prescribing of drug A specifically is rising.
   `alpha_it` does not absorb this. Test with a physician-drug level event study
   and pre-trend test, reusing the existing diagnostics.
2. **Spillover.** Two kinds, and the second is not hypothetical. If engagement
   about drug A suppresses drug B, then B is not a clean control and tau is
   overstated; testable by checking whether payments about drug A move the same
   manufacturer's unrelated products. Separately, Agha and Zeltzer (2022, AEJ:
   Policy) measured **peer** spillovers in anticoagulants and found they account
   for roughly a quarter of the total effect. Untreated physicians are therefore
   not a clean control either, and a design ignoring this is mis-specified.
3. **Annual periodicity.** Part D is annual, so payments and prescribing within
   the same year cannot be ordered. This genuinely weakens short-run event
   studies and cannot be fixed with these data.
4. **Crosswalk error.** Mismatched drug names attenuate tau toward zero. Match
   rates get reported, not buried.

## Analysis window

Part D by Provider and Drug covers data years 2013 to 2024. The Open Payments
API metastore exposes general payment data for **2019 to 2025 only**, so the
usable overlap is **2019 to 2024, six annual periods**. Earlier Open Payments
years exist historically but are not retrievable through this API, so treatment
cannot be defined before 2019 and Part D years before 2019 are not usable here.

## Therapeutic class

Diabetes agents (SGLT2 inhibitors, DPP-4 inhibitors, GLP-1 receptor agonists),
29 brand families, of which 15 clear the volume threshold used in estimation. Chosen over anticoagulants because six manufacturers compete at
meaningful volume (Lilly, Novo Nordisk, AstraZeneca, Boehringer Ingelheim,
Merck, Janssen) rather than two, which yields more independent
payment-to-product relationships. Anticoagulants reduce in practice to Eliquis
against Xarelto, with Pradaxa and Savaysa returning 20 and 5 Indiana rows
respectively.

Filtering to one class also makes the data tractable: Part D by Provider and
Drug is 28.0 million rows per year unfiltered and roughly 900,000 for this
class.

## Pipeline

| Stage | Work | Status |
| --- | --- | --- |
| 1. Scope | Select therapeutic class | Done: diabetes |
| 2. Acquire Part D | Brand families, 2019 to 2024 | Done: 3,715,046 rows |
| 3. Acquire payments | Stream and filter by brand family | Done: about 5.49M records |
| 4. Link | NPI join | Done |
| 5. Crosswalk | Brand families on both sides | Done: 28 of 29 families matched |
| 6. Panel | Physician x family x year | Done: 7,654,946 cells |
| 7. Estimate | Three-way fixed effects | Done: [results.md](results.md) |
| 8. Validate | Simulated ground truth, four regimes | Done |

Stage 8 came before believing any real-data number, and then had to be revisited after it. The same
discipline applied to the existing estimators applies here: prove the design
recovers a known effect on simulated data where truth exists, then point it at
CMS data where it does not.
