# Related work, and what is not new here

This document exists because an earlier version of this repository implied the
design here was an advance on the literature. It is not. The core identification
strategy was published years ago, in better journals, on better data. Stating
that plainly is more useful than the alternative.

## The two literatures

Work on industry payments and prescribing splits into two bodies that rarely
cite each other, and conflating them produces a badly wrong picture.

**Clinical and health services research.** This is what the 2021 *Annals of
Internal Medicine* systematic review surveyed, finding 21 of 36 studies at
serious risk of bias and calling for stronger causal designs. That assessment is
accurate for the literature it covers.

**Health economics.** This literature had already answered that call. It is
small, technical, and published in economics journals, which is likely why the
clinical review did not reflect it.

Positioning a project against only the first literature, as this one originally
did, overstates the gap by a wide margin.

## Carey, Lieber and Miller (2021)

*Drug Firms' Payments and Physicians' Prescribing Behavior in Medicare Part D.*
Journal of Public Economics 197. NBER Working Paper 26751.

This paper already does what `rxinc.estimators.triple_diff` does, and does it
better:

- **Fixed effects for each physician-drug combination**, which is the same
  identification logic as the triple difference here.
- **Monthly** granularity, linking Open Payments 2013 to 2015 to prescription
  data for a large panel of Part D enrollees. They identify abrupt changes in
  prescribing in the months right after a payment.
- They explicitly motivate the design by the same problem this repository
  motivates it by: physicians who receive payments tend to be ex ante
  higher-volume prescribers.
- They **already report the pre-trend result**, finding no evidence of
  differential trends between paid and unpaid physicians before payment.

Their monthly data also dissolves the limitation this repository cannot escape.
Part D public files are annual, so payments and prescribing within a year cannot
be ordered. Carey et al. order them to the month.

Scale figures worth knowing from that paper: more than 85% of drug firms'
marketing expenditure targets physicians, more than one fifth of branded Part D
expenditure comes from a physician who recently received a payment for that
drug, and 29% of Part D physicians are paid for at least one drug over the
sample period.

They also go further than identification. Five case studies of patent expiry
show physicians receiving payments transition patients to generics just as
quickly as those who do not, and hand-collected efficacy data across three
therapeutic classes shows prescribed drug quality is largely unaffected by
payments.

## Agha and Zeltzer (2022)

*Drug Diffusion Through Peer Networks: The Influence of Industry Payments.*
American Economic Journal: Economic Policy 14(2), 1 to 33. NBER Working Paper
26338.

Relevant for two reasons.

First, it studies **anticoagulants**, the therapeutic class this project
originally selected before switching to diabetes. Over 2014 to 2016, payments
associated with anticoagulant marketing increased prescription volume by **23
percent**.

Second, it quantifies a channel this repository listed as an unaddressed threat.
**Peer spillovers contribute roughly a quarter of that increase**: prescribing
rises for the paid physician and for that physician's peers. The paper also finds
payments increase prescribing to both recommended and contraindicated patients.

Spillover is therefore not a hypothetical caveat to note in a limitations
section. It is a measured, first-order part of the effect, and a design that
ignores it is mis-specified.

## So what, if anything, is left

Three things, all narrower than a research contribution.

**1. Public-data reproducibility.** Carey et al. use prescription data for a
panel of Part D enrollees, which requires a CMS data use agreement. This pipeline
runs entirely on public provider-level aggregates that anyone can re-execute
today with no application and no credentials. That lowers the barrier to
replication and teaching. It also means weaker data: provider-year aggregates
rather than enrollee-month claims.

**2. An estimator-bias benchmark.** The published papers argue their design
handles targeting. This repository *measures* what each design returns under
each selection regime against known ground truth, and reports interval coverage
rather than point estimates alone. That artifact does not exist in this
literature and is useful for calibrating how much to trust a given design. It is
a methods and teaching contribution, not a finding.

**3. A genuinely open empirical question.** Carey et al. end in 2015. Agha and
Zeltzer end in 2016. The data assembled here runs **2019 to 2024** and covers the
GLP-1 era, which is the largest promotional event in diabetes in decades and
falls entirely after both studies. Ozempic, Rybelsus and Mounjaro alone account
for roughly 845,000 of the 998,693 matched payment records in 2023. Whether the
published effect sizes hold in that period, for drugs with unprecedented demand
and persistent supply constraints, is not something either paper can speak to.

That third point is the only part of this project that could produce a result
the literature does not already have.

**Update, after estimation.** It did. On 2019 to 2024 data the effect among
physicians is +4.0% (conservatively +2.5%), in line with Carey et al., so published magnitudes
roughly hold into the GLP-1 era. A fourth point also emerged: on annual public
data, the published physician-by-drug specification fails a falsification test
(next year's payment predicts this year's prescribing, t = 16), and adding
physician-by-year effects removes most of the failure and lowers the estimate by
about 45%. That is a statement about annual data, not a criticism of the
original monthly analysis, whose pre-trends are flat. See
[results.md](results.md).

## How to read this repository

As a reproducible open-data replication and extension to a period the
published work does not cover, with an estimator-bias benchmark. Not as a new
identification strategy. The physician-by-drug effects are Carey, Lieber and
Miller's; adding physician-by-year effects is a modest refinement that annual
data turn out to need.
