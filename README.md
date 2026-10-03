# rx-incrementality

**A reproducible, public-data harness for measuring which causal designs survive
the way pharmaceutical manufacturers actually target physicians, and which
silently do not.**

Not a new identification strategy. The design used here was published by Carey,
Lieber and Miller (2021) on better data. What this adds is a benchmark of
estimator bias against known ground truth, an end-to-end pipeline anyone can run
without a data use agreement, and data covering a period the published studies
predate. See [What is new here, and what is not](#what-is-new-here-and-what-is-not).

## The problem

Manufacturers reported **$14.67 billion** in payments to clinicians in Open
Payments Program Year 2025 across 17.07 million records. More than 85% of drug
firms' marketing expenditure targets physicians, and 29% of Part D physicians
are paid for at least one drug over a typical sample period.

Measuring whether that spend *changes* prescribing is hard for a structural
reason: manufacturers select physicians on the very behaviour being measured.
Exposure is correlated with the outcome by construction.

## Finding 1: standard designs fail, and say nothing about it

Simulate a physician panel where the true effect is known to be exactly **+0.05
log points**. Vary only how manufacturers choose whom to pay.

| Estimator | Random | Targets **high** prescribers | Targets **rising** prescribers |
| --- | ---: | ---: | ---: |
| Naive pooled OLS | 0.092 | 0.606 | 0.612 |
| Two-way fixed effects | 0.048 | **0.050** | 0.098 |
| DiD vs never-treated | 0.048 | **0.048** | -0.007 |
| ITS level change | 0.048 | **0.051** | 0.110 |
| *95% CI coverage* | *88 to 93%* | *95 to 100%* | ***0%*** |

If selection is on **levels**, panel methods work. If selection is on
**momentum**, they break incoherently: fixed effects and interrupted time series
roughly double the effect, differences-in-differences erases it by anchoring on
a transitory pre-onset peak. **Coverage falls to 0% for all four.** Each reports
a tight interval that excludes the truth.

## Finding 2: comparing drugs instead of physicians survives it

Open Payments attributes each payment to a specific product and Part D reports
prescribing by drug, so the comparison can be made **within a physician-year,
across drugs**:

```
log(claims_ijt) = tau * Paid_ijt + alpha_it + delta_jt + e_ijt
```

The physician-by-year effect absorbs that physician's entire trajectory by
construction. Monte Carlo, 40 replications, true effect 0.050:

| Selection regime | Physician-drug + period FE | Triple difference |
| --- | --- | --- |
| Random | 0.049, coverage 97.5% | 0.056, coverage 97.5% |
| **Physician momentum** | **0.128, coverage 0%** | **0.056, coverage 97.5%** |
| Drug-specific momentum | 0.205, coverage 0% | 0.094, coverage 47.5% |

The conventional design overstates by 2.6x with zero coverage; this one is
unbiased with correct coverage.

**The bottom row is a real limitation.** Drug-specific targeting does break the
design, because that confounder varies within the physician-year, which is the
dimension identification relies on. The diagnostic separates the cases:
pre-trends pass under physician momentum and reject at p ~ 1e-88 under drug
momentum. Report the test with the estimate.

## What is new here, and what is not

**Not new.** The identification strategy. [Carey, Lieber and Miller
(2021)](https://doi.org/10.1016/j.jpubeco.2021.104402) use fixed effects for each
physician-drug combination, at **monthly** granularity, on Open Payments 2013 to
2015 linked to Part D enrollee claims. They motivate it with the same targeting
problem and already report the pre-trend result. Their monthly data also
dissolves this project's annual-periodicity limitation. [Agha and Zeltzer
(2022)](https://doi.org/10.1257/pol.20200044) estimate a **23% volume increase**
in anticoagulants over 2014 to 2016, and find **peer spillovers account for about
a quarter of it**.

**Arguably useful.**

1. **Public-data reproducibility.** The published work uses enrollee-level Part D
   data requiring a CMS data use agreement. This runs on public provider-level
   aggregates with no credentials. Lower barrier, weaker data.
2. **An estimator-bias benchmark.** The papers argue their design handles
   targeting. This measures what each design returns under each regime against
   known truth, reporting interval coverage rather than point estimates alone.
   That artifact does not exist in this literature.
3. **An uncovered period.** Carey et al. end in 2015, Agha and Zeltzer in 2016.
   The data here runs **2019 to 2024** and covers the GLP-1 era, the largest
   promotional event in diabetes in decades. Ozempic, Rybelsus and Mounjaro alone
   account for roughly 845,000 of the 998,693 matched 2023 payment records.
   Whether published effect sizes hold there is genuinely open.

Full accounting in [docs/related-work.md](docs/related-work.md).

## Data acquired

Diabetes agents (SGLT2, DPP-4, GLP-1), 23 brands, 2019 to 2024.

| Source | Grain | Scale |
| --- | --- | --- |
| Part D by Provider and Drug | physician-drug-year | 3,491,826 rows, 142,022 to 268,694 prescribers per year |
| Open Payments general | payment record | ~5.47M records matched across 6 years |

Both public, key-free, and free of protected health information.

Feasibility measured before committing to the design: **47.3%** of Indiana
anticoagulant prescribers use two or more competing drugs; **59.6%** of
physicians are paid about two or more diabetes brands. Payment sizes rule out a
binary treatment, with a median physician-brand-year total of **$33.32** and 97%
of records being food and beverage, so treatment is `log(1 + amount)`.

## Quickstart

```bash
pip install -e ".[dev]"
make demo      # reproduce the tables above
make test      # 70 tests
```

```python
from rxinc import simulate_drug_panel, triple_diff
from rxinc.estimators import drug_event_study
from rxinc.diagnostics import pretrend_test

panel = simulate_drug_panel()                   # physician momentum selection
print(triple_diff(panel))                       # recovers the truth
print(pretrend_test(drug_event_study(panel)))   # and confirms it is identified
```

```bash
rxinc catalog
python3 scripts/fetch_class_panel.py --class diabetes --start 2019 --end 2024
python3 scripts/stream_open_payments.py --year 2023 --class diabetes
```

## Layout

```
src/rxinc/
  simulate.py     single-outcome panel DGP, three selection regimes
  drugpanel.py    physician-drug-year DGP for validating the triple difference
  estimators.py   naive OLS, two-way FE, DiD, event study, ITS, triple diff
  diagnostics.py  pre-trend test, placebo, balance, Monte Carlo
  datasets.py     CMS Part D and Open Payments clients
  linkage.py      NPI-keyed linkage with name fallback
  panel.py        CMS data to estimator-ready panel
scripts/          data acquisition
docs/             problem statement, approach, related work, methodology, data, findings
```

Fixed effects use alternating projections, checked against explicit dummy
regression on balanced and unbalanced panels. Errors cluster on physician.

## Honest limitations

- **No real-data estimate has been produced yet.** Data is acquired and the
  estimator is validated on simulation, but the physician-drug-year panel is not
  built and no CMS number is reported here.
- **Peer spillover is unhandled.** Agha and Zeltzer measured it at roughly a
  quarter of the total effect, which means untreated physicians are not a clean
  control. A design ignoring this is mis-specified, and this one currently does.
- Part D public files are annual, so payments and prescribing within a year
  cannot be ordered. The published work uses monthly data and does not have this
  problem.
- Part D covers Medicare beneficiaries only, which matters especially for GLP-1s
  given substantial commercial and cash-pay use.
- Simulation parameters are set to make each failure visible, not calibrated to
  observed manufacturer behaviour. The *direction* of each failure is a property
  of the design; the *magnitudes* are not estimates of real-world bias.
- The triple difference does not solve drug-specific targeting. It makes that
  case detectable, which is not the same as solving it.

## Sources

- [Carey, Lieber and Miller (2021), Drug Firms' Payments and Physicians' Prescribing Behavior in Medicare Part D, Journal of Public Economics 197](https://doi.org/10.1016/j.jpubeco.2021.104402) ([NBER WP 26751](https://www.nber.org/papers/w26751))
- [Agha and Zeltzer (2022), Drug Diffusion Through Peer Networks, American Economic Journal: Economic Policy 14(2)](https://doi.org/10.1257/pol.20200044)
- [Are Financial Payments From the Pharmaceutical Industry Associated With Physician Prescribing? A Systematic Review, Annals of Internal Medicine 174(3)](https://doi.org/10.7326/M20-5665)
- [Open Payments, CMS](https://www.cms.gov/priorities/key-initiatives/open-payments)
- [Medicare Part D Prescribers by Provider and Drug, CMS](https://data.cms.gov/provider-summary-by-type-of-service/medicare-part-d-prescribers/medicare-part-d-prescribers-by-provider-and-drug)

## License

MIT
