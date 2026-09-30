# rx-incrementality

**Standard causal designs cannot measure whether pharmaceutical payments change
prescribing, because manufacturers target physicians on the very behaviour
being measured. This repository shows how badly those designs fail, gives a
design that survives the realistic case, and builds the reproducible pipeline
over public CMS data.**

## Why this matters

Manufacturers reported **$14.67 billion** in payments to clinicians in Open
Payments Program Year 2025, across 17.07 million records covering roughly 1.08
million physicians. It is the highest annual total since the Sunshine Act began
collecting in 2013.

Whether that spend *changes* prescribing, rather than merely tracking it,
decides how it should be read by regulators, by health systems, and by the
commercial analytics teams who plan it. The evidence base is weaker than the
volume of literature suggests: a systematic review in *Annals of Internal
Medicine* found **21 of 36 studies at serious risk of bias**, named the
mechanism, that dose-response patterns "may also reflect residual confounding
if industry targets clinicians who already have higher baseline prescribing
volumes," and called for designs that move past association.

## Finding 1: the standard designs fail, silently

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
**momentum**, which is what dynamic-cohort commercial targeting actually does,
they break incoherently: fixed effects and interrupted time series roughly
double the effect, differences-in-differences erases it by anchoring on a
transitory pre-onset peak. **Coverage falls to 0% for all four.** Each reports a
tight interval that excludes the truth.

## Finding 2: comparing drugs instead of physicians survives it

Open Payments attributes each payment to a specific product, and Part D reports
prescribing by drug. So instead of comparing paid physicians to unpaid
physicians, compare **drugs within the same physician-year**:

```
log(claims_ijt) = tau * Paid_ijt + alpha_it + delta_jt + e_ijt
```

The physician-by-year effect `alpha_it` absorbs that physician's entire
trajectory by construction. Monte Carlo, 40 replications, true effect 0.050:

| Selection regime | Physician-drug + period FE | Triple difference |
| --- | --- | --- |
| Random | 0.049, coverage 97.5% | 0.056, coverage 97.5% |
| **Physician momentum** | **0.128, coverage 0%** | **0.056, coverage 97.5%** |
| Drug-specific momentum | 0.205, coverage 0% | 0.094, coverage 47.5% |

The middle row is the result. The conventional design overstates by 2.6x with
zero coverage; the triple difference is unbiased with correct coverage.

**The bottom row is a real limitation, not a footnote.** Drug-specific targeting
does break the design, because that confounder varies within the physician-year,
which is the dimension identification relies on. The diagnostic separates the
cases: pre-trends pass under physician momentum and reject at p ~ 1e-88 under
drug momentum. Report the test with the estimate.

## Verified feasibility on real data

Both required margins of variation were measured on CMS data before the design
was committed to.

- **47.3%** of Indiana anticoagulant prescribers use two or more competing drugs
  (DY2024, 5,695 physicians).
- **59.6%** of physicians are paid about two or more diabetes brands (2023,
  132,756 physicians).

Payment sizes rule out a binary treatment: the median physician-brand-year total
is **$33.32** and 97% of records are "Food and Beverage". Treatment is specified
as `log(1 + amount)`.

## Data acquired

Diabetes agents (SGLT2, DPP-4, GLP-1), 23 brands, 2019 to 2024.

| Source | Grain | Scale |
| --- | --- | --- |
| Part D by Provider and Drug | physician-drug-year | 3,491,826 rows, 142,022 to 268,694 prescribers per year |
| Open Payments general | payment record | ~5.47M records matched across 6 years |

Both public, key-free, and free of protected health information. See
[docs/data-sources.md](docs/data-sources.md) for endpoints, limits and the
acquisition constraints that shaped the tooling.

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

panel = simulate_drug_panel()            # physician momentum selection
print(triple_diff(panel))                # recovers the truth
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
docs/             problem statement, approach, methodology, data, findings
```

Fixed effects use alternating projections, checked against explicit dummy
regression on balanced and unbalanced panels. Errors cluster on physician.

## Honest limitations

- **No real-data estimate has been produced yet.** Data is acquired and the
  estimator is validated on simulation, but the physician-drug-year panel is not
  yet built and no CMS number is reported here.
- Simulation parameters are set to make each failure visible, not calibrated to
  observed manufacturer behaviour. The *direction* of each failure is a property
  of the design; the *magnitudes* are not estimates of real-world bias.
- Part D is annual, so payments and prescribing within a year cannot be ordered.
  This genuinely weakens short-run event studies and cannot be fixed with these
  data.
- The Open Payments API exposes 2019 onward only, so the usable window is six
  annual periods.
- The triple difference does not solve drug-specific targeting. It makes that
  case detectable, which is not the same as solving it.

## Sources

- [Open Payments, CMS](https://www.cms.gov/priorities/key-initiatives/open-payments)
- [Medicare Part D Prescribers by Provider and Drug, CMS](https://data.cms.gov/provider-summary-by-type-of-service/medicare-part-d-prescribers/medicare-part-d-prescribers-by-provider-and-drug)
- [Are Financial Payments From the Pharmaceutical Industry Associated With Physician Prescribing? A Systematic Review, *Annals of Internal Medicine*](https://doi.org/10.7326/M20-5665)
- [Association between industry payments and prescribing costly medications, *BMC Health Services Research*](https://link.springer.com/article/10.1186/s12913-018-3043-8)
- [Comparison of antibiotic prescriptions in IQVIA Xponent and Medicare Part D, PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC9972535/)

## License

MIT
