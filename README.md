# rx-incrementality

**When industry targets physicians whose prescribing is already rising, every
standard causal design gets the answer wrong — and none of them says so.**

Pharmaceutical manufacturers reported **$14.67 billion** in payments to
clinicians in Open Payments Program Year 2025 across 17.07 million records, the
highest annual total since the Sunshine Act began collecting in 2013. Whether
those payments *change* prescribing, rather than merely tracking it, is the
question that decides how the spend should be read — by regulators, by health
systems, and by the commercial analytics teams who plan it.

The evidence base is weaker than the volume of literature suggests. A
systematic review in *Annals of Internal Medicine* found **21 of 36 studies at
serious risk of bias**, and named the mechanism: dose-response patterns "may
also reflect residual confounding if industry targets clinicians who already
have higher baseline prescribing volumes." It closed by calling for designs —
instrumental variables, interrupted time series, policy natural experiments —
that move past association.

This repository takes that concern seriously enough to measure it.

## The result

Simulate a physician panel where the true effect of a payment relationship is
known to be exactly **+0.05 log points** (~5% more prescribing). Vary only how
manufacturers choose whom to pay. Run the standard estimators.

| Estimator | Random assignment | Targets **high** prescribers | Targets **rising** prescribers |
| --- | ---: | ---: | ---: |
| Naive pooled OLS | 0.092 | 0.606 | 0.612 |
| Two-way fixed effects | 0.048 | **0.050** | 0.098 |
| DiD vs never-treated | 0.048 | **0.048** | −0.007 |
| ITS level change | 0.048 | **0.051** | 0.110 |
| *95% CI coverage* | *88–93%* | *95–100%* | ***0%*** |

*True effect = 0.050. Mean over 40 replications, 1,200 physicians × 16 quarters.*

Read the last two columns together. When targeting keys on **levels**, panel
methods work — they recover 0.050 and their intervals cover the truth 95–100%
of the time. When targeting keys on **momentum**, they break, and they break
*incoherently*:

- Two-way fixed effects **doubles** the effect (0.098).
- Interrupted time series does the same (0.110).
- DiD against never-treated **erases** it (−0.007), because it anchors on the
  period before onset — a transitory peak — and mean reversion cancels the real
  effect. An Ashenfelter dip, inverted.

**Coverage falls to 0% for all four.** Each reports a tight interval that
excludes the truth. A study using any one of them would look precise and be
wrong, and a literature using several would disagree without converging.

## The part that is actionable

The failure is invisible in the estimates and obvious in the diagnostics.

| Regime | Pre-trend test | Pre-period growth gap |
| --- | --- | ---: |
| Random | PASS (p = 0.13) | ~0.000 |
| Targets high prescribers | PASS (p = 0.17) | +0.004 |
| Targets rising prescribers | **FAIL (p = 4×10⁻¹⁰⁷)** | **+0.022** |

A joint Wald test on pre-onset event-study coefficients separates the
identified case from the unidentified one decisively. The balance table shows
the mechanism: physicians who later receive payments were already growing at
0.032 log points per quarter beforehand, against 0.010 for those who never do.

So the recommendation is narrow and cheap: **publish the pre-trend test and a
placebo alongside any estimate.** They run on the same data, cost nothing, and
distinguish the regime where your number means something from the regime where
it does not.

## Quickstart

```bash
pip install -e ".[dev]"
make demo      # reproduce the tables above
make test      # 53 tests
```

```bash
rxinc simulate --targeting dynamic --out panel.csv
rxinc estimate --panel panel.csv --true-effect 0.05
```

```python
from rxinc import simulate_panel, twoway_fe
from rxinc.diagnostics import pretrend_test
from rxinc.estimators import event_study

panel = simulate_panel()
print(twoway_fe(panel))                    # 0.0984 — nearly double the truth
print(pretrend_test(event_study(panel)))   # FAIL — and it tells you why
```

## Real data

The same estimators run unchanged on CMS data, because
`rxinc.panel.build_panel` emits the identical schema to the simulator.

```bash
rxinc catalog                                  # 25 Part D years, 2016–2025 payments
rxinc fetch-partd --state IN --max-rows 5000
```

Both sources are public, key-free, and contain **no protected health
information** — they are provider-level aggregates.

- **Medicare Part D Prescribers** — actual adjudicated claim counts, ~1.42M
  prescribers per year. Not projections from a pharmacy sample, which is why
  validation studies use Part D as the benchmark.
- **Open Payments** — every reportable manufacturer payment, bulk CSV per year.

**A correction worth stating.** Secondary sources commonly assert that Open
Payments carries no NPI and that linkage requires name-and-address matching.
That is outdated: the current detailed general-payments schema includes
`Covered_Recipient_NPI`, verified against the published 91-field data
dictionary. `rxinc.linkage` keys on NPI — validated against its Luhn check
digit — and keeps a surname-plus-state fallback for older years, reporting how
many records took each route. Ambiguous blocks are left unmatched rather than
resolved arbitrarily, because a wrong link fabricates a treatment assignment.

## Layout

```
src/rxinc/
  simulate.py     panel DGP with three targeting regimes
  estimators.py   naive OLS, two-way FE, DiD, event study, ITS
  diagnostics.py  pre-trend test, placebo, balance, Monte Carlo
  datasets.py     CMS Part D and Open Payments clients
  linkage.py      NPI-keyed linkage with name fallback
  panel.py        CMS data to estimator-ready panel
docs/
  methodology.md  estimands, identification, inference
  data-sources.md what the data is and is not
  findings.md     full results and limitations
```

Fixed effects use alternating projections rather than explicit dummies, and the
transform is checked against a dummy-variable regression on both balanced and
unbalanced panels. Standard errors cluster on physician throughout.

## Honest limitations

- Dynamic-targeting parameters are set to make the problem visible, not
  calibrated to observed manufacturer behaviour. The *direction* of each
  failure is a property of the design; the *magnitudes* are not estimates of
  real-world bias.
- The CMS pipeline is verified against live endpoints and unit-tested on
  synthetic frames, but a full Open Payments × Part D panel has not yet been
  built and analysed. That is the next step.
- **No estimator here solves the dynamic case.** Instrumental variables and
  policy natural experiments are not implemented. This repository establishes
  that the problem is real and detectable; it does not claim to fix it.

## Sources

- [Open Payments — CMS](https://www.cms.gov/priorities/key-initiatives/open-payments)
- [Medicare Part D Prescribers by Provider — CMS](https://data.cms.gov/provider-summary-by-type-of-service/medicare-part-d-prescribers/medicare-part-d-prescribers-by-provider)
- [Are Financial Payments From the Pharmaceutical Industry Associated With Physician Prescribing? A Systematic Review — *Annals of Internal Medicine*](https://doi.org/10.7326/M20-5665)
- [Association between industry payments and prescribing costly medications — *BMC Health Services Research*](https://link.springer.com/article/10.1186/s12913-018-3043-8)
- [Comparison of antibiotic prescriptions in IQVIA Xponent and Medicare Part D — PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC9972535/)

## License

MIT
