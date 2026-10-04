# rx-incrementality

**Do pharmaceutical payments to physicians still move prescribing in the GLP-1
era, and how much of the apparent effect is really targeting? Estimated on
public CMS data for 2019 to 2024, with every design first validated against
simulated ground truth.**

## Results

A payment relationship about a diabetes drug is associated with **+3.9%** more
Medicare claims for that drug (95% CI 3.6 to 4.2%), and a **+2.7 point** higher
chance of prescribing it at all, among 89,011 established prescribers. The
effect builds from +2.3% in the year of first payment to +6.0% one to two years
later, with flat pre-trends. A conservative reading that controls for
anticipatory targeting gives **+2.6%**.

That is in line with Carey, Lieber and Miller's +2.2% to +5.2% for 2013 to 2015.
**The published effect sizes roughly hold into the GLP-1 era**: +4.5% for GLP-1
agents, +3.7% for SGLT2 inhibitors, and a precise zero for the mature DPP-4
class.

## Most of the apparent effect is targeting

| Fixed effects | Effect on claims | Future payment coefficient (should be zero) |
| --- | ---: | ---: |
| Physician-year + drug-year | +23.3% | +0.122 |
| Physician-drug + drug-year (the published design) | +6.8% | +0.026 |
| **Physician-year + drug-year + physician-drug** | **+3.9%** | **+0.005** |

A payment next year cannot cause prescribing this year, so the right-hand column
is a falsification test. Without physician-drug effects the estimate is six
times too large, because representatives target physicians who already favour
their drug. On annual data, the published design also fails the test, and adding
physician-year effects cuts the estimate by a further 40%.

That does not contradict the published work, whose monthly data show flat
pre-trends; at monthly granularity a lead is weeks, not a year. It does mean
annual public data need all three sets of fixed effects. Full results, event
study and caveats: [docs/results.md](docs/results.md).

## Why the design was chosen: simulated ground truth

Before touching CMS data, each design was run on simulated panels where the true
effect is known to be exactly **0.050**, under four targeting regimes. Thirty
replications each:

| How representatives choose physicians | 2-way | 3-way |
| --- | --- | --- |
| At random | 0.048, coverage 93% | 0.048, coverage 97% |
| Physicians whose prescribing is rising | 0.047, coverage 87% | 0.048, coverage 100% |
| **Physicians who already favour the drug** | **0.686, coverage 0%** | **0.046, coverage 97%** |
| Physicians whose use of that drug is rising | 0.091, coverage 37% | 0.097, coverage 0% |

The 3-way design survives the two realistic forms of targeting and is three times
more precise. **The last row is a real limitation**: targeting on a physician's
rising use of a specific drug breaks every design here. The event-study pre-trend
test is how to detect it, and on the real data it passes.

An earlier single-outcome simulation shows the standard physician-level designs
(fixed effects, differences-in-differences, interrupted time series) collapsing
to 0% coverage under momentum targeting: [docs/findings.md](docs/findings.md).

## What is new here, and what is not

**Not new.** Physician-by-drug fixed effects come from [Carey, Lieber and Miller
(2021)](https://doi.org/10.1016/j.jpubeco.2021.104402), on monthly enrollee data.
[Agha and Zeltzer (2022)](https://doi.org/10.1257/pol.20200044) estimate effects
for anticoagulants and show peer spillovers matter.

**What this adds.**

1. **The GLP-1 era.** Both studies end by 2016. This covers 2019 to 2024.
2. **Evidence that, on annual data, the published design needs physician-year
   effects**, demonstrated with a falsification test rather than asserted.
3. **A bias benchmark**: each design measured against known truth under each
   targeting regime, with interval coverage.
4. **Fully public and reproducible.** No data use agreement. Anyone can re-run
   it.

Full accounting in [docs/related-work.md](docs/related-work.md).

## Data

| Source | Grain | Scale |
| --- | --- | --- |
| Part D by Provider and Drug | physician-drug-year | 3,715,046 rows, 2019 to 2024 |
| Open Payments, general payments | payment record | about 5.49M records matched |

Public, key-free, no protected health information. Joined on NPI and on **brand
family**, because Part D splits a brand across pack sizes and devices (`Victoza
2-Pak`, `Victoza 3-Pak`) while Open Payments reports the brand (`VICTOZA`). An
exact-name join silently lost 1.58M Victoza claims in 2019; the crosswalk now
matches 28 of 29 families on both sides. See
[docs/data-sources.md](docs/data-sources.md).

## Quickstart

```bash
pip install -e ".[dev]"
make test      # 115 tests
make demo      # simulation tables
python3 scripts/fetch_class_panel.py --class diabetes --start 2019 --end 2024
for y in 2019 2020 2021 2022 2023 2024; do python3 scripts/stream_open_payments.py --year $y; done
python3 scripts/estimate_real.py
```

```python
from rxinc import simulate_drug_panel, triple_diff
from rxinc.estimators import drug_event_study
from rxinc.diagnostics import pretrend_test

THREE_WAY = ("physician_year", "drug_year", "physician_drug")
panel = simulate_drug_panel()
print(triple_diff(panel, absorb=THREE_WAY))
print(pretrend_test(drug_event_study(panel, absorb=THREE_WAY)))
```

## Layout

```
src/rxinc/
  simulate.py     single-outcome panel DGP
  drugpanel.py    physician-drug-year DGP, four targeting regimes
  estimators.py   OLS, two-way FE, DiD, event studies, ITS, k-way triple difference
  diagnostics.py  pre-trend test, placebo, balance, Monte Carlo
  crosswalk.py    brand-family mapping between Part D and Open Payments
  drugdata.py     CMS extracts to the estimation panel
  datasets.py     CMS clients, robust to catalog layout changes
scripts/          acquisition and estimation
results/          estimates.json
docs/             results, problem, approach, related work, methods, data, simulation
```

Fixed effects use alternating projections, checked against explicit dummy
regression. Standard errors cluster on physician, with nested fixed effects
excluded from the small-sample correction.

## Limitations

- **Not proof of causation.** Pre-trends pass and the falsification test nearly
  passes, but the residual future-payment coefficient (+0.005, SE 0.0018) shows a
  little anticipatory targeting remains. The honest range is about 2.6 to 3.9%.
- **Medicare only**, which matters for GLP-1 agents given commercial and cash-pay
  use.
- **Established prescribers.** The panel is physicians observed all six years;
  31% of payment dollars go to others.
- **Suppression.** Cells under 11 claims are not reported, so the intensive margin
  is conditional on clearing that floor.
- **Annual timing** cannot order payment and prescribing within a year.
- **Peers** are proxied by city, much cruder than the shared-patient networks in
  the published work.
- **The GLP-1 era split** (+0.7% before 2022, +5.0% after) includes the Mounjaro
  launch, exactly where drug-specific targeting is most plausible.

## Sources

- [Carey, Lieber and Miller (2021), Journal of Public Economics 197](https://doi.org/10.1016/j.jpubeco.2021.104402) ([NBER WP 26751](https://www.nber.org/papers/w26751))
- [Agha and Zeltzer (2022), American Economic Journal: Economic Policy 14(2)](https://doi.org/10.1257/pol.20200044)
- [Annals of Internal Medicine systematic review (2021)](https://doi.org/10.7326/M20-5665)
- [Open Payments, CMS](https://www.cms.gov/priorities/key-initiatives/open-payments)
- [Medicare Part D Prescribers by Provider and Drug, CMS](https://data.cms.gov/provider-summary-by-type-of-service/medicare-part-d-prescribers/medicare-part-d-prescribers-by-provider-and-drug)

## License

MIT
