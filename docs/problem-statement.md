# Problem statement

## Context

Pharmaceutical manufacturers maintain large, continuous financial relationships
with prescribing clinicians. In Open Payments Program Year 2025, manufacturers
reported $14.67 billion across 17.07 million records, the highest annual total
since disclosure began in 2013, covering roughly 1.08 million physicians.

Two audiences need to know what that spend does:

- **Commercial analytics teams** allocate promotional budget across physicians
  and channels, on an assumed return that is rarely measured causally.
- **Regulators, payers and health systems** need to know whether disclosure
  addresses a real behavioral effect or documents an association.

Both work from the same weak evidence base.

## The problem

**There is no reliable, reproducible estimate of how much prescribing
pharmaceutical promotional spend actually causes, and the standard designs used
to produce such estimates fail in ways their own output does not reveal.**

The obstruction is structural rather than statistical noise. Manufacturers
select which physicians to engage on the basis of prescribing behavior, so
exposure is correlated with the outcome by construction. Modern commercial
targeting compounds this: propensity and switching-likelihood models assign
engagement using dynamic cohorts that update as new prescription data arrives,
meaning selection keys on prescribing trajectory, not merely level.

The two cases have different consequences:

- If selection is on **level** (persistently high prescribers), within-physician
  designs difference it away and published estimates are defensible.
- If selection is on **trajectory** (physicians already rising), within-physician
  designs do not fix it, because treated physicians were already diverging
  before exposure.

Which case holds is an empirical question the literature has largely not
separated.

## Why it is unsolved

A 2021 systematic review in *Annals of Internal Medicine* assessed 36 studies of
payments and prescribing and found 21 at serious risk of bias. It named the
mechanism explicitly, that dose-response patterns "may also reflect residual
confounding if industry targets clinicians who already have higher baseline
prescribing volumes," and closed by calling for instrumental variables,
interrupted time series and policy natural experiments.

That call has largely not been answered, for three reasons:

1. **The counterfactual is unobservable.** Randomising promotional exposure is
   neither legal nor commercially plausible, so there is no experimental
   benchmark.
2. **Industry-standard exposure and outcome data are proprietary and measured
   with error.** Commercial prescriber panels project national estimates from
   roughly 93% retail sampling using geospatial assumptions known to be false;
   validation work found per-physician projection errors near -8.4% and -5.6%
   in rural areas. That error is rarely propagated downstream.
3. **Analysis is fragmented and unreproducible.** Commercial work is
   proprietary, academic work is one-off, and there is no shared open tooling
   for the linkage.

## Why it is tractable anyway

The necessary data is public and free of protected health information:

- **Medicare Part D Prescribers**, with actual adjudicated claim counts rather
  than sample projections, which is why validation studies use it as the
  benchmark.
- **Open Payments**, with every reportable manufacturer payment, attributed to
  a specific product.

Both are provider-level aggregates. See [data-sources.md](data-sources.md) for
verified scale, coverage and limits.

## Scope

**In scope:** measuring how much standard estimators misstate the effect under
realistic selection; making that failure detectable from observed data;
reproducible tooling over public CMS data.

**Out of scope:** establishing the true magnitude of the effect; claims about
individual physicians or manufacturers; anything requiring proprietary
prescriber panels or patient-level data.

## What success looks like

1. Quantified, replicated bias for each standard design under each selection
   regime.
2. A diagnostic that reliably separates the identified case from the
   unidentified one.
3. An end-to-end reproducible pipeline over public data.
4. An honest statement of what remains unresolved.

## Key risk

The project can characterise and detect the problem. If the deliverable is read
as "here is the causal effect," it overclaims. The framing must remain: here is
why existing estimates cannot be trusted, and here is how to tell.
