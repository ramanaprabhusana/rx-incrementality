"""The numbers quoted in the README and results page must match the results files.

Documentation drifts: an estimate is re-run, a table is edited by hand, and a
figure that was right last month is quietly wrong. Each claim below is derived
from ``results/estimates.json`` or ``results/sensitivity.json``, formatted the
way the documents quote it, and looked up in the text. Changing a result
without updating the prose, or the prose without the result, fails here.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
README, RESULTS = "README.md", "docs/results.md"
BOTH = (README, RESULTS)


@pytest.fixture(scope="module")
def est():
    return json.loads((ROOT / "results/estimates.json").read_text())


@pytest.fixture(scope="module")
def sens():
    return json.loads((ROOT / "results/sensitivity.json").read_text())


@pytest.fixture(scope="module")
def docs():
    # Collapse whitespace so a number is found even where the prose wraps mid-phrase.
    return {name: re.sub(r"\s+", " ", (ROOT / name).read_text()) for name in (*BOTH, "docs/data-sources.md")}


def pct(b: float) -> float:
    return 100.0 * math.expm1(b)


def spec(est, label: str, population: str = "physician") -> dict:
    rows = [r for r in est["estimates"] if r["spec"] == label and r["population"] == population]
    assert len(rows) == 1, f"expected exactly one estimate labelled {label!r}, found {len(rows)}"
    return rows[0]


def effect(row: dict) -> str:
    return f"{pct(row['coef']):+.1f}%"


def interval(row: dict) -> str:
    lo, hi = pct(row["ci_low"]), pct(row["ci_high"])
    return f"{lo:.1f} to {hi:.1f}%" if lo > 0 else f"{lo:.1f} to {hi:+.1f}%"


def require(docs, text: str, where) -> None:
    for name in (where,) if isinstance(where, str) else where:
        assert text in docs[name], f"{name} does not quote {text!r}"


def test_sample_sizes(est, docs):
    phys, npp, unbal = (est["panels"][k] for k in ("physician", "npp", "physician_unbalanced"))
    require(docs, f"{phys['n_physicians']:,}", BOTH)
    for value in (phys["n_rows"], phys["n_observed"], phys["n_clean_onset_pairs"], phys["n_never_paid_pairs"],
                  npp["n_physicians"], npp["n_rows_exposure_known"], unbal["n_physicians"]):
        require(docs, f"{value:,}", RESULTS)
    require(docs, f"{100 * phys['share_paid_cells']:.1f}%", RESULTS)


def test_headline_estimates(est, docs):
    main = spec(est, "intensive: log claims ~ paid")
    require(docs, f"**{effect(main)}**", BOTH)
    require(docs, interval(main), RESULTS)
    lo, hi = pct(main["ci_low"]), pct(main["ci_high"])
    assert f"95% CI {lo:.1f} to {hi:.1f}%" in docs[README]

    lead = spec(est, "3-way: all three + next year's payment")
    require(docs, f"**{effect(lead)}**", README)
    require(docs, f"| {effect(lead)} | {interval(lead)} |", RESULTS)

    floor = spec(est, "intensive, pairs above the floor in every year")
    require(docs, f"({effect(floor)})", BOTH)
    require(docs, f"| {effect(floor)} | {interval(floor)} |", RESULTS)

    dollars = spec(est, "intensive: log claims ~ log(1 + $)")
    require(docs, f"| {effect(dollars)} | {interval(dollars)} |", RESULTS)

    adopt = spec(est, "extensive: prescribes at >= 11 claims ~ paid")
    require(docs, f"+{100 * adopt['coef']:.1f} points | {100 * adopt['ci_low']:.1f} to {100 * adopt['ci_high']:.1f} |",
            RESULTS)
    require(docs, f"+{100 * adopt['coef']:.1f} point estimate", README)


def test_specification_comparison(est, docs):
    rows = [("2-way: physician-year + drug-year", "Physician-year + drug-year"),
            ("published: physician-drug + drug-year", "Physician-drug + drug-year (Carey, Lieber and Miller)"),
            ("3-way: all three", "**All three**")]
    for label, name in rows:
        base = spec(est, label)
        extra = spec(est, f"{label} + next year's payment")["extra"]
        lead, se = extra["pay_any_lead"], extra["pay_any_lead_se"]
        bold = "**" if name.startswith("**") else ""
        require(docs, f"| {name} | {bold}{effect(base)}{bold} | {bold}{lead:+.3f} (SE {se:.3f}){bold} |", RESULTS)
        # The figures' alt text quotes both series without signs.
        require(docs, f"{pct(base['coef']):.1f}", BOTH)
        require(docs, f"{pct(lead):.1f}", BOTH)
    naive = spec(est, "2-way: physician-year + drug-year")
    require(docs, effect(naive), README)
    residual = spec(est, "3-way: all three + next year's payment")["extra"]
    require(docs, f"({residual['pay_any_lead']:+.3f}, SE {residual['pay_any_lead_se']:.3f})", BOTH)


def test_event_study_paths(est, docs):
    studies = est["event_studies"]
    claims, pooled = studies["physician"]["cohort (Sun-Abraham)"], studies["physician"]["pooled dummies"]
    adoption = studies["physician_extensive"]["cohort (Sun-Abraham)"]
    pooled_by_period = dict(zip(pooled["rel_periods"], pooled["coefs"], strict=True))
    for rel, c_claims, c_adopt in zip(claims["rel_periods"], claims["coefs"], adoption["coefs"], strict=True):
        middle = f"{pct(pooled_by_period[rel]):+.1f}%" if rel in pooled_by_period else ""
        row = f"| {rel:+d} |" if rel > 0 else f"| {rel} |"
        found = re.search(re.escape(row) + r" \**" + re.escape(f"{pct(c_claims):+.1f}%") + r"\** \| "
                          + re.escape(middle) + r" ?\| " + re.escape(f"{100 * c_adopt:+.2f} pt") + r" \|",
                          docs[RESULTS])
        assert found, f"event-study row for year {rel} does not match the results file"
    require(docs, f"p = {claims['pretrend_p']:.2f}", BOTH)
    require(docs, f"p = {pooled['pretrend_p']:.2f}", RESULTS)
    assert adoption["pretrend_p"] < 0.001
    by_period = dict(zip(claims["rel_periods"], claims["coefs"], strict=True))
    require(docs, f"{pct(by_period[0]):+.1f}% in the year of first payment to {pct(by_period[2]):+.1f}% two years",
            README)
    npp = studies["npp"]["cohort (Sun-Abraham)"]
    require(docs, f"(p = {npp['pretrend_p']:.2f})", RESULTS)


def test_trend_adjusted_estimates(est, docs):
    claims = est["event_studies"]["physician"]["cohort (Sun-Abraham)"]["detrended"]
    adoption = est["event_studies"]["physician_extensive"]["cohort (Sun-Abraham)"]["detrended"]
    # Standard error carried to the percent scale by the delta method.
    se_pct = 100 * math.exp(claims["post_average"]) * claims["post_average_se"]
    require(docs, f"**{pct(claims['post_average']):+.1f}%** (SE {se_pct:.1f})", RESULTS)
    require(docs, f"**{100 * adoption['post_average']:+.2f} pt** (SE {100 * adoption['post_average_se']:.2f})", RESULTS)
    require(docs, f"{100 * adoption['slope']:+.2f} pt (SE {100 * adoption['slope_se']:.2f})", RESULTS)


def test_heterogeneity_and_robustness(est, docs):
    table = {
        "class: GLP-1": "GLP-1 agents", "class: SGLT2": "SGLT2 inhibitors", "class: DPP-4": "DPP-4 inhibitors",
        "GLP-1, 2019-2021": "GLP-1, 2019 to 2021", "GLP-1, 2022-2024": "GLP-1, 2022 to 2024",
        "GLP-1, 2022-2024, excluding Mounjaro": "GLP-1, 2022 to 2024, excluding Mounjaro",
        "specialty: Family or general practice": "Family or general practice",
        "specialty: Internal medicine": "Internal medicine", "specialty: Endocrinology": "Endocrinology",
        "specialty: Cardiology": "Cardiology", "specialty: Nephrology": "Nephrology",
        "paid last year": "Paid last year rather than this year",
        "payment >= $25": "Payment of at least $25", "payment >= $100": "Payment of at least $100",
        "paid + same-city peer share paid": "Controlling for same-city peer exposure",
        "paid + paid about another drug by same manufacturer":
            "Controlling for payments about the same manufacturer's other drugs",
        "all classes excluding Mounjaro (launch)": "Excluding Mounjaro",
        "food only vs any non-food (speaker, consulting, travel)": "Food and beverage only",
    }
    for label, name in table.items():
        row = spec(est, label)
        require(docs, f"| {name} | {effect(row)} | {interval(row)} |", RESULTS)
    for label, name in {"specialty: Nurse practitioner": "Nurse practitioners, 2021 to 2024",
                        "specialty: Physician assistant": "Physician assistants, 2021 to 2024"}.items():
        row = spec(est, label, "npp")
        require(docs, f"| {name} | {effect(row)} | {interval(row)} |", RESULTS)

    for label in ("class: GLP-1", "class: SGLT2", "GLP-1, 2022-2024, excluding Mounjaro"):
        require(docs, effect(spec(est, label)), README)
    npp = spec(est, "non-physician practitioners, 2021-2024: intensive", "npp")
    require(docs, f"({effect(npp)})", README)
    require(docs, f"{effect(npp)} overall", RESULTS)

    food = spec(est, "food only vs any non-food (speaker, consulting, travel)")
    nonfood = f"{pct(food['extra']['pay_nonfood_any']):+.1f}%"
    require(docs, f"({nonfood} against {effect(food)})", README)
    require(docs, f"| Any speaker, consulting or travel payment | {nonfood} |", RESULTS)

    specialties = [pct(r["coef"]) for r in est["estimates"]
                   if r["spec"].startswith("specialty:") and r["population"] == "physician"]
    require(docs, f"({min(specialties):.1f}% to {max(specialties):.1f}%)", README)
    require(docs, f"between +{min(specialties):.1f}% and +{max(specialties):.1f}%", RESULTS)

    other = spec(est, "paid + paid about another drug by same manufacturer")["extra"]
    require(docs, f"by {abs(pct(other['pay_same_mfr_other'])):.1f}% (SE {100 * other['pay_same_mfr_other_se']:.1f}%)",
            RESULTS)
    unbalanced = spec(est, "physicians, unbalanced panel (all years present not required)", "physician_unbalanced")
    require(docs, f"| {effect(unbalanced)} | {interval(unbalanced)} |", RESULTS)


def test_sensitivity_breakdowns(sens, docs):
    studies = sens["studies"]
    claims, adoption, npp = studies["physician"], studies["physician_extensive"], studies["npp"]

    def m(study, target="first_two_years"):
        return f"{100 * study['targets'][target]['breakdown_m']:.2f}"

    require(docs, f"| Breakdown, bends either way | {m(claims)} | {m(adoption)} | {m(npp)} |", RESULTS)
    require(docs, f"M = {m(claims)}", BOTH)
    require(docs, f"M = {m(adoption)}", BOTH)

    names = {"t+0": "Onset year (t0)", "t+1": "Year after (t+1)", "t+2": "Two years after (t+2)",
             "t+3": "Three years after (t+3)", "all_years": "All five years averaged"}
    for key, name in names.items():
        c, a = claims["targets"][key], adoption["targets"][key]
        found = re.search(re.escape(f"| {name} | {pct(c['estimate_m0']):+.1f}% | {m(claims, key)} | "
                                    f"{100 * a['estimate_m0']:+.2f} pt") + r"[^|]*\| " + re.escape(m(adoption, key)),
                          docs[RESULTS])
        assert found, f"breakdown row {name!r} does not match the results file"
    c, a = claims["targets"]["first_two_years"], adoption["targets"]["first_two_years"]
    require(docs, f"**{pct(c['estimate_m0']):+.1f}%** | **{m(claims)}** | **{100 * a['estimate_m0']:+.2f} pt** | "
                  f"**{m(adoption)}**", RESULTS)


def test_one_sided_restrictions(sens, docs):
    studies = sens["studies"]
    claims, adoption, npp = (studies[k]["restricted"] for k in ("physician", "physician_extensive", "npp"))
    assert adoption["concave"]["verdict"] == "supported"
    assert claims["concave"]["verdict"] == "mixed"
    for study in (claims, adoption, npp):
        assert math.isinf(study["concave"]["breakdown_first_two_years"])
    floor = 100 * adoption["concave"]["lowest_lower_limit_on_grid"]
    require(docs, f"+{floor:.2f} points", BOTH)
    require(docs, f"| Lowest lower limit at any M, can only flatten | "
                  f"{pct(claims['concave']['lowest_lower_limit_on_grid']):+.2f}% | **{floor:+.2f} pt** | "
                  f"{pct(npp['concave']['lowest_lower_limit_on_grid']):+.2f}% |", RESULTS)
    for key in ("increasing", "positive_bias"):
        row = " | ".join(f"{100 * s[key]['breakdown_first_two_years']:.2f}" for s in (claims, adoption, npp))
        assert f"| {row} |" in docs[RESULTS], f"{key} breakdowns do not match the results file"


def test_pre_period_bending(sens, docs):
    claims = {b["centred_on"]: b for b in sens["studies"]["physician"]["pre_period_bending"]}
    adoption = {b["centred_on"]: b for b in sens["studies"]["physician_extensive"]["pre_period_bending"]}
    words = {-2: "two", -3: "three", -4: "four"}
    for centre, word in words.items():
        c, a = claims[centre], adoption[centre]
        found = re.search(re.escape(f"| Centred {word} years before onset | {100 * c['value']:+.2f} "
                                    f"(SE {100 * c['se']:.2f}) | ") + r"\**"
                          + re.escape(f"{100 * a['value']:+.2f} (SE {100 * a['se']:.2f})"), docs[RESULTS])
        assert found, f"bending centred {word} years before onset does not match the results file"


def test_unlinkable_records(est, docs):
    quality = est["data_quality"]
    payments = quality["open_payments"]
    assert quality["part_d"]["rows_dropped_invalid_npi"] == 0 and payments["rows_invalid_npi"] == 0
    rows = f"{payments['rows_blank_npi']:,}"
    require(docs, rows, (README, "docs/data-sources.md"))
    require(docs, f"{payments['rows']:,}", (README, "docs/data-sources.md"))
    share_rows = f"{100 * payments['rows_blank_npi'] / payments['rows']:.2f}%"
    share_dollars = f"{100 * payments['dollars_blank_npi'] / payments['dollars']:.2f}%"
    for name in (README, "docs/data-sources.md"):
        assert share_rows in docs[name] and share_dollars in docs[name], f"{name}: unlinkable shares are out of date"


def test_provenance_describes_the_committed_results(est, docs):
    provenance = json.loads((ROOT / "results/provenance.json").read_text())
    for name in ("estimates", "sensitivity"):
        digest = hashlib.sha256((ROOT / f"results/{name}.json").read_bytes()).hexdigest()
        assert provenance[f"{name}_sha256"] == digest, f"results/{name}.json changed after `make provenance` last ran"

    payments = [f for f in provenance["files"] if f["dataset"].startswith("Open Payments")]
    assert len(payments) == 6
    for record in payments:
        assert f"| {record['year']} |" in docs["docs/data-sources.md"]
        assert f"| {int(record['rows']):,} |" in docs["docs/data-sources.md"], f"{record['year']}: matched rows differ"
    assert sum(int(f["rows"]) for f in payments) == est["data_quality"]["open_payments"]["rows"]
    claims = [f for f in provenance["files"] if f["dataset"].startswith("Medicare Part D")]
    require(docs, f"{sum(int(f['rows']) for f in claims):,} rows", README)
