"""Build a physician-by-brand-family-by-year panel from the CMS extracts.

The output feeds :func:`rxinc.estimators.triple_diff` and
:func:`rxinc.estimators.drug_event_study` directly. Three constraints of the
public data shape every choice here, and each is stated rather than hidden:

**Suppression.** Part D omits any physician-drug-year with fewer than 11
claims. A missing cell means "under 11", not zero. So the panel carries two
outcomes: ``any_rx``, whether the physician prescribed the family at the
reporting threshold, defined on the full grid; and ``log_rx``, log claims,
defined only where observed. Neither alone is the full effect.

**Left-censored treatment.** Open Payments here starts in 2019. A pair paid in
2019 may also have been paid in 2018, so its true onset is unknown. Such pairs
are flagged ``left_censored`` and excluded from event studies, which need a
known onset. Pooled contemporaneous estimates can still use them.

**Coverage.** Open Payments did not report payments to nurse practitioners,
physician assistants and other non-physician practitioners until program year
2021, under the SUPPORT Act; the extracts contain zero such records for 2019 and
2020 and about 360,000 a year after. Their earlier exposure is unknown, not
zero. Pharmacists and physicians in training are never covered recipients, so
their exposure is structurally zero and carries no information. The builder
therefore classifies every prescriber, marks exposure missing before coverage
begins, and moves the left-censoring year to the start of coverage.

**Attribution.** A payment record can name up to five products. Its amount is
split equally across every product it names, including out-of-class ones, so a
lunch about two drugs is not counted twice. Record counts are kept unsplit as an
alternative intensity measure.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from rxinc.crosswalk import FAMILY_CLASS, FAMILY_MANUFACTURER, brand_family, in_class

# Open Payments covered-recipient status by Part D prescriber type. Physicians
# (including dentists, podiatrists and optometrists) are covered throughout the
# window; these non-physician practitioners only from program year 2021.
NPP_TYPES = {
    "Nurse Practitioner", "Physician Assistant", "Certified Clinical Nurse Specialist",
    "Clinical Nurse Specialist", "Certified Nurse Midwife", "Midwife",
    "Certified Registered Nurse Anesthetist (CRNA)", "Anesthesiologist Assistant",
}
_NOT_COVERED_KEYS = (
    "pharmacist", "registered nurse", "licensed practical", "student",
    "psycholog", "social worker", "counselor", "case manager",
    "specialist/technologist", "dietitian", "nutrition",
)
COVERAGE_START = {"physician": 2013, "npp": 2021}


def coverage_class(prescriber_type: object) -> str:
    """Open Payments coverage of a Part D prescriber type.

    Returns ``"physician"`` (covered throughout), ``"npp"`` (non-physician
    practitioner, covered from 2021) or ``"not_covered"`` (never a covered
    recipient, so payment exposure is structurally zero).
    """
    t = "" if prescriber_type is None or (isinstance(prescriber_type, float)
                                         and pd.isna(prescriber_type)) else str(prescriber_type)
    if t in NPP_TYPES:
        return "npp"
    if any(k in t.lower() for k in _NOT_COVERED_KEYS):
        return "not_covered"
    return "physician"


def specialty_group(prescriber_type: object) -> str:
    """Collapse Part D prescriber types into groups large enough to compare."""
    t = str(prescriber_type or "").lower()
    if t in ("family practice", "general practice", "family medicine"):
        return "Family or general practice"
    if t == "internal medicine":
        return "Internal medicine"
    if t == "endocrinology":
        return "Endocrinology"
    if "cardiology" in t or "heart failure" in t:
        return "Cardiology"
    if t == "nephrology":
        return "Nephrology"
    if t == "nurse practitioner":
        return "Nurse practitioner"
    if t == "physician assistant":
        return "Physician assistant"
    return "Other"


PRODUCT_COLS = [
    f"Name_of_Drug_or_Biological_or_Device_or_Medical_Supply_{i}" for i in range(1, 6)
]
SUPPRESSION_FLOOR = 11


def _family_lookup(values: pd.Series) -> pd.Series:
    """Vectorised brand_family over a Series, computed once per distinct value."""
    uniq = values.dropna().unique()
    mapping = {v: brand_family(v) for v in uniq}
    return values.map(mapping)


def load_part_d(paths: Iterable[str | Path]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read Part D class extracts and collapse to brand family.

    Args:
        paths: CSVs written by ``scripts/fetch_class_panel.py``.

    Returns:
        ``(claims, attrs)``. ``claims`` has one row per (npi, family, year)
        with summed ``claims``. ``attrs`` has one row per (npi, year) with
        ``specialty``, ``city`` and ``state`` for peer-group construction.
    """
    use = ["Prscrbr_NPI", "Brnd_Name", "Tot_Clms", "Prscrbr_Type",
           "Prscrbr_City", "Prscrbr_State_Abrvtn", "data_year"]
    frames = [
        pd.read_csv(p, usecols=lambda c: c in set(use), dtype={"Prscrbr_NPI": str},
                    low_memory=False)
        for p in paths
    ]
    raw = pd.concat(frames, ignore_index=True)
    raw["family"] = _family_lookup(raw["Brnd_Name"])
    raw = raw[raw["family"].map(in_class)]
    raw = raw.rename(columns={"Prscrbr_NPI": "npi", "data_year": "year"})

    claims = (
        raw.groupby(["npi", "family", "year"], as_index=False)["Tot_Clms"].sum()
        .rename(columns={"Tot_Clms": "claims"})
    )
    attrs = (
        raw.sort_values("Tot_Clms", ascending=False)
        .groupby(["npi", "year"], as_index=False)
        .agg(specialty=("Prscrbr_Type", "first"),
             city=("Prscrbr_City", "first"),
             state=("Prscrbr_State_Abrvtn", "first"))
    )
    return claims, attrs


def load_payments(paths_by_year: Mapping[int, str | Path]) -> pd.DataFrame:
    """Read Open Payments class extracts and attribute amounts to families.

    Args:
        paths_by_year: Program year to CSV written by
            ``scripts/stream_open_payments.py``.

    Returns:
        One row per (npi, family, year) with ``pay_amount`` (split equally
        across all products a record names), ``pay_records`` (records naming
        the family, unsplit), and the same two restricted to payments other
        than food and beverage, which are mostly speaker fees: under 2% of
        records but about half of all dollars.
    """
    use = ["Covered_Recipient_NPI", "Total_Amount_of_Payment_USDollars",
           "Nature_of_Payment_or_Transfer_of_Value", *PRODUCT_COLS]
    out = []
    for year, path in paths_by_year.items():
        df = pd.read_csv(path, usecols=lambda c: c in set(use), dtype=str, low_memory=False)
        if "Nature_of_Payment_or_Transfer_of_Value" not in df.columns:
            raise ValueError(
                f"{path}: missing Nature_of_Payment_or_Transfer_of_Value, needed to "
                "separate food from speaker and consulting payments"
            )
        df = df[df["Covered_Recipient_NPI"].fillna("").str.len() > 0].copy()
        df["amount"] = pd.to_numeric(df["Total_Amount_of_Payment_USDollars"], errors="coerce").fillna(0.0)
        df["record"] = np.arange(len(df))
        df["food"] = df["Nature_of_Payment_or_Transfer_of_Value"].fillna("").eq("Food and Beverage")
        slots = [c for c in PRODUCT_COLS if c in df.columns]
        df["n_listed"] = df[slots].notna().sum(axis=1).clip(lower=1)
        long = df.melt(id_vars=["record", "Covered_Recipient_NPI", "amount", "n_listed", "food"],
                       value_vars=slots, value_name="product").dropna(subset=["product"])
        long["family"] = _family_lookup(long["product"])
        long = long[long["family"].map(in_class)]
        long["share"] = long["amount"] / long["n_listed"]
        # The same family named twice on one record counts once.
        long = long.drop_duplicates(["record", "family"])
        long["share_nonfood"] = np.where(long["food"], 0.0, long["share"])
        long["record_nonfood"] = np.where(long["food"], np.nan, long["record"])
        agg = long.groupby(["Covered_Recipient_NPI", "family"], as_index=False).agg(
            pay_amount=("share", "sum"), pay_records=("record", "nunique"),
            pay_amount_nonfood=("share_nonfood", "sum"),
            pay_records_nonfood=("record_nonfood", "nunique"))
        agg["year"] = int(year)
        out.append(agg)
    pay = pd.concat(out, ignore_index=True).rename(columns={"Covered_Recipient_NPI": "npi"})
    return pay


@dataclass
class DrugPanelReport:
    """What survived construction of the real-data panel."""

    population: str
    n_physicians: int
    n_families: int
    years: list[int]
    n_rows: int
    n_rows_exposure_known: int
    n_observed: int
    families: list[str]
    dropped_families: list[str]
    share_paid_cells: float
    n_left_censored_pairs: int
    n_clean_onset_pairs: int
    n_never_paid_pairs: int
    payment_match_rate: float
    excluded_by_coverage: dict[str, int]

    def __str__(self) -> str:
        excl = ", ".join(f"{k} {v:,}" for k, v in self.excluded_by_coverage.items()) or "none"
        return (
            f"[{self.population}] {self.n_physicians:,} prescribers x {self.n_families} families x "
            f"{len(self.years)} years = {self.n_rows:,} cells, {self.n_rows_exposure_known:,} with "
            f"payment exposure observed, {self.n_observed:,} at >= {SUPPRESSION_FLOOR} claims. "
            f"Paid cells {self.share_paid_cells:.1%}. Pairs: {self.n_clean_onset_pairs:,} clean onset, "
            f"{self.n_left_censored_pairs:,} left-censored, {self.n_never_paid_pairs:,} never paid. "
            f"{self.payment_match_rate:.1%} of in-class payment dollars land on the panel. "
            f"Excluded by coverage: {excl}."
        )


POPULATIONS = {"physician": {"physician"}, "npp": {"npp"}, "covered": {"physician", "npp"}}
_PAY_COLS = ["pay_amount", "pay_records", "pay_amount_nonfood", "pay_records_nonfood"]


def build_drug_panel(
    claims: pd.DataFrame,
    attrs: pd.DataFrame,
    payments: pd.DataFrame,
    min_family_peak_claims: float = 100_000,
    require_all_years: bool = True,
    families: list[str] | None = None,
    population: str = "physician",
    full_grid: bool = True,
) -> tuple[pd.DataFrame, DrugPanelReport]:
    """Assemble the estimator-ready grid.

    The prescriber set is fixed across years so that physician-by-year effects
    are compared on a stable population. Families with small national volume
    are dropped because their cells are almost all below the suppression floor.

    Args:
        claims: From :func:`load_part_d`.
        attrs: From :func:`load_part_d`.
        payments: From :func:`load_payments`.
        min_family_peak_claims: Keep families whose largest annual total reaches
            this, unless ``families`` is given.
        require_all_years: Keep only prescribers observed in every year.
        families: Explicit family list, overriding the threshold. Pass this when
            building on a sample.
        population: ``"physician"`` (covered throughout, the primary
            population), ``"npp"`` (non-physician practitioners, exposure
            observed from 2021) or ``"covered"`` (both). Never-covered types are
            always excluded.
        full_grid: Build every prescriber x family x year cell, needed for the
            extensive margin. ``False`` keeps observed cells only, which is
            enough for the intensive margin and far smaller.

    Returns:
        ``(panel, report)``. Estimation frames must filter on
        ``exposure_known``: before a prescriber's coverage starts, payment
        columns are missing rather than zero.
    """
    if population not in POPULATIONS:
        raise ValueError(f"population must be one of {sorted(POPULATIONS)}, got {population!r}")
    years = sorted(claims["year"].unique().tolist())
    first_year = years[0]

    peak = claims.groupby(["family", "year"])["claims"].sum().groupby("family").max()
    if families is not None:
        keep_fams = sorted(set(families) & set(peak.index))
    else:
        keep_fams = sorted(peak[peak >= min_family_peak_claims].index)
    dropped = sorted(set(peak.index) - set(keep_fams))
    c = claims[claims["family"].isin(keep_fams)]

    # Classify each prescriber once, by their most common type across years.
    modal = (attrs.dropna(subset=["specialty"])
             .groupby(["npi", "specialty"]).size().rename("n").reset_index()
             .sort_values(["npi", "n"], ascending=[True, False])
             .drop_duplicates("npi").set_index("npi")["specialty"])
    cov = modal.map(coverage_class)
    wanted = POPULATIONS[population]
    in_scope = c["npi"].map(cov)
    excluded = (c.loc[~in_scope.isin(wanted), "npi"].drop_duplicates().map(cov)
                .value_counts().to_dict())
    c = c[in_scope.isin(wanted)]

    if require_all_years:
        n_years = c.groupby("npi")["year"].nunique()
        npis = sorted(n_years[n_years == len(years)].index)
    else:
        npis = sorted(c["npi"].unique())
    c = c[c["npi"].isin(set(npis))]

    if full_grid:
        fam_years = c[["family", "year"]].drop_duplicates()
        grid = pd.DataFrame({"npi": npis}).merge(fam_years, how="cross")
        grid = grid.merge(c, on=["npi", "family", "year"], how="left")
    else:
        grid = c.copy()
    grid["any_rx"] = grid["claims"].notna().astype(np.int8)
    grid["log_rx"] = np.log(grid["claims"])
    grid["claims"] = grid["claims"].fillna(0.0)

    pay = payments[payments["family"].isin(keep_fams)]
    total_dollars = float(pay["pay_amount"].sum())
    grid = grid.merge(pay, on=["npi", "family", "year"], how="left")
    for col in _PAY_COLS:
        if col not in grid:
            grid[col] = 0.0
        grid[col] = grid[col].fillna(0.0)
    landed = float(grid["pay_amount"].sum())

    # Exposure is only observed once the prescriber is a covered recipient.
    grid["coverage"] = grid["npi"].map(cov)
    obs_start = grid["coverage"].map(COVERAGE_START).clip(lower=first_year)
    grid["exposure_start"] = obs_start.astype(int)
    grid["exposure_known"] = grid["year"] >= grid["exposure_start"]
    unknown = ~grid["exposure_known"]
    grid.loc[unknown, _PAY_COLS] = np.nan

    def flag(cond: pd.Series) -> pd.Series:
        return cond.astype(float).where(grid["exposure_known"])

    grid["pay_any"] = flag(grid["pay_records"] > 0)
    grid["pay_log"] = np.log1p(grid["pay_amount"])
    grid["pay_nonfood_any"] = flag(grid["pay_records_nonfood"] > 0)
    grid["pay_food_only"] = flag((grid["pay_records"] > 0) & ~(grid["pay_records_nonfood"] > 0))
    grid["pay_ge25"] = flag(grid["pay_amount"] >= 25)
    grid["pay_ge100"] = flag(grid["pay_amount"] >= 100)

    grid = grid.sort_values(["npi", "family", "year"]).reset_index(drop=True)
    by_pair = grid.groupby(["npi", "family"], sort=False)
    # Above the suppression floor in every year the family exists. Treatment
    # cannot change whether these pairs are observed, so an intensive-margin
    # estimate on them is free of selection through the floor.
    grid["pair_always_observed"] = by_pair["any_rx"].transform("min").astype(bool)
    grid["pay_any_lag"] = by_pair["pay_any"].shift(1)
    grid["pay_any_lead"] = by_pair["pay_any"].shift(-1)

    paid = grid.loc[grid["pay_any"] == 1].groupby(["npi", "family"])["year"].min()
    grid = grid.merge(paid.rename("first_paid_year"), on=["npi", "family"], how="left")
    grid["left_censored"] = grid["first_paid_year"] == grid["exposure_start"]
    clean = grid["first_paid_year"].notna() & ~grid["left_censored"]
    grid["first_treat_period"] = np.where(clean, grid["first_paid_year"], -1).astype(int)
    grid["event_time"] = np.where(clean, grid["year"] - grid["first_paid_year"], np.nan)
    grid["treated"] = (clean & (grid["year"] >= grid["first_paid_year"])).astype(np.int8)
    grid["ever_treated"] = clean

    grid = grid.merge(attrs[["npi", "year", "city", "state"]], on=["npi", "year"], how="left")
    grid["specialty"] = grid["npi"].map(modal)
    grid["specialty_group"] = grid["specialty"].map(specialty_group)
    grid["manufacturer"] = grid["family"].map(FAMILY_MANUFACTURER)
    grid["drug_class"] = grid["family"].map(FAMILY_CLASS)

    own = grid["pay_any"].fillna(0.0)
    own_n = grid["pay_any"].notna().astype(float)
    grid["peer_key"] = grid["city"].fillna("").str.upper() + "|" + grid["state"].fillna("")
    g = grid.groupby(["peer_key", "family", "year"])["pay_any"]
    denom = g.transform("count") - own_n
    grid["peer_pay_share"] = np.where(denom > 0, (g.transform("sum") - own) / denom.where(denom > 0), 0.0)
    grid.loc[unknown, "peer_pay_share"] = np.nan
    grid["has_peers"] = (denom > 0).astype(np.int8)

    m = grid.groupby(["npi", "year", "manufacturer"])["pay_any"].transform("sum")
    grid["pay_same_mfr_other"] = flag((m - own) > 0)

    year_idx = grid["year"] - first_year
    grid["physician_id"] = grid["npi"].map({n: i for i, n in enumerate(npis)}).astype(np.int64)
    grid["drug_id"] = grid["family"].map({f: i for i, f in enumerate(keep_fams)}).astype(np.int64)
    grid["period"] = year_idx.astype(int)
    grid["physician_year"] = grid["physician_id"] * len(years) + year_idx
    grid["drug_year"] = grid["drug_id"] * len(years) + year_idx
    grid["physician_drug"] = grid["physician_id"] * len(keep_fams) + grid["drug_id"]
    grid.attrs["true_effect"] = None

    pairs = grid.drop_duplicates(["npi", "family"])
    known = grid["exposure_known"]
    report = DrugPanelReport(
        population=population,
        n_physicians=len(npis),
        n_families=len(keep_fams),
        years=years,
        n_rows=len(grid),
        n_rows_exposure_known=int(known.sum()),
        n_observed=int(grid["any_rx"].sum()),
        families=keep_fams,
        dropped_families=dropped,
        share_paid_cells=float(grid.loc[known, "pay_any"].mean()),
        n_left_censored_pairs=int(pairs["left_censored"].sum()),
        n_clean_onset_pairs=int(pairs["ever_treated"].sum()),
        n_never_paid_pairs=int(pairs["first_paid_year"].isna().sum()),
        payment_match_rate=landed / total_dollars if total_dollars else 0.0,
        excluded_by_coverage={k: int(v) for k, v in excluded.items()},
    )
    return grid, report
