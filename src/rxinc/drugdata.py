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
        across all products a record names) and ``pay_records`` (records
        naming the family, unsplit).
    """
    use = ["Covered_Recipient_NPI", "Total_Amount_of_Payment_USDollars", *PRODUCT_COLS]
    out = []
    for year, path in paths_by_year.items():
        df = pd.read_csv(path, usecols=lambda c: c in set(use), dtype=str, low_memory=False)
        df = df[df["Covered_Recipient_NPI"].fillna("").str.len() > 0].copy()
        df["amount"] = pd.to_numeric(df["Total_Amount_of_Payment_USDollars"], errors="coerce").fillna(0.0)
        df["record"] = np.arange(len(df))
        slots = [c for c in PRODUCT_COLS if c in df.columns]
        df["n_listed"] = df[slots].notna().sum(axis=1).clip(lower=1)
        long = df.melt(id_vars=["record", "Covered_Recipient_NPI", "amount", "n_listed"],
                       value_vars=slots, value_name="product").dropna(subset=["product"])
        long["family"] = _family_lookup(long["product"])
        long = long[long["family"].map(in_class)]
        long["share"] = long["amount"] / long["n_listed"]
        # The same family named twice on one record counts once.
        long = long.drop_duplicates(["record", "family"])
        agg = long.groupby(["Covered_Recipient_NPI", "family"], as_index=False).agg(
            pay_amount=("share", "sum"), pay_records=("record", "nunique"))
        agg["year"] = int(year)
        out.append(agg)
    pay = pd.concat(out, ignore_index=True).rename(columns={"Covered_Recipient_NPI": "npi"})
    return pay


@dataclass
class DrugPanelReport:
    """What survived construction of the real-data panel."""

    n_physicians: int
    n_families: int
    years: list[int]
    n_rows: int
    n_observed: int
    families: list[str]
    dropped_families: list[str]
    share_paid_cells: float
    n_left_censored_pairs: int
    n_clean_onset_pairs: int
    n_never_paid_pairs: int
    payment_match_rate: float

    def __str__(self) -> str:
        return (
            f"{self.n_physicians:,} physicians x {self.n_families} families x "
            f"{len(self.years)} years = {self.n_rows:,} cells "
            f"({self.n_observed:,} observed at >= {SUPPRESSION_FLOOR} claims). "
            f"Paid cells {self.share_paid_cells:.1%}. Pairs: {self.n_clean_onset_pairs:,} "
            f"clean onset, {self.n_left_censored_pairs:,} left-censored, "
            f"{self.n_never_paid_pairs:,} never paid. "
            f"{self.payment_match_rate:.1%} of in-class payment dollars land on panel physicians."
        )


def build_drug_panel(
    claims: pd.DataFrame,
    attrs: pd.DataFrame,
    payments: pd.DataFrame,
    min_family_peak_claims: float = 100_000,
    require_all_years: bool = True,
    families: list[str] | None = None,
) -> tuple[pd.DataFrame, DrugPanelReport]:
    """Assemble the estimator-ready grid.

    The physician set is fixed across years so that physician-by-year effects
    are compared on a stable population rather than one whose entry and exit
    could itself respond to promotion. Families with small national volume are
    dropped because they add cells that are almost all below the suppression
    floor and carry little information.

    Args:
        claims: From :func:`load_part_d`.
        attrs: From :func:`load_part_d`.
        payments: From :func:`load_payments`.
        min_family_peak_claims: Keep families whose largest annual total, summed
            over observed cells, reaches this.
        require_all_years: Keep only physicians observed in every year.
        families: Explicit family list, overriding the volume threshold. Pass
            this when building on a physician sample, since a threshold applied
            to a sample drops families that clear it easily nationally.

    Returns:
        ``(panel, report)``. Key columns: ``physician_id`` and ``drug_id``
        (integer codes), ``year``, ``family``, ``drug_class``, ``manufacturer``,
        outcomes ``claims``, ``any_rx``, ``log_rx``; exposures ``pay_amount``,
        ``pay_records``, ``pay_any``, ``pay_log``, ``pay_any_lag``,
        ``pay_any_lead``; event-study fields ``first_treat_period``,
        ``event_time``, ``treated``, ``left_censored``; peer exposure
        ``peer_pay_share``; and integer composite keys ``physician_year``,
        ``drug_year`` and ``physician_drug``.
    """
    years = sorted(claims["year"].unique().tolist())
    first_year = years[0]

    peak = claims.groupby(["family", "year"])["claims"].sum().groupby("family").max()
    if families is not None:
        keep_fams = sorted(set(families) & set(peak.index))
    else:
        keep_fams = sorted(peak[peak >= min_family_peak_claims].index)
    dropped = sorted(set(peak.index) - set(keep_fams))
    c = claims[claims["family"].isin(keep_fams)]

    if require_all_years:
        n_years = c.groupby("npi")["year"].nunique()
        npis = sorted(n_years[n_years == len(years)].index)
    else:
        npis = sorted(c["npi"].unique())
    c = c[c["npi"].isin(set(npis))]

    # Family-years that exist at all (Mounjaro has no 2019-2021 cells).
    fam_years = c[["family", "year"]].drop_duplicates()
    phys = pd.DataFrame({"npi": npis})
    grid = phys.merge(fam_years, how="cross")

    grid = grid.merge(c, on=["npi", "family", "year"], how="left")
    grid["any_rx"] = grid["claims"].notna().astype(np.int8)
    grid["log_rx"] = np.log(grid["claims"])
    grid["claims"] = grid["claims"].fillna(0.0)

    pay = payments[payments["family"].isin(keep_fams)]
    total_dollars = float(pay["pay_amount"].sum())
    grid = grid.merge(pay, on=["npi", "family", "year"], how="left")
    landed = float(grid["pay_amount"].sum())
    grid["pay_amount"] = grid["pay_amount"].fillna(0.0)
    grid["pay_records"] = grid["pay_records"].fillna(0).astype(int)
    grid["pay_any"] = (grid["pay_records"] > 0).astype(np.int8)
    grid["pay_log"] = np.log1p(grid["pay_amount"])

    grid = grid.sort_values(["npi", "family", "year"]).reset_index(drop=True)
    by_pair = grid.groupby(["npi", "family"], sort=False)
    grid["pay_any_lag"] = by_pair["pay_any"].shift(1)
    grid["pay_any_lead"] = by_pair["pay_any"].shift(-1)

    paid_years = grid.loc[grid["pay_any"] == 1].groupby(["npi", "family"])["year"].min()
    paid_years = paid_years.rename("first_paid_year")
    grid = grid.merge(paid_years, on=["npi", "family"], how="left")
    grid["left_censored"] = grid["first_paid_year"] == first_year
    clean = grid["first_paid_year"].notna() & ~grid["left_censored"]
    grid["first_treat_period"] = np.where(clean, grid["first_paid_year"], -1).astype(int)
    grid["event_time"] = np.where(clean, grid["year"] - grid["first_paid_year"], np.nan)
    grid["treated"] = (clean & (grid["year"] >= grid["first_paid_year"])).astype(np.int8)
    grid["ever_treated"] = clean

    grid = grid.merge(attrs, on=["npi", "year"], how="left")
    grid["peer_key"] = grid["city"].fillna("").str.upper() + "|" + grid["state"].fillna("")
    g = grid.groupby(["peer_key", "family", "year"])["pay_any"]
    grp_sum, grp_n = g.transform("sum"), g.transform("size")
    # Leave-one-out share of same-city panel physicians paid about this family.
    grid["peer_pay_share"] = np.where(grp_n > 1, (grp_sum - grid["pay_any"]) / (grp_n - 1), 0.0)
    grid["has_peers"] = (grp_n > 1).astype(np.int8)

    grid["family_code"] = grid["family"].map({f: i for i, f in enumerate(keep_fams)})
    year_idx = grid["year"] - first_year
    grid["physician_id"] = grid["npi"].map({n: i for i, n in enumerate(npis)}).astype(np.int64)
    grid["drug_id"] = grid["family_code"].astype(np.int64)
    grid["period"] = year_idx.astype(int)
    grid["physician_year"] = grid["physician_id"] * len(years) + year_idx
    grid["drug_year"] = grid["drug_id"] * len(years) + year_idx
    grid["physician_drug"] = grid["physician_id"] * len(keep_fams) + grid["drug_id"]
    grid["drug_class"] = grid["family"].map(FAMILY_CLASS)
    grid["manufacturer"] = grid["family"].map(FAMILY_MANUFACTURER)
    grid.attrs["true_effect"] = None

    pairs = grid.drop_duplicates(["npi", "family"])
    report = DrugPanelReport(
        n_physicians=len(npis),
        n_families=len(keep_fams),
        years=years,
        n_rows=len(grid),
        n_observed=int(grid["any_rx"].sum()),
        families=keep_fams,
        dropped_families=dropped,
        share_paid_cells=float(grid["pay_any"].mean()),
        n_left_censored_pairs=int(pairs["left_censored"].sum()),
        n_clean_onset_pairs=int(pairs["ever_treated"].sum()),
        n_never_paid_pairs=int(pairs["first_paid_year"].isna().sum()),
        payment_match_rate=landed / total_dollars if total_dollars else 0.0,
    )
    return grid, report
