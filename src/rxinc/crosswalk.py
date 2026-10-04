"""Map Part D brand strings and Open Payments product names to brand families.

The two sources name products differently, and the difference is not cosmetic.
Open Payments reports what a payment was *about* at the brand level, so a lunch
promoting the Janumet franchise is recorded as ``JANUMET``. Part D splits that
same franchise across formulations and devices: ``Janumet`` and ``Janumet Xr``,
``Victoza 2-Pak`` and ``Victoza 3-Pak``, ``Bydureon Pen`` and ``Bydureon Bcise``.

Joining on raw strings therefore silently drops prescribing. An earlier version
of this pipeline matched on exact uppercase names and lost Victoza entirely,
1.58 million claims in 2019, more than Ozempic that year.

The fix is to join at the **brand family**: the unit a payment can be about.
Part D claims are summed within a family. :func:`crosswalk_report` makes every
mapping, and every name that fails to map, visible for review, because a
crosswalk error does not raise; it just biases the estimate toward zero.
"""

from __future__ import annotations

import re
from collections.abc import Iterable

import pandas as pd

# Tokens that distinguish formulations, pack sizes or devices within a family.
_VARIANT_TOKENS = {"XR", "ER", "PEN", "BCISE", "2-PAK", "3-PAK", "PAK"}
_PACK = re.compile(r"\b\d+-PAK\b")
# Strengths appear as "100-33" in Part D and "100/33" in Open Payments.
_STRENGTH = re.compile(r"\b\d+(?:\.\d+)?(?:[-/]\d+(?:\.\d+)?)?\b")
_PUNCT_ONLY = re.compile(r"^[^A-Z0-9]+$")

FAMILY_CLASS: dict[str, str] = {
    # SGLT2 inhibitors and their fixed-dose combinations
    "JARDIANCE": "SGLT2", "FARXIGA": "SGLT2", "INVOKANA": "SGLT2",
    "STEGLATRO": "SGLT2", "SYNJARDY": "SGLT2", "XIGDUO": "SGLT2",
    "GLYXAMBI": "SGLT2", "INVOKAMET": "SGLT2", "QTERN": "SGLT2",
    "TRIJARDY": "SGLT2",
    "SEGLUROMET": "SGLT2", "STEGLUJAN": "SGLT2",
    # DPP-4 inhibitors and their fixed-dose combinations
    "JANUVIA": "DPP-4", "JANUMET": "DPP-4", "TRADJENTA": "DPP-4",
    "ONGLYZA": "DPP-4", "JENTADUETO": "DPP-4", "KOMBIGLYZE": "DPP-4",
    "NESINA": "DPP-4", "KAZANO": "DPP-4", "OSENI": "DPP-4",
    # GLP-1 receptor agonists, including the dual GIP/GLP-1 agonist
    "OZEMPIC": "GLP-1", "TRULICITY": "GLP-1", "RYBELSUS": "GLP-1",
    "VICTOZA": "GLP-1", "MOUNJARO": "GLP-1", "BYETTA": "GLP-1",
    "BYDUREON": "GLP-1", "ADLYXIN": "GLP-1",
}

FAMILY_MANUFACTURER: dict[str, str] = {
    "JARDIANCE": "BI/Lilly", "SYNJARDY": "BI/Lilly", "GLYXAMBI": "BI/Lilly",
    "TRADJENTA": "BI/Lilly", "JENTADUETO": "BI/Lilly", "TRIJARDY": "BI/Lilly",
    "FARXIGA": "AstraZeneca", "XIGDUO": "AstraZeneca", "QTERN": "AstraZeneca",
    "ONGLYZA": "AstraZeneca", "KOMBIGLYZE": "AstraZeneca",
    "BYDUREON": "AstraZeneca", "BYETTA": "AstraZeneca",
    "INVOKANA": "Janssen", "INVOKAMET": "Janssen",
    "JANUVIA": "Merck", "JANUMET": "Merck", "STEGLATRO": "Merck",
    "SEGLUROMET": "Merck", "STEGLUJAN": "Merck",
    "OZEMPIC": "Novo Nordisk", "RYBELSUS": "Novo Nordisk", "VICTOZA": "Novo Nordisk",
    "TRULICITY": "Lilly", "MOUNJARO": "Lilly",
    "NESINA": "Takeda", "KAZANO": "Takeda", "OSENI": "Takeda",
    "ADLYXIN": "Sanofi",
}

# Families deliberately outside the comparison set, kept so they are reported
# as excluded rather than silently unmatched.
EXCLUDED_FAMILIES: dict[str, str] = {
    "SOLIQUA": "insulin combination",
    "XULTOPHY": "insulin combination",
    "WEGOVY": "obesity indication",
    "SAXENDA": "obesity indication",
    "ZEPBOUND": "obesity indication",
}


def brand_family(name: object) -> str | None:
    """Normalise a raw product or brand string to its brand family.

    Upper-cases, drops dosage strengths such as ``100-33``, and drops variant
    tokens that distinguish formulations, pack sizes or devices. Returns
    ``None`` for missing input.

    Args:
        name: A Part D ``Brnd_Name`` or an Open Payments product name.

    Returns:
        The family key, e.g. ``"VICTOZA"`` for ``"Victoza 3-Pak"``, or ``None``.
    """
    if name is None or (isinstance(name, float) and pd.isna(name)):
        return None
    text = str(name).upper().strip()
    if not text:
        return None
    # Pack sizes go first: stripping strengths first would eat the "2" of
    # "2-PAK" and leave "-PAK" behind, which is exactly how Victoza was lost.
    text = _PACK.sub(" ", text)
    text = _STRENGTH.sub(" ", text)
    tokens = [
        t for t in re.split(r"\s+", text)
        if t and t not in _VARIANT_TOKENS and not _PUNCT_ONLY.match(t)
    ]
    return " ".join(tokens) or None


def in_class(family: str | None) -> bool:
    """Whether a family belongs to the analysed comparison set."""
    return family is not None and family in FAMILY_CLASS


def crosswalk_report(
    part_d_names: Iterable[str], payment_names: Iterable[str]
) -> pd.DataFrame:
    """Tabulate every raw name on both sides against its family.

    Args:
        part_d_names: Distinct Part D ``Brnd_Name`` values.
        payment_names: Distinct Open Payments product names.

    Returns:
        One row per (source, raw name) with ``family``, ``in_class``, and
        ``status``: ``"matched"`` when the family appears on both sides,
        ``"part_d_only"`` or ``"payments_only"`` when it appears on one,
        ``"excluded"`` for deliberate exclusions, and ``"out_of_class"``
        otherwise. One-sided families are the thing to inspect: a family with
        prescribing but no payments may be genuinely unpromoted, or may be a
        naming mismatch that would bias the estimate.
    """
    rows = []
    for source, names in (("part_d", part_d_names), ("payments", payment_names)):
        for raw in sorted({str(n) for n in names if n is not None}):
            rows.append({"source": source, "raw_name": raw, "family": brand_family(raw)})
    report = pd.DataFrame(rows, columns=["source", "raw_name", "family"])
    if report.empty:
        return report.assign(in_class=[], status=[])

    fams_d = set(report.loc[report.source == "part_d", "family"].dropna())
    fams_p = set(report.loc[report.source == "payments", "family"].dropna())

    def status(fam: str | None) -> str:
        if fam in EXCLUDED_FAMILIES:
            return "excluded"
        if not in_class(fam):
            return "out_of_class"
        if fam in fams_d and fam in fams_p:
            return "matched"
        return "part_d_only" if fam in fams_d else "payments_only"

    report["in_class"] = report["family"].map(in_class)
    report["status"] = report["family"].map(status)
    return report
