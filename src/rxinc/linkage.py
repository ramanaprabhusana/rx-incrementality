"""Link Open Payments records to Part D prescribers.

Current Open Payments program years carry ``Covered_Recipient_NPI``, so the
join is a direct key match on NPI and the hard entity-resolution problem
mostly disappears. Two caveats keep a fallback necessary:

* Older program years predate the NPI field entirely, and much published work
  on this linkage was done by matching names and addresses.
* Even in recent years the field can be blank for some covered recipients.

So :func:`link_payments_to_prescribers` joins on NPI where it is present and
usable, falls back to blocking on normalised surname and state, and reports
exactly how many records were matched by each route. The report matters more
than the match: a linkage whose fallback rate is high should not be treated as
interchangeable with one that keyed cleanly on NPI.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

import numpy as np
import pandas as pd

_SUFFIXES = {
    "JR", "SR", "II", "III", "IV", "V",
    "MD", "DO", "PHD", "DDS", "DMD", "NP", "PA", "RN", "DPM", "OD",
}
_NON_ALPHA = re.compile(r"[^A-Z]")
_NPI_PREFIX = "80840"


def normalize_name(value: object) -> str:
    """Upper-case a name and strip punctuation and credential suffixes.

    Args:
        value: Raw name, possibly missing.

    Returns:
        Normalised name, or ``""`` if the input was missing or had no letters.
    """
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return ""
    text = str(value).upper().strip()
    parts = [p for p in re.split(r"[\s,.]+", text) if p]
    kept = [_NON_ALPHA.sub("", p) for p in parts]
    kept = [p for p in kept if p and p not in _SUFFIXES]
    return " ".join(kept)


def is_valid_npi(value: object) -> bool:
    """Check a National Provider Identifier against its Luhn check digit.

    An NPI is ten digits whose final digit is a Luhn check computed over the
    constant prefix ``80840`` followed by the first nine digits. Validating
    catches transcription errors and placeholder values such as ``0000000000``
    that would otherwise produce confident, wrong joins.

    Args:
        value: Candidate NPI.

    Returns:
        ``True`` if the value is a structurally valid NPI.
    """
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return False
    digits = re.sub(r"\D", "", str(value))
    if len(digits) != 10:
        return False

    payload = _NPI_PREFIX + digits[:9]
    total = 0
    # Luhn: double every second digit counting from the right of the payload.
    for index, char in enumerate(reversed(payload)):
        digit = int(char)
        if index % 2 == 0:
            digit *= 2
            if digit > 9:
                digit -= 9
        total += digit
    return (10 - total % 10) % 10 == int(digits[9])


@dataclass
class LinkageReport:
    """How a linkage went, broken down by route.

    Attributes:
        n_payment_rows: Payment records considered.
        n_prescribers: Prescriber records available to match against.
        matched_by_npi: Payment records joined on a valid NPI.
        matched_by_name: Payment records joined by the surname/state fallback.
        unmatched: Payment records left unjoined.
        invalid_npi: Records whose NPI failed the check-digit test.
        ambiguous_name: Records whose fallback block held more than one
            candidate prescriber and so were left unmatched.
    """

    n_payment_rows: int
    n_prescribers: int
    matched_by_npi: int
    matched_by_name: int
    unmatched: int
    invalid_npi: int
    ambiguous_name: int

    @property
    def match_rate(self) -> float:
        """Share of payment records linked by any route."""
        if self.n_payment_rows == 0:
            return 0.0
        return (self.matched_by_npi + self.matched_by_name) / self.n_payment_rows

    @property
    def fallback_share(self) -> float:
        """Share of matches that relied on name blocking rather than NPI."""
        total = self.matched_by_npi + self.matched_by_name
        return self.matched_by_name / total if total else 0.0

    def __str__(self) -> str:
        return (
            f"Linked {self.match_rate:.1%} of {self.n_payment_rows:,} payment rows "
            f"({self.matched_by_npi:,} by NPI, {self.matched_by_name:,} by name; "
            f"{self.unmatched:,} unmatched, {self.invalid_npi:,} invalid NPI, "
            f"{self.ambiguous_name:,} ambiguous). "
            f"Fallback share {self.fallback_share:.1%}."
        )


def link_payments_to_prescribers(
    payments: pd.DataFrame,
    prescribers: pd.DataFrame,
    payment_npi_col: str = "Covered_Recipient_NPI",
    payment_last_col: str = "Covered_Recipient_Last_Name",
    payment_first_col: str = "Covered_Recipient_First_Name",
    payment_state_col: str = "Recipient_State",
    prescriber_npi_col: str = "Prscrbr_NPI",
    prescriber_last_col: str = "Prscrbr_Last_Org_Name",
    prescriber_first_col: str = "Prscrbr_First_Name",
    prescriber_state_col: str = "Prscrbr_State_Abrvtn",
) -> tuple[pd.DataFrame, LinkageReport]:
    """Attach a prescriber NPI to each payment record.

    NPI is used where present and check-digit valid. Remaining records are
    blocked on normalised surname plus state; a block is accepted only when it
    resolves to exactly one prescriber whose normalised first name matches, or
    shares a first initial when no exact first-name match exists. Blocks with
    several surviving candidates are deliberately left unmatched rather than
    resolved arbitrarily, because a wrong link is worse than a missing one --
    it fabricates a treatment assignment.

    Args:
        payments: Open Payments records.
        prescribers: Part D prescriber records.
        payment_npi_col: NPI column in ``payments``.
        payment_last_col: Surname column in ``payments``.
        payment_first_col: First-name column in ``payments``.
        payment_state_col: State column in ``payments``.
        prescriber_npi_col: NPI column in ``prescribers``.
        prescriber_last_col: Surname column in ``prescribers``.
        prescriber_first_col: First-name column in ``prescribers``.
        prescriber_state_col: State column in ``prescribers``.

    Returns:
        ``(linked, report)`` where ``linked`` is ``payments`` plus
        ``matched_npi`` (nullable) and ``match_method`` (``"npi"``,
        ``"name"`` or ``"unmatched"``).
    """
    work = payments.copy()
    n_rows = len(work)

    known_npis = set(
        str(v) for v in prescribers[prescriber_npi_col].astype(str) if is_valid_npi(v)
    )

    raw_npi = work[payment_npi_col].astype(str).str.replace(r"\D", "", regex=True)
    valid = raw_npi.map(is_valid_npi)
    in_part_d = raw_npi.isin(known_npis)

    matched_npi = pd.Series(pd.NA, index=work.index, dtype="object")
    method = pd.Series("unmatched", index=work.index, dtype="object")

    npi_hit = valid & in_part_d
    matched_npi[npi_hit] = raw_npi[npi_hit]
    method[npi_hit] = "npi"

    has_npi_value = raw_npi.str.len() > 0
    invalid_npi = int((has_npi_value & ~valid).sum())

    # Fallback: block on normalised surname + state.
    need_fallback = ~npi_hit
    ambiguous = 0
    if need_fallback.any():
        ref = prescribers.assign(
            _last=prescribers[prescriber_last_col].map(normalize_name),
            _first=prescribers[prescriber_first_col].map(normalize_name),
            _state=prescribers[prescriber_state_col].astype(str).str.upper().str.strip(),
            _npi=prescribers[prescriber_npi_col].astype(str),
        )
        ref = ref[ref["_npi"].map(is_valid_npi)]
        blocks: dict[tuple[str, str], list[tuple[str, str]]] = {}
        for last, first, state, npi in zip(
            ref["_last"], ref["_first"], ref["_state"], ref["_npi"]
        ):
            if last and state:
                blocks.setdefault((last, state), []).append((first, npi))

        sub = work.loc[need_fallback]
        cand_last = sub[payment_last_col].map(normalize_name)
        cand_first = sub[payment_first_col].map(normalize_name)
        cand_state = sub[payment_state_col].astype(str).str.upper().str.strip()

        for idx, last, first, state in zip(
            sub.index, cand_last, cand_first, cand_state
        ):
            candidates = blocks.get((last, state))
            if not candidates:
                continue
            exact = [npi for cand_first_name, npi in candidates if cand_first_name == first]
            if len(exact) == 1:
                matched_npi[idx] = exact[0]
                method[idx] = "name"
                continue
            if len(exact) > 1:
                ambiguous += 1
                continue
            if first:
                initial = [
                    npi
                    for cand_first_name, npi in candidates
                    if cand_first_name[:1] == first[:1]
                ]
                if len(initial) == 1:
                    matched_npi[idx] = initial[0]
                    method[idx] = "name"
                elif len(initial) > 1:
                    ambiguous += 1

    work["matched_npi"] = matched_npi
    work["match_method"] = method

    report = LinkageReport(
        n_payment_rows=n_rows,
        n_prescribers=len(prescribers),
        matched_by_npi=int((method == "npi").sum()),
        matched_by_name=int((method == "name").sum()),
        unmatched=int((method == "unmatched").sum()),
        invalid_npi=invalid_npi,
        ambiguous_name=ambiguous,
    )
    return work, report
