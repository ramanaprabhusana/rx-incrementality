"""The Open Payments stream must never leave a truncated extract that looks complete."""

from __future__ import annotations

import importlib.util
import io
import sys
from pathlib import Path

import pandas as pd
import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "stream_open_payments.py"
PRODUCT = "Name_of_Drug_or_Biological_or_Device_or_Medical_Supply_1"


@pytest.fixture
def stream(monkeypatch):
    spec = importlib.util.spec_from_file_location("stream_open_payments", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    monkeypatch.setattr(module, "open_payments_general_distribution",
                        lambda year: type("D", (), {"download_url": f"https://example.invalid/PGYR{year}.csv"})())
    return module


def _serve(monkeypatch, module, body: bytes, declared_length: int):
    class Response(io.BytesIO):
        headers = {"Content-Length": str(declared_length)}

    monkeypatch.setattr(module.urllib.request, "urlopen", lambda request, timeout=None: Response(body))


def _body() -> bytes:
    rows = pd.DataFrame({
        "Covered_Recipient_NPI": ["1234567893", "1234567893", "1234567893"],
        "Total_Amount_of_Payment_USDollars": ["20", "30", "40"],
        "Nature_of_Payment_or_Transfer_of_Value": ["Food and Beverage"] * 3,
        PRODUCT: ["JARDIANCE", "LIPITOR", "Ozempic"],
    })
    return rows.to_csv(index=False).encode()


def _run(monkeypatch, module, tmp_path) -> int:
    monkeypatch.setattr(sys, "argv", ["stream_open_payments.py", "--year", "2023", "--out-dir", str(tmp_path)])
    return module.main()


def test_complete_download_is_renamed_into_place(stream, monkeypatch, tmp_path):
    body = _body()
    _serve(monkeypatch, stream, body, len(body))
    assert _run(monkeypatch, stream, tmp_path) == 0
    out = tmp_path / "openpayments_diabetes_2023.csv"
    kept = pd.read_csv(out)
    assert sorted(kept[PRODUCT]) == ["JARDIANCE", "Ozempic"]
    assert not out.with_name(out.name + ".part").exists()


def test_truncated_download_fails_and_leaves_no_final_file(stream, monkeypatch, tmp_path, capsys):
    body = _body()
    _serve(monkeypatch, stream, body, len(body) + 1000)   # the server promised more than arrived
    assert _run(monkeypatch, stream, tmp_path) == 1
    assert not (tmp_path / "openpayments_diabetes_2023.csv").exists()
    assert "INCOMPLETE" in capsys.readouterr().err


def test_a_failed_run_does_not_disturb_an_earlier_complete_file(stream, monkeypatch, tmp_path):
    body = _body()
    _serve(monkeypatch, stream, body, len(body))
    assert _run(monkeypatch, stream, tmp_path) == 0
    out = tmp_path / "openpayments_diabetes_2023.csv"
    before = out.read_bytes()
    _serve(monkeypatch, stream, body, len(body) + 1000)
    assert _run(monkeypatch, stream, tmp_path) == 1
    assert out.read_bytes() == before
