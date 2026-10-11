"""The README's Python example must run exactly as printed."""

from __future__ import annotations

import re
from pathlib import Path

README = Path(__file__).resolve().parents[1] / "README.md"


def test_readme_python_blocks_run(capsys):
    blocks = re.findall(r"```python\n(.*?)```", README.read_text(), flags=re.S)
    assert blocks, "the README no longer has a Python example"
    for block in blocks:
        exec(compile(block, str(README), "exec"), {})
    printed = capsys.readouterr().out
    assert "Triple diff" in printed and "Pre-trend test" in printed
