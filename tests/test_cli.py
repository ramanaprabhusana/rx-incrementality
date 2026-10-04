"""Smoke tests: the entry points a first-time user touches must import and parse."""

from __future__ import annotations

import pytest

from rxinc.cli import build_parser


def test_demo_lives_in_the_package():
    """The demo used to import from examples/, which an installed package lacks."""
    from rxinc.demo import main  # noqa: F401


@pytest.mark.parametrize(
    "argv",
    [["demo"], ["simulate", "--targeting", "static"], ["catalog"],
     ["estimate", "--panel", "p.csv"], ["fetch-partd", "--state", "IN"]],
)
def test_every_subcommand_parses(argv):
    args = build_parser().parse_args(argv)
    assert callable(args.func)
