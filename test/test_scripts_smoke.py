"""Smoke tests for the ``scripts/`` and ``launchers/`` entry points.

The launchers are "edit the variables at the top, then run" scripts that do
real work at module top level, so they cannot be imported safely.  Instead we
byte-compile every entry point (guaranteeing it parses and imports its
dependencies' names) and run ``--help`` on the argparse-based CLI scripts.
"""

import os
import py_compile
import subprocess
import sys

import pytest

from conftest import PACKAGE_ROOT

SCRIPTS_DIR = PACKAGE_ROOT / "scripts"
LAUNCHERS_DIR = PACKAGE_ROOT / "launchers"

ALL_ENTRY_POINTS = sorted(SCRIPTS_DIR.glob("*.py")) + sorted(LAUNCHERS_DIR.glob("*.py"))

# CLI scripts whose argparse help can be exercised without side effects.
ARGPARSE_SCRIPTS = [
    "lelantos_pixel.py",
    "lelantos_prop.py",
    "lelantos_cut_dla_catalog.py",
    "lelantos_cut_qso_catalog.py",
    "lelantos_show_qso_catalog.py",
    "lelantos_void_statistics.py",
]


def test_entry_points_discovered():
    # 10 scripts + 17 launchers are shipped; guard against silent globbing loss
    assert len(ALL_ENTRY_POINTS) >= 25


@pytest.mark.parametrize(
    "path", ALL_ENTRY_POINTS, ids=[p.name for p in ALL_ENTRY_POINTS]
)
def test_entry_point_compiles(path):
    """Every script/launcher byte-compiles (no syntax errors)."""
    py_compile.compile(str(path), doraise=True)


@pytest.mark.parametrize("script", ARGPARSE_SCRIPTS)
def test_argparse_script_help(script):
    """Each argparse CLI prints usage and exits 0 for ``--help``."""
    env = dict(os.environ, PYTHONPATH=str(PACKAGE_ROOT))
    result = subprocess.run(
        [sys.executable, str(SCRIPTS_DIR / script), "--help"],
        cwd=str(PACKAGE_ROOT),
        env=env,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert "usage" in result.stdout.lower()
