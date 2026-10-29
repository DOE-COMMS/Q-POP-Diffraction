"""
The converted legacy examples in examples/ (see examples/convert_legacy_examples.py).

Every example must run, and the disordered polar vortex must reproduce
reference_legacy.h5, the intensity of a legacy-format run on the same inputs.

    pytest tests/test_examples.py
"""

import shutil
import sys
from pathlib import Path

import h5py
import numpy as np
import pytest

from qpop_diffraction.cli import main

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"
CASES = sorted(p.parent.relative_to(EXAMPLES) for p in EXAMPLES.glob("**/system_setting.json"))


def run(case, tmp_path, monkeypatch):
    work = tmp_path / case
    shutil.copytree(EXAMPLES / case, work)
    monkeypatch.chdir(work)
    monkeypatch.setattr(sys, "argv", ["qpop-diffraction"])
    assert main() == 0
    with h5py.File(work / "diffraction.h5", "r") as fh:
        return fh["intensity"][:]


@pytest.mark.parametrize("case", CASES, ids=str)
def test_example_runs(case, tmp_path, monkeypatch):
    I = run(case, tmp_path, monkeypatch)
    assert np.isfinite(I).all() and I.max() > 1e-30


def test_polar_vortex_matches_legacy(tmp_path, monkeypatch):
    case = Path("4_polar_vortex_PbTiO3/disordered")
    I = run(case, tmp_path, monkeypatch)
    with h5py.File(EXAMPLES / case / "reference_legacy.h5", "r") as fh:
        ref = fh["intensity"][:]
    strong = ref > 1e-8 * ref.max()
    np.testing.assert_allclose(I[strong], ref[strong], rtol=1e-4)
