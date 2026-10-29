"""
The multi-GPU program qpop_diffraction_mgpu vs the Python path on the planted runs of test_planted.py.

Needs the multi-GPU build (bin/qpop_diffraction_mgpu) and a launch command
inside a GPU allocation, e.g.

    QPOP_MGPU_LAUNCH="srun -n 3 --gpus-per-node=4 --gpu-bind=none" pytest tests/test_mgpu.py

Three ranks on NY = 32 rows exercise an uneven split and the wrapped y-fftshift.
Skipped when QPOP_MGPU_LAUNCH is not set or the program is not built.
"""

import os
import shlex
import subprocess
import sys

import h5py
import numpy as np
import pytest

from qpop_diffraction.cli import MGPU_PROGRAM, main
from test_planted import DX, NX, NY, NZ, write_run, zeros

LAUNCH = os.environ.get("QPOP_MGPU_LAUNCH")
pytestmark = pytest.mark.skipif(not LAUNCH or not MGPU_PROGRAM.exists(),
                                reason="set QPOP_MGPU_LAUNCH and build bin/qpop_diffraction_mgpu")


def intensity(tmp_path, monkeypatch, mgpu, limits=None):
    """Run the CLI (and qpop_diffraction_mgpu for mgpu); return /intensity and the q axes."""
    if limits is not None:
        cfg = (tmp_path / "input.json").read_text().replace('"output_file"', f'{limits}, "output_file"')
        (tmp_path / "input.json").write_text(cfg)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", ["qpop-diffraction"] + (["--mgpu"] if mgpu else []))
    assert main() == 0
    if mgpu:
        subprocess.run(shlex.split(LAUNCH) + [str(MGPU_PROGRAM), str(tmp_path / "diffraction_params.h5")], check=True)
    with h5py.File(tmp_path / "diffraction.h5", "r") as fh:
        return fh["intensity"][:], [fh[f"q_{c}"][:] for c in "xyz"]


def compare(tmp_path, monkeypatch, P, u, limits=None, **run):
    ref, mgpu = tmp_path / "ref", tmp_path / "mgpu"
    for d in (ref, mgpu):
        d.mkdir()
        write_run(d, P, u, **run)
    I_ref, q_ref = intensity(ref, monkeypatch, False, limits)
    I, q = intensity(mgpu, monkeypatch, True, limits)
    for a, b in zip(q, q_ref):
        np.testing.assert_allclose(a, b, rtol=1e-14)
    np.testing.assert_allclose(I, I_ref, rtol=1e-9, atol=1e-12 * I_ref.max())


def stripes():
    P = zeros()
    P[0] = 0.1 * np.sin(2 * np.pi * np.arange(NX) / 16)[None, None, :]
    P[2] = np.where((np.arange(NX) // 8) % 2 == 0, 0.5, -0.5)[None, None, :]
    return P


def test_bulk(tmp_path, monkeypatch):
    u = zeros()
    u[2] = 0.01 * DX * np.arange(NZ)[:, None, None]
    u[0] = 0.003 * DX * np.cos(2 * np.pi * np.arange(NZ) / NZ)[:, None, None]
    compare(tmp_path, monkeypatch, stripes(), u)


def test_film_on_substrate(tmp_path, monkeypatch):
    # Non-scattering substrate, u_z reference ramp from c_film != c_sub
    compare(tmp_path, monkeypatch, stripes(), zeros(), geometry="Thin-film",
            layers=([0, 1], [8, 12]), sub_a=3.9e-10, film_a=4.05e-10)


def test_cropped(tmp_path, monkeypatch):
    # The ramp is summed over the full z before the crop
    compare(tmp_path, monkeypatch, stripes(), zeros(), limits='"x_limits": [5, 36], "y_limits": [3, 27], "z_limits": [4, 17]',
            geometry="Thin-film", layers=([0, 1], [8, 12]), sub_a=3.9e-10, film_a=4.05e-10)


def test_phase_and_region_files(tmp_path, monkeypatch):
    # Fractional phase and a smooth region weight, read as gzip float32
    z, y, x = np.meshgrid(np.arange(NZ), np.arange(NY), np.arange(NX), indexing="ij")
    phase = 0.5 + 0.4 * np.cos(2 * np.pi * x / NX)
    region = np.exp(-((y - NY / 2) / (NY / 4)) ** 2)
    compare(tmp_path, monkeypatch, stripes(), zeros(), phase=phase, region=region)
