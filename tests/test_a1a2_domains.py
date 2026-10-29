"""
a1/a2 domains in a PbTiO3 film: in-plane P along x or y in stripes whose walls
are (110) planes, on a non-scattering substrate with air above, seen at 2 0 0.

A cell (i, j) is a1 when (i + j) mod p < p / 2, so the domains repeat every
p * dx / sqrt(2) along [110] and each is w = p * dx / (2 sqrt 2) wide.  With
F1, F2 the structure factors of the two variants (see test_single_domain),
the kinematic answers are:

  - thickness T = t * dz: along q_z through G the film is a uniform slab, so
        I(G + k_z) / I(G) = sinc^2(k_z T / 2) |F(G + k_z)|^2 / |F(G)|^2,  F = (F1 + F2) / 2,
    with zeros at k_z = 2 pi m / T;
  - domain width: the first satellites sit at k1 = (2 pi / (p dx)) (1, 1, 0), with
        I(G + k1) / I(G) = S(k1)^2 |c1|^2 |F1(G + k1) - F2(G + k1)|^2 / |F(G)|^2,
    S the voxel shape factor and c1 the first Fourier coefficient of the square wave;
    nothing appears along (1, -1, 0), and the even orders vanish for equal widths.

The code expands exp(-i k . d) to first order in the ionic shifts d, which
costs O((k d)^2) in the satellite intensity; the rod has k . d = 0 and is exact.

    pytest tests/test_a1a2_domains.py
"""

import numpy as np
import pytest

from test_planted import A, at, atom_in, run, write_run
from test_single_domain import PS, structure_factor

NX, NY, NZ = 32, 48, 32
SUBSTRATE = 8
G = np.array([2 * np.pi * 2 / A, 0.0, 0.0])
F1 = lambda q: structure_factor(q, (PS, 0.0, 0.0))
F2 = lambda q: structure_factor(q, (0.0, PS, 0.0))


def stripes(tmp_path, monkeypatch, t, p):
    """a1/a2 stripes of period p (cells along x and y) in a film t cells thick."""
    i = np.arange(NX)[None, None, :]
    j = np.arange(NY)[None, :, None]
    a1 = np.broadcast_to(((i + j) % p) < p // 2, (NZ, NY, NX))
    P = np.zeros((3, NZ, NY, NX))
    P[0] = PS * a1
    P[1] = PS * ~a1
    write_run(tmp_path, P, np.zeros_like(P), geometry="Thin-film",
              layers=([0, 1, 0], [SUBSTRATE, t, NZ - SUBSTRATE - t]), atom=atom_in(G))
    return run(tmp_path, monkeypatch)


@pytest.mark.parametrize("p", [8, 16])
@pytest.mark.parametrize("t", [8, 16])
def test_thickness_fringes(tmp_path, monkeypatch, t, p):
    I, q = stripes(tmp_path, monkeypatch, t, p)
    F = lambda q: (F1(q) + F2(q)) / 2
    I0 = at(I, q, G)
    for kz in 2 * np.pi * np.fft.fftfreq(NZ, 1.0):
        qz = G + np.array([0.0, 0.0, kz])
        expected = np.sinc(kz * t / (2 * np.pi)) ** 2 * abs(F(qz)) ** 2 / abs(F(G)) ** 2
        if expected > 1e-8:
            np.testing.assert_allclose(at(I, q, qz) / I0, expected, rtol=1e-6)
        elif expected < 1e-20:
            assert at(I, q, qz) < 1e-20 * I0, kz


@pytest.mark.parametrize("p", [8, 16])
@pytest.mark.parametrize("t", [8, 16])
def test_domain_satellites(tmp_path, monkeypatch, t, p):
    I, q = stripes(tmp_path, monkeypatch, t, p)
    kappa = 2 * np.pi / p
    k1 = kappa * np.array([1.0, 1.0, 0.0])
    c1 = np.mean(np.where(np.arange(p) < p // 2, 1.0, 0.0) * np.exp(-2j * np.pi * np.arange(p) / p))
    S = np.sinc(kappa / (2 * np.pi)) ** 2
    expected = S ** 2 * abs(c1) ** 2 * abs(F1(G + k1) - F2(G + k1)) ** 2 / abs((F1(G) + F2(G)) / 2) ** 2

    I0 = at(I, q, G)
    np.testing.assert_allclose(at(I, q, G + k1) / I0, expected, rtol=1e-2)
    np.testing.assert_allclose(at(I, q, G - k1) / I0, abs(c1) ** 2 * S ** 2 * abs(F1(G - k1) - F2(G - k1)) ** 2
                               / abs((F1(G) + F2(G)) / 2) ** 2, rtol=1e-2)
    assert at(I, q, G + kappa * np.array([1.0, -1.0, 0.0])) < 1e-20 * I0
    assert at(I, q, G + 2 * k1) < 1e-20 * I0
