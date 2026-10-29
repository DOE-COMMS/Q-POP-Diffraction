"""
Single-domain PbTiO3: uniform P along each of the six tetragonal variants.

With no strain, the intensity at a reciprocal lattice point G is the structure
factor of the cell with its P-induced ionic shifts d_a,

    I(G) = |sum_a f_a exp(-i G . (R_a + d_a))|^2 / V^2,

computed here independently of the code.  A spontaneous strain e along P
moves the reflections with a component along P to G / (1 + e) and leaves the
others at G (to first order in e, as the code's phase exp(-i G . u) is).

    pytest tests/test_single_domain.py
"""

import numpy as np
import pytest

from test_planted import DX, A, FORM, POS, SHIFT, at, atom_in, peak, run, write_run

NX, NY, NZ = 40, 30, 20
PS = 0.75  # C/m^2
VARIANTS = {"c+": (0, 0, 1), "c-": (0, 0, -1), "a1+": (1, 0, 0), "a1-": (-1, 0, 0), "a2+": (0, 1, 0), "a2-": (0, -1, 0)}
REFLECTIONS = [(0, 0, 2), (2, 0, 0), (0, 2, 0), (1, 1, 1)]


def reciprocal(hkl):
    return 2 * np.pi * np.array(hkl, dtype=float) / A


def structure_factor(q, P):
    """sum_a f_a exp(-i q . (R_a + d_a)), positions in nm, P in C/m^2."""
    r = POS * A + SHIFT * np.asarray(P) / 1e-9
    return np.sum(FORM * np.exp(-1j * (r @ q)))


def single_domain(tmp_path, monkeypatch, P, hkl, strain=None):
    """Uniform P (C/m^2) around reflection hkl; strain (axis, e) adds u_axis = e * r_axis."""
    fields = np.zeros((3, NZ, NY, NX))
    for d in range(3):
        fields[d] = P[d]
    u = np.zeros((3, NZ, NY, NX))
    if strain is not None:
        d, e = strain
        r = np.arange((NX, NY, NZ)[d]) * DX
        u[d] = e * r.reshape([-1 if i == 2 - d else 1 for i in range(3)])
    write_run(tmp_path, fields, u, atom=atom_in(reciprocal(hkl)))
    return run(tmp_path, monkeypatch)


@pytest.mark.parametrize("hkl", REFLECTIONS, ids=lambda h: "".join(map(str, h)))
@pytest.mark.parametrize("variant", VARIANTS)
def test_structure_factor(tmp_path, monkeypatch, variant, hkl):
    P = PS * np.array(VARIANTS[variant])
    G = reciprocal(hkl)
    I, q = single_domain(tmp_path, monkeypatch, P, hkl)
    np.testing.assert_allclose(peak(I, q), G, atol=1e-9)
    np.testing.assert_allclose(at(I, q, G), abs(structure_factor(G, P)) ** 2 / A ** 6, rtol=1e-9)


def test_variant_symmetry(tmp_path, monkeypatch):
    # The cell is cubic: P along d seen at the 2 0 0 reflection along d is the same for every d.
    # A shift perpendicular to G is invisible.  The imaginary form factors break Friedel's law,
    # so P and -P differ along G.
    I = {}
    for name, P in list(VARIANTS.items()) + [("0", (0, 0, 0))]:
        for hkl in [(0, 0, 2), (2, 0, 0), (0, 2, 0)]:
            work = tmp_path / f"{name}_{''.join(map(str, hkl))}"
            work.mkdir()
            out, q = single_domain(work, monkeypatch, PS * np.array(P), hkl)
            I[name, hkl] = at(out, q, reciprocal(hkl))

    for sign in "+-":
        parallel = [I["c" + sign, (0, 0, 2)], I["a1" + sign, (2, 0, 0)], I["a2" + sign, (0, 2, 0)]]
        np.testing.assert_allclose(parallel, parallel[0], rtol=1e-12)
        perpendicular = [I["a1" + sign, (0, 0, 2)], I["a2" + sign, (0, 0, 2)], I["c" + sign, (2, 0, 0)]]
        np.testing.assert_allclose(perpendicular, I["0", (0, 0, 2)], rtol=1e-12)
    assert abs(I["c+", (0, 0, 2)] / I["c-", (0, 0, 2)] - 1) > 1e-3


@pytest.mark.parametrize("variant", ["c+", "a1+", "a2+"])
def test_spontaneous_strain(tmp_path, monkeypatch, variant):
    # e * L * |G_200| / 2 pi is an integer on every axis (0.1 * 40, 30, 20), so u = e r is periodic in
    # exp(-i G . u) and the shifted peak lies on the grid.  The code is first order in u, so the peak
    # is at G (1 - e); e is small enough that G e^2 is well inside half a grid step of G / (1 + e).
    e = 0.02
    axis = int(np.argmax(np.abs(VARIANTS[variant])))
    for hkl in [(0, 0, 2), (2, 0, 0), (0, 2, 0)]:
        work = tmp_path / "".join(map(str, hkl))
        work.mkdir()
        I, q = single_domain(work, monkeypatch, PS * np.array(VARIANTS[variant]), hkl, strain=(axis, e))
        G = reciprocal(hkl)
        expected = G.copy()
        expected[axis] = G[axis] / (1 + e)
        dq = 2 * np.pi / (np.array([NX, NY, NZ]) * DX * 1e9)
        assert np.all(np.abs(peak(I, q) - expected) < dq / 2), (hkl, peak(I, q), expected)
