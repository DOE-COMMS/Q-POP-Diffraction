"""
Planted-signal tests: synthetic FerroDyn runs with a known answer.

Each test writes a small run directory in FerroDyn layout (datasets shaped
(nz, ny, nx), SI units), runs the qpop-diffraction CLI on it and checks where
the intensity lands.  The grid is non-cubic so a wrong axis order is caught.

    pytest tests/
"""

import json
import sys

import h5py
import numpy as np
import pytest

from qpop_diffraction import SizeContext, load_diffraction_hdf5
from qpop_diffraction.cli import main

NX, NY, NZ = 48, 32, 20
DX = 1e-9
A = 0.4              # lattice constant in parameter.atom.in (nm)
G = 2 * np.pi * 2 / A  # (002), nm^-1

# PbTiO3-like 5-atom cell, one phase (examples/4_polar_vortex_PbTiO3): Pb, Ti, and the O on the
# x, y and z faces of the Ti octahedron.  SHIFT[a, d] is atom a's displacement along d (m) per C/m^2
# of P along d.
POS = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]])
FORM = np.array([19.3678 + 0.746466j, 40.5697 + 1.67439j, 8.035 + 0.021j, 8.035 + 0.021j, 8.035 + 0.021j])
SHIFT = np.array([[1.01, 1.01, 1.01], [-0.56, -0.56, -0.56],
                  [-3.66, -3.15, -3.15], [-3.15, -3.66, -3.15], [-3.15, -3.15, -3.66]]) * 1e-11


def atom_in(q00):
    """parameter.atom.in of the PbTiO3 cell, with the diffraction pattern computed around q00 (nm^-1)."""
    lines = [" ".join(str(c) for c in q00), f"{A} 0.0 0.0", f"0.0 {A} 0.0", f"0.0 0.0 {A}", str(len(POS))]
    lines += [" ".join(str(c) for c in x) for x in POS]
    lines += ["1"] + [f"{f.real} {f.imag}" for f in FORM]
    for s in SHIFT:
        lines += [" ".join(str(s[k] if d == k else 0.0) for d in range(3)) for k in range(3)]
    return "\n".join(lines) + "\n"


ATOM_IN = atom_in((0.0, 0.0, G))


def material(name, is_fe, a):
    return {
        "Name": name,
        "Ferroelectric": {"Is_FE": is_fe},
        "Mechanical": {
            "LatticeConstants": {"a": a, "b": a, "c": a, "ReferenceTemperature": 300.0},
            "ThermalExpansion": [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
        },
    }


def write_run(tmp_path, P, u, geometry="Bulk", layers=None, sub_a=4e-10, film_a=4e-10, phase=None, region=None,
              atom=ATOM_IN):
    """
    Write system_setting.json, materials/, data/ and the CLI inputs; P and u are (3, nz, ny, nx),
    which sets the grid, phase (nz, ny, nx) the fraction of the one phase and region (nz, ny, nx)
    the scattering weight.
    """
    names = ["Substrate", "Film"]
    nz, ny, nx = P.shape[1:]
    layer_types, layer_sizes = layers or ([1], [nz])
    setting = {
        "Mesh": {"NX": nx, "NY": ny, "NZ": nz, "DX": DX, "DY": DX, "DZ": DX},
        "Rotation": {"Enabled": False, "EulerAngles": [0.0, 0.0, 0.0]},
        "BoundaryConditions": {"GeometryType": geometry},
        "SimulationControls": {"Temperature": {"Initial": 300.0}},
        "Materials": {
            "NumberOfMaterials": 2,
            "MaterialNames": names,
            "NumberOfLayers": len(layer_types),
            "LayerTypes": layer_types,
            "LayerSizes": layer_sizes,
        },
        "Solvers": {},
    }
    (tmp_path / "system_setting.json").write_text(json.dumps(setting))
    (tmp_path / "materials").mkdir()
    (tmp_path / "materials" / "Substrate.json").write_text(json.dumps(material("Substrate", False, sub_a)))
    (tmp_path / "materials" / "Film.json").write_text(json.dumps(material("Film", True, film_a)))
    (tmp_path / "data").mkdir()

    mat = np.repeat(np.repeat(layer_types, layer_sizes), ny * nx).reshape(nz, ny, nx).astype(np.int32)
    with h5py.File(tmp_path / "data" / "materialField_0.h5", "w") as fh:
        fh["materialField"] = mat
    with h5py.File(tmp_path / "data" / "P_Lattice_0.h5", "w") as fh:
        fh.attrs["Timestep"] = 0
        for c, comp in zip("xyz", P):
            fh[f"P_Lattice_{c}"] = comp
    with h5py.File(tmp_path / "data" / "Displacement_0.h5", "w") as fh:
        for c, comp in zip("xyz", u):
            fh[f"Displacement_{c}"] = comp

    cfg = {
        "platt_file": "data/P_Lattice_0.h5",
        "displacement_file": "data/Displacement_0.h5",
        "output_file": "diffraction",
    }
    if phase is not None:
        with h5py.File(tmp_path / "data" / "phase.h5", "w") as fh:
            fh["phase_1"] = phase.astype(np.float32)
        cfg["phase_file"] = "data/phase.h5"
    if region is not None:
        with h5py.File(tmp_path / "data" / "region.h5", "w") as fh:
            fh.create_dataset("region", data=region.astype(np.float32), compression="gzip", shuffle=True)
        cfg["region_file"] = "data/region.h5"

    (tmp_path / "parameter.atom.in").write_text(atom)
    (tmp_path / "input.json").write_text(json.dumps(cfg))


def run(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", ["qpop-diffraction"])
    assert main() == 0
    out = load_diffraction_hdf5(tmp_path / "diffraction.h5")
    return out["intensity"][0], out["q_vector"]


def at(I, q, target):
    """Intensity at the grid point nearest to the q target."""
    d = sum((q[i] - target[i]) ** 2 for i in range(3))
    return I.flat[np.argmin(d)]


def peak(I, q):
    i = np.argmax(I)
    return np.array([q[0].flat[i], q[1].flat[i], q[2].flat[i]])


def zeros():
    return np.zeros((3, NZ, NY, NX))


def test_uniform_p_peaks_at_g(tmp_path, monkeypatch):
    P = zeros()
    P[2] = 0.5
    write_run(tmp_path, P, zeros())
    I, q = run(tmp_path, monkeypatch)
    np.testing.assert_allclose(peak(I, q), [0.0, 0.0, G], atol=1e-9)


def test_uz_ramp_shifts_peak(tmp_path, monkeypatch):
    # u_z = e * z in metres: the 002 peak moves to G (1 - e).  A missing m -> nm
    # conversion would leave it at G.
    e = 0.02
    u = zeros()
    u[2] = e * DX * np.arange(NZ)[:, None, None]
    write_run(tmp_path, zeros(), u)
    I, q = run(tmp_path, monkeypatch)
    dqz = 2 * np.pi / (NZ * DX * 1e9)
    assert abs(peak(I, q)[2] - G * (1 - e)) < dqz


def test_stripe_satellites_along_x_only(tmp_path, monkeypatch):
    # 180-degree Pz stripes with period L along x give satellites at G +- 2 pi / L
    # on q_x; nothing appears at the same offset on q_y.
    L = 16
    P = zeros()
    P[2] = np.where((np.arange(NX) // (L // 2)) % 2 == 0, 0.5, -0.5)[None, None, :]
    write_run(tmp_path, P, zeros())
    I, q = run(tmp_path, monkeypatch)
    k = 2 * np.pi / (L * DX * 1e9)
    on_x = at(I, q, (k, 0.0, G))
    on_y = at(I, q, (0.0, k, G))
    assert on_x > 1e-6 * I.max()
    assert on_x > 1e6 * on_y


def test_substrate_misfit(tmp_path):
    # Thin film on a substrate: the in-plane lattice is the substrate's,
    # referred to aC; the substrate is not a phase.
    write_run(tmp_path, zeros(), zeros(), geometry="Thin-film", layers=([0, 1], [8, 12]), sub_a=3.9e-10)
    params = SizeContext()
    assert params.read_json(tmp_path / "system_setting.json") == [1]
    dc = params.reference_strain(np.array([0.41, 0.41, 0.41]), 300.0)
    np.testing.assert_allclose(params.strainAvg, [0.39 / 0.41 - 1, 0.39 / 0.41 - 1, 0, 0, 0, 0])
    np.testing.assert_allclose(dc, [0.39 / 0.41 - 1, 0.4 / 0.41 - 1])


def test_film_lattice_differs_from_ac(tmp_path, monkeypatch):
    # FerroDyn strain is measured from the film's own lattice: with u = 0 the
    # 002 peak must sit at 4 pi / c_film, not at G of the atom file's aC.
    c = 0.392
    write_run(tmp_path, zeros(), zeros(), film_a=c * 1e-9)
    I, q = run(tmp_path, monkeypatch)
    dqz = 2 * np.pi / (NZ * DX * 1e9)
    assert abs(peak(I, q)[2] - 4 * np.pi / c) < dqz


def test_temperature_ramp(tmp_path):
    write_run(tmp_path, zeros(), zeros())
    params = SizeContext()
    params.read_json(tmp_path / "system_setting.json")
    params.setting["SimulationControls"]["Temperature"] = {
        "Initial": 300.0, "Transient": True, "Final": 700.0, "Interval": 4000, "Step": 25}
    assert [params.temperature(s) for s in (0, 3999, 4000, 7500, 67500, 99999)] == [300, 300, 325, 325, 700, 700]


def test_phase_and_region_weights(tmp_path, monkeypatch):
    # The weight of a cell is phase fraction x region: uniform 0.6 and 0.5 scale I by 0.3^2.
    # A region that is zero on half the x range halves the amplitude at G.
    P = zeros()
    P[2] = 0.5
    runs = {"plain": {}, "uniform": {"phase": np.full((NZ, NY, NX), 0.6), "region": np.full((NZ, NY, NX), 0.5)},
            "half": {"region": (np.arange(NX) < NX // 2)[None, None, :] * np.ones((NZ, NY, NX))}}
    I = {}
    for name, weights in runs.items():
        (tmp_path / name).mkdir()
        write_run(tmp_path / name, P, zeros(), **weights)
        I[name], q = run(tmp_path / name, monkeypatch)
    np.testing.assert_allclose(I["uniform"], 0.09 * I["plain"], rtol=1e-6)
    np.testing.assert_allclose(at(I["half"], q, (0.0, 0.0, G)), 0.25 * at(I["plain"], q, (0.0, 0.0, G)), rtol=1e-6)


def test_rotated_run_rejected(tmp_path):
    write_run(tmp_path, zeros(), zeros())
    setting = json.loads((tmp_path / "system_setting.json").read_text())
    setting["Rotation"]["Enabled"] = True
    (tmp_path / "system_setting.json").write_text(json.dumps(setting))
    with pytest.raises(ValueError):
        SizeContext().read_json(tmp_path / "system_setting.json")
