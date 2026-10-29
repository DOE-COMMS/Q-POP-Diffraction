#!/usr/bin/env python3
"""
Rebuild the four legacy examples (Examples.tar.gz, text inputs of the old C++
code) as runs that qpop-diffraction reads directly.

Each case gets the FerroDyn layout: system_setting.json, materials/,
data/materialField_0.h5, data/P_Lattice_0.h5, data/Displacement_0.h5,
data/phase.h5, data/region.h5, parameter.atom.in and input.json.  Fields are
float32 with shuffle + gzip; reference intensities are float64.

How the legacy inputs map:
  - displace.in (nm)       -> Displacement_x/y/z in m
  - strucOrd.in (C/m^2)    -> P_Lattice_x/y/z; zero when there is none
  - phaseFra.in            -> phase.h5:/phase_<id>
  - region.in              -> region.h5:/region; without one, the legacy default
                              tanh window, so results stay comparable
  - parameter.system.in    -> Mesh; one Is_FE material per phase whose lattice is
                              diag(aC), so the average strain is zero as before
  - parameter.atom.in      -> copied; with nStruc = 0, rewritten with zero
                              P -> shift coefficients for the 3 P components
  - I.00000000.dat         -> reference_legacy.h5:/intensity (nz, ny, nx), a regression
                              baseline from a legacy-format run (not in the tarball;
                              taken from --legacy-outputs, the extracted Examples/)

    python examples/convert_legacy_examples.py Examples.tar.gz [--legacy-outputs Examples]
"""

import argparse
import io
import json
import tarfile
from pathlib import Path

import h5py
import numpy as np

from qpop_diffraction import AtomList

CASES = {
    "1. Two-phase mixture/isotropic": "1_two_phase_mixture/isotropic",
    "1. Two-phase mixture/anisotropic": "1_two_phase_mixture/anisotropic",
    "2. Simple two-variant twinned crystal": "2_two_variant_twin",
    "3. Three-variant twinned hexagonal crystal": "3_three_variant_hexagonal_twin",
    "4. Polar vortex in PbTiO3/ordered": "4_polar_vortex_PbTiO3/ordered",
    "4. Polar vortex in PbTiO3/disordered": "4_polar_vortex_PbTiO3/disordered",
}
FIELD = dict(compression="gzip", compression_opts=9, shuffle=True)


def numbers(text):
    """Numbers of a legacy parameter file, comments after '!' removed."""
    return [float(t) for line in text.splitlines() for t in line.split("!")[0].split()]


def read_system(text):
    v = numbers(text)
    return {"l": v[0:3], "n": [int(x) for x in v[3:6]], "ns": int(v[6]), "nf": int(v[7]),
            "nPhase": int(v[8]), "nStruc": int(v[9]), "strain": v[10:16]}


def read_field(text, ncomp, n):
    """Legacy 4-D text field -> (ncomp, nx, ny, nz)."""
    data = np.loadtxt(io.StringIO(text), skiprows=1)
    i, j, k = (data[:, d].astype(int) - 1 for d in range(3))
    field = np.zeros((ncomp, *n))
    field[:, i, j, k] = data[:, 3:3 + ncomp].T
    return field


def legacy_window(n):
    """Default region of the legacy code (nf = 0): a tanh window per axis, normalised to 1."""
    w = np.ones(n)
    for d, m in enumerate(n):
        c1, c2, t = (m + 1) / 2 - m / np.pi, (m + 1) / 2 + m / np.pi, m / 8
        i = np.arange(1, m + 1)
        w = w * (np.tanh((i - c1) / t) - np.tanh((i - c2) / t)).reshape([-1 if e == d else 1 for e in range(3)])
    return w / w.max()


def write_fields(path, fields, dtype=np.float32, **attrs):
    """Datasets in FerroDyn layout (nz, ny, nx); fields are (nx, ny, nz)."""
    with h5py.File(path, "w") as fh:
        for key, value in attrs.items():
            fh.attrs[key] = value
        for name, f in fields.items():
            fh.create_dataset(name, data=np.ascontiguousarray(f.transpose(2, 1, 0), dtype=dtype),
                              **(FIELD if dtype != np.float64 else dict(compression="gzip", compression_opts=9)))


def atom_file(text, path, nPhase, nStruc, title):
    """Copy parameter.atom.in, adding zero P -> shift coefficients when the legacy case had none."""
    if nStruc == 3:
        path.write_text(text)
        return
    atoms = AtomList()
    path.write_text(text)
    assert atoms.read(str(path), nPhase, nStruc)
    out = [f"! {title}",
           f"{atoms.q00[0]} {atoms.q00[1]} {atoms.q00[2]}    !reciprocal lattice point (nm^-1)"]
    out += [f"{a[0]} {a[1]} {a[2]}    !lattice vector (nm)" for a in atoms.aC]
    out += [f"{atoms.nAtom}    !# of atoms in a unit cell"]
    out += [f"{x[0]} {x[1]} {x[2]}    !reduced coordinates of atom #{i + 1}" for i, x in enumerate(atoms.xAtom)]
    for p in range(nPhase):
        out += [f"{p + 1}    !phase ID"]
        out += [f"{atoms.fAtomR[p, i]} {atoms.fAtomI[p, i]}    !form factor of atom #{i + 1}" for i in range(atoms.nAtom)]
        out += ["0.0 0.0 0.0    !P -> shift coefficient (none in this example)"] * (3 * atoms.nAtom)
    path.write_text("\n".join(out) + "\n")


def lattice(path, nPhase):
    """aC (nm) of a converted parameter.atom.in."""
    atoms = AtomList()
    assert atoms.read(str(path), nPhase, 3)
    return atoms.aC


def convert(tar, source, case, legacy_outputs=None):
    files = {Path(m.name).name: m for m in tar.getmembers()
             if m.isfile() and Path(m.name).parent.as_posix().endswith(source) and "Zone.Identifier" not in m.name}
    text = {name: tar.extractfile(m).read().decode() for name, m in files.items()}
    for name in ("I.00000000.dat", "region.00000000.dat"):
        if legacy_outputs is not None and (legacy_outputs / source / name).is_file():
            text[name] = (legacy_outputs / source / name).read_text()
    system = read_system(text["parameter.system.in"])
    if system["ns"] or system["nf"] or any(system["strain"]):
        raise ValueError(f"{source}: films and average strain are not handled by this converter")
    n, nPhase, title = system["n"], system["nPhase"], text["parameter.system.in"].splitlines()[0].lstrip("! ")
    print(f"{source} -> {case}: grid {n}, {nPhase} phase(s)")

    (case / "materials").mkdir(parents=True, exist_ok=True)
    (case / "data").mkdir(exist_ok=True)
    atom_file(text["parameter.atom.in"], case / "parameter.atom.in", nPhase, system["nStruc"], title)
    aC = np.diag(lattice(case / "parameter.atom.in", nPhase))

    names = ["Film"] if nPhase == 1 else [f"Phase{p + 1}" for p in range(nPhase)]
    for name in names:
        (case / "materials" / f"{name}.json").write_text(json.dumps({
            "Name": name,
            "Ferroelectric": {"Is_FE": True},
            "Mechanical": {
                "LatticeConstants": {"a": aC[0] * 1e-9, "b": aC[1] * 1e-9, "c": aC[2] * 1e-9, "ReferenceTemperature": 300.0},
                "ThermalExpansion": [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
            },
        }, indent=4) + "\n")
    (case / "system_setting.json").write_text(json.dumps({
        "Mesh": {"NX": n[0], "NY": n[1], "NZ": n[2],
                 "DX": system["l"][0] / n[0] * 1e-9, "DY": system["l"][1] / n[1] * 1e-9, "DZ": system["l"][2] / n[2] * 1e-9},
        "Rotation": {"Enabled": False},
        "BoundaryConditions": {"GeometryType": "Bulk"},
        "SimulationControls": {"Temperature": {"Initial": 300.0}},
        "Materials": {"NumberOfMaterials": nPhase, "MaterialNames": names,
                      "NumberOfLayers": 1, "LayerTypes": [0], "LayerSizes": [n[2]]},
        "Solvers": {},
    }, indent=4) + "\n")

    write_fields(case / "data" / "materialField_0.h5", {"materialField": np.zeros(n)}, dtype=np.int8)
    P = read_field(text["strucOrd.in"], 3, n) if "strucOrd.in" in text else np.zeros((3, *n))
    write_fields(case / "data" / "P_Lattice_0.h5", {f"P_Lattice_{c}": P[d] for d, c in enumerate("xyz")}, Timestep=0)
    cfg = {"platt_file": "data/P_Lattice_0.h5", "displacement_file": None, "phase_file": None,
           "region_file": "data/region.h5", "x_limits": None, "y_limits": None, "z_limits": None,
           "output_file": "diffraction"}
    if "displace.in" in text:
        u = read_field(text["displace.in"], 3, n) * 1e-9
        write_fields(case / "data" / "Displacement_0.h5", {f"Displacement_{c}": u[d] for d, c in enumerate("xyz")})
        cfg["displacement_file"] = "data/Displacement_0.h5"
    if "phaseFra.in" in text:
        f = read_field(text["phaseFra.in"], nPhase, n)
        write_fields(case / "data" / "phase.h5", {f"phase_{p + 1}": f[p] for p in range(nPhase)})
        cfg["phase_file"] = "data/phase.h5"
    region = read_field(text["region.in"], 1, n)[0] if "region.in" in text else legacy_window(n)
    write_fields(case / "data" / "region.h5", {"region": region})
    if "region.00000000.dat" in text:
        written = read_field(text["region.00000000.dat"], 1, n)[0]
        print(f"  window vs legacy region.00000000.dat: max |diff| = {np.abs(region - written).max():.1e}")
    (case / "input.json").write_text(json.dumps(cfg, indent=4) + "\n")

    if "I.00000000.dat" in text:
        write_fields(case / "reference_legacy.h5", {"intensity": read_field(text["I.00000000.dat"], 1, n)[0]},
                     dtype=np.float64)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("tarball", type=Path, help="Examples.tar.gz of the legacy code")
    parser.add_argument("--output", "-o", type=Path, default=Path(__file__).resolve().parent,
                        help="Directory to write the cases into (default: this examples/ directory)")
    parser.add_argument("--legacy-outputs", type=Path, default=None,
                        help="Extracted Examples/ holding legacy I.00000000.dat / region.00000000.dat to keep as references")
    args = parser.parse_args()
    with tarfile.open(args.tarball) as tar:
        for source, case in CASES.items():
            convert(tar, source, args.output / case, args.legacy_outputs)


if __name__ == "__main__":
    main()
