# Inputs

A run reads three files that you pass on the command line, plus the FerroDyn files they point to:

```bash
qpop-diffraction -i input.json -s system_setting.json -a parameter.atom.in
```

These are also the defaults, so `qpop-diffraction` alone works in a directory that holds all three.

## input.json

```json
{
    "platt_file": "data/P_Lattice_15000.h5",
    "displacement_file": "data/Displacement_15000.h5",
    "x_limits": null,
    "y_limits": null,
    "z_limits": null,
    "phase_file": null,
    "region_file": null,
    "output_file": "diffraction_15000"
}
```

- `platt_file` (required): the polarization snapshot.
- `displacement_file`: the elastic displacement of the same step. `null` means zero displacement. FerroDyn can save the elastic fields at a different interval from P, so pick a step that has both.
- `x_limits`, `y_limits`, `z_limits`: `[start, end]` cell indices, 0-based and inclusive, to compute a sub-volume. `null` means the whole axis.
- `phase_file`: optional phase fractions (see below). `null` means one phase per ferroelectric material, taken from the material field.
- `region_file`: optional scattering weight (see below). `null` means 1.
- `output_file`: base name of the result (default `diffraction`).

Relative paths are relative to the directory you run in. Other keys, such as `_comment`, are ignored.

Leave `z_limits` at `null` for films. The FFT is periodic, so cropping z to the film wraps its top onto its bottom. The `u_z` ramp then makes a phase jump there, which smears intensity along q_z and shifts the peaks. The non-scattering substrate and air above the film act as padding.

## system_setting.json and the FerroDyn files

The FerroDyn input of the run. Read from it:

- `Mesh`: `NX, NY, NZ` and `DX, DY, DZ` (m)
- `Materials`: `MaterialNames` and `LayerTypes`
- `BoundaryConditions.GeometryType`: `Thin-film` or bulk
- `SimulationControls.Temperature`: `Initial`, and `Transient, Step, Interval, Final` for a ramp
- `Solvers.ElasticStrain.InitialConditions`: `MisfitStrain` and `AppliedStrain`
- `Rotation.Enabled`: must be false; rotated runs are rejected

Read from beside it:

- `materials/<name>.json`, for every material: `Ferroelectric.Is_FE` (ferroelectric materials scatter, the others do not), `Mechanical.LatticeConstants` (`a, b, c` in m, and `ReferenceTemperature`) and `Mechanical.ThermalExpansion` (3×3, the diagonal is used).
- `data/materialField_0.h5`: the integer material index of every cell, dataset `materialField`.

The field files:

- `P_Lattice_<step>.h5`: datasets `P_Lattice_x/y/z` (or the older `PLatt_x/y/z`) in C/m², and the attribute `Timestep`, which sets the temperature.
- `Displacement_<step>.h5`: datasets `Displacement_x/y/z` in m.

## HDF5 field layout

Every field dataset has the shape `(nz, ny, nx)` of FerroDyn, with x varying fastest, and covers the full grid of `Mesh`. Any other shape is rejected. Files written before 2026-05-26 are `(nx, ny, nz)`; convert them with `qpop_monitor.migrate_hdf5_layout`. Any floating-point type and any standard HDF5 compression are accepted, and values are converted to double on read. Float32 with shuffle + gzip keeps example data small.

## Phase fractions: phase_file

By default every material with `Is_FE` true is one phase, in `MaterialNames` order, and a cell belongs entirely to the phase of its material. A `phase_file` replaces this with fractions. It holds datasets `phase_1` to `phase_<nPhase>`, one per phase ID of `parameter.atom.in`, laid out like the fields. Cells of non-ferroelectric materials still do not scatter.

## Scattering weight: region_file

A `region_file` holds one dataset `region`, laid out like the fields, which multiplies the scattering of every phase in each cell. Use it to restrict scattering to a grain or a probed volume, or to apply a smooth window that suppresses truncation streaks from the box edges.

## parameter.atom.in

The crystal: lattice, atoms, form factors and how P moves the atoms. The file is read as a stream of numbers. Anything after `!` on a line is a comment, and line breaks do not matter. The numbers, in order:

- `q00`: 3 numbers, the point in reciprocal space (nm⁻¹) to compute around. The pattern is centred on the reciprocal lattice point nearest to it, so e.g. `0 0 31.416` selects 002 of a 0.4 nm cubic cell.
- `aC`: 9 numbers, the lattice vectors a, b, c (nm), one per row. Choose them close to the film's average lattice (see [Method](method.md)).
- `nAtom`: the number of atoms in the unit cell.
- `nAtom` × 3 numbers: the reduced coordinates of each atom.
- One block per phase, `nPhase` blocks in all, where `nPhase` is the number of `Is_FE` materials.

Each phase block holds:

- the phase ID: 1, 2, ... for the ferroelectric materials in `MaterialNames` order;
- `nAtom` × 2 numbers: the real and imaginary parts of each atom's form factor;
- `nAtom` × 3 × 3 numbers: for each atom, three lines (for P_x, P_y, P_z), each giving the atom's displacement vector (u1 u2 u3) per unit polarization, in m per C/m² (m³/C).

In the PbTiO3 example, the Ti line for P_z is `0.0 0.0 -0.56E-11`. At P_z = 0.75 C/m², Ti is therefore displaced by −4.2 pm along z on top of the cell's displacement field. See `examples/4_polar_vortex_PbTiO3/*/parameter.atom.in` for a complete five-atom file.
