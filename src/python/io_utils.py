"""
I/O utilities: Q-POP-FerroDyn HDF5 input (fields, phase fractions, scattering
weight), diffraction HDF5/XDMF output and the parameters of the multi-GPU program.
"""

from __future__ import annotations

import os

import h5py
import numpy as np

from .constants import L0


# ── Q-POP-FerroDyn HDF5 input ───────────────────────────────────────────────
# FerroDyn writes one dataset per component, shaped (nz, ny, nx) with x fastest,
# in SI units.  Everything here returns (x, y, z) order, C-flattened, as the
# diffraction code expects.


def crop_slices(nx, ny, nz, x_lim, y_lim, z_lim):
    """0-indexed inclusive limits (None → full axis) as (x, y, z) slices."""
    return tuple(
        slice(lim[0], lim[1] + 1) if lim is not None else slice(0, n)
        for lim, n in ((x_lim, nx), (y_lim, ny), (z_lim, nz))
    )


def _load_field(fh, name, nx, ny, nz, sl):
    """Read one (nz, ny, nx) FerroDyn dataset and return it cropped in (x, y, z) order."""
    field = fh[name]
    if field.shape != (nz, ny, nx):
        raise ValueError(
            f"Dataset '{name}' in {fh.filename} has shape {field.shape}, expected "
            f"(nz, ny, nx) = ({nz}, {ny}, {nx}) from system_setting.json.  Files "
            f"written before 2026-05-26 are (nx, ny, nz): convert them with "
            f"qpop_monitor.migrate_hdf5_layout."
        )
    return np.asarray(field[sl[2], sl[1], sl[0]], dtype=np.float64).transpose(2, 1, 0)


def _load_vector_hdf5(filepath, prefixes, nx, ny, nz, x_lim, y_lim, z_lim):
    """
    Load ``{prefix}_x/y/z`` from a FerroDyn HDF5 file, trying each prefix in
    order, and return the cropped field as (3, nx_s, ny_s, nz_s).
    """
    sl = crop_slices(nx, ny, nz, x_lim, y_lim, z_lim)
    with h5py.File(filepath, "r") as fh:
        prefix = next(p for p in prefixes if f"{p}_x" in fh)
        return np.stack([_load_field(fh, f"{prefix}_{c}", nx, ny, nz, sl) for c in "xyz"])


def load_platt_hdf5(filepath, nx, ny, nz, x_lim=None, y_lim=None, z_lim=None):
    """
    Load a FerroDyn ``P_Lattice_<step>.h5`` (or older ``PLatt_<step>.h5``).

    Returns P in C/m², flattened in C order: [Px | Py | Pz].
    """
    return _load_vector_hdf5(filepath, ("P_Lattice", "PLatt"), nx, ny, nz, x_lim, y_lim, z_lim).ravel()


def load_displacement_hdf5(filepath, nx, ny, nz, x_lim=None, y_lim=None, z_lim=None):
    """
    Load a FerroDyn ``Displacement_<step>.h5``.

    FerroDyn writes u in metres; returns u in nm, flattened in C order: [ux | uy | uz].
    """
    return _load_vector_hdf5(filepath, ("Displacement",), nx, ny, nz, x_lim, y_lim, z_lim).ravel() / L0


def load_phase_hdf5(filepath, nPhase, nx, ny, nz, x_lim=None, y_lim=None, z_lim=None):
    """
    Load phase fractions ``phase_1`` .. ``phase_<nPhase>`` (phase IDs of
    parameter.atom.in), cropped, as (nPhase, nx_s, ny_s, nz_s).
    """
    sl = crop_slices(nx, ny, nz, x_lim, y_lim, z_lim)
    with h5py.File(filepath, "r") as fh:
        return np.stack([_load_field(fh, f"phase_{p + 1}", nx, ny, nz, sl) for p in range(nPhase)])


def load_region_hdf5(filepath, nx, ny, nz, x_lim=None, y_lim=None, z_lim=None):
    """Load the scattering weight ``region``, cropped, as (nx_s, ny_s, nz_s)."""
    sl = crop_slices(nx, ny, nz, x_lim, y_lim, z_lim)
    with h5py.File(filepath, "r") as fh:
        return _load_field(fh, "region", nx, ny, nz, sl)


def load_material_field(filepath, nx, ny, nz):
    """Load FerroDyn's ``materialField_0.h5`` uncropped, as integer material indices in (x, y, z) order."""
    with h5py.File(filepath, "r") as fh:
        return _load_field(fh, "materialField", nx, ny, nz, crop_slices(nx, ny, nz, None, None, None)).astype(int)


def write_diffraction_output(base_name, I, qx, qy, qz, step):
    """
    Write ``{base_name}.h5`` and the ParaView sidecar ``{base_name}.xdmf``.

    The HDF5 holds ``/intensity`` in (nz, ny, nx) order, as FerroDyn stores
    fields, and the 1-D axes ``/q_x``, ``/q_y``, ``/q_z`` in nm^-1.  The q grid
    is rectilinear, so the XDMF is a 3DRectMesh; take log10 in ParaView.
    """
    write_diffraction_axes(base_name, qx, qy, qz, step)
    with h5py.File(f"{base_name}.h5", "a") as fh:
        fh["intensity"] = I.transpose(2, 1, 0)


def write_diffraction_axes(base_name, qx, qy, qz, step):
    """
    Write ``{base_name}.h5`` without ``/intensity`` (attributes and q axes)
    and the ParaView sidecar ``{base_name}.xdmf``.  The multi-GPU program adds
    ``/intensity`` to this file.
    """
    nx, ny, nz = len(qx), len(qy), len(qz)
    print(f"Writing diffraction output → {base_name}.h5")
    with h5py.File(f"{base_name}.h5", "w") as fh:
        fh.attrs["Timestep"] = step
        fh.attrs["q_units"] = "nm^-1"
        fh["q_x"], fh["q_y"], fh["q_z"] = qx, qy, qz

    h5 = os.path.basename(f"{base_name}.h5")
    item = '<DataItem Dimensions="{}" NumberType="Float" Precision="8" Format="HDF">' + h5 + ':/{}</DataItem>'
    with open(f"{base_name}.xdmf", "w", encoding="utf-8") as fh:
        fh.write(f"""<?xml version="1.0" ?>
<Xdmf Version="3.0">
  <Domain>
    <Grid Name="DiffractionPattern" GridType="Uniform">
      <Topology TopologyType="3DRectMesh" Dimensions="{nz} {ny} {nx}"/>
      <Geometry GeometryType="VXVYVZ">
        {item.format(nx, "q_x")}
        {item.format(ny, "q_y")}
        {item.format(nz, "q_z")}
      </Geometry>
      <Attribute Name="Intensity" AttributeType="Scalar" Center="Node">
        {item.format(f"{nz} {ny} {nx}", "intensity")}
      </Attribute>
    </Grid>
  </Domain>
</Xdmf>
""")
    print(f"Writing ParaView sidecar  → {base_name}.xdmf")


def write_mgpu_params(path, full, limits, spacing, G, qC, V, R, fAtom, su, phase_of, dc,
                      platt_file, displacement_file, material_file, phase_file, region_file, output_file):
    """
    Write the inputs of the multi-GPU program qpop_diffraction_mgpu.

    Args:
        full:    (nx, ny, nz) of the FerroDyn grid
        limits:  (x0, x1, y0, y1, z0, z1), 0-indexed inclusive crop
        spacing: (dx, dy, dz) in nm
        G, qC:   reciprocal lattice point and strain-corrected centre, nm^-1
        V:       unit-cell volume, nm^3
        R:       (nAtom, 3) atom positions in nm
        fAtom:   (nPhase, nAtom) complex form factors
        su:      (3, nStruc, nPhase, nAtom) ionic shifts in nm per C/m^2
        phase_of: phase index of each material, -1 for non-scattering ones
        dc:      c_m(T) / aC_z - 1 of each material
        *_file:  inputs and the output .h5; displacement_file, phase_file and region_file may be None
    """
    with h5py.File(path, "w") as fh:
        fh["grid"] = np.asarray(full, dtype=np.int64)
        fh["crop"] = np.asarray(limits, dtype=np.int64)
        fh["spacing"] = np.asarray(spacing, dtype=np.float64)
        fh["G"] = np.asarray(G, dtype=np.float64)
        fh["qC"] = np.asarray(qC, dtype=np.float64)
        fh["V"] = np.array([V], dtype=np.float64)
        fh["R"] = np.asarray(R, dtype=np.float64)
        fh["fAtom_re"] = np.ascontiguousarray(np.real(fAtom), dtype=np.float64)
        fh["fAtom_im"] = np.ascontiguousarray(np.imag(fAtom), dtype=np.float64)
        fh["su"] = np.asarray(su, dtype=np.float64)
        fh["phase_of"] = np.asarray(phase_of, dtype=np.int32)
        fh["dc"] = np.asarray(dc, dtype=np.float64)
        fh["u_scale"] = np.array([1.0 / L0])
        with h5py.File(platt_file, "r") as fp:
            platt_prefix = next(p for p in ("P_Lattice", "PLatt") if f"{p}_x" in fp)
        for key, value in (("platt_file", platt_file), ("platt_prefix", platt_prefix),
                           ("displacement_file", displacement_file), ("material_file", material_file),
                           ("phase_file", phase_file), ("region_file", region_file), ("output_file", output_file)):
            if value is not None:
                fh.attrs[key] = np.bytes_(str(value))
    print(f"Writing multi-GPU parameters → {path}")


def load_diffraction_hdf5(filepath) -> dict:
    """
    Load the output of :func:`write_diffraction_output`.

    Returns a dict with ``intensity`` and ``log10_intensity`` shaped
    (1, nx, ny, nz), the 1-D axes ``q_x``, ``q_y``, ``q_z``, ``q_vector``
    (3, nx, ny, nz) broadcast from them, ``dims`` and ``q_units``.
    """
    with h5py.File(filepath, "r") as fh:
        I = fh["intensity"][:].transpose(2, 1, 0)
        q = [fh[f"q_{c}"][:] for c in "xyz"]
        q_units = str(fh.attrs["q_units"])
    shape = I.shape
    return {
        "dims": shape,
        "intensity": I[np.newaxis],
        "log10_intensity": np.log10(I)[np.newaxis],
        "q_x": q[0], "q_y": q[1], "q_z": q[2],
        "q_vector": np.stack([np.broadcast_to(qd.reshape([-1 if d == i else 1 for i in range(3)]), shape)
                              for d, qd in enumerate(q)]),
        "q_units": q_units,
    }
