#!/usr/bin/env python3
"""
Q-POP Diffraction Simulation for Q-POP-FerroDyn output

Usage
-----
    qpop-diffraction [--input input.json]
                     [--system system_setting.json]
                     [--atom parameter.atom.in]
                     [--mgpu]

Without --mgpu the whole calculation runs here, on one GPU (CuPy) or the CPU.
With --mgpu this only prepares a run of the multi-GPU program
(qpop_diffraction_mgpu, from the multi-GPU build): it writes {output_file}_params.h5 and the
output files without /intensity, and prints the srun line that fills it in.

Required input files
--------------------
system_setting.json   FerroDyn input of the run: grid (Mesh), materials,
                      temperature and misfit/applied strain.  The run's
                      materials/ directory and data/materialField_0.h5 are
                      read from beside it.
parameter.atom.in     Reciprocal lattice point, lattice vectors, atomic
                      structure, form factors, displacement coefficients.
                      One phase block per Is_FE material.
input.json            Paths to the P_Lattice / Displacement HDF5 files,
                      optional phase fractions and scattering weight, and
                      optional x/y/z slice limits.  See json_input.py.

Output
------
{output_file}.h5      /intensity in (nz, ny, nx) order, as FerroDyn stores
                      fields, and the 1-D q axes /q_x, /q_y, /q_z in nm⁻¹.
{output_file}.xdmf    XDMF 3.0 sidecar for direct ParaView 3-D
                      visualisation (open this file in ParaView).
{output_file}_params.h5  Inputs of the multi-GPU program (--mgpu only).
"""

import argparse
import math
import shutil
import sys
from pathlib import Path

import h5py
import numpy as np

from .constants import OS0
from .size_context import SizeContext
from .atom_list import AtomList
from .diffraction import diffraction, reciprocal_centre, q_axes
from .io_utils import (
    load_platt_hdf5,
    load_displacement_hdf5,
    load_material_field,
    load_phase_hdf5,
    load_region_hdf5,
    crop_slices,
    write_diffraction_output,
    write_diffraction_axes,
    write_mgpu_params,
)
from .json_input import load_input_json


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Q-POP Diffraction Simulation for Q-POP-FerroDyn output",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input", "-i",
        default="input.json",
        metavar="JSON",
        help="Path to input.json (P_Lattice/Displacement paths, slice limits, output name).",
    )
    parser.add_argument(
        "--system", "-s",
        default="system_setting.json",
        metavar="JSON",
        help="Path to the FerroDyn system_setting.json of the run.",
    )
    parser.add_argument(
        "--atom", "-a",
        default="parameter.atom.in",
        metavar="FILE",
        help="Path to parameter.atom.in.",
    )
    parser.add_argument(
        "--mgpu",
        action="store_true",
        help="Prepare a run of the multi-GPU program qpop_diffraction_mgpu instead of computing here (for grids larger than one GPU).",
    )
    return parser.parse_args()


# Device bytes per cell in qpop_diffraction_mgpu: P and u (6 doubles), phase index, amp,
# the FFT buffer and HeFFTe's workspace (about three complex arrays); plus one
# double for region_file and one per phase for phase_file.
MGPU_BYTES_PER_CELL = 6 * 8 + 1 + 2 * 16 + 3 * 16
MGPU_PROGRAM = Path(shutil.which("qpop_diffraction_mgpu")
                    or Path(__file__).resolve().parents[2] / "bin" / "qpop_diffraction_mgpu")


def _prepare_mgpu(args, cfg, params, fe_ids, atoms, lims):
    """Write the inputs of qpop_diffraction_mgpu and the output skeleton; print the srun line."""
    full = (params.nx, params.ny, params.nz)
    with h5py.File(cfg["platt_file"], "r") as fh:
        step = int(fh.attrs["Timestep"])

    print(f"\nTimestep {step}:")
    dc = params.reference_strain(np.diag(atoms.aC), params.temperature(step))

    print("\nApplying domain slice:")
    limits = params.apply_slice(*lims)
    params.setup()

    G, qC = reciprocal_centre(params, atoms)
    phase_of = np.full(len(params.names), -1)
    phase_of[fe_ids] = np.arange(len(fe_ids))

    out = Path(f"{cfg['output_file']}.h5").resolve()
    params_file = Path(f"{cfg['output_file']}_params.h5").resolve()
    write_diffraction_axes(cfg["output_file"], *q_axes(params, qC), step)
    write_mgpu_params(
        params_file, full, limits, (params.dx, params.dy, params.dz),
        G, qC, np.linalg.det(atoms.aC), atoms.xAtom @ atoms.aC, atoms.fAtom,
        np.stack([atoms.su1, atoms.su2, atoms.su3]) / OS0, phase_of, dc,
        Path(cfg["platt_file"]).resolve(),
        Path(cfg["displacement_file"]).resolve() if cfg["displacement_file"] is not None else None,
        (Path(args.system).parent / "data" / "materialField_0.h5").resolve(),
        Path(cfg["phase_file"]).resolve() if cfg["phase_file"] is not None else None,
        Path(cfg["region_file"]).resolve() if cfg["region_file"] is not None else None,
        out,
    )

    per_cell = MGPU_BYTES_PER_CELL + 8 * (cfg["region_file"] is not None) \
        + 8 * params.nPhase * (cfg["phase_file"] is not None)
    gb = params.n * per_cell / 1e9
    gpus = max(1, math.ceil(gb / 36.0))
    nodes = math.ceil(gpus / 4)
    print(f"\nqpop_diffraction_mgpu needs about {gb:.0f} GB of GPU memory in total, "
          f"e.g. {gpus} GPU(s) of 40 GB or {max(1, math.ceil(gb / 76.0))} of 80 GB.")
    print("Run it with one MPI rank per GPU, e.g. on nodes with 4 GPUs of 40 GB:")
    print(f"  srun -N {nodes} --ntasks-per-node=4 --gpus-per-node=4 --gpu-bind=none {MGPU_PROGRAM} {params_file}")
    return 0


def main() -> int:
    """Main simulation function."""
    args = _parse_args()

    # ── 1. Parse input.json ──────────────────────────────────────────────────
    cfg = load_input_json(args.input)
    lims = (cfg["x_limits"], cfg["y_limits"], cfg["z_limits"])

    # ── 2. Read system_setting.json (full-domain dimensions) ─────────────────
    params = SizeContext()
    fe_ids = params.read_json(args.system)
    full = (params.nx, params.ny, params.nz)
    dz = params.lz / params.nz
    sl = crop_slices(*full, *lims)

    # ── 3. Read parameter.atom.in ────────────────────────────────────────────
    atoms = AtomList()
    if not atoms.read(args.atom, params.nPhase, params.nStruc):
        return 1
    if not atoms.setup():
        return 1

    if args.mgpu:
        return _prepare_mgpu(args, cfg, params, fe_ids, atoms, lims)

    # ── 4. Load P_Lattice → oStruc, Displacement → u (cropped, x,y,z order) ─
    print(f"Loading P_Lattice file: {cfg['platt_file']}")
    oStruc = load_platt_hdf5(cfg["platt_file"], *full, *lims) / OS0
    with h5py.File(cfg["platt_file"], "r") as fh:
        step = int(fh.attrs["Timestep"])

    u = None
    if cfg["displacement_file"] is not None:
        print(f"Loading Displacement file: {cfg['displacement_file']}")
        u = load_displacement_hdf5(cfg["displacement_file"], *full, *lims)
    else:
        print("No displacement_file provided — using zero displacement field.")

    # ── 5. Phases and region: non-FE cells of materialField do not scatter ───
    #   phase_file replaces the one-hot phases, region_file weights the cells.
    mat_file = Path(args.system).parent / "data" / "materialField_0.h5"
    print(f"Loading material field: {mat_file}")
    mat = load_material_field(mat_file, *full)
    oPhase = np.stack([mat[sl] == m for m in fe_ids]).astype(np.float64)
    region = oPhase.sum(axis=0).ravel()
    if cfg["phase_file"] is not None:
        print(f"Loading phase fractions: {cfg['phase_file']}")
        oPhase = load_phase_hdf5(cfg["phase_file"], params.nPhase, *full, *lims)
    oPhase = oPhase.ravel()
    if cfg["region_file"] is not None:
        print(f"Loading scattering weight: {cfg['region_file']}")
        region *= load_region_hdf5(cfg["region_file"], *full, *lims).ravel()

    # ── 6. Strains referred to aC at the snapshot temperature ────────────────
    #   FerroDyn measures strain from each material's own lattice.  The
    #   out-of-plane mismatch c_m/aC_z - 1 becomes a u_z ramp summed up the
    #   full stack (cell centres), so each layer diffracts at its own c.
    print(f"\nTimestep {step}:")
    dc = params.reference_strain(np.diag(atoms.aC), params.temperature(step))[mat]
    uz_ref = ((np.cumsum(dc, axis=2) - 0.5 * dc) * dz)[sl].ravel()

    # ── 7. Apply slice to SizeContext ────────────────────────────────────────
    print("\nApplying domain slice:")
    params.apply_slice(*lims)
    params.setup()

    if u is None:
        u = np.zeros(3 * params.n)
    u[2 * params.n:] += uz_ref

    # ── 8. Diffraction and output ────────────────────────────────────────────
    print("\nComputing diffraction")
    I, qx, qy, qz = diffraction(params, atoms, oPhase, oStruc, u, region)
    write_diffraction_output(cfg["output_file"], I, qx, qy, qz, step)

    print("\nDiffraction simulation completed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
