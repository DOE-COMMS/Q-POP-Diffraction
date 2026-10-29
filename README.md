# Q-POP-Diffraction

Kinematic X-ray diffraction patterns computed directly from [Q-POP-FerroDyn](https://github.com/DOE-COMMS/Q-POP-FerroDyn) phase-field output.

Given a FerroDyn snapshot (polarization `P_Lattice_<step>.h5`, displacement `Displacement_<step>.h5`) and a description of the crystal (`parameter.atom.in`), it computes the 3-D intensity around a chosen reciprocal lattice point. The output is a reciprocal-space volume you can open in ParaView or slice into reciprocal-space maps. The grid, materials, temperature and substrate strain are read from the run's own `system_setting.json`.

## Builds

There are two builds; choose one when you install.

| Build | Computes on | For | Needs |
|---|---|---|---|
| Standard | one GPU (CuPy) or the CPU (NumPy) | grids that fit one GPU, e.g. 512×512×104 in ~12 s on an A100 | Python, numpy, h5py; CuPy for the GPU |
| Multi-GPU | many GPUs, with the MPI program `qpop_diffraction_mgpu` | larger grids, e.g. 4000×4000×68 in 1.5 min on 4 A100s | the standard build plus CUDA, GPU-aware MPI, [HeFFTe](https://github.com/icl-utk-edu/heffte) 2.4.1 and parallel HDF5 |

Both read the same inputs and write the same result.

### Standard build

```bash
git clone https://github.com/DOE-COMMS/Q-POP-Diffraction.git
cd Q-POP-Diffraction
pip install -r requirements.txt   # numpy, h5py, matplotlib, pytest
pip install -e .                  # the qpop_diffraction package and its commands
pip install cupy                  # optional: GPU support
```

### Multi-GPU build

Do the standard build first. The multi-GPU build also needs:

- the CUDA toolkit;
- an MPI library with CUDA-aware support;
- parallel HDF5;
- [HeFFTe](https://github.com/icl-utk-edu/heffte) 2.4.1, built with the cuFFT backend and GPU-aware MPI.

HeFFTe is not part of this repository and must be installed first; [docs/installation.md](docs/installation.md) shows how. Then compile the multi-GPU program from the repository root:

```bash
make HEFFTE_DIR=/path/to/heffte HDF5_DIR=/path/to/parallel/hdf5 MPI_DIR=/path/to/mpi
```

This compiles the C++/CUDA sources in `src/heffte/` with nvcc. It links the result against HeFFTe, cuFFT, parallel HDF5 and MPI, and writes the executable `bin/qpop_diffraction_mgpu`. `qpop-diffraction --mgpu` looks for the program on `PATH`, then in `bin/`. The build is tested with CUDA 12.9 on A100 GPUs, GCC 14, Cray MPICH 9.0.1, HeFFTe 2.4.1 and HDF5 1.14.3.

## Quick start

```bash
cd examples/4_polar_vortex_PbTiO3/disordered
qpop-diffraction                  # reads input.json, system_setting.json, parameter.atom.in
qpop-diffraction-rsm -d diffraction.h5
```

This writes `diffraction.h5` and `diffraction.xdmf`, which you can open in ParaView, and three 2-D maps in `reciprocal_space_maps/`.

## Using it on a FerroDyn run

In the run directory, write an `input.json` (template: [examples/input.json](examples/input.json)):

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

Then add a `parameter.atom.in` for your crystal, with one phase block per ferroelectric material. Run:

```bash
qpop-diffraction -i input.json -s system_setting.json -a parameter.atom.in
```

With the multi-GPU build, run larger grids in two steps. `--mgpu` prepares the run without loading the fields, and prints an `srun` line sized to your grid. That line launches `qpop_diffraction_mgpu`, one MPI rank per GPU:

```bash
qpop-diffraction -i input.json -s system_setting.json -a parameter.atom.in --mgpu
srun -N 2 --ntasks-per-node=4 --gpus-per-node=4 --gpu-bind=none \
    /path/to/Q-POP-Diffraction/bin/qpop_diffraction_mgpu diffraction_15000_params.h5
```

[examples/run_mgpu.sbatch](examples/run_mgpu.sbatch) is a Slurm job script to adapt to your system.

## Visualizing the output

A run writes `{output_file}.h5` with the intensity `/intensity`, shaped `(nz, ny, nx)`, and the q axes `/q_x`, `/q_y`, `/q_z` in nm⁻¹. It also writes `{output_file}.xdmf`, which describes the same data for ParaView. The intensity spans many decades, so always look at it on a log scale.

**Reciprocal-space maps.** Project the volume onto the qx–qy, qx–qz and qy–qz planes:

```bash
qpop-diffraction-rsm -d diffraction.h5 -o maps
```

This writes `maps/qx_qy.png`, `maps/qx_qz.png` and `maps/qy_qz.png` (log10 intensity), plus `.npz` files with the maps and their axes for your own plots. `--aggregation sum` (the default) integrates over the third axis like an area detector; `mean` averages instead.

**ParaView (3-D).** Open `diffraction.xdmf` and click Apply. Add a Calculator filter with the expression `log10(Intensity)`, then:

- use Slice for a plane through the peak, or
- use Contour for isosurfaces of the peak and its satellites.

The axes are q in nm⁻¹.

**Python.** A q_x–q_z cut through the strongest peak:

```python
import h5py
import matplotlib.pyplot as plt
import numpy as np

with h5py.File("diffraction.h5", "r") as f:
    I = f["intensity"][:]                                   # (nz, ny, nx)
    qx, qy, qz = f["q_x"][:], f["q_y"][:], f["q_z"][:]

j = np.unravel_index(np.argmax(I), I.shape)[1]              # q_y index of the peak
plt.pcolormesh(qx, qz, np.log10(I[:, j, :]), shading="auto", cmap="jet")
plt.xlabel("$q_x$ (nm$^{-1}$)")
plt.ylabel("$q_z$ (nm$^{-1}$)")
plt.colorbar(label="log$_{10}$ I")
plt.savefig("qx_qz.png", dpi=200)
```

`qpop_diffraction.load_diffraction_hdf5("diffraction.h5")` returns the same data in `(x, y, z)` order, with a broadcast `q_vector`. See [docs/outputs.md](docs/outputs.md).

## For AI agents

[skills/qpop-diffraction/SKILL.md](skills/qpop-diffraction/SKILL.md) is a skill for agentic AI tools, such as Q-POP-Manager, Claude Code, Codex and others that read `SKILL.md` skills. It covers setting up a run, choosing a build, reading the output, interpreting peaks and the known traps. Link it into your agent's skill directory, for example:

```bash
ln -s "$PWD/skills/qpop-diffraction" ~/.claude/skills/qpop-diffraction
```

[AGENTS.md](AGENTS.md) points agents that read that file to the same skill.

## Documentation

- [Installation](docs/installation.md): the two builds, GPU support, installing HeFFTe, compiling `qpop_diffraction_mgpu`
- [Inputs](docs/inputs.md): `input.json`, `parameter.atom.in`, the FerroDyn files, phase and region files
- [Method](docs/method.md): the diffraction model, the strain reference, temperature, limits
- [Multi-GPU runs](docs/multi-gpu.md): workflow, sizing, launch, performance
- [Outputs](docs/outputs.md): the HDF5/XDMF result and reciprocal-space maps
- [Testing](docs/testing.md): the test suite and the known answers it checks
- [Diffraction flow](docs/diffraction-flow.md): pipeline flowchart and the q-pop-monitor contract

## Repository layout

```
src/python/    the qpop_diffraction Python package: inputs, single-GPU/CPU calculation, outputs, reciprocal-space maps
src/heffte/    C++/CUDA source of the multi-GPU program qpop_diffraction_mgpu (MPI, HeFFTe, parallel HDF5)
Makefile       builds bin/qpop_diffraction_mgpu (multi-GPU build)
examples/      ready-to-run examples, input.json template, Slurm job script for multi-GPU runs
tests/         pytest suite
docs/          documentation (published to the project website)
```

## Examples

[examples/](examples/README.md) holds four ready-to-run cases: a two-phase mixture, a two-variant twin, a three-variant hexagonal twin, and a polar vortex array in PbTiO3.

## Tests

```bash
pytest tests/
```

The suite checks the code against analytic answers: structure factors of single-domain PbTiO3 variants, thickness fringes and domain satellites of a1/a2 films, strain shifts and planted modulations. It also covers the examples, and `qpop_diffraction_mgpu` against the Python path. See [docs/testing.md](docs/testing.md).

## Acknowledgement

This development of the software was supported as part of the Computational Materials Sciences Program funded by the U.S. Department of Energy, Office of Science, Basic Energy Sciences, under Award No. DE-SC0020145.
