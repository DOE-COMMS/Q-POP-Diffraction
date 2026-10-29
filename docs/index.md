# Q-POP-Diffraction

Q-POP-Diffraction computes kinematic X-ray diffraction from Q-POP-FerroDyn phase-field output. It takes one snapshot of the polarization and displacement fields plus a description of the crystal. From these it builds the complex amplitude scattered by every cell, Fourier-transforms it, and writes the 3-D intensity around a chosen reciprocal lattice point.

The intensity is written in reciprocal space (nm⁻¹) on the FFT grid of the simulation, centred on the strain-corrected Bragg position. From there you can view it in ParaView, project it into reciprocal-space maps, or compare it with measured RSMs.

## Two ways to run

- `qpop-diffraction`: the whole calculation in Python, on one GPU with CuPy or on the CPU with NumPy. Use it whenever the grid fits in one GPU; 512×512×104 takes about 12 s on an A100.
- `qpop-diffraction --mgpu` followed by `qpop_diffraction_mgpu` (multi-GPU build): Python prepares the run, and a compiled C++/CUDA MPI program computes it on many GPUs, with the FFTs done by HeFFTe. Use it for grids larger than one GPU, e.g. 4000×4000×68.

Both read the same inputs and write the same output file.

## Pages

- [Installation](installation.md): the standard and multi-GPU builds, installing HeFFTe, compiling `qpop_diffraction_mgpu`
- [Inputs](inputs.md): every input file and key, including the `parameter.atom.in` format
- [Method](method.md): the diffraction model, the strain reference, temperature, limits and approximations
- [Multi-GPU runs](multi-gpu.md): workflow, sizing, launch, performance, troubleshooting
- [Outputs](outputs.md): the result file and reciprocal-space maps
- [Testing](testing.md): the test suite and the analytic answers it checks
- [Diffraction flow](diffraction-flow.md): pipeline flowchart and the contract for q-pop-monitor

## Quick start

```bash
pip install -r requirements.txt
pip install -e .
cd examples/4_polar_vortex_PbTiO3/disordered
qpop-diffraction
qpop-diffraction-rsm -d diffraction.h5
```
