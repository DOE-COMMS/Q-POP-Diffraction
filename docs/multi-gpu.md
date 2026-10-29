# Multi-GPU runs

With the multi-GPU build, grids that do not fit in one GPU are run in two steps. `qpop-diffraction --mgpu` does the setup in Python. `qpop_diffraction_mgpu`, a compiled C++/CUDA MPI program (`bin/qpop_diffraction_mgpu`), then does every per-cell step on the GPUs, with the FFTs distributed by HeFFTe. The result is the same file the standard build writes. See [Installation](installation.md) for how to compile the program.

## Workflow

```bash
qpop-diffraction -i input.json -s system_setting.json -a parameter.atom.in --mgpu
srun -N 2 --ntasks-per-node=4 --gpus-per-node=4 --gpu-bind=none \
    /path/to/Q-POP-Diffraction/bin/qpop_diffraction_mgpu diffraction_15000_params.h5
```

The `srun` line is for Slurm with 4 GPUs per node; with another launcher, start one MPI rank per GPU (e.g. `mpirun -n 8 ...`).

The first step is cheap and never loads the fields, so it can run on a login node or at the start of the job. It does three things:

- It reads `input.json`, `system_setting.json` and `parameter.atom.in`, works out the temperature and strain reference, and finds G and qC.
- It writes `{output_file}.h5`, holding the q axes and attributes but not `/intensity`, and `{output_file}.xdmf`.
- It writes `{output_file}_params.h5` for `qpop_diffraction_mgpu` and prints an `srun` line sized to the grid.

`qpop_diffraction_mgpu` then runs on one MPI rank per GPU:

- Each rank reads its slab of P, u, the material field and, if given, the phase and region files, with collective parallel HDF5.
- On its GPU it builds the u_z reference ramp and the phase of every cell, then for every atom four source fields, four HeFFTe forward transforms, and the amplitude update.
- It computes the intensity, applies the fftshift, and writes its slab of `/intensity` into `{output_file}.h5`.

`examples/run_mgpu.sbatch` is a Slurm job script that does both steps; adapt its account, partition and environment to your system.

## Decomposition

The grid is split into slabs of y rows, one per rank, as evenly as possible. Every rank holds all x and all z of its rows, in real and in Fourier space. That matches FerroDyn's `(nz, ny, nx)` files, so reads and writes are a few large contiguous blocks. There can't be more ranks than rows in y.

## Sizing

`qpop_diffraction_mgpu` needs about 130 bytes of GPU memory per cell:

- P and u: 48 bytes;
- the amplitude and the FFT buffer: 32 bytes;
- HeFFTe's workspace: about 48 bytes;
- the phase index: 1 byte.

A `region_file` adds 8 bytes per cell, and a `phase_file` adds 8 bytes per phase.

- 4000×4000×68 (1.09e9 cells) needs about 140 GB: 8 GPUs of 40 GB (two nodes of 4), or 4 GPUs of 80 GB (one node).
- 512×512×104 (2.7e7 cells) fits one GPU; the standard build is simpler there.

The `--mgpu` step prints its estimate. Rank 0 prints the largest device memory actually used on any rank.

## Launch

- One rank per GPU, with `--gpu-bind=none`. Each rank picks the GPU matching its rank on the node.
- GPU-aware MPI must be enabled; for Cray MPICH, `export MPICH_GPU_SUPPORT_ENABLED=1`.
- Run `qpop_diffraction_mgpu` from anywhere: the params file holds absolute paths to all inputs and to the output.

## Performance

Rank 0 prints the time spent reading, computing the amplitude (and how much of that was FFTs), and writing. At this size the FFTs themselves are a small part. Measured on NVIDIA A100 nodes with 4 GPUs each, one 1024³ complex transform takes about 0.3 s with HeFFTe on 8 GPUs, and a 5-atom cell needs 20 of them. Reading tens of GB of float64 fields usually dominates. Off-node transposes are bandwidth-bound, so prefer the fewest GPUs that fit, e.g. 80 GB nodes over twice as many 40 GB GPUs.

## Troubleshooting

- "GPU_SUPPORT_ENABLED is requested, but GTL library is not linked" (Cray MPICH): the program was built without Cray MPICH's GPU transport layer. Rebuild with the Cray MPICH module loaded, so the Makefile detects it.
- Jobs fail at launch with "Permission denied" (exit code 126): the binary was overwritten while they were queued. Don't rebuild with jobs pending.
- "has shape (...), expected (nz, ny, nx)": the field files are in the pre-2026-05-26 layout; convert them with `qpop_monitor.migrate_hdf5_layout`.
- "More ranks than rows in y": use fewer ranks, or don't crop y so tightly.
- Out of device memory: add GPUs or use 80 GB nodes; see Sizing.

## Validation

`tests/test_mgpu.py` runs `qpop_diffraction_mgpu` and the Python path on the same small runs and requires the same intensity:

- bulk, film on substrate, cropped, and phase/region-weighted cases;
- 3 ranks over 32 rows, an uneven split that also exercises the wrapped fftshift.

Run it inside a GPU allocation after building ([Testing](testing.md)).
