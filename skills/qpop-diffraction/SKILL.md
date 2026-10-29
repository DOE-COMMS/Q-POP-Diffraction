---
name: qpop-diffraction
description: Compute, visualise and interpret kinematic X-ray diffraction (reciprocal-space volumes, RSMs, rocking curves, satellites, thickness fringes) from Q-POP-FerroDyn phase-field output with Q-POP-Diffraction. Use whenever someone asks for a diffraction pattern, XRD, an RSM or a peak position from a FerroDyn run; when writing input.json or parameter.atom.in; when choosing between the standard and multi-GPU builds; or when reading diffraction.h5.
---

# Q-POP-Diffraction

Q-POP-Diffraction turns one FerroDyn snapshot into the 3-D kinematic X-ray intensity around one reciprocal lattice point. The inputs are P (`P_Lattice_<step>.h5`), u (`Displacement_<step>.h5`), the run's `system_setting.json` and a crystal description (`parameter.atom.in`). The output is `{output_file}.h5`: `/intensity` shaped `(nz, ny, nx)` and the q axes `/q_x`, `/q_y`, `/q_z` in nm⁻¹. Repository docs: `docs/inputs.md`, `docs/method.md`, `docs/outputs.md`.

## 1. Set up a run

Work in the FerroDyn run directory, which holds `system_setting.json`, `materials/` and `data/materialField_0.h5`.

- **input.json:**
  - `platt_file` (required);
  - `displacement_file`: the same step as P, which matters because the elastic save interval can differ;
  - `x/y/z_limits`, `phase_file`, `region_file`, `output_file`.
  - Start from `examples/input.json`.
- **parameter.atom.in:** read as a stream of numbers, with `!` starting a comment:
  - `q00` (3 numbers), then `aC` (9, rows a, b, c in nm), then `nAtom`, then the reduced coordinates;
  - then one block per `Is_FE` material, in `MaterialNames` order. Each block has the phase ID, `nAtom` form factors (real, imaginary), and `nAtom` × 3 lines of P→shift vectors (m per C/m², one line each for P_x, P_y, P_z).
  - Every phase block must have all 3 × nAtom shift lines, zeros if P does not move the atoms. A missing line silently shifts every later number: the next phase ID is read as a coefficient.
  - `q00` picks the reflection: the nearest reciprocal lattice point of `aC` is used, so e.g. `0 0 31.416` with a = 0.4 nm gives 002.
  - Choose `aC` close to the film's average lattice (see section 4).
- **Check before running:**
  - field datasets are `(nz, ny, nx)` over the full `Mesh`;
  - `Rotation.Enabled` is false (rotated runs are rejected);
  - the number of `Is_FE` materials equals the number of phase blocks.

## 2. Choose the build and run

- **Standard build** (Python; CuPy GPU or NumPy CPU). Needs about 250 bytes per cell on one GPU, so about 150 M cells on 40 GB. Run:
  ```bash
  qpop-diffraction -i input.json -s system_setting.json -a parameter.atom.in
  ```
  `QPOP_DIFFRACTION_DEVICE=cpu` forces NumPy. Its memory is large (about 400 GB for 1.1e9 cells) and it is slow.
- **Multi-GPU build** (`bin/qpop_diffraction_mgpu`: MPI, HeFFTe, parallel HDF5). Needs about 130 bytes per cell in total, spread over the GPUs:
  - `qpop-diffraction ... --mgpu` writes `{output_file}_params.h5` and the output without `/intensity`, and prints an `srun` line. It never loads the fields.
  - Then launch `qpop_diffraction_mgpu {output_file}_params.h5` with one MPI rank per GPU and `--gpu-bind=none`.
  - There can't be more ranks than y rows.
  - Rank 0 prints the read, amplitude, FFT and write times and the peak device memory.
  - A 4000×4000×68 grid needs about 140 GB: 8 × 40 GB or 4 × 80 GB.
- Ask the user which machine and queue to use. Do not rebuild `bin/qpop_diffraction_mgpu` while jobs that use it are queued: they fail at launch with "Permission denied".

## 3. Read and show the output

- `qpop-diffraction-rsm -d diffraction.h5 -o maps` writes the qx–qy, qx–qz and qy–qz projections (`.png` and `.npz`), summed over the third axis.
- In Python, `h5py`: `I = f["intensity"][:]` is `(nz, ny, nx)`; pair `I[k, j, i]` with `(q_x[i], q_y[j], q_z[k])`.
  - `qpop_diffraction.load_diffraction_hdf5` instead returns `(x, y, z)` order, with a leading length-1 axis, plus a broadcast `q_vector`.
  - Mixing the two orders silently swaps axes.
- ParaView: open `diffraction.xdmf`, then add a Calculator with `log10(Intensity)`, and a Slice or Contour.
- Always plot log10 intensity. Values are floored at 1e-50, so a map full of −50 means there is no signal there, not that something failed.
- To compare with an experiment, render the same observable: same reflection, same projection, axis ranges and orientation. Do this before claiming agreement.

## 4. Interpret

- q = qC + k. qC is G shifted by the in-plane average strain: the substrate's lattice for a film, or the misfit. Out-of-plane, each material's c(T)/aC_z − 1 enters as a u_z ramp summed up the stack, so each layer's 00l peak sits at its own c.
- Temperature comes from the P file's `Timestep` attribute and FerroDyn's `Temperature` ramp. A run restarted with a different `Initial` gets the wrong T.
- Only cells of `Is_FE` materials scatter, so there is no substrate peak in a film.
- Known answers (all checked by the tests in `tests/`):
  - single domain: I(G) = |Σ f e^{−iG·(R+δ)}|²/V²;
  - a shift perpendicular to G is invisible;
  - P and −P differ when the form factors have imaginary parts;
  - a film of thickness T: Laue fringes sinc²(k_z T/2), zeros at 2πm/T;
  - stripe domains of period Λ: satellites at 2π/Λ along the stripe wavevector; with equal widths, the even orders vanish.
- Resolution is Δq = 2π/L per axis, and the range is ±π/dx. A small box cannot resolve close satellites.
- The model is first order in k·u. A strain e moves a peak to G(1 − e), not G/(1 + e); the error is G·e². Keep `aC` close to the lattice so that e is small.
- It is kinematic: no dynamical effects, absorption or instrument resolution.

## 5. Traps

- Never crop z to the film (`z_limits`). The FFT is periodic, so the film's top wraps onto its bottom and the u_z ramp makes a phase jump there. That smears and shifts the 00l peaks (one case moved from 47.86 to 48.07 nm⁻¹). Keep the substrate and air as padding.
- Field files written before 2026-05-26 are `(nx, ny, nz)` and are rejected. Convert them with `qpop_monitor.migrate_hdf5_layout`; do not transpose by hand.
- `phaseFra_file` (text) is gone. Use `phase_file` (HDF5 datasets `phase_1..N`) and `region_file` (`region`), laid out like the fields. Any float type and compression are fine.
- `data/materialField_0.h5` is read from beside `system_setting.json`. If the run's `materials/` or `data/` are symlinks, check that they resolve.
- The `Bulk-stress-free` macrostrain is not stored in the input, so it is not included.
- Intensities are in arbitrary units (|F|²/V² per cell, averaged over the box), not normalised to a beam. Compare shapes and positions, or ratios within one run.

## 6. Verify

- After installing: `pytest tests/` (planted signals, single-domain PbTiO3 variants, a1/a2 films, examples). With the multi-GPU build, also `QPOP_MGPU_LAUNCH="srun -n 3 ..." pytest tests/test_mgpu.py` inside a GPU allocation.
- For a new material or setup, first run a uniform single domain and check that the peak sits at 2π·h/d for its lattice. Only then run the real snapshot.
- Report what ran and what did not. A pattern from an untested setup is not a result.
