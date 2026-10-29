# Testing

```bash
pip install -r requirements.txt     # includes pytest
pytest tests/
```

The tests build small synthetic runs in FerroDyn layout, run `qpop-diffraction` on them, and compare the intensity with answers derived independently of the code. Every grid is non-cubic, so a swapped axis fails. They run on the CPU in seconds, or on the GPU when CuPy is installed.

## Single-domain PbTiO3: test_single_domain.py

Uniform P = 0.75 C/m² along each of the six tetragonal variants c±, a1± and a2±.

- Structure factor: at the 002, 200, 020 and 111 reflections the peak sits at G, and I(G) = |Σ_a f_a exp(-i G·(R_a + d_a))|² / V². Here d_a is the P-induced shift, and the test computes the sum directly in numpy.
- Cubic symmetry: c at 002, a1 at 200 and a2 at 020 give the same intensity. Mixing up P components or axes breaks this.
- A shift perpendicular to G is invisible: an a domain at 002 scatters like P = 0.
- Friedel's law is broken: P and −P differ at 002, because the form factors have imaginary parts.
- Spontaneous strain: a strain e along P moves only the reflections along P, to G/(1 + e). The test uses e = 2%, small enough for the first-order model (see [Method](method.md)).

## a1/a2 domains in a film: test_a1a2_domains.py

Stripes of a1 and a2 domains with (110) walls and period p cells (p = 8 and 16). They sit in films t = 8 and 16 nm thick on a non-scattering substrate with air above, and are seen at the 200 reflection.

- Film thickness: along q_z through G the intensity follows the slab Laue function sinc²(k_z T/2) × |F(G + k_z)|² / |F(G)|², with zeros at k_z = 2πm/T. This is exact.
- Domain width: the first satellites sit at k = 2π/(p dx) (1, 1, 0) with the square-wave intensity S²|c₁|²|F_a1 − F_a2|² / |F|², to 1%. The difference is the first-order expansion in k·d.
- Orientation: nothing appears along (1, −1, 0).
- Equal widths: the second order is absent.

## Planted signals: test_planted.py

- uniform P peaks at G;
- a u_z = e·z ramp moves the 002 peak to G(1 − e), which catches a missing m → nm conversion;
- 180° stripes along x give satellites on q_x only;
- a substrate gives the expected in-plane strain relative to `aC`, and is not a phase;
- a film whose lattice differs from `aC` peaks at 4π/c_film;
- the temperature follows FerroDyn's ramp;
- rotated runs are rejected;
- `phase_file` and `region_file` weight the cells.

## Examples: test_examples.py

Runs every case in `examples/`. The disordered polar vortex is also compared with `reference_legacy.h5`, the intensity of a legacy-format run on the same inputs.

## GPU: test_gpu.py

Checks CuPy's FFT against NumPy. Skipped without CuPy or a GPU.

## Multi-GPU build: test_mgpu.py

Runs `qpop_diffraction_mgpu` and the Python path on the same planted runs and requires the same intensity:

- bulk, film on substrate, cropped, and phase/region-weighted cases;
- 3 ranks over 32 rows, an uneven split.

It needs the multi-GPU build and a launch command inside a GPU allocation, and is skipped otherwise:

```bash
QPOP_MGPU_LAUNCH="srun -n 3 --gpus-per-node=4 --gpu-bind=none" pytest tests/test_mgpu.py
```
