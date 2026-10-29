# Installation

## Choose a build

Q-POP-Diffraction comes in two builds. Decide which one you need before installing.

- Standard build: the Python package only. It computes on one GPU (with CuPy) or on the CPU. Choose it when your grids fit in the memory of one GPU. It needs roughly 250 bytes per cell, so about 150 million cells on a 40 GB GPU; 512×512×104 (27 million cells) takes about 12 s on an NVIDIA A100.
- Multi-GPU build: the standard build plus `qpop_diffraction_mgpu`, a compiled C++/CUDA MPI program that spreads the calculation over many GPUs. Choose it for grids larger than one GPU, e.g. 4000×4000×68. It needs CUDA, CUDA-aware MPI, HeFFTe and parallel HDF5.

Both builds read the same inputs and write the same output. The multi-GPU build only adds the option `qpop-diffraction --mgpu` and the program that carries it out.

## Standard build

Python 3.9 or newer. From the repository root:

```bash
pip install -r requirements.txt
pip install -e .
```

- The first command installs the Python packages: numpy and h5py (required), matplotlib (for reciprocal-space maps) and pytest (for the tests).
- The second installs the `qpop_diffraction` package from this directory, with its commands `qpop-diffraction` and `qpop-diffraction-rsm`.

For GPU support, also install CuPy, following its installation guide (https://docs.cupy.dev/en/stable/install.html).

When CuPy imports, the whole calculation runs on the GPU and the code prints `CuPy detected - computing on the GPU`. Without CuPy it runs with NumPy on the CPU. Set `QPOP_DIFFRACTION_DEVICE=cpu` to force the CPU even when CuPy is installed. `pytest tests/test_gpu.py` checks that CuPy works.

## Multi-GPU build

Do the standard build first, then install the requirements below and compile the program.

### Requirements

- NVIDIA GPUs and the CUDA toolkit (nvcc and cuFFT).
- An MPI library with CUDA-aware (GPU-aware) support.
- HeFFTe 2.4.1, built with the cuFFT backend and GPU-aware MPI (see below).
- HDF5 built with parallel (MPI-IO) support.
- A C++17 compiler.

The multi-GPU build is tested with:

- CUDA 12.9 on NVIDIA A100 GPUs
- GCC 14
- Cray MPICH 9.0.1
- HeFFTe 2.4.1
- HDF5 1.14.3.9 (parallel)

Other versions may work but have not been tried.

### Installing HeFFTe

HeFFTe (Highly Efficient FFT for Exascale, https://github.com/icl-utk-edu/heffte) provides the distributed 3-D FFTs. It is not part of this repository and must be installed first. Build version 2.4.1 with CMake, using the same compiler, CUDA and MPI you will use for the multi-GPU build:

```bash
git clone --branch v2.4.1 https://github.com/icl-utk-edu/heffte.git
cmake -S heffte -B heffte/build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_COMPILER=mpicxx \
    -DCMAKE_CUDA_ARCHITECTURES=80 \
    -DHeffte_ENABLE_CUDA=ON \
    -DHeffte_ENABLE_GPU_AWARE_MPI=ON \
    -DBUILD_SHARED_LIBS=ON \
    -DCMAKE_INSTALL_PREFIX=$HOME/heffte \
    -DCMAKE_INSTALL_RPATH=$HOME/heffte/lib64 \
    -DCMAKE_BUILD_WITH_INSTALL_RPATH=ON
cmake --build heffte/build -j 16
cmake --install heffte/build
```

What the options do:

- `CMAKE_CXX_COMPILER=mpicxx` compiles HeFFTe with MPI. Use your system's MPI C++ wrapper, e.g. `CC` on HPE Cray systems.
- `CMAKE_CUDA_ARCHITECTURES` is the GPU architecture: 80 for A100, 90 for H100.
- `Heffte_ENABLE_CUDA=ON` selects the cuFFT backend, the one `qpop_diffraction_mgpu` uses.
- `Heffte_ENABLE_GPU_AWARE_MPI=ON` lets MPI move GPU buffers directly.
- The two `RPATH` options make the installed library find its dependencies at run time.

If nvcc cannot find `mpi.h` (common when the MPI wrapper is not `mpicxx`), add `-DCMAKE_CUDA_FLAGS="-I/path/to/mpi/include"`.

The result is installed under `$HOME/heffte` (`include/` and `lib64/` or `lib/`).

### Compiling qpop_diffraction_mgpu

With the same compiler, CUDA and MPI in your environment:

From the repository root:

```bash
make HEFFTE_DIR=$HOME/heffte HDF5_DIR=/path/to/parallel/hdf5 MPI_DIR=/path/to/mpi
```

`make` runs the `Makefile` at the repository root. It compiles `src/heffte/qpop_diffraction_mgpu.cu`, together with the headers next to it, using nvcc. It links the result against HeFFTe, cuFFT, parallel HDF5 and MPI, and writes the executable `bin/qpop_diffraction_mgpu`. `qpop-diffraction --mgpu` and the tests look for the program on `PATH` first, then in `bin/`, so you can also copy it to a directory on your `PATH`.

The variables:

- `HEFFTE_DIR`: the HeFFTe install prefix. Defaults to `$Heffte_ROOT` if set.
- `HDF5_DIR`: the parallel HDF5 prefix, with `include/` and `lib/`. Defaults to `$HDF5_DIR` if set, as many HDF5 modules do.
- `MPI_DIR`: the MPI prefix, with `include/mpi.h` and `lib/libmpi`. Defaults to `$MPICH_DIR` if set. nvcc does not compile through the MPI compiler wrapper, so MPI is passed explicitly.
- `CUDA_ARCH`: the GPU architecture, default `sm_80` (A100); e.g. `sm_90` for H100.
- `CUDA_HOME`: the CUDA toolkit, by default found from the location of `nvcc`.
- `CUFFT_DIR`: the directory holding `libcufft.so`. It is found in `$CUDA_HOME/lib64`, or in `math_libs/<version>/lib64` of the NVIDIA HPC SDK, where cuFFT is installed apart from the compiler.

On HPE Cray systems, the Makefile links Cray MPICH's GPU transport layer (GTL) automatically when `CRAY_MPICH_ROOTDIR` is set. Without it, Cray MPICH aborts CUDA programs with "GPU_SUPPORT_ENABLED is requested, but GTL library is not linked". There, loading the programming environment, CUDA, MPICH and parallel HDF5 modules sets `MPICH_DIR` and `HDF5_DIR`, so only `HEFFTE_DIR` needs to be passed.

At run time, enable GPU-aware MPI. For Cray MPICH, that is `export MPICH_GPU_SUPPORT_ENABLED=1`. To check the build, run `tests/test_mgpu.py` ([Testing](testing.md)).

Do not recompile while jobs that use the program are queued or running. Overwriting the executable makes those jobs fail at launch with "Permission denied".
