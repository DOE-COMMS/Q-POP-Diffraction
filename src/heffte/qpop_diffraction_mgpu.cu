/*
 * Multi-GPU kinematic diffraction for grids larger than one GPU.
 *
 *     qpop-diffraction --mgpu ...                      (writes <output>_params.h5)
 *     srun ... --gpu-bind=none qpop_diffraction_mgpu <output>_params.h5
 *
 * Same calculation as qpop_diffraction/diffraction.py, with the grid split into
 * y-slabs, one per rank and GPU, and the FFTs done by HeFFTe:
 *
 *     Amp(k) = S(k) sum_a exp(-i (k + G) . R_a) / V [FFT(fExpN_a) - i k . FFT(fExpN_a u_a)]
 *
 * Each rank reads its slab of P, u, the material field and, when given, the phase
 * fractions and scattering weight with parallel HDF5, and
 * writes its slab of /intensity into the output file prepared by the Python side.
 */

#include <algorithm>
#include <cstdio>
#include <string>
#include <vector>

#include <cuda_runtime.h>
#include <heffte.h>
#include <mpi.h>

#include "util.hpp"
#include "params.hpp"
#include "hdf5_io.hpp"
#include "kernels.cuh"

constexpr int THREADS = 256;

inline unsigned blocks(long n) { return unsigned((n + THREADS - 1) / THREADS); }

// Everything that holds MPI or GPU resources lives here, so it is released before MPI_Finalize.
static void run(const char* paramsFile) {
    int rank, nRanks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nRanks);

    // One GPU per rank on each node
    MPI_Comm node;
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL, &node);
    int localRank, nDevices;
    MPI_Comm_rank(node, &localRank);
    CUDA_CHECK(cudaGetDeviceCount(&nDevices));
    CUDA_CHECK(cudaSetDevice(localRank % nDevices));

    mgpuParams prm;
    prm.read(paramsFile);

    slab s;
    s.nx = prm.nx();
    s.ny = prm.ny();
    s.nz = prm.nz();
    s.dx = prm.spacing[0];
    s.dy = prm.spacing[1];
    s.dz = prm.spacing[2];
    if (nRanks > s.ny) fail("More ranks (" + std::to_string(nRanks) + ") than rows in y (" + std::to_string(s.ny) + ").");
    s.yLow = rank * (s.ny / nRanks) + std::min<long>(rank, s.ny % nRanks);
    s.nyl = s.ny / nRanks + (rank < s.ny % nRanks ? 1 : 0);
    long n = s.n();

    if (rank == 0) {
        std::printf("Grid %ld x %ld x %ld (crop of %ld x %ld x %ld), %d ranks, %d GPU(s) per node\n",
                    s.nx, s.ny, s.nz, prm.full[0], prm.full[1], prm.full[2], nRanks, nDevices);
        std::printf("%d atom(s), %d phase(s)\n", prm.nAtom, prm.nPhase);
    }
    double t0 = MPI_Wtime();

    // ── Read P, u and the material field of this slab ──────────────────────────
    std::array<hsize_t, 3> start = {hsize_t(prm.crop[4]), hsize_t(prm.crop[2] + s.yLow), hsize_t(prm.crop[0])};
    std::array<hsize_t, 3> count = {hsize_t(s.nz), hsize_t(s.nyl), hsize_t(s.nx)};
    std::vector<double> host;

    double *P, *U;
    CUDA_CHECK(cudaMalloc(&P, 3 * n * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&U, 3 * n * sizeof(double)));

    hid_t file = openParallel(prm.plattFile, H5F_ACC_RDONLY, MPI_COMM_WORLD);
    for (int d = 0; d < 3; d++) {
        readBlock(file, prm.plattFile, prm.plattPrefix + "_" + "xyz"[d], H5T_NATIVE_DOUBLE, prm.full, start, count, host);
        CUDA_CHECK(cudaMemcpy(P + d * n, host.data(), n * sizeof(double), cudaMemcpyHostToDevice));
    }
    H5Fclose(file);

    if (prm.displacementFile.empty()) {
        CUDA_CHECK(cudaMemset(U, 0, 3 * n * sizeof(double)));
    } else {
        file = openParallel(prm.displacementFile, H5F_ACC_RDONLY, MPI_COMM_WORLD);
        for (int d = 0; d < 3; d++) {
            readBlock(file, prm.displacementFile, std::string("Displacement_") + "xyz"[d], H5T_NATIVE_DOUBLE, prm.full, start, count, host);
            CUDA_CHECK(cudaMemcpy(U + d * n, host.data(), n * sizeof(double), cudaMemcpyHostToDevice));
        }
        H5Fclose(file);
        scale<<<blocks(3 * n), THREADS>>>(U, 3 * n, prm.uScale);
    }

    double *fraction = nullptr, *region = nullptr;
    if (!prm.phaseFile.empty()) {
        CUDA_CHECK(cudaMalloc(&fraction, prm.nPhase * n * sizeof(double)));
        file = openParallel(prm.phaseFile, H5F_ACC_RDONLY, MPI_COMM_WORLD);
        for (int p = 0; p < prm.nPhase; p++) {
            readBlock(file, prm.phaseFile, "phase_" + std::to_string(p + 1), H5T_NATIVE_DOUBLE, prm.full, start, count, host);
            CUDA_CHECK(cudaMemcpy(fraction + p * n, host.data(), n * sizeof(double), cudaMemcpyHostToDevice));
        }
        H5Fclose(file);
    }
    if (!prm.regionFile.empty()) {
        CUDA_CHECK(cudaMalloc(&region, n * sizeof(double)));
        file = openParallel(prm.regionFile, H5F_ACC_RDONLY, MPI_COMM_WORLD);
        readBlock(file, prm.regionFile, "region", H5T_NATIVE_DOUBLE, prm.full, start, count, host);
        CUDA_CHECK(cudaMemcpy(region, host.data(), n * sizeof(double), cudaMemcpyHostToDevice));
        H5Fclose(file);
    }
    std::vector<double>().swap(host);

    std::vector<int> mat;
    file = openParallel(prm.materialFile, H5F_ACC_RDONLY, MPI_COMM_WORLD);
    readBlock(file, prm.materialFile, "materialField", H5T_NATIVE_INT, prm.full,
              {0, start[1], start[2]}, {hsize_t(prm.full[2]), count[1], count[2]}, mat);
    H5Fclose(file);
    for (int m : mat)
        if (m < 0 || m >= prm.nMaterial) fail("materialField holds material " + std::to_string(m) + ", outside the run's MaterialNames.");

    // ── u_z reference ramp and the phase of each cell ──────────────────────────
    int *dMat, *dPhaseOf;
    double* dDc;
    signed char* phase;
    CUDA_CHECK(cudaMalloc(&dMat, mat.size() * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&dPhaseOf, prm.nMaterial * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&dDc, prm.nMaterial * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&phase, n * sizeof(signed char)));
    CUDA_CHECK(cudaMemcpy(dMat, mat.data(), mat.size() * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(dPhaseOf, prm.phaseOf.data(), prm.nMaterial * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(dDc, prm.dc.data(), prm.nMaterial * sizeof(double), cudaMemcpyHostToDevice));
    columns<<<blocks(s.nx * s.nyl), THREADS>>>(dMat, prm.full[2], prm.crop[4], dDc, dPhaseOf, U + 2 * n, phase, s);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaFree(dMat));
    std::vector<int>().swap(mat);

    double t1 = MPI_Wtime();

    // ── Amplitude: per atom, FFT of fExpN and of fExpN u_x, u_y, u_z ───────────
    double *su, *fRe, *fIm;
    CUDA_CHECK(cudaMalloc(&su, prm.su.size() * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&fRe, prm.fRe.size() * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&fIm, prm.fIm.size() * sizeof(double)));
    CUDA_CHECK(cudaMemcpy(su, prm.su.data(), prm.su.size() * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(fRe, prm.fRe.data(), prm.fRe.size() * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(fIm, prm.fIm.data(), prm.fIm.size() * sizeof(double), cudaMemcpyHostToDevice));

    heffte::box3d<> box({0, int(s.yLow), 0}, {int(s.nx - 1), int(s.yLow + s.nyl - 1), int(s.nz - 1)});
    auto options = heffte::default_options<heffte::backend::cufft>();
    options.algorithm = heffte::reshape_algorithm::p2p_plined;
    heffte::fft3d<heffte::backend::cufft> fft(box, box, MPI_COMM_WORLD, options);

    cufftDoubleComplex *amp, *buf, *work;
    CUDA_CHECK(cudaMalloc(&amp, n * sizeof(cufftDoubleComplex)));
    CUDA_CHECK(cudaMalloc(&buf, n * sizeof(cufftDoubleComplex)));
    CUDA_CHECK(cudaMalloc(&work, fft.size_workspace() * sizeof(cufftDoubleComplex)));
    CUDA_CHECK(cudaMemset(amp, 0, n * sizeof(cufftDoubleComplex)));
    double memory = usedDeviceMemoryGB();

    double3 qC = make_double3(prm.qC[0], prm.qC[1], prm.qC[2]);
    double3 G = make_double3(prm.G[0], prm.G[1], prm.G[2]);
    double tFFT = 0.0;
    for (int a = 0; a < prm.nAtom; a++) {
        double3 R = make_double3(prm.R[3 * a], prm.R[3 * a + 1], prm.R[3 * a + 2]);
        for (int pass = 0; pass < 4; pass++) {
            source<<<blocks(n), THREADS>>>(buf, P, U, phase, fraction, region, su, fRe, fIm, qC, prm.nStruc, prm.nPhase, prm.nAtom, a, pass, n);
            CUDA_CHECK(cudaGetLastError());
            CUDA_CHECK(cudaDeviceSynchronize());
            double tf = MPI_Wtime();
            fft.forward(buf, buf, work, heffte::scale::full);
            CUDA_CHECK(cudaDeviceSynchronize());
            tFFT += MPI_Wtime() - tf;
            accumulate<<<blocks(n), THREADS>>>(amp, buf, R, G, prm.V, pass, s);
            CUDA_CHECK(cudaGetLastError());
        }
    }

    // ── Intensity, fftshifted, written in place of the FFT buffer ──────────────
    double* I = reinterpret_cast<double*>(buf);
    intensity<<<blocks(n), THREADS>>>(I, amp, s);
    CUDA_CHECK(cudaGetLastError());
    std::vector<double> out(n);
    CUDA_CHECK(cudaMemcpy(out.data(), I, n * sizeof(double), cudaMemcpyDeviceToHost));
    double t2 = MPI_Wtime();

    writeIntensity(prm.outputFile, MPI_COMM_WORLD, s.nx, s.ny, s.nz, s.yLow, s.nyl, out);
    double t3 = MPI_Wtime();

    double maxMemory;
    MPI_Reduce(&memory, &maxMemory, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        std::printf("Read %.1f s, amplitude %.1f s (FFT %.1f s, %d transforms), write %.1f s, total %.1f s\n",
                    t1 - t0, t2 - t1, tFFT, 4 * prm.nAtom, t3 - t2, t3 - t0);
        std::printf("Device memory used, max over ranks: %.1f GB\n", maxMemory);
        std::printf("Wrote /intensity to %s\n", prm.outputFile.c_str());
    }

    for (void* ptr : {(void*)P, (void*)U, (void*)phase, (void*)dPhaseOf, (void*)dDc, (void*)su, (void*)fRe, (void*)fIm,
                      (void*)amp, (void*)buf, (void*)work, (void*)fraction, (void*)region})
        CUDA_CHECK(cudaFree(ptr));
    MPI_Comm_free(&node);
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    bool help = argc == 2 && (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help");
    if (argc != 2 || help) {
        if (rank == 0)
            std::printf("Usage: srun ... --gpu-bind=none qpop_diffraction_mgpu <output>_params.h5\n\n"
                        "Multi-GPU kinematic diffraction. Prepare <output>_params.h5 and <output>.h5 with\n"
                        "qpop-diffraction --mgpu; this fills in <output>.h5:/intensity. One rank per GPU.\n");
        MPI_Finalize();
        return help ? 0 : 1;
    }

    run(argv[1]);
    MPI_Finalize();
    return 0;
}
