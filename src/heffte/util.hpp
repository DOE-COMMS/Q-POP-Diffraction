#ifndef UTIL_HPP_
#define UTIL_HPP_

#include <cstdlib>
#include <iostream>
#include <string>

#include <cuda_runtime.h>
#include <mpi.h>

// Any rank can fail; the others may be waiting in a collective, so abort them all.
[[noreturn]] inline void fail(const std::string& message) {
    std::cerr << "qpop_diffraction_mgpu: " << message << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
    std::exit(1);
}

#define CUDA_CHECK(call)                                                                  \
    do {                                                                                  \
        cudaError_t err = (call);                                                         \
        if (err != cudaSuccess) fail(std::string(#call) + ": " + cudaGetErrorString(err)); \
    } while (0)

inline double usedDeviceMemoryGB() {
    size_t free, total;
    CUDA_CHECK(cudaMemGetInfo(&free, &total));
    return (total - free) / 1e9;
}

#endif
