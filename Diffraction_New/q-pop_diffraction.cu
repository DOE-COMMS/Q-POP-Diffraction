#include <array>
#include <complex>
#include <iostream>
#include <random>
#include <vector>
#include <cufftXt.h>

#include "src/utils.hpp"
#include "src/sizeContext.hpp"
#include "src/readArrays.hpp"
#include "src/atomList.hpp"
#include "src/constants.hpp"

__global__
void scaling_kernel(cufftComplex* data, int element_count, float scale) {
    const int tid = threadIdx.x + blockIdx.x * blockDim.x;
    const int stride = blockDim.x * gridDim.x;
    for (auto i = tid; i < element_count; i += stride) {
        data[i].x *= scale;
        data[i].y *= scale;
    }
}

using cpudata_t = std::vector<std::complex<float>>;
using gpus_t = std::vector<int>;
using dim_t = std::array<size_t, 3>;

int main(int argc, char* argv[]) {
    int readErrors = 0;

    sizeContext params;
    atomList atoms;

    params.read("parameter.system.in");
    params.setup();

    // oPhase from phaseFra.in
    size_t oPhase_size = size_t(params.nPhase * params.nx * params.ny * params.nz);
    std::vector<double> oPhase(oPhase_size, 0.0);

    readErrors = read4D("phaseFra.in", oPhase, params.nPhase, params.nx, params.ny, params.nz);
    if (-1 == readErrors) {
        std::cout << "File phaseFra.in not provided. Using a pure phase 1." << std::endl;
        for (int i = 0; i < params.nx; i++) {
            for (int j = 0; j < params.ny; j++) {
                for (int k = 0; k < params.nz; k++) {
                    oPhase[i * params.ny * params.nz + j * params.nz + k] = 1.0;
                }
            }
        }
    }
    else if (0 == readErrors) {
        std::cerr << "Error reading phaseFra.in. Exiting." << std::endl;
        return EXIT_FAILURE;
    }
    else if (1 == readErrors) {
        std::cout << "Successfully read phaseFra.in." << std::endl;

        for (int p = 0; p < params.nPhase; p++) {
            for (int i = 0; i < params.nx; i++) {
                for (int j = 0; j < params.ny; j++) {
                    for (int k = 0; k < params.nz; k++) {
                        if (k < params.k1 - 2 || k > params.k2)
                            oPhase[p * params.n + i * params.ny * params.nz + j * params.nz + k] = 1.0;
                    }
                }
            }
        }
    }

    // oStruc from strucOrd.in
    size_t oStruc_size = size_t(params.nStruc * params.nx * params.ny * params.nz);
    std::vector<double> oStruc(oStruc_size, 0.0);

    readErrors = read4D("strucOrd.in", oStruc, params.nStruc, params.nx, params.ny, params.nz);
    if (-1 == readErrors) {
        if (params.nStruc > 0) {
            std::cout << "File strucOrd.in not provided. Skipping." << std::endl;
            return EXIT_FAILURE;
        }
    }
    else if (0 == readErrors) {
        std::cerr << "Error reading strucOrd.in. Exiting." << std::endl;
        return EXIT_FAILURE;
    }
    else if (1 == readErrors && params.nStruc > 0) {
        std::cout << "Successfully read strucOrd.in." << std::endl;

        for (auto& i : oStruc) {
            i /= os0;
        }
    }


    // u (displacement field) from displace.in
    size_t u_size = size_t(3 * params.nx * params.ny * params.nz);
    std::vector<double> u(u_size, 0.0);

    readErrors = read4D("displace.in", u, 3, params.nx, params.ny, params.nz);
    if (-1 == readErrors) {
        std::cout << "File displace.in not provided. Using a zero displacement field." << std::endl;
        for (auto &i : u) {
            i = 0.0;
        }
    }
    else if (0 == readErrors) {
        std::cerr << "Error reading displace.in. Exiting." << std::endl;
        return EXIT_FAILURE;
    }
    else if (1 == readErrors) {
        std::cout << "Successfully read displace.in." << std::endl;
    }

    atoms.read("parameter.atom.in", params.nPhase, params.nStruc);
    atoms.setup();

    // Region
    size_t region_size = size_t(params.nx * params.ny * params.nz);
    std::vector<double> region(region_size, 0.0);
    readErrors = read4D("region.in", region, 1, params.nx, params.ny, params.nz);
    if (-1 == readErrors) {
        std::cout << "File region.in not provided. Using a default region." << std::endl;

        double ix1 = (params.nx + 1) / 2.0 - params.nx / M_PI;
        double iy1 = (params.ny + 1) / 2.0 - params.ny / M_PI;
        double iz1 = (params.nz + 1) / 2.0 - params.nz / M_PI;

        double ix2 = (params.nx + 1) / 2.0 + params.nx / M_PI;
        double iy2 = (params.ny + 1) / 2.0 + params.ny / M_PI;
        double iz2 = (params.nz + 1) / 2.0 + params.nz / M_PI;

        double tx = params.nx / 8.0;
        double ty = params.ny / 8.0;
        double tz = params.nz / 8.0;

        std::fill(region.begin(), region.end(), 1.0);

        if (params.nf == 0) {
            for (int i = 0; i < params.nx; i++) {
                for (int j = 0; j < params.ny; j++) {
                    for (int k = 0; k < params.nz; k++) {
                        int idx = i * params.ny * params.nz + j * params.nz + k;
                        region[idx] = tanh((k - 1 - iz1) / tz) - tanh((k - 1 - iz2) / tz);
                        region[idx] *= tanh((j - 1 - iy1) / ty) - tanh((j - 1 - iy2) / ty);
                        region[idx] *= tanh((i - 1 - ix1) / tx) - tanh((i - 1 - ix2) / tx);
                    }
                }
            }
        }

        double region_max = *std::max_element(region.begin(), region.end());
        for (auto &i : region) {
            i /= region_max;
        }
    }
    else if (0 == readErrors) {
        std::cerr << "Error reading region.in. Exiting." << std::endl;
        return EXIT_FAILURE;
    }
    else if (1 == readErrors) {
        std::cout << "Successfully read region.in." << std::endl;
    }

    return EXIT_SUCCESS;
};