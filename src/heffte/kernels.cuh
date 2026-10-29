#ifndef KERNELS_CUH_
#define KERNELS_CUH_

#include <cufft.h>

// The rank's block of the cropped grid: all x, rows yLow .. yLow + nyl - 1, all z.
// Cell i = x + nx * (j + nyl * z), the same layout in real and Fourier space.
struct slab {
    long nx, ny, nz, nyl, yLow;
    double dx, dy, dz;
    long n() const { return nx * nyl * nz; }
};

// 2 pi * numpy.fft.fftfreq(n, h)[i]
__device__ inline double waveNumber(long i, long n, double h) {
    long m = i <= (n - 1) / 2 ? i : i - n;
    return 2.0 * M_PI * double(m) / (double(n) * h);
}

__device__ inline void cellIndex(long i, const slab& s, long& x, long& j, long& z) {
    x = i % s.nx;
    j = (i / s.nx) % s.nyl;
    z = i / (s.nx * s.nyl);
}

__global__ void scale(double* a, long n, double factor) {
    long i = blockIdx.x * long(blockDim.x) + threadIdx.x;
    if (i < n) a[i] *= factor;
}

// One thread per (x, y) column of the full-z material block (nzFull, nyl, nx):
// u_z += (cumsum_z dc[mat] - dc[mat] / 2) * dz over the full stack, then cropped to z0 .. z0 + nz - 1;
// phase = phaseOf[mat] on the crop.
__global__ void columns(const int* mat, long nzFull, long z0, const double* dc, const int* phaseOf,
                        double* uz, signed char* phase, slab s) {
    long c = blockIdx.x * long(blockDim.x) + threadIdx.x;
    long nc = s.nx * s.nyl;
    if (c >= nc) return;

    double sum = 0.0;
    for (long z = 0; z < nzFull; z++) {
        int m = mat[c + nc * z];
        sum += dc[m];
        if (z >= z0 && z < z0 + s.nz) {
            long i = c + nc * (z - z0);
            uz[i] += (sum - 0.5 * dc[m]) * s.dz;
            phase[i] = phaseOf[m];
        }
    }
}

// buf = sum_p w_p f[p, a] exp(-i qC . u_ap) * (1 for pass 0, u_ap[pass - 1] otherwise), u_ap = u + su[:, :, p, a] . P.
// w_p is 1 on the cell's phase, or the phase fraction when given, times the region weight when given;
// cells that do not scatter (phase -1) give 0.
__global__ void source(cufftDoubleComplex* buf, const double* P, const double* U, const signed char* phase,
                       const double* fraction, const double* region, const double* su, const double* fRe,
                       const double* fIm, double3 qC, int nStruc, int nPhase, int nAtom, int a, int pass, long n) {
    long i = blockIdx.x * long(blockDim.x) + threadIdx.x;
    if (i >= n) return;

    double re = 0.0, im = 0.0;
    int cellPhase = phase[i];
    if (cellPhase >= 0) {
        double r = region ? region[i] : 1.0;
        int pLow = fraction ? 0 : cellPhase;
        int pHigh = fraction ? nPhase : cellPhase + 1;
        for (int p = pLow; p < pHigh; p++) {
            double w = r * (fraction ? fraction[p * n + i] : 1.0);
            if (w == 0.0) continue;

            double ua[3];
            for (int d = 0; d < 3; d++) {
                ua[d] = U[d * n + i];
                for (int k = 0; k < nStruc; k++) ua[d] += su[((d * nStruc + k) * nPhase + p) * nAtom + a] * P[k * n + i];
            }

            double sn, cs;
            sincos(-(qC.x * ua[0] + qC.y * ua[1] + qC.z * ua[2]), &sn, &cs);
            double fr = fRe[p * nAtom + a], fi = fIm[p * nAtom + a];
            if (pass > 0) w *= ua[pass - 1];
            re += w * (fr * cs - fi * sn);
            im += w * (fr * sn + fi * cs);
        }
    }
    buf[i] = make_cuDoubleComplex(re, im);
}

// amp += exp(-i (k + G) . R_a) / V * buf for pass 0, and the same times -i k[pass - 1] otherwise.
__global__ void accumulate(cufftDoubleComplex* amp, const cufftDoubleComplex* buf, double3 R, double3 G,
                           double V, int pass, slab s) {
    long i = blockIdx.x * long(blockDim.x) + threadIdx.x;
    if (i >= s.n()) return;

    long x, j, z;
    cellIndex(i, s, x, j, z);
    double k[3] = {waveNumber(x, s.nx, s.dx), waveNumber(s.yLow + j, s.ny, s.dy), waveNumber(z, s.nz, s.dz)};

    double sn, cs;
    sincos(-((k[0] + G.x) * R.x + (k[1] + G.y) * R.y + (k[2] + G.z) * R.z), &sn, &cs);
    double re = cs / V, im = sn / V;
    if (pass > 0) {
        double kc = k[pass - 1];
        double t = re;
        re = kc * im;
        im = -kc * t;
    }

    cufftDoubleComplex b = buf[i];
    amp[i].x += re * b.x - im * b.y;
    amp[i].y += re * b.y + im * b.x;
}

// I = max(|S amp|^2, 1e-50), S = prod_d sinc(k_d h_d / 2), fftshifted along all three axes.
// Rows are stored sorted by the shifted y' (see writeIntensity): y' < y0 first, then y' >= y0.
__global__ void intensity(double* out, const cufftDoubleComplex* amp, slab s) {
    long i = blockIdx.x * long(blockDim.x) + threadIdx.x;
    if (i >= s.n()) return;

    long x, j, z;
    cellIndex(i, s, x, j, z);
    double k[3] = {waveNumber(x, s.nx, s.dx), waveNumber(s.yLow + j, s.ny, s.dy), waveNumber(z, s.nz, s.dz)};
    double h[3] = {s.dx, s.dy, s.dz};

    double S = 1.0;
    for (int d = 0; d < 3; d++) {
        double t = 0.5 * k[d] * h[d];
        if (t != 0.0) S *= sin(t) / t;
    }
    double I = fmax(S * S * (amp[i].x * amp[i].x + amp[i].y * amp[i].y), 1e-50);

    long y0 = (s.yLow + s.ny / 2) % s.ny;
    long wrapped = s.nyl - min(s.nyl, s.ny - y0);
    long ys = (s.yLow + j + s.ny / 2) % s.ny;
    long row = ys < y0 ? ys : ys - y0 + wrapped;
    out[(x + s.nx / 2) % s.nx + s.nx * (row + s.nyl * ((z + s.nz / 2) % s.nz))] = I;
}

#endif
