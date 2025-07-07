#ifndef DIFFRACTION_SETUP_HPP_
#define DIFFRACTION_SETUP_HPP_

#include <cufft.h>
#include <complex>
#include <array>
#include <vector>

#include "writeArrays.hpp"

struct diffraction {
    double VCell;
    double bC[3][3];
    int hM1, hM2, hM3;
    double G10, G20, G30;
    double qC10, qC20, qC30;
    double qC1, qC2, qC3;

    std::vector<double> mqk1, mqk2, mqk3;
    std::vector<double> mqk1_out, mqk2_out, mqk3_out;
    std::vector<double> mqkA1, mqkA2, mqkA3;
    std::vector<double> mqkB1, mqkB2, mqkB3;

    std::vector<std::complex<double>> fA, fB;
    std::vector<std::complex<double>> kfExpA, kfExpB;

    std::vector<double> dmqk1, dmqk2, dmqk3;

    std::vector<double> DRAtom;
    int nAtom;

    void diffraction_setup(sizeContext* params, atomList* atoms);

    void diffraction_calc(sizeContext* params, atomList* atoms,
        std::vector<double>& IDiffr, std::vector<double>& region,
        std::vector<double>& QCenter, std::vector<double>& oPhase,
        std::vector<double>& oStruc, std::vector<double>& u);

    void ArrayFourierToRegular(const std::vector<double> A, const std::vector<double> B, std::vector<double>& out, sizeContext* params);
};

void diffraction::diffraction_setup(sizeContext* params, atomList* atoms) {
    VCell = static_cast<double>(atoms->aC[0][0] * atoms->aC[1][1] * atoms->aC[2][2]) + \
        static_cast<double>(atoms->aC[0][1] * atoms->aC[1][2] * atoms->aC[2][0]) + \
        static_cast<double>(atoms->aC[0][2] * atoms->aC[1][0] * atoms->aC[2][1]) - \
        static_cast<double>(atoms->aC[0][0] * atoms->aC[1][2] * atoms->aC[2][1]) - \
        static_cast<double>(atoms->aC[0][1] * atoms->aC[1][0] * atoms->aC[2][2]) - \
        static_cast<double>(atoms->aC[0][2] * atoms->aC[1][1] * atoms->aC[2][0]);

    bC[0][0] = 2.0 * M_PI / VCell * static_cast<double>(atoms->aC[1][1] * atoms->aC[2][2] - atoms->aC[1][2] * atoms->aC[2][1]);
    bC[0][1] = 2.0 * M_PI / VCell * static_cast<double>(atoms->aC[1][2] * atoms->aC[2][0] - atoms->aC[1][0] * atoms->aC[2][2]);
    bC[0][2] = 2.0 * M_PI / VCell * static_cast<double>(atoms->aC[1][0] * atoms->aC[2][1] - atoms->aC[1][1] * atoms->aC[2][0]);

    bC[1][0] = 2.0 * M_PI / VCell * static_cast<double>(atoms->aC[0][2] * atoms->aC[2][1] - atoms->aC[0][1] * atoms->aC[2][2]);
    bC[1][1] = 2.0 * M_PI / VCell * static_cast<double>(atoms->aC[0][0] * atoms->aC[2][2] - atoms->aC[0][2] * atoms->aC[2][0]);
    bC[1][2] = 2.0 * M_PI / VCell * static_cast<double>(atoms->aC[0][1] * atoms->aC[2][0] - atoms->aC[0][0] * atoms->aC[2][1]);

    bC[2][0] = 2.0 * M_PI / VCell * static_cast<double>(atoms->aC[0][1] * atoms->aC[1][2] - atoms->aC[0][2] * atoms->aC[1][1]);
    bC[2][1] = 2.0 * M_PI / VCell * static_cast<double>(atoms->aC[0][2] * atoms->aC[1][0] - atoms->aC[0][0] * atoms->aC[1][2]);
    bC[2][2] = 2.0 * M_PI / VCell * static_cast<double>(atoms->aC[0][0] * atoms->aC[1][1] - atoms->aC[0][1] * atoms->aC[1][0]);

    hM1 = std::lround((atoms->q00[0] * atoms->aC[0][0] + atoms->q00[1] * atoms->aC[0][1] + atoms->q00[2] * atoms->aC[0][2]) / (2 * M_PI));
    hM2 = std::lround((atoms->q00[0] * atoms->aC[1][0] + atoms->q00[1] * atoms->aC[1][1] + atoms->q00[2] * atoms->aC[1][2]) / (2 * M_PI));
    hM3 = std::lround((atoms->q00[0] * atoms->aC[2][0] + atoms->q00[1] * atoms->aC[2][1] + atoms->q00[2] * atoms->aC[2][2]) / (2 * M_PI));

    G10 = hM1 * bC[0][0] + hM2 * bC[1][0] + hM3 * bC[2][0];
    G20 = hM1 * bC[0][1] + hM2 * bC[1][1] + hM3 * bC[2][1];
    G30 = hM1 * bC[0][2] + hM2 * bC[1][2] + hM3 * bC[2][2];

    qC10 = G10;
    qC20 = G20;
    qC30 = G30;

    mqk1.resize(params->Cn1 * params->Cn2 * params->Cn3);
    mqk2.resize(params->Cn1 * params->Cn2 * params->Cn3);
    mqk3.resize(params->Cn1 * params->Cn2 * params->Cn3);

    mqkA1.resize(params->Cn1 * params->Cn2 * params->Cn3);
    mqkA2.resize(params->Cn1 * params->Cn2 * params->Cn3);
    mqkA3.resize(params->Cn1 * params->Cn2 * params->Cn3);

    mqkB1.resize(params->Cn1 * params->Cn2 * params->Cn3);
    mqkB2.resize(params->Cn1 * params->Cn2 * params->Cn3);
    mqkB3.resize(params->Cn1 * params->Cn2 * params->Cn3);

    for (int i = 0; i < params->Cn1; i++) {
        for (int j = 0; j < params->Cn2; j++) {
            for (int k = 0; k < params->Cn3; k++) {
                int idx = i * params->Cn2 * params->Cn3 + j * params->Cn3 + k;

                mqk1[idx] = static_cast<double>(i);
                if (i >= params->nx / 2) mqk1[idx] = static_cast<double>(i - params->nx);

                mqk1[idx] *= 2.0 * M_PI / static_cast<double>(params->lx);

                mqkA1[idx] = mqk1[idx] + qC10;
                mqkB1[idx] = -mqk1[idx] + qC10;

                mqk2[idx] = static_cast<double>(j);
                if (j >= params->ny / 2) mqk2[idx] = static_cast<double>(j - params->ny);

                mqk2[idx] *= 2.0 * M_PI / static_cast<double>(params->ly);

                mqkA2[idx] = mqk2[idx] + qC20;
                mqkB2[idx] = -mqk2[idx] + qC20;

                mqk3[idx] = static_cast<double>(k);
                if (k >= params->nz / 2) mqk3[idx] = static_cast<double>(k - params->nz);

                mqk3[idx] *= 2.0 * M_PI / static_cast<double>(params->lz);

                mqkA3[idx] = mqk3[idx] + qC30;
                mqkB3[idx] = -mqk3[idx] + qC30;

                if (params->nx % 2 == 0 && i == params->nx / 2) {
                    mqkB1[idx] = mqk1[idx] + qC10;
                }

                if (params->ny % 2 == 0 && j == params->ny / 2) {
                    mqkB2[idx] = mqk2[idx] + qC20;
                }

                if (params->nz % 2 == 0 && k == params->nz / 2) {
                    mqkB3[idx] = mqk3[idx] + qC30;
                }
            }
        }

    }

    kfExpA.resize(params->Cn1 * params->Cn2 * params->Cn3);
    kfExpB.resize(params->Cn1 * params->Cn2 * params->Cn3);

    for (int i = 0; i < params->Cn1; i++) {
        for (int j = 0; j < params->Cn2; j++) {
            for (int k = 0; k < params->Cn3; k++) {
                int idx = i * params->Cn2 * params->Cn3 + j * params->Cn3 + k;

                std::complex kTemp = { 0.0, (mqkA1[idx] - G10) * params->dx };
                std::complex temp1 = { kTemp.real() / 2.0, kTemp.imag() / 2.0 };
                std::complex temp2 = { -kTemp.real() / 2.0, -kTemp.imag() / 2.0 };
                std::complex<double> kfTemp1 = 1.0 / kTemp * (exp(temp1) - exp(temp2));
                if (std::abs(kTemp) <= 1e-8) kfTemp1 = { 1.0, 0.0 };

                kTemp = { 0.0, (mqkA2[idx] - G20) * params->dy };
                temp1 = { kTemp.real() / 2.0, kTemp.imag() / 2.0 };
                temp2 = { -kTemp.real() / 2.0, -kTemp.imag() / 2.0 };
                std::complex<double> kfTemp2 = 1.0 / kTemp * (exp(temp1) - exp(temp2));
                if (std::abs(kTemp) <= 1e-8) kfTemp2 = { 1.0, 0.0 };

                kTemp = { 0.0, (mqkA3[idx] - G30) * params->dz };
                temp1 = { kTemp.real() / 2.0, kTemp.imag() / 2.0 };
                temp2 = { -kTemp.real() / 2.0, -kTemp.imag() / 2.0 };
                std::complex<double> kfTemp3 = 1.0 / kTemp * (exp(temp1) - exp(temp2));
                if (std::abs(kTemp) <= 1e-8) kfTemp3 = { 1.0, 0.0 };

                kfExpA[idx] = kfTemp1 * kfTemp2 * kfTemp3;

                kTemp = { 0.0, (mqkB1[idx] - G10) * params->dx };
                temp1 = { kTemp.real() / 2.0, kTemp.imag() / 2.0 };
                temp2 = { -kTemp.real() / 2.0, -kTemp.imag() / 2.0 };
                kfTemp1 = 1.0 / kTemp * (exp(temp1) - exp(temp2));
                if (std::abs(kTemp) <= 1e-8) kfTemp1 = { 1.0, 0.0 };

                kTemp = { 0.0, (mqkB2[idx] - G20) * params->dy };
                temp1 = { kTemp.real() / 2.0, kTemp.imag() / 2.0 };
                temp2 = { -kTemp.real() / 2.0, -kTemp.imag() / 2.0 };
                kfTemp2 = 1.0 / kTemp * (exp(temp1) - exp(temp2));
                if (std::abs(kTemp) <= 1e-8) kfTemp2 = { 1.0, 0.0 };

                kTemp = { 0.0, (mqkB3[idx] - G30) * params->dz };
                temp1 = { kTemp.real() / 2.0, kTemp.imag() / 2.0 };
                temp2 = { -kTemp.real() / 2.0, -kTemp.imag() / 2.0 };
                kfTemp3 = 1.0 / kTemp * (exp(temp1) - exp(temp2));
                if (std::abs(kTemp) <= 1e-8) kfTemp3 = { 1.0, 0.0 };

                kfExpB[idx] = kfTemp1 * kfTemp2 * kfTemp3;
            }
        }
    }

    this->nAtom = atoms->nAtom;

    fA.resize(params->Cn1 * params->Cn2 * params->Cn3 * this->nAtom);
    fB.resize(params->Cn1 * params->Cn2 * params->Cn3 * this->nAtom);

    DRAtom.resize(3 * this->nAtom);

    for (int i = 0; i < this->nAtom; i++) {

        DRAtom[i * 3 + 0] = static_cast<double>(atoms->xAtom[i][0] * atoms->aC[0][0]) + \
            static_cast<double>(atoms->xAtom[i][1] * atoms->aC[1][0]) + \
            static_cast<double>(atoms->xAtom[i][2] * atoms->aC[2][0]);
        DRAtom[i * 3 + 1] = static_cast<double>(atoms->xAtom[i][0] * atoms->aC[0][1]) + \
            static_cast<double>(atoms->xAtom[i][1] * atoms->aC[1][1]) + \
            static_cast<double>(atoms->xAtom[i][2] * atoms->aC[2][1]);
        DRAtom[i * 3 + 2] = static_cast<double>(atoms->xAtom[i][0] * atoms->aC[0][2]) + \
            static_cast<double>(atoms->xAtom[i][1] * atoms->aC[1][2]) + \
            static_cast<double>(atoms->xAtom[i][2] * atoms->aC[2][2]);

        for (int x = 0; x < params->Cn1; x++) {
            for (int y = 0; y < params->Cn2; y++) {
                for (int z = 0; z < params->Cn3; z++) {
                    int idx = x * params->Cn2 * params->Cn3 + y * params->Cn3 + z;
                    int aidx = idx + i * (params->Cn1 * params->Cn2 * params->Cn3);

                    std::complex<double> temp = { 0.0, -mqkA1[idx] * DRAtom[i * 3] - mqkA2[idx] * DRAtom[i * 3 + 1] - mqkA3[idx] * DRAtom[i * 3 + 2] };
                    fA[aidx] = exp(temp) / VCell;

                    temp = { 0.0, -mqkB1[idx] * DRAtom[i * 3] - mqkB2[idx] * DRAtom[i * 3 + 1] - mqkB3[idx] * DRAtom[i * 3 + 2] };
                    fB[aidx] = exp(temp) / VCell;
                }
            }
        }
    }

    mqk1_out.resize(params->n);
    mqk2_out.resize(params->n);
    mqk3_out.resize(params->n);

    ArrayFourierToRegular(mqkA1, mqkB1, mqk1_out, params);
    ArrayFourierToRegular(mqkA2, mqkB2, mqk2_out, params);
    ArrayFourierToRegular(mqkA3, mqkB3, mqk3_out, params);

    for (int i = 0; i < params->Cn1 * params->Cn2 * params->Cn3; i++) {
        mqkA1[i] = mqkA1[i] - qC10;
        mqkB1[i] = mqkB1[i] - qC10;

        mqkA2[i] = mqkA2[i] - qC20;
        mqkB2[i] = mqkB2[i] - qC20;

        mqkA3[i] = mqkA3[i] - qC30;
        mqkB3[i] = mqkB3[i] - qC30;
    }
}

void diffraction::diffraction_calc(sizeContext* params, atomList* atoms,
    std::vector<double>& IDiffr, std::vector<double>& region,
    std::vector<double>& QCenter, std::vector<double>& oPhase,
    std::vector<double>& oStruc, std::vector<double>& u) {
    double G1 = G10 - qC10 * params->strainAvg[0] - qC20 * params->strainAvg[5] - qC30 * params->strainAvg[4];
    double G2 = G20 - qC10 * params->strainAvg[5] - qC20 * params->strainAvg[1] - qC30 * params->strainAvg[3];
    double G3 = G30 - qC10 * params->strainAvg[4] - qC20 * params->strainAvg[3] - qC30 * params->strainAvg[2];

    qC1 = qC10 + G1 - G10;
    qC2 = qC20 + G2 - G20;
    qC3 = qC30 + G3 - G30;

    std::vector<double> q0Gr(params->n, 0.0);
    for (int i = 0; i < params->nx; i++) {
        for (int j = 0; j < params->ny; j++) {
            for (int k = 0; k < params->nz; k++) {
                int idx = i * params->ny * params->nz + j * params->nz + k;
                q0Gr[idx] += (qC1 - G1) * static_cast<double>(i + 1) * static_cast<double>(params->dx);
                q0Gr[idx] += (qC2 - G2) * static_cast<double>(j + 1) * static_cast<double>(params->dy);
                q0Gr[idx] += (qC3 - G3) * static_cast<double>(k + 1) * static_cast<double>(params->dz);
            }
        }
    }

    std::vector<double> fExpN_r(params->n * atoms->nAtom, 0.0);
    std::vector<double> fU1ExpN_r(params->n * atoms->nAtom, 0.0);
    std::vector<double> fU2ExpN_r(params->n * atoms->nAtom, 0.0);
    std::vector<double> fU3ExpN_r(params->n * atoms->nAtom, 0.0);

    std::vector<double> fExpN_i(params->n * atoms->nAtom,  0.0);
    std::vector<double> fU1ExpN_i(params->n * atoms->nAtom, 0.0);
    std::vector<double> fU2ExpN_i(params->n * atoms->nAtom, 0.0);
    std::vector<double> fU3ExpN_i(params->n * atoms->nAtom, 0.0);

    for (int n = 0; n < atoms->nAtom; n++) {
        for (int m = 0; m < params->nPhase; m++) {
            std::vector<double> uNM1(params->n, 0.0);
            std::vector<double> uNM2(params->n, 0.0);
            std::vector<double> uNM3(params->n, 0.0);

            if (u.size() == static_cast<size_t>(3 * params->n)) {
                std::copy(u.begin(), u.begin() + static_cast<size_t>(params->n), uNM1.begin());
                std::copy(u.begin() + static_cast<size_t>(params->n), u.begin() + static_cast<size_t>(2 * params->n), uNM2.begin());
                std::copy(u.begin() + static_cast<size_t>(2 * params->n), u.begin() + static_cast<size_t>(3 * params->n), uNM3.begin());
            }

            for (int k = 0; k < params->nStruc; k++) {
                for (int i = 0; i < params->n; i++) {
                    uNM1[i] += static_cast<double>(atoms->su1[k][m][n]) * oStruc[k * params->n + i];
                    uNM2[i] += static_cast<double>(atoms->su2[k][m][n]) * oStruc[k * params->n + i];
                    uNM3[i] += static_cast<double>(atoms->su3[k][m][n]) * oStruc[k * params->n + i];
                }
            }

            for (int i = 0; i < params->n; i++) {
                std::complex<double> exponent = { 0.0, -(qC1 * uNM1[i] + qC2 * uNM2[i] + qC3 * uNM3[i] + q0Gr[i]) };
                std::complex<double> exp_term = std::exp(exponent);
                std::complex<double> fExpNM = exp_term * static_cast<std::complex<double>>(atoms->fAtom[m][n]) * oPhase[m * params->n + i] * region[i];

                int idx = n * params->n + i;
                fExpN_r[idx] += fExpNM.real();
                fU1ExpN_r[idx] += fExpNM.real() * uNM1[i];
                fU2ExpN_r[idx] += fExpNM.real() * uNM2[i];
                fU3ExpN_r[idx] += fExpNM.real() * uNM3[i];

                fExpN_i[idx] += fExpNM.imag();
                fU1ExpN_i[idx] += fExpNM.imag() * uNM1[i];
                fU2ExpN_i[idx] += fExpNM.imag() * uNM2[i];
                fU3ExpN_i[idx] += fExpNM.imag() * uNM3[i];
            }
        }
    }

    cufftDoubleReal *d_input;
    cufftDoubleComplex *d_output;

    cudaMalloc(&d_input, static_cast<size_t>(params->n * atoms->nAtom) * sizeof(cufftDoubleReal));
    cudaMalloc(&d_output, static_cast<size_t>(params->Cn1 * params->Cn2 * params->Cn3 * atoms->nAtom) * sizeof(cufftDoubleComplex));

    cufftHandle plan;
    cufftCreate(&plan);

    std::array<int, 3> fft = { params->Rn1, params->Rn2, params->Rn3 };
    int batch_size = atoms->nAtom;
    int fft_size = params->n;

    cufftPlanMany(&plan, 3, fft.data(), 
                         nullptr, 1, fft_size, // *inembed, istride, idist
                         nullptr, 1, params->Cn1 * params->Cn2 * params->Cn3, // *onembed, ostride, odist
                         CUFFT_D2Z, batch_size);

    double normalization = 1.0 / static_cast<double>(params->n);

    cudaMemcpy(d_input, fExpN_r.data(), static_cast<size_t>(params->n * atoms->nAtom) * sizeof(cufftDoubleReal), cudaMemcpyHostToDevice);
    fExpN_r.resize(0);
    cufftExecD2Z(plan, d_input, d_output);
    std::vector<std::complex<double>> fExpN_rk(static_cast<size_t>(params->Cn1 * params->Cn2 * params->Cn3 * atoms->nAtom));
    cudaMemcpy(fExpN_rk.data(), d_output, static_cast<size_t>(params->Cn1 * params->Cn2 * params->Cn3 * atoms->nAtom) * sizeof(cufftDoubleComplex), cudaMemcpyDeviceToHost);

    cudaMemcpy(d_input, fExpN_i.data(), params->n * atoms->nAtom * sizeof(cufftDoubleReal), cudaMemcpyHostToDevice);
    fExpN_i.resize(0);
    cufftExecD2Z(plan, d_input, d_output);
    std::vector<std::complex<double>> fExpN_ik(params->Cn1 * params->Cn2 * params->Cn3 * atoms->nAtom);
    cudaMemcpy(fExpN_ik.data(), d_output, params->Cn1 * params->Cn2 * params->Cn3 * atoms->nAtom * sizeof(cufftDoubleComplex), cudaMemcpyDeviceToHost);   

    cudaMemcpy(d_input, fU1ExpN_r.data(), params->n * atoms->nAtom * sizeof(cufftDoubleReal), cudaMemcpyHostToDevice);
    fU1ExpN_r.resize(0);
    cufftExecD2Z(plan, d_input, d_output);
    std::vector<std::complex<double>> fU1ExpN_rk(params->Cn1 * params->Cn2 * params->Cn3 * atoms->nAtom);
    cudaMemcpy(fU1ExpN_rk.data(), d_output, params->Cn1 * params->Cn2 * params->Cn3 * atoms->nAtom * sizeof(cufftDoubleComplex), cudaMemcpyDeviceToHost);

    cudaMemcpy(d_input, fU1ExpN_i.data(), params->n * atoms->nAtom * sizeof(cufftDoubleReal), cudaMemcpyHostToDevice);
    fU1ExpN_i.resize(0);
    cufftExecD2Z(plan, d_input, d_output);
    std::vector<std::complex<double>> fU1ExpN_ik(params->Cn1 * params->Cn2 * params->Cn3 * atoms->nAtom);
    cudaMemcpy(fU1ExpN_ik.data(), d_output, params->Cn1* params->Cn2* params->Cn3* atoms->nAtom * sizeof(cufftDoubleComplex), cudaMemcpyDeviceToHost);

    cudaMemcpy(d_input, fU2ExpN_r.data(), params->n * atoms->nAtom * sizeof(cufftDoubleReal), cudaMemcpyHostToDevice);
    fU2ExpN_r.resize(0);
    cufftExecD2Z(plan, d_input, d_output);
    std::vector<std::complex<double>> fU2ExpN_rk(params->Cn1 * params->Cn2 * params->Cn3 * atoms->nAtom);
    cudaMemcpy(fU2ExpN_rk.data(), d_output, params->Cn1 * params->Cn2 * params->Cn3 * atoms->nAtom * sizeof(cufftDoubleComplex), cudaMemcpyDeviceToHost);

    cudaMemcpy(d_input, fU2ExpN_i.data(), params->n * atoms->nAtom * sizeof(cufftDoubleReal), cudaMemcpyHostToDevice);
    fU2ExpN_i.resize(0);
    cufftExecD2Z(plan, d_input, d_output);
    std::vector<std::complex<double>> fU2ExpN_ik(params->Cn1 * params->Cn2 * params->Cn3 * atoms->nAtom);
    cudaMemcpy(fU2ExpN_ik.data(), d_output, params->Cn1 * params->Cn2 * params->Cn3 * atoms->nAtom * sizeof(cufftDoubleComplex), cudaMemcpyDeviceToHost);

    cudaMemcpy(d_input, fU3ExpN_r.data(), params->n * atoms->nAtom * sizeof(cufftDoubleReal), cudaMemcpyHostToDevice);
    fU3ExpN_r.resize(0);
    cufftExecD2Z(plan, d_input, d_output);
    std::vector<std::complex<double>> fU3ExpN_rk(params->Cn1 * params->Cn2 * params->Cn3 * atoms->nAtom);
    cudaMemcpy(fU3ExpN_rk.data(), d_output, params->Cn1 * params->Cn2 * params->Cn3 * atoms->nAtom * sizeof(cufftDoubleComplex), cudaMemcpyDeviceToHost);

    cudaMemcpy(d_input, fU3ExpN_i.data(), params->n * atoms->nAtom * sizeof(cufftDoubleReal), cudaMemcpyHostToDevice);
    fU3ExpN_i.resize(0);
    cufftExecD2Z(plan, d_input, d_output);
    std::vector<std::complex<double>> fU3ExpN_ik(params->Cn1 * params->Cn2 * params->Cn3 * atoms->nAtom);
    cudaMemcpy(fU3ExpN_ik.data(), d_output, params->Cn1* params->Cn2* params->Cn3* atoms->nAtom * sizeof(cufftDoubleComplex), cudaMemcpyDeviceToHost);

    cufftDestroy(plan);
    cudaFree(d_input);
    cudaFree(d_output);

    for (int i = 0; i < params->Cn1 * params->Cn2 * params->Cn3 * atoms->nAtom; i++) {
        fExpN_rk[i] *= normalization;
        fU1ExpN_rk[i] *= normalization;
        fU2ExpN_rk[i] *= normalization;
        fU3ExpN_rk[i] *= normalization;

        fExpN_ik[i] *= normalization;
        fU1ExpN_ik[i] *= normalization;
        fU2ExpN_ik[i] *= normalization;
        fU3ExpN_ik[i] *= normalization;
    }

    // Calculate intensities
    std::vector<double> IA(params->Cn1 * params->Cn2 * params->Cn3, 0.0);
    std::vector<double> IB(params->Cn1 * params->Cn2 * params->Cn3, 0.0);

    std::vector<std::complex<double>> AmpA(params->Cn1 * params->Cn2 * params->Cn3, { 0.0, 0.0 });
    std::vector<std::complex<double>> AmpB(params->Cn1 * params->Cn2 * params->Cn3, { 0.0, 0.0 });

    for (int i = 0; i < params->Cn1; i++) {
        for (int j = 0; j < params->Cn2; j++) {
            for (int k = 0; k < params->Cn3; k++) {
                int idx = i * params->Cn2 * params->Cn3 + j * params->Cn3 + k;

                for (int a = 0; a < atoms->nAtom; a++) {
                    int aidx = a * params->Cn1 * params->Cn2 * params->Cn3 + idx;

                    std::complex<double> temp_A = fExpN_rk[aidx] + std::complex<double>{0.0, 1.0} *fExpN_ik[aidx];
                    temp_A -= std::complex<double>{0.0, 1.0} *
                        (mqkA1[idx] * fU1ExpN_rk[aidx] + mqkA2[idx] * fU2ExpN_rk[aidx] + mqkA3[idx] * fU3ExpN_rk[aidx]);
                    temp_A += (mqkA1[idx] * fU1ExpN_ik[aidx] + mqkA2[idx] * fU2ExpN_ik[aidx] + mqkA3[idx] * fU3ExpN_ik[aidx]);

                    std::complex<double> temp_B = std::conj(fExpN_rk[aidx]) + std::complex<double>{0.0, 1.0} * std::conj(fExpN_ik[aidx]);
                    temp_B -= std::complex<double>{0.0, 1.0} *
                        (mqkB1[idx] * std::conj(fU1ExpN_rk[aidx]) + mqkB2[idx] * std::conj(fU2ExpN_rk[aidx]) + mqkB3[idx] * std::conj(fU3ExpN_rk[aidx]));
                    temp_B += (mqkB1[idx] * std::conj(fU1ExpN_ik[aidx]) + mqkB2[idx] * std::conj(fU2ExpN_ik[aidx]) + mqkB3[idx] * std::conj(fU3ExpN_ik[aidx]));

                    AmpA[idx] += temp_A * fA[aidx];
                    AmpB[idx] += temp_B * fB[aidx];
                }

                AmpA[idx] *= kfExpA[idx];
                AmpB[idx] *= kfExpB[idx];
            }
        }
    }

    for (int i = 0; i < params->Cn1 * params->Cn2 * params->Cn3; i++) {
        IA[i] = std::max(std::abs(AmpA[i]) * std::abs(AmpA[i]), 1e-50);
        IB[i] = std::max(std::abs(AmpB[i]) * std::abs(AmpB[i]), 1e-50);
    }

    ArrayFourierToRegular(IA, IB, IDiffr, params);
}

void diffraction::ArrayFourierToRegular(const std::vector<double> A, const std::vector<double> B, std::vector<double> &out, sizeContext* params) {

    std::vector<double> temp(params->n, 0.0);
    for (int i = 0; i < params->Cn1; i++) {
        for (int j = 0; j < params->Cn2; j++) {
            for (int k = 0; k < params->Cn3; k++) {
                
                int kk = ((params->nz - k) % params->nz);
                int jj = ((params->ny - j) % params->ny);
                int ii = ((params->nx - i) % params->nx);

                int temp_idx = i * params->Rn2 * params->Rn3 + j * params->Rn3 + k;
                int fourier_idx = i * params->Cn2 * params->Cn3 + j * params->Cn3 + k;
                
                temp[temp_idx] = A[fourier_idx];
                
                if (kk >= params->Cn3 && kk < params->nz && jj < params->ny && ii < params->nx) {
                    int b_idx = ii * params->Rn2 * params->Rn3 + jj * params->Rn3 + kk;
                    temp[b_idx] = B[fourier_idx];
                }
            }
        }
    }

    for (int i = 0; i < params->Rn1; i++) {
        for (int j = 0; j < params->Rn2; j++) {
            for (int k = 0; k < params->Rn3; k++) {
                int ii = ((i + params->nx/2) % params->nx);
                int jj = ((j + params->ny/2) % params->ny);
                int kk = ((k + params->nz/2) % params->nz);

                int out_idx = ii * params->Rn2 * params->Rn3 + jj * params->Rn3 + kk;
                int temp_idx = i * params->Rn2 * params->Rn3 + j * params->Rn3 + k;
                
                out[out_idx] = temp[temp_idx];
            }
        }
    }
}

#endif