#ifndef DIFFRACTION_SETUP_HPP_
#define DIFFRACTION_SETUP_HPP_

#include <cufft.h>
#include <complex>

struct diffraction {
    double VCell;
    double bC[3][3];
    int hM1, hM2, hM3;
    double G10, G20, G30;
    double qC10, qC20, qC30;

    double* mqk1_h, * mqk2_h, * mqk3_h;
    double* mqk1_3, * mqk2_3, * mqk3_3;

    double* mqkA1_h, * mqkA2_h, * mqkA3_h;
    double* mqkA1, * mqkA2, * mqkA3;

    double* mqkB1_h, * mqkB2_h, * mqkB3_h;
    double* mqkB1, * mqkB2, * mqkB3;

    cufftDoubleComplex** fA_h, **fB_h, **fA, **fB;
    cufftDoubleComplex* kfExpA, * kfExpB, * kfExpA_d, * kfExpB_d;

    double *dmqk1, *dmqk2, *dmqk3;

    double** DRAtom;
    int nAtom;

    void diffraction_setup(sizeContext* params, atomList* atoms);

    void diffraction_calc(sizeContext* params, atomList* atoms,
                          std::vector<double>& IDiffr, std::vector<double>& region,
                          std::vector<double>& QCenter, std::vector<double>& oPhase,
                          std::vector<double>& oStruc, std::vector<double>& u);

    void cleanup() {
        delete[] mqkA1_h;
        delete[] mqkA2_h;
        delete[] mqkA3_h;

        delete[] mqkB1_h;
        delete[] mqkB2_h;
        delete[] mqkB3_h;

        delete[] mqk1_h;
        delete[] mqk2_h;
        delete[] mqk3_h;

        cudaFree(mqkA1);
        cudaFree(mqkA2);
        cudaFree(mqkA3);

        cudaFree(mqkB1);
        cudaFree(mqkB2);
        cudaFree(mqkB3);

        cudaFree(mqk1_3);
        cudaFree(mqk2_3);
        cudaFree(mqk3_3);

        delete[] kfExpA;
        delete[] kfExpB;

        cudaFree(kfExpA_d);
        cudaFree(kfExpB_d);

        for (int i = 0; i < nAtom; i++) {
            delete[] DRAtom[i];
            delete[] fA_h[i];
            delete[] fB_h[i];
            cudaFree(fA);
            cudaFree(fB);
        }

        delete[] DRAtom;
        delete[] fA_h;
        delete[] fB_h;
        delete[] fA;
        delete[] fB;
    }
};

void diffraction::diffraction_setup(sizeContext* params, atomList* atoms) {
    VCell = atoms->aC[0][0] * atoms->aC[1][1] * atoms->aC[2][2] + \
        atoms->aC[0][1] * atoms->aC[1][2] * atoms->aC[2][0] + \
        atoms->aC[0][2] * atoms->aC[1][0] * atoms->aC[2][1] - \
        atoms->aC[0][0] * atoms->aC[1][2] * atoms->aC[2][1] - \
        atoms->aC[0][1] * atoms->aC[1][0] * atoms->aC[2][2] - \
        atoms->aC[0][2] * atoms->aC[1][1] * atoms->aC[2][0];


    bC[0][0] = 2 * M_PI / VCell * (atoms->aC[1][1] * atoms->aC[2][2] - atoms->aC[1][2] * atoms->aC[2][1]);
    bC[0][1] = 2 * M_PI / VCell * (atoms->aC[1][2] * atoms->aC[2][0] - atoms->aC[1][0] * atoms->aC[2][2]);
    bC[0][2] = 2 * M_PI / VCell * (atoms->aC[1][0] * atoms->aC[2][1] - atoms->aC[1][1] * atoms->aC[2][0]);

    bC[1][0] = 2 * M_PI / VCell * (atoms->aC[0][2] * atoms->aC[2][1] - atoms->aC[0][1] * atoms->aC[2][2]);
    bC[1][1] = 2 * M_PI / VCell * (atoms->aC[0][0] * atoms->aC[2][2] - atoms->aC[0][2] * atoms->aC[2][0]);
    bC[1][2] = 2 * M_PI / VCell * (atoms->aC[0][1] * atoms->aC[2][0] - atoms->aC[0][0] * atoms->aC[2][1]);

    bC[2][0] = 2 * M_PI / VCell * (atoms->aC[0][1] * atoms->aC[1][2] - atoms->aC[0][2] * atoms->aC[1][1]);
    bC[2][1] = 2 * M_PI / VCell * (atoms->aC[0][2] * atoms->aC[1][0] - atoms->aC[0][0] * atoms->aC[1][2]);
    bC[2][2] = 2 * M_PI / VCell * (atoms->aC[0][0] * atoms->aC[1][1] - atoms->aC[0][1] * atoms->aC[1][0]);

    hM1 = std::lround((atoms->q00[0] * atoms->aC[0][0] + atoms->q00[1] * atoms->aC[0][1] + atoms->q00[2] * atoms->aC[0][2]) / (2 * M_PI));
    hM2 = std::lround((atoms->q00[0] * atoms->aC[1][0] + atoms->q00[1] * atoms->aC[1][1] + atoms->q00[2] * atoms->aC[1][2]) / (2 * M_PI));
    hM3 = std::lround((atoms->q00[0] * atoms->aC[2][0] + atoms->q00[1] * atoms->aC[2][1] + atoms->q00[2] * atoms->aC[2][2]) / (2 * M_PI));

    G10 = hM1 * bC[0][0] + hM2 * bC[1][0] + hM3 * bC[2][0];
    G20 = hM1 * bC[0][1] + hM2 * bC[1][1] + hM3 * bC[2][1];
    G30 = hM1 * bC[0][2] + hM2 * bC[1][2] + hM3 * bC[2][2];

    qC10 = G10;
    qC20 = G20;
    qC30 = G30;

    mqk1_h = new double[params->nx];
    mqk2_h = new double[params->ny];
    mqk3_h = new double[params->nz];

    mqkA1_h = new double[params->nx];
    mqkA2_h = new double[params->ny];
    mqkA3_h = new double[params->nz];

    mqkB1_h = new double[params->nx];
    mqkB2_h = new double[params->ny];
    mqkB3_h = new double[params->nz];

    for (int i = 0; i < params->nx; i++) {
        if (i < params->nx / 2) {
            mqk1_h[i] = (i - params->nx);
        }
        else {
            mqk1_h[i] = i;
        }

        mqk1_h[i] *= 2.0 * M_PI / params->lx;

        mqkA1_h[i] = mqk1_h[i] + qC10;
        mqkB1_h[i] = -mqk1_h[i] + qC10;
    }

    for (int i = 0; i < params->ny; i++) {
        if (i < params->ny / 2) {
            mqk2_h[i] = (i - params->ny);
        }
        else {
            mqk2_h[i] = i;
        }

        mqk2_h[i] *= 2.0 * M_PI / params->ly;

        mqkA2_h[i] = mqk2_h[i] + qC20;
        mqkB2_h[i] = -mqk2_h[i] + qC20;
    }

    for (int i = 0; i < params->nz; i++) {
        if (i < params->nz / 2) {
            mqk3_h[i] = (i - params->nz);
        }
        else {
            mqk3_h[i] = i;
        }

        mqk3_h[i] *= 2.0 * M_PI / params->lz;

        mqkA3_h[i] = mqk3_h[i] + qC30;
        mqkB3_h[i] = -mqk3_h[i] + qC30;
    }

    if (params->nx % 2 == 0) {
        mqkB1_h[params->nx / 2] = mqk3_h[params->nx / 2] + qC10;
    }

    if (params->ny % 2 == 0) {
        mqkB2_h[params->ny / 2] = mqk3_h[params->ny / 2] + qC20;
    }

    if (params->nz % 2 == 0) {
        mqkB3_h[params->nz / 2] = mqk3_h[params->nz / 2] + qC30;
    }

    kfExpA = new cufftDoubleComplex[params->nx * params->ny * params->nz];
    kfExpB = new cufftDoubleComplex[params->nx * params->ny * params->nz];

    cudaMalloc((void**)&kfExpA_d, params->nx * params->ny * params->nz * sizeof(cufftDoubleComplex));
    cudaMalloc((void**)&kfExpB_d, params->nx * params->ny * params->nz * sizeof(cufftDoubleComplex));

    for (int i = 0; i < params->nx; i++) {
        for (int j = 0; j < params->ny; j++) {
            for (int k = 0; k < params->nz; k++) {
                int idx = i * params->ny * params->nz + j * params->nz + k;

                std::complex kTemp = { 0.0, (mqkA1_h[i] - G10) * params->dx };
                std::complex temp1 = { kTemp.real() / 2.0, kTemp.imag() / 2.0 };
                std::complex temp2 = { -kTemp.real() / 2.0, -kTemp.imag() / 2.0 };
                std::complex<double> kfTemp1 = 1.0 / kTemp * (exp(temp1) - exp(temp2));
                if (std::abs(kTemp) <= 1e-8) kfTemp1 = { 1.0, 0.0 };

                kTemp = { 0.0, (mqkA2_h[j] - G20) * params->dy };
                temp1 = { kTemp.real() / 2.0, kTemp.imag() / 2.0 };
                temp2 = { -kTemp.real() / 2.0, -kTemp.imag() / 2.0 };
                std::complex<double> kfTemp2 = 1.0 / kTemp * (exp(temp1) - exp(temp2));
                if (std::abs(kTemp) <= 1e-8) kfTemp2 = { 1.0, 0.0 };

                kTemp = { 0.0, (mqkA3_h[k] - G30) * params->dz };
                temp1 = { kTemp.real() / 2.0, kTemp.imag() / 2.0 };
                temp2 = { -kTemp.real() / 2.0, -kTemp.imag() / 2.0 };
                std::complex<double> kfTemp3 = 1.0 / kTemp * (exp(temp1) - exp(temp2));
                if (std::abs(kTemp) <= 1e-8) kfTemp3 = { 1.0, 0.0 };

                std::complex<double> kfTemp = kfTemp1 * kfTemp2 * kfTemp3;
                kfExpA[idx] = { kfTemp.real(), kfTemp.imag() };

                kTemp = { 0.0, (mqkB1_h[i] - G10) * params->dx };
                temp1 = { kTemp.real() / 2.0, kTemp.imag() / 2.0 };
                temp2 = { -kTemp.real() / 2.0, -kTemp.imag() / 2.0 };
                kfTemp1 = 1.0 / kTemp * (exp(temp1) - exp(temp2));
                if (std::abs(kTemp) <= 1e-8) kfTemp1 = { 1.0, 0.0 };

                kTemp = { 0.0, (mqkB2_h[j] - G20) * params->dy };
                temp1 = { kTemp.real() / 2.0, kTemp.imag() / 2.0 };
                temp2 = { -kTemp.real() / 2.0, -kTemp.imag() / 2.0 };
                kfTemp2 = 1.0 / kTemp * (exp(temp1) - exp(temp2));
                if (std::abs(kTemp) <= 1e-8) kfTemp2 = { 1.0, 0.0 };

                kTemp = { 0.0, (mqkB3_h[k] - G30) * params->dz };
                temp1 = { kTemp.real() / 2.0, kTemp.imag() / 2.0 };
                temp2 = { -kTemp.real() / 2.0, -kTemp.imag() / 2.0 };
                kfTemp3 = 1.0 / kTemp * (exp(temp1) - exp(temp2));
                if (std::abs(kTemp) <= 1e-8) kfTemp3 = { 1.0, 0.0 };

                kfTemp = kfTemp1 * kfTemp2 * kfTemp3;
                kfExpB[idx] = { kfTemp.real(), kfTemp.imag() };
            }
        }
    }

    cudaMemcpy(kfExpA_d, kfExpA, params->nx * params->ny * params->nz * sizeof(cufftDoubleComplex), cudaMemcpyHostToDevice);
    cudaMemcpy(kfExpB_d, kfExpB, params->nx * params->ny * params->nz * sizeof(cufftDoubleComplex), cudaMemcpyHostToDevice);

    this->nAtom = atoms->nAtom;

    fA_h = new cufftDoubleComplex * [this->nAtom];
    fB_h = new cufftDoubleComplex * [this->nAtom];

    fA = new cufftDoubleComplex * [this->nAtom];
    fB = new cufftDoubleComplex * [this->nAtom];

    DRAtom = new double* [this->nAtom];
    for (int i = 0; i < this->nAtom; i++) {
        DRAtom[i] = new double[3];

        DRAtom[i][0] = atoms->xAtom[i][0] * atoms->aC[0][0] + \
            atoms->xAtom[i][1] * atoms->aC[1][0] + \
            atoms->xAtom[i][2] * atoms->aC[2][0];
        DRAtom[i][1] = atoms->xAtom[i][0] * atoms->aC[0][1] + \
            atoms->xAtom[i][1] * atoms->aC[1][1] + \
            atoms->xAtom[i][2] * atoms->aC[2][1];
        DRAtom[i][2] = atoms->xAtom[i][0] * atoms->aC[0][2] + \
            atoms->xAtom[i][1] * atoms->aC[1][2] + \
            atoms->xAtom[i][2] * atoms->aC[2][2];

        fA_h[i] = new cufftDoubleComplex[params->nx * params->ny * params->nz];
        fB_h[i] = new cufftDoubleComplex[params->nx * params->ny * params->nz];

        cudaMalloc((void**)&fA[i], params->nx * params->ny * params->nz * sizeof(cufftDoubleComplex));
        cudaMalloc((void**)&fB[i], params->nx * params->ny * params->nz * sizeof(cufftDoubleComplex));

        for (int x = 0; x < params->nx; x++) {
            for (int y = 0; y < params->ny; y++) {
                for (int z = 0; z < params->nz; z++) {
                    int idx = x * params->ny * params->nz + y * params->nz + z;

                    std::complex<double> temp = { 0.0, -mqkA1_h[x] * DRAtom[i][0] - mqkA2_h[y] * DRAtom[i][1] - mqkA3_h[z] * DRAtom[i][2] };
                    std::complex<double> exptemp = exp(temp) / VCell;
                    fA_h[i][idx] = { exptemp.real(), exptemp.imag() };

                    temp = { 0.0, -mqkB1_h[x] * DRAtom[i][0] - mqkB2_h[y] * DRAtom[i][1] - mqkB3_h[z] * DRAtom[i][2] };
                    exptemp = exp(temp) / VCell;
                    fB_h[i][idx] = { exptemp.real(), exptemp.imag() };
                }
            }
        }

        cudaMemcpy(fA[i], fA_h[i], params->nx * params->ny * params->nz * sizeof(cufftDoubleComplex), cudaMemcpyHostToDevice);
        cudaMemcpy(fB[i], fB_h[i], params->nx * params->ny * params->nz * sizeof(cufftDoubleComplex), cudaMemcpyHostToDevice);
    }

    // call ArrayFourierToRegular(mqkA1,mqkB1,mqk1_out,trans)
    // call ArrayFourierToRegular(mqkA2,mqkB2,mqk2_out,trans)
    // call ArrayFourierToRegular(mqkA3,mqkB3,mqk3_out,trans)

    // DQ
    // dmqk1 = mqk1_out - qC10
    // dmqk2 = mqk2_out - qC20
    // dmqk3 = mqk3_out - qC30   !dmqk is also /Delta q

    for (int i = 0; i < params->nx; i++) {
        mqkA1_h[i] = mqkA1_h[i] - qC10;
        mqkB1_h[i] = mqkB1_h[i] - qC10;
    }

    for (int i = 0; i < params->ny; i++) {
        mqkA2_h[i] = mqkA2_h[i] - qC20;
        mqkB2_h[i] = mqkB2_h[i] - qC20;
    }

    for (int i = 0; i < params->nz; i++) {
        mqkA3_h[i] = mqkA3_h[i] - qC30;
        mqkB3_h[i] = mqkB3_h[i] - qC30;
    }

    cudaMalloc((void**)&mqk1_3, params->nx * sizeof(double));
    cudaMalloc((void**)&mqk2_3, params->ny * sizeof(double));
    cudaMalloc((void**)&mqk3_3, params->nz * sizeof(double));

    cudaMalloc((void**)&mqkA1, params->nx * sizeof(double));
    cudaMalloc((void**)&mqkA2, params->ny * sizeof(double));
    cudaMalloc((void**)&mqkA3, params->nz * sizeof(double));

    cudaMalloc((void**)&mqkB1, params->nx * sizeof(double));
    cudaMalloc((void**)&mqkB2, params->ny * sizeof(double));
    cudaMalloc((void**)&mqkB3, params->nz * sizeof(double));

    cudaMemcpy(mqk1_3, mqk1_h, params->nx * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(mqk2_3, mqk2_h, params->ny * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(mqk3_3, mqk3_h, params->nz * sizeof(double), cudaMemcpyHostToDevice);

    cudaMemcpy(mqkA1, mqkA1_h, params->nx * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(mqkA2, mqkA2_h, params->ny * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(mqkA3, mqkA3_h, params->nz * sizeof(double), cudaMemcpyHostToDevice);

    cudaMemcpy(mqkB1, mqkB1_h, params->nx * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(mqkB2, mqkB2_h, params->ny * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(mqkB3, mqkB3_h, params->nz * sizeof(double), cudaMemcpyHostToDevice);
}

void diffraction::diffraction_calc(sizeContext* params, atomList* atoms,
                          std::vector<double>& IDiffr, std::vector<double>& region,
                          std::vector<double>& QCenter, std::vector<double>& oPhase,
                          std::vector<double>& oStruc, std::vector<double>& u) {
    double G1 = G10 - qC10 * params->strainAvg[0] - qC20 * params->strainAvg[5] - qC30 * params->strainAvg[4];
    double G2 = G20 - qC10 * params->strainAvg[5] - qC20 * params->strainAvg[1] - qC30 * params->strainAvg[3];
    double G3 = G30 - qC10 * params->strainAvg[4] - qC20 * params->strainAvg[3] - qC30 * params->strainAvg[2];

    double qC1 = qC10 + G1 - G10;
    double qC2 = qC20 + G2 - G20;
    double qC3 = qC30 + G3 - G30;

    // mqk1_out = dmqk1 + qC1
    // mqk2_out = dmqk2 + qC2
    // mqk3_out = dmqk3 + qC3    

    std::vector<double> q0Gr(params->nx * params->ny * params->nz, 0.0);
    for (int i = 0; i < params->nx; i++) {
        for (int j = 0; j < params->ny; j++) {
            for (int k = 0; k < params->nz; k++) {
                int idx = i * params->ny * params->nz + j * params->nz + k;
                q0Gr[idx] += (qC1 - G1) * i * params->dx;
                q0Gr[idx] += (qC2 - G2) * j * params->dy;
                q0Gr[idx] += (qC3 - G3) * k * params->dz;
            }
        }
    }

    for (int n = 0; n < atoms->nAtom; n++) {
        for (int m = 0; m < params->nPhase; m++) {
    // fExpN = 0.d0
    // fU1ExpN = 0.d0
    // fU2ExpN = 0.d0
    // fU3ExpN = 0.d0

    // do n = 1,nAtom
    //   do m = 1,nPhase
    //     uNM1(:,:,:,n) = u1(:,:,:)
    //     uNM2(:,:,:,n) = u2(:,:,:)
    //     uNM3(:,:,:,n) = u3(:,:,:)
    //     do k = 1,nStruc
    //       uNM1(:,:,:,n) = uNM1(:,:,:,n) + su1(k,m,n)*os(k,:,:,:)
    //       uNM2(:,:,:,n) = uNM2(:,:,:,n) + su2(k,m,n)*os(k,:,:,:)
    //       uNM3(:,:,:,n) = uNM3(:,:,:,n) + su3(k,m,n)*os(k,:,:,:)
    //     enddo
    //     fExpNM(:,:,:,n) = cdexp( -(0,1)* (qC1*uNM1(:,:,:,n) + qC2*uNM2(:,:,:,n) + qC3*uNM3(:,:,:,n) + q0Gr) ) *fAtom(m,n)*op(m,:,:,:)*region
    //       fExpN(:,:,:,n) =   fExpN(:,:,:,n) + fExpNM(:,:,:,n)
    //     fU1ExpN(:,:,:,n) = fU1ExpN(:,:,:,n) + uNM1(:,:,:,n) * fExpNM(:,:,:,n)
    //     fU2ExpN(:,:,:,n) = fU2ExpN(:,:,:,n) + uNM2(:,:,:,n) * fExpNM(:,:,:,n)
    //     fU3ExpN(:,:,:,n) = fU3ExpN(:,:,:,n) + uNM3(:,:,:,n) * fExpNM(:,:,:,n)
    //   enddo
    // enddo

    // do n = 1,nAtom
    //   tempR = real(  fExpN(:,:,:,n));  call mupro_fft_forward(tempR,tempC);    fExpNRk(:,:,:,n) = tempC/(nx*ny*nz)
    //   tempR = imag(  fExpN(:,:,:,n));  call mupro_fft_forward(tempR,tempC);    fExpNIk(:,:,:,n) = tempC/(nx*ny*nz)
    //   tempR = real(fU1ExpN(:,:,:,n));  call mupro_fft_forward(tempR,tempC);  fU1ExpNRk(:,:,:,n) = tempC/(nx*ny*nz)
    //   tempR = imag(fU1ExpN(:,:,:,n));  call mupro_fft_forward(tempR,tempC);  fU1ExpNIk(:,:,:,n) = tempC/(nx*ny*nz)
    //   tempR = real(fU2ExpN(:,:,:,n));  call mupro_fft_forward(tempR,tempC);  fU2ExpNRk(:,:,:,n) = tempC/(nx*ny*nz)
    //   tempR = imag(fU2ExpN(:,:,:,n));  call mupro_fft_forward(tempR,tempC);  fU2ExpNIk(:,:,:,n) = tempC/(nx*ny*nz)
    //   tempR = real(fU3ExpN(:,:,:,n));  call mupro_fft_forward(tempR,tempC);  fU3ExpNRk(:,:,:,n) = tempC/(nx*ny*nz)
    //   tempR = imag(fU3ExpN(:,:,:,n));  call mupro_fft_forward(tempR,tempC);  fU3ExpNIk(:,:,:,n) = tempC/(nx*ny*nz)
    // enddo

    // AmpA = (0.,0.)
    // AmpB = (0.,0.)
    // do n=1,nAtom
    //   AmpA = AmpA + fA(:,:,:,n) * (                          fExpNRk(:,:,:,n)  + (0,1)*          fExpNIk(:,:,:,n)   &
    //                                 -(0,1)* (mqkA1*       fU1ExpNRk(:,:,:,n)  + mqkA2*       fU2ExpNRk(:,:,:,n)  + mqkA3*       fU3ExpNRk(:,:,:,n)) &
    //                                 +       (mqkA1*       fU1ExpNIk(:,:,:,n)  + mqkA2*       fU2ExpNIk(:,:,:,n)  + mqkA3*       fU3ExpNIk(:,:,:,n))  )
    //   AmpB = AmpB + fB(:,:,:,n) * (                 dconjg(  fExpNRk(:,:,:,n)) + (0,1)* dconjg(  fExpNIk(:,:,:,n))  &
    //                                 -(0,1)* (mqkB1*dconjg(fU1ExpNRk(:,:,:,n)) + mqkB2*dconjg(fU2ExpNRk(:,:,:,n)) + mqkB3*dconjg(fU3ExpNRk(:,:,:,n)))&
    //                                 +       (mqkB1*dconjg(fU1ExpNIk(:,:,:,n)) + mqkB2*dconjg(fU2ExpNIk(:,:,:,n)) + mqkB3*dconjg(fU3ExpNIk(:,:,:,n))) )
    // enddo

    // AmpA = AmpA * kfExpA
    // AmpB = AmpB * kfExpB

    // IA = max(cdabs(AmpA)**2,1.d-50)
    // IB = max(cdabs(AmpB)**2,1.d-50)

    // call ArrayFourierToRegular(IA,IB,I_out,trans)            
        }
    }
}

#endif