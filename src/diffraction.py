"""
Diffraction calculation class
"""

import numpy as np
import math

# Try to import CuPy for GPU acceleration
try:
    import cupy as cp
    from cupyx.scipy.fft import rfftn as cp_rfftn
    HAS_CUPY = True
    print("CuPy detected - using GPU acceleration for FFTs")
except ImportError:
    cp = np
    HAS_CUPY = False
    from numpy.fft import rfftn as cp_rfftn
    print("CuPy not available - using CPU (NumPy) for FFTs")


class Diffraction:
    """Container for diffraction calculation parameters and methods"""
    
    def __init__(self):
        self.VCell = 0.0
        self.bC = np.zeros((3, 3))
        
        self.hM1 = 0
        self.hM2 = 0
        self.hM3 = 0
        
        self.G10 = 0.0
        self.G20 = 0.0
        self.G30 = 0.0
        
        self.qC10 = 0.0
        self.qC20 = 0.0
        self.qC30 = 0.0
        
        self.qC1 = 0.0
        self.qC2 = 0.0
        self.qC3 = 0.0
        
        # Wave vectors
        self.mqk1 = None
        self.mqk2 = None
        self.mqk3 = None
        
        self.mqk1_out = None
        self.mqk2_out = None
        self.mqk3_out = None
        
        self.mqkA1 = None
        self.mqkA2 = None
        self.mqkA3 = None
        
        self.mqkB1 = None
        self.mqkB2 = None
        self.mqkB3 = None
        
        # Form factors and exponentials
        self.fA = None
        self.fB = None
        self.kfExpA = None
        self.kfExpB = None
        
        # Displacement arrays
        self.dmqk1 = None
        self.dmqk2 = None
        self.dmqk3 = None
        
        self.DRAtom = None
        self.nAtom = 0
    
    def diffraction_setup(self, params, atoms):
        """
        Setup diffraction calculation parameters.
        
        Args:
            params: SizeContext object
            atoms: AtomList object
        """
        # Calculate unit cell volume
        self.VCell = (atoms.aC[0, 0] * atoms.aC[1, 1] * atoms.aC[2, 2] +
                      atoms.aC[0, 1] * atoms.aC[1, 2] * atoms.aC[2, 0] +
                      atoms.aC[0, 2] * atoms.aC[1, 0] * atoms.aC[2, 1] -
                      atoms.aC[0, 0] * atoms.aC[1, 2] * atoms.aC[2, 1] -
                      atoms.aC[0, 1] * atoms.aC[1, 0] * atoms.aC[2, 2] -
                      atoms.aC[0, 2] * atoms.aC[1, 1] * atoms.aC[2, 0])
        
        # Calculate reciprocal lattice vectors
        self.bC[0, 0] = 2.0 * math.pi / self.VCell * (atoms.aC[1, 1] * atoms.aC[2, 2] - atoms.aC[1, 2] * atoms.aC[2, 1])
        self.bC[0, 1] = 2.0 * math.pi / self.VCell * (atoms.aC[1, 2] * atoms.aC[2, 0] - atoms.aC[1, 0] * atoms.aC[2, 2])
        self.bC[0, 2] = 2.0 * math.pi / self.VCell * (atoms.aC[1, 0] * atoms.aC[2, 1] - atoms.aC[1, 1] * atoms.aC[2, 0])
        
        self.bC[1, 0] = 2.0 * math.pi / self.VCell * (atoms.aC[0, 2] * atoms.aC[2, 1] - atoms.aC[0, 1] * atoms.aC[2, 2])
        self.bC[1, 1] = 2.0 * math.pi / self.VCell * (atoms.aC[0, 0] * atoms.aC[2, 2] - atoms.aC[0, 2] * atoms.aC[2, 0])
        self.bC[1, 2] = 2.0 * math.pi / self.VCell * (atoms.aC[0, 1] * atoms.aC[2, 0] - atoms.aC[0, 0] * atoms.aC[2, 1])
        
        self.bC[2, 0] = 2.0 * math.pi / self.VCell * (atoms.aC[0, 1] * atoms.aC[1, 2] - atoms.aC[0, 2] * atoms.aC[1, 1])
        self.bC[2, 1] = 2.0 * math.pi / self.VCell * (atoms.aC[0, 2] * atoms.aC[1, 0] - atoms.aC[0, 0] * atoms.aC[1, 2])
        self.bC[2, 2] = 2.0 * math.pi / self.VCell * (atoms.aC[0, 0] * atoms.aC[1, 1] - atoms.aC[0, 1] * atoms.aC[1, 0])
        
        # Calculate Miller indices
        self.hM1 = round((atoms.q00[0] * atoms.aC[0, 0] + atoms.q00[1] * atoms.aC[0, 1] + atoms.q00[2] * atoms.aC[0, 2]) / (2 * math.pi))
        self.hM2 = round((atoms.q00[0] * atoms.aC[1, 0] + atoms.q00[1] * atoms.aC[1, 1] + atoms.q00[2] * atoms.aC[1, 2]) / (2 * math.pi))
        self.hM3 = round((atoms.q00[0] * atoms.aC[2, 0] + atoms.q00[1] * atoms.aC[2, 1] + atoms.q00[2] * atoms.aC[2, 2]) / (2 * math.pi))
        
        # Calculate G vectors
        self.G10 = self.hM1 * self.bC[0, 0] + self.hM2 * self.bC[1, 0] + self.hM3 * self.bC[2, 0]
        self.G20 = self.hM1 * self.bC[0, 1] + self.hM2 * self.bC[1, 1] + self.hM3 * self.bC[2, 1]
        self.G30 = self.hM1 * self.bC[0, 2] + self.hM2 * self.bC[1, 2] + self.hM3 * self.bC[2, 2]
        
        self.qC10 = self.G10
        self.qC20 = self.G20
        self.qC30 = self.G30
        
        # Initialize wave vector arrays
        self.mqk1 = np.zeros(params.Cn1 * params.Cn2 * params.Cn3)
        self.mqk2 = np.zeros(params.Cn1 * params.Cn2 * params.Cn3)
        self.mqk3 = np.zeros(params.Cn1 * params.Cn2 * params.Cn3)
        
        self.mqkA1 = np.zeros(params.Cn1 * params.Cn2 * params.Cn3)
        self.mqkA2 = np.zeros(params.Cn1 * params.Cn2 * params.Cn3)
        self.mqkA3 = np.zeros(params.Cn1 * params.Cn2 * params.Cn3)
        
        self.mqkB1 = np.zeros(params.Cn1 * params.Cn2 * params.Cn3)
        self.mqkB2 = np.zeros(params.Cn1 * params.Cn2 * params.Cn3)
        self.mqkB3 = np.zeros(params.Cn1 * params.Cn2 * params.Cn3)
        
        # Setup wave vectors in Fourier space
        for i in range(params.Cn1):
            for j in range(params.Cn2):
                for k in range(params.Cn3):
                    idx = i * params.Cn2 * params.Cn3 + j * params.Cn3 + k
                    
                    self.mqk1[idx] = float(i)
                    if i >= params.nx // 2:
                        self.mqk1[idx] = float(i - params.nx)
                    self.mqk1[idx] *= 2.0 * math.pi / float(params.lx)
                    
                    self.mqkA1[idx] = self.mqk1[idx] + self.qC10
                    self.mqkB1[idx] = -self.mqk1[idx] + self.qC10
                    
                    self.mqk2[idx] = float(j)
                    if j >= params.ny // 2:
                        self.mqk2[idx] = float(j - params.ny)
                    self.mqk2[idx] *= 2.0 * math.pi / float(params.ly)
                    
                    self.mqkA2[idx] = self.mqk2[idx] + self.qC20
                    self.mqkB2[idx] = -self.mqk2[idx] + self.qC20
                    
                    self.mqk3[idx] = float(k)
                    if k >= params.nz // 2:
                        self.mqk3[idx] = float(k - params.nz)
                    self.mqk3[idx] *= 2.0 * math.pi / float(params.lz)
                    
                    self.mqkA3[idx] = self.mqk3[idx] + self.qC30
                    self.mqkB3[idx] = -self.mqk3[idx] + self.qC30
                    
                    # Handle Nyquist frequencies
                    if params.nx % 2 == 0 and i == params.nx // 2:
                        self.mqkB1[idx] = self.mqk1[idx] + self.qC10
                    
                    if params.ny % 2 == 0 and j == params.ny // 2:
                        self.mqkB2[idx] = self.mqk2[idx] + self.qC20
                    
                    if params.nz % 2 == 0 and k == params.nz // 2:
                        self.mqkB3[idx] = self.mqk3[idx] + self.qC30
        
        # Setup form factor exponentials
        self.kfExpA = np.zeros(params.Cn1 * params.Cn2 * params.Cn3, dtype=complex)
        self.kfExpB = np.zeros(params.Cn1 * params.Cn2 * params.Cn3, dtype=complex)
        
        for i in range(params.Cn1):
            for j in range(params.Cn2):
                for k in range(params.Cn3):
                    idx = i * params.Cn2 * params.Cn3 + j * params.Cn3 + k
                    
                    # Calculate kfExpA
                    kTemp = 1j * (self.mqkA1[idx] - self.G10) * params.dx
                    temp1 = kTemp / 2.0
                    temp2 = -kTemp / 2.0
                    kfTemp1 = (np.exp(temp1) - np.exp(temp2)) / kTemp if abs(kTemp) > 1e-8 else 1.0 + 0j
                    
                    kTemp = 1j * (self.mqkA2[idx] - self.G20) * params.dy
                    temp1 = kTemp / 2.0
                    temp2 = -kTemp / 2.0
                    kfTemp2 = (np.exp(temp1) - np.exp(temp2)) / kTemp if abs(kTemp) > 1e-8 else 1.0 + 0j
                    
                    kTemp = 1j * (self.mqkA3[idx] - self.G30) * params.dz
                    temp1 = kTemp / 2.0
                    temp2 = -kTemp / 2.0
                    kfTemp3 = (np.exp(temp1) - np.exp(temp2)) / kTemp if abs(kTemp) > 1e-8 else 1.0 + 0j
                    
                    self.kfExpA[idx] = kfTemp1 * kfTemp2 * kfTemp3
                    
                    # Calculate kfExpB
                    kTemp = 1j * (self.mqkB1[idx] - self.G10) * params.dx
                    temp1 = kTemp / 2.0
                    temp2 = -kTemp / 2.0
                    kfTemp1 = (np.exp(temp1) - np.exp(temp2)) / kTemp if abs(kTemp) > 1e-8 else 1.0 + 0j
                    
                    kTemp = 1j * (self.mqkB2[idx] - self.G20) * params.dy
                    temp1 = kTemp / 2.0
                    temp2 = -kTemp / 2.0
                    kfTemp2 = (np.exp(temp1) - np.exp(temp2)) / kTemp if abs(kTemp) > 1e-8 else 1.0 + 0j
                    
                    kTemp = 1j * (self.mqkB3[idx] - self.G30) * params.dz
                    temp1 = kTemp / 2.0
                    temp2 = -kTemp / 2.0
                    kfTemp3 = (np.exp(temp1) - np.exp(temp2)) / kTemp if abs(kTemp) > 1e-8 else 1.0 + 0j
                    
                    self.kfExpB[idx] = kfTemp1 * kfTemp2 * kfTemp3
        
        # Setup atom form factors
        self.nAtom = atoms.nAtom
        self.fA = np.zeros((params.Cn1 * params.Cn2 * params.Cn3 * self.nAtom), dtype=complex)
        self.fB = np.zeros((params.Cn1 * params.Cn2 * params.Cn3 * self.nAtom), dtype=complex)
        
        self.DRAtom = np.zeros(3 * self.nAtom)
        
        for i in range(self.nAtom):
            # Calculate real space atomic positions
            self.DRAtom[i * 3 + 0] = (atoms.xAtom[i, 0] * atoms.aC[0, 0] +
                                       atoms.xAtom[i, 1] * atoms.aC[1, 0] +
                                       atoms.xAtom[i, 2] * atoms.aC[2, 0])
            self.DRAtom[i * 3 + 1] = (atoms.xAtom[i, 0] * atoms.aC[0, 1] +
                                       atoms.xAtom[i, 1] * atoms.aC[1, 1] +
                                       atoms.xAtom[i, 2] * atoms.aC[2, 1])
            self.DRAtom[i * 3 + 2] = (atoms.xAtom[i, 0] * atoms.aC[0, 2] +
                                       atoms.xAtom[i, 1] * atoms.aC[1, 2] +
                                       atoms.xAtom[i, 2] * atoms.aC[2, 2])
            
            for x in range(params.Cn1):
                for y in range(params.Cn2):
                    for z in range(params.Cn3):
                        idx = x * params.Cn2 * params.Cn3 + y * params.Cn3 + z
                        aidx = idx + i * (params.Cn1 * params.Cn2 * params.Cn3)
                        
                        temp = 1j * (-self.mqkA1[idx] * self.DRAtom[i * 3] -
                                     self.mqkA2[idx] * self.DRAtom[i * 3 + 1] -
                                     self.mqkA3[idx] * self.DRAtom[i * 3 + 2])
                        self.fA[aidx] = np.exp(temp) / self.VCell
                        
                        temp = 1j * (-self.mqkB1[idx] * self.DRAtom[i * 3] -
                                     self.mqkB2[idx] * self.DRAtom[i * 3 + 1] -
                                     self.mqkB3[idx] * self.DRAtom[i * 3 + 2])
                        self.fB[aidx] = np.exp(temp) / self.VCell
        
        # Setup output wave vectors
        self.mqk1_out = np.zeros(params.n)
        self.mqk2_out = np.zeros(params.n)
        self.mqk3_out = np.zeros(params.n)
        
        self.ArrayFourierToRegular(self.mqkA1, self.mqkB1, self.mqk1_out, params)
        self.ArrayFourierToRegular(self.mqkA2, self.mqkB2, self.mqk2_out, params)
        self.ArrayFourierToRegular(self.mqkA3, self.mqkB3, self.mqk3_out, params)
        
        # Adjust wave vectors relative to center
        for i in range(params.Cn1 * params.Cn2 * params.Cn3):
            self.mqkA1[i] = self.mqkA1[i] - self.qC10
            self.mqkB1[i] = self.mqkB1[i] - self.qC10
            
            self.mqkA2[i] = self.mqkA2[i] - self.qC20
            self.mqkB2[i] = self.mqkB2[i] - self.qC20
            
            self.mqkA3[i] = self.mqkA3[i] - self.qC30
            self.mqkB3[i] = self.mqkB3[i] - self.qC30
    
    def diffraction_calc(self, params, atoms, IDiffr, region, QCenter, oPhase, oStruc, u):
        """
        Calculate diffraction intensities.
        
        Args:
            params: SizeContext object
            atoms: AtomList object
            IDiffr: Output intensity array
            region: Region mask array
            QCenter: Q-space center (not used in current implementation)
            oPhase: Phase fraction array
            oStruc: Structural order parameter array
            u: Displacement field array
        """
        # Calculate adjusted G vectors with strain
        G1 = self.G10 - self.qC10 * params.strainAvg[0] - self.qC20 * params.strainAvg[5] - self.qC30 * params.strainAvg[4]
        G2 = self.G20 - self.qC10 * params.strainAvg[5] - self.qC20 * params.strainAvg[1] - self.qC30 * params.strainAvg[3]
        G3 = self.G30 - self.qC10 * params.strainAvg[4] - self.qC20 * params.strainAvg[3] - self.qC30 * params.strainAvg[2]
        
        self.qC1 = self.qC10 + G1 - self.G10
        self.qC2 = self.qC20 + G2 - self.G20
        self.qC3 = self.qC30 + G3 - self.G30
        
        # Calculate phase factor q0·r
        q0Gr = np.zeros(params.n)
        for i in range(params.nx):
            for j in range(params.ny):
                for k in range(params.nz):
                    idx = i * params.ny * params.nz + j * params.nz + k
                    q0Gr[idx] += (self.qC1 - G1) * float(i + 1) * float(params.dx)
                    q0Gr[idx] += (self.qC2 - G2) * float(j + 1) * float(params.dy)
                    q0Gr[idx] += (self.qC3 - G3) * float(k + 1) * float(params.dz)
        
        # Initialize form factor arrays (real and imaginary parts)
        fExpN_r = np.zeros(params.n * atoms.nAtom)
        fU1ExpN_r = np.zeros(params.n * atoms.nAtom)
        fU2ExpN_r = np.zeros(params.n * atoms.nAtom)
        fU3ExpN_r = np.zeros(params.n * atoms.nAtom)
        
        fExpN_i = np.zeros(params.n * atoms.nAtom)
        fU1ExpN_i = np.zeros(params.n * atoms.nAtom)
        fU2ExpN_i = np.zeros(params.n * atoms.nAtom)
        fU3ExpN_i = np.zeros(params.n * atoms.nAtom)
        
        # Calculate form factors with displacements
        for n in range(atoms.nAtom):
            for m in range(params.nPhase):
                uNM1 = np.zeros(params.n)
                uNM2 = np.zeros(params.n)
                uNM3 = np.zeros(params.n)
                
                if u.size == 3 * params.n:
                    uNM1[:] = u[0:params.n]
                    uNM2[:] = u[params.n:2*params.n]
                    uNM3[:] = u[2*params.n:3*params.n]
                
                for k in range(params.nStruc):
                    for i in range(params.n):
                        uNM1[i] += atoms.su1[k, m, n] * oStruc[k * params.n + i]
                        uNM2[i] += atoms.su2[k, m, n] * oStruc[k * params.n + i]
                        uNM3[i] += atoms.su3[k, m, n] * oStruc[k * params.n + i]
                
                for i in range(params.n):
                    exponent = 1j * (-(self.qC1 * uNM1[i] + self.qC2 * uNM2[i] + self.qC3 * uNM3[i] + q0Gr[i]))
                    exp_term = np.exp(exponent)
                    fExpNM = exp_term * atoms.fAtom[m, n] * oPhase[m * params.n + i] * region[i]
                    
                    idx = n * params.n + i
                    fExpN_r[idx] += fExpNM.real
                    fU1ExpN_r[idx] += fExpNM.real * uNM1[i]
                    fU2ExpN_r[idx] += fExpNM.real * uNM2[i]
                    fU3ExpN_r[idx] += fExpNM.real * uNM3[i]
                    
                    fExpN_i[idx] += fExpNM.imag
                    fU1ExpN_i[idx] += fExpNM.imag * uNM1[i]
                    fU2ExpN_i[idx] += fExpNM.imag * uNM2[i]
                    fU3ExpN_i[idx] += fExpNM.imag * uNM3[i]
        
        # Perform FFTs using numpy
        # Reshape arrays for FFT
        normalization = 1.0 / float(params.n)
        
        # Transfer data to GPU if CuPy is available
        if HAS_CUPY:
            fExpN_r_gpu = cp.asarray(fExpN_r)
            fExpN_i_gpu = cp.asarray(fExpN_i)
            fU1ExpN_r_gpu = cp.asarray(fU1ExpN_r)
            fU1ExpN_i_gpu = cp.asarray(fU1ExpN_i)
            fU2ExpN_r_gpu = cp.asarray(fU2ExpN_r)
            fU2ExpN_i_gpu = cp.asarray(fU2ExpN_i)
            fU3ExpN_r_gpu = cp.asarray(fU3ExpN_r)
            fU3ExpN_i_gpu = cp.asarray(fU3ExpN_i)
        else:
            fExpN_r_gpu = fExpN_r
            fExpN_i_gpu = fExpN_i
            fU1ExpN_r_gpu = fU1ExpN_r
            fU1ExpN_i_gpu = fU1ExpN_i
            fU2ExpN_r_gpu = fU2ExpN_r
            fU2ExpN_i_gpu = fU2ExpN_i
            fU3ExpN_r_gpu = fU3ExpN_r
            fU3ExpN_i_gpu = fU3ExpN_i
        
        # FFT of each atom's contribution (on GPU if available)
        fExpN_rk = cp.zeros((atoms.nAtom, params.Cn1, params.Cn2, params.Cn3), dtype=complex)
        fExpN_ik = cp.zeros((atoms.nAtom, params.Cn1, params.Cn2, params.Cn3), dtype=complex)
        fU1ExpN_rk = cp.zeros((atoms.nAtom, params.Cn1, params.Cn2, params.Cn3), dtype=complex)
        fU1ExpN_ik = cp.zeros((atoms.nAtom, params.Cn1, params.Cn2, params.Cn3), dtype=complex)
        fU2ExpN_rk = cp.zeros((atoms.nAtom, params.Cn1, params.Cn2, params.Cn3), dtype=complex)
        fU2ExpN_ik = cp.zeros((atoms.nAtom, params.Cn1, params.Cn2, params.Cn3), dtype=complex)
        fU3ExpN_rk = cp.zeros((atoms.nAtom, params.Cn1, params.Cn2, params.Cn3), dtype=complex)
        fU3ExpN_ik = cp.zeros((atoms.nAtom, params.Cn1, params.Cn2, params.Cn3), dtype=complex)
        
        for a in range(atoms.nAtom):
            # Reshape to 3D for FFT
            start_idx = a * params.n
            end_idx = (a + 1) * params.n
            
            fExpN_rk[a] = cp_rfftn(fExpN_r_gpu[start_idx:end_idx].reshape(params.nx, params.ny, params.nz)) * normalization
            fExpN_ik[a] = cp_rfftn(fExpN_i_gpu[start_idx:end_idx].reshape(params.nx, params.ny, params.nz)) * normalization
            fU1ExpN_rk[a] = cp_rfftn(fU1ExpN_r_gpu[start_idx:end_idx].reshape(params.nx, params.ny, params.nz)) * normalization
            fU1ExpN_ik[a] = cp_rfftn(fU1ExpN_i_gpu[start_idx:end_idx].reshape(params.nx, params.ny, params.nz)) * normalization
            fU2ExpN_rk[a] = cp_rfftn(fU2ExpN_r_gpu[start_idx:end_idx].reshape(params.nx, params.ny, params.nz)) * normalization
            fU2ExpN_ik[a] = cp_rfftn(fU2ExpN_i_gpu[start_idx:end_idx].reshape(params.nx, params.ny, params.nz)) * normalization
            fU3ExpN_rk[a] = cp_rfftn(fU3ExpN_r_gpu[start_idx:end_idx].reshape(params.nx, params.ny, params.nz)) * normalization
            fU3ExpN_ik[a] = cp_rfftn(fU3ExpN_i_gpu[start_idx:end_idx].reshape(params.nx, params.ny, params.nz)) * normalization
        
        # Calculate intensities (on GPU if available)
        IA = cp.zeros(params.Cn1 * params.Cn2 * params.Cn3)
        IB = cp.zeros(params.Cn1 * params.Cn2 * params.Cn3)
        
        AmpA = cp.zeros(params.Cn1 * params.Cn2 * params.Cn3, dtype=complex)
        AmpB = cp.zeros(params.Cn1 * params.Cn2 * params.Cn3, dtype=complex)
        
        # Transfer lookup arrays to GPU if needed
        if HAS_CUPY:
            mqkA1_gpu = cp.asarray(self.mqkA1)
            mqkA2_gpu = cp.asarray(self.mqkA2)
            mqkA3_gpu = cp.asarray(self.mqkA3)
            mqkB1_gpu = cp.asarray(self.mqkB1)
            mqkB2_gpu = cp.asarray(self.mqkB2)
            mqkB3_gpu = cp.asarray(self.mqkB3)
            fA_gpu = cp.asarray(self.fA)
            fB_gpu = cp.asarray(self.fB)
            kfExpA_gpu = cp.asarray(self.kfExpA)
            kfExpB_gpu = cp.asarray(self.kfExpB)
        else:
            mqkA1_gpu = self.mqkA1
            mqkA2_gpu = self.mqkA2
            mqkA3_gpu = self.mqkA3
            mqkB1_gpu = self.mqkB1
            mqkB2_gpu = self.mqkB2
            mqkB3_gpu = self.mqkB3
            fA_gpu = self.fA
            fB_gpu = self.fB
            kfExpA_gpu = self.kfExpA
            kfExpB_gpu = self.kfExpB
        
        for i in range(params.Cn1):
            for j in range(params.Cn2):
                for k in range(params.Cn3):
                    idx = i * params.Cn2 * params.Cn3 + j * params.Cn3 + k
                    
                    for a in range(atoms.nAtom):
                        aidx = a * params.Cn1 * params.Cn2 * params.Cn3 + idx
                        
                        temp_A = fExpN_rk[a, i, j, k] + 1j * fExpN_ik[a, i, j, k]
                        temp_A -= 1j * (mqkA1_gpu[idx] * fU1ExpN_rk[a, i, j, k] +
                                        mqkA2_gpu[idx] * fU2ExpN_rk[a, i, j, k] +
                                        mqkA3_gpu[idx] * fU3ExpN_rk[a, i, j, k])
                        temp_A += (mqkA1_gpu[idx] * fU1ExpN_ik[a, i, j, k] +
                                   mqkA2_gpu[idx] * fU2ExpN_ik[a, i, j, k] +
                                   mqkA3_gpu[idx] * fU3ExpN_ik[a, i, j, k])
                        
                        temp_B = cp.conj(fExpN_rk[a, i, j, k]) + 1j * cp.conj(fExpN_ik[a, i, j, k])
                        temp_B -= 1j * (mqkB1_gpu[idx] * cp.conj(fU1ExpN_rk[a, i, j, k]) +
                                        mqkB2_gpu[idx] * cp.conj(fU2ExpN_rk[a, i, j, k]) +
                                        mqkB3_gpu[idx] * cp.conj(fU3ExpN_rk[a, i, j, k]))
                        temp_B += (mqkB1_gpu[idx] * cp.conj(fU1ExpN_ik[a, i, j, k]) +
                                   mqkB2_gpu[idx] * cp.conj(fU2ExpN_ik[a, i, j, k]) +
                                   mqkB3_gpu[idx] * cp.conj(fU3ExpN_ik[a, i, j, k]))
                        
                        # Get fA and fB for this atom and position
                        fA_val = fA_gpu[aidx]
                        fB_val = fB_gpu[aidx]
                        
                        AmpA[idx] += temp_A * fA_val
                        AmpB[idx] += temp_B * fB_val
                    
                    AmpA[idx] *= kfExpA_gpu[idx]
                    AmpB[idx] *= kfExpB_gpu[idx]
        
        for i in range(params.Cn1 * params.Cn2 * params.Cn3):
            IA[i] = max(abs(AmpA[i]) ** 2, 1e-50)
            IB[i] = max(abs(AmpB[i]) ** 2, 1e-50)
        
        # Transfer results back to CPU if using GPU
        if HAS_CUPY:
            IA = cp.asnumpy(IA)
            IB = cp.asnumpy(IB)
        
        self.ArrayFourierToRegular(IA, IB, IDiffr, params)
    
    def ArrayFourierToRegular(self, A, B, out, params):
        """
        Transform from Fourier space to regular space arrangement.
        
        Args:
            A: Array A in Fourier space
            B: Array B in Fourier space (conjugate part)
            out: Output array in regular space
            params: SizeContext object
        """
        temp = np.zeros(params.n)
        
        for i in range(params.Cn1):
            for j in range(params.Cn2):
                for k in range(params.Cn3):
                    kk = (params.nz - k) % params.nz
                    jj = (params.ny - j) % params.ny
                    ii = (params.nx - i) % params.nx
                    
                    temp_idx = i * params.Rn2 * params.Rn3 + j * params.Rn3 + k
                    fourier_idx = i * params.Cn2 * params.Cn3 + j * params.Cn3 + k
                    
                    temp[temp_idx] = A[fourier_idx]
                    
                    if kk >= params.Cn3 and kk < params.nz and jj < params.ny and ii < params.nx:
                        b_idx = ii * params.Rn2 * params.Rn3 + jj * params.Rn3 + kk
                        temp[b_idx] = B[fourier_idx]
        
        # Shift to center the output
        for i in range(params.Rn1):
            for j in range(params.Rn2):
                for k in range(params.Rn3):
                    ii = (i + params.nx // 2) % params.nx
                    jj = (j + params.ny // 2) % params.ny
                    kk = (k + params.nz // 2) % params.nz
                    
                    out_idx = ii * params.Rn2 * params.Rn3 + jj * params.Rn3 + kk
                    temp_idx = i * params.Rn2 * params.Rn3 + j * params.Rn3 + k
                    
                    out[out_idx] = temp[temp_idx]
