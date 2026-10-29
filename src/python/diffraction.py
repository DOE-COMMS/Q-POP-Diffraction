"""
Kinematic diffraction, evaluated on the GPU with CuPy when it is installed.

For each atom a of the unit cell, with the atom displacement
u_a(r) = u(r) + su_a · P(r) in nm:

    fExpN_a(r) = sum_phases f_a · phase(r) · region(r) · exp(-i qC · u_a(r))
    Amp(k)     = S(k) · sum_a exp(-i (k + G) · R_a) / V · [FFT(fExpN_a) - i k · FFT(fExpN_a u_a)]

k is the FFT wavevector, G the reciprocal lattice point nearest q00, qC the
same point corrected by the average strain, R_a the atom position in the cell,
V the cell volume and S(k) = prod_d sinc(k_d d_d / 2) the voxel shape factor.
The intensity |Amp|^2 is centred on qC: q = qC + k.

Set QPOP_DIFFRACTION_DEVICE=cpu to compute with NumPy even when CuPy is installed,
e.g. on a CPU node or for a grid larger than GPU memory.
"""

import math
import os

import numpy as np

if os.environ.get("QPOP_DIFFRACTION_DEVICE", "").lower() == "cpu":
    xp = np
    print("QPOP_DIFFRACTION_DEVICE=cpu - computing on the CPU (NumPy)")
else:
    try:
        import cupy as xp
        print("CuPy detected - computing on the GPU")
    except ImportError:
        xp = np
        print("CuPy not available - computing on the CPU (NumPy)")


def reciprocal_centre(params, atoms):
    """
    Returns:
        (G, qC): the reciprocal lattice point nearest q00 and the same point
        corrected by the average strain, both in nm^-1
    """
    aC = atoms.aC
    G = np.round(aC @ atoms.q00 / (2 * math.pi)) @ (2 * math.pi * np.linalg.inv(aC).T)
    sA = params.strainAvg
    eps = np.array([[sA[0], sA[5], sA[4]], [sA[5], sA[1], sA[3]], [sA[4], sA[3], sA[2]]])
    qC = G - eps @ G
    print(f"  G = {G}, qC = {qC}")
    return G, qC


def wavevectors(params):
    """FFT wavevectors along x, y, z in nm^-1, in FFT order."""
    shape = (params.nx, params.ny, params.nz)
    return [2 * math.pi * np.fft.fftfreq(m, d) for m, d in zip(shape, (params.dx, params.dy, params.dz))]


def q_axes(params, qC):
    """1-D q axes in nm^-1 of the fftshifted intensity."""
    return tuple(np.fft.fftshift(kd) + qC[d] for d, kd in enumerate(wavevectors(params)))


def diffraction(params, atoms, oPhase, oStruc, u, region):
    """
    Args:
        params: SizeContext after setup() (grid in nm, strainAvg)
        atoms:  AtomList after setup() (aC, q00, xAtom, fAtom; su in nm per C/m^2)
        oPhase: (nPhase * n,) phase fractions
        oStruc: (nStruc * n,) order parameters (P in C/m^2)
        u:      (3 * n,) displacement in nm
        region: (n,) scattering weight

    Returns:
        (I, qx, qy, qz): intensity (nx, ny, nz), fftshifted, and its 1-D q axes in nm^-1
    """
    shape = (params.nx, params.ny, params.nz)
    n = params.n

    aC = atoms.aC
    V = np.linalg.det(aC)
    G, qC = reciprocal_centre(params, atoms)

    # FFT wavevectors along each axis and the voxel shape factor
    k = wavevectors(params)
    kg = [xp.asarray(kd).reshape([-1 if d == i else 1 for i in range(3)]) for d, kd in enumerate(k)]
    S = kg[0] * 0 + 1.0
    for d, h in enumerate((params.dx, params.dy, params.dz)):
        S = S * xp.sinc(kg[d] * h / (2 * math.pi))

    P = xp.asarray(oStruc.reshape(params.nStruc, n))
    U = xp.asarray(u.reshape(3, n))
    w = xp.asarray(oPhase.reshape(params.nPhase, n)) * xp.asarray(region)
    qCg = xp.asarray(qC)
    su = np.stack([atoms.su1, atoms.su2, atoms.su3])     # (direction, nStruc, nPhase, nAtom)
    R = atoms.xAtom @ aC                                  # (nAtom, 3) atom positions, nm

    amp = xp.zeros(shape, dtype=complex)
    for a in range(atoms.nAtom):
        fE = xp.zeros(n, dtype=complex)
        fU = xp.zeros((3, n), dtype=complex)
        for p in range(params.nPhase):
            ua = U + xp.asarray(su[:, :, p, a]) @ P
            e = atoms.fAtom[p, a] * w[p] * xp.exp(-1j * (qCg @ ua))
            fE += e
            fU += e * ua
        FE = xp.fft.fftn(fE.reshape(shape)) / n
        FU = xp.fft.fftn(fU.reshape(3, *shape), axes=(1, 2, 3)) / n
        phase = xp.exp(-1j * (kg[0] + G[0]) * R[a, 0]) * xp.exp(-1j * (kg[1] + G[1]) * R[a, 1]) \
            * xp.exp(-1j * (kg[2] + G[2]) * R[a, 2]) / V
        amp += phase * (FE - 1j * (kg[0] * FU[0] + kg[1] * FU[1] + kg[2] * FU[2]))

    I = xp.fft.fftshift(xp.maximum(xp.abs(S * amp) ** 2, 1e-50))
    I = I.get() if xp is not np else I
    return (I,) + q_axes(params, qC)
