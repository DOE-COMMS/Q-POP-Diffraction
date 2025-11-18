#!/usr/bin/env python3
"""
Q-POP Diffraction Simulation
Python conversion of the C++/CUDA implementation
"""

import numpy as np
import sys
import math

from src.constants import OS0
from src.size_context import SizeContext
from src.atom_list import AtomList
from src.diffraction import Diffraction
from src.io_utils import read4D, write4D


def main():
    """Main simulation function"""
    
    # Initialize contexts
    params = SizeContext()
    atoms = AtomList()
    
    # Read and setup parameters
    if not params.read("parameter.system.in"):
        return 1
    if not params.setup():
        return 1
    
    # Read phase fractions from phaseFra.in
    oPhase_size = params.nPhase * params.nx * params.ny * params.nz
    oPhase = np.zeros(oPhase_size)
    
    readErrors = read4D("phaseFra.in", oPhase, params.nPhase, params.nx, params.ny, params.nz)
    
    if readErrors == -1:
        print("File phaseFra.in not provided. Using a pure phase 1.")
        for i in range(params.nx):
            for j in range(params.ny):
                for k in range(params.nz):
                    oPhase[i * params.ny * params.nz + j * params.nz + k] = 1.0
    elif readErrors == 0:
        print("Error reading phaseFra.in. Exiting.", file=sys.stderr)
        return 1
    elif readErrors == 1:
        print("Successfully read phaseFra.in.")
        
        for p in range(params.nPhase):
            for i in range(params.nx):
                for j in range(params.ny):
                    for k in range(params.nz):
                        if k < params.k1 - 2 or k > params.k2:
                            oPhase[p * params.n + i * params.ny * params.nz + j * params.nz + k] = 1.0
    
    # Read structural order parameters from strucOrd.in
    oStruc_size = params.nStruc * params.nx * params.ny * params.nz
    oStruc = np.zeros(oStruc_size)
    
    readErrors = read4D("strucOrd.in", oStruc, params.nStruc, params.nx, params.ny, params.nz)
    
    if readErrors == -1:
        if params.nStruc > 0:
            print("File strucOrd.in not provided. Skipping.")
            return 1
    elif readErrors == 0:
        print("Error reading strucOrd.in. Exiting.", file=sys.stderr)
        return 1
    elif readErrors == 1 and params.nStruc > 0:
        print("Successfully read strucOrd.in.")
        oStruc /= OS0
    
    # Read displacement field from displace.in
    u_size = 3 * params.nx * params.ny * params.nz
    u = np.zeros(u_size)
    
    readErrors = read4D("displace.in", u, 3, params.nx, params.ny, params.nz)
    
    if readErrors == -1:
        print("File displace.in not provided. Using a zero displacement field.")
        u.fill(0.0)
    elif readErrors == 0:
        print("Error reading displace.in. Exiting.", file=sys.stderr)
        return 1
    elif readErrors == 1:
        print("Successfully read displace.in.")
    
    # Read atom parameters
    if not atoms.read("parameter.atom.in", params.nPhase, params.nStruc):
        return 1
    if not atoms.setup():
        return 1
    
    # Read region from region.in
    region_size = params.nx * params.ny * params.nz
    region = np.zeros(region_size)
    
    readErrors = read4D("region.in", region, 1, params.nx, params.ny, params.nz)
    
    if readErrors == -1:
        print("File region.in not provided. Using a default region.")
        
        ix1 = (params.nx + 1) / 2.0 - params.nx / math.pi
        ix2 = (params.nx + 1) / 2.0 + params.nx / math.pi
        
        iy1 = (params.ny + 1) / 2.0 - params.ny / math.pi
        iy2 = (params.ny + 1) / 2.0 + params.ny / math.pi
        
        iz1 = (params.nz + 1) / 2.0 - params.nz / math.pi
        iz2 = (params.nz + 1) / 2.0 + params.nz / math.pi
        
        tx = params.nx / 8.0
        ty = params.ny / 8.0
        tz = params.nz / 8.0
        
        region.fill(1.0)
        
        for i in range(params.nx):
            for j in range(params.ny):
                for k in range(params.nz):
                    idx = i * params.ny * params.nz + j * params.nz + k
                    if params.nf == 0:
                        region[idx] = math.tanh((k + 1 - iz1) / tz) - math.tanh((k + 1 - iz2) / tz)
                    region[idx] *= math.tanh((j + 1 - iy1) / ty) - math.tanh((j + 1 - iy2) / ty)
                    region[idx] *= math.tanh((i + 1 - ix1) / tx) - math.tanh((i + 1 - ix2) / tx)
        
        region_max = np.max(region)
        if region_max == 0.0:
            print("Region max is zero. Exiting.", file=sys.stderr)
            return 1
        region /= region_max
        
    elif readErrors == 0:
        print("Error reading region.in. Exiting.", file=sys.stderr)
        return 1
    elif readErrors == 1:
        print("Successfully read region.in.")
    
    write4D("region.00000000.dat", region, 1, params.nx, params.ny, params.nz)
    
    # Initialize diffraction arrays
    IDiffr = np.full(params.n, 3.0)
    DQ = np.zeros(3 * params.n)
    QCenter = np.zeros(3)
    
    print("\nSetting up diffraction\n")
    
    # Setup and calculate diffraction
    diffContext = Diffraction()
    diffContext.diffraction_setup(params, atoms)
    diffContext.diffraction_calc(params, atoms, IDiffr, region, QCenter, oPhase, oStruc, u)
    
    # Write outputs
    write4D("I.00000000.dat", IDiffr, 1, params.nx, params.ny, params.nz)
    
    # Write log10 of intensity
    IDiffr_log = np.log10(IDiffr)
    write4D("lg_{10}I.00000000.dat", IDiffr_log, 1, params.nx, params.ny, params.nz)
    
    # Write q-vectors
    for idx in range(params.n):
        DQ[0 * params.n + idx] = diffContext.mqk1_out[idx] + diffContext.qC1 - diffContext.qC10
        DQ[1 * params.n + idx] = diffContext.mqk2_out[idx] + diffContext.qC2 - diffContext.qC20
        DQ[2 * params.n + idx] = diffContext.mqk3_out[idx] + diffContext.qC3 - diffContext.qC30
    
    write4D("qVector.00000000.dat", DQ, 3, params.nx, params.ny, params.nz)
    
    print("\nDiffraction simulation completed.")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
