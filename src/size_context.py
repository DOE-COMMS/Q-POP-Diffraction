"""
Size context class for handling simulation parameters and grid setup
"""

import numpy as np
import math


class SizeContext:
    """Container for simulation size parameters and grid configuration"""
    
    def __init__(self):
        # System dimensions
        self.lx = 0.0
        self.ly = 0.0
        self.lz = 0.0
        
        # Grid points
        self.nx = 0
        self.ny = 0
        self.nz = 0
        
        # Film parameters
        self.ns = 0  # Substrate thickness
        self.nf = 0  # Film thickness
        
        # Phase and structure parameters
        self.nPhase = 0
        self.nStruc = 0
        
        # Strain average
        self.strainAvg = np.zeros(6)
        
        # Computed parameters
        self.n = 0
        self.dx = 0.0
        self.dy = 0.0
        self.dz = 0.0
        
        self.k1 = 0
        self.k2 = 0
        self.kt = 0.0
        
        self.ns1 = 0
        self.k0 = 0
        
        self.hf = 0.0
        self.hs = 0.0
        self.h1 = 0.0
        self.h2 = 0.0
        self.hd = 0.0
        self.h = 0.0
        self.hstep = 0.0
        
        # FFT dimensions
        self.Rn1 = 0  # Real space
        self.Rn2 = 0
        self.Rn3 = 0
        
        self.Cn1 = 0  # Complex/Fourier space
        self.Cn2 = 0
        self.Cn3 = 0
    
    def read(self, filename):
        """
        Read parameters from input file.
        
        Args:
            filename: Path to parameter file
            
        Returns:
            bool: True if successful, False otherwise
        """
        try:
            with open(filename, 'r') as f:
                lines = f.readlines()
        except IOError:
            print(f"Cannot open file: {filename}")
            return False
        
        numbers = []
        for line in lines:
            # Remove comments (after '!')
            if '!' in line:
                line = line[:line.index('!')]
            
            # Strip whitespace
            line = line.strip()
            if not line:
                continue
            
            # Extract numbers
            tokens = line.split()
            for token in tokens:
                try:
                    numbers.append(float(token))
                except ValueError:
                    continue
        
        if len(numbers) < 16:
            print("Not enough parameters in input file.")
            return False
        
        # Parse parameters
        idx = 0
        self.lx = numbers[idx]; idx += 1
        self.ly = numbers[idx]; idx += 1
        self.lz = numbers[idx]; idx += 1
        self.nx = int(numbers[idx]); idx += 1
        self.ny = int(numbers[idx]); idx += 1
        self.nz = int(numbers[idx]); idx += 1
        self.ns = int(numbers[idx]); idx += 1
        self.nf = int(numbers[idx]); idx += 1
        self.nPhase = int(numbers[idx]); idx += 1
        self.nStruc = int(numbers[idx]); idx += 1
        
        for i in range(6):
            self.strainAvg[i] = numbers[idx]
            idx += 1
        
        print("\nRead parameters:")
        print(f"  lx = {self.lx}, ly = {self.ly}, lz = {self.lz}")
        print(f"  nx = {self.nx}, ny = {self.ny}, nz = {self.nz}")
        print(f"  ns = {self.ns}, nf = {self.nf}")
        print(f"  nPhase = {self.nPhase}, nStruc = {self.nStruc}")
        print(f"  strainAvg = {self.strainAvg}\n")
        
        # Validation
        if self.nPhase < 1 or self.nPhase > 12:
            print(f"This program allows 1~12 phases only. You cannot claim {self.nPhase} phases.")
            return False
        
        if self.nStruc < 0 or self.nStruc > 12:
            print(f"This program allows 0~12 structural order parameters only. You cannot claim {self.nStruc} structural order parameters.")
            return False
        
        return True
    
    def setup(self):
        """
        Setup computed parameters based on read values.
        
        Returns:
            bool: True if successful
        """
        self.n = self.nx * self.ny * self.nz
        
        self.Rn1 = self.nx
        self.Rn2 = self.ny
        self.Rn3 = self.nz
        
        self.Cn1 = self.nx
        self.Cn2 = self.ny
        self.Cn3 = math.floor(self.nz / 2) + 1
        
        self.dx = self.lx / self.nx
        self.dy = self.ly / self.ny
        self.dz = self.lz / self.nz
        
        self.k1 = self.ns + 1
        self.k2 = self.ns + self.nf
        
        if self.nf == 0:
            self.k1 = 1
            self.k2 = self.nz
        
        self.ns1 = 0
        self.k0 = self.ns1 + 1
        
        if self.nf != 0:
            self.k1 = self.ns + 1  # film bottom
            self.k2 = self.ns + self.nf  # film surface
            
            if self.k2 > self.nz:
                self.k2 = self.nz
        else:
            self.hf = self.nf * self.dz
            self.hs = (self.ns - self.ns1 - 1) * self.dz
            self.h1 = self.dz * (self.k1 - 1)
            self.h2 = self.dz * (self.k2 - 1)
            self.hd = self.h2 - self.h1
            
            self.h = max(self.hs, self.hf)
            self.hstep = self.dz
        
        return True
