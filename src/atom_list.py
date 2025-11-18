"""
Atom list class for handling atomic structure and form factors
"""

import numpy as np
from .constants import L0, OS0


class AtomList:
    """Container for atomic structure information and form factors"""
    
    def __init__(self):
        self.nPhase = 0
        self.nStruc = 0
        
        # Reciprocal lattice vector
        self.q00 = np.zeros(3)
        
        # Crystal lattice vectors
        self.aC = np.zeros((3, 3))
        
        # Atomic positions (fractional coordinates)
        self.xAtom = None
        
        # Form factors (complex)
        self.fAtom = None
        self.fAtomR = None
        self.fAtomI = None
        
        # Structural displacements
        self.su1 = None
        self.su2 = None
        self.su3 = None
        
        self.nAtom = 0
    
    def read(self, filename, nPhase, nStruc):
        """
        Read atomic parameters from input file.
        
        Args:
            filename: Path to parameter file
            nPhase: Number of phases
            nStruc: Number of structural order parameters
            
        Returns:
            bool: True if successful, False otherwise
        """
        self.nPhase = nPhase
        self.nStruc = nStruc
        
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
        
        idx = 0
        
        # Read q00
        self.q00[0] = numbers[idx]; idx += 1
        self.q00[1] = numbers[idx]; idx += 1
        self.q00[2] = numbers[idx]; idx += 1
        
        # Read crystal lattice vectors
        for i in range(3):
            for j in range(3):
                self.aC[i, j] = numbers[idx]
                idx += 1
        
        # Read number of atoms
        self.nAtom = int(numbers[idx])
        idx += 1
        
        if self.nAtom <= 0:
            print(f"Invalid number of atoms: {self.nAtom}")
            return False
        
        # Read atomic coordinates
        print("Reading atom coordinates")
        self.xAtom = np.zeros((self.nAtom, 3))
        for i in range(self.nAtom):
            self.xAtom[i, 0] = numbers[idx]; idx += 1
            self.xAtom[i, 1] = numbers[idx]; idx += 1
            self.xAtom[i, 2] = numbers[idx]; idx += 1
        
        # Read atom form factors
        print("Reading atom form factors")
        self.fAtom = np.zeros((nPhase, self.nAtom), dtype=complex)
        self.fAtomR = np.zeros((nPhase, self.nAtom))
        self.fAtomI = np.zeros((nPhase, self.nAtom))
        
        # Initialize structural displacements
        self.su1 = np.zeros((nStruc, nPhase, self.nAtom))
        self.su2 = np.zeros((nStruc, nPhase, self.nAtom))
        self.su3 = np.zeros((nStruc, nPhase, self.nAtom))
        
        for i in range(nPhase):
            m = int(numbers[idx]) - 1  # Convert to 0-indexed
            idx += 1
            print(f"Read input parameters for phase # {m + 1}")
            
            # Read form factors
            for j in range(self.nAtom):
                self.fAtomR[m, j] = numbers[idx]; idx += 1
                self.fAtomI[m, j] = numbers[idx]; idx += 1
                self.fAtom[m, j] = complex(self.fAtomR[m, j], self.fAtomI[m, j])
            
            # Read structural displacements
            for j in range(self.nAtom):
                for k in range(nStruc):
                    self.su1[k, m, j] = numbers[idx]; idx += 1
                    self.su2[k, m, j] = numbers[idx]; idx += 1
                    self.su3[k, m, j] = numbers[idx]; idx += 1
        
        return True
    
    def validate(self):
        """
        Validate and print atom parameters.
        
        Returns:
            bool: True if successful
        """
        print("Validating read atom parameters")
        print(f"q00: {self.q00}")
        
        print("aC:")
        for i in range(3):
            print(f"  {self.aC[i, :]}")
        
        print(f"nAtom: {self.nAtom}")
        
        print("xAtom:")
        for i in range(self.nAtom):
            print(f"  {self.xAtom[i, :]}")
        
        print("fAtom (per phase):")
        for i in range(self.nPhase):
            print(f"  Atom {i}: ", end="")
            for p in range(self.nAtom):
                print(f"({self.fAtomR[i, p]},{self.fAtomI[i, p]}) ", end="")
            print()
        
        print("su1, su2, su3 (per structure, phase, atom):")
        for i in range(self.nPhase):
            print(f"  Phase {i}:")
            for s in range(self.nAtom):
                print(f"    Atom {s}: ")
                for p in range(self.nStruc):
                    print(f"      {self.su1[p, i, s]} {self.su2[p, i, s]} {self.su3[p, i, s]}")
                print()
        
        return True
    
    def setup(self):
        """
        Setup computed parameters.
        
        Returns:
            bool: True if successful
        """
        print("Setting up atom parameters")
        
        # Normalize structural displacements
        for i in range(self.nStruc):
            for j in range(self.nPhase):
                for k in range(self.nAtom):
                    self.su1[i, j, k] /= (L0 / OS0)
                    self.su2[i, j, k] /= (L0 / OS0)
                    self.su3[i, j, k] /= (L0 / OS0)
        
        return True
