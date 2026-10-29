"""
Size context class for handling simulation parameters and grid setup
"""

import json
import math
from pathlib import Path

import numpy as np

from .constants import L0


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
    
    def read_json(self, filename):
        """
        Read the grid and phases from a Q-POP-FerroDyn system_setting.json
        and the materials/<name>.json files beside it.

        Phases are the Is_FE materials, in MaterialNames order, so
        parameter.atom.in needs one phase block per FE material.  The average
        strain is set later by reference_strain(), once aC and the snapshot
        temperature are known.

        Returns:
            list: material indices of the phases, in phase order
        """
        path = Path(filename)
        self.setting = json.loads(path.read_text())
        if self.setting.get("Rotation", {}).get("Enabled", False):
            raise ValueError("Rotated runs (Rotation.Enabled) are not supported.")

        mesh = self.setting["Mesh"]
        self.nx, self.ny, self.nz = mesh["NX"], mesh["NY"], mesh["NZ"]
        self.lx = self.nx * mesh["DX"] / L0
        self.ly = self.ny * mesh["DY"] / L0
        self.lz = self.nz * mesh["DZ"] / L0
        self.ns = 0
        self.nf = 0
        self.nStruc = 3

        self.names = self.setting["Materials"]["MaterialNames"]
        self.materials = [json.loads((path.parent / "materials" / f"{n}.json").read_text()) for n in self.names]
        fe_ids = [i for i, m in enumerate(self.materials) if m.get("Ferroelectric", {}).get("Is_FE", False)]
        self.nPhase = len(fe_ids)
        self.film = fe_ids[0]

        print("\nRead parameters:")
        print(f"  lx = {self.lx}, ly = {self.ly}, lz = {self.lz}")
        print(f"  nx = {self.nx}, ny = {self.ny}, nz = {self.nz}")
        print(f"  nPhase = {self.nPhase} ({[self.names[i] for i in fe_ids]}), nStruc = {self.nStruc}\n")

        return fe_ids

    def temperature(self, step):
        """Temperature at a timestep, following FerroDyn's ramp (timer.hpp updateTemperature)."""
        t = self.setting["SimulationControls"]["Temperature"]
        if not t.get("Transient", False):
            return t["Initial"]
        T = t["Initial"] + t["Step"] * (step // t["Interval"])
        return min(T, t["Final"]) if t["Step"] > 0 else max(T, t["Final"])

    def lattice(self, m, T):
        """Lattice constants (a, b, c) of material m at temperature T in nm, as in FerroDyn materials.cpp."""
        mech = self.materials[m]["Mechanical"]
        lc = mech["LatticeConstants"]
        alpha = mech["ThermalExpansion"]
        return np.array([lc[k] * (1.0 + alpha[i][i] * (T - lc["ReferenceTemperature"])) for i, k in enumerate("abc")]) / L0

    def reference_strain(self, aC, T):
        """
        Refer FerroDyn's strains, which are measured from each material's own
        lattice, to the diagonal lattice aC (nm) of parameter.atom.in.

        Sets strainAvg to the in-plane strain of the film relative to aC:
        a thin film on a substrate takes the substrate's in-plane lattice, a
        film with MisfitStrain takes its own lattice times (1 + misfit), as in
        elastic_solver_base.cpp; otherwise the film keeps its own lattice.
        AppliedStrain on the film is added.

        Returns:
            np.ndarray: c_m(T) / aC_z - 1 per material, the out-of-plane
            reference mismatch to be carried by the displacement field
        """
        elastic = self.setting["Solvers"].get("ElasticStrain", {}).get("InitialConditions", {})
        misfit = elastic.get("MisfitStrain", {})
        applied = elastic.get("AppliedStrain", {})
        layers = self.setting["Materials"]["LayerTypes"]
        air = [i for i, n in enumerate(self.names) if n in ("Air", "Vacuum")]
        sub = layers[0] if layers[0] not in air and any(t not in [layers[0]] + air for t in layers) else None

        a = self.lattice(self.film, T)
        self.strainAvg = np.zeros(6)
        if self.setting["BoundaryConditions"].get("GeometryType", "Bulk").lower() == "thin-film":
            if misfit.get("Enabled", False):
                a[:2] *= 1.0 + np.array(misfit["Strain"][:2])
                self.strainAvg[5] = misfit["Strain"][2]
            elif sub is not None:
                a[:2] = self.lattice(sub, T)[:2]
        self.strainAvg[:2] = a[:2] / aC[:2] - 1.0
        if applied.get("Enabled", False) and str(self.film) in applied:
            self.strainAvg += applied[str(self.film)]

        print(f"  T = {T} K, strainAvg = {self.strainAvg}")
        return np.array([self.lattice(m, T)[2] for m in range(len(self.names))]) / aC[2] - 1.0

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

    def apply_slice(self, x_limits, y_limits, z_limits):
        """
        Restrict the simulation domain to a sub-region defined by index limits.

        Must be called after :meth:`read` and before array allocation.  The
        method saves the original grid spacing, updates nx/ny/nz and lx/ly/lz
        to match the slice, and re-runs :meth:`setup` to refresh all derived
        attributes (n, FFT dims, k1/k2, …).

        ns and nf are reset to 0 so that the slice is treated as a uniform
        bulk domain (no film/substrate phase-forcing is applied to the sliced
        data).

        Parameters
        ----------
        x_limits : [int, int] or None
            0-indexed inclusive [start, end] along x; None → full axis.
        y_limits : [int, int] or None
            0-indexed inclusive [start, end] along y; None → full axis.
        z_limits : [int, int] or None
            0-indexed inclusive [start, end] along z; None → full axis.

        Returns
        -------
        tuple (x0, x1, y0, y1, z0, z1)
            The resolved integer bounds actually applied.

        Raises
        ------
        IndexError
            If any limit exceeds the corresponding full-domain axis size.
        """
        # Compute original grid spacing (must be called after read())
        dx_orig = self.lx / self.nx
        dy_orig = self.ly / self.ny
        dz_orig = self.lz / self.nz

        # Resolve defaults
        x0 = x_limits[0] if x_limits is not None else 0
        x1 = x_limits[1] if x_limits is not None else self.nx - 1
        y0 = y_limits[0] if y_limits is not None else 0
        y1 = y_limits[1] if y_limits is not None else self.ny - 1
        z0 = z_limits[0] if z_limits is not None else 0
        z1 = z_limits[1] if z_limits is not None else self.nz - 1

        # Bounds check
        if x1 >= self.nx:
            raise IndexError(f"x_limits end {x1} is out of range for nx={self.nx}.")
        if y1 >= self.ny:
            raise IndexError(f"y_limits end {y1} is out of range for ny={self.ny}.")
        if z1 >= self.nz:
            raise IndexError(f"z_limits end {z1} is out of range for nz={self.nz}.")

        # Update grid dimensions to slice size
        self.nx = x1 - x0 + 1
        self.ny = y1 - y0 + 1
        self.nz = z1 - z0 + 1

        # Update physical sizes proportionally
        self.lx = self.nx * dx_orig
        self.ly = self.ny * dy_orig
        self.lz = self.nz * dz_orig

        # Treat the slice as a uniform bulk domain (no film/substrate logic)
        self.ns = 0
        self.nf = 0

        # Recompute all derived attributes
        self.setup()

        print(f"  Slice applied  : x=[{x0},{x1}], y=[{y0},{y1}], z=[{z0},{z1}]")
        print(f"  Slice grid     : nx={self.nx}, ny={self.ny}, nz={self.nz}")
        print(f"  Slice domain   : lx={self.lx:.4f} nm, ly={self.ly:.4f} nm, lz={self.lz:.4f} nm")

        return x0, x1, y0, y1, z0, z1
