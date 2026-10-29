# Method

## Kinematic amplitude

Each cell r of the simulation grid (spacing dx, dy, dz) is treated as a uniform block of crystal. A phase p fills a fraction w_p(r) of it. The cell carries the displacement u(r) and the polarization P(r). Atom a of phase p sits at R_a + u_ap(r), where R_a is its position in the unit cell and

```
u_ap(r) = u(r) + su_pa · P(r)
```

adds the P-induced ionic shift from `parameter.atom.in`. The amplitude around q = qC + k is

```
A(k) = S(k) Σ_a exp(-i (k + G) · R_a) / V · [ FFT(fE_a)(k) - i k · FFT(fE_a u_a)(k) ]

fE_a(r) = Σ_p w_p(r) f_pa exp(-i qC · u_ap(r))
```

and the intensity is I(k) = |A(k)|², floored at 1e-50.

- G is the reciprocal lattice point of `aC` nearest to `q00`. qC is G shifted by the average strain, see below. V is the unit-cell volume.
- FFT is normalised by the number of cells. A perfect single domain therefore gives I(G) = |F(G)|²/V², the squared structure factor per cell volume, whatever the grid size.
- S(k) = Π_d sinc(k_d h_d / 2) is the transform of one cell. Each cell is a continuous block, so a film of t cells has the Laue fringes of a slab of thickness t·dz.
- The term in brackets is exp(-i (qC + k)·u) expanded to first order in k·u. It is exact at the Bragg position k = 0, and accurate when k·u stays small over the region that scatters.

That is four FFTs per atom: fE_a and the three components of fE_a u_a.

## Reciprocal-space grid

The intensity is computed on the FFT grid of the (possibly cropped) simulation box:

- q_d = qC_d + k_d, with k_d = 2π m / (N_d h_d);
- the resolution along axis d is Δq_d = 2π / (N_d h_d), so a larger box gives finer q;
- the range is ±π / h_d around qC.

A 400 nm box with 0.4 nm cells resolves 0.016 nm⁻¹ over ±7.9 nm⁻¹.

## Which cells scatter

The weight of phase p in a cell is its fraction × the region weight. The fraction is one-hot from the material field unless a `phase_file` gives it, and the region weight is 1 unless a `region_file` gives it. Cells whose material is not ferroelectric (substrate, air, vacuum, electrodes) do not scatter. In a thin film the substrate therefore produces no peak, only the film does.

## Strain reference

FerroDyn measures strain from each material's own lattice at the current temperature. The diffraction model uses one lattice, `aC` of `parameter.atom.in`, so every strain is referred to `aC`.

In-plane, the average strain `strainAvg` moves the scattering centre to qC = G − ε·G. It follows FerroDyn's `elastic_solver_base.cpp`:

- a thin film on a substrate takes the substrate's in-plane lattice a(T), b(T);
- a thin film with `MisfitStrain.Enabled` takes its own lattice × (1 + misfit);
- otherwise the film keeps its own lattice;
- `AppliedStrain` on the film is added.

Out-of-plane, each material m adds a ramp to u_z, summed up the full stack from the bottom of the box to the cell's centre:

```
u_z(z) += Σ_{cells z' ≤ z} dz (c_m(z')(T) / aC_z - 1) - dz (c_m(z)(T) / aC_z - 1) / 2
```

Each layer then diffracts at its own c_m(T)·(1 + ∂u_z/∂z). This is what places the 00l peaks of a superlattice such as PbTiO3/SrTiO3 correctly with a single `aC`. The ramp is built over the full z extent before any crop.

Choose `aC` close to the film's average lattice. The ramp then stays small, where the first-order expansion in k·u is most accurate.

## Temperature

The temperature of a snapshot comes from its `Timestep` attribute and FerroDyn's ramp in `SimulationControls.Temperature`:

- with `Transient` true, the temperature is `Initial + Step · floor(step / Interval)`, capped at `Final`, as in FerroDyn's `timer.hpp`; otherwise it is `Initial`;
- it sets each lattice through `a(T) = a · (1 + alpha · (T − ReferenceTemperature))`.

A run that was restarted with a different `Initial` gets the wrong temperature.

## Approximations and limits

- Kinematic: single scattering, no absorption, extinction or instrumental resolution.
- First order in k·u. A homogeneous strain e entered through u puts the peak at G(1 − e) instead of G/(1 + e), an error of G·e². Referring the lattice to `aC` keeps e small.
- Periodic box: the FFT treats the box as periodic. Keep the substrate and air as padding instead of cropping z to the film.
- The order parameter is the lattice polarization P (three components). Other order parameters (e.g. octahedral tilts) are not included.
- Not supported: rotated runs (`Rotation.Enabled`, rejected), and the `Bulk-stress-free` macrostrain, which FerroDyn computes at run time and does not store in the input.
