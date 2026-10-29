# Examples

Each directory is a complete run in the Q-POP-FerroDyn layout. To run one:

```bash
cd examples/4_polar_vortex_PbTiO3/disordered
qpop-diffraction                                  # reads input.json, system_setting.json, parameter.atom.in
qpop-diffraction-rsm -d diffraction.h5            # 2-D projections
```

| Example | Grid | Fields | What it shows |
|---|---|---|---|
| `1_two_phase_mixture/isotropic`, `anisotropic` | 800×800×1 | `phase.h5` (two phase fractions) | diffuse scattering of a two-phase microstructure |
| `2_two_variant_twin` | 2000×20×1 | `Displacement_0.h5` | peak splitting by two twin variants (rotated `aC`) |
| `3_three_variant_hexagonal_twin` | 1024×1024×1 | `Displacement_0.h5`, `region.h5` (circular grain) | three-variant twinning in a hexagonal lattice |
| `4_polar_vortex_PbTiO3/ordered`, `disordered` | 1×400×400 | `P_Lattice_0.h5` | satellites of a polar vortex array in PbTiO3 |

`examples/input.json` is a template for real FerroDyn runs, and `run_mgpu.sbatch` is a Slurm job script for the multi-GPU program.

## Where they come from

They are the four examples of the original C++ code (`Examples.tar.gz`, text inputs), rebuilt by `convert_legacy_examples.py`. The fields are stored as float32 with shuffle + gzip. The conversion:

- **Units:** `displace.in` (nm) becomes `Displacement_x/y/z` (m), and `strucOrd.in` becomes `P_Lattice_x/y/z` (C/m²). Examples without P get a zero `P_Lattice_0.h5`, and their `parameter.atom.in` gets zero P→shift coefficients for the three P components.
- **Phases and region:** `phaseFra.in` becomes `phase.h5:/phase_<id>`, and `region.in` becomes `region.h5:/region`. The examples that had no `region.in` get the legacy code's default tanh window as `region.h5`, so their patterns stay comparable with the original ones.
- **Lattice:** each phase is an `Is_FE` material whose lattice is `diag(aC)` of `parameter.atom.in`. That makes the average strain zero, as in the legacy inputs. These examples are synthetic, so their material files only serve this bookkeeping.

`4_polar_vortex_PbTiO3/disordered/reference_legacy.h5` holds the intensity of a legacy-format run on the same inputs. `tests/test_examples.py` checks the current code against it, and runs every example.

Rebuild from the tarball with:
```bash
python examples/convert_legacy_examples.py /path/to/Examples.tar.gz --legacy-outputs /path/to/extracted/Examples
```
