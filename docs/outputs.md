# Outputs

## The result file

`{output_file}.h5`, written by both the standard build and `qpop_diffraction_mgpu`:

- `/intensity`: float64, shape `(nz, ny, nx)` like the FerroDyn fields, fftshifted so the Bragg position is in the middle of the box;
- `/q_x`, `/q_y`, `/q_z`: the 1-D q axes in nm⁻¹, with q = qC + k;
- attribute `Timestep`: the FerroDyn step of the snapshot;
- attribute `q_units`: `nm^-1`.

The q grid is rectilinear, so `/intensity[k, j, i]` sits at `(q_x[i], q_y[j], q_z[k])`.

`{output_file}.xdmf` describes the same data for ParaView as a 3DRectMesh on the q axes. Open it in ParaView and apply `log10(Intensity)` with the Calculator filter to see the weak features.

## Reading it in Python

```python
from qpop_diffraction import load_diffraction_hdf5

out = load_diffraction_hdf5("diffraction.h5")
I = out["intensity"][0]          # (nx, ny, nz)
qx, qy, qz = out["q_x"], out["q_y"], out["q_z"]
```

The dictionary also holds `log10_intensity` and `q_vector`, the q axes broadcast to `(3, nx, ny, nz)`. Note the `(x, y, z)` order, the reverse of the file's.

## Reciprocal-space maps

```bash
qpop-diffraction-rsm -d diffraction.h5 [-o reciprocal_space_maps] [--bins N] [--aggregation sum|mean]
```

It projects the 3-D intensity onto the qx–qy, qx–qz and qy–qz planes. For each plane it writes:

- `<plane>.npz`, with the map and its log10 on the q axes;
- `<plane>.png`, with log10 intensity on a colour scale.

`--aggregation sum` (the default) integrates over the third axis, like a detector that collects all of it. `mean` averages instead. `--bins` rebins the axes.

## Multi-GPU parameters

With `--mgpu` there is also `{output_file}_params.h5`, the input of `qpop_diffraction_mgpu`:

- datasets `grid`, `crop`, `spacing`, `G`, `qC`, `V`, `R`, `fAtom_re`, `fAtom_im`, `su`, `phase_of`, `dc` and `u_scale`;
- attributes holding the absolute paths of the field files and of the output.

It can be deleted after `qpop_diffraction_mgpu` has run.
