# Notes for AI agents

To compute, visualise or interpret diffraction from Q-POP-FerroDyn output with this repository, read [skills/qpop-diffraction/SKILL.md](skills/qpop-diffraction/SKILL.md) first. It covers setting up `input.json` and `parameter.atom.in`, choosing the standard or multi-GPU build, reading the output, interpreting peaks, and the traps that fail silently.

Repository map:

- `src/python/`: the `qpop_diffraction` Python package (`qpop-diffraction`, `qpop-diffraction-rsm`)
- `src/heffte/`: C++/CUDA source of the multi-GPU program, built into `bin/qpop_diffraction_mgpu` by the root `Makefile`
- `tests/`: pytest suite with analytic answers; run `pytest tests/` after any change to the calculation
- `docs/`: user documentation; `docs/method.md` has the model and its approximations
- `examples/`: ready-to-run cases
