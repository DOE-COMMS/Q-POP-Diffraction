"""
JSON input file parser for Q-POP-Diffraction.

Reads input.json which specifies:
  - platt_file        : path to P_Lattice_*.h5 / PLatt_*.h5 (required)
  - displacement_file : path to Displacement_*.h5 (optional; null → zero displacements)
  - x_limits          : [start, end]  0-indexed, inclusive (optional; null → full axis)
  - y_limits          : [start, end]  0-indexed, inclusive (optional; null → full axis)
  - z_limits          : [start, end]  0-indexed, inclusive (optional; null → full axis)
  - phase_file        : HDF5 of phase fractions phase_1 .. phase_<nPhase> (optional; null → one-hot
                        phases from materialField)
  - region_file       : HDF5 with a scattering weight `region` (optional; null → 1)
  - output_file       : base name for output files, no extension (optional; default "diffraction")
"""

from __future__ import annotations

import json
from pathlib import Path


def load_input_json(path: str | Path = "input.json") -> dict:
    """
    Parse and validate input.json.

    Returns a dict with all keys present and defaults applied:
        {
            "platt_file":        Path,
            "displacement_file": Path | None,
            "phase_file":        Path | None,
            "region_file":       Path | None,
            "x_limits":          [int, int] | None,
            "y_limits":          [int, int] | None,
            "z_limits":          [int, int] | None,
            "output_file":       str,
        }

    Raises
    ------
    FileNotFoundError
        If a required or explicitly-specified file does not exist.
    ValueError
        If the JSON is malformed or a parameter value is out of range.
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Input JSON not found: {path}")

    with path.open("r", encoding="utf-8") as fh:
        try:
            data = json.load(fh)
        except json.JSONDecodeError as exc:
            raise ValueError(f"Failed to parse {path}: {exc}") from exc

    # ── required ──────────────────────────────────────────────────────────────
    if "platt_file" not in data:
        raise ValueError("input.json must contain 'platt_file'.")

    platt_path = Path(data["platt_file"])
    if not platt_path.exists():
        raise FileNotFoundError(f"platt_file not found: {platt_path}")

    # ── optional displacement ─────────────────────────────────────────────────
    disp_path: Path | None = None
    if data.get("displacement_file") is not None:
        disp_path = Path(data["displacement_file"])
        if not disp_path.exists():
            raise FileNotFoundError(f"displacement_file not found: {disp_path}")

    # ── optional phase fractions and scattering weight ────────────────────────
    def _optional_file(key: str) -> Path | None:
        if data.get(key) is None:
            return None
        file = Path(data[key])
        if not file.is_file():
            raise FileNotFoundError(f"{key} not found: {file}")
        return file

    if data.get("phaseFra_file") is not None:
        raise ValueError("'phaseFra_file' (text) was replaced by 'phase_file' (HDF5 datasets phase_1 .. phase_<nPhase>).")
    phase_path = _optional_file("phase_file")
    region_path = _optional_file("region_file")

    # ── x/y/z limits ─────────────────────────────────────────────────────────
    def _validate_limits(key: str, raw) -> list[int] | None:
        if raw is None:
            return None
        if not isinstance(raw, (list, tuple)) or len(raw) != 2:
            raise ValueError(f"'{key}' must be a 2-element [start, end] array, got: {raw!r}")
        start, end = int(raw[0]), int(raw[1])
        if start < 0:
            raise ValueError(f"'{key}' start index must be >= 0, got {start}.")
        if end < start:
            raise ValueError(f"'{key}' end ({end}) must be >= start ({start}).")
        return [start, end]

    x_limits = _validate_limits("x_limits", data.get("x_limits"))
    y_limits = _validate_limits("y_limits", data.get("y_limits"))
    z_limits = _validate_limits("z_limits", data.get("z_limits"))

    # ── output file base name ─────────────────────────────────────────────────
    output_file = str(data.get("output_file", "diffraction"))

    return {
        "platt_file": platt_path,
        "displacement_file": disp_path,
        "phase_file": phase_path,
        "region_file": region_path,
        "x_limits": x_limits,
        "y_limits": y_limits,
        "z_limits": z_limits,
        "output_file": output_file,
    }
