#!/usr/bin/env python3
"""
Reciprocal-space maps from a Q-POP-Diffraction output file.

Projects the 3-D intensity of diffraction.h5 onto the qx-qy, qx-qz and qy-qz
planes and writes each map as .npz (intensity and log10 intensity on its q
axes) and .png.

    qpop-diffraction-rsm -d diffraction.h5 [-o reciprocal_space_maps]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .io_utils import load_diffraction_hdf5


PAIR_DEFINITIONS = (
    (0, 1, "qx_qy"),
    (0, 2, "qx_qz"),
    (1, 2, "qy_qz"),
)
AXIS_LABELS = ("Qx (nm$^{-1}$)", "Qy (nm$^{-1}$)", "Qz (nm$^{-1}$)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--diffraction", "-d",
        type=Path,
        required=True,
        metavar="HDF5",
        help="Diffraction HDF5 written by qpop-diffraction or qpop_diffraction_mgpu (/intensity, /q_x, /q_y, /q_z).",
    )
    parser.add_argument(
        "--output-dir", "-o",
        type=Path,
        default=Path("reciprocal_space_maps"),
        help="Directory for generated maps and figures.",
    )
    parser.add_argument(
        "--bins",
        type=int,
        default=None,
        help=(
            "Override the number of bins per axis used in the 2D histograms. "
            "By default the script uses the number of unique q samples per component."
        ),
    )
    parser.add_argument(
        "--aggregation",
        choices=("sum", "mean"),
        default="sum",
        help="How to accumulate intensity when multiple voxels fall into the same reciprocal-space bin.",
    )
    return parser.parse_args()


def build_edges(values: np.ndarray, n_bins_override: int | None) -> tuple[np.ndarray, np.ndarray]:
    if n_bins_override is not None:
        edges = np.linspace(values.min(), values.max(), n_bins_override + 1)
        centers = 0.5 * (edges[:-1] + edges[1:])
        return edges, centers

    unique_values = np.unique(np.round(values, decimals=10))
    if unique_values.size == 1:
        delta = 0.5
        edges = np.array([unique_values[0] - delta, unique_values[0] + delta])
        centers = unique_values.copy()
        return edges, centers

    centers = unique_values
    midpoints = 0.5 * (centers[:-1] + centers[1:])
    first_edge = centers[0] - 0.5 * (centers[1] - centers[0])
    last_edge = centers[-1] + 0.5 * (centers[-1] - centers[-2])
    edges = np.concatenate(([first_edge], midpoints, [last_edge]))
    return edges, centers


def build_map(
    q_a: np.ndarray,
    q_b: np.ndarray,
    intensity: np.ndarray,
    bins_override: int | None,
    aggregation: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    edges_a, centers_a = build_edges(q_a, bins_override)
    edges_b, centers_b = build_edges(q_b, bins_override)

    weighted_sum, _, _ = np.histogram2d(q_a, q_b, bins=(edges_a, edges_b), weights=intensity)
    if aggregation == "sum":
        return centers_a, centers_b, weighted_sum

    counts, _, _ = np.histogram2d(q_a, q_b, bins=(edges_a, edges_b))
    mean_map = np.divide(weighted_sum, counts, out=np.zeros_like(weighted_sum), where=counts > 0)
    return centers_a, centers_b, mean_map


def save_map_data(output_path: Path, q_a: np.ndarray, q_b: np.ndarray, map_data: np.ndarray) -> None:
    log_map = np.log10(np.maximum(map_data, np.finfo(np.float64).tiny))
    np.savez(
        output_path,
        q_a=q_a,
        q_b=q_b,
        q_units="nm^-1",
        intensity=map_data,
        log10_intensity=log_map,
    )


def save_map_figure(
    output_path: Path,
    q_a: np.ndarray,
    q_b: np.ndarray,
    map_data: np.ndarray,
    label_a: str,
    label_b: str,
    aggregation: str,
) -> None:
    display_data = np.log10(np.maximum(map_data, np.finfo(np.float64).tiny))

    fig, ax = plt.subplots(figsize=(8, 6))
    image = ax.pcolormesh(q_a, q_b, display_data.T, shading="auto", cmap="turbo")
    ax.set_xlabel(label_a)
    ax.set_ylabel(label_b)
    ax.set_title(f"Reciprocal-space map: {label_a} vs {label_b} ({aggregation})")
    colorbar = fig.colorbar(image, ax=ax)
    colorbar.set_label("log10 intensity")
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def main() -> int:
    args = parse_args()

    data = load_diffraction_hdf5(args.diffraction)
    print(f"Loaded diffraction HDF5: {args.diffraction}")

    args.output_dir.mkdir(parents=True, exist_ok=True)

    intensity = data["intensity"][0].reshape(-1)
    q_components = [data["q_vector"][index].reshape(-1) for index in range(3)]

    for axis_a, axis_b, stem in PAIR_DEFINITIONS:
        q_a = q_components[axis_a]
        q_b = q_components[axis_b]
        centers_a, centers_b, map_data = build_map(
            q_a, q_b, intensity, args.bins, args.aggregation
        )

        save_map_data(args.output_dir / f"{stem}.npz", centers_a, centers_b, map_data)
        save_map_figure(
            args.output_dir / f"{stem}.png",
            centers_a,
            centers_b,
            map_data,
            AXIS_LABELS[axis_a],
            AXIS_LABELS[axis_b],
            args.aggregation,
        )
        print(
            f"Wrote {stem}.npz and {stem}.png with shape {map_data.shape} "
            f"to {args.output_dir}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())