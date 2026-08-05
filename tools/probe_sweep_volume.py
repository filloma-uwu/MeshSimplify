"""Evaluate axis-sweep envelopes by added visual-hull volume."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from scipy import ndimage

from probe_projection_topology import rasterize, read_obj


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--size", type=int, default=128)
    parser.add_argument("--closing-ratio", type=float, default=0.03)
    args = parser.parse_args()
    vertices, faces = read_obj(args.input)
    xy = rasterize(vertices, faces, (0, 1), args.size)
    xz = rasterize(vertices, faces, (0, 2), args.size)
    yz = rasterize(vertices, faces, (1, 2), args.size)
    visual_hull = (
        xy[None, :, :]
        & xz[:, None, :]
        & yz[:, :, None]
    )  # z, y, x
    radius = max(1, round(args.closing_ratio * args.size))
    structure_2d = np.ones((2 * radius + 1, 2 * radius + 1), dtype=bool)
    structure_3d = ndimage.generate_binary_structure(3, 1)
    projections = {
        "x": (yz, lambda mask: np.broadcast_to(mask[:, :, None], visual_hull.shape)),
        "y": (xz, lambda mask: np.broadcast_to(mask[:, None, :], visual_hull.shape)),
        "z": (xy, lambda mask: np.broadcast_to(mask[None, :, :], visual_hull.shape)),
    }
    for axis, (projection, expand) in projections.items():
        closed = ndimage.binary_fill_holes(projection)
        closed |= ndimage.binary_closing(closed, structure=structure_2d)
        candidate = expand(closed)
        added = candidate & ~visual_hull
        labels, count = ndimage.label(added, structure=structure_3d)
        sizes = np.bincount(labels.ravel())[1:]
        ratios = sorted((float(value / visual_hull.size) for value in sizes), reverse=True)
        print(axis, "candidate_ratio", float(candidate.mean()),
              "added_ratio", float(added.mean()), "components", count,
              "largest", ratios[:12])


if __name__ == "__main__":
    main()
