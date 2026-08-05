"""Probe scale-space cavity candidates on arbitrary triangle soups."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import scipy.ndimage as ndi
import trimesh


@dataclass
class Candidate:
    radius: int
    voxels: int
    volume: float
    bounds_voxels: tuple[int, int, int]
    box_ratio: float
    equivalent_depth: float


def probe(path: Path, resolution: int, maximum_radius: int) -> list[Candidate]:
    mesh = trimesh.load(path, force="mesh", process=False)
    pitch = float(np.max(mesh.extents)) / resolution
    surface = mesh.voxelized(pitch=pitch, method="subdivide")
    matrix = np.pad(surface.matrix.astype(bool), maximum_radius + 3)
    distance = ndi.distance_transform_edt(~matrix)
    structure = ndi.generate_binary_structure(3, 1)

    # Radius-zero exterior is the reference envelope of the triangle soup.
    def envelope(radius: int) -> np.ndarray:
        centers = distance > radius + 0.25
        seeds = np.zeros_like(centers)
        seeds[0, :, :] = centers[0, :, :]
        seeds[-1, :, :] = centers[-1, :, :]
        seeds[:, 0, :] = centers[:, 0, :]
        seeds[:, -1, :] = centers[:, -1, :]
        seeds[:, :, 0] = centers[:, :, 0]
        seeds[:, :, -1] = centers[:, :, -1]
        reachable = ndi.binary_propagation(seeds, structure=structure, mask=centers)
        if radius:
            reachable = ndi.binary_dilation(reachable, structure=ball(radius))
        return ~reachable

    baseline = envelope(0)
    previous = baseline
    result: list[Candidate] = []
    for radius in range(1, maximum_radius + 1):
        current = envelope(radius)
        added = current & ~previous
        labels, count = ndi.label(added, structure=structure)
        objects = ndi.find_objects(labels)
        for label, slices in enumerate(objects, 1):
            if slices is None:
                continue
            component = labels[slices] == label
            voxels = int(component.sum())
            if voxels < 4:
                continue
            dimensions = tuple(int(value.stop - value.start) for value in slices)
            box_volume = int(np.prod(dimensions))
            boundary = component & ~ndi.binary_erosion(component, structure=structure)
            mouth_area = max(int(boundary.sum()), 1) * pitch * pitch
            result.append(Candidate(
                radius=radius,
                voxels=voxels,
                volume=voxels * pitch**3,
                bounds_voxels=dimensions,
                box_ratio=voxels / box_volume,
                equivalent_depth=voxels * pitch**3 / mouth_area,
            ))
        previous = current
    print(f"model={path.stem} pitch={pitch:.6g} grid={matrix.shape} surface={matrix.sum()}")
    for item in sorted(result, key=lambda value: value.volume, reverse=True)[:80]:
        print(item)
    return result


def ball(radius: int) -> np.ndarray:
    coordinates = np.arange(-radius, radius + 1)
    x, y, z = np.meshgrid(coordinates, coordinates, coordinates, indexing="ij")
    return x * x + y * y + z * z <= radius * radius


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--resolution", type=int, default=256)
    parser.add_argument("--maximum-radius", type=int, default=12)
    args = parser.parse_args()
    probe(args.input, args.resolution, args.maximum_radius)


if __name__ == "__main__":
    main()
