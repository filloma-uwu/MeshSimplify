"""Probe projected holes and boundary concavities before the C++ implementation."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw
from scipy import ndimage
from scipy.spatial import ConvexHull


def read_obj(path: Path) -> tuple[np.ndarray, np.ndarray]:
    vertices: list[list[float]] = []
    faces: list[list[int]] = []
    for raw in path.read_text(encoding="utf-8-sig").splitlines():
        fields = raw.split()
        if fields[:1] == ["v"]:
            vertices.append([float(value) for value in fields[1:4]])
        elif fields[:1] == ["f"]:
            faces.append([int(value.split("/")[0]) - 1 for value in fields[1:4]])
    return np.asarray(vertices), np.asarray(faces)


def rasterize(vertices: np.ndarray, faces: np.ndarray, axes: tuple[int, int], size: int) -> np.ndarray:
    points = vertices[:, axes]
    lower = points.min(axis=0)
    extent = np.maximum(points.max(axis=0) - lower, 1.0e-12)
    margin = 4
    scale = (size - 2 * margin - 1) / max(extent)
    projected = (points - lower) * scale + margin
    image = Image.new("1", (size, size), 0)
    draw = ImageDraw.Draw(image)
    for face in faces:
        polygon = [(float(projected[index, 0]), float(projected[index, 1])) for index in face]
        draw.polygon(polygon, fill=1)
    return np.asarray(image, dtype=bool)


def classify_background(mask: np.ndarray) -> tuple[np.ndarray, np.ndarray, list[dict[str, float]]]:
    background = ~mask
    labels, count = ndimage.label(background)
    border_labels = set(np.unique(np.concatenate([
        labels[0], labels[-1], labels[:, 0], labels[:, -1]
    ])))
    holes = np.zeros_like(mask)
    records: list[dict[str, float]] = []
    exterior = np.isin(labels, list(border_labels))
    for label in range(1, count + 1):
        component = labels == label
        if label in border_labels:
            continue
        holes |= component
        distance = ndimage.distance_transform_edt(component)
        records.append({"pixels": float(component.sum()), "radius": float(distance.max())})
    return holes, exterior, sorted(records, key=lambda item: item["pixels"], reverse=True)


def classify_concavities(mask: np.ndarray) -> list[dict[str, float]]:
    occupied = np.argwhere(mask)
    points = occupied[:, [1, 0]]
    hull_points = points[ConvexHull(points).vertices]
    hull_image = Image.new("1", (mask.shape[1], mask.shape[0]), 0)
    ImageDraw.Draw(hull_image).polygon([tuple(point) for point in hull_points], fill=1)
    hull = np.asarray(hull_image, dtype=bool)
    missing = hull & ~mask
    labels, count = ndimage.label(missing)
    hull_boundary = hull & ~ndimage.binary_erosion(hull)
    distance_to_hull = ndimage.distance_transform_edt(~hull_boundary)
    records: list[dict[str, float]] = []
    for label in range(1, count + 1):
        component = labels == label
        touches_mouth = bool(np.any(ndimage.binary_dilation(component) & hull_boundary))
        records.append({
            "kind": "concavity" if touches_mouth else "hole",
            "pixels": float(component.sum()),
            "depth": float(distance_to_hull[component].max()),
            "depth_ratio": float(distance_to_hull[component].max() / min(mask.shape)),
        })
    return sorted(records, key=lambda item: item["pixels"], reverse=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--size", type=int, default=1024)
    args = parser.parse_args()
    vertices, faces = read_obj(args.input)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for name, axes in (("xy", (0, 1)), ("xz", (0, 2)), ("yz", (1, 2))):
        mask = rasterize(vertices, faces, axes, args.size)
        holes, exterior, records = classify_background(mask)
        concavities = classify_concavities(mask)
        visualization = np.zeros((*mask.shape, 3), dtype=np.uint8)
        visualization[exterior] = (25, 28, 31)
        visualization[mask] = (190, 198, 202)
        visualization[holes] = (232, 91, 76)
        Image.fromarray(visualization).save(args.output_dir / f"{name}.png")
        print(name, "closed_holes", len(records), "largest", records[:12])
        print(name, "hull_missing", concavities[:16])


if __name__ == "__main__":
    main()
