"""Print connected-component OBB diagnostics for an OBJ mesh."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import trimesh


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("mesh", type=Path)
    parser.add_argument("--limit", type=int, default=40)
    parser.add_argument("--cluster-gap-ratio", type=float)
    parser.add_argument("--negative-protrusion", nargs=2, metavar=("AXIS", "PLANE"))
    args = parser.parse_args()

    mesh = trimesh.load_mesh(args.mesh, process=True, maintain_order=True)
    components = mesh.split(only_watertight=False)
    print(f"faces={len(mesh.faces)} vertices={len(mesh.vertices)} components={len(components)}")
    rows = []
    for component in components:
        obb = component.bounding_box_oriented
        source_area = float(component.area)
        obb_area = float(obb.area)
        rows.append(
            (
                len(component.faces),
                source_area,
                max(obb_area - source_area, 0.0) / max(source_area, 1.0e-30),
                float(obb.volume),
                component.bounds.mean(axis=0),
                obb.extents,
            )
        )
    rows.sort(reverse=True, key=lambda row: row[0])
    for index, (faces, area, excess, volume, center, extents) in enumerate(rows[: args.limit]):
        print(
            f"{index:3d} faces={faces:5d} area={area:12.1f} "
            f"area_excess={excess:9.6f} obb_volume={volume:12.1f} "
            f"center={center.round(1)} extents={extents.round(1)}"
        )

    if args.cluster_gap_ratio is not None:
        bounds = np.asarray([component.bounds for component in components])
        extents = bounds[:, 1] - bounds[:, 0]
        model_extents = mesh.bounds[1] - mesh.bounds[0]
        model_diagonal = float(np.linalg.norm(model_extents))
        areas = np.asarray([float(component.area) for component in components])
        dominant = (areas >= 0.03 * float(mesh.area)) | (
            (extents >= 0.3 * model_extents).sum(axis=1) >= 2
        )
        parent = np.arange(len(components))

        def find(index: int) -> int:
            while parent[index] != index:
                parent[index] = parent[parent[index]]
                index = int(parent[index])
            return index

        def union(first: int, second: int) -> None:
            first_root = find(first)
            second_root = find(second)
            if first_root != second_root:
                parent[second_root] = first_root

        maximum_gap = args.cluster_gap_ratio * model_diagonal
        active = np.flatnonzero(~dominant)
        for offset, first in enumerate(active):
            gap = np.maximum(
                np.maximum(bounds[first, 0] - bounds[active[offset + 1 :], 1],
                           bounds[active[offset + 1 :], 0] - bounds[first, 1]),
                0.0,
            )
            distances = np.linalg.norm(gap, axis=1)
            for second in active[offset + 1 :][distances <= maximum_gap]:
                union(int(first), int(second))
        groups: dict[int, list[int]] = {}
        for index in active:
            groups.setdefault(find(int(index)), []).append(int(index))
        group_rows = []
        for group in groups.values():
            union_mesh = trimesh.util.concatenate([components[index] for index in group])
            obb = union_mesh.bounding_box_oriented
            group_rows.append((
                len(group), len(union_mesh.faces), float(union_mesh.area),
                float(obb.area), union_mesh.bounds.mean(axis=0), obb.extents,
            ))
        group_rows.sort(reverse=True, key=lambda row: row[1])
        print(f"cluster_gap={maximum_gap:.3f} dominant={int(dominant.sum())} groups={len(group_rows)}")
        for index, (parts, faces, area, obb_area, center, group_extents) in enumerate(group_rows[: args.limit]):
            excess = max(obb_area - area, 0.0) / max(area, 1.0e-30)
            print(
                f"group {index:3d} parts={parts:4d} faces={faces:5d} "
                f"area_excess={excess:9.6f} center={center.round(1)} "
                f"extents={group_extents.round(1)}"
            )

    if args.negative_protrusion is not None:
        axis_name, plane_text = args.negative_protrusion
        axis = {"x": 0, "y": 1, "z": 2}[axis_name.lower()]
        plane = float(plane_text)
        bounds = np.asarray([component.bounds for component in components])
        selected = np.flatnonzero(bounds[:, 0, axis] < plane - 1.0e-6)
        projection_axes = [value for value in range(3) if value != axis]
        parent = np.arange(len(components))

        def find_projected(index: int) -> int:
            while parent[index] != index:
                parent[index] = parent[parent[index]]
                index = int(parent[index])
            return index

        def union_projected(first: int, second: int) -> None:
            first_root = find_projected(first)
            second_root = find_projected(second)
            if first_root != second_root:
                parent[second_root] = first_root

        model_diagonal = float(np.linalg.norm(mesh.bounds[1] - mesh.bounds[0]))
        maximum_gap = 0.01 * model_diagonal
        for offset, first in enumerate(selected):
            others = selected[offset + 1 :]
            gap = np.maximum(
                np.maximum(
                    bounds[first, 0, projection_axes] - bounds[others, 1][:, projection_axes],
                    bounds[others, 0][:, projection_axes] - bounds[first, 1, projection_axes],
                ),
                0.0,
            )
            distances = np.linalg.norm(gap, axis=1)
            for second in others[distances <= maximum_gap]:
                union_projected(int(first), int(second))
        groups: dict[int, list[int]] = {}
        for index in selected:
            groups.setdefault(find_projected(int(index)), []).append(int(index))
        protrusion_rows = []
        for group in groups.values():
            union_mesh = trimesh.util.concatenate([components[index] for index in group])
            obb = union_mesh.bounding_box_oriented
            protrusion_rows.append((
                len(group), len(union_mesh.faces), float(union_mesh.area),
                float(obb.area), union_mesh.bounds,
            ))
        protrusion_rows.sort(reverse=True, key=lambda row: row[1])
        print(f"negative_protrusions={len(selected)} groups={len(protrusion_rows)} gap={maximum_gap:.3f}")
        for index, (parts, faces, area, obb_area, group_bounds) in enumerate(protrusion_rows[: args.limit]):
            excess = max(obb_area - area, 0.0) / max(area, 1.0e-30)
            print(
                f"protrusion {index:3d} parts={parts:4d} faces={faces:5d} "
                f"area_excess={excess:9.6f} "
                f"center={group_bounds.mean(axis=0).round(1)} "
                f"extents={(group_bounds[1] - group_bounds[0]).round(1)}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
