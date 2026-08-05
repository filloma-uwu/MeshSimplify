"""Apply the final C++ triangle export invariants to an existing viewer batch."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def clean_proxy(path: Path) -> tuple[int, dict[int, int]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    vertices = [tuple(map(float, fields[1:4]))
                for line in lines
                if (fields := line.split()) and fields[0] == "v"]
    lower = [min(vertex[axis] for vertex in vertices) for axis in range(3)]
    upper = [max(vertex[axis] for vertex in vertices) for axis in range(3)]
    diagonal = math.sqrt(sum((upper[axis] - lower[axis]) ** 2
                             for axis in range(3)))
    area_epsilon_squared = max(diagonal**4 * 1.0e-28, 1.0e-48)
    quantum = max(diagonal * 1.0e-10, 1.0e-12)
    output: list[str] = []
    seen: set[tuple[tuple[float, float, float], ...]] = set()
    group_id: int | None = None
    group_counts: dict[int, int] = {}
    triangle_count = 0
    for line in lines:
        fields = line.split()
        if not fields:
            output.append(line)
            continue
        if fields[0] == "v" and len(fields) >= 4:
            output.append(line)
            continue
        if fields[0] == "g":
            name = fields[1] if len(fields) > 1 else ""
            group_id = (int(name.split("_", 2)[1])
                        if name.startswith("primitive_") else None)
            output.append(line)
            continue
        if fields[0] != "f" or len(fields) != 4:
            output.append(line)
            continue
        indices = [int(field.split("/", 1)[0]) - 1 for field in fields[1:]]
        points = [vertices[index] for index in indices]
        first = tuple(points[1][axis] - points[0][axis] for axis in range(3))
        second = tuple(points[2][axis] - points[0][axis] for axis in range(3))
        cross = (first[1] * second[2] - first[2] * second[1],
                 first[2] * second[0] - first[0] * second[2],
                 first[0] * second[1] - first[1] * second[0])
        if sum(value * value for value in cross) <= area_epsilon_squared:
            continue
        key = tuple(sorted(tuple(int(round(value / quantum)) for value in point)
                           for point in points))
        if key in seen:
            continue
        seen.add(key)
        output.append(line)
        triangle_count += 1
        if group_id is not None:
            group_counts[group_id] = group_counts.get(group_id, 0) + 1
    path.write_text("\n".join(output) + "\n", encoding="utf-8")
    return triangle_count, group_counts


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    args = parser.parse_args()
    manifest_path = args.manifest.resolve()
    manifest = json.loads(manifest_path.read_text(encoding="ascii"))
    for model in manifest["models"]:
        metadata_path = (manifest_path.parent / model["metadata"]).resolve()
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        proxy_path = metadata_path.parent / metadata["proxy"]
        triangle_count, group_counts = clean_proxy(proxy_path)
        stats = metadata["stats"]
        stats["triangulated_proxy_triangles"] = triangle_count
        stats["proxy_triangles"] = triangle_count
        for component in metadata.get("proxy_components", []):
            component["triangulated_face_count"] = group_counts.get(
                int(component["id"]), 0)
        metadata_path.write_text(
            json.dumps(metadata, separators=(",", ":"), ensure_ascii=False) + "\n",
            encoding="utf-8")
        print(f"model={model['id']} proxy_triangles={triangle_count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
