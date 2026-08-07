"""Report overlap between coplanar semantic surface primitives.

Every plane bucket is projected through one shared orthonormal frame.  The old
diagnostic independently dropped the dominant coordinate of each primitive;
near a 45-degree normal, floating-point ties could make two primitives choose
different coordinates and manufacture a false overlap.
"""

from collections import defaultdict
from pathlib import Path
import sys

import numpy as np
from shapely.geometry import Polygon
from shapely.ops import unary_union


path = Path(sys.argv[1])
vertices: list[tuple[float, float, float]] = []
groups: list[dict[str, object]] = []
current: dict[str, object] | None = None
for line in path.read_text().splitlines():
    fields = line.split()
    if not fields:
        continue
    if fields[0] == "v":
        vertices.append(tuple(map(float, fields[1:4])))
    elif fields[0] == "g":
        name = fields[1]
        parts = name.split("_")
        current = {"id": int(parts[1]), "kind": parts[-1], "faces": []}
        groups.append(current)
    elif fields[0] == "f" and current is not None:
        current["faces"].append(
            [int(value.split("/", 1)[0]) - 1 for value in fields[1:]])


def canonical_plane(group: dict[str, object]):
    # Bands are not planar; treating their final triangle as a polygon was the
    # other source of misleading overlap reports.
    if group["kind"] in {"cylindricalband", "conicalband"}:
        return None
    faces = group["faces"]
    if not faces:
        return None
    normal = np.zeros(3)
    point = None
    for face in faces:
        if len(face) < 3:
            continue
        first = np.asarray(vertices[face[0]])
        for index in range(1, len(face) - 1):
            second = np.asarray(vertices[face[index]])
            third = np.asarray(vertices[face[index + 1]])
            normal += np.cross(second - first, third - first)
        point = first
    length = np.linalg.norm(normal)
    if point is None or length <= 1.0e-14:
        return None
    normal /= length
    dominant = int(np.argmax(np.abs(normal)))
    if normal[dominant] < 0.0:
        normal = -normal
    distance = float(np.dot(normal, point))
    key = tuple(np.round(normal, 7)) + (round(distance, 5),)
    return key, normal, distance


plane_groups: dict[tuple[float, ...], list[dict[str, object]]] = defaultdict(list)
for group in groups:
    plane = canonical_plane(group)
    if plane is None:
        continue
    key, normal, distance = plane
    group["normal"] = normal
    group["distance"] = distance
    plane_groups[key].append(group)


def plane_frame(normal: np.ndarray):
    helper = np.array([1.0, 0.0, 0.0])
    if abs(float(np.dot(helper, normal))) > 0.9:
        helper = np.array([0.0, 1.0, 0.0])
    first = np.cross(normal, helper)
    first /= np.linalg.norm(first)
    second = np.cross(normal, first)
    return first, second


results = []
projected_groups: dict[int, tuple[Polygon, np.ndarray, float]] = {}
for items in plane_groups.values():
    normal = items[0]["normal"]
    first_axis, second_axis = plane_frame(normal)
    projected = []
    for group in items:
        polygons = []
        for face in group["faces"]:
            points = np.asarray([vertices[index] for index in face])
            coordinates = [
                (float(np.dot(point, first_axis)),
                 float(np.dot(point, second_axis)))
                for point in points
            ]
            polygon = Polygon(coordinates).buffer(0)
            if not polygon.is_empty and polygon.area > 1.0e-12:
                polygons.append(polygon)
        polygon = unary_union(polygons).buffer(0)
        projected.append((group["id"], polygon))
        projected_groups[group["id"]] = (
            polygon, normal, float(group["distance"]))
    for first in range(len(projected)):
        for second in range(first + 1, len(projected)):
            area = projected[first][1].intersection(projected[second][1]).area
            if area > 1.0e-5:
                results.append((area, projected[first][0], projected[second][0],
                                projected[first][1].area,
                                projected[second][1].area))

for result in sorted(results, reverse=True):
    print(result)
print("overlap_pairs", len(results),
      "pairwise_overlap_area", sum(item[0] for item in results))

if len(sys.argv) > 2:
    selected = int(sys.argv[2])
    source = projected_groups[selected]
    for item_id, (polygon, normal, distance) in projected_groups.items():
        if item_id == selected or abs(float(np.dot(normal, source[1]))) < 1.0 - 1.0e-7:
            continue
        # Reprojecting non-coplanar polygons is intentionally omitted here;
        # this optional mode is only a plane-separation diagnostic.
        separation = abs(distance - source[2])
        if separation <= 20.0:
            print("parallel", selected, item_id, "separation", separation)
