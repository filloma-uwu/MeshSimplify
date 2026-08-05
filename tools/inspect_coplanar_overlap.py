from collections import defaultdict
from pathlib import Path
import sys

import numpy as np
from shapely.geometry import Polygon


path = Path(sys.argv[1])
vertices = []
groups = []
current = None
for line in path.read_text().splitlines():
    fields = line.split()
    if not fields:
        continue
    if fields[0] == "v":
        vertices.append(tuple(map(float, fields[1:4])))
    elif fields[0] == "g":
        current = {"id": int(fields[1].split("_")[1]), "indices": []}
        groups.append(current)
    elif fields[0] == "f":
        current["indices"] = [int(value) - 1 for value in fields[1:]]

planes = defaultdict(list)
for group in groups:
    points = np.array([vertices[index] for index in group["indices"]])
    normal = np.zeros(3)
    for first, second in zip(points, np.roll(points, -1, axis=0)):
        normal += np.cross(first, second)
    normal /= np.linalg.norm(normal)
    dominant = int(np.argmax(np.abs(normal)))
    if normal[dominant] < 0:
        normal = -normal
    distance = float(np.dot(normal, points[0]))
    key = tuple(np.round(normal, 7)) + (round(distance, 5),)
    axes = [axis for axis in range(3) if axis != dominant]
    polygon = Polygon(points[:, axes]).buffer(0)
    planes[key].append((group["id"], polygon, dominant, distance))

results = []
for items in planes.values():
    for first in range(len(items)):
        for second in range(first + 1, len(items)):
            area = items[first][1].intersection(items[second][1]).area
            if area > 1.0e-5:
                results.append((area, items[first][0], items[second][0],
                                items[first][2], items[first][3],
                                items[first][1].area, items[second][1].area))
for result in sorted(results, reverse=True):
    print(result)
print("overlap_pairs", len(results), "pairwise_overlap_area", sum(item[0] for item in results))

if len(sys.argv) > 3 and sys.argv[2] == "--plane-z":
    requested_z = float(sys.argv[3])
    for items in planes.values():
        for item in items:
            if item[2] == 2 and abs(item[3] - requested_z) <= 1.0e-4:
                print("plane_z", requested_z, "id", item[0], "area", item[1].area,
                      "vertices", len(item[1].exterior.coords) - 1)
    raise SystemExit

if len(sys.argv) > 2:
    selected = int(sys.argv[2])
    source = next(item for items in planes.values() for item in items if item[0] == selected)
    for items in planes.values():
        for item in items:
            if item[0] == selected or item[2] != source[2]:
                continue
            overlap = source[1].intersection(item[1]).area
            if overlap > 1.0e-5:
                print("projected", selected, item[0], "separation", abs(source[3] - item[3]),
                      "overlap", overlap, "selected_area", source[1].area,
                      "other_area", item[1].area)
    outward = [item for items in planes.values() for item in items
               if item[2] == source[2] and item[3] < source[3] and
               source[3] - item[3] <= 20.0]
    if outward:
        from shapely.ops import unary_union
        covered = source[1].intersection(unary_union([item[1] for item in outward])).area
        print("nearer_outer_union", "covered", covered, "selected_area", source[1].area,
              "ratio", covered / source[1].area,
              "ids", [item[0] for item in outward])
