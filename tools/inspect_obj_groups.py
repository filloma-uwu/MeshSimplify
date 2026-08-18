from pathlib import Path
import numpy as np
import sys

path = Path(sys.argv[1])
vertices = []
groups = []
current = None
for line in path.read_text().splitlines():
    fields = line.split()
    if not fields:
        continue
    if fields[0] == "v":
        vertices.append(np.asarray(list(map(float, fields[1:4]))))
    elif fields[0] == "g":
        current = {"name": fields[1], "faces": []}
        groups.append(current)
    elif fields[0] == "f":
        current["faces"].append([int(value.split("/")[0]) - 1 for value in fields[1:]])

for index, group in enumerate(groups):
    ids = sorted({vertex for face in group["faces"] for vertex in face})
    points = np.asarray([vertices[vertex] for vertex in ids])
    normal = np.zeros(3)
    for face in group["faces"]:
        normal += np.cross(points[0] * 0 + vertices[face[1]] - vertices[face[0]],
                           vertices[face[2]] - vertices[face[0]])
    length = np.linalg.norm(normal)
    print(index, group["name"], "bounds", points.min(0).round(1),
          points.max(0).round(1), "normal", (normal / length).round(3) if length else normal)
