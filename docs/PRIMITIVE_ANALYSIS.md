# Conservative Primitive Analysis

This document describes the offline primitive-analysis stage. Collision-tree
construction and traversal are deliberately outside its current scope.

## Contract

For source triangle set `F`, the analyzer returns primitive regions `P_i` and
an explicitly excluded redundant-surface set `E` satisfying

```text
F = E disjoint_union (disjoint_union_i F(P_i)).
```

Every non-excluded source triangle occurs exactly once. Adjacent regions may
share boundary vertices and edges, but they cannot share a triangle. Surface
replacement is conservative: fitted planes are moved to their outward support
position, outer circular radii expand, and annular inner radii shrink. Hole
filling may only add collision surface.

## Representation

The semantic output types are:

- planar polygon;
- disk;
- annulus;
- cylindrical band;
- conical band.

They are surface patches, not closed physics-engine colliders. A separate stage
triangulates them to the ordinary OBJ consumed by PQSS. A simple n-gon produces
exactly n-2 triangles without Steiner vertices. Circular surfaces use the
configured circumferential segment count.

## Analysis Algorithm

1. Remove numerical zero-area faces and weld vertices by exact position.
2. Form bounded normal domains from stable, high-area faces, then split them by
   fitted plane distance and edge connectivity. This recovers one CAD surface
   despite decimal-coordinate noise without merging arbitrary nearby faces.
3. Classify circular planar domains as disks or annuli. Fit connected lateral
   domains as cylindrical or conical bands and verify their side normals.
4. Evaluate every inner boundary independently. Its conservative loop area
   times model depth is divided by the whole-model AABB volume. Small holes may
   close; a large body cavity remains.
5. Build a model-level hierarchy of surface and enclosing candidates. A merge
   is admissible only when the maximum sampled directed distance from its proxy
   surface to the source open surface does not exceed
   `maximum_open_error_distance`. The deterministic MeshLab-style sampling
   covers vertices, edges, and area-weighted face interiors; every sampled
   point uses an exact closest-triangle query. Added volume is not a merge
   criterion.
6. Among admissible candidates, minimize triangulated proxy workload, then
   canonicalize exact coplanar coverage and triangulate the semantic surfaces.

The approach combines the face-based hierarchical framework of Attene,
Falcidieno, and Spagnuolo (2006), *Hierarchical Mesh Segmentation Based on
Fitting Primitives*, with the adjacency collapse and added-volume principle of
Knodt and Gao (2026), *Convex Primitive Decomposition for Collision Detection*.

No watertightness, manifoldness, consistent winding, or self-intersection
property is required.

## Output

Each analyzed model directory contains:

- `source.obj`: source surface after numerical-degenerate filtering;
- `regions.obj`: source triangles grouped by final responsibility;
- `primitives.obj`: semantic surface groups for inspection;
- `proxy.obj`: the triangulated collision input;
- `model.json`: type counts, triangle counts, fill diagnostics, and timing.

The output root contains `viewer_manifest.json`. Group names use matching IDs:

```text
region_00012_annulus
primitive_00012_annulus
```

This makes region and primitive colors stable in the viewer.

## Current Limitations

- Candidate circular axes come from PCA and stable normal domains; the fit is
  not a global continuous optimum.
- Curved semantic surfaces are sampled in `proxy.obj`; PQSS currently consumes
  those triangles rather than analytic surface equations.
- Non-circular holes that fail the per-cavity fill test fall back to exact
  polygons/triangles rather than a general polygon-with-holes semantic type.
