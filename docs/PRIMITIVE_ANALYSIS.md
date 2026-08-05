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
   fitted plane distance and edge connectivity. Single-boundary approximate
   planes become outward support polygons; this recovers one CAD surface despite
   decimal-coordinate noise without merging arbitrary nearby faces.
3. Classify circular planar domains as disks or annuli. A candidate takes
   responsibility only for source triangles inside its actual radial domain;
   islands inside an annular opening remain available to later analysis. Fit
   connected lateral domains as cylindrical or conical bands and verify their
   side normals. Their circumscribed tessellation is audited directly as a
   surface, including complete 360-degree bands.
4. Evaluate every inner boundary independently. Its conservative loop area
   times model depth is divided by the whole-model AABB volume. Small holes may
   close; a large body cavity remains. Repeated closed components with matching
   cross sections are analyzed in the same model frame; each gap is assigned
   its own volume ratio and may be bridged without consuming a global budget.
5. Merge every adjacent patch that has the same certified analytic surface;
   clipped patches retain their parameter-domain boundary, while a certified
   original complete surface may be restored.
6. Remove exact overlap and certified internal surfaces. Stage 3 is a
   surface-candidate search only: a planar, circular, cylindrical, or conical
   patch may merge only into another supported surface patch. A fitted box,
   cylinder solid, frustum solid, or the boundary of such a fitted volume is not
   an eligible simplification candidate. Volume certificates may prove that an
   existing surface is occluded, but cannot create replacement geometry.
   Every approximate surface replacement must pass the sampled directed
   simplified-to-hole-filled distance limit `maximum_open_error_distance`.
   Hole filling itself is not charged because phase 1 is the distance reference.
   PQSS workload and a target primitive count are not part of the geometry
   decision. Source faces still failing the final conservative audit remain exact
   surface triangles; repair never introduces a fitted solid.
7. Canonicalize the final outer surface and triangulate each semantic patch.

The approach combines the face-based hierarchical framework of Attene,
Falcidieno, and Spagnuolo (2006), *Hierarchical Mesh Segmentation Based on
Fitting Primitives*, with the adjacency collapse and added-volume principle of
Knodt and Gao (2026), *Convex Primitive Decomposition for Collision Detection*.

No watertightness, manifoldness, consistent winding, or self-intersection
property is required.

## Output

Each analyzed model directory contains:

- `source.obj`: source surface after numerical-degenerate filtering;
- `phase1_hole_filled.obj`: accepted holes sealed, their hidden inner walls
  removed, and over-budget cavities retained;
- `phase2_recognized_surfaces.obj`: adjacent certified analytic surfaces after
  exact same-surface merging;
- `primitives.obj`: phase 3 after overlap/internal removal and iterative
  Hausdorff-limited merging;
- `proxy.obj`: phase 4, the ordinary triangulated collision input;
- `regions.obj`: source triangles grouped by final responsibility;
- `model.json`: type counts, triangle counts, fill diagnostics, and timing.
- `intercomponent_gap_profile.json`: accepted per-gap bridge bounds and volume
  ratios when repeated closed components contain fillable spaces;
- `spatial_group_fixed_point_profile.json`: records that stage 3 is operating in
  surface-only mode and the requested directed-error limit;
- `coverage_audit_pre_repair.json`: exported workload before conservative
  fallback repair, the exact source-face IDs, and the reduced repair workload;
- `final_occlusion_certificate_profile.json`: historical versus active closed
  certificates used by final overlap removal;
- `stage_error_profile.jsonl`: primitive workload and measured maximum error at
  phase 2, phase 3, final cleanup, and conservative repair.

The maximum and mean simplification errors in `model.json` and the viewer are
directed from phase 3/4 to phase 1. Strict conservative coverage is still
audited independently against `source.obj`. Exact source-triangle safety
repairs are not approximation error and remain exact surface triangles.

Audit a completed batch, including finite coordinates, zero-area triangles,
exact duplicate triangles, containment, error limits, timings, and peak memory:

```powershell
python tools/audit_staged_surface_outputs.py `
  outputs/surface_primitives_stage3_complete_v4/viewer_manifest.json
```

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
