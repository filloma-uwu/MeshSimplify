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
   original complete surface may be restored. The planar layer pass accepts
   only same-facing parallel patches. It may flatten a shallow parallel step
   within the directed-error limit, but it never projects differently oriented
   faces onto one support plane.
6. Build the stage-3 adjacency graph from both source-responsibility edges and
   the actual boundaries of the current stage-2 patches. Boundary contact is
   tested between line segments, so T junctions, partially overlapping edges,
   and an endpoint lying inside another edge do not disappear merely because
   the two OBJ patches lack an identical vertex. For every adjacent pair, fit
   supported surface replacements from the union of their unique
   responsibility vertices. Coplanar work remains a single polygon surface.
   Certified disks, annuli, cylindrical bands, and conical bands are atomic at
   this stage. A complete model or
   geometric component is first tested for a conservative revolved envelope.
   This requires either certified circular end rings or a circular projected
   convex silhouette; rectangular silhouettes therefore cannot become round.
   The candidate is always a closed side-plus-two-caps surface and must reduce
   triangle work and pass the user's directed-distance limit.

   There is no generic non-planar fallback. If adjacent patches cannot be
   certified as one of the supported analytic surface types, they remain as
   their current surface patches. In particular, neither an oriented box nor a
   three-dimensional convex hull may replace merely complex or slightly curved
   geometry. Pairwise analytic growth is reserved for at most 512 active groups
   to bound failed-fit work. There is no model identifier, filename, or
   coordinate-specific path.

   Stage 3 exports surface candidates only. Revolved envelopes are emitted as
   ordinary disk and band surfaces; there is no box, convex hull, or other
   volumetric primitive in the output interface. Candidate distance is
   certified first with the 1-Lipschitz property of point-to-mesh distance; an
   unresolved certificate falls back to the full deterministic surface
   sampler. Hole filling itself is not charged because phase 1 is the distance
   reference. PQSS workload and a target primitive count are not part of the
   geometry decision.

   Coverage auditing retains the pre-canonicalization owner of every source
   face. This is necessary because coplanar union may rewrite the final
   `source_faces` map even though the conservative surface union is unchanged.
   Source faces still failing that independent audit remain exact surface
   triangles; repair never introduces a fitted solid.
7. Remove composite-proxy interfaces from the actual collision shell. Accepted
   support protrusions emit only their exposed five faces; the side faces stop
   at the support plane, and their contact footprints are subtracted from the
   supporting polygon. These hidden support owners remain audit-only coverage
   certificates and are never exported to PQSS. Contact cutouts cannot be
   refilled by the error merge, and multi-contact supports are decomposed into
   local strips to avoid long, BVH-unfriendly triangles. The candidate remains
   local by construction: its longest fitted edge is limited to 18% of the
   whole-model AABB diagonal. Larger candidates stay in the ordinary polygon or
   certified round-surface path instead of becoming a box-shaped shell.
   Final volume-occlusion certificates are rebuilt only from explicit
   `enclosure_group` shells. Ungrouped planar patches are not speculatively
   combined into extrusions; doing so is quadratic on large CAD soups and can
   erase surfaces without a closed-envelope ownership proof.
8. Canonicalize the final outer surface and triangulate each semantic patch.
   Boolean-clipping fragments are deduplicated at model-relative tolerance;
   loops with fewer than three unique non-collinear points are zero-area
   artifacts and are discarded before occlusion. The final source-coverage
   audit remains authoritative after this cleanup.

For oversized phase-2 sets, stage 3 now changes search direction before the
quadratic adjacency loop. It evaluates a spatial hierarchy from the root down;
each node competes as a set of analytic surface patches, not as a volumetric
primitive. The current families are a conservative planar box shell (six
polygons) and a conservative surface of revolution (one cylindrical/conical
band plus two disks). Every candidate must reduce triangulated query work by at
least four times, pass the same directed proxy-to-phase-1 distance limit, and
retain a source-coverage certificate. If the node fails, it is split and the
same competition is repeated on its children. Candidate sets below the fixed
workload guard continue through the ordinary adjacency fixed point unchanged;
this preserves already good planar, disk, annulus, cylindrical-band, and
conical-band recognition. Production candidates never invoke QuickHull.

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
- `adjacent_envelope_group_profile.json`: adjacency count, accepted certified
  analytic surface merges, distance work, and final group workload;
- `top_down_surface_profile.json`: workload-guard decision plus planar,
  revolved-surface, and six-polygon shell candidate counts for the bounded
  top-down decomposition;
- `coverage_audit_pre_repair.json`: exported workload before conservative
  fallback repair, the exact source-face IDs, and the reduced repair workload;
- `final_occlusion_certificate_profile.json`: historical versus active closed
  certificates used by final overlap removal;
- `stage_error_profile.jsonl`: reserved stage-profile stream. Full directed
  distance is measured once by the final audit; intermediate full-surface
  scans are deliberately omitted because they do not participate in candidate
  acceptance and dominate runtime on large CAD meshes.

The maximum and mean simplification errors in `model.json` and the viewer are
directed from phase 3/4 to phase 1. Strict conservative coverage is still
audited independently against `source.obj`. Exact source-triangle safety
repairs are not approximation error and remain exact surface triangles.

Coverage audit is geometric first. A source triangle with no surviving
per-face owner is not repaired immediately: the audit first checks whether the
current proxy surfaces or active closed-volume certificates already cover it.
Only genuinely uncovered triangles enter the exact repair path. This keeps
top-down simplification from being undone by bookkeeping-only owner loss.

Audit a completed batch, including finite coordinates, zero-area triangles,
exact duplicate triangles, containment, error limits, timings, and peak memory:

```powershell
python tools/audit_staged_surface_outputs.py `
  outputs/surface_adjacent_fixed_point_full_v4/viewer_manifest.json
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
