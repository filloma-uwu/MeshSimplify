# Unified v165 near-rectangle report

## Change

The final planar stage now proposes the minimum-area enclosing rectangle for a
polygon only when all of the following model-independent certificates pass:

- the rectangle contains every original boundary vertex;
- it reduces triangulated workload;
- added area fits the model-global tiny-detail budget;
- local added area is at most 3% of the polygon area;
- every original boundary vertex is within 1% of the rectangle diagonal from
  the rectangle boundary;
- the polygon is not marked as a preserved cavity opening.

This converts shallow bumps around an otherwise rectangular CAD silhouette into
a conservative four-corner rectangle, while retaining structural hexagons and
deep or narrow recesses.

Manifest:

`outputs/surface_primitives_v165_near_rectangle_regularized/viewer_manifest.json`

## Result

- Model 13 primitive 0: 13 boundary vertices / 11 triangles to 4 vertices / 2 triangles.
- Model 13 total: 114 to 105 proxy triangles.
- Model 12: two shallow six-vertex outlines each reduce from 4 to 2 triangles.
- All other evaluated models retain their v164 primitive and triangle counts.
- All 13 models: `unassigned_source_faces=0`, `failed_source_faces=0`.
- Maximum measured process memory: 229.8 MB on model 16.

## Verification

- C++ primitive analyzer regression: passed, including paired shallow-rectangle
  acceptance and structural-hexagon rejection cases.
- Python suites: 25 passed, 3 subtests passed.
- Playwright viewer suite: 3 passed.
- Complete 13-model rebuild and coverage audit: passed.
