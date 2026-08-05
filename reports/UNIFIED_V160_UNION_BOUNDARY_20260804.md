# Unified v160 union-boundary proxy report

## Result

The production proxy now emits the exterior boundary of the union of certified
closed envelopes. A planar region strictly crossed by another certified volume
is clipped even when the resulting boundary needs more triangles. Coplanar
union then removes coincident output patches.

The same Release binary and identical options generated all 13 evaluated
models. No branch uses a model ID, filename, scene pose, or collision result.

Manifest:

`outputs/surface_primitives_v160_union_boundary_complete/viewer_manifest.json`

## Selection policy

- Added volume remains normalized by the complete model-scale volume.
- Proxy workload uses the triangle count after containment removal, overlap
  clipping, and coplanar canonicalization.
- Added volume and the largest connected added component have bounded workload
  penalties. Closed-volume responsibility has a bounded workload credit, so
  overlap removal may buy some triangles but can never select an arbitrarily
  fragmented zero-error candidate.
- Model 5 therefore selects the zero-added-volume frontier instead of the old
  12-triangle candidate with 28.86% added volume.

## Complete batch

| Model | Source tris | Primitives | Proxy tris | Time (s) | Peak MB | Added volume | Clipped patches | Removed internal area |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 2 | 62,888 | 6 | 12 | 6.31 | 50.3 | 0.000% | 0 | 0.0 |
| 3 | 16,976 | 15 | 50 | 4.50 | 23.1 | 2.709% | 5 | 198,400.0 |
| 4 | 65,538 | 11 | 22 | 5.79 | 52.5 | 0.000% | 0 | 0.0 |
| 5 | 72,304 | 12 | 36 | 9.76 | 55.8 | 0.000% | 6 | 1,031,174.0 |
| 12 | 11,326 | 46 | 165 | 8.86 | 36.3 | 4.963% | 30 | 97,495.4 |
| 13 | 45,496 | 16 | 94 | 13.97 | 80.8 | 0.000% | 4 | 34,760.3 |
| 14 | 41,981 | 30 | 172 | 22.48 | 69.2 | 25.632% | 19 | 219,264.5 |
| 15 | 25,030 | 16 | 32 | 2.74 | 43.9 | 1.875% | 6 | 56,033.4 |
| 16 | 154,877 | 54 | 141 | 47.37 | 220.2 | 22.301% | 31 | 1,584,118.3 |
| 17 | 24,914 | 11 | 22 | 5.66 | 32.4 | 5.526% | 4 | 16,951.7 |
| 18 | 46,065 | 11 | 22 | 12.96 | 60.9 | 3.866% | 4 | 233,944.5 |
| 19 | 26,488 | 3 | 96 | 2.83 | 40.1 | 0.000% | 0 | 0.0 |
| 20 | 11,480 | 6 | 12 | 6.07 | 28.9 | 15.161% | 0 | 0.0 |

All models have `unassigned_source_faces=0` and `failed_source_faces=0`.

## Memory diagnosis

Model 14 originally exceeded 8 GB during overlap clipping. Per-candidate stage
profiling showed the growth occurred inside Clipper with only 753 primitives:
one rectangular solid had been recognized as the same extrusion along three
axes, producing many coincident plane sections. Certificates are now deduplicated
by their quantized 3D vertex sets and projected sections are deduplicated before
Boolean union.

After the fix, model 14 peaks at 69.2 MB. The complete four-process batch peaks
at 220.2 MB for any one process. The generator records `peak_memory_mb.txt` per
model and aborts the complete batch if one process exceeds 8 GB or system free
memory falls below 4 GB.

## Verification

- C++ primitive analyzer regression: passed.
- Python primitive-analysis suite: 19 passed.
- Playwright viewer suite: 3 passed.
- All 13 viewer models: nonblank canvas, selectable last primitive, and opacity
  reset after selection/model changes.
- Screenshots: `viewer/screenshots/v160`.
