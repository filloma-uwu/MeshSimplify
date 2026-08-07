# PQSSProxyMesh

`PQSSProxyMesh` is an offline CPU project for generating a pool of independent
visual collision proxy meshes optimized for PQSS `BuildStrategy::Optimized`
BVHs.

The first supported mode is `ConservativeOuter`: the generated proxy is a
certified union of supported analytic surface patches. It may produce
false-positive collision results, but must not introduce false negatives on
the validated collision workload. Unsupported non-planar geometry remains a
surface mesh; it is never replaced by a generic box or convex solid.

## Fixed Decisions

- The simplifier input is a group of models matching one PQSS `ModelPool`.
- The original path keeps one physical proxy per logical model. The
  independent-primitive path maps one logical model to one or more physical
  PQSS models, one RSS BVH per primitive. Primitives never cross logical model
  boundaries.
- Candidate pools use PQSS's normal model-pool subdivision policy.
- The caller specifies a maximum PQSS Optimized BVH depth applied to every
  model in the pool.
- Build time is outside the optimization objective.
- Pool-level query cost and BVH quality take priority over minimum face count.

PQSS itself is the source of truth for subdivision. A candidate is loaded as a
fresh complete pool, subdivision remains enabled, and `EndPool()` computes the
shared threshold and builds every model. The simplifier does not maintain a
second implementation of that policy.

See [docs/DESIGN.md](docs/DESIGN.md) for the current design contract.

## Implemented Generator

The original adaptive outer generator uses 16-bucket world-axis SAH to
partition source triangles, then fits boxes, prisms, and local hulls. It is
retained as a baseline; its depth-derived face budget and minimum-volume
selection are not the target algorithm because they produced loose boxes and
poor scoped query work.

The new experimental `pqss-proxy-primitives` path implements the first stages
of PQSS-targeted conservative primitive decomposition: dense SAH
over-segmentation, area-normal quadrics, bottom-up spatial merging, certified
OBB/prism/hull fitting, and removal of internal primitive faces by CSG union.
Its implementation is intentionally allowed to be slow. Fixed part counts are
currently hierarchy checkpoints, not the final selection objective.

The complete target design and research basis are in
[docs/PQSS_TARGETED_ALGORITHM.md](docs/PQSS_TARGETED_ALGORITHM.md). Remaining
work is the PQSS static RSS/overlap inspector and pool-level Pareto selection;
until those are complete, the primitive path is experimental.

Models that already contain at most 12 faces are copied exactly. This avoids
turning open planes such as real-scene `21.obj` into large artificial solids;
exact pass-through is conservative but has zero outward clearance.

```powershell
$env:PYTHONPATH=(Resolve-Path python).Path
python -m pqss_proxy_mesh.cli `
  --input-dir test_data/real_scene/source_pool `
  --output-dir outputs/real_scene/conservative_outer_adaptive_depth8_strict `
  --max-depth 8 `
  --depth-safety-margin 4
```

Experimental primitive decomposition:

```powershell
$env:PYTHONPATH=(Resolve-Path python).Path
python -m pqss_proxy_mesh.primitive_cli `
  --input-dir test_data/real_scene/source_pool `
  --output-dir outputs/primitive_candidate
```

Historical independent primitive BLAS ablation (the directory name
`strength8` records an older fixed-eight checkpoint; it is not the current
`analysis-strength` parameter):

```powershell
$env:PYTHONPATH=(Resolve-Path python).Path
python -m pqss_proxy_mesh.primitive_cli `
  --input-dir test_data/real_scene/source_pool `
  --output-dir outputs/real_scene/independent_primitives_strength8 `
  --analysis-strength 0.5 `
  --independent-parts
```

This writes `primitive_pool_manifest.tsv`. The C++ loader expands the manifest
into one large PQSS `ModelPool`: the historical fixed-eight real-scene pool has 114
base physical models and 13 appended original workpieces, for 127 independent
RSS BVHs. PQSS computes one shared subdivision threshold over all 114 base
physical models before building each BVH.

Python 3.11+, NumPy, and SciPy are declared in `pyproject.toml`.

## Conservative Primitive Analysis

The C++ collision-proxy path separates semantic surface analysis from
tessellation. Analysis emits planar polygons, disks, annuli, cylindrical bands,
and conical bands. These are surfaces rather than closed physics colliders. It
writes `primitives.obj`, where an n-sided planar region remains one OBJ polygon,
then writes `proxy.obj` in a distinct triangulation stage for direct PQSS
consumption. Each simple n-gon produces exactly n-2 triangles; circular-surface
density is controlled by `--round-surface-segments`.

```powershell
.\build\Release\pqss-primitive-mesh-analyze.exe `
  --input test_data/real_scene/source_pool/4.obj `
  --output-dir outputs/model4 `
  --primitive-types polygon,surface `
  --round-surface-segments 24 `
  --maximum-open-error-distance 80 `
  --maximum-process-memory-gb 2
```

`--maximum-open-error-distance` is a sampled proxy-to-source distance limit in
the OBJ's model units. If omitted, the limit is 100 model units for every model
in the pool; it is not scaled by each model's AABB. Like MeshLab's surface
Hausdorff workflow, the C++ audit samples the
emitted proxy surface and computes every sample's exact closest point and
closest source triangle through a triangle BVH. Sampling is deterministic and
includes triangle vertices, uniformly spaced edge points, and area-weighted
face-interior points. A source triangle is excluded from the error statistic
only when rays in both normal directions are blocked by other source triangles,
which classifies it as an internal surface. `model.json` records the mean and
maximum accepted sample distance, the maximum point pair, sample count, and
containment result. The point-to-mesh calculation is exact; the maximum over
the continuous proxy surface remains a sampling estimate. PyMeshLab is not a
runtime dependency.

On Windows the analyzer installs a Job Object before loading the OBJ. Its
default process-commit limit is 2 GB; `--maximum-process-memory-gb` may lower or
retain that limit but cannot raise it. Allocation beyond the limit fails inside the
analyzer instead of exhausting system memory. The uniform batch launcher also
defaults to one model at a time, monitors both private bytes and working set,
and writes `peak_private_memory_mb.txt` and `peak_memory_mb.txt` for every
completed model.

Circular holes are evaluated before a surface is fragmented. Each hole's
conservative swept-volume upper bound is independently normalized by the
whole-model AABB volume. Accepted holes produce a disk; a retained concentric
hole produces an annulus. CAD coordinate quantization is handled by a
complete-surface best-fit plane certificate rather than by relaxing arbitrary
triangle pairs.

```powershell
$env:PYTHONPATH=(Resolve-Path python).Path
python -m pqss_proxy_mesh.primitive_analysis_cli `
  --input-dir test_data/real_scene/source_pool `
  --output-dir outputs/primitive_analysis `
  --models 1,21 `
  --primitive-types sphere,capsule,rss,triangle `
  --analysis-strength 0.5
```

Inspect the source mesh, its triangle-region partition, and the fitted
primitives with the standalone viewer:

```powershell
python tools/serve_primitive_viewer.py `
  --manifest outputs/spatial_group_gap_validation/viewer_manifest.json `
  --port 8091
```

Long full-pool runs can be resumed without recomputing completed models:

```powershell
python -m pqss_proxy_mesh.primitive_analysis_cli `
  --input-dir test_data/real_scene/source_pool `
  --output-dir outputs/primitive_analysis `
  --models 2,3,4,5,12,13,14,15,16,17,18,19,20 `
  --analysis-strength 0.5 `
  --resume
```

Open the URL printed by the server. The viewer provides source, regions,
primitives, overlay, and split modes; primitive-type filters; individual region
selection; opacity and wireframe controls; and sampled outward-deviation
statistics. Entering a non-negative maximum error and choosing `重新生成` runs
the C++ analyzer again, caches the result, and reloads it without changing the
camera. Switching display modes also preserves the current camera. See
[docs/PRIMITIVE_ANALYSIS.md](docs/PRIMITIVE_ANALYSIS.md) for
the algorithm and data contract.

## Real-Scene Validation

`pqss-real-scene-compare` compiles the actual PQSS CPU implementation. It
builds original and proxy pools with `BuildStrategy::Optimized`, appends the
original workpieces, runs the supplied real-scene pose file, and writes aggregate,
per-model, and per-model-pair statistics.

`pqss-independent-real-scene-compare` keeps the original and official data on
the reference one-model/one-BVH query path and uses the independent scheduler
only for the generated manifest. It reports PQSS internal BV tests, primitive
root RSS tests, their sum, triangle tests, containment-cache skips, and the
confusion matrix separately. The project compiles local copies of the PQSS
files that it modifies; the sibling `PQSS` repository is not changed.

`pqss-primitive-bvh-real-scene-compare` is the current candidate path. Each
analyzed logical model has one BVH whose leaves are analytic primitives.
Internal RSS nodes are fitted with PQSS Optimized directly to the union of the
source triangles assigned to descendant leaves. The fitted RSS is certified
but not inflated. Sphere and capsule leaves use degenerate RSS lengths, and
the resulting tree is installed into the normal PQSS `ModelPool`. Queries use
PQSS `Query`/`QueryRecurse` directly; leaf overlap is the final conservative
collision answer for analyzed models.

The expensive custom BVHs are cached one file per analyzed model under
`outputs/primitive_bvh_official_scope/bvh_cache`. Cache keys include the model
TSV contents, source OBJ contents, build strategy, scalar size, and cache format
version, so changing one model rebuilds only that model.

```powershell
build/Release/pqss-primitive-bvh-real-scene-compare.exe `
  test_data/real_scene/source_pool `
  ../PQSS/test_data/real_scene/obj_simplified `
  ../PQSS/test_data/real_scene/obj `
  outputs/primitive_bvh_official_scope/primitive_bvh_pool.tsv `
  test_data/real_scene/official_simplified_scope/collision_pairs.txt `
  reports/real_scene/primitive_leaf_bvh_official_scope `
  1 `
  original_91194221491100 `
  official_simplified_v1
```

The first fixed-eight result is in
`reports/real_scene/independent_primitives_strength8/comparison_report.md`.
It has zero false negatives, but direct primitive-pair scheduling performs
16,494,180 root RSS tests. This is a measured negative result: independent
BLASes remove centroid interleaving inside a monolithic BVH, but an unindexed
M-by-N root scan is too expensive and must gain a primitive-level broad phase.

The canonical comparison scope replaces only the 13 models changed by the
dataset's official MeshLab output: `2-5` and `12-20`. All other base models and
all appended workpieces retain original geometry. The fixed query file keeps
only poses involving at least one of those 13 models; it is generated once,
not filtered by the benchmark executable on every run:

```text
test_data/real_scene/official_simplified_scope/collision_pairs.txt
```

See `reports/real_scene/SCOPED_COMPARISON.md` for the current baseline and
`test_data/real_scene/official_simplified_scope/scope_manifest.json` for the
model set, counts, and content hashes. Older reports over all 882,008 poses are
historical and must not be used to rank new candidates.

The current canonical result and reproducibility details are in
[reports/real_scene/SCOPED_COMPARISON.md](reports/real_scene/SCOPED_COMPARISON.md).

## Build

```powershell
cmake -S . -B build
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

Both C++ and Python tests are registered with CTest.
