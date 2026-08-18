# PQSS-Targeted Conservative Primitive Decomposition

> **Historical design only.** Its PQSS-workload objective is superseded by
> repository rule 15 and `PHASE2_SIMPLIFICATION_RESEARCH.md`: the active
> objective is minimum final triangulated OBJ face count under hard conservative
> coverage and directed Hausdorff constraints. PQSS workload is post-generation
> diagnostic data and must not drive candidate selection.

## Objective

The generator runs offline before query poses are known. It must not execute
or optimize against collision queries. Real-scene BV and triangle test counts
are held-out acceptance measurements, not generator inputs.

The output for one source model is a triangle soup representing a union of
independent convex primitives. Primitives may overlap. The complete output is
not required to be watertight, manifold, consistently connected, or free of
self-intersections.

Hard constraints are:

1. every source triangle is assigned to at least one output primitive;
2. the assigned convex primitive contains all three triangle vertices;
3. containment is certified from the primitive half-spaces with positive
   numerical clearance;
4. models outside the selected simplification set pass through exactly;
5. an optional PQSS depth limit rejects a pool but is not an optimization
   objective.

The primary objective is low PQSS query work. Since poses are unavailable,
the generator minimizes intrinsic surrogates derived from the production PQSS
builder. Face count, depth, volume, and build time are not primary objectives.

## Research Basis

### Convex Primitive Decomposition

Knodt and Gao, *Convex Primitive Decomposition for Collision Detection*
(2026), arXiv:2602.07369, initializes a covering primitive per source face and
greedily merges adjacent primitives. Each fitted primitive expands to enclose
all assigned vertices. The method directly supports non-manifold,
non-watertight, multi-component meshes and uses OBBs, spheres, cylinders,
capsules, frustums, and trapezoidal prisms. Its merge cost is added primitive
volume, inspired by QEM. This supplies the coverage invariant and bottom-up
primitive hierarchy used here.

### Split, Merge, and Refine

Park and Sung, *Split, Merge, and Refine: Fitting Tight Bounding Boxes via
Over-Segmentation and Iterative Search* (3DV 2024),
DOI:10.1109/3DV62453.2024.00146, demonstrates that tight complete coverage
benefits from over-segmentation followed by hierarchical merging and local
search. This supplies the search structure and motivates retaining multiple
candidate decompositions instead of committing to one greedy path.

### Error-Bounded Simplification

Cohen et al., *Simplification Envelopes* (SIGGRAPH 1996),
DOI:10.1145/237170.237220, constrains simplification between geometric
envelopes and measures error in both directions. The present algorithm keeps
the conservative outer-coverage side as a hard invariant and uses one-way
proxy-to-source distance/empty volume to reject loose candidates.

### Bounding Proxies

Calderon and Boubekeur, *Bounding Proxies for Shape Approximation* (TOG
2017), DOI:10.1145/3072959.3073714, argues that complex, multi-component input
must be regularized while preserving strict enclosure before coarse bounding
optimization. The present algorithm avoids their global closing because PQSS
does not require a single watertight proxy, but retains their strict bounding
contract.

### PQSS Construction

PQSS `BuildStrategy::Optimized` partitions triangle centroids with a
world-axis, 16-bin SAH. Its split cost is child AABB surface area multiplied by
child triangle count. Query traversal expands the non-leaf RSS with larger
`BV::Size()`, where:

```text
RSSSize = length0 * length1
        + pi * radius * (length0 + length1 + 2 * radius)
```

These exact definitions, rather than a generic face-count objective, define
the static workload surrogate below.

## Algorithm

### 1. Weld And Over-Segment

Vertices are welded by a scale-relative coordinate tolerance solely to recover
face adjacency. Source coordinates are never moved. Degenerate faces remain
covered but do not contribute normals.

Each welded connected component is over-segmented into compact, normal-
coherent patches. Candidate splits include:

- topology boundaries and disconnected components;
- large spatial gaps;
- normal discontinuities;
- the production 16-bin world-axis SAH split;
- local principal-axis splits for curved or diagonal parts.

Splitting continues well past the eventual output complexity. This prevents a
large early box from bridging holes or distant parts.

### 2. Fit Conservative Primitive Families

For every patch, an area-weighted face-normal quadric supplies an orientation.
The following candidates are fitted and expanded to enclose every assigned
source vertex:

- OBB;
- oriented 4/6/8/12-sided capped prisms;
- circumscribed cylinder and capsule meshes;
- oriented 14/18/26-DOPs;
- a bounded-face local convex hull.

Candidate certification evaluates all assigned vertices against every convex
half-space. A primitive that fails positive-clearance certification is
discarded. Since a convex primitive contains all three vertices, it contains
the complete source triangle.

The original experiment merged independent shells into one PQSS model, with
an optional CSG union. Both forms still produced poor BV traversal on complex
`16.obj`: separate shells interleaved their centroids, while the CSG union did
not sufficiently reduce sibling overlap. That path remains a negative
baseline.

The current experiment instead sends every fitted primitive to PQSS as a
separate physical model with its own RSS BVH. Logical model boundaries are
maintained by a manifest and query scheduler. Cross-BLAS RSS containment is
certified at build time and used to reuse query-local RSS separation results.
This removes centroid interleaving between primitives, but introduces up to
`M*N` root RSS tests for logical models with `M` and `N` primitives.

### 3. Bottom-Up Merge Hierarchy

The initial graph contains topologically adjacent patches plus a small number
of spatial-neighbor edges. A spatial edge is allowed only when its gap and
added empty volume are bounded, preventing arbitrary merges between remote
components.

Merging two nodes unions their assigned faces and normal quadrics, refits all
primitive families, and adds the best certified candidates to a Pareto set.
The merge hierarchy uses a beam rather than one irrevocable greedy path.

No recognizer owns source triangles. A cylinder detector, coplanar detector,
box fitter, or cavity proposal can only add a certified candidate covering a
face set. Candidate families never run as sequential cleanup passes and never
erase responsibility needed by later families. The model-level optimizer is
the only code allowed to select candidates.

### 4. CPU Watertight Error Reference

The input and final proxy remain arbitrary triangle soups. Watertightness is
used only inside the offline error estimator. Following PaMO stage 1, source
triangles are conservatively rasterized into a padded isotropic grid, expanded
by a three-cell narrow band, and closed by exterior flood fill. Isotropic cells
give the closing band one model-space scale, so repeated gaps such as sleepers
are not accidentally protected only because they lie on a short AABB axis.
Three-axis paired ray intervals supplement the flood fill for thin or
disconnected CAD exports.
This is a CPU implementation and does not depend on PaMO's CUDA/Warp runtime.

The grid may certify that an already existing surface is hidden by a closed
region, but it never creates a closed box, cylinder, or frustum candidate.
Candidate geometry is surface-only.  A box-shaped result, when it arises from
the input geometry, is represented by six independently certified polygon
faces; a round result is represented by disk/annulus faces and cylindrical or
conical side faces.  The final OBJ comes only from triangulating those surfaces,
not from marching cubes or from decomposing a fitted solid.

### 4.1 Analytic Sweep Atoms

Approximate planar analysis proposes circular disk or annulus surfaces.
Non-planar analysis proposes cylindrical or conical side surfaces directly.
Compatible disconnected CAD patches may be united only when they certify the
same axis, radius law, axial domain, and complementary angular domain.  The
candidate remains one side surface; end disks are separate candidates with
their own responsibility.  No sweep volume exists in the candidate pool.

### 5. PQSS Static Work Vector

Every candidate state is scored by a vector, not a single weighted face or
volume number:

```text
W = (
    predicted post-subdivision leaf triangles,
    PQSS 16-bin SAH sum,
    depth-weighted RSS Size sum,
    sibling RSS overlap,
    proxy-to-source one-way distance,
    added empty volume
)
```

The first four terms estimate triangle and BV work. The last two prevent the
optimizer from obtaining a cheap BVH by filling holes or distant gaps.
Dominated states are removed. Fixed primitive counts are hierarchy checkpoints,
not objectives. Depth is reported and may reject a state when a
caller supplies a hard limit, but shallower depth does not outrank lower
predicted work.

All error coordinates in the frontier are normalized by the complete input
model. A candidate-local empty fraction is forbidden: it makes an identical
geometric error change meaning when an unrelated hierarchy split changes the
candidate bounds, and in practice causes open CAD soups to fall back to tens
of thousands of triangle primitives.

Sibling overlap is mandatory, not optional. Experiments on `16.obj` showed
that lower root and total RSS size can still produce substantially more BV
tests when triangle centroids interleave and both children overlap. Candidate
selection must therefore reproduce the production centroid partition and
measure its child overlap rather than infer quality from primitive count,
volume, depth, or aggregate RSS alone.

### 6. Pool-Level Search

PQSS subdivision uses one area threshold shared by the full model pool.
Candidate selection is therefore performed on complete pools. At hierarchy
checkpoints the production CPU `Optimized` builder is run without queries to
obtain exact subdivision counts, RSS statistics, and depth. Coordinate descent
then replaces one model candidate at a time while retaining the pool-level
Pareto frontier.

No real-scene pose or collision result is available to this search.

### 7. Held-Out Acceptance

After a candidate pool is frozen, the fixed official-simplification query set
is run once. Ranking uses, in order:

1. zero false negatives;
2. total BV tests;
3. total triangle tests;
4. false positives;
5. query wall time.

Build/generation time is reported but never used for acceptance.

## Implementation Plan

1. Primitive fitting and half-space certification.
2. Welded adjacency and over-segmentation.
3. Bottom-up merge hierarchy with Pareto/beam retention.
4. PQSS static pool inspector for SAH, RSS, overlap, and subdivision data.
5. Pool-level candidate search.
6. Fixed-scope real-scene acceptance run.

## Independent BLAS Ablation

The first fixed-scope run was a historical fixed-eight checkpoint: it forced
eight physical primitives for each of the 13 selected models. This predates the
current `[0,1]` boundary-sensitivity meaning of `analysis_strength` and must not
be interpreted as the current analyzer's behavior. Unchanged base models and
all workpieces remained one physical model each. The resulting pool contained
127 physical models and 28,391 RSS nodes.

| Metric | Official simplified | Independent strength 8 |
|---|---:|---:|
| PQSS internal BV tests | 256,818 | 592,024 |
| Root/scheduler RSS tests | 741,326 | 16,494,180 |
| Triangle tests | 3,575 | 109,792 |
| Containment-cache skips | 0 | 3,722,102 |
| False positives | 48 | 9,231 |
| False negatives | 0 | 0 |

The result rejects direct all-pairs primitive scheduling as the final query
architecture. Node containment caching is correct and measurable, but most
root pairs are disjoint rather than nested. A primitive-level TLAS or another
sublinear broad phase is required before higher analysis strengths can be
competitive. The result is recorded in
`reports/real_scene/independent_primitives_strength8/comparison_report.md`.
