# PQSSProxyMesh Algorithm Rules

These rules are mandatory for every simplification change in this repository.

1. Never branch on a model ID, file name, known test-scene role, or hand-picked
   coordinate. A fix must be expressed as a geometry rule and run unchanged on
   every input OBJ.
2. Phase 1 takes the original arbitrary OBJ triangle soup and emits
   `phase1_halfedge.bin` as its only formal geometry artifact. Any occupancy
   grid used by phase 1 is temporary internal state, never a persisted phase
   result. Phase 2 takes the frozen halfedge file as its sole geometry input;
   it must not rerun phase 1 or use a display OBJ, a voxel file, or a
   reconstructed mesh in its place.
3. False positives are allowed; false negatives are not. A source triangle may
   be removed only when a geometry certificate proves that the retained proxy
   union still covers its responsibility.
4. Primitive analysis emits surface primitives: planar polygons, disks, annuli,
   cylindrical bands, and conical bands. A box is six polygons, not a separate
   output primitive. A distinct triangulation stage converts all surface
   primitives to the ordinary triangle OBJ consumed by PQSS; an n-vertex simple
   polygon must produce exactly n-2 triangles without Steiner vertices.
5. Every candidate and the final proxy must satisfy the same user-provided
   directed proxy-to-phase1 Hausdorff limit. The limit is a feasibility
   threshold, not a local, cumulative, normalized, or global error budget.
6. Geometrically coplanar output patches must be unioned before triangulation.
   Internal edges and overlapping coplanar faces must not survive merely because
   the source OBJ uses different vertex indices or disconnected topology.
7. Never delete a proxy because it merely looks internal or unimportant. Delete
   it only after union-coverage certification; record the reason in statistics.
8. Offline generation must not use scene poses or collision-query outcomes.
   Geometry recognition and simplification granularity are decided only by
   geometric type certificates, conservative coverage, and the user-provided
   directed proxy-to-source Hausdorff limit. PQSS workload must not prune,
   rank, or select geometry candidates.
9. BV tests, triangle tests, BVH depth, and primitive count are post-generation
   evaluation diagnostics. They may motivate a later algorithm redesign, but
   they never override a geometric merge decision in the current run.
10. Every algorithm change must add a model-independent regression, run the full
    test suite, regenerate every evaluated model with identical options, and
    update one non-mixed visualization manifest. Do not compare outputs produced
    by different algorithm versions as if they were one uniform result.
11. Primitive IDs and triangulation diagonals are not semantic identities and
    may change after a planar union. Validate covered regions, overlap area,
    conservative source coverage, and query workload instead of preserving a
    screenshot's numeric primitive ID.
12. Geometry-specific recognizers may only propose certified candidates. They
    must never greedily consume faces, mutate shared responsibility, or choose
    the final result in a family-specific pass. Exact fallback, planar, curved,
    convex, and enclosing candidates must compete in one model-level search.
13. A result is not eligible for the current viewer or collision benchmark
    because one named model looks correct. The same binary and identical
    options must pass the complete selected-model batch. Any output with
    per-face fallback at industrial scale (hundreds or thousands of isolated
    triangle primitives) is an optimizer failure, not a valid fine result.
14. A watertight mesh may be built as an offline occupancy/error reference, but
    it is never a required input property or a final-output property. The final
    conservative OBJ must remain certified directly against source triangles.
15. The optimization objective is the minimum actual triangle count in the
    final OBJ after coplanar union, exposed-boundary extraction, and semantic
    primitive triangulation, subject to every hard coverage and directed
    Hausdorff constraint. The Hausdorff limit is a feasibility threshold, not
    a local, cumulative, or global error budget. Target reduction ratios, QEM,
    Chamfer distance, primitive count, and PQSS workload may not replace this
    objective or authorize/reject a geometric merge.
