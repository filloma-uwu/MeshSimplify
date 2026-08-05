# PQSSProxyMesh Design

## Goal

Generate a pool of independent visual proxy meshes that:

1. preserves logical model order and boundaries while allowing a logical
   model to own multiple physical primitive models;
2. conservatively contains each corresponding source model;
3. builds high-quality PQSS Optimized RSS BVHs using normal pool subdivision;
4. keeps every physical primitive BVH at or below an optional requested
   maximum depth;
5. runs offline on the CPU, with build time excluded from the objective.

## Subdivision Policy

Subdivision must behave exactly as it does when the production application
builds a normal PQSS `ModelPool`. It is enabled by default and PQSS itself
performs the operation; the simplifier must not duplicate the rule.

For a fresh complete candidate pool, PQSS currently:

1. computes `area_metric = length_squared(cross(p2-p1, p3-p1))` for every
   proxy triangle in every model;
2. sets one shared threshold to the arithmetic mean of that metric over the
   whole pool;
3. recursively splits a triangle into four midpoint triangles while its area
   metric is greater than the shared threshold;
4. independently builds one Optimized BVH per model from the resulting
   triangles.

If a model starts with `F` proxy faces and PQSS performs `S` recursive split
operations, its built triangle/leaf count is `F + 3*S`, and a non-empty binary
BVH has `2*(F + 3*S) - 1` nodes. These values and the actual depth must still be
read from PQSS diagnostics, rather than inferred in production code.

The entire base candidate pool, including all physical primitives, must be
created before calling `EndPool()`. Appending
to an already processed PQSS pool reuses the existing threshold and does not
rebuild older models, so it is not equivalent to evaluating the pool as a
group.

Changing one proxy can change the shared threshold and therefore the leaf
counts and depths of every other proxy. Candidate evaluation and acceptance
are consequently atomic at pool level.

## Independent Primitive BLAS Path

The independent path stores every fitted primitive as a separate OBJ and a
separate physical model in one PQSS `ModelPool`. A manifest maps each logical
model ID to its physical model IDs. A query between logical models with `M`
and `N` primitives considers at most `M*N` root pairs and stops globally on
the first triangle pair within tolerance.

After PQSS builds the physical BVHs, the project-local adapter restores every
parent-relative RSS node to model coordinates. For each node it records:

- its ordinary ancestors in the same BLAS, based on triangle-set inclusion;
- RSS nodes in other BLASes of the same logical model whose swept volume
  completely contains it.

During one query, a cache stores only node pairs rejected by the RSS distance
test. A current pair `(u,v)` can be skipped when the cache contains `(a,b)`
with `u` covered by `a` and `v` covered by `b`. The cache is pair-specific and
query-local; a negative triangle traversal is never reused as a pruning
certificate. This preserves early-stop behavior and does not infer separation
from overlapping triangle shells.

The strength-8 real-scene experiment constructs 127 physical models and
records 13,024 cross-tree containment links. It safely skips 3,722,102 node
pairs and has zero false negatives, but still executes 16,494,180 primitive
root RSS tests. The next design requirement is therefore a primitive-level
broad phase; increasing primitive count without one is not viable.

## PQSS-Aware Quality Rules

Each proxy should preserve spatial partitionability rather than optimize face
count alone:

- do not bridge disconnected components or large empty gaps;
- bound triangle diameter, maximum edge length, and aspect ratio separately;
- keep triangle area and sampling density reasonably uniform so PCA-based RSS
  fitting is not biased by tessellation density;
- prefer compact, normal-coherent patches that produce low-radius RSS nodes;
- reject operations that substantially increase local SAH cost, RSS size, or
  sibling overlap.

## Conservative Adaptive Outer Algorithm

The implemented path does not use one convex hull for an entire model.

1. Source triangles are recursively partitioned with the same world-axis,
   16-bucket SAH cost form used by PQSS.
2. Every partition evaluates expanded OBB, oriented 4/6/8-sided prism, and
   local convex-hull candidates.
3. A dynamic program chooses one candidate per partition under the face budget
   derived from the requested depth.
4. Every source triangle remains assigned to one convex proxy part. Checking
   its three vertices against that part's half-spaces certifies the complete
   triangle because the proxy part is convex.
5. The complete proxy pool is built by the real PQSS Optimized builder. Its
   measured depth, not the face-budget estimate, is the acceptance criterion.

Already-simple models are exact pass-through candidates. This is important for
open planes: adding artificial thickness can dramatically increase root
overlap even though it is nominally an outer approximation.

## Candidate Evaluation

> Superseded for new development by `PQSS_TARGETED_ALGORITHM.md`. In
> particular, measured query traces are held-out acceptance data and are not
> available to the offline generator, and minimum volume is not the primary
> selection objective.

Every complete candidate pool is rebuilt by the production PQSS Optimized
builder with subdivision enabled. Per-model containment and maximum depth are
hard constraints. Candidate ranking uses, in priority order:

1. measured PQSS query time on representative traces;
2. mean and P99 BV test counts;
3. query-pair false-positive rate and root candidate rate;
4. per-model root RSS size, internal RSS size, and sibling overlap;
5. post-subdivision triangle test count and memory;
6. total proxy face count.

The initial static score, before representative query traces are available, is:

```text
pool_score =
    sum_over_models(
        w_root * root_rss_size
      + w_internal * depth_weighted_internal_rss_size
      + w_overlap * sibling_overlap
      + w_nodes * post_subdivision_node_count)
```

## Depth Contract

The public contract is:

```text
for every model: actual_max_depth <= requested_max_depth
```

Exact equality is not required. The actual depth and subdivision statistics of
every model are reported in `ProxyPoolMetrics`.

## Implementation Status

1. Model-pool OBJ input/output: complete.
2. PQSS CPU Optimized adapter and build diagnostics: complete.
3. Adaptive local outer proxy generation: complete.
4. Per-part conservative containment certification: complete.
5. Real-scene full-pose comparison and false-positive statistics: complete.
6. Automatic iterative feedback from measured PQSS depth into generation:
   currently performed by the depth safety margin and external validation;
   direct in-process feedback remains future work.
