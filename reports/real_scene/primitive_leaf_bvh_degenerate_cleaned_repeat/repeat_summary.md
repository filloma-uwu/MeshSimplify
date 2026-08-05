# Degenerate-cleaned Primitive BVH Repeat Summary

- Queries per run: 741326
- Threads: 1
- Primitive BVH cache format: v2
- Primitive BVH cache hits per run: 13/13
- Original and official simplified baseline: standard PQSS Optimized builder/query
- Candidate: analytic primitive leaves with descendant-responsibility-triangle Optimized RSS internal nodes

The local PQSS copy includes the zero-length RSS core-interval fix: overlapping
one-dimensional side constraints are collapsed before corner classification.
Every generated internal RSS passed the existing source-vertex containment
certification during the full v2 cache rebuild.

| Run | Original ms | Official simplified ms | Primitive BVH ms |
|---:|---:|---:|---:|
|1|95.3605|91.3985|86.6153|
|2|93.2842|92.3685|88.4191|
|3|92.8125|91.7905|88.6976|
|4|93.2936|95.2211|86.9541|
|5|94.7777|92.6156|86.7955|
|Median|93.2936|92.3685|86.9541|

The primitive BVH median is 5.86% lower than the official simplified median
and 6.80% lower than the original-geometry median.

| Metric | Original | Official simplified | Primitive BVH |
|---|---:|---:|---:|
|Total BV tests|1082278|998144|1006532|
|Triangle/analytic-leaf tests|7880|3575|1975|
|False positives|0|48|94|
|False negatives|0|0|0|

Relative to original geometry, the primitive BVH reduces total BV tests by
7.00% and terminal tests by 74.94%. Relative to the official simplification,
it performs 0.84% more total BV tests but 44.76% fewer terminal tests.

The one-time uncached v2 primitive BVH rebuild took 253240 ms. Build time is
reported only and is not part of the optimization objective.
