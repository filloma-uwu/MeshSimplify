# Primitive-leaf BVH Real-scene Comparison

- Queries: 741326
- Threads: 1
- Truth: original source geometry with the official PQSS Optimized builder
- Official simplified: dataset `obj_simplified`, standard PQSS Optimized builder and query
- Candidate: one BVH per logical model, analytic primitives as leaves
- Internal nodes: PQSS Optimized RSS fitted to descendant responsibility triangles
- Build time is reported only and is not an optimization objective.

## Build

| Metric | Value |
|---|---:|
| Original pool build ms |2431.670900|
| Official simplified pool build ms |377.332100|
| Reference import ms |61.068800|
| Analyzed primitive BVH build ms |71.121000|
| Primitive BVH cache hits |13|
| Primitive BVH cache misses |0|
| Primitive BVH cache load ms |1.500800|
| Primitive BVH uncached build ms |0.000000|
| Primitive BVH cache save ms |0.000000|
| Install into PQSS ModelPool ms |4.922400|
| Reference pool memory bytes |371230920|
| Official simplified pool memory bytes |56422656|

## Query Work

PQSS's counter omits the model-pair root test, so one root test per query is shown separately.

| Metric | Original geometry | Official simplified | Primitive BVH |
|---|---:|---:|---:|
| Root BV tests |741326|741326|741326|
| Internal BV tests |340952|256818|265206|
| Total BV tests |1082278|998144|1006532|
| Triangle/leaf-pair tests |7880|3575|1975|
| Positives |1881|1929|1975|
| Parallel query wall ms |93.284200|92.368500|88.419100|

## Accuracy Against Original Geometry

| Candidate | TP | TN | FP | FN | False-positive rate |
|---|---:|---:|---:|---:|---:|
| Official simplified |1881|739397|48|0|0.000065|
| Primitive BVH |1881|739351|94|0|0.000127|

A valid conservative result requires `FN = 0`.

The candidate uses the same `ModelPool::Query` and `QueryRecurse` implementation as PQSS. Its leaf-test counter records analytic RSS leaf pairs rather than triangle pairs.
