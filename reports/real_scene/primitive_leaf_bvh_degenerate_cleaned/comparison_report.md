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
| Original pool build ms |2413.951400|
| Official simplified pool build ms |379.541000|
| Reference import ms |58.405600|
| Analyzed primitive BVH build ms |253239.877200|
| Primitive BVH cache hits |2|
| Primitive BVH cache misses |11|
| Primitive BVH cache load ms |0.725100|
| Primitive BVH uncached build ms |253150.373500|
| Primitive BVH cache save ms |18.163200|
| Install into PQSS ModelPool ms |4.710900|
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
| Parallel query wall ms |107.074600|91.953300|85.763600|

## Accuracy Against Original Geometry

| Candidate | TP | TN | FP | FN | False-positive rate |
|---|---:|---:|---:|---:|---:|
| Official simplified |1881|739397|48|0|0.000065|
| Primitive BVH |1881|739351|94|0|0.000127|

A valid conservative result requires `FN = 0`.

The candidate uses the same `ModelPool::Query` and `QueryRecurse` implementation as PQSS. Its leaf-test counter records analytic RSS leaf pairs rather than triangle pairs.
