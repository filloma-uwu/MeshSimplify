# Primitive-leaf BVH Real-scene Comparison

- Queries: 741326
- Threads: 1
- Truth: original source geometry with the official PQSS Optimized builder
- Candidate: one BVH per logical model, analytic primitives as leaves
- Internal nodes: PQSS Optimized RSS fitted to descendant responsibility triangles
- Build time is reported only and is not an optimization objective.

## Build

| Metric | Value |
|---|---:|
| Reference pool build ms |2424.141800|
| Reference import ms |70.764400|
| Analyzed primitive BVH build ms |96.707300|
| Primitive BVH cache hits |13|
| Primitive BVH cache misses |0|
| Primitive BVH cache load ms |4.278200|
| Primitive BVH uncached build ms |0.000000|
| Primitive BVH cache save ms |0.000000|
| Install into PQSS ModelPool ms |5.740900|
| Reference pool memory bytes |371230920|

## Query Work

PQSS's counter omits the model-pair root test, so one root test per reference query is shown separately.

| Metric | Original geometry | Primitive BVH |
|---|---:|---:|
| Root BV tests |741326|741326|
| Internal BV tests |340952|275046|
| Total BV tests |1082278|1016372|
| Triangle/leaf-pair tests |7880|1975|
| Positives |1881|1975|
| Parallel query wall ms |96.268500|96.747100|

## Accuracy Against Original Geometry

| TP | TN | FP | FN | False-positive rate |
|---:|---:|---:|---:|---:|
|1881|739351|94|0|0.000127|

A valid conservative result requires `FN = 0`.

The candidate uses the same `ModelPool::Query` and `QueryRecurse` implementation as PQSS. Its leaf-test counter records analytic RSS leaf pairs rather than triangle pairs.
