# PQSS Real Scene Conservative Proxy Comparison

- Build strategy: `Optimized`
- Threads: 1
- Queries: 741326
- Base models: proxy `1-23` versus original `1-23`
- Appended workpieces: original geometry in both pools

- Requested proxy max depth: 8
- Actual proxy max depth: 19
- Depth constraint satisfied: `false`

## Build

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Base OBJ load ms |932.459800|137.760100|
| Base EndPool ms |1384.061600|192.740500|
| Base total ms |2316.521400|330.500600|
| Total including appends ms |2358.765600|367.800200|
| Pool memory bytes |371230896|56422632|
| Area threshold |3188826352605639680.000000|22164839914464608256.000000|
| Cache hits |36|36|
| Cache misses |0|0|

- Original cache tag: `original_91194221491100`
- Proxy cache tag: `official_simplified_v1`

### Cache Phases

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Eligible models |36|36|
| Area threshold ms |3.298800|0.338000|
| Cache load ms |1396.680000|213.992700|
| Cache build ms |0.000000|0.000000|
| Cache save ms |0.000000|0.000000|

### Appended Workpieces

| Model | Original append ms | Proxy-pool append ms |
|---:|---:|---:|
|24|9.288700|8.926700|
|26|9.142600|8.792800|
|27|3.267200|2.930800|
|28|4.876900|1.254300|
|29|0.184600|0.137800|
|30|0.137300|0.111900|
|31|0.150600|0.135000|
|32|0.124500|0.097500|
|33|2.945300|2.858600|
|34|0.127100|0.116100|
|35|3.206500|2.767500|
|36|0.196200|0.112300|
|37|8.579900|9.046700|

## Base Model And BVH Data

| Metric | Original `1-23` | Proxy `1-23` |
|---|---:|---:|
| Input triangles |627190|90233|
| Built triangles |629236|90743|
| Subdivisions |682|170|
| BV nodes |1258449|181463|
| Leaf nodes |629236|90743|
| Internal nodes |629213|90720|
| Maximum depth |59|19|
| Root size sum |1313317799907.206787|1313318275682.316895|
| Total size |13819613236236.109375|12357103958880.779297|
| Leaf size sum |2000290757799.567383|2000181976429.811035|
| Internal size sum |11819322478436.542969|10356921982451.007812|

## Queries

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Parallel wall ms |87.213800|84.243000|
| Summed call ms |87.198000|84.225400|
| BV tests |340952|256818|
| Triangle tests |7880|3575|
| Positive results |1881|1929|

## Correctness

| TP | TN | FP | FN | False-positive rate |
|---:|---:|---:|---:|---:|
|1881|739397|48|0|0.000065|

The false-positive denominator is `FP + TN`. A valid conservative result requires `FN = 0`.

Per-model geometry and BVH data are in `model_bvh_stats.csv`; query-pair data are in `query_pair_stats.csv`.
