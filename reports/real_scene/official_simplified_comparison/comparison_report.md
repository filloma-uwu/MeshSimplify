# PQSS Real Scene Conservative Proxy Comparison

- Build strategy: `Optimized`
- Threads: 1
- Queries: 882008
- Base models: proxy `1-23` versus original `1-23`
- Appended workpieces: original geometry in both pools

- Requested proxy max depth: 64
- Actual proxy max depth: 19
- Depth constraint satisfied: `true`

## Build

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Base OBJ load ms |934.159800|140.823500|
| Base EndPool ms |1385.504400|148085.095600|
| Base total ms |2319.664300|148225.919100|
| Total including appends ms |2360.676800|155176.062000|
| Pool memory bytes |371230896|56422632|
| Area threshold |3188826352605639680.000000|22164839914464608256.000000|
| Cache hits |36|8|
| Cache misses |0|28|

- Original cache tag: `original_91194221491100`
- Proxy cache tag: `official_simplified_v1`

### Cache Phases

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Eligible models |36|36|
| Area threshold ms |3.152100|0.330900|
| Cache load ms |1397.215200|29.290300|
| Cache build ms |0.000000|154536.731100|
| Cache save ms |0.000000|334.945700|

### Appended Workpieces

| Model | Original append ms | Proxy-pool append ms |
|---:|---:|---:|
|24|8.794400|3147.253900|
|26|8.932700|8.948900|
|27|3.001500|1248.558900|
|28|4.691800|1.302100|
|29|0.156000|0.145300|
|30|0.114300|0.106200|
|31|0.141700|2.972700|
|32|0.102800|0.172900|
|33|2.855900|1253.415100|
|34|0.121100|0.191400|
|35|2.989200|1277.579700|
|36|0.137100|0.210400|
|37|8.953400|9.242800|

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
| Parallel wall ms |251.139100|235.545500|
| Summed call ms |251.122700|235.528000|
| BV tests |2594432|2510298|
| Triangle tests |140359|136054|
| Positive results |16557|16605|

## Correctness

| TP | TN | FP | FN | False-positive rate |
|---:|---:|---:|---:|---:|
|16557|865403|48|0|0.000055|

The false-positive denominator is `FP + TN`. A valid conservative result requires `FN = 0`.

Per-model geometry and BVH data are in `model_bvh_stats.csv`; query-pair data are in `query_pair_stats.csv`.
