# PQSS Real Scene Conservative Proxy Comparison

- Build strategy: `Optimized`
- Threads: 1
- Queries: 741326
- Base models: proxy `1-23` versus original `1-23`
- Appended workpieces: original geometry in both pools

- Requested proxy max depth: 64
- Actual proxy max depth: 19
- Depth constraint satisfied: `true`

## Build

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Base OBJ load ms |941.375200|129.597800|
| Base EndPool ms |1398.269600|1559.176100|
| Base total ms |2339.644900|1688.773900|
| Total including appends ms |2380.589500|1727.355600|
| Pool memory bytes |371230896|53058936|
| Area threshold |3188826352605639680.000000|23996064645422276608.000000|
| Cache hits |36|35|
| Cache misses |0|1|

- Original cache tag: `original_91194221491100`
- Proxy cache tag: `official_simplified_v1`

### Cache Phases

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Eligible models |36|36|
| Area threshold ms |2.231500|0.319400|
| Cache load ms |1411.295200|203.449700|
| Cache build ms |0.000000|1371.437100|
| Cache save ms |0.000000|5.342300|

### Appended Workpieces

| Model | Original append ms | Proxy-pool append ms |
|---:|---:|---:|
|24|8.846800|9.236500|
|26|8.727900|8.704500|
|27|2.899800|3.082700|
|28|5.042100|1.589000|
|29|0.210300|0.180500|
|30|0.165900|0.150000|
|31|0.131300|0.137500|
|32|0.112200|0.133700|
|33|2.863400|3.039900|
|34|0.131500|0.170600|
|35|2.851600|3.030500|
|36|0.126800|0.130100|
|37|8.821500|8.981100|

## Base Model And BVH Data

| Metric | Original `1-23` | Proxy `1-23` |
|---|---:|---:|
| Input triangles |627190|83347|
| Built triangles |629236|83857|
| Subdivisions |682|170|
| BV nodes |1258449|167691|
| Leaf nodes |629236|83857|
| Internal nodes |629213|83834|
| Maximum depth |59|19|
| Root size sum |1313317799907.206787|1313318655363.675781|
| Total size |13819613236236.109375|12357105415468.939453|
| Leaf size sum |2000290757799.567383|2000182857132.899414|
| Internal size sum |11819322478436.542969|10356922558336.078125|

## Queries

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Parallel wall ms |92.617500|85.947100|
| Summed call ms |92.597000|85.932900|
| BV tests |340952|262508|
| Triangle tests |7880|3575|
| Positive results |1881|1929|

## Correctness

| TP | TN | FP | FN | False-positive rate |
|---:|---:|---:|---:|---:|
|1881|739397|48|0|0.000065|

The false-positive denominator is `FP + TN`. A valid conservative result requires `FN = 0`.

Per-model geometry and BVH data are in `model_bvh_stats.csv`; query-pair data are in `query_pair_stats.csv`.
