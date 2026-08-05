# PQSS Real Scene Conservative Proxy Comparison

- Build strategy: `Optimized`
- Threads: 1
- Queries: 741326
- Base models: proxy `1-23` versus original `1-23`
- Appended workpieces: original geometry in both pools

- Requested proxy max depth: 64
- Actual proxy max depth: 49
- Depth constraint satisfied: `true`

## Build

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Base OBJ load ms |942.351400|140.077700|
| Base EndPool ms |1387.321700|160253.540300|
| Base total ms |2329.673100|160393.618100|
| Total including appends ms |2370.748600|168008.699000|
| Pool memory bytes |371230896|56422632|
| Area threshold |3188826352605639680.000000|22164839914464608256.000000|
| Cache hits |36|8|
| Cache misses |0|28|

- Original cache tag: `original_91194221491100`
- Proxy cache tag: `official_simplified_proxy_bvh_v2`

### Cache Phases

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Eligible models |36|36|
| Area threshold ms |2.144900|0.422700|
| Cache load ms |1400.342800|29.151200|
| Cache build ms |0.000000|167367.640600|
| Cache save ms |0.000000|336.648000|

### Appended Workpieces

| Model | Original append ms | Proxy-pool append ms |
|---:|---:|---:|
|24|8.839300|3605.275100|
|26|8.892800|8.782400|
|27|2.847000|1329.283300|
|28|4.810500|1.346300|
|29|0.191500|0.174500|
|30|0.121700|0.126000|
|31|0.131900|3.178300|
|32|0.115400|0.184900|
|33|3.140700|1309.170400|
|34|0.192900|0.179700|
|35|2.883900|1348.289300|
|36|0.132100|0.196000|
|37|8.761700|8.867800|

## Base Model And BVH Data

| Metric | Original `1-23` | Proxy `1-23` |
|---|---:|---:|
| Input triangles |627190|90233|
| Built triangles |629236|90743|
| Subdivisions |682|170|
| BV nodes |1258449|181463|
| Leaf nodes |629236|90743|
| Internal nodes |629213|90720|
| Maximum depth |59|49|
| Root size sum |1313317799907.206787|1313318275682.316895|
| Total size |13819613236236.109375|12326234945818.734375|
| Leaf size sum |2000290757799.567383|2000181976429.811035|
| Internal size sum |11819322478436.542969|10326052969388.958984|

## Queries

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Parallel wall ms |99.423800|90.349000|
| Summed call ms |99.401800|90.334000|
| BV tests |340952|406212|
| Triangle tests |7880|3144|
| Positive results |1881|1929|

## Correctness

| TP | TN | FP | FN | False-positive rate |
|---:|---:|---:|---:|---:|
|1881|739397|48|0|0.000065|

The false-positive denominator is `FP + TN`. A valid conservative result requires `FN = 0`.

Per-model geometry and BVH data are in `model_bvh_stats.csv`; query-pair data are in `query_pair_stats.csv`.
