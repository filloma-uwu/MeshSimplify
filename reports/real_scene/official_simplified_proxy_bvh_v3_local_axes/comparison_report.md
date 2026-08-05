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
| Base OBJ load ms |951.456800|141.632400|
| Base EndPool ms |1402.028800|147195.838200|
| Base total ms |2353.485600|147337.470700|
| Total including appends ms |2394.059400|154305.379300|
| Pool memory bytes |371230896|56422632|
| Area threshold |3188826352605639680.000000|22164839914464608256.000000|
| Cache hits |36|8|
| Cache misses |0|28|

- Original cache tag: `original_91194221491100`
- Proxy cache tag: `official_simplified_proxy_bvh_v3_local_axes`

### Cache Phases

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Eligible models |36|36|
| Area threshold ms |2.171800|0.379400|
| Cache load ms |1415.042200|29.683000|
| Cache build ms |0.000000|153649.405000|
| Cache save ms |0.000000|348.260700|

### Appended Workpieces

| Model | Original append ms | Proxy-pool append ms |
|---:|---:|---:|
|24|8.719500|3147.574300|
|26|9.003000|9.055000|
|27|2.847100|1262.783300|
|28|4.603500|1.437400|
|29|0.152000|0.146000|
|30|0.136600|0.109500|
|31|0.140100|3.506300|
|32|0.100500|0.165100|
|33|3.343600|1241.398400|
|34|0.189800|0.190700|
|35|2.842200|1292.502800|
|36|0.132400|0.228000|
|37|8.346800|8.786900|

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
| Total size |13819613236236.109375|12357103276996.212891|
| Leaf size sum |2000290757799.567383|2000181976429.811035|
| Internal size sum |11819322478436.542969|10356921300566.439453|

## Queries

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Parallel wall ms |99.434600|85.492200|
| Summed call ms |99.413500|85.427700|
| BV tests |340952|256080|
| Triangle tests |7880|3145|
| Positive results |1881|1929|

## Correctness

| TP | TN | FP | FN | False-positive rate |
|---:|---:|---:|---:|---:|
|1881|739397|48|0|0.000065|

The false-positive denominator is `FP + TN`. A valid conservative result requires `FN = 0`.

Per-model geometry and BVH data are in `model_bvh_stats.csv`; query-pair data are in `query_pair_stats.csv`.
