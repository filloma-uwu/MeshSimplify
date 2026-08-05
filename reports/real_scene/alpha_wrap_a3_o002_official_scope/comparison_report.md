# PQSS Real Scene Conservative Proxy Comparison

- Build strategy: `Optimized`
- Threads: 1
- Queries: 741326
- Base models: proxy `1-23` versus original `1-23`
- Appended workpieces: original geometry in both pools

- Requested proxy max depth: 64
- Actual proxy max depth: 21
- Depth constraint satisfied: `true`

## Build

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Base OBJ load ms |934.289800|97.719200|
| Base EndPool ms |1411.847200|181172.762700|
| Base total ms |2346.137000|181270.481900|
| Total including appends ms |2387.522000|188222.133900|
| Pool memory bytes |371230896|66505800|
| Area threshold |3188826352605639680.000000|18731327208216809472.000000|
| Cache hits |36|8|
| Cache misses |0|28|

- Original cache tag: `original_91194221491100`
- Proxy cache tag: `alpha_wrap_a3_o002_v1`

### Cache Phases

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Eligible models |36|36|
| Area threshold ms |2.534800|0.401500|
| Cache load ms |1424.292400|24.875100|
| Cache build ms |0.000000|187558.912100|
| Cache save ms |0.000000|380.145600|

### Appended Workpieces

| Model | Original append ms | Proxy-pool append ms |
|---:|---:|---:|
|24|9.182500|3156.595400|
|26|8.932100|8.830300|
|27|2.906600|1259.802200|
|28|4.959300|1.359200|
|29|0.198600|0.164800|
|30|0.130700|0.147200|
|31|0.137200|3.117500|
|32|0.107300|0.160500|
|33|2.967600|1237.775600|
|34|0.123100|0.228300|
|35|2.872600|1274.166700|
|36|0.126900|0.198500|
|37|8.728000|9.080100|

## Base Model And BVH Data

| Metric | Original `1-23` | Proxy `1-23` |
|---|---:|---:|
| Input triangles |627190|106773|
| Built triangles |629236|107283|
| Subdivisions |682|170|
| BV nodes |1258449|214543|
| Leaf nodes |629236|107283|
| Internal nodes |629213|107260|
| Maximum depth |59|21|
| Root size sum |1313317799907.206787|1313317926707.854492|
| Total size |13819613236236.109375|12357017546320.994141|
| Leaf size sum |2000290757799.567383|2000173749672.437988|
| Internal size sum |11819322478436.542969|10356843796648.595703|

## Queries

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Parallel wall ms |101.691200|88.251000|
| Summed call ms |101.676800|88.240000|
| BV tests |340952|261416|
| Triangle tests |7880|3681|
| Positive results |1881|1905|

## Correctness

| TP | TN | FP | FN | False-positive rate |
|---:|---:|---:|---:|---:|
|1881|739421|24|0|0.000032|

The false-positive denominator is `FP + TN`. A valid conservative result requires `FN = 0`.

Per-model geometry and BVH data are in `model_bvh_stats.csv`; query-pair data are in `query_pair_stats.csv`.
