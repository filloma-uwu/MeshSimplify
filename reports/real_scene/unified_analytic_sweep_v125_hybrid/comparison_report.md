# PQSS Real Scene Conservative Proxy Comparison

- Build strategy: `Optimized`
- BVH builder: `PQSS standard`
- Threads: 1
- Queries: 741326
- Base models: proxy `1-23` versus original `1-23`
- Appended workpieces: original geometry in both pools

- Requested proxy max depth: 100
- Actual proxy max depth: 18
- Depth constraint satisfied: `true`

## Build

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Base OBJ load ms |930.690900|38.094800|
| Base EndPool ms |1438.095700|946.511600|
| Base total ms |2368.786700|984.606500|
| Total including appends ms |2411.586500|1031.229500|
| Pool memory bytes |371230920|19306008|
| Area threshold |3188826352605639680.000000|82318077071649538048.000000|
| Cache hits |36|26|
| Cache misses |0|10|

- Original cache tag: `original_91194221491100`
- Proxy cache tag: `official_simplified_v1`

### Cache Phases

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Eligible models |36|36|
| Area threshold ms |2.171100|0.132000|
| Cache load ms |1452.205200|106.600900|
| Cache build ms |0.000000|854.043100|
| Cache save ms |0.000000|16.433200|

### Appended Workpieces

| Model | Original append ms | Proxy-pool append ms |
|---:|---:|---:|
|24|9.582900|11.504300|
|26|9.166100|9.168400|
|27|3.070700|4.731700|
|28|4.999300|1.485800|
|29|0.170800|0.240600|
|30|0.112600|0.143000|
|31|0.121300|0.478100|
|32|0.100700|0.178200|
|33|2.954300|4.705200|
|34|0.124900|0.164000|
|35|2.900000|4.453600|
|36|0.467900|0.169900|
|37|9.010200|9.179000|

## Base Model And BVH Data

| Metric | Original `1-23` | Proxy `1-23` |
|---|---:|---:|
| Input triangles |627190|24296|
| Built triangles |629236|24806|
| Subdivisions |682|170|
| BV nodes |1258449|49589|
| Leaf nodes |629236|24806|
| Internal nodes |629213|24783|
| Maximum depth |59|18|
| Root size sum |1313317799907.206787|1313320817794.130127|
| Total size |13819613236236.109375|12356367863071.453125|
| Leaf size sum |2000290757799.567383|2000184570673.971680|
| Internal size sum |11819322478436.542969|10356183292397.519531|

## Queries

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Parallel wall ms |90.262200|197.931500|
| Summed call ms |90.242900|197.912500|
| BV tests |340952|2260210|
| Triangle tests |7880|899107|
| Positive results |1881|2432|

## Correctness

| TP | TN | FP | FN | False-positive rate |
|---:|---:|---:|---:|---:|
|1881|738894|551|0|0.000745|

The false-positive denominator is `FP + TN`. A valid conservative result requires `FN = 0`.

Per-model geometry and BVH data are in `model_bvh_stats.csv`; query-pair data are in `query_pair_stats.csv`.
