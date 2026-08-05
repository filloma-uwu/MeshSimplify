# PQSS Real Scene Conservative Proxy Comparison

- Build strategy: `Optimized`
- BVH builder: `PQSS standard`
- Threads: 1
- Queries: 741326
- Base models: proxy `1-23` versus original `1-23`
- Appended workpieces: original geometry in both pools

- Requested proxy max depth: 100
- Actual proxy max depth: 32
- Depth constraint satisfied: `true`

## Build

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Base OBJ load ms |928.860900|339.143500|
| Base EndPool ms |1423.649100|414005.686400|
| Base total ms |2352.510100|414344.829900|
| Total including appends ms |2395.223800|420808.199400|
| Pool memory bytes |371230920|139399992|
| Area threshold |3188826352605639680.000000|8881783464497157120.000000|
| Cache hits |36|8|
| Cache misses |0|28|

- Original cache tag: `original_91194221491100`
- Proxy cache tag: `unified_analytic_sweep_v125`

### Cache Phases

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Eligible models |36|36|
| Area threshold ms |2.167600|0.825000|
| Cache load ms |1437.958700|48.968000|
| Cache build ms |0.000000|419189.881300|
| Cache save ms |0.000000|868.937800|

### Appended Workpieces

| Model | Original append ms | Proxy-pool append ms |
|---:|---:|---:|
|24|9.278100|2964.411400|
|26|9.217300|9.092400|
|27|3.107600|1119.770100|
|28|4.865000|4.915300|
|29|0.205500|0.163600|
|30|0.123900|0.261800|
|31|0.120300|3.176500|
|32|0.109600|0.153700|
|33|2.881100|1144.350200|
|34|0.153800|0.202400|
|35|3.235700|1207.576600|
|36|0.206100|0.200800|
|37|9.193400|9.068900|

## Base Model And BVH Data

| Metric | Original `1-23` | Proxy `1-23` |
|---|---:|---:|
| Input triangles |627190|225180|
| Built triangles |629236|227226|
| Subdivisions |682|682|
| BV nodes |1258449|454429|
| Leaf nodes |629236|227226|
| Internal nodes |629213|227203|
| Maximum depth |59|32|
| Root size sum |1313317799907.206787|1313320782720.861084|
| Total size |13819613236236.109375|13818424286026.982422|
| Leaf size sum |2000290757799.567383|2000193788399.974854|
| Internal size sum |11819322478436.542969|11818230497627.009766|

## Queries

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Parallel wall ms |95.198400|665.687400|
| Summed call ms |95.134000|665.628300|
| BV tests |340952|10930138|
| Triangle tests |7880|5125352|
| Positive results |1881|2432|

## Correctness

| TP | TN | FP | FN | False-positive rate |
|---:|---:|---:|---:|---:|
|1881|738894|551|0|0.000745|

The false-positive denominator is `FP + TN`. A valid conservative result requires `FN = 0`.

Per-model geometry and BVH data are in `model_bvh_stats.csv`; query-pair data are in `query_pair_stats.csv`.
