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
| Base OBJ load ms |963.639300|1114.795100|
| Base EndPool ms |1431.270400|589674.204800|
| Base total ms |2394.909700|590788.999900|
| Total including appends ms |2438.118700|597475.417600|
| Pool memory bytes |371230920|202028328|
| Area threshold |3188826352605639680.000000|6214229343389104128.000000|
| Cache hits |36|8|
| Cache misses |0|28|

- Original cache tag: `original_91194221491100`
- Proxy cache tag: `primitive_mesh_obj_s041_standard_batched_v1`

### Cache Phases

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Eligible models |36|36|
| Area threshold ms |2.246600|1.135300|
| Cache load ms |1445.789600|111.284800|
| Cache build ms |0.000000|594517.351700|
| Cache save ms |0.000000|1232.485700|

### Appended Workpieces

| Model | Original append ms | Proxy-pool append ms |
|---:|---:|---:|
|24|9.286100|3094.944600|
|26|10.053700|9.797200|
|27|3.144800|1165.458300|
|28|4.874600|5.283600|
|29|0.164500|0.237300|
|30|0.111600|0.132200|
|31|0.130200|2.974700|
|32|0.104300|0.174400|
|33|2.997300|1189.345000|
|34|0.127200|0.347400|
|35|3.119300|1208.312400|
|36|0.202600|0.214600|
|37|8.872500|9.157600|

## Base Model And BVH Data

| Metric | Original `1-23` | Proxy `1-23` |
|---|---:|---:|
| Input triangles |627190|321842|
| Built triangles |629236|323888|
| Subdivisions |682|682|
| BV nodes |1258449|647753|
| Leaf nodes |629236|323888|
| Internal nodes |629213|323865|
| Maximum depth |59|32|
| Root size sum |1313317799907.206787|1313318713146.486084|
| Total size |13819613236236.109375|13818970930579.156250|
| Leaf size sum |2000290757799.567383|2000297382320.368896|
| Internal size sum |11819322478436.542969|11818673548258.789062|

## Queries

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Parallel wall ms |89.669000|93.798500|
| Summed call ms |89.650900|93.782200|
| BV tests |340952|351208|
| Triangle tests |7880|7869|
| Positive results |1881|1881|

## Correctness

| TP | TN | FP | FN | False-positive rate |
|---:|---:|---:|---:|---:|
|1881|739445|0|0|0.000000|

The false-positive denominator is `FP + TN`. A valid conservative result requires `FN = 0`.

Per-model geometry and BVH data are in `model_bvh_stats.csv`; query-pair data are in `query_pair_stats.csv`.
