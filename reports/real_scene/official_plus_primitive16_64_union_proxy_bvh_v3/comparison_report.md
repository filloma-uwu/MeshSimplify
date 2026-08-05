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
| Base OBJ load ms |943.683000|128.199100|
| Base EndPool ms |1399.712000|1530.399400|
| Base total ms |2343.395100|1658.598500|
| Total including appends ms |2384.777500|1697.531400|
| Pool memory bytes |371230896|53058936|
| Area threshold |3188826352605639680.000000|23996064645422276608.000000|
| Cache hits |36|35|
| Cache misses |0|1|

- Original cache tag: `original_91194221491100`
- Proxy cache tag: `official_simplified_proxy_bvh_v3_local_axes`

### Cache Phases

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Eligible models |36|36|
| Area threshold ms |2.175700|0.321800|
| Cache load ms |1413.289900|202.847900|
| Cache build ms |0.000000|1342.811500|
| Cache save ms |0.000000|5.143200|

### Appended Workpieces

| Model | Original append ms | Proxy-pool append ms |
|---:|---:|---:|
|24|9.525900|9.130000|
|26|8.667900|9.497100|
|27|2.950700|3.039600|
|28|5.061400|1.239500|
|29|0.195400|0.148400|
|30|0.117600|0.110800|
|31|0.124700|0.120200|
|32|0.107000|0.097900|
|33|2.866300|2.945900|
|34|0.126600|0.140100|
|35|2.853100|3.120500|
|36|0.126300|0.124000|
|37|8.644900|9.207000|

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
| Total size |13819613236236.109375|12357104769062.152344|
| Leaf size sum |2000290757799.567383|2000182857132.899414|
| Internal size sum |11819322478436.542969|10356921911929.292969|

## Queries

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Parallel wall ms |89.710300|86.574800|
| Summed call ms |89.652000|86.560700|
| BV tests |340952|262112|
| Triangle tests |7880|3145|
| Positive results |1881|1929|

## Correctness

| TP | TN | FP | FN | False-positive rate |
|---:|---:|---:|---:|---:|
|1881|739397|48|0|0.000065|

The false-positive denominator is `FP + TN`. A valid conservative result requires `FN = 0`.

Per-model geometry and BVH data are in `model_bvh_stats.csv`; query-pair data are in `query_pair_stats.csv`.
