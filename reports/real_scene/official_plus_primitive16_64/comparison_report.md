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
| Base OBJ load ms |948.238700|128.642800|
| Base EndPool ms |1409.705200|1131.893000|
| Base total ms |2357.943900|1260.535800|
| Total including appends ms |2398.552300|1298.659400|
| Pool memory bytes |371230896|52799880|
| Area threshold |3188826352605639680.000000|24085648566329577472.000000|
| Cache hits |36|35|
| Cache misses |0|1|

- Original cache tag: `original_91194221491100`
- Proxy cache tag: `official_simplified_v1`

### Cache Phases

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Eligible models |36|36|
| Area threshold ms |2.604400|0.307100|
| Cache load ms |1422.014000|202.804100|
| Cache build ms |0.000000|945.584900|
| Cache save ms |0.000000|4.372900|

### Appended Workpieces

| Model | Original append ms | Proxy-pool append ms |
|---:|---:|---:|
|24|8.892600|8.942500|
|26|8.691200|9.122600|
|27|2.944800|2.896200|
|28|4.742700|1.252500|
|29|0.188300|0.136800|
|30|0.121000|0.106200|
|31|0.124300|0.134400|
|32|0.104700|0.106000|
|33|3.058300|3.356200|
|34|0.154400|0.177700|
|35|2.860300|3.055200|
|36|0.130800|0.270900|
|37|8.580000|8.553800|

## Base Model And BVH Data

| Metric | Original `1-23` | Proxy `1-23` |
|---|---:|---:|
| Input triangles |627190|83037|
| Built triangles |629236|83547|
| Subdivisions |682|170|
| BV nodes |1258449|167071|
| Leaf nodes |629236|83547|
| Internal nodes |629213|83524|
| Maximum depth |59|19|
| Root size sum |1313317799907.206787|1313318649467.802734|
| Total size |13819613236236.109375|12357112972599.771484|
| Leaf size sum |2000290757799.567383|2000185046966.212158|
| Internal size sum |11819322478436.542969|10356927925633.597656|

## Queries

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Parallel wall ms |86.793300|86.814900|
| Summed call ms |86.730600|86.797700|
| BV tests |340952|262790|
| Triangle tests |7880|3575|
| Positive results |1881|1929|

## Correctness

| TP | TN | FP | FN | False-positive rate |
|---:|---:|---:|---:|---:|
|1881|739397|48|0|0.000065|

The false-positive denominator is `FP + TN`. A valid conservative result requires `FN = 0`.

Per-model geometry and BVH data are in `model_bvh_stats.csv`; query-pair data are in `query_pair_stats.csv`.
