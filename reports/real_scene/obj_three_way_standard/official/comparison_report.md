# PQSS Real Scene Conservative Proxy Comparison

- Build strategy: `Optimized`
- BVH builder: `PQSS standard`
- Threads: 1
- Queries: 741326
- Base models: proxy `1-23` versus original `1-23`
- Appended workpieces: original geometry in both pools

- Requested proxy max depth: 100
- Actual proxy max depth: 19
- Depth constraint satisfied: `true`

## Build

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Base OBJ load ms |965.402800|147.108900|
| Base EndPool ms |1441.007900|147526.606500|
| Base total ms |2406.410700|147673.715400|
| Total including appends ms |2449.649500|154430.060800|
| Pool memory bytes |371230920|56422656|
| Area threshold |3188826352605639680.000000|22164839914464608256.000000|
| Cache hits |36|8|
| Cache misses |0|28|

- Original cache tag: `original_91194221491100`
- Proxy cache tag: `official_simplified_standard_obj_v1`

### Cache Phases

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Eligible models |36|36|
| Area threshold ms |2.301000|0.422300|
| Cache load ms |1455.311600|31.504600|
| Cache build ms |0.000000|153765.688500|
| Cache save ms |0.000000|347.478400|

### Appended Workpieces

| Model | Original append ms | Proxy-pool append ms |
|---:|---:|---:|
|24|9.410400|3138.230000|
|26|9.078500|9.160500|
|27|3.004100|1167.901300|
|28|5.204500|1.419000|
|29|0.250600|0.167300|
|30|0.134500|0.116900|
|31|0.168000|3.325500|
|32|0.155900|0.197800|
|33|3.019100|1191.207500|
|34|0.135600|0.184200|
|35|3.059400|1234.780700|
|36|0.146300|0.222500|
|37|9.449600|9.399100|

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
| Root size sum |1313317799907.206787|1313318277783.705322|
| Total size |13819613236236.109375|12357104329873.414062|
| Leaf size sum |2000290757799.567383|2000181976429.811035|
| Internal size sum |11819322478436.542969|10356922353443.642578|

## Queries

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Parallel wall ms |93.336500|91.652600|
| Summed call ms |93.312000|91.633300|
| BV tests |340952|255900|
| Triangle tests |7880|3533|
| Positive results |1881|1929|

## Correctness

| TP | TN | FP | FN | False-positive rate |
|---:|---:|---:|---:|---:|
|1881|739397|48|0|0.000065|

The false-positive denominator is `FP + TN`. A valid conservative result requires `FN = 0`.

Per-model geometry and BVH data are in `model_bvh_stats.csv`; query-pair data are in `query_pair_stats.csv`.
