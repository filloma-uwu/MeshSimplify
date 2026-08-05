# PQSS Real Scene Conservative Proxy Comparison

- Build strategy: `Optimized`
- Threads: 1
- Queries: 882008
- Base models: proxy `1-23` versus original `1-23`
- Appended workpieces: original geometry in both pools

- Requested proxy max depth: 8
- Actual proxy max depth: 8
- Depth constraint satisfied: `true`

## Build

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Base OBJ load ms |960.791000|1.279200|
| Base EndPool ms |1403.275600|2.489900|
| Base total ms |2364.066700|3.769200|
| Total including appends ms |2404.986600|41.878800|
| Pool memory bytes |371230896|5932056|
| Area threshold |3188826352605639680.000000|8097165996745788751872.000000|
| Cache hits |36|36|
| Cache misses |0|0|

- Original cache tag: `original_91194221491100`
- Proxy cache tag: `adaptive_depth8_strict_20260723`

### Cache Phases

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Eligible models |36|36|
| Area threshold ms |2.181400|0.001900|
| Cache load ms |1416.354000|25.939000|
| Cache build ms |0.000000|0.000000|
| Cache save ms |0.000000|0.000000|

### Appended Workpieces

| Model | Original append ms | Proxy-pool append ms |
|---:|---:|---:|
|24|8.984500|9.548000|
|26|8.803200|9.261800|
|27|2.896500|3.038600|
|28|4.737800|0.206800|
|29|0.201500|0.128400|
|30|0.114300|0.111500|
|31|0.128300|0.139500|
|32|0.104200|0.105700|
|33|2.882000|3.123700|
|34|0.120800|0.182200|
|35|2.998700|3.144800|
|36|0.168700|0.143000|
|37|8.758300|8.960900|

## Base Model And BVH Data

| Metric | Original `1-23` | Proxy `1-23` |
|---|---:|---:|
| Input triangles |627190|247|
| Built triangles |629236|277|
| Subdivisions |682|10|
| BV nodes |1258449|531|
| Leaf nodes |629236|277|
| Internal nodes |629213|254|
| Maximum depth |59|8|
| Root size sum |1313317799907.206787|1313336005584.172852|
| Total size |13819613236236.109375|7627639749885.146484|
| Leaf size sum |2000290757799.567383|2000241704474.182373|
| Internal size sum |11819322478436.542969|5627398045410.963867|

## Queries

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Parallel wall ms |242.681300|221.939500|
| Summed call ms |242.665500|221.921800|
| BV tests |2594432|2075940|
| Triangle tests |140359|182707|
| Positive results |16557|49285|

## Correctness

| TP | TN | FP | FN | False-positive rate |
|---:|---:|---:|---:|---:|
|16557|832723|32728|0|0.037816|

The false-positive denominator is `FP + TN`. A valid conservative result requires `FN = 0`.

Per-model geometry and BVH data are in `model_bvh_stats.csv`; query-pair data are in `query_pair_stats.csv`.
