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
| Base OBJ load ms |948.001900|129.621200|
| Base EndPool ms |1403.545300|4285.640800|
| Base total ms |2351.547300|4415.262000|
| Total including appends ms |2392.597700|4454.029300|
| Pool memory bytes |371230896|54032808|
| Area threshold |3188826352605639680.000000|23465089812654706688.000000|
| Cache hits |36|35|
| Cache misses |0|1|

- Original cache tag: `original_91194221491100`
- Proxy cache tag: `official_simplified_v1`

### Cache Phases

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Eligible models |36|36|
| Area threshold ms |2.599100|0.319800|
| Cache load ms |1415.984500|201.918300|
| Cache build ms |0.000000|4088.582000|
| Cache save ms |0.000000|13.732300|

### Appended Workpieces

| Model | Original append ms | Proxy-pool append ms |
|---:|---:|---:|
|24|8.902900|9.198800|
|26|8.919700|9.315800|
|27|3.022400|3.056700|
|28|4.675300|1.282600|
|29|0.154700|0.144400|
|30|0.113700|0.120500|
|31|0.133900|0.123400|
|32|0.100000|0.108000|
|33|2.866700|2.920100|
|34|0.117300|0.119600|
|35|3.027600|2.876000|
|36|0.464800|0.118300|
|37|8.536800|9.368900|

## Base Model And BVH Data

| Metric | Original `1-23` | Proxy `1-23` |
|---|---:|---:|
| Input triangles |627190|85233|
| Built triangles |629236|85743|
| Subdivisions |682|170|
| BV nodes |1258449|171463|
| Leaf nodes |629236|85743|
| Internal nodes |629213|85720|
| Maximum depth |59|21|
| Root size sum |1313317799907.206787|1313318422688.375732|
| Total size |13819613236236.109375|12357101573688.660156|
| Leaf size sum |2000290757799.567383|2000182217337.612793|
| Internal size sum |11819322478436.542969|10356919356351.085938|

## Queries

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Parallel wall ms |100.570500|85.491200|
| Summed call ms |100.553700|85.454200|
| BV tests |340952|315902|
| Triangle tests |7880|3575|
| Positive results |1881|1929|

## Correctness

| TP | TN | FP | FN | False-positive rate |
|---:|---:|---:|---:|---:|
|1881|739397|48|0|0.000065|

The false-positive denominator is `FP + TN`. A valid conservative result requires `FN = 0`.

Per-model geometry and BVH data are in `model_bvh_stats.csv`; query-pair data are in `query_pair_stats.csv`.
