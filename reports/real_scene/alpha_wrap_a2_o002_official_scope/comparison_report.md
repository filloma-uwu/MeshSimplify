# PQSS Real Scene Conservative Proxy Comparison

- Build strategy: `Optimized`
- Threads: 1
- Queries: 741326
- Base models: proxy `1-23` versus original `1-23`
- Appended workpieces: original geometry in both pools

- Requested proxy max depth: 64
- Actual proxy max depth: 22
- Depth constraint satisfied: `true`

## Build

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Base OBJ load ms |936.016000|167.648700|
| Base EndPool ms |1387.008200|340928.072300|
| Base total ms |2323.024200|341095.721000|
| Total including appends ms |2366.398700|348062.567000|
| Pool memory bytes |371230896|114557448|
| Area threshold |3188826352605639680.000000|10648436543311988736.000000|
| Cache hits |36|8|
| Cache misses |0|28|

- Original cache tag: `original_91194221491100`
- Proxy cache tag: `alpha_wrap_a2_o002_v1`

### Cache Phases

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Eligible models |36|36|
| Area threshold ms |2.472700|0.699800|
| Cache load ms |1401.205900|34.329600|
| Cache build ms |0.000000|346941.559700|
| Cache save ms |0.000000|644.581400|

### Appended Workpieces

| Model | Original append ms | Proxy-pool append ms |
|---:|---:|---:|
|24|10.062200|3138.429000|
|26|9.268800|9.166700|
|27|3.052000|1270.895300|
|28|4.744900|4.979500|
|29|0.162800|0.193900|
|30|0.115800|0.114800|
|31|0.130400|3.213500|
|32|0.114100|0.183200|
|33|3.010900|1262.430000|
|34|0.277200|0.263200|
|35|3.261400|1267.715700|
|36|0.228100|0.199200|
|37|8.926900|9.034900|

## Base Model And BVH Data

| Metric | Original `1-23` | Proxy `1-23` |
|---|---:|---:|
| Input triangles |627190|187821|
| Built triangles |629236|189867|
| Subdivisions |682|682|
| BV nodes |1258449|379711|
| Leaf nodes |629236|189867|
| Internal nodes |629213|189844|
| Maximum depth |59|22|
| Root size sum |1313317799907.206787|1313317958837.806396|
| Total size |13819613236236.109375|13819157910903.843750|
| Leaf size sum |2000290757799.567383|2000180303610.733398|
| Internal size sum |11819322478436.542969|11818977607293.111328|

## Queries

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Parallel wall ms |103.199800|87.085500|
| Summed call ms |103.182500|87.071700|
| BV tests |340952|269388|
| Triangle tests |7880|3469|
| Positive results |1881|1893|

## Correctness

| TP | TN | FP | FN | False-positive rate |
|---:|---:|---:|---:|---:|
|1881|739433|12|0|0.000016|

The false-positive denominator is `FP + TN`. A valid conservative result requires `FN = 0`.

Per-model geometry and BVH data are in `model_bvh_stats.csv`; query-pair data are in `query_pair_stats.csv`.
