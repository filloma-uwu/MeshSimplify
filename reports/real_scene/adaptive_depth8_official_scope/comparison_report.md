# PQSS Real Scene Conservative Proxy Comparison

- Build strategy: `Optimized`
- Threads: 1
- Queries: 741326
- Base models: proxy `1-23` versus original `1-23`
- Appended workpieces: original geometry in both pools

- Requested proxy max depth: 8
- Actual proxy max depth: 15
- Depth constraint satisfied: `false`

## Build

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Base OBJ load ms |933.271300|3.482100|
| Base EndPool ms |1384.400200|2440.393600|
| Base total ms |2317.671500|2443.875700|
| Total including appends ms |2359.551300|2480.766800|
| Pool memory bytes |371230896|7143240|
| Area threshold |3188826352605639680.000000|997506235010558459904.000000|
| Cache hits |36|28|
| Cache misses |0|8|

- Original cache tag: `original_91194221491100`
- Proxy cache tag: `adaptive_depth8_strict_20260723`

### Cache Phases

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Eligible models |36|36|
| Area threshold ms |2.259100|0.007900|
| Cache load ms |1398.742800|25.955000|
| Cache build ms |0.000000|2417.416200|
| Cache save ms |0.000000|17.534400|

### Appended Workpieces

| Model | Original append ms | Proxy-pool append ms |
|---:|---:|---:|
|24|8.562900|9.046700|
|26|9.305400|9.040500|
|27|2.960100|2.944400|
|28|4.646800|0.410000|
|29|0.167700|0.131600|
|30|0.139200|0.108700|
|31|0.139200|0.118400|
|32|0.111700|0.119500|
|33|3.134700|2.838400|
|34|0.179600|0.128000|
|35|3.262400|2.914200|
|36|0.169900|0.132100|
|37|9.082500|8.945500|

## Base Model And BVH Data

| Metric | Original `1-23` | Proxy `1-23` |
|---|---:|---:|
| Input triangles |627190|2005|
| Built triangles |629236|2131|
| Subdivisions |682|42|
| BV nodes |1258449|4239|
| Leaf nodes |629236|2131|
| Internal nodes |629213|2108|
| Maximum depth |59|15|
| Root size sum |1313317799907.206787|1313335925031.191895|
| Total size |13819613236236.109375|9628061963310.171875|
| Leaf size sum |2000290757799.567383|2000241286739.673340|
| Internal size sum |11819322478436.542969|7627820676570.504883|

## Queries

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Parallel wall ms |89.575500|116.233400|
| Summed call ms |89.553600|116.214700|
| BV tests |340952|931862|
| Triangle tests |7880|165173|
| Positive results |1881|19874|

## Correctness

| TP | TN | FP | FN | False-positive rate |
|---:|---:|---:|---:|---:|
|1881|721452|17993|0|0.024333|

The false-positive denominator is `FP + TN`. A valid conservative result requires `FN = 0`.

Per-model geometry and BVH data are in `model_bvh_stats.csv`; query-pair data are in `query_pair_stats.csv`.
