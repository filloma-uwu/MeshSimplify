# PQSS Real Scene Conservative Proxy Comparison

- Build strategy: `Optimized`
- Threads: 1
- Queries: 882008
- Base models: proxy `1-23` versus original `1-23`
- Appended workpieces: original geometry in both pools

## Build

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Base OBJ load ms |935.483900|1.679700|
| Base EndPool ms |1366.882500|934.350600|
| Base total ms |2302.366400|936.030300|
| Total including appends ms |2343.268600|7946.527900|
| Pool memory bytes |371230896|6439368|
| Area threshold |3188826352605639680.000000|2072538860453929156608.000000|
| Cache hits |36|5|
| Cache misses |0|31|

- Original cache tag: `original_91194221491100`
- Proxy cache tag: `adaptive_v2_20260723`

### Appended Workpieces

| Model | Original append ms | Proxy-pool append ms |
|---:|---:|---:|
|24|8.956500|3147.549400|
|26|8.909700|8.783000|
|27|3.005500|1256.795200|
|28|4.556600|72.674300|
|29|0.162900|3.410400|
|30|0.115800|1.064300|
|31|0.130600|3.269800|
|32|0.107300|0.185900|
|33|2.848600|1240.876700|
|34|0.126700|0.172600|
|35|3.214700|1266.614200|
|36|0.173200|0.221000|
|37|8.582300|8.855600|

## Queries

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Parallel wall ms |248.756500|226.454100|
| Summed call ms |248.742800|226.441000|
| BV tests |2594432|2404452|
| Triangle tests |140359|185598|
| Positive results |16557|27298|

## Correctness

| TP | TN | FP | FN | False-positive rate |
|---:|---:|---:|---:|---:|
|16557|854710|10741|0|0.012411|

The false-positive denominator is `FP + TN`. A valid conservative result requires `FN = 0`.

Per-model geometry and BVH data are in `model_bvh_stats.csv`; query-pair data are in `query_pair_stats.csv`.
