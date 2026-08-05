# PQSS Real Scene Conservative Proxy Comparison

- Build strategy: `Optimized`
- Threads: 1
- Queries: 882008
- Base models: proxy `1-23` versus original `1-23`
- Appended workpieces: original geometry in both pools

## Build

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Base OBJ load ms |957.376800|1.056100|
| Base EndPool ms |1374.833500|217.115700|
| Base total ms |2332.210300|218.171900|
| Total including appends ms |2373.664200|7210.383500|
| Pool memory bytes |371230896|5932056|
| Area threshold |3188826352605639680.000000|8097165996745788751872.000000|
| Cache hits |36|5|
| Cache misses |0|31|

- Original cache tag: `original_91194221491100`
- Proxy cache tag: `adaptive_depth8_strict_20260723`

### Appended Workpieces

| Model | Original append ms | Proxy-pool append ms |
|---:|---:|---:|
|24|9.246100|3162.702200|
|26|8.604600|8.847900|
|27|3.377200|1256.130400|
|28|4.774400|15.575700|
|29|0.185300|3.338900|
|30|0.113700|1.217700|
|31|0.121800|3.731900|
|32|0.108200|0.191200|
|33|2.858400|1242.690800|
|34|0.129600|0.192700|
|35|2.829200|1287.571300|
|36|0.123200|0.218000|
|37|8.968300|9.771100|

## Queries

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Parallel wall ms |241.926200|223.746100|
| Summed call ms |241.863300|223.728800|
| BV tests |2594432|2075940|
| Triangle tests |140359|182707|
| Positive results |16557|49285|

## Correctness

| TP | TN | FP | FN | False-positive rate |
|---:|---:|---:|---:|---:|
|16557|832723|32728|0|0.037816|

The false-positive denominator is `FP + TN`. A valid conservative result requires `FN = 0`.

Per-model geometry and BVH data are in `model_bvh_stats.csv`; query-pair data are in `query_pair_stats.csv`.
