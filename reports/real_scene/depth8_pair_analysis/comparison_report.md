# PQSS Real Scene Conservative Proxy Comparison

- Build strategy: `Optimized`
- Threads: 1
- Queries: 882008
- Base models: proxy `1-23` versus original `1-23`
- Appended workpieces: original geometry in both pools

## Build

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Base OBJ load ms |938.373100|1.968900|
| Base EndPool ms |1371.771000|4.852100|
| Base total ms |2310.144200|6.821000|
| Total including appends ms |2351.054000|43.513300|
| Pool memory bytes |371230896|6505536|
| Area threshold |3188826352605639680.000000|14711067562543117500416.000000|
| Cache hits |36|36|
| Cache misses |0|0|

- Original cache tag: `original_91194221491100`
- Proxy cache tag: `proxy_92385015748300`

### Appended Workpieces

| Model | Original append ms | Proxy-pool append ms |
|---:|---:|---:|
|24|8.803400|9.253700|
|26|8.732300|8.885500|
|27|3.213400|2.959400|
|28|4.551400|0.209400|
|29|0.176600|0.150400|
|30|0.128500|0.119800|
|31|0.131100|0.132300|
|32|0.125000|0.102900|
|33|2.830400|3.034900|
|34|0.130300|0.167100|
|35|3.251600|2.832100|
|36|0.181300|0.130700|
|37|8.640500|8.701500|

## Queries

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Parallel wall ms |237.377000|336.161800|
| Summed call ms |237.318100|336.148500|
| BV tests |2594432|6290438|
| Triangle tests |140359|254644|
| Positive results |16557|32252|

## Correctness

| TP | TN | FP | FN | False-positive rate |
|---:|---:|---:|---:|---:|
|16557|849756|15695|0|0.018135|

The false-positive denominator is `FP + TN`. A valid conservative result requires `FN = 0`.

Per-model geometry and BVH data are in `model_bvh_stats.csv`; query-pair data are in `query_pair_stats.csv`.
