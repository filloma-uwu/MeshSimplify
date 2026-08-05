# PQSS Real Scene Conservative Proxy Comparison

- Build strategy: `Optimized`
- Threads: 1
- Queries: 882008
- Base models: proxy `1-23` versus original `1-23`
- Appended workpieces: original geometry in both pools

## Build

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Base OBJ load ms |937.616800|1.846700|
| Base EndPool ms |1182713.745100|1239.692900|
| Base total ms |1183651.361900|1241.539600|
| Total including appends ms |1190610.140000|8220.971000|
| Pool memory bytes |371230896|6505536|
| Area threshold |3188826352605639680.000000|14711067562543117500416.000000|
| Cache hits |8|5|
| Cache misses |28|31|

### Appended Workpieces

| Model | Original append ms | Proxy-pool append ms |
|---:|---:|---:|
|24|3134.414300|3144.239000|
|26|9.096800|8.993500|
|27|1271.343700|1264.757200|
|28|5.119900|15.319100|
|29|0.267200|3.275700|
|30|0.130400|1.559700|
|31|3.255500|3.091600|
|32|0.190400|0.175000|
|33|1265.048900|1244.228300|
|34|0.226900|0.180500|
|35|1260.645500|1284.201900|
|36|0.183600|0.177400|
|37|8.812300|9.165700|

## Queries

| Metric | Original pool | Proxy pool |
|---|---:|---:|
| Parallel wall ms |249.094300|339.838800|
| Summed call ms |249.078500|339.826700|
| BV tests |2594432|6290438|
| Triangle tests |140359|254644|
| Positive results |16557|32252|

## Correctness

| TP | TN | FP | FN | False-positive rate |
|---:|---:|---:|---:|---:|
|16557|849756|15695|0|0.018135|

The false-positive denominator is `FP + TN`. A valid conservative result requires `FN = 0`.

Per-model geometry and BVH data are in `model_bvh_stats.csv`.
