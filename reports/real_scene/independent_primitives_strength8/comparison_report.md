# Independent Primitive RSS BVH Comparison

- Queries: 741326
- Threads: 1
- Truth: original geometry, original one-model/one-BVH query
- Official: official simplified geometry, original one-model/one-BVH query
- Candidate: independent primitive BLASes, cross-tree RSS containment pruning
- Candidate builder: project-local six-axis Optimized builder

## Build

| Metric | Original | Official | Candidate |
|---|---:|---:|---:|
| Base build ms |2354.261700|337.071500|26.582500|
| Total build ms |2395.899600|374.591200|105.765400|
| Containment analysis ms |0|0|42.378700|
| PQSS pool memory bytes |371230896|56422632|8304456|
| Physical models |36|36|127|
| RSS BVH nodes |-| -|28391|
| Tree ancestor links |0|0|264776|
| Cross-tree containment links |0|0|13024|

## Query Work

The original PQSS counter excludes each model-pair root RSS test. Both views are reported.

| Metric | Original truth | Official simplified | Candidate |
|---|---:|---:|---:|
| PQSS internal BV tests |340952|256818|592024|
| Root/scheduler RSS tests |741326|741326|16494180|
| All RSS distance tests |1082278|998144|17086204|
| Triangle tests |7880|3575|109792|
| Containment-cache skips |0|0|3722102|
| Positive results |1881|1929|11112|
| Parallel wall ms |97.079600|91.565100|2423.704600|

## Correctness Against Original Geometry

| Candidate | TP | TN | FP | FN | False-positive rate |
|---|---:|---:|---:|---:|---:|
| Official simplified |1881|739397|48|0|0.000065|
| Independent primitives |1881|730214|9231|0|0.012484|

A conservative candidate must have `FN = 0`. Detailed data are in `model_bvh_stats.csv` and `query_pair_stats.csv`.
