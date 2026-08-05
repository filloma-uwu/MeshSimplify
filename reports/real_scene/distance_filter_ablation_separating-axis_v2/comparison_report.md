# Independent Primitive RSS BVH Comparison

- Queries: 741326
- Threads: 1
- Truth: original geometry, original one-model/one-BVH query
- Official: official simplified geometry, original one-model/one-BVH query
- Candidate: independent primitive BLASes, cross-tree RSS containment pruning
- Candidate builder: project-local six-axis Optimized builder

- Distance filter: `separating-axis`

## Build

| Metric | Original | Official | Candidate |
|---|---:|---:|---:|
| Base build ms |2272.993900|328.594800|26.548000|
| Total build ms |2314.145200|365.409800|104.391400|
| Containment analysis ms |0|0|42.032900|
| PQSS pool memory bytes |371230896|56422632|8304456|
| Physical models |36|36|127|
| RSS BVH nodes |-| -|28391|
| Tree ancestor links |0|0|264776|
| Cross-tree containment links |0|0|13024|

## Query Work

The original PQSS counter excludes each model-pair root RSS test. Both views are reported.

| Metric | Original truth | Official simplified | Candidate |
|---|---:|---:|---:|
| PQSS internal BV tests |340952|256818|533825|
| Root/scheduler RSS tests |741326|741326|880973|
| All RSS distance tests |1082278|998144|1414798|
| Triangle tests |7880|3575|109792|
| Distance-filter skips |0|0|19393508|
| Containment checks |0|0|0|
| Containment skips |0|0|0|
| Separating-axis checks |0|0|21240366|
| Separating-axis skips |0|0|19393508|
| Separating-axis generation checks |0|0|1734898|
| Positive results |1881|1929|11112|
| Parallel wall ms |94.796000|101.549400|1439.017600|

## Correctness Against Original Geometry

| Candidate | TP | TN | FP | FN | False-positive rate |
|---|---:|---:|---:|---:|---:|
| Official simplified |1881|739397|48|0|0.000065|
| Independent primitives |1881|730214|9231|0|0.012484|

A conservative candidate must have `FN = 0`. Detailed data are in `model_bvh_stats.csv` and `query_pair_stats.csv`.
