# Three-way Repeat Summary

- Queries per run: 741326
- Threads: 1
- Builder/query for original and official simplified: standard PQSS Optimized `ModelPool`
- Primitive candidate query: the same PQSS `ModelPool::Query` after replacing 13 model BVHs

| Run | Original ms | Official simplified ms | Primitive BVH ms |
|---:|---:|---:|---:|
|1|95.8692|92.8363|90.8889|
|2|96.1568|94.6803|89.6101|
|3|96.1405|93.2673|89.7740|
|4|95.5614|94.4445|88.4961|
|5|95.6729|98.8975|89.4204|
|Median|95.8692|94.4445|89.6101|

The primitive BVH median is 5.12% lower than the official simplified median
and 6.53% lower than the original-geometry median.

| Metric | Original | Official simplified | Primitive BVH |
|---|---:|---:|---:|
|Total BV tests|1082278|998144|1016372|
|Triangle/analytic-leaf tests|7880|3575|1975|
|False positives|0|48|94|
|False negatives|0|0|0|
