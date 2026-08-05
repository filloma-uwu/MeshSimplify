# Primitive-analysis v2 BVH Repeat Summary

- Queries per run: 741326
- Threads: 1
- Candidate analysis: 13 models, 605628 active source triangles, 797 analytic primitives
- Candidate cache: 13/13 hits in every timed run
- Baselines: standard PQSS Optimized builder and query

| Run | Original ms | Official simplified ms | Primitive BVH v2 ms |
|---:|---:|---:|---:|
|1|96.2784|91.2754|87.3949|
|2|98.8795|96.7448|88.9805|
|3|93.8354|94.7029|88.6101|
|4|96.0272|93.2234|85.6420|
|5|96.1897|90.0074|88.9174|
|Median|96.1897|93.2234|88.6101|

The v2 primitive BVH median is 7.88% lower than original geometry and 4.95%
lower than the official simplification.

## Query Work

The root test is unavoidable and contributes one test for every query. Internal
BV tests therefore show the change attributable to BVH traversal more directly.

| Metric | Original | Official simplified | Primitive BVH v2 |
|---|---:|---:|---:|
|Root BV tests|741326|741326|741326|
|Internal BV tests|340952|256818|265614|
|Total BV tests|1082278|998144|1006940|
|Triangle/analytic-leaf tests|7880|3575|1975|
|False positives|0|48|94|
|False negatives|0|0|0|

Relative to original geometry, v2 reduces internal BV traversal by 22.10%,
total BV tests by 6.96%, and terminal tests by 74.94%. Relative to the official
simplification, it performs 3.43% more internal BV tests but 44.76% fewer
terminal tests. The more expensive terminal work reduction is sufficient to
produce the lower measured query time.

The cold candidate BVH build took 259614 ms and produced 13 persistent cache
files. Build time is reported only and is excluded from the query objective.
