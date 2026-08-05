# PQSSProxyMesh Real-Scene Final Report

> Historical full-query result. The canonical benchmark now uses only poses
> involving at least one officially simplified model. See
> `SCOPED_COMPARISON.md`; do not compare new candidates against the totals in
> this report.

## Test Contract

- CPU PQSS with `BuildStrategy::Optimized` and one query thread.
- All 882,008 poses from the PQSS real-scene `collision_pairs.txt` file.
- Original baseline: original `1.obj` through `23.obj`.
- Proxy candidate: conservative proxy `1.obj` through `23.obj`.
- Both pools append the original workpieces `24, 26-37` in the production order.
- Requested maximum proxy BVH depth: 8.
- False-positive rate: `FP / (FP + TN)`.

## Result

The depth-8 candidate passed both hard constraints:

- actual maximum PQSS BVH depth: 8;
- false negatives: 0.

| Metric | Original | Depth-8 proxy | Change |
|---|---:|---:|---:|
| Base input triangles | 627,190 | 247 | -99.9606% |
| Base built triangles | 629,236 | 277 | -99.9560% |
| Base BV nodes | 1,258,449 | 531 | -99.9578% |
| Pool memory | 371,230,896 B | 5,932,056 B | -98.4021% |
| Maximum base-model depth | 59 | 8 | satisfied |
| Query wall time | 242.6813 ms | 221.9395 ms | -8.5469% |
| BV tests | 2,594,432 | 2,075,940 | -19.9848% |
| Triangle tests | 140,359 | 182,707 | +30.1712% |
| Positive results | 16,557 | 49,285 | +32,728 |

Correctness matrix:

| TP | TN | FP | FN | False-positive rate |
|---:|---:|---:|---:|---:|
| 16,557 | 832,723 | 32,728 | 0 | 3.7816% |

The higher triangle-test count comes from conservative false positives. The
shallower, much smaller BVHs still reduce total query time.

## Construction Time

The original fresh-tag Optimized construction was allowed to finish on this
machine. There were no cache files from an earlier run under that tag; models
with identical content may still reuse a file created earlier in the same run.
Construction and query timings are intentionally separated.

| Metric | Original | Depth-8 proxy |
|---|---:|---:|
| Base OBJ loading | 937.6168 ms | 1.0561 ms |
| Base `EndPool()` | 1,182,713.7451 ms | 217.1157 ms |
| Base total | 1,183,651.3619 ms | 218.1719 ms |
| Total including appended workpieces | 1,190,610.1400 ms | 7,210.3835 ms |

The base construction speedup is about 5,425x. Build time is not part of the
query objective, but is reported as requested. Repeated validation loads the
same BVHs from tagged caches and is not used for this table.

## Pareto Alternative

The adaptive depth-10 candidate uses 965 input triangles. It measured 226.4541
ms in its fresh comparison run, had no false negatives, and reduced the
false-positive rate to 1.2411%. The depth-8 candidate is faster in the latest
cached validation run but has more false positives. Both are retained because
they represent different depth/accuracy choices.

## Artifacts

- `outputs/real_scene/conservative_outer_adaptive_depth8_strict/`: final proxy OBJ pool.
- `reports/real_scene/adaptive_depth8_validated/comparison_report.md`: complete aggregate report.
- `reports/real_scene/adaptive_depth8_validated/model_bvh_stats.csv`: per-model geometry and BVH data.
- `reports/real_scene/adaptive_depth8_validated/query_pair_stats.csv`: per-model-pair query and correctness data.
- `reports/real_scene/adaptive_v2_depth8_comparison/`: depth-10 Pareto candidate measured with a requested budget of 8 but rejected by actual-depth validation.
