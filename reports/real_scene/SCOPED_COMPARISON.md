# Official-Simplification-Scope Comparison

## Canonical Scope

- Replacement models: `2, 3, 4, 5, 12, 13, 14, 15, 16, 17, 18, 19, 20`.
- Models `1, 6-11, 21-23` remain original in every candidate base pool.
- Appended workpieces `24, 26-37` remain original.
- A query is retained when at least one endpoint is a replacement model.
- Fixed query file: `test_data/real_scene/official_simplified_scope/collision_pairs.txt`.

The fixed file retains 741,326 of 882,008 queries and excludes 140,682. Its
SHA-256 is `cf21de26fb2a5020cc77d8f413e5c417328741f9937a5545df5dcee81a8429ae`.
Future real-scene comparisons use this file directly; the comparator does not
filter queries at runtime.

## Result

All results use PQSS CPU `Optimized`, one query thread, the fixed query file,
and original geometry for models outside the replacement set.

| Metric | Original | Official simplified | Current strict outer |
|---|---:|---:|---:|
| Replaced-model maximum depth | 59 | 19 | 7 |
| Whole base-pool maximum depth | 59 | 19 | 15 |
| BV tests | 340,952 | 256,818 | 931,862 |
| Triangle tests | 7,880 | 3,575 | 165,173 |
| Query wall time | 87.2138 ms | 84.2430 ms | 116.2334 ms |
| Positive results | 1,881 | 1,929 | 19,874 |
| False positives | 0 | 48 | 17,993 |
| False negatives | 0 | 0 | 0 |
| False-positive rate | 0 | 0.00649% | 2.4333% |

The whole-pool depth of the current hybrid is 15 because unchanged original
`6.obj` has depth 15. The simplification depth contract is evaluated only on
the 13 replacement models; their current maximum is 7 and satisfies depth 8.

The current strict outer result does not beat the official simplified data.
Its large boxes create broad-phase overlap and turn many originally rejected
pairs into primitive tests. The largest regressions are `6 x 18`, `15 x 5`,
`37 x 18`, `16 x 5`, `6 x 13`, and `19 x 17`.

The new acceptance targets are therefore:

1. certified strict containment and zero false negatives;
2. fewer than 256,818 BV tests;
3. fewer than 3,575 triangle tests;
4. replacement-model maximum PQSS depth at most the caller's optional cap.

BV and triangle work decide quality. Depth is only a reported value or a hard
caller constraint; a shallower BVH does not outrank a lower-work candidate.

The detailed reports are in `official_simplified_scoped/` and
`adaptive_depth8_official_scope/`.
