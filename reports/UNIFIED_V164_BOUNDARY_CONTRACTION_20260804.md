# Unified v164 boundary-contraction report

## Result

The final proxy is clipped against the complete closed-shell certificates saved
before contained faces are removed. This matters because containment cleanup can
open a shell whose volume is still needed to remove overlap from a neighboring
proxy. One-sided owner sections retain true boundaries, owner crossing sections
remove surfaces strictly inside their own certified solid, and opposite-sided
owner/other sections contract shared interfaces.

No rule depends on a model ID, filename, coordinate, scene pose, or collision
query. All 13 models use one Release binary and identical options.

Manifest:

`outputs/surface_primitives_v164_union_boundary_contracted/viewer_manifest.json`

## Complete batch

| Model | Source tris | Primitives | Proxy tris | Time (s) | Peak MB | Added volume | Certificates | Clipped patches | Removed internal area |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 2 | 62,888 | 6 | 12 | 5.95 | 53.7 | 0.000% | 2 | 0 | 0.0 |
| 3 | 16,976 | 16 | 52 | 4.18 | 22.4 | 2.709% | 6 | 8 | 355,720.0 |
| 4 | 65,538 | 18 | 28 | 5.50 | 53.6 | 0.000% | 3 | 1 | 522,000.0 |
| 5 | 72,304 | 12 | 36 | 6.93 | 54.2 | 0.000% | 2 | 6 | 1,031,174.0 |
| 12 | 11,326 | 39 | 147 | 9.63 | 36.4 | 4.963% | 10 | 32 | 127,449.3 |
| 13 | 45,496 | 36 | 114 | 15.57 | 80.4 | 0.000% | 3 | 6 | 63,625.3 |
| 14 | 41,981 | 29 | 188 | 24.88 | 69.1 | 25.632% | 7 | 23 | 350,864.4 |
| 15 | 25,030 | 30 | 44 | 2.70 | 44.8 | 1.875% | 3 | 10 | 271,354.1 |
| 16 | 154,877 | 54 | 148 | 47.88 | 218.7 | 22.301% | 6 | 31 | 1,595,264.5 |
| 17 | 24,914 | 18 | 28 | 6.04 | 31.5 | 5.526% | 2 | 5 | 644,282.3 |
| 18 | 46,065 | 18 | 28 | 13.46 | 60.5 | 3.866% | 2 | 6 | 547,194.4 |
| 19 | 26,488 | 3 | 96 | 2.72 | 40.9 | 0.000% | 0 | 0 | 0.0 |
| 20 | 11,480 | 6 | 12 | 6.29 | 28.6 | 15.161% | 1 | 0 | 0.0 |

All models have `unassigned_source_faces=0` and `failed_source_faces=0`.
The largest single-process peak is 218.7 MB on model 16.

## Model 3 contraction

The previously uncontracted synthetic face spanned
`x=[2280,2780], z=[3485.75,3675.75]`. Complete certificates prove that its two
16-unit side strips and 20-unit lower strip lie inside adjacent envelopes. The
retained exposed boundary is `x=[2296,2764], z=[3505.75,3675.75]`.

Model 3 now recognizes six complete envelope certificates, clips eight boundary
patches, and removes 355,720 square units of internal overlap while preserving
zero failed source coverage.

## Verification

- C++ primitive analyzer regression: passed.
- Python suites: 25 passed, 3 subtests passed.
- Playwright viewer suite: 3 passed.
- Complete 13-model coverage audit: passed.
- Complete batch memory guard: passed; no process exceeded 218.7 MB.
