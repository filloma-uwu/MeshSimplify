# Algorithm Direction

The project is a constrained minimum-complexity outer approximation problem.

## Feasible Set

For source mesh \(M\), threshold \(\varepsilon\), and candidate proxy \(P\), a
candidate is feasible only if:

- it conservatively covers the source according to the project's geometric
  coverage certificate; and
- its directed proxy-to-source error is at most \(\varepsilon\).

## Ranking

For feasible candidates:

1. fewer final OBJ triangles is better;
2. ties may use deterministic geometric tie-breakers;
3. PQSS/BVH/query metrics are never tie-breakers.

## Current Research Focus

Stage 2 and Stage 3 should evolve toward a model-level candidate search rather
than a collection of unrelated local heuristics. Useful candidate generators
may include planar union/replacement, analytic surface recognition, conservative
hole/cavity treatment, and enclosing surface constructions, but every candidate
must enter the same feasibility-and-face-count comparison.

## Certification Gap

The current code contains deterministic surface sampling for directed
proxy-to-source distance. Sampling is useful validation, but a sampled maximum
does not prove the continuous one-sided Hausdorff bound.

A future proof-level certificate can exploit the 1-Lipschitz property of
point-to-mesh distance with adaptive subdivision/bounds, exact distance extrema
for restricted primitive families, or another conservative upper-bounding
scheme. Until then documentation and output metadata must distinguish
"measured" from "certified".
