# Unified cavity filling for conservative collision proxies

## Problem

The gaps between the sleepers in model 2 and the planar holes in model 3 must
not be handled by model-specific rules. Both are empty-space accessibility
features: an outside probe of a chosen scale cannot enter them through their
mouths. A rail-side recess and model 3's main cavity are also the same class of
feature, but differ in geometric depth and therefore in outward-envelope error.

## Relevant work

Portaneri, Rouxel-Labbe, Hemmer, Cohen-Steiner, and Alliez, *Alpha Wrapping with
an Offset*, ACM Transactions on Graphics 41(4), 2022,
DOI `10.1145/3528223.3530152`, supplies the closest algorithmic foundation.
It wraps from outside to inside through a Delaunay complex. A gate is traversable
only when its empty-ball clearance is at least `alpha`. Consequently cavities
and holes with smaller mouths are inaccessible and remain inside the wrap. The
algorithm strictly encloses arbitrary triangle soups. Its independent `offset`
parameter places output vertices on an unsigned-distance level set and controls
tightness.

CGAL's 3D Alpha Wrapping manual documents the proof-level properties and the
outside flood-fill interpretation:
`https://doc.cgal.org/latest/Alpha_wrap_3/index.html`.

Jiang et al., *Feature-Preserving Shrink Wrapping with Adaptive Alpha*, Computer
Aided Geometric Design 113, 2024, DOI `10.1016/j.cagd.2024.102321`, extends the
idea to a spatially adaptive alpha field. This is relevant when one global
feature scale is insufficient, but adaptation must be driven by geometry rather
than model identity.

Classical binary morphological closing, `(S xor B_r) ominus B_r`, is the regular
set analogue: it fills gaps and concavities that a radius-`r` structuring ball
cannot enter. A voxel implementation is only an approximation unless its cell
error is included in the bound.

## Required unified decision

1. Construct an outside empty-space decomposition (Delaunay cells or an
   adaptive octree with certified triangle distance bounds).
2. Flood from infinity while recording each gate's clearance. This produces a
   maximin accessibility value for every empty-space cell: the radius of the
   largest probe that can reach it from outside.
3. For a scale `alpha`, cells with accessibility below `alpha` are sealed by one
   rule. This handles both model 3's planar holes and model 2's sleeper gaps.
4. Exterior-connected recesses that remain reachable are separate candidate
   basins. Generate a cap from their mouth and measure the one-sided outward
   deviation from the cap to the source surface. Fill only if the certified
   maximum is below `epsilon_fill` (optionally also constrain added area or
   volume). This fills shallow rail-side recesses while retaining model 3's deep
   main cavity.
5. Remove every source or proxy face that becomes internal after sealing. Emit
   only the boundary, then classify boundary patches as rectangle, frustum, or
   triangle.

`alpha` therefore controls accessibility scale; `epsilon_fill` controls maximum
outward-envelope error. Neither parameter identifies a model or a semantic part.

## Implemented first stage

The analyzer now constructs an axis-envelope candidate for every axis; there is
no aspect-ratio or model-ID path. Closed projection components and exterior
concavity components are evaluated independently by their conservative prism
volume before being filled. Each complete candidate is then measured against a
three-projection visual hull. A candidate is usable only below the total added
volume limit, and is selected over the structural result only when it reduces
the proxy primitive count by the configured factor. Remaining projected cavity
count is the primary tie-breaker between candidates.

This is a conservative raster approximation of the intended empty-space method,
not yet an exact Delaunay accessibility implementation. Its discretization error
is one projection cell and is exposed through the selected/candidate statistics.
