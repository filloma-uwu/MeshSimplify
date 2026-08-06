# Algorithm Invariants

The repository-level rules in `AGENTS.md` are the acceptance contract for the
proxy generator. In particular, a visually observed defect is evidence for a
missing general geometric operation, not authorization for a model-specific
path.

For planar output, the canonical operation is:

```text
geometrically coplanar proxy faces
    -> project into one 2D plane
    -> compute the exact covered-region union
    -> remove internal and duplicate edges
    -> emit each union component as one semantic Polygon
    -> recognize analytic surface families where certified
    -> triangulate surface primitives in a separate collision stage
```

Semantic primitive analysis emits `Polygon`, `Disk`, `Annulus`,
`CylindricalBand`, and `ConicalBand`. These are surfaces, not closed physics
colliders. A cylinder is two disks and one cylindrical band; an annular flange
uses annuli plus inner and outer bands. `primitives.obj` preserves this semantic
grouping for inspection. `proxy.obj` is the separate PQSS input: a simple
n-vertex polygon becomes exactly n-2 triangles without Steiner vertices, while
analytic circular surfaces use the configured circumferential segment count.

Hole decisions happen on the complete fitted surface domain, before fallback
fragmentation. Each inner loop is independently evaluated by a conservative
swept-volume upper bound divided by the whole-model scale volume. "Global"
refers to this whole-model denominator; holes do not consume a shared budget.
An unfilled concentric circular loop produces an `Annulus`; it must not force
the complete assembly into one closed frustum.

Topology in the OBJ is only evidence. Geometric coincidence and intersection
must still be handled when faces do not share vertex indices. A redundant proxy
is removed only when the retained proxy union covers its complete geometric
support; vertex-only or sparse-sample checks are not sufficient.

The collision OBJ contains only the exposed boundary of an accepted composite
proxy.  When an outward protrusion is replaced by five box-shaped surface
patches, all four side patches are clipped at the exact support plane and the
covered contact footprint is subtracted from the support polygon.  Omitting the
box back face alone is insufficient: leaving the support face underneath would
create an internal overlapping surface and enlarge RSS nodes.  Contact cutouts
are topologically protected from the later Hausdorff merge, because filling
them has zero directed surface error but would recreate query-only internal
geometry.  A support polygon with several contacts is partitioned into local
strips rather than triangulated across all holes, avoiding long triangles whose
RSS bounds span distant attachments.

Two tempting shortcuts are explicitly invalid:

- projected overlap between parallel planes is not three-dimensional coverage;
  a lower planar proxy cannot delete a raised patch merely because their 2D
  projections overlap;
- a small globally normalized area increase does not authorize replacing an
  arbitrary planar silhouette with its bounding rectangle. Structural polygon
  boundaries, including concave or chamfered side faces, remain exact unless a
  complete conservative 3D envelope is emitted and certified.

A local box simplification therefore operates on a complete three-dimensional
responsibility region. It fits the full enclosing box, evaluates added visual-
hull volume against the whole-model volume and added shell area against the
whole-model surface, then emits all six Polygon faces. Shallow exterior feet are
part of that region; they must never be discarded as inward blind-cavity walls.

Parallel-shell processing has two ordered phases. First, genuinely adjacent
thin layers are collapsed only when bridge faces certify their connection; this
phase cannot cross the object's local depth interval. Exact coplanar union then
creates stable planar regions. Second, deeper terraced recesses are considered
using a layered projection integral: each depth band contributes only the union
area actually present in that band, divided by whole-model scale volume. Two
complete opposing exterior faces are never collapsed by this cavity phase.

Minor silhouette details may be removed only by a candidate that reduces the
later triangle count, exactly contains the current polygon, and passes the
whole-model added-area budget. A short bevel is replaced by extending its two
adjacent main boundary lines to their intersection. The algorithm must not
invent an axis-aligned bounding-box corner when retaining a slanted edge has
equal or lower collision complexity.

When two opposing conservative polygons have been synchronized to the same
silhouette, stale bridge faces are not automatically retained.  Their original
triangle responsibility is projected into the cap plane and certified inside
the polygonal extrusion and between its two cap planes.  Only then may the
bridge shell be rebuilt from one quad per cap edge.  Acceptance requires a
strict reduction of the implied swept cross-section and permits at most two
additional proxy triangles.  Thus an inward-sloping side is introduced because
it removes globally unnecessary volume at low collision cost, not merely to
make the exported surface watertight.

Opposing-silhouette synchronization and extrusion-side rebuilding are event
driven.  Their cap candidates are only polygons moved by the immediately
preceding parallel-shell phases (`shallow_shell_coalesced`); unchanged polygons
cannot have acquired stale bridge geometry and must not enter the cubic
all-polygon search.  The geometric certificate remains unchanged after this
broad-phase restriction.

The planar Boolean operation uses the vendored Clipper2 C++ implementation
under the Boost Software License 1.0. This dependency permits commercial and
closed-source use and imposes no watertightness or manifoldness requirement on
the input triangle soup.

Planar-detail absorption is transactional: a filled-notch or tiny-detail
candidate replaces the current boundary only when it remains a positive-area
simple polygon and still triangulates to exactly n-2 faces. Decimal-coordinate
triples whose altitude is below 1e-12 of their longest edge are treated as
numerical zero-area input slivers and removed before region analysis; zero-area
slivers produced by planar Boolean triangulation are likewise not emitted as
semantic primitives.

Local coplanar rectangle merging uses a stable-index candidate priority queue.
The initial pair set is evaluated once; after a merge, only candidates touching
the updated rectangle are recomputed.  This preserves the same pairwise
coplanarity, rectangular-hull, and whole-model added-area tests while avoiding
the former full all-pairs rescan after every accepted merge.
