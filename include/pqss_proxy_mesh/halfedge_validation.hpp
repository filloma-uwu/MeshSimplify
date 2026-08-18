#pragma once

#include "pqss_proxy_mesh/topology_fill.hpp"

#include <cstddef>

namespace pqss_proxy_mesh
{

struct HalfedgeTopologyValidationReport
{
    std::size_t vertices = 0;
    std::size_t triangles = 0;
    std::size_t halfedges = 0;
    std::size_t boundary_halfedges = 0;
    std::size_t face_components = 0;
    int euler_characteristic = 0;
    double signed_volume = 0.0;
    // Geometric triangle-triangle intersections are a separate, more
    // expensive certificate and are deliberately not implied by this report.
    bool geometric_intersections_tested = false;
};

// Offline acceptance check for a phase-1 halfedge artifact. This performs a
// complete combinatorial traversal and throws when the mesh is not a closed,
// connected, consistently oriented genus-zero triangular 2-manifold. It does
// not reject coincident geometric positions, degenerate geometric triangles,
// or a negative signed volume after a surface projection.
// Generation and production loading do not call this function implicitly.
[[nodiscard]] HalfedgeTopologyValidationReport validatePhase1HalfedgeTopology(
    const OrientedSurfaceMesh& mesh);

// Validates one or more disjoint closed, consistently oriented genus-zero
// components. This is the phase-2 acceptance gate for multi-part proxies.
[[nodiscard]] HalfedgeTopologyValidationReport validateClosedHalfedgeTopology(
    const OrientedSurfaceMesh& mesh);

} // namespace pqss_proxy_mesh
