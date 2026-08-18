#pragma once

#include "pqss_proxy_mesh/topology_fill.hpp"

#include <cstddef>
#include <cstdint>
#include <filesystem>

namespace pqss_proxy_mesh
{

// Phase-2 surface analysis driven only by the frozen phase-1 geometry and
// halfedge opposite adjacency. Adjacent regions repeatedly compete for a
// conservative support-plane polygon merge until no certified merge remains.
struct HalfedgeSurfaceAnalysisOptions
{
    // One-sided proxy-to-phase1 distance limit in model units.
    double maximum_directed_hausdorff = 100.0;
    std::uint32_t maximum_certificate_depth = 18;
    double numerical_tolerance = 1.0e-9;
};

struct HalfedgeSurfaceAnalysisStats
{
    std::size_t phase1_vertices = 0;
    std::size_t phase1_triangles = 0;
    // Kept as an API name; this is the number of final support-plane or exact regions.
    std::size_t planar_regions = 0;
    std::size_t emitted_planar_polygons = 0;
    std::size_t fallback_triangles = 0;
    std::size_t accepted_region_merges = 0;
    std::size_t rejected_orientation_merges = 0;
    std::size_t rejected_triangle_count_merges = 0;
    std::size_t rejected_hausdorff_merges = 0;
    std::size_t final_vertices = 0;
    std::size_t final_triangles = 0;
    // Every final region is either an exact source triangle or a support-plane
    // convex hull whose projection contains all responsibility triangles.
    bool coverage_certified = false;
    bool global_hausdorff_certified = false;
    double maximum_directed_hausdorff = 0.0;
    double certified_directed_hausdorff_upper_bound = 0.0;
    std::size_t global_hausdorff_reference_queries = 0;
    double elapsed_seconds = 0.0;
};

// Simplifies a phase-1 halfedge mesh into a conservative surface proxy and
// writes proxy.obj (ordinary triangle OBJ) and model.json under output_directory.
[[nodiscard]] HalfedgeSurfaceAnalysisStats analyzeHalfedgeSurface(
    const OrientedSurfaceMesh& mesh,
    const std::filesystem::path& output_directory,
    const HalfedgeSurfaceAnalysisOptions& options);

} // namespace pqss_proxy_mesh
