#pragma once

#include "pqss_proxy_mesh/mesh_pool.hpp"

#include <cstddef>
#include <filesystem>
#include <vector>

namespace pqss_proxy_mesh
{

struct Phase0OpenHoleOptions
{
    // Numerical tolerance used by the 2D parameter-domain triangulation. This
    // is not a simplification error and does not require the 3D loop to be planar.
    double projection_relative_tolerance = 1.0e-12;
};

struct Phase0OpenHoleCap
{
    std::size_t boundary_vertices = 0;
    std::size_t first_cap_triangle = 0;
    std::size_t cap_triangles = 0;
    double area = 0.0;
    double maximum_planarity_error = 0.0;
};

struct Phase0OpenHoleStats
{
    std::size_t source_vertices = 0;
    std::size_t source_triangles = 0;
    std::size_t analysis_vertices = 0;
    std::size_t analysis_triangles = 0;
    std::size_t reality_boundary_edges = 0;
    std::size_t closed_boundary_loops = 0;
    std::size_t rejected_branched_boundaries = 0;
    std::size_t nonplanar_loops = 0;
    std::size_t rejected_nonsimple_loops = 0;
    std::size_t rejected_triangulation_failures = 0;
    std::size_t capped_loops = 0;
    std::size_t cap_triangles = 0;
    double elapsed_seconds = 0.0;
};

struct Phase0OpenHoleResult
{
    // The original source triangles are retained byte-for-byte at the front.
    MeshModel covered_mesh;
    MeshModel caps;
    std::vector<Phase0OpenHoleCap> cap_records;
    Phase0OpenHoleStats stats;
};

// Detects physical, single-use boundary segments after geometric welding and
// T-junction conformance. A cap is emitted for every closed, unbranched reality
// boundary loop whose 3D cycle has a simple 2D parameterization. The loop is
// not required to be planar. Occupancy grids and model-specific rules are not used.
[[nodiscard]] Phase0OpenHoleResult coverPlanarOpenHoles(
    const MeshModel& source,
    const Phase0OpenHoleOptions& options = {});

// Generates an independent phase-0 inspection result. The combined OBJ keeps
// the original source faces and appends caps; phase0_caps.obj contains only the
// synthetic faces so the viewer can overlay them distinctly.
[[nodiscard]] Phase0OpenHoleStats generatePhase0OpenHoleModel(
    const std::filesystem::path& source_obj,
    const std::filesystem::path& output_directory,
    const Phase0OpenHoleOptions& options = {});

} // namespace pqss_proxy_mesh
