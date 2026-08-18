#pragma once

#include "pqss_proxy_mesh/hausdorff_simplifier.hpp"
#include "pqss_proxy_mesh/topology_fill.hpp"

#include <cstddef>
#include <filesystem>
#include <string>

namespace pqss_proxy_mesh
{

struct WorkingHalfedgeValidationReport
{
    std::size_t vertices = 0;
    std::size_t triangles = 0;
    std::size_t halfedges = 0;
    std::size_t paired_edges = 0;
    std::size_t boundary_halfedges = 0;
    std::size_t face_components = 0;
};

struct Phase2HalfedgeOptions
{
    double maximum_directed_hausdorff = 100.0;
    std::uint32_t maximum_certificate_depth = 18;
    double numerical_tolerance = 1.0e-9;
    double coplanar_relative_tolerance = 1.0e-9;
    std::string model_id;
};

struct Phase2HalfedgeStats
{
    std::size_t input_vertices = 0;
    std::size_t input_triangles = 0;
    std::size_t output_vertices = 0;
    std::size_t output_triangles = 0;
    std::size_t output_halfedges = 0;
    std::size_t paired_edges = 0;
    std::size_t boundary_halfedges = 0;
    std::size_t face_components = 0;
    std::size_t dropped_degenerate_triangles = 0;
    std::size_t dropped_duplicate_triangles = 0;
    std::size_t nonmanifold_edge_groups = 0;
    std::size_t inconsistent_orientation_edges = 0;
    std::size_t candidate_faces_considered = 0;
    std::size_t candidate_faces_accepted = 0;
    std::size_t candidate_faces_rejected_disconnected = 0;
    std::size_t source_faces_retriangulated = 0;
    std::size_t closed_boundary_loops = 0;
    std::size_t boundary_cap_faces = 0;
    std::size_t outward_vertices_moved = 0;
    std::size_t outward_support_planes = 0;
    bool source_coverage_certified = false;
    std::string selected_candidate;
    DirectedHausdorffCertificate directed_hausdorff;
    double elapsed_seconds = 0.0;
};

// Validates a working halfedge complex. Boundaries and spatial self-intersection
// are permitted; every stored opposite must still be the reverse use of the
// same geometric segment and every topology vertex must own one face fan.
[[nodiscard]] WorkingHalfedgeValidationReport validateWorkingHalfedgeTopology(
    const OrientedSurfaceMesh& mesh,
    double relative_tolerance = 1.0e-12);

// Phase 2: conservatively relax the frozen phase-1 surface outward onto simpler
// local support structure while preserving a valid halfedge topology. Phase 1 is
// the sole geometry input. Face-count minimization and coplanar union belong to
// later stages; phase 2 writes phase2_halfedge.bin as its formal artifact.
[[nodiscard]] Phase2HalfedgeStats generatePhase2Halfedge(
    const std::filesystem::path& phase1_halfedge,
    const std::filesystem::path& output_directory,
    const Phase2HalfedgeOptions& options = {});

} // namespace pqss_proxy_mesh
