#pragma once

#include "pqss_proxy_mesh/mesh_pool.hpp"
#include "pqss_proxy_mesh/topology_fill.hpp"

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

namespace pqss_proxy_mesh
{

struct DirectedHausdorffCertificate
{
    bool passed = false;
    double lower_bound = 0.0;
    double upper_bound = 0.0;
    std::size_t reference_queries = 0;
    std::size_t subdivision_nodes = 0;
    std::array<double, 3> maximum_proxy_point{};
    std::array<double, 3> nearest_reference_point{};
    std::size_t maximum_proxy_face = 0;
    std::string failure_reason;
};

struct HausdorffSimplificationOptions
{
    // Metadata label only. Geometry generation never reads or branches on it.
    std::string model_id;
    double maximum_directed_hausdorff = 100.0;
    std::uint32_t maximum_certificate_depth = 18;
    double numerical_tolerance = 1.0e-9;
    bool enable_axis_aligned_box = true;
    bool enable_convex_hull = true;
    bool enable_adaptive_convex_cover = true;
    bool enable_bounded_box_cover = true;
    bool enable_discrete_orientation_polytopes = true;
    bool enable_exact_coplanar_union = true;
    bool enable_planar_component_convexification = true;
};

struct HausdorffCandidateStats
{
    std::string name;
    std::size_t triangles = 0;
    bool conservative_coverage = false;
    DirectedHausdorffCertificate hausdorff;
    bool selected = false;
};

struct HausdorffSimplificationStats
{
    std::size_t phase1_vertices = 0;
    std::size_t phase1_triangles = 0;
    std::size_t final_vertices = 0;
    std::size_t final_triangles = 0;
    double maximum_directed_hausdorff = 0.0;
    double elapsed_seconds = 0.0;
    std::string selected_candidate;
    std::vector<HausdorffCandidateStats> candidates;
};

// Production reader: checks the binary framing and payload length only.
// Offline topology acceptance is an explicit halfedge_validation operation.
[[nodiscard]] OrientedSurfaceMesh readAnalysisHalfedgeMesh(
    const std::filesystem::path& path);

[[nodiscard]] DirectedHausdorffCertificate certifyDirectedHausdorff(
    const MeshModel& proxy,
    const MeshModel& reference,
    double maximum_distance,
    std::uint32_t maximum_depth = 18,
    double numerical_tolerance = 1.0e-9);

[[nodiscard]] MeshModel convexifyCertifiedPlanarComponents(
    const MeshModel& proxy,
    const MeshModel& reference,
    double maximum_distance,
    std::uint32_t maximum_depth = 18,
    double numerical_tolerance = 1.0e-9);

// Reuses one immutable reference-mesh BVH across many candidate certificates.
// The reference mesh must outlive the certifier.
class DirectedHausdorffCertifier
{
public:
    explicit DirectedHausdorffCertifier(const MeshModel& reference);
    ~DirectedHausdorffCertifier();
    DirectedHausdorffCertifier(DirectedHausdorffCertifier&&) noexcept;
    DirectedHausdorffCertifier& operator=(DirectedHausdorffCertifier&&) noexcept;
    DirectedHausdorffCertifier(const DirectedHausdorffCertifier&) = delete;
    DirectedHausdorffCertifier& operator=(const DirectedHausdorffCertifier&) = delete;

    [[nodiscard]] DirectedHausdorffCertificate certify(
        const MeshModel& proxy,
        double maximum_distance,
        std::uint32_t maximum_depth = 18,
        double numerical_tolerance = 1.0e-9) const;

private:
    struct Implementation;
    std::unique_ptr<Implementation> implementation_;
};

[[nodiscard]] HausdorffSimplificationStats simplifyPhase1Halfedge(
    const std::filesystem::path& phase1_halfedge,
    const std::filesystem::path& source_obj,
    const std::filesystem::path& output_directory,
    const HausdorffSimplificationOptions& options = {});

} // namespace pqss_proxy_mesh
