#pragma once

#include "pqss_proxy_mesh/mesh_pool.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <vector>

namespace pqss_proxy_mesh
{

struct VoxelGrid
{
    std::array<std::uint32_t, 3> shape{};
    double pitch = 1.0;
    Position3 origin{};
    std::vector<std::uint8_t> occupancy;

    [[nodiscard]] std::size_t index(
        std::uint32_t x, std::uint32_t y, std::uint32_t z) const;
    [[nodiscard]] bool occupied(
        std::uint32_t x, std::uint32_t y, std::uint32_t z) const;
};

struct TriangleVoxelProvenance
{
    // surface_voxels[face] contains every grid cell rasterized by that source
    // triangle. A cell may occur in multiple face lists.
    std::vector<std::vector<std::uint32_t>> surface_voxels;
};

inline constexpr std::uint32_t invalid_surface_index =
    std::numeric_limits<std::uint32_t>::max();

struct SurfaceHalfedge
{
    std::uint32_t origin = invalid_surface_index;
    std::uint32_t face = invalid_surface_index;
    std::uint32_t next = invalid_surface_index;
    std::uint32_t opposite = invalid_surface_index;
};

struct VoxelBoundaryCrossing
{
    std::array<std::uint32_t, 3> inside{};
    std::array<std::uint32_t, 3> outside{};
};

// This is the phase-1 boundary representation used by later C++ stages.  The
// OBJ is only a visualization export; adjacency and solid-boundary provenance
// live here.
struct OrientedSurfaceMesh
{
    MeshModel geometry;
    std::vector<SurfaceHalfedge> halfedges;
    std::vector<std::uint32_t> vertex_halfedges;
    std::vector<std::uint32_t> face_halfedges;
    std::vector<VoxelBoundaryCrossing> boundary_crossings;
    std::vector<std::uint32_t> face_boundary_crossings;
    // Set for source-constrained reconstruction. Invalid denotes a newly
    // generated fill surface rather than a face supported by source geometry.
    std::vector<std::uint32_t> face_source_faces;
    int euler_characteristic = 0;
    double signed_volume = 0.0;
};

struct Phase1Solid
{
    VoxelGrid occupancy;
    OrientedSurfaceMesh boundary;
};

struct AnalysisHalfedgeStats
{
    std::size_t dropped_degenerate_triangles = 0;
    std::size_t paired_edges = 0;
    std::size_t boundary_halfedges = 0;
    std::size_t nonmanifold_edge_groups = 0;
    std::size_t inconsistent_orientation_edges = 0;
    std::size_t face_components = 0;
    std::size_t dropped_duplicate_triangles = 0;
};

// Builds an analysis halfedge mesh without imposing watertightness or spatial
// non-intersection. Exactly two oppositely directed uses form an opposite pair;
// boundary and conflicting/non-manifold edge uses remain unpaired.
[[nodiscard]] OrientedSurfaceMesh buildAnalysisHalfedgeMesh(
    const MeshModel& mesh,
    AnalysisHalfedgeStats* stats = nullptr);
void writeAnalysisHalfedgeMesh(
    const std::filesystem::path& path,
    const OrientedSurfaceMesh& mesh);

// A source-constrained volumetric complex is the reconstruction domain for
// phase 1.  Unlike the old dual grid, its geometry is cut by source triangles;
// the frozen voxels label cells but do not place boundary vertices.
struct SourceConstrainedCell
{
    std::array<std::uint32_t, 4> vertices{};
    std::array<std::uint32_t, 4> neighbors{
        invalid_surface_index, invalid_surface_index,
        invalid_surface_index, invalid_surface_index};
    // Source triangle constraining the face opposite each tetrahedron vertex.
    // Invalid means that the face is free to connect adjacent volume cells.
    std::array<std::uint32_t, 4> source_faces{
        invalid_surface_index, invalid_surface_index,
        invalid_surface_index, invalid_surface_index};
    bool occupied = false;
};

struct SourceConstrainedCellComplex
{
    std::vector<Position3> vertices;
    std::vector<SourceConstrainedCell> cells;
};

// Extracts the oriented interface between occupied and empty tetrahedra.  A
// production complex is generated from the original triangle soup; the public
// operation is kept independent so topology, provenance, and exact source
// planes can be regression-tested without a meshing dependency.
[[nodiscard]] OrientedSurfaceMesh extractCellComplexBoundary(
    const SourceConstrainedCellComplex& complex);
[[nodiscard]] MeshModel canonicalizeCoplanarTriangleSoup(
    const MeshModel& source,
    double relative_tolerance = 1.0e-9);

struct BettiNumbers
{
    int beta0 = 0;
    int beta1 = 0;
    int beta2 = 0;

    auto operator<=>(const BettiNumbers&) const = default;
};

struct TopologyFillStep
{
    int axis = 0;
    int radius = 0;
    std::size_t added_voxels = 0;
    BettiNumbers topology_before{};
    BettiNumbers topology_after{};
};

struct TopologyFillStats
{
    BettiNumbers input_betti{};
    BettiNumbers output_betti{};
    std::size_t source_triangles = 0;
    std::size_t input_voxels = 0;
    std::size_t output_voxels = 0;
    std::size_t added_voxels = 0;
    std::size_t cavity_fill_voxels = 0;
    std::size_t removed_voxels = 0;
    std::size_t exposed_voxel_faces = 0;
    std::size_t mesh_triangles = 0;
    std::size_t planar_regions = 0;
    std::size_t mesh_vertices = 0;
    std::size_t retained_source_triangles = 0;
    std::size_t removed_internal_triangles = 0;
    std::size_t ambiguous_source_triangles = 0;
    std::size_t hole_boundary_loops = 0;
    std::size_t cap_triangles = 0;
    std::size_t halfedge_count = 0;
    std::size_t paired_halfedge_edges = 0;
    std::size_t boundary_halfedges = 0;
    std::size_t nonmanifold_edge_groups = 0;
    std::size_t inconsistent_orientation_edges = 0;
    std::size_t dropped_degenerate_triangles = 0;
    std::size_t halfedge_face_components = 0;
    std::size_t dropped_duplicate_triangles = 0;
    bool mesh_watertight = false;
    bool mesh_oriented = false;
    bool mesh_manifold = false;
    bool mesh_connected = false;
    bool mesh_has_only_boundary_faces = false;
    int mesh_euler_characteristic = 0;
    double mesh_signed_volume = 0.0;
    double pitch = 0.0;
    std::array<std::uint32_t, 3> grid_shape{};
    double elapsed_seconds = 0.0;
    std::vector<TopologyFillStep> steps;
};

struct TopologyFillOptions
{
    std::size_t maximum_grid_voxels = 12'000'000;
    std::uint32_t padding = 4;
    std::uint32_t maximum_steps = 32;
};

[[nodiscard]] BettiNumbers voxelBettiNumbers(const VoxelGrid& grid);
[[nodiscard]] VoxelGrid enclosingTopologyFill(
    const VoxelGrid& source,
    std::uint32_t maximum_steps,
    TopologyFillStats* stats = nullptr,
    VoxelGrid* cavity_fill_labels = nullptr);
[[nodiscard]] VoxelGrid voxelizeTriangleSoup(
    const MeshModel& mesh,
    std::size_t maximum_grid_voxels,
    std::uint32_t padding,
    TriangleVoxelProvenance* provenance = nullptr);
// Rasterizes and flood-classifies the source on an existing grid. This does
// not run topology repair; it is used to compare a frozen filled solid with
// the original triangle soup on exactly the same cells.
[[nodiscard]] VoxelGrid voxelizeTriangleSoupOnGrid(
    const MeshModel& mesh,
    const VoxelGrid& grid_template,
    TriangleVoxelProvenance* provenance = nullptr);

struct OriginalObjHoleSurgeryStats
{
    std::size_t source_voxels = 0;
    std::size_t filled_voxels = 0;
    std::size_t added_voxels = 0;
    std::size_t source_triangles = 0;
    std::size_t retained_triangles = 0;
    std::size_t removed_internal_triangles = 0;
    std::size_t removed_closed_cavity_triangles = 0;
    std::size_t removed_closed_cavity_components = 0;
    std::size_t closed_cavity_voxels = 0;
    std::size_t ambiguous_triangles = 0;
    std::size_t interface_edges = 0;
    std::size_t closed_loops = 0;
    std::size_t cap_triangles = 0;
};

// Uses the frozen filled occupancy only as an inside/outside certificate.
// Existing boundary geometry comes directly from the input OBJ; only
// certified internal walls are removed and their source-edge loops are capped.
[[nodiscard]] MeshModel buildOriginalObjHoleSurgery(
    const VoxelGrid& filled,
    const VoxelGrid& cavity_fill_labels,
    const VoxelGrid& original,
    const TriangleVoxelProvenance& provenance,
    const MeshModel& source,
    OriginalObjHoleSurgeryStats* stats = nullptr);
[[nodiscard]] Phase1Solid buildPhase1Solid(
    VoxelGrid occupancy,
    const MeshModel& source);
[[nodiscard]] Phase1Solid buildPhase1Solid(VoxelGrid occupancy);
[[nodiscard]] OrientedSurfaceMesh projectPhase1BoundaryToSource(
    const OrientedSurfaceMesh& boundary,
    const VoxelGrid& filled_occupancy,
    const VoxelGrid& source_occupancy,
    const MeshModel& source);
// Visual inspection mesh that keeps certified source-boundary triangles exact
// and uses the frozen topology surface only where the source has no support.
[[nodiscard]] MeshModel buildSourceDominatedPhase1Preview(
    const VoxelGrid& occupancy,
    const MeshModel& source);
[[nodiscard]] MeshModel readTriangleSoupObj(const std::filesystem::path& path);
void writePhase1BoundaryObj(
    const std::filesystem::path& path,
    const OrientedSurfaceMesh& surface);
void writeVoxelGrid(const std::filesystem::path& path, const VoxelGrid& grid);
[[nodiscard]] VoxelGrid readVoxelGrid(const std::filesystem::path& path);
[[nodiscard]] TopologyFillStats generateTopologyFillModel(
    const std::filesystem::path& source_obj,
    const std::filesystem::path& output_directory,
    const TopologyFillOptions& options);
[[nodiscard]] std::vector<TopologyFillStats> generateTopologyFillBatch(
    const std::filesystem::path& source_directory,
    const std::vector<int>& model_ids,
    const std::filesystem::path& output_root,
    const TopologyFillOptions& options);

} // namespace pqss_proxy_mesh
