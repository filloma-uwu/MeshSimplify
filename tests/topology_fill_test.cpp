#include "pqss_proxy_mesh/topology_fill.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>

namespace
{

void require(const bool condition, const std::string& message)
{
    if (!condition) throw std::runtime_error(message);
}

pqss_proxy_mesh::VoxelGrid makeGrid(const std::uint32_t size)
{
    pqss_proxy_mesh::VoxelGrid grid;
    grid.shape = {size, size, size};
    grid.occupancy.assign(static_cast<std::size_t>(size) * size * size, 0);
    return grid;
}

void testSingleCube()
{
    auto grid = makeGrid(5);
    grid.occupancy[grid.index(2, 2, 2)] = 1;
    require(pqss_proxy_mesh::voxelBettiNumbers(grid) ==
                pqss_proxy_mesh::BettiNumbers{1, 0, 0},
            "one voxel must be one contractible component");
}

void testCavity()
{
    auto grid = makeGrid(9);
    for (std::uint32_t x = 2; x <= 6; ++x)
        for (std::uint32_t y = 2; y <= 6; ++y)
            for (std::uint32_t z = 2; z <= 6; ++z)
                if (x == 2 || x == 6 || y == 2 || y == 6 || z == 2 || z == 6)
                    grid.occupancy[grid.index(x, y, z)] = 1;
    require(pqss_proxy_mesh::voxelBettiNumbers(grid) ==
                pqss_proxy_mesh::BettiNumbers{1, 0, 1},
            "closed voxel shell must contain one cavity");
}

void testDigitalConnectivityPair()
{
    auto grid = makeGrid(7);
    grid.occupancy[grid.index(2, 2, 2)] = 1;
    grid.occupancy[grid.index(3, 3, 2)] = 1;
    require(pqss_proxy_mesh::voxelBettiNumbers(grid) ==
                pqss_proxy_mesh::BettiNumbers{2, 0, 0},
            "edge-touching voxels must remain separate under foreground 6-connectivity");
    grid.occupancy.assign(grid.occupancy.size(), 0);
    grid.occupancy[grid.index(2, 2, 2)] = 1;
    grid.occupancy[grid.index(3, 3, 3)] = 1;
    require(pqss_proxy_mesh::voxelBettiNumbers(grid) ==
                pqss_proxy_mesh::BettiNumbers{2, 0, 0},
            "vertex-touching voxels must remain separate under foreground 6-connectivity");
}

void testTorusFill()
{
    constexpr std::uint32_t size = 32;
    auto grid = makeGrid(size);
    const double center = 0.5 * (size - 1);
    const double major = size * 0.24;
    const double minor = size * 0.10;
    for (std::uint32_t x = 0; x < size; ++x)
        for (std::uint32_t y = 0; y < size; ++y)
            for (std::uint32_t z = 0; z < size; ++z)
            {
                const double dx = x - center;
                const double dy = y - center;
                const double dz = z - center;
                const double radial = std::sqrt(dx * dx + dy * dy) - major;
                if (radial * radial + dz * dz <= minor * minor)
                    grid.occupancy[grid.index(x, y, z)] = 1;
            }
    require(pqss_proxy_mesh::voxelBettiNumbers(grid) ==
                pqss_proxy_mesh::BettiNumbers{1, 1, 0},
            "voxel torus must contain one handle");
    pqss_proxy_mesh::TopologyFillStats stats;
    const auto filled = pqss_proxy_mesh::enclosingTopologyFill(grid, 32, &stats);
    require(stats.output_betti == pqss_proxy_mesh::BettiNumbers{1, 0, 0},
            "enclosing fill must remove the torus handle");
    require(stats.removed_voxels == 0 && stats.added_voxels > 0,
            "enclosing fill must be a strict superset");
    for (std::size_t index = 0; index < grid.occupancy.size(); ++index)
        require(!grid.occupancy[index] || filled.occupancy[index],
                "enclosing fill lost an occupied voxel");
}

void testBeta2IsFilledOnlyAfterHandleClosing()
{
    auto grid = makeGrid(13);
    // A hollow box whose cavity is connected to the exterior by a one-voxel
    // axial mouth. It begins with beta2=0; closing the mouth creates beta2,
    // and the final single cavity-fill pass must record that entire region.
    for (std::uint32_t x = 3; x <= 9; ++x)
        for (std::uint32_t y = 3; y <= 9; ++y)
            for (std::uint32_t z = 3; z <= 9; ++z)
                if (x == 3 || x == 9 || y == 3 || y == 9 || z == 3 || z == 9)
                    grid.occupancy[grid.index(x, y, z)] = 1;
    for (std::uint32_t x = 3; x <= 9; ++x)
        grid.occupancy[grid.index(x, 6, 6)] = 0;
    require(pqss_proxy_mesh::voxelBettiNumbers(grid).beta2 == 0,
            "open cavity must not be classified as beta2 before its mouth closes");
    pqss_proxy_mesh::TopologyFillStats stats;
    pqss_proxy_mesh::VoxelGrid labels;
    const auto filled = pqss_proxy_mesh::enclosingTopologyFill(
        grid, 32, &stats, &labels);
    require(stats.output_betti == pqss_proxy_mesh::BettiNumbers{1, 0, 0},
            "two-stage fill must end at beta=(1,0,0)");
    require(stats.cavity_fill_voxels == static_cast<std::size_t>(std::count(
                    labels.occupancy.begin(), labels.occupancy.end(), std::uint8_t{1})),
            "final beta2 fill must retain its exact post-responsibility difference");
    for (std::size_t index = 0; index < labels.occupancy.size(); ++index)
        require(!labels.occupancy[index] ||
                    (filled.occupancy[index] && !grid.occupancy[index]),
                "beta2 label must be a newly added final occupied voxel");
}

void testHandleRepairDoesNotFillUnrelatedExteriorNotch()
{
    auto grid = makeGrid(25);
    // A solid block with a through tunnel contributes beta1=1.  A separate
    // notch cut down from the exterior contributes no Betti defect and must
    // not be retained merely because axial closing spans both gaps.
    for (std::uint32_t x = 4; x <= 20; ++x)
        for (std::uint32_t y = 4; y <= 20; ++y)
            for (std::uint32_t z = 4; z <= 20; ++z)
                grid.occupancy[grid.index(x, y, z)] = 1;
    for (std::uint32_t x = 4; x <= 20; ++x)
        for (std::uint32_t y = 9; y <= 11; ++y)
            for (std::uint32_t z = 9; z <= 11; ++z)
                grid.occupancy[grid.index(x, y, z)] = 0;
    for (std::uint32_t x = 8; x <= 10; ++x)
        for (std::uint32_t y = 4; y <= 8; ++y)
            for (std::uint32_t z = 17; z <= 20; ++z)
                grid.occupancy[grid.index(x, y, z)] = 0;
    require(pqss_proxy_mesh::voxelBettiNumbers(grid) ==
                pqss_proxy_mesh::BettiNumbers{1, 1, 0},
            "through tunnel plus exterior notch must have exactly one handle");
    const auto filled = pqss_proxy_mesh::enclosingTopologyFill(grid, 32);
    require(pqss_proxy_mesh::voxelBettiNumbers(filled) ==
                pqss_proxy_mesh::BettiNumbers{1, 0, 0},
            "handle repair must remove the through tunnel");
    require(!filled.occupied(9, 5, 19),
            "handle repair retained a beta1-irrelevant exterior notch fill");
}

void testIncrementalClosingProfileConnectsSeparateComponents()
{
    auto grid = makeGrid(17);
    for (std::uint32_t x = 3; x <= 5; ++x)
        for (std::uint32_t y = 6; y <= 10; ++y)
            for (std::uint32_t z = 6; z <= 10; ++z)
                grid.occupancy[grid.index(x, y, z)] = 1;
    for (std::uint32_t x = 11; x <= 13; ++x)
        for (std::uint32_t y = 6; y <= 10; ++y)
            for (std::uint32_t z = 6; z <= 10; ++z)
                grid.occupancy[grid.index(x, y, z)] = 1;
    require(pqss_proxy_mesh::voxelBettiNumbers(grid) ==
                pqss_proxy_mesh::BettiNumbers{2, 0, 0},
            "disconnected closing fixture must begin with two components");
    pqss_proxy_mesh::TopologyFillStats stats;
    const auto filled = pqss_proxy_mesh::enclosingTopologyFill(grid, 8, &stats);
    require(stats.output_betti == pqss_proxy_mesh::BettiNumbers{1, 0, 0},
            "incremental closing must connect separate foreground components");
    for (std::uint32_t x = 6; x <= 10; ++x)
        require(filled.occupied(x, 8, 8),
                "incremental closing did not retain the required bridge");
    for (std::size_t index = 0; index < grid.occupancy.size(); ++index)
        require(!grid.occupancy[index] || filled.occupancy[index],
                "incremental closing bridge lost a source voxel");
}

void testAnalysisHalfedgeSplitsConflictingFans()
{
    pqss_proxy_mesh::MeshModel mesh;
    mesh.name = "halfedge_conflicts";
    mesh.vertices = {
        {0,0,0}, {1,0,0}, {0,1,0}, {1,1,0},
        {0,0,0}, {0,-1,0}, {0,0,1}};
    mesh.triangles = {
        {{0,1,2}}, {{1,0,3}},
        {{4,1,5}}, {{0,1,6}}};
    pqss_proxy_mesh::AnalysisHalfedgeStats stats;
    const auto result = pqss_proxy_mesh::buildAnalysisHalfedgeMesh(mesh, &stats);
    require(result.halfedges.size() == result.geometry.triangles.size() * 3,
            "analysis halfedge must create three halfedges per face");
    require(stats.nonmanifold_edge_groups == 1,
            "three or more uses of a geometric edge must remain unpaired");
    for (std::uint32_t face = 0; face < result.geometry.triangles.size(); ++face)
    {
        const std::uint32_t first = result.face_halfedges[face];
        const std::uint32_t second = result.halfedges[first].next;
        const std::uint32_t third = result.halfedges[second].next;
        require(result.halfedges[third].next == first,
                "analysis halfedge face must be a closed three-cycle");
    }
    for (std::uint32_t edge = 0; edge < result.halfedges.size(); ++edge)
        if (result.halfedges[edge].opposite !=
            pqss_proxy_mesh::invalid_surface_index)
            require(result.halfedges[result.halfedges[edge].opposite].opposite == edge,
                    "analysis opposite relation must be symmetric");
}

void testPhase1SolidBoundary()
{
    auto grid = makeGrid(7);
    grid.occupancy[grid.index(3, 3, 3)] = 1;
    const auto solid = pqss_proxy_mesh::buildPhase1Solid(grid);
    const auto& boundary = solid.boundary;
    require(!boundary.geometry.triangles.empty(), "phase-1 boundary must not be empty");
    require(boundary.halfedges.size() == boundary.geometry.triangles.size() * 3,
            "each boundary triangle must own three halfedges");
    require(boundary.face_boundary_crossings.size() == boundary.geometry.triangles.size(),
            "each face must retain its voxel-boundary certificate");
    require(boundary.euler_characteristic == 2 && boundary.signed_volume > 0.0,
            "single-solid boundary must be an outward genus-zero surface");
    for (std::uint32_t index = 0; index < boundary.halfedges.size(); ++index)
    {
        const auto& halfedge = boundary.halfedges[index];
        require(halfedge.opposite < boundary.halfedges.size(),
                "closed boundary halfedge must have an opposite");
        require(boundary.halfedges[halfedge.opposite].opposite == index,
                "halfedge opposite relation must be symmetric");
        require(halfedge.next < boundary.halfedges.size() &&
                    boundary.halfedges[halfedge.next].face == halfedge.face,
                "halfedge next relation must stay within its face");
    }
}

void testDualVertexDisconnectedFansAreSplit()
{
    auto grid = makeGrid(9);
    std::uint32_t state = 144;
    std::array<std::uint32_t, 3> cursor{4,4,4};
    for (int step = 0; step < 32; ++step)
    {
        grid.occupancy[grid.index(cursor[0], cursor[1], cursor[2])] = 1;
        state = state * 1664525u + 1013904223u;
        const int axis = static_cast<int>((state >> 16) % 3);
        const int direction = (state & 1u) ? 1 : -1;
        cursor[axis] = static_cast<std::uint32_t>(std::clamp(
            static_cast<int>(cursor[axis]) + direction, 2, 6));
    }
    require(pqss_proxy_mesh::voxelBettiNumbers(grid) ==
                pqss_proxy_mesh::BettiNumbers{1, 0, 0},
            "fan-split fixture must be one contractible occupied set");
    const auto solid = pqss_proxy_mesh::buildPhase1Solid(grid);
    std::set<std::array<double, 3>> positions;
    bool split_position = false;
    for (const auto point : solid.boundary.geometry.vertices)
        if (!positions.insert({point.x, point.y, point.z}).second)
            split_position = true;
    require(split_position,
            "fan-split fixture did not exercise duplicated topology vertices");
    require(solid.boundary.euler_characteristic == 2 &&
                solid.boundary.signed_volume > 0.0,
            "fan splitting must preserve the closed genus-zero boundary");
    for (std::uint32_t edge = 0; edge < solid.boundary.halfedges.size(); ++edge)
    {
        const auto opposite = solid.boundary.halfedges[edge].opposite;
        require(opposite < solid.boundary.halfedges.size() &&
                    solid.boundary.halfedges[opposite].opposite == edge,
                "fan-split boundary must retain symmetric opposite halfedges");
    }
}

void testDualCheckerboardEdgeIsSplitByOccupiedSector()
{
    auto grid = makeGrid(9);
    const std::array<std::array<std::uint32_t, 3>, 5> path{{
        {{4,3,3}}, {{3,3,3}}, {{3,4,3}}, {{3,4,4}}, {{4,4,4}}}};
    for (const auto& cell : path)
        grid.occupancy[grid.index(cell[0], cell[1], cell[2])] = 1;
    require(pqss_proxy_mesh::voxelBettiNumbers(grid) ==
                pqss_proxy_mesh::BettiNumbers{1, 0, 0},
            "checkerboard dual-edge fixture must be contractible");
    const auto solid = pqss_proxy_mesh::buildPhase1Solid(grid);
    require(solid.boundary.euler_characteristic == 2 &&
                solid.boundary.signed_volume > 0.0,
            "checkerboard dual edge must extract as a closed genus-zero surface");
    for (std::uint32_t edge = 0; edge < solid.boundary.halfedges.size(); ++edge)
    {
        const auto opposite = solid.boundary.halfedges[edge].opposite;
        require(opposite < solid.boundary.halfedges.size() &&
                    solid.boundary.halfedges[opposite].opposite == edge,
                "checkerboard dual edge left an unpaired halfedge");
    }
}

void testVoxelGridRoundTrip()
{
    auto grid = makeGrid(7);
    grid.pitch = 0.125;
    grid.origin = {-1.0, 2.0, 3.5};
    grid.occupancy[grid.index(2, 3, 4)] = 1;
    grid.occupancy[grid.index(5, 1, 2)] = 1;
    const auto path = std::filesystem::temp_directory_path() /
        "pqss_topology_fill_round_trip.vox";
    std::error_code ignored;
    std::filesystem::remove(path, ignored);
    pqss_proxy_mesh::writeVoxelGrid(path, grid);
    const auto loaded = pqss_proxy_mesh::readVoxelGrid(path);
    std::filesystem::remove(path, ignored);
    require(loaded.shape == grid.shape && loaded.pitch == grid.pitch &&
                loaded.origin.x == grid.origin.x && loaded.origin.y == grid.origin.y &&
                loaded.origin.z == grid.origin.z && loaded.occupancy == grid.occupancy,
            "PQSSVOX1 round trip must preserve the complete solid occupancy");
}

void testOpenSoupEdgeCannotAttractBoundaryVertices()
{
    auto grid = makeGrid(7);
    grid.occupancy[grid.index(3, 3, 3)] = 1;
    pqss_proxy_mesh::MeshModel open_soup_triangle;
    open_soup_triangle.vertices = {
        {1.0, 2.2, 2.5}, {5.0, 2.2, 2.5}, {1.0, 6.2, 2.5}};
    open_soup_triangle.triangles = {{{0, 1, 2}}};
    const auto solid = pqss_proxy_mesh::buildPhase1Solid(grid, open_soup_triangle);
    for (const auto vertex : solid.boundary.geometry.vertices)
    {
        require(std::abs(vertex.y - 2.2) > 1.0e-12,
                "an uncertified open soup edge attracted a boundary vertex");
    }
}

void testNearCoplanarSourcePatchesDoNotCreateSawTeeth()
{
    auto grid = makeGrid(7);
    grid.occupancy[grid.index(3, 3, 3)] = 1;
    pqss_proxy_mesh::MeshModel source;
    source.vertices = {
        {1.0, 1.0, 2.5}, {6.0, 1.0, 2.5}, {1.0, 6.0, 2.5},
        {2.0, 2.0, 2.51}, {4.0, 2.0, 2.51}, {2.0, 4.0, 2.51}};
    source.triangles = {{{0, 1, 2}}, {{3, 4, 5}}};
    const auto solid = pqss_proxy_mesh::buildPhase1Solid(grid, source);
    for (const auto vertex : solid.boundary.geometry.vertices)
        if (vertex.z < 3.0)
            require(std::abs(vertex.z - 2.501379310344828) < 1.0e-9,
                "near-coplanar source patches produced alternating boundary heights");
}

void testSourceProjectionStaysInsideActiveCell()
{
    auto grid = makeGrid(7);
    grid.occupancy[grid.index(3, 3, 3)] = 1;
    pqss_proxy_mesh::MeshModel source;
    source.vertices = {
        {0.0, 0.0, 1.0}, {7.0, 0.0, 1.0}, {0.0, 7.0, 1.0},
        {0.0, 0.0, 5.0}, {0.0, 7.0, 5.0}, {7.0, 0.0, 5.0}};
    source.triangles = {{{0, 1, 2}}, {{3, 4, 5}}};

    const auto solid = pqss_proxy_mesh::buildPhase1Solid(grid, source);
    const auto& boundary = solid.boundary;
    require(boundary.euler_characteristic == 2 && boundary.signed_volume > 0.0,
            "cell-bounded source projection must preserve the closed halfedge surface");
    for (const auto vertex : boundary.geometry.vertices)
        for (const double coordinate : {vertex.x, vertex.y, vertex.z})
            require(coordinate >= 2.0 && coordinate <= 4.0,
                    "source projection moved a dual vertex outside its active cell union");
    for (std::uint32_t index = 0; index < boundary.halfedges.size(); ++index)
    {
        const auto opposite = boundary.halfedges[index].opposite;
        require(opposite < boundary.halfedges.size() &&
                    boundary.halfedges[opposite].opposite == index,
                "cell-bounded source projection damaged halfedge opposites");
    }
}

void testNearestSourceProjectionPreservesHalfedgeTopology()
{
    auto grid = makeGrid(7);
    grid.occupancy[grid.index(3, 3, 3)] = 1;
    const auto solid = pqss_proxy_mesh::buildPhase1Solid(grid);

    pqss_proxy_mesh::MeshModel source;
    source.vertices = {
        {2.0, 2.0, 2.0}, {4.0, 2.0, 2.0}, {4.0, 4.0, 2.0}, {2.0, 4.0, 2.0},
        {2.0, 2.0, 4.0}, {4.0, 2.0, 4.0}, {4.0, 4.0, 4.0}, {2.0, 4.0, 4.0}};
    source.triangles = {
        {{0,2,1}},{{0,3,2}},{{4,5,6}},{{4,6,7}},
        {{0,1,5}},{{0,5,4}},{{1,2,6}},{{1,6,5}},
        {{2,3,7}},{{2,7,6}},{{3,0,4}},{{3,4,7}}};

    const auto projected = pqss_proxy_mesh::projectPhase1BoundaryToSource(
        solid.boundary, grid, grid, source);
    require(projected.geometry.triangles.size() <=
                solid.boundary.geometry.triangles.size() &&
                projected.euler_characteristic == 2 && projected.signed_volume > 0.0,
            "nearest-source projection did not preserve a closed genus-zero surface");
    for (std::uint32_t index = 0; index < projected.halfedges.size(); ++index)
    {
        const auto opposite = projected.halfedges[index].opposite;
        require(opposite < projected.halfedges.size() &&
                    projected.halfedges[opposite].opposite == index,
                "nearest-source projection produced an invalid opposite relation");
    }
    for (const auto vertex : projected.geometry.vertices)
    {
        const double source_plane_distance = std::min({
            std::abs(vertex.x - 2.0), std::abs(vertex.x - 4.0),
            std::abs(vertex.y - 2.0), std::abs(vertex.y - 4.0),
            std::abs(vertex.z - 2.0), std::abs(vertex.z - 4.0)});
        require(source_plane_distance <= 1.0e-4,
            "projected vertex exceeds the bounded source-plane residual");
    }
}

void testSourceConstrainedBoundaryKeepsExactPlane()
{
    // Two tetrahedra share z=0.  The occupied tetrahedron is below it and the
    // empty one above it, so extraction must reuse the exact constrained face
    // instead of creating voxel- or QEF-positioned vertices.
    pqss_proxy_mesh::SourceConstrainedCellComplex complex;
    complex.vertices = {
        {0.0, 0.0, 0.0}, {4.0, 0.0, 0.0}, {0.0, 3.0, 0.0},
        {0.0, 0.0, -2.0}, {0.0, 0.0, 2.0}};
    pqss_proxy_mesh::SourceConstrainedCell below;
    below.vertices = {0, 2, 1, 3};
    below.neighbors[3] = 1;
    below.source_faces[3] = 0;
    below.occupied = true;
    pqss_proxy_mesh::SourceConstrainedCell above;
    above.vertices = {0, 1, 2, 4};
    above.neighbors[3] = 0;
    above.source_faces[3] = 0;
    complex.cells = {below, above};

    const auto boundary = pqss_proxy_mesh::extractCellComplexBoundary(complex);
    require(boundary.face_source_faces.size() == boundary.geometry.triangles.size(),
            "cell-complex boundary must retain source-face provenance");
    bool found_source_face = false;
    for (std::size_t face_id = 0; face_id < boundary.geometry.triangles.size(); ++face_id)
    {
        const auto triangle = boundary.geometry.triangles[face_id];
        bool is_source_face = true;
        for (const auto vertex : triangle)
            is_source_face &= std::abs(boundary.geometry.vertices[vertex].z) < 1.0e-15;
        if (is_source_face)
            require(boundary.face_source_faces[face_id] == 0,
                    "exact source interface lost its source-face certificate");
        found_source_face |= is_source_face;
    }
    require(found_source_face,
            "cell-complex extraction must retain the exact source-constrained face");
    for (const auto vertex : boundary.geometry.vertices)
        require(vertex.z == 0.0 || vertex.z == -2.0 || vertex.z == 2.0,
                "cell-complex extraction introduced a projected vertex");
}

void testCoplanarSoupCanonicalization()
{
    pqss_proxy_mesh::MeshModel source;
    source.name = "split_rectangle";
    source.vertices = {
        {0,0,0}, {2,0,0}, {4,0,0}, {0,3,0}, {2,3,0}, {4,3,0}};
    source.triangles = {{{0,1,4}}, {{0,4,3}}, {{1,2,5}}, {{1,5,4}}};
    const auto result = pqss_proxy_mesh::canonicalizeCoplanarTriangleSoup(source);
    require(result.triangles.size() == 2,
            "a split planar rectangle must canonicalize to exactly two triangles");
    double area = 0.0;
    for (const auto face : result.triangles)
    {
        const auto a = result.vertices[face[0]];
        const auto b = result.vertices[face[1]];
        const auto c = result.vertices[face[2]];
        const double ux = b.x - a.x, uy = b.y - a.y;
        const double vx = c.x - a.x, vy = c.y - a.y;
        area += 0.5 * std::abs(ux * vy - uy * vx);
        require(a.z == 0.0 && b.z == 0.0 && c.z == 0.0,
                "canonical planar vertices must remain on the exact source plane");
    }
    require(std::abs(area - 12.0) < 1.0e-9,
            "coplanar union must preserve the covered source area");
}

void testCoplanarSoupCanonicalizationPreservesHole()
{
    pqss_proxy_mesh::MeshModel source;
    source.name = "planar_ring";
    source.vertices = {
        {0,0,0}, {4,0,0}, {4,4,0}, {0,4,0},
        {1,1,0}, {3,1,0}, {3,3,0}, {1,3,0}};
    source.triangles = {
        {{0,1,5}}, {{0,5,4}}, {{1,2,6}}, {{1,6,5}},
        {{2,3,7}}, {{2,7,6}}, {{3,0,4}}, {{3,4,7}}};
    const auto result = pqss_proxy_mesh::canonicalizeCoplanarTriangleSoup(source);
    double area = 0.0;
    for (const auto face : result.triangles)
    {
        const auto a = result.vertices[face[0]];
        const auto b = result.vertices[face[1]];
        const auto c = result.vertices[face[2]];
        area += 0.5 * std::abs(
            (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x));
        const double center_x = (a.x + b.x + c.x) / 3.0;
        const double center_y = (a.y + b.y + c.y) / 3.0;
        require(center_x <= 1.0 || center_x >= 3.0 ||
                    center_y <= 1.0 || center_y >= 3.0,
                "coplanar union filled a source hole");
    }
    require(std::abs(area - 12.0) < 1.0e-9,
            "coplanar union must preserve a polygon-with-hole area");
}

} // namespace

int main()
{
    try
    {
        testSingleCube();
        testCavity();
        testDigitalConnectivityPair();
        testTorusFill();
        testBeta2IsFilledOnlyAfterHandleClosing();
        testHandleRepairDoesNotFillUnrelatedExteriorNotch();
        testIncrementalClosingProfileConnectsSeparateComponents();
        testAnalysisHalfedgeSplitsConflictingFans();
        testPhase1SolidBoundary();
        testDualVertexDisconnectedFansAreSplit();
        testDualCheckerboardEdgeIsSplitByOccupiedSector();
        testVoxelGridRoundTrip();
        testOpenSoupEdgeCannotAttractBoundaryVertices();
        testNearCoplanarSourcePatchesDoNotCreateSawTeeth();
        testSourceProjectionStaysInsideActiveCell();
        testNearestSourceProjectionPreservesHalfedgeTopology();
        testSourceConstrainedBoundaryKeepsExactPlane();
        testCoplanarSoupCanonicalization();
        testCoplanarSoupCanonicalizationPreservesHole();
        std::cout << "topology fill tests passed\n";
        return EXIT_SUCCESS;
    }
    catch (const std::exception& exception)
    {
        std::cerr << exception.what() << '\n';
        return EXIT_FAILURE;
    }
}
