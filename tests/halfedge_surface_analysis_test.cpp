#include "pqss_proxy_mesh/halfedge_surface_analysis.hpp"
#include "pqss_proxy_mesh/hausdorff_simplifier.hpp"
#include "pqss_proxy_mesh/halfedge_validation.hpp"

#include <array>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace
{
void require(const bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}

// A planar N x N quad grid in the z = 0 plane, triangulated into 2 * N * N
// triangles. Coplanar union must reduce it to exactly n-2 = 2 triangles.
pqss_proxy_mesh::MeshModel makeSubdividedQuad(const std::size_t cells)
{
    pqss_proxy_mesh::MeshModel mesh;
    mesh.name = "subdivided_quad";
    for (std::size_t y = 0; y <= cells; ++y)
        for (std::size_t x = 0; x <= cells; ++x)
            mesh.vertices.push_back({2.0 * x / cells - 1.0, 2.0 * y / cells - 1.0, 0.0});
    const auto index = [&](std::size_t x, std::size_t y)
    {
        return static_cast<std::uint32_t>(y * (cells + 1) + x);
    };
    for (std::size_t y = 0; y < cells; ++y)
        for (std::size_t x = 0; x < cells; ++x)
        {
            mesh.triangles.push_back({index(x, y), index(x + 1, y), index(x + 1, y + 1)});
            mesh.triangles.push_back({index(x, y), index(x + 1, y + 1), index(x, y + 1)});
        }
    return mesh;
}

pqss_proxy_mesh::MeshModel makeCube()
{
    pqss_proxy_mesh::MeshModel cube;
    cube.name = "cube";
    cube.vertices = {{0,0,0},{0,0,1},{0,1,0},{0,1,1},
                     {1,0,0},{1,0,1},{1,1,0},{1,1,1}};
    cube.triangles = {{{0,2,3}},{{0,3,1}},{{4,5,7}},{{4,7,6}},
                      {{0,1,5}},{{0,5,4}},{{2,6,7}},{{2,7,3}},
                      {{0,4,6}},{{0,6,2}},{{1,3,7}},{{1,7,5}}};
    return cube;
}

// Two adjacent planar strips meet at a shallow ridge. The four source
// triangles have a two-triangle support-plane hull at a sufficiently loose
// directed distance limit.
pqss_proxy_mesh::MeshModel makeShallowFold()
{
    pqss_proxy_mesh::MeshModel mesh;
    mesh.name = "shallow_fold";
    mesh.vertices = {{0,-1,0},{0,1,0},{1,-1,0.05},{1,1,0.05},{2,-1,0},{2,1,0}};
    mesh.triangles = {{{0,2,3}},{{0,3,1}},{{2,4,5}},{{2,5,3}}};
    return mesh;
}

// An L-shaped planar patch. Its full convex hull fills the missing quadrant,
// so a tight Hausdorff limit must reject that non-convex-region merge.
pqss_proxy_mesh::MeshModel makeNonConvexPatch()
{
    pqss_proxy_mesh::MeshModel mesh;
    mesh.name = "non_convex_patch";
    mesh.vertices = {{0,0,0},{1,0,0},{2,0,0},
                     {0,1,0},{1,1,0},{2,1,0},
                     {0,2,0},{1,2,0}};
    mesh.triangles = {{{0,1,4}},{{0,4,3}},{{1,2,5}},{{1,5,4}},
                      {{3,4,7}},{{3,7,6}}};
    return mesh;
}

// The triangles use opposite directions across their paired edge and have
// opposing normals. They must not acquire one averaged support plane.
pqss_proxy_mesh::MeshModel makeOpposingAdjacentFaces()
{
    pqss_proxy_mesh::MeshModel mesh;
    mesh.name = "opposing_adjacent_faces";
    mesh.vertices = {{0,0,0},{1,0,0},{0,1,0},{0,1,0.1}};
    mesh.triangles = {{{0,1,2}},{{1,0,3}}};
    return mesh;
}

} // namespace

int main()
{
    try
    {
        const auto temporary = std::filesystem::temp_directory_path() /
            "pqss_halfedge_surface_analysis_test";
        std::filesystem::create_directories(temporary);

        // Subdivided coplanar quad: 2 * N * N triangles must collapse to 2.
        const pqss_proxy_mesh::MeshModel quad = makeSubdividedQuad(4);
        const pqss_proxy_mesh::OrientedSurfaceMesh quad_halfedge =
            pqss_proxy_mesh::buildAnalysisHalfedgeMesh(quad);
        pqss_proxy_mesh::HalfedgeSurfaceAnalysisOptions options;
        options.maximum_directed_hausdorff = 0.01;
        const auto quad_stats = pqss_proxy_mesh::analyzeHalfedgeSurface(
            quad_halfedge, temporary / "quad", options);
        require(quad_stats.final_triangles == 2,
                "subdivided coplanar quad must reduce to exactly two triangles");
        require(quad_stats.coverage_certified,
                "zero-error coplanar union must pass the area coverage certificate");
        require(quad_stats.global_hausdorff_certified &&
                    quad_stats.certified_directed_hausdorff_upper_bound <=
                        options.maximum_directed_hausdorff,
                "assembled coplanar union must pass its global Hausdorff certificate");

        // Closed cube remains conservatively covered. Support-plane triangles
        // may also bridge a convex edge when the directed limit certifies it.
        const pqss_proxy_mesh::MeshModel cube = makeCube();
        const pqss_proxy_mesh::OrientedSurfaceMesh cube_halfedge =
            pqss_proxy_mesh::buildAnalysisHalfedgeMesh(cube);
        const auto cube_stats = pqss_proxy_mesh::analyzeHalfedgeSurface(
            cube_halfedge, temporary / "cube", options);
        require(cube_stats.final_triangles > 0 && cube_stats.final_triangles <= 12,
                "closed cube support-plane analysis must not increase triangle count");
        require(cube_stats.coverage_certified,
                "cube coplanar union must pass the area coverage certificate");

        const pqss_proxy_mesh::MeshModel fold = makeShallowFold();
        const pqss_proxy_mesh::OrientedSurfaceMesh fold_halfedge =
            pqss_proxy_mesh::buildAnalysisHalfedgeMesh(fold);
        pqss_proxy_mesh::HalfedgeSurfaceAnalysisOptions strict_options;
        strict_options.maximum_directed_hausdorff = 0.005;
        const auto strict_stats = pqss_proxy_mesh::analyzeHalfedgeSurface(
            fold_halfedge, temporary / "fold_strict", strict_options);
        pqss_proxy_mesh::HalfedgeSurfaceAnalysisOptions loose_options;
        loose_options.maximum_directed_hausdorff = 0.1;
        const auto loose_stats = pqss_proxy_mesh::analyzeHalfedgeSurface(
            fold_halfedge, temporary / "fold_loose", loose_options);
        require(strict_stats.final_triangles == fold.triangles.size(),
                "strict limit must retain the shallow fold");
        require(loose_stats.final_triangles < strict_stats.final_triangles,
                "loose limit must simplify adjacent shallow-fold regions");
        require(strict_stats.coverage_certified && loose_stats.coverage_certified,
                "both shallow-fold results must retain conservative coverage");
        require(strict_stats.global_hausdorff_certified &&
                    loose_stats.global_hausdorff_certified,
                "both assembled shallow-fold proxies must pass global certification");

        const pqss_proxy_mesh::MeshModel non_convex = makeNonConvexPatch();
        const auto non_convex_stats = pqss_proxy_mesh::analyzeHalfedgeSurface(
            pqss_proxy_mesh::buildAnalysisHalfedgeMesh(non_convex),
            temporary / "non_convex", strict_options);
        require(non_convex_stats.final_triangles > 2,
                "tight limit must not replace an L-shaped patch by its full convex hull");
        require(non_convex_stats.rejected_hausdorff_merges > 0,
                "non-convex hull fill must be rejected by the Hausdorff certificate");
        require(non_convex_stats.coverage_certified,
                "non-convex fallback children must retain coverage");

        const pqss_proxy_mesh::MeshModel opposing = makeOpposingAdjacentFaces();
        const auto opposing_stats = pqss_proxy_mesh::analyzeHalfedgeSurface(
            pqss_proxy_mesh::buildAnalysisHalfedgeMesh(opposing),
            temporary / "opposing", loose_options);
        require(opposing_stats.final_triangles == opposing.triangles.size(),
                "opposing face directions must not merge");
        require(opposing_stats.rejected_orientation_merges > 0,
                "opposing face directions must fail orientation consistency");
        require(opposing_stats.coverage_certified,
                "opposing exact fallbacks must retain coverage");

        std::cout << "halfedge surface analysis tests passed\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
