#include "pqss_proxy_mesh/phase0_open_hole.hpp"

#include <cstdlib>
#include <iostream>
#include <stdexcept>

namespace
{

void require(const bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}

void testPlanarBoundaryLoopIsCoveredWithoutChangingSourceFaces()
{
    pqss_proxy_mesh::MeshModel open_box;
    open_box.name = "open_box";
    open_box.vertices = {
        {0,0,0}, {2,0,0}, {2,2,0}, {0,2,0},
        {0,0,1}, {2,0,1}, {2,2,1}, {0,2,1}};
    open_box.triangles = {
        {{0,2,1}}, {{0,3,2}},
        {{0,1,5}}, {{0,5,4}},
        {{1,2,6}}, {{1,6,5}},
        {{2,3,7}}, {{2,7,6}},
        {{3,0,4}}, {{3,4,7}}};
    const auto result = pqss_proxy_mesh::coverPlanarOpenHoles(open_box);
    require(result.stats.closed_boundary_loops == 1, "top opening was not recognized");
    require(result.stats.capped_loops == 1 && result.stats.cap_triangles == 2,
            "four-edge opening did not produce exactly two cap triangles");
    require(result.covered_mesh.triangles.size() == open_box.triangles.size() + 2,
            "phase 0 changed source triangle responsibility");
    for (std::size_t face = 0; face < open_box.triangles.size(); ++face)
        require(result.covered_mesh.triangles[face] == open_box.triangles[face],
                "phase 0 reordered or replaced an original source triangle");
}

void testOpenSheetOutlineIsNotMistakenForClosedRealityLoopAfterBranching()
{
    pqss_proxy_mesh::MeshModel branched;
    branched.name = "branched_boundary";
    branched.vertices = {
        {0,0,0}, {1,0,0}, {0,1,0}, {1,1,0}, {2,0,0}};
    branched.triangles = {
        {{0,1,2}}, {{1,3,2}}, {{1,4,3}}};
    const auto result = pqss_proxy_mesh::coverPlanarOpenHoles(branched);
    require(result.stats.capped_loops == 1,
            "a simple planar open-sheet boundary must remain a recognized boundary loop");
    require(result.stats.cap_triangles == 3,
            "five-edge planar boundary must triangulate to n-2 faces");
}

} // namespace

int main()
{
    try
    {
        testPlanarBoundaryLoopIsCoveredWithoutChangingSourceFaces();
        testOpenSheetOutlineIsNotMistakenForClosedRealityLoopAfterBranching();
        std::cout << "phase0 open-hole tests passed\n";
        return EXIT_SUCCESS;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
