#include "pqss_proxy_mesh/hausdorff_simplifier.hpp"
#include "pqss_proxy_mesh/halfedge_validation.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace
{
void require(const bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}
}

int main()
{
    try
    {
        pqss_proxy_mesh::MeshModel reference;
        reference.vertices = {{-2,-2,0},{2,-2,0},{2,2,0},{-2,2,0}};
        reference.triangles = {{{0,1,2}},{{0,2,3}}};
        auto exact = reference;
        const auto exact_certificate = pqss_proxy_mesh::certifyDirectedHausdorff(
            exact, reference, 0.0, 12, 1.0e-10);
        require(exact_certificate.passed, "identical meshes must have zero directed Hausdorff distance");

        auto offset = reference;
        for (auto& point : offset.vertices) point.z = 0.5;
        const auto loose = pqss_proxy_mesh::certifyDirectedHausdorff(
            offset, reference, 0.51, 12, 1.0e-10);
        require(loose.passed && loose.upper_bound <= 0.51 + 1.0e-8,
                "parallel offset inside the limit must be certified");
        const auto strict = pqss_proxy_mesh::certifyDirectedHausdorff(
            offset, reference, 0.49, 12, 1.0e-10);
        require(!strict.passed && strict.lower_bound >= 0.5 - 1.0e-10,
                "parallel offset beyond the limit must be rejected");

        pqss_proxy_mesh::MeshModel cube;
        cube.vertices = {{0,0,0},{0,0,1},{0,1,0},{0,1,1},
                         {1,0,0},{1,0,1},{1,1,0},{1,1,1}};
        cube.triangles = {{{0,2,3}},{{0,3,1}},{{4,5,7}},{{4,7,6}},
                          {{0,1,5}},{{0,5,4}},{{2,6,7}},{{2,7,3}},
                          {{0,4,6}},{{0,6,2}},{{1,3,7}},{{1,7,5}}};
        for (auto& triangle : cube.triangles)
            std::swap(triangle[1], triangle[2]);
        pqss_proxy_mesh::AnalysisHalfedgeStats halfedge_stats;
        const auto cube_halfedge = pqss_proxy_mesh::buildAnalysisHalfedgeMesh(cube, &halfedge_stats);
        require(halfedge_stats.boundary_halfedges == 0,
                "cube regression must form a closed halfedge mesh");
        const auto validation =
            pqss_proxy_mesh::validatePhase1HalfedgeTopology(cube_halfedge);
        require(validation.euler_characteristic == 2 &&
                    validation.face_components == 1,
                "explicit halfedge validation must accept the closed cube");
        const auto temporary = std::filesystem::temp_directory_path() /
            "pqss_hausdorff_simplifier_cube";
        std::filesystem::create_directories(temporary);
        pqss_proxy_mesh::writeAnalysisHalfedgeMesh(temporary / "phase1.bin", cube_halfedge);
        auto damaged = cube_halfedge;
        damaged.halfedges.front().opposite = pqss_proxy_mesh::invalid_surface_index;
        pqss_proxy_mesh::writeAnalysisHalfedgeMesh(
            temporary / "unchecked_phase1.bin", damaged);
        const auto unchecked = pqss_proxy_mesh::readAnalysisHalfedgeMesh(
            temporary / "unchecked_phase1.bin");
        require(unchecked.halfedges.front().opposite ==
                    pqss_proxy_mesh::invalid_surface_index,
                "production halfedge loading must not run topology validation");
        bool rejected = false;
        try
        {
            (void)pqss_proxy_mesh::validatePhase1HalfedgeTopology(unchecked);
        }
        catch (const std::runtime_error&)
        {
            rejected = true;
        }
        require(rejected,
                "explicit halfedge validation must reject a damaged opposite relation");
        {
            std::ofstream source(temporary / "source.obj");
            source << "o cube\n";
        }
        pqss_proxy_mesh::HausdorffSimplificationOptions cube_options;
        cube_options.maximum_directed_hausdorff = 0.0;
        const auto cube_result = pqss_proxy_mesh::simplifyPhase1Halfedge(
            temporary / "phase1.bin", temporary / "source.obj",
            temporary / "output", cube_options);
        require(cube_result.final_triangles == 12,
                "zero-error cube must retain the minimum six triangulated polygon faces");
        std::cout << "hausdorff simplifier tests passed\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
