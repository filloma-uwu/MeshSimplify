#include "pqss_proxy_mesh/hausdorff_simplifier.hpp"
#include "pqss_proxy_mesh/phase2_halfedge.hpp"

#include <cstdlib>
#include <array>
#include <filesystem>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>

namespace
{

void require(const bool condition, const std::string& message)
{
    if (!condition) throw std::runtime_error(message);
}

pqss_proxy_mesh::MeshModel makeSubdividedClosedBox(const std::uint32_t cells)
{
    pqss_proxy_mesh::MeshModel box;
    box.name = "subdivided_closed_box";
    std::map<std::array<double, 3>, std::uint32_t> vertices;
    const auto vertex = [&](const double x, const double y, const double z)
    {
        const std::array<double, 3> key{x, y, z};
        const auto [iterator, inserted] = vertices.try_emplace(
            key, static_cast<std::uint32_t>(box.vertices.size()));
        if (inserted) box.vertices.push_back({x, y, z});
        return iterator->second;
    };
    const auto addFace = [&](const std::array<double, 3> origin,
                             const std::array<double, 3> u,
                             const std::array<double, 3> v)
    {
        for (std::uint32_t y = 0; y < cells; ++y)
            for (std::uint32_t x = 0; x < cells; ++x)
            {
                const auto point = [&](const std::uint32_t px, const std::uint32_t py)
                {
                    const double sx = static_cast<double>(px) / cells;
                    const double sy = static_cast<double>(py) / cells;
                    return vertex(origin[0] + u[0] * sx + v[0] * sy,
                                  origin[1] + u[1] * sx + v[1] * sy,
                                  origin[2] + u[2] * sx + v[2] * sy);
                };
                const std::uint32_t a = point(x, y);
                const std::uint32_t b = point(x + 1, y);
                const std::uint32_t c = point(x + 1, y + 1);
                const std::uint32_t d = point(x, y + 1);
                box.triangles.push_back({a, b, c});
                box.triangles.push_back({a, c, d});
            }
    };

    addFace({0,0,0}, {0,1,0}, {1,0,0});
    addFace({0,0,1}, {1,0,0}, {0,1,0});
    addFace({0,0,0}, {0,0,1}, {0,1,0});
    addFace({1,0,0}, {0,1,0}, {0,0,1});
    addFace({0,0,0}, {1,0,0}, {0,0,1});
    addFace({0,1,0}, {0,0,1}, {1,0,0});
    return box;
}

void testExactPhase1FallbackIsAValidPhase2Artifact()
{
    pqss_proxy_mesh::MeshModel rectangle;
    rectangle.name = "tetrahedron_with_near_duplicate";
    rectangle.vertices = {
        {0,0,0}, {4,0,0}, {0,4,0}, {0,0,4}, {1.0e-10,0,0}};
    rectangle.triangles = {
        {{4,2,1}}, {{0,1,3}}, {{1,2,3}}, {{2,0,3}}};
    pqss_proxy_mesh::AnalysisHalfedgeStats input_stats;
    auto input = pqss_proxy_mesh::buildAnalysisHalfedgeMesh(rectangle, &input_stats);
    const std::filesystem::path root = std::filesystem::temp_directory_path() /
        "pqss_phase2_halfedge_regression";
    std::error_code ignored;
    std::filesystem::remove_all(root, ignored);
    std::filesystem::create_directories(root);
    const std::filesystem::path input_path = root / "phase1_halfedge.bin";
    pqss_proxy_mesh::writeAnalysisHalfedgeMesh(input_path, input);

    pqss_proxy_mesh::Phase2HalfedgeOptions options;
    options.model_id = "rectangle";
    options.maximum_directed_hausdorff = 0.01;
    const auto stats = pqss_proxy_mesh::generatePhase2Halfedge(
        input_path, root / "output", options);
    require(stats.output_triangles == 4 && stats.boundary_halfedges == 0 &&
                stats.face_components == 1,
            "phase-2 exact fallback did not produce one closed component");
    const auto output = pqss_proxy_mesh::readAnalysisHalfedgeMesh(
        root / "output" / "phase2_halfedge.bin");
    const auto validation = pqss_proxy_mesh::validateWorkingHalfedgeTopology(output);
    require(validation.triangles == 4 && validation.boundary_halfedges == 0 &&
                validation.face_components == 1,
            "serialized phase-2 halfedge failed validation");
    std::filesystem::remove_all(root, ignored);
}

void testMismatchedOppositeIsRejected()
{
    pqss_proxy_mesh::MeshModel tetrahedron;
    tetrahedron.vertices = {{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}};
    tetrahedron.triangles = {
        {{0,2,1}}, {{0,1,3}}, {{1,2,3}}, {{2,0,3}}};
    auto mesh = pqss_proxy_mesh::buildAnalysisHalfedgeMesh(tetrahedron);
    mesh.halfedges[0].opposite = 1;
    bool rejected = false;
    try
    {
        (void)pqss_proxy_mesh::validateWorkingHalfedgeTopology(mesh);
    }
    catch (const std::exception&)
    {
        rejected = true;
    }
    require(rejected, "geometrically mismatched opposite was accepted");
}

void testOpenInputIsRejectedWithoutFabricatedOpposite()
{
    pqss_proxy_mesh::MeshModel open_tetrahedron;
    open_tetrahedron.name = "open_tetrahedron";
    open_tetrahedron.vertices = {{0,0,0}, {3,0,0}, {0,3,0}, {0,0,3}};
    open_tetrahedron.triangles = {
        {{0,2,1}}, {{0,1,3}}, {{2,0,3}}};
    const auto input = pqss_proxy_mesh::buildAnalysisHalfedgeMesh(open_tetrahedron);
    const std::filesystem::path root = std::filesystem::temp_directory_path() /
        "pqss_phase2_boundary_closure_regression";
    std::error_code ignored;
    std::filesystem::remove_all(root, ignored);
    std::filesystem::create_directories(root);
    pqss_proxy_mesh::writeAnalysisHalfedgeMesh(root / "phase1_halfedge.bin", input);
    bool rejected = false;
    try
    {
        (void)pqss_proxy_mesh::generatePhase2Halfedge(
            root / "phase1_halfedge.bin", root / "output");
    }
    catch (const std::exception&)
    {
        rejected = true;
    }
    require(rejected,
            "open input was silently closed with a fabricated opposite");
    std::filesystem::remove_all(root, ignored);
}

void testGlobalCoplanarUnionHandlesSubdividedClosedPlanes()
{
    const auto box = makeSubdividedClosedBox(2);
    require(box.triangles.size() == 48, "subdivided closed box fixture changed");
    const auto input = pqss_proxy_mesh::buildAnalysisHalfedgeMesh(box);
    const std::filesystem::path root = std::filesystem::temp_directory_path() /
        "pqss_phase2_global_coplanar_union_regression";
    std::error_code ignored;
    std::filesystem::remove_all(root, ignored);
    std::filesystem::create_directories(root);
    pqss_proxy_mesh::writeAnalysisHalfedgeMesh(root / "phase1_halfedge.bin", input);

    pqss_proxy_mesh::Phase2HalfedgeOptions options;
    options.model_id = "subdivided_box";
    options.maximum_directed_hausdorff = 0.0;
    const auto stats = pqss_proxy_mesh::generatePhase2Halfedge(
        root / "phase1_halfedge.bin", root / "output", options);
    require(stats.output_triangles == 12,
            "global coplanar union must collapse each box side to two triangles");
    require(stats.boundary_halfedges == 0 && stats.face_components == 1,
            "global coplanar union must remain one closed halfedge component");
    const auto output = pqss_proxy_mesh::readAnalysisHalfedgeMesh(
        root / "output" / "phase2_halfedge.bin");
    const auto validation = pqss_proxy_mesh::validateWorkingHalfedgeTopology(output);
    require(validation.triangles == 12 && validation.boundary_halfedges == 0,
            "global coplanar union output failed phase-2 halfedge validation");
    std::filesystem::remove_all(root, ignored);
}

} // namespace

int main()
{
    try
    {
        testExactPhase1FallbackIsAValidPhase2Artifact();
        testGlobalCoplanarUnionHandlesSubdividedClosedPlanes();
        testMismatchedOppositeIsRejected();
        testOpenInputIsRejectedWithoutFabricatedOpposite();
        std::cout << "phase2 halfedge tests passed\n";
        return EXIT_SUCCESS;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
