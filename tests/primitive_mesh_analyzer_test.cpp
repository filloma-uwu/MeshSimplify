#include "pqss_proxy_mesh/primitive_mesh_analyzer.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>

namespace
{

void writeText(const std::filesystem::path& path, const std::string& text)
{
    std::ofstream stream(path);
    if (!stream) throw std::runtime_error("cannot create test OBJ");
    stream << text;
}

void require(const bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}

std::string readText(const std::filesystem::path& path)
{
    std::ifstream stream(path);
    return {std::istreambuf_iterator<char>(stream),
            std::istreambuf_iterator<char>()};
}

} // namespace

int main()
{
    try
    {
        const auto root = std::filesystem::temp_directory_path() /
            "pqss_primitive_mesh_analyzer_test";
        std::filesystem::remove_all(root);
        std::filesystem::create_directories(root);

        const auto box = root / "box.obj";
        writeText(box,
            "v 0 0 0\nv 1 0 0\nv 1 1 0\nv 0 1 0\n"
            "v 0 0 1\nv 1 0 1\nv 1 1 1\nv 0 1 1\n"
            "f 1 2 3\nf 1 3 4\nf 5 7 6\nf 5 8 7\n"
            "f 1 5 6\nf 1 6 2\nf 2 6 7\nf 2 7 3\n"
            "f 3 7 8\nf 3 8 4\nf 4 8 5\nf 4 5 1\n");
        pqss_proxy_mesh::PrimitiveMeshAnalysisOptions exact;
        exact.allow_frustum = false;
        exact.maximum_open_error_distance = 0.0;
        const auto box_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
            box, root / "box", exact);
        require(box_stats.containment_validation_passed,
                "exact box must pass conservative coverage");
        require(box_stats.primitive_count == 6 &&
                    box_stats.proxy_triangles == 12,
                "box must become six polygon surfaces and twelve triangles");
        require(box_stats.open_max_distance <= 1.0e-9,
                "exact box must have zero phase3-to-phase1 surface error");
        for (const char* file : {
                 "phase1_hole_filled.obj",
                 "phase2_recognized_surfaces.obj",
                 "primitives.obj",
                 "proxy.obj"})
            require(std::filesystem::is_regular_file(root / "box" / file),
                    "every pipeline phase must have a viewer asset");
        const std::string box_metadata = readText(root / "box" / "model.json");
        require(box_metadata.find(
                    "sampled_directed_phase3_to_phase1_surface_distance") !=
                    std::string::npos &&
                box_metadata.find(
                    "\"reference\":\"phase1_hole_filled.obj\"") !=
                    std::string::npos,
                "metadata must state that simplification error uses phase 1");

        // Three disconnected coplanar rectangles need two accepted merges.
        // The filled gaps are 0.1 wide, so the directed maximum is 0.05.
        const auto gaps = root / "coplanar_gaps.obj";
        writeText(gaps,
            "v 0 0 0\nv 1 0 0\nv 1 1 0\nv 0 1 0\n"
            "v 1.1 0 0\nv 2.1 0 0\nv 2.1 1 0\nv 1.1 1 0\n"
            "v 2.2 0 0\nv 3.2 0 0\nv 3.2 1 0\nv 2.2 1 0\n"
            "f 1 2 3\nf 1 3 4\nf 5 6 7\nf 5 7 8\n"
            "f 9 10 11\nf 9 11 12\n");

        auto loose = exact;
        loose.maximum_open_error_distance = 0.051;
        const auto loose_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
            gaps, root / "gaps_loose", loose);
        require(loose_stats.containment_validation_passed,
                "merged coplanar proxy must remain conservative");
        require(loose_stats.primitive_count == 1 &&
                    loose_stats.proxy_triangles == 2 &&
                    loose_stats.merged_local_planar_primitives == 2,
                "Hausdorff-compliant merging must iterate to a fixed point");
        require(loose_stats.open_max_distance <= 0.051 + 1.0e-9,
                "accepted fixed-point merge must respect the user limit");

        auto strict = exact;
        strict.maximum_open_error_distance = 0.049;
        const auto strict_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
            gaps, root / "gaps_strict", strict);
        require(strict_stats.primitive_count == 3 &&
                    strict_stats.proxy_triangles == 6 &&
                    strict_stats.merged_local_planar_primitives == 0,
                "a merge beyond the directed Hausdorff limit must be rejected");

        std::cout << "primitive mesh staged-pipeline tests passed\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
