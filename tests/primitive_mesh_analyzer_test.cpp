#include "pqss_proxy_mesh/primitive_mesh_analyzer.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
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

void appendBoxObj(std::ostringstream& obj, std::size_t& next_vertex,
                  const double x0, const double y0, const double z0,
                  const double x1, const double y1, const double z1)
{
    const std::size_t base = next_vertex;
    obj << "v " << x0 << ' ' << y0 << ' ' << z0 << '\n'
        << "v " << x1 << ' ' << y0 << ' ' << z0 << '\n'
        << "v " << x1 << ' ' << y1 << ' ' << z0 << '\n'
        << "v " << x0 << ' ' << y1 << ' ' << z0 << '\n'
        << "v " << x0 << ' ' << y0 << ' ' << z1 << '\n'
        << "v " << x1 << ' ' << y0 << ' ' << z1 << '\n'
        << "v " << x1 << ' ' << y1 << ' ' << z1 << '\n'
        << "v " << x0 << ' ' << y1 << ' ' << z1 << '\n';
    const auto face = [&](const std::size_t a, const std::size_t b,
                          const std::size_t c)
    {
        obj << "f " << base + a << ' ' << base + b << ' ' << base + c
            << '\n';
    };
    face(0, 1, 2); face(0, 2, 3);
    face(4, 6, 5); face(4, 7, 6);
    face(0, 4, 5); face(0, 5, 1);
    face(1, 5, 6); face(1, 6, 2);
    face(2, 6, 7); face(2, 7, 3);
    face(3, 7, 4); face(3, 4, 0);
    next_vertex += 8;
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

        // A spatial group is not restricted to coplanar surfaces. A loose
        // directed-distance limit may replace this non-convex L assembly by
        // one conservative six-face box; an exact limit must reject it.
        const auto spatial_group = root / "spatial_group.obj";
        writeText(spatial_group,
            "v 0 0 0\nv 2 0 0\nv 2 1 0\nv 1 1 0\nv 1 3 0\nv 0 3 0\n"
            "v 0 0 1\nv 2 0 1\nv 2 1 1\nv 1 1 1\nv 1 3 1\nv 0 3 1\n"
            "f 1 3 2\nf 1 4 3\nf 1 5 4\nf 1 6 5\n"
            "f 7 8 9\nf 7 9 10\nf 7 10 11\nf 7 11 12\n"
            "f 1 2 8\nf 1 8 7\nf 2 3 9\nf 2 9 8\n"
            "f 3 4 10\nf 3 10 9\nf 4 5 11\nf 4 11 10\n"
            "f 5 6 12\nf 5 12 11\nf 6 1 7\nf 6 7 12\n");

        auto group_loose = exact;
        group_loose.maximum_cavity_added_volume_ratio = 0.0;
        group_loose.maximum_open_error_distance = 10.0;
        const auto group_loose_stats =
            pqss_proxy_mesh::analyzePrimitiveMeshObj(
                spatial_group, root / "spatial_group_loose", group_loose);
        require(group_loose_stats.containment_validation_passed &&
                    group_loose_stats.primitive_count == 6 &&
                    group_loose_stats.proxy_triangles == 12 &&
                    group_loose_stats.merged_spatial_primitive_groups > 0,
                "loose error must permit a noncoplanar spatial group box");

        auto group_strict = group_loose;
        group_strict.maximum_open_error_distance = 0.0;
        const auto group_strict_stats =
            pqss_proxy_mesh::analyzePrimitiveMeshObj(
                spatial_group, root / "spatial_group_strict", group_strict);
        require(group_strict_stats.merged_spatial_primitive_groups == 0 &&
                    group_strict_stats.primitive_count > 6,
                "exact error must reject the noncoplanar spatial group box");

        // Three equal closed components contain two independent narrow gaps.
        // Each gap is assessed against the model volume on its own.
        const auto narrow_gaps = root / "narrow_component_gaps.obj";
        std::ostringstream narrow_gap_obj;
        std::size_t next_vertex = 1;
        appendBoxObj(narrow_gap_obj, next_vertex, 0.0, 0, 0, 1.0, 1, 1);
        appendBoxObj(narrow_gap_obj, next_vertex, 1.1, 0, 0, 2.1, 1, 1);
        appendBoxObj(narrow_gap_obj, next_vertex, 2.2, 0, 0, 3.2, 1, 1);
        writeText(narrow_gaps, narrow_gap_obj.str());
        auto gap_fill = exact;
        gap_fill.maximum_open_error_distance = 0.0;
        gap_fill.maximum_cavity_added_volume_ratio = 0.04;
        const auto gap_fill_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
            narrow_gaps, root / "narrow_component_gaps", gap_fill);
        require(gap_fill_stats.filled_intercomponent_gaps == 2,
                "two narrow component gaps must be filled independently");

        // A large empty region is not a narrow gap and must survive the same
        // per-gap budget.
        const auto large_gap = root / "large_component_gap.obj";
        std::ostringstream large_gap_obj;
        next_vertex = 1;
        appendBoxObj(large_gap_obj, next_vertex, 0, 0, 0, 1, 1, 1);
        appendBoxObj(large_gap_obj, next_vertex, 3, 0, 0, 4, 1, 1);
        writeText(large_gap, large_gap_obj.str());
        const auto large_gap_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
            large_gap, root / "large_component_gap", gap_fill);
        require(large_gap_stats.filled_intercomponent_gaps == 0,
                "a component gap beyond the per-gap budget must be preserved");

        std::cout << "primitive mesh staged-pipeline tests passed\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
