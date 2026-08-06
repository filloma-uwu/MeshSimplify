#include "pqss_proxy_mesh/primitive_mesh_analyzer.hpp"

#include <cmath>
#include <array>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

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

void appendCylinderBandObj(std::ostringstream& obj, std::size_t& next_vertex,
                           const std::size_t segments,
                           const double radius, const double height)
{
    const std::size_t base = next_vertex;
    constexpr double pi = 3.14159265358979323846;
    for (std::size_t index = 0; index < segments; ++index)
    {
        const double angle = 2.0 * pi * static_cast<double>(index) /
                             static_cast<double>(segments);
        const double x = radius * std::cos(angle);
        const double y = radius * std::sin(angle);
        obj << "v " << x << ' ' << y << " 0\n"
            << "v " << x << ' ' << y << ' ' << height << '\n';
    }
    for (std::size_t index = 0; index < segments; ++index)
    {
        const std::size_t next = (index + 1) % segments;
        const std::size_t lower = base + 2 * index;
        const std::size_t upper = lower + 1;
        const std::size_t next_lower = base + 2 * next;
        const std::size_t next_upper = next_lower + 1;
        obj << "f " << lower << ' ' << next_lower << ' ' << next_upper << '\n'
            << "f " << lower << ' ' << next_upper << ' ' << upper << '\n';
    }
    next_vertex += 2 * segments;
}

bool hasLargeFaceOnX(const std::filesystem::path& path, const double x)
{
    std::ifstream stream(path);
    std::vector<std::array<double, 3>> vertices;
    std::string line;
    while (std::getline(stream, line))
    {
        std::istringstream fields(line);
        std::string tag;
        fields >> tag;
        if (tag == "v")
        {
            std::array<double, 3> vertex{};
            fields >> vertex[0] >> vertex[1] >> vertex[2];
            vertices.push_back(vertex);
            continue;
        }
        if (tag != "f") continue;
        std::vector<std::size_t> face;
        std::string token;
        while (fields >> token)
            face.push_back(static_cast<std::size_t>(
                std::stoull(token.substr(0, token.find('/'))) - 1));
        if (face.size() < 3) continue;
        double lower_y = 1.0e30;
        double upper_y = -1.0e30;
        double lower_z = 1.0e30;
        double upper_z = -1.0e30;
        bool on_plane = true;
        for (const auto index : face)
        {
            if (index >= vertices.size() ||
                std::abs(vertices[index][0] - x) > 1.0e-8)
            {
                on_plane = false;
                break;
            }
            lower_y = std::min(lower_y, vertices[index][1]);
            upper_y = std::max(upper_y, vertices[index][1]);
            lower_z = std::min(lower_z, vertices[index][2]);
            upper_z = std::max(upper_z, vertices[index][2]);
        }
        if (on_plane && upper_y - lower_y > 0.9 &&
            upper_z - lower_z > 0.9)
            return true;
    }
    return false;
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
        exact.allow_round_surfaces = false;
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
                "metadata must use the phase-1 error reference");

        // A complete analytic side surface is a surface primitive, not a solid
        // cylinder candidate.  Its outward circumscribed tessellation must pass
        // coverage directly without falling back to the source triangles.
        const auto cylinder_band = root / "cylinder_band.obj";
        std::ostringstream cylinder_band_obj;
        std::size_t cylinder_vertex = 1;
        appendCylinderBandObj(
            cylinder_band_obj, cylinder_vertex, 24, 1.0, 2.0);
        writeText(cylinder_band, cylinder_band_obj.str());
        auto analytic = exact;
        analytic.allow_round_surfaces = true;
        const auto cylinder_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
            cylinder_band, root / "cylinder_band", analytic);
        require(cylinder_stats.containment_validation_passed &&
                    cylinder_stats.cylindrical_band_count == 1 &&
                    cylinder_stats.primitive_count == 1 &&
                    cylinder_stats.proxy_triangles == 48,
                "a complete cylinder side must remain one certified surface");

        // Error tolerance does not create adjacency. Three disconnected
        // coplanar rectangles remain independent stage-3 surfaces.
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
                "disconnected coplanar proxy must remain conservative");
        require(loose_stats.primitive_count == 3 &&
                    loose_stats.proxy_triangles == 6 &&
                    loose_stats.merged_local_planar_primitives == 0,
                "maximum error must not turn separated surfaces into neighbors");
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

        // Stage 3 contains surface candidates only.  A large error limit must
        // not introduce a fitted box solid around a non-convex surface set.
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
                    group_loose_stats.primitive_count > 0 &&
                    group_loose_stats.merged_spatial_primitive_groups == 0,
                "a loose error must retain surface candidates and never create a solid candidate");
        require(readText(root / "spatial_group_loose" /
                         "surface_merge_profile.json").find(
                    "\"remaining_acceptable_candidates\":0") !=
                    std::string::npos,
                "adjacent surface merging must terminate at a verified fixed point");
        require(readText(root / "spatial_group_loose" /
                         "surface_merge_profile.json").find(
                    "\"strategy\":\"streaming_greedy_adjacent_fixed_point\"") !=
                    std::string::npos,
                "adjacent pairs must be evaluated and accepted immediately instead of pre-fitted globally");
        require(readText(root / "spatial_group_loose" /
                         "coverage_audit_pre_repair.json").find(
                    "\"repair_face_count\":0") != std::string::npos,
                "pre-canonicalization responsibility certificates must avoid safety repair");

        // Geometric adjacency is edge-to-edge, not vertex-to-vertex. The
        // vertical rectangle starts in the interior of the horizontal
        // rectangle's boundary edge, so the pair forms a T junction without a
        // shared OBJ vertex. A loose error limit must still give this adjacent
        // pair a stage-3 merge attempt.
        const auto t_junction = root / "surface_t_junction.obj";
        writeText(t_junction,
            "v 0 0 0\nv 2 0 0\nv 2 1 0\nv 0 1 0\n"
            "v 1 1 0\nv 1 2 0\nv 1 2 1\nv 1 1 1\n"
            "f 1 2 3\nf 1 3 4\nf 5 6 7\nf 5 7 8\n");
        const auto t_junction_stats =
            pqss_proxy_mesh::analyzePrimitiveMeshObj(
                t_junction, root / "surface_t_junction", group_loose);
        require(t_junction_stats.containment_validation_passed &&
                    t_junction_stats.primitive_count > 0,
                "an isolated non-coplanar T junction must remain conservatively enclosed after merging");
        require(readText(root / "surface_t_junction" /
                         "surface_merge_profile.json").find(
                    "\"segment_contact_adjacencies\":1") !=
                    std::string::npos,
                "T-junction adjacency must come from segment contact");
        require(readText(root / "surface_t_junction" /
                         "surface_merge_profile.json").find(
                    "\"connectivity_rejections\":") !=
                    std::string::npos,
                "surface merge profiling must retain connectivity diagnostics");

        auto group_strict = group_loose;
        group_strict.maximum_open_error_distance = 0.0;
        const auto group_strict_stats =
            pqss_proxy_mesh::analyzePrimitiveMeshObj(
                spatial_group, root / "spatial_group_strict", group_strict);
        require(group_strict_stats.merged_spatial_primitive_groups == 0 &&
                    group_strict_stats.primitive_count > 6,
                "exact error must also preserve the surface-only candidate set");

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

        // Thin outer plates can make their cavity-facing walls look like
        // redundant parallel skins.  The open channel between them is much
        // larger than the fill budget, so both inner walls must survive.
        const auto open_channel = root / "open_channel.obj";
        std::ostringstream open_channel_obj;
        next_vertex = 1;
        appendBoxObj(open_channel_obj, next_vertex,
                     0.0, 0.0, 0.0, 0.02, 1.0, 1.0);
        appendBoxObj(open_channel_obj, next_vertex,
                     0.98, 0.0, 0.0, 1.0, 1.0, 1.0);
        appendBoxObj(open_channel_obj, next_vertex,
                     0.02, 0.0, 0.0, 0.98, 1.0, 0.02);
        writeText(open_channel, open_channel_obj.str());
        auto preserve_channel = exact;
        preserve_channel.maximum_cavity_added_volume_ratio = 0.08;
        const auto channel_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
            open_channel, root / "open_channel", preserve_channel);
        require(channel_stats.containment_validation_passed &&
                    hasLargeFaceOnX(
                        root / "open_channel" / "primitives.obj", 0.02) &&
                    hasLargeFaceOnX(
                        root / "open_channel" / "primitives.obj", 0.98),
                "over-budget open-channel inner walls must not be occluded");

        std::cout << "primitive mesh staged-pipeline tests passed\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
