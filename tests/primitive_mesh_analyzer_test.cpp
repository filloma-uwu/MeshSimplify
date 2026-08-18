#include "pqss_proxy_mesh/primitive_mesh_analyzer.hpp"

#include <cmath>
#include <array>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <set>
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

std::set<std::uint64_t> positiveEnclosureGroups(const std::string& metadata)
{
    constexpr std::string_view key = "\"enclosure_group\":";
    std::set<std::uint64_t> groups;
    std::size_t offset = 0;
    while ((offset = metadata.find(key, offset)) != std::string::npos)
    {
        offset += key.size();
        const std::uint64_t group = std::stoull(metadata.substr(offset));
        if (group != 0) groups.insert(group);
    }
    return groups;
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

void appendExtrudedPolygonObj(
    std::ostringstream& obj, std::size_t& next_vertex,
    const std::vector<std::array<double, 2>>& polygon,
    const double z0, const double z1)
{
    const std::size_t base = next_vertex;
    for (const double z : {z0, z1})
        for (const auto& point : polygon)
            obj << "v " << point[0] << ' ' << point[1] << ' ' << z << '\n';
    for (std::size_t corner = 1; corner + 1 < polygon.size(); ++corner)
    {
        obj << "f " << base << ' ' << base + corner + 1 << ' '
            << base + corner << '\n';
        obj << "f " << base + polygon.size() << ' '
            << base + polygon.size() + corner << ' '
            << base + polygon.size() + corner + 1 << '\n';
    }
    for (std::size_t edge = 0; edge < polygon.size(); ++edge)
    {
        const std::size_t next = (edge + 1) % polygon.size();
        obj << "f " << base + edge << ' ' << base + next << ' '
            << base + polygon.size() + next << '\n'
            << "f " << base + edge << ' '
            << base + polygon.size() + next << ' '
            << base + polygon.size() + edge << '\n';
    }
    next_vertex += 2 * polygon.size();
}

void appendRectangularRingPrismObj(
    std::ostringstream& obj, std::size_t& next_vertex,
    const double outer_min, const double outer_max,
    const double inner_min, const double inner_max,
    const double z0, const double z1)
{
    const std::size_t base = next_vertex;
    const std::array<std::array<double, 2>, 8> rings{{
        {{outer_min, outer_min}}, {{outer_max, outer_min}},
        {{outer_max, outer_max}}, {{outer_min, outer_max}},
        {{inner_min, inner_min}}, {{inner_max, inner_min}},
        {{inner_max, inner_max}}, {{inner_min, inner_max}},
    }};
    for (const double z : {z0, z1})
        for (const auto& point : rings)
            obj << "v " << point[0] << ' ' << point[1] << ' ' << z << '\n';
    const auto quad = [&](const std::size_t a, const std::size_t b,
                          const std::size_t c, const std::size_t d)
    {
        obj << "f " << base + a << ' ' << base + b << ' ' << base + c << '\n'
            << "f " << base + a << ' ' << base + c << ' ' << base + d << '\n';
    };
    for (std::size_t edge = 0; edge < 4; ++edge)
    {
        const std::size_t next = (edge + 1) % 4;
        quad(edge, next, 4 + next, 4 + edge);
        quad(8 + edge, 12 + edge, 12 + next, 8 + next);
        quad(edge, 8 + edge, 8 + next, next);
        quad(4 + edge, 4 + next, 12 + next, 12 + edge);
    }
    next_vertex += 16;
}

void appendQuadGridObj(
    std::ostringstream& obj, std::size_t& next_vertex,
    const std::array<double, 3>& origin,
    const std::array<double, 3>& first,
    const std::array<double, 3>& second,
    const std::size_t subdivisions)
{
    const auto point = [&](const double first_parameter,
                           const double second_parameter)
    {
        std::array<double, 3> result{};
        for (int axis = 0; axis < 3; ++axis)
            result[axis] = origin[axis] +
                first_parameter * first[axis] +
                second_parameter * second[axis];
        return result;
    };
    for (std::size_t row = 0; row < subdivisions; ++row)
        for (std::size_t column = 0; column < subdivisions; ++column)
        {
            const double u0 = static_cast<double>(column) / subdivisions;
            const double u1 = static_cast<double>(column + 1) / subdivisions;
            const double v0 = static_cast<double>(row) / subdivisions;
            const double v1 = static_cast<double>(row + 1) / subdivisions;
            const std::array corners{
                point(u0, v0), point(u1, v0),
                point(u1, v1), point(u0, v1)};
            const std::size_t base = next_vertex;
            for (const auto& vertex : corners)
                obj << "v " << vertex[0] << ' ' << vertex[1] << ' '
                    << vertex[2] << '\n';
            obj << "f " << base << ' ' << base + 1 << ' ' << base + 2 << '\n'
                << "f " << base << ' ' << base + 2 << ' ' << base + 3 << '\n';
            next_vertex += 4;
        }
}

void appendSubdividedBoxObj(
    std::ostringstream& obj, std::size_t& next_vertex,
    const double x0, const double y0, const double z0,
    const double x1, const double y1, const double z1,
    const std::size_t subdivisions)
{
    appendQuadGridObj(obj, next_vertex, {x0, y0, z0},
                      {0, y1 - y0, 0}, {x1 - x0, 0, 0}, subdivisions);
    appendQuadGridObj(obj, next_vertex, {x0, y0, z1},
                      {x1 - x0, 0, 0}, {0, y1 - y0, 0}, subdivisions);
    appendQuadGridObj(obj, next_vertex, {x0, y0, z0},
                      {x1 - x0, 0, 0}, {0, 0, z1 - z0}, subdivisions);
    appendQuadGridObj(obj, next_vertex, {x0, y1, z0},
                      {0, 0, z1 - z0}, {x1 - x0, 0, 0}, subdivisions);
    appendQuadGridObj(obj, next_vertex, {x0, y0, z0},
                      {0, 0, z1 - z0}, {0, y1 - y0, 0}, subdivisions);
    appendQuadGridObj(obj, next_vertex, {x1, y0, z0},
                      {0, y1 - y0, 0}, {0, 0, z1 - z0}, subdivisions);
}

void appendNegativeYAttachmentObj(
    std::ostringstream& obj, std::size_t& next_vertex,
    const double x0, const double z0, const double y_front,
    const double x1, const double z1, const std::size_t subdivisions)
{
    appendQuadGridObj(obj, next_vertex, {x0, y_front, z0},
                      {x1 - x0, 0, 0}, {0, 0, z1 - z0}, subdivisions);
    appendQuadGridObj(obj, next_vertex, {x0, y_front, z0},
                      {0, 0, z1 - z0}, {0, -y_front, 0}, subdivisions);
    appendQuadGridObj(obj, next_vertex, {x1, y_front, z0},
                      {0, -y_front, 0}, {0, 0, z1 - z0}, subdivisions);
    appendQuadGridObj(obj, next_vertex, {x0, y_front, z0},
                      {0, -y_front, 0}, {x1 - x0, 0, 0}, subdivisions);
    appendQuadGridObj(obj, next_vertex, {x0, y_front, z1},
                      {x1 - x0, 0, 0}, {0, -y_front, 0}, subdivisions);
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

void appendCylinderArcObj(std::ostringstream& obj, std::size_t& next_vertex,
                          const std::size_t full_segments,
                          const std::size_t arc_segments,
                          const double radius, const double height)
{
    const std::size_t base = next_vertex;
    constexpr double pi = 3.14159265358979323846;
    for (std::size_t index = 0; index <= arc_segments; ++index)
    {
        const double angle = 2.0 * pi * static_cast<double>(index) /
                             static_cast<double>(full_segments);
        const double x = radius * std::cos(angle);
        const double y = radius * std::sin(angle);
        obj << "v " << x << ' ' << y << " 0\n"
            << "v " << x << ' ' << y << ' ' << height << '\n';
    }
    for (std::size_t index = 0; index < arc_segments; ++index)
    {
        const std::size_t lower = base + 2 * index;
        const std::size_t upper = lower + 1;
        const std::size_t next_lower = lower + 2;
        const std::size_t next_upper = next_lower + 1;
        obj << "f " << lower << ' ' << next_lower << ' ' << next_upper << '\n'
            << "f " << lower << ' ' << next_upper << ' ' << upper << '\n';
    }
    next_vertex += 2 * (arc_segments + 1);
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

bool hasHorizontalSurfaceAt(const std::filesystem::path& path,
                            const double x, const double y, const double z)
{
    std::ifstream stream(path);
    std::vector<std::array<double, 3>> vertices;
    std::string line;
    const auto signedArea = [](const std::array<double, 3>& first,
                               const std::array<double, 3>& second,
                               const double px, const double py)
    {
        return (second[0] - first[0]) * (py - first[1]) -
               (second[1] - first[1]) * (px - first[0]);
    };
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
        for (std::size_t corner = 1; corner + 1 < face.size(); ++corner)
        {
            const auto& first = vertices[face[0]];
            const auto& second = vertices[face[corner]];
            const auto& third = vertices[face[corner + 1]];
            if (std::abs(first[2] - z) > 1.0e-8 ||
                std::abs(second[2] - z) > 1.0e-8 ||
                std::abs(third[2] - z) > 1.0e-8)
                continue;
            const double first_side = signedArea(first, second, x, y);
            const double second_side = signedArea(second, third, x, y);
            const double third_side = signedArea(third, first, x, y);
            const bool nonnegative = first_side >= -1.0e-8 &&
                second_side >= -1.0e-8 && third_side >= -1.0e-8;
            const bool nonpositive = first_side <= 1.0e-8 &&
                second_side <= 1.0e-8 && third_side <= 1.0e-8;
            if (nonnegative || nonpositive) return true;
        }
    }
    return false;
}

bool hasAxisAlignedSurfaceAt(const std::filesystem::path& path,
                             const int normal_axis,
                             const std::array<double, 3>& point)
{
    std::ifstream stream(path);
    std::vector<std::array<double, 3>> vertices;
    std::string line;
    const int first_axis = (normal_axis + 1) % 3;
    const int second_axis = (normal_axis + 2) % 3;
    const auto signedArea = [&](const std::array<double, 3>& first,
                                const std::array<double, 3>& second)
    {
        return (second[first_axis] - first[first_axis]) *
                   (point[second_axis] - first[second_axis]) -
               (second[second_axis] - first[second_axis]) *
                   (point[first_axis] - first[first_axis]);
    };
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
        for (std::size_t corner = 1; corner + 1 < face.size(); ++corner)
        {
            const auto& first = vertices[face[0]];
            const auto& second = vertices[face[corner]];
            const auto& third = vertices[face[corner + 1]];
            if (std::abs(first[normal_axis] - point[normal_axis]) > 1.0e-8 ||
                std::abs(second[normal_axis] - point[normal_axis]) > 1.0e-8 ||
                std::abs(third[normal_axis] - point[normal_axis]) > 1.0e-8)
                continue;
            const double first_side = signedArea(first, second);
            const double second_side = signedArea(second, third);
            const double third_side = signedArea(third, first);
            if ((first_side >= -1.0e-8 && second_side >= -1.0e-8 &&
                 third_side >= -1.0e-8) ||
                (first_side <= 1.0e-8 && second_side <= 1.0e-8 &&
                 third_side <= 1.0e-8))
                return true;
        }
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
                "metadata must use the filled phase-1 error reference");

        // A closed genus-one component has no boundary hole. Until a concrete
        // added volume is accepted by the topology pass, the legacy planar
        // loop cleanup must leave both its outer shell and tunnel untouched.
        const auto closed_ring = root / "closed_ring_prism.obj";
        std::ostringstream closed_ring_obj;
        std::size_t closed_ring_vertex = 1;
        appendRectangularRingPrismObj(
            closed_ring_obj, closed_ring_vertex,
            0.0, 4.0, 1.0, 3.0, 0.0, 1.0);
        writeText(closed_ring, closed_ring_obj.str());
        const auto closed_ring_stats =
            pqss_proxy_mesh::analyzePrimitiveMeshObj(
                closed_ring, root / "closed_ring_prism", exact);
        const auto closed_ring_phase1 = root / "closed_ring_prism" /
            "phase1_hole_filled.obj";
        require(closed_ring_stats.containment_validation_passed &&
                    closed_ring_stats.filled_planar_holes == 0 &&
                    closed_ring_stats.removed_sealed_void_wall_primitives == 0 &&
                    hasAxisAlignedSurfaceAt(
                        closed_ring_phase1, 0, {0.0, 2.0, 0.5}) &&
                    hasAxisAlignedSurfaceAt(
                        closed_ring_phase1, 0, {1.0, 2.0, 0.5}),
                "legacy hole cleanup must not consume a closed manifold component");

        // Independent closed components are searched separately. Enclosure IDs
        // produced by those searches must be made model-global before later
        // occlusion passes use them as semantic ownership keys.
        const auto detached_boxes = root / "detached_closed_boxes.obj";
        std::ostringstream detached_boxes_obj;
        std::size_t detached_box_vertex = 1;
        appendBoxObj(detached_boxes_obj, detached_box_vertex,
                     0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        appendBoxObj(detached_boxes_obj, detached_box_vertex,
                     3.0, 0.0, 0.0, 4.0, 1.0, 1.0);
        writeText(detached_boxes, detached_boxes_obj.str());
        const auto detached_box_stats =
            pqss_proxy_mesh::analyzePrimitiveMeshObj(
                detached_boxes, root / "detached_closed_boxes", exact);
        const std::string detached_box_metadata = readText(
            root / "detached_closed_boxes" / "model.json");
        require(detached_box_stats.containment_validation_passed &&
                    positiveEnclosureGroups(detached_box_metadata).size() >= 2 &&
                    readText(root / "detached_closed_boxes" /
                             "top_down_surface_profile.json").find(
                        "\"closed_component_partitioned\":true") !=
                        std::string::npos,
                "closed component searches must use globally unique enclosure IDs");

        // Hole filling may remove an inner boundary loop, but it must not turn
        // an outer concavity into a free convex-hull rectangle before the
        // Hausdorff-controlled merge stage.  This U-shaped planar patch has a
        // rectangular convex hull and therefore catches that distinction.
        const auto concave_patch = root / "concave_planar_patch.obj";
        writeText(concave_patch,
            "v 0 0 0\nv 3 0 0\nv 3 3 0\nv 2 3 0\n"
            "v 2 1 0\nv 1 1 0\nv 1 3 0\nv 0 3 0\n"
            "f 1 2 5\nf 1 5 6\nf 1 6 7\nf 1 7 8\n"
            "f 2 3 4\nf 2 4 5\n");
        const auto concave_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
            concave_patch, root / "concave_planar_patch", exact);
        require(concave_stats.containment_validation_passed &&
                    concave_stats.primitive_count == 1 &&
                    concave_stats.proxy_triangles == 6,
                "phase 1 must preserve a concave outer planar boundary");
        require(readText(root / "concave_planar_patch" /
                         "phase1_hole_filled.obj").find("_polygon") !=
                    std::string::npos,
                "a concave phase-1 patch must remain a polygon, not a rectangle hull");

        // Partially overlapping coplanar source patches are one collision
        // surface. Coverage repair and final canonicalization must emit their
        // geometric union, not retain two stacked rectangles merely because
        // they came from disconnected OBJ components.
        const auto overlapping_patches = root / "overlapping_patches.obj";
        writeText(overlapping_patches,
            "v 0 0 0\nv 2 0 0\nv 2 1 0\nv 0 1 0\n"
            "v 1 0 0\nv 3 0 0\nv 3 1 0\nv 1 1 0\n"
            "f 1 2 3\nf 1 3 4\nf 5 6 7\nf 5 7 8\n");
        const auto overlapping_stats =
            pqss_proxy_mesh::analyzePrimitiveMeshObj(
                overlapping_patches, root / "overlapping_patches", exact);
        require(overlapping_stats.containment_validation_passed &&
                    overlapping_stats.primitive_count == 1 &&
                    overlapping_stats.proxy_triangles == 2,
                "overlapping coplanar patches must become one non-overlapping union");

        // Filling repeated cross-members closes the intervals between matching
        // boxes. Their horizontal top responsibility must remain, while the
        // bridge provides a horizontal surface over each former opening.
        const auto repeated_cross_members = root / "repeated_cross_members.obj";
        std::ostringstream cross_member_obj;
        std::size_t cross_member_vertex = 1;
        appendBoxObj(cross_member_obj, cross_member_vertex,
                     0.0, 0.0, 0.5, 0.5, 6.0, 1.0);
        appendBoxObj(cross_member_obj, cross_member_vertex,
                     3.5, 0.0, 0.5, 4.0, 6.0, 1.0);
        for (const double center : {1.2, 3.0, 4.8})
            appendBoxObj(cross_member_obj, cross_member_vertex,
                         0.0, center - 0.2, 0.0,
                         4.0, center + 0.2, 0.5);
        writeText(repeated_cross_members, cross_member_obj.str());
        auto cross_member_options = exact;
        cross_member_options.maximum_open_error_distance = 0.1;
        const auto cross_member_stats =
            pqss_proxy_mesh::analyzePrimitiveMeshObj(
                repeated_cross_members,
                root / "repeated_cross_members", cross_member_options);
        const auto cross_member_phase1 = root / "repeated_cross_members" /
            "phase1_hole_filled.obj";
        require(cross_member_stats.containment_validation_passed &&
                    hasHorizontalSurfaceAt(
                        cross_member_phase1, 2.0, 1.2, 0.5) &&
                    hasHorizontalSurfaceAt(
                        cross_member_phase1, 2.0, 3.0, 0.5) &&
                    hasHorizontalSurfaceAt(
                        cross_member_phase1, 2.0, 4.8, 0.5) &&
                    hasHorizontalSurfaceAt(
                        cross_member_phase1, 2.0, 2.0, 0.5),
                "phase-1 filling must close repeated-member gaps without losing their top surfaces");

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
        analytic.maximum_open_error_distance = 0.01;
        const auto cylinder_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
            cylinder_band, root / "cylinder_band", analytic);
        require(cylinder_stats.containment_validation_passed &&
                    cylinder_stats.cylindrical_band_count == 1 &&
                    cylinder_stats.primitive_count == 1 &&
                    cylinder_stats.proxy_triangles == 48,
                "a complete cylinder side must remain one certified surface");

        // An exposed cylindrical arc is still one analytic surface. It must
        // not be completed into a closed cylinder or replaced by rectangular
        // box faces merely because caps or neighboring walls are absent.
        const auto cylinder_arc = root / "cylinder_arc.obj";
        std::ostringstream cylinder_arc_obj;
        std::size_t cylinder_arc_vertex = 1;
        appendCylinderArcObj(
            cylinder_arc_obj, cylinder_arc_vertex, 96, 24, 1.0, 2.0);
        writeText(cylinder_arc, cylinder_arc_obj.str());
        const auto cylinder_arc_stats =
            pqss_proxy_mesh::analyzePrimitiveMeshObj(
                cylinder_arc, root / "cylinder_arc", analytic);
        require(cylinder_arc_stats.containment_validation_passed &&
                    cylinder_arc_stats.cylindrical_band_count == 1 &&
                    cylinder_arc_stats.primitive_count == 1 &&
                    cylinder_arc_stats.proxy_triangles < 48,
                "a partial cylinder must remain one trimmed analytic surface");

        // Equal corner radii do not make a circular cross-section.  This is a
        // rectangular prism with short corner chamfers: the long straight
        // edges create nonuniform angular support and must prevent promotion to
        // disks or a cylindrical band.
        const auto chamfered_prism = root / "chamfered_rectangular_prism.obj";
        std::ostringstream chamfered_obj;
        const std::array<std::array<double, 2>, 8> chamfered_ring{{
            {{-230.0, -220.0}}, {{-210.0, -240.0}},
            {{ 210.0, -240.0}}, {{ 230.0, -220.0}},
            {{ 230.0,  220.0}}, {{ 210.0,  240.0}},
            {{-210.0,  240.0}}, {{-230.0,  220.0}},
        }};
        for (const double x : {0.0, 10.0})
            for (const auto& point : chamfered_ring)
                chamfered_obj << "v " << x << ' ' << point[0] << ' '
                              << point[1] << '\n';
        for (std::size_t corner = 1; corner + 1 < chamfered_ring.size(); ++corner)
            chamfered_obj << "f 1 " << corner + 1 << ' ' << corner + 2 << '\n'
                          << "f 9 " << 9 + corner + 1 << ' '
                          << 9 + corner << '\n';
        for (std::size_t edge = 0; edge < chamfered_ring.size(); ++edge)
        {
            const std::size_t next = (edge + 1) % chamfered_ring.size();
            chamfered_obj << "f " << edge + 1 << ' ' << next + 1 << ' '
                          << next + 9 << '\n'
                          << "f " << edge + 1 << ' ' << next + 9 << ' '
                          << edge + 9 << '\n';
        }
        writeText(chamfered_prism, chamfered_obj.str());
        const auto chamfered_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
            chamfered_prism, root / "chamfered_rectangular_prism", analytic);
        require(chamfered_stats.containment_validation_passed &&
                    chamfered_stats.disk_count == 0 &&
                    chamfered_stats.cylindrical_band_count == 0 &&
                    chamfered_stats.conical_band_count == 0 &&
                    chamfered_stats.primitive_count == 10 &&
                    chamfered_stats.proxy_triangles == 28,
                "a chamfered rectangular prism must not be promoted to a round surface");

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

        // A loose directed-distance limit does not make a concave L extrusion
        // a member of the box surface family. Its inset side faces do not lie
        // on any of the six OBB support planes.
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
                "a loose error may still simplify the L shell through local adjacent merges");
        const std::string top_down_profile = readText(
            root / "spatial_group_loose" /
                "top_down_surface_profile.json");
        require(top_down_profile.find("\"surface_candidates_only\":true") !=
                    std::string::npos &&
                    top_down_profile.find("\"skipped\":true") ==
                    std::string::npos &&
                    top_down_profile.find("\"accepted_box_surface_sets\":0") ==
                        std::string::npos,
                "top-down candidate search must still accept unprotected low-error boxes");
        const std::string adjacent_envelope_profile = readText(
            root / "spatial_group_loose" /
                "adjacent_envelope_group_profile.json");
        require(adjacent_envelope_profile.find("\"complete\":true") !=
                    std::string::npos &&
                    adjacent_envelope_profile.find(
                        "\"quickhull_kernel_seconds\":0") !=
                    std::string::npos,
                "top-down output must still reach the adjacent-envelope fixed point");
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

        // The top face of this deep narrow slot is cheap in directed
        // proxy-to-source distance, even with a very loose limit. The slot
        // bottom and side walls are nevertheless real boundary surfaces and
        // disqualify the enclosing OBB as a box-family candidate.
        const auto deep_slot = root / "deep_slot_prism.obj";
        std::ostringstream deep_slot_obj;
        std::size_t deep_slot_next_vertex = 1;
        appendExtrudedPolygonObj(
            deep_slot_obj, deep_slot_next_vertex,
            {{{0.0, 0.0}}, {{300.0, 0.0}}, {{300.0, 500.0}},
             {{155.0, 500.0}}, {{155.0, 200.0}}, {{145.0, 200.0}},
             {{145.0, 500.0}}, {{0.0, 500.0}}},
            0.0, 300.0);
        writeText(deep_slot, deep_slot_obj.str());
        auto deep_slot_policy = exact;
        deep_slot_policy.maximum_cavity_added_volume_ratio = 0.0;
        deep_slot_policy.maximum_open_error_distance = 100.0;
        const auto deep_slot_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
            deep_slot, root / "deep_slot_prism", deep_slot_policy);
        const std::string deep_slot_top_down_profile = readText(
            root / "deep_slot_prism" / "top_down_surface_profile.json");
        const std::string deep_slot_adjacent_profile = readText(
            root / "deep_slot_prism" /
                "adjacent_envelope_group_profile.json");
        require(deep_slot_stats.containment_validation_passed &&
                deep_slot_top_down_profile.find(
                    "\"boundary_responsibility_rejections\":0") ==
                    std::string::npos &&
                deep_slot_top_down_profile.find(
                    "\"unowned_intrusion_rejections\":0") ==
                    std::string::npos &&
                deep_slot_adjacent_profile.find(
                    "\"shared_responsibility_rejections\":0") ==
                    std::string::npos &&
                hasAxisAlignedSurfaceAt(
                    root / "deep_slot_prism" / "primitives.obj", 1,
                    {150.0, 200.0, 150.0}),
                "neither top-down nor pairwise enclosure growth may consume a deeper structured responsibility");

        // A detached surface can overlap a candidate hull's AABB while lying
        // wholly outside one of its oblique supporting planes.  The stage-3
        // intrusion index must prove that separation in its hierarchy and
        // leave the detached responsibility untouched.
        const auto spatial_group_with_external =
            root / "spatial_group_with_external.obj";
        writeText(spatial_group_with_external,
            readText(spatial_group) +
            "v 1.85 2.80 0.40\n"
            "v 1.95 2.80 0.40\n"
            "v 1.90 2.90 0.40\n"
            "f 13 14 15\n");
        const auto spatial_external_stats =
            pqss_proxy_mesh::analyzePrimitiveMeshObj(
                spatial_group_with_external,
                root / "spatial_group_with_external", group_loose);
        require(spatial_external_stats.containment_validation_passed &&
                    spatial_external_stats.coverage_failed_source_faces == 0,
                "hierarchical halfspace pruning must preserve detached source responsibility");

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
        auto t_junction_policy = group_loose;
        t_junction_policy.maximum_open_error_distance = 0.0;
        const auto t_junction_stats =
            pqss_proxy_mesh::analyzePrimitiveMeshObj(
                t_junction, root / "surface_t_junction", t_junction_policy);
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

        // Collision workload is a post-generation diagnostic, never an
        // acceptance condition. These six nearly coplanar patches meet along
        // offset T junctions. Whether top-down or adjacent selection handles
        // them first, neither candidate stage may reject geometry because of
        // its intermediate triangle count.
        const auto workload_independent_merge =
            root / "workload_independent_merge.obj";
        writeText(workload_independent_merge,
            "v 0 0 -1\n"
            "v 0 0 0\nv 1 0 0\nv 1 1 0\nv 0 1 0\n"
            "v 1 0 0.00000001\nv 2 0 0.00000001\n"
            "v 2 1 0.00000001\nv 1 1 0.00000001\n"
            "v 0 1 0\nv 2 1 0\nv 2 2 0.00000002\n"
            "v 0 2 0.00000002\n"
            "v 0 0 0\nv 2 0 0\nv 2 -1 0.00000002\n"
            "v 0 -1 0.00000002\n"
            "v 0 0 0\nv 0 1 0\nv -1 1 0.00000002\n"
            "v -1 0 0.00000002\n"
            "v 2 0 0\nv 2 1 0\nv 3 1 0.00000002\n"
            "v 3 0 0.00000002\n"
            "f 2 3 4\nf 2 4 5\nf 6 7 8\nf 6 8 9\n"
            "f 10 11 12\nf 10 12 13\nf 14 16 15\nf 14 17 16\n"
            "f 18 19 20\nf 18 20 21\nf 22 24 23\nf 22 25 24\n");
        const auto workload_independent_stats =
            pqss_proxy_mesh::analyzePrimitiveMeshObj(
                workload_independent_merge,
                root / "workload_independent_merge", group_loose);
        const std::string workload_independent_profile = readText(
            root / "workload_independent_merge" /
                "surface_merge_profile.json");
        const std::string workload_independent_top_down_profile = readText(
            root / "workload_independent_merge" /
                "top_down_surface_profile.json");
        require(workload_independent_stats.containment_validation_passed &&
                    (workload_independent_stats.merged_local_planar_primitives > 0 ||
                     workload_independent_stats.merged_spatial_primitive_groups > 0) &&
                    workload_independent_profile.find(
                        "\"workload_rejections\":0") !=
                        std::string::npos &&
                    workload_independent_top_down_profile.find(
                        "\"workload_rejections\":0") !=
                        std::string::npos,
                "an error-bounded merge must not be rejected or ranked by triangle workload");

        auto group_strict = group_loose;
        group_strict.maximum_open_error_distance = 0.0;
        const auto group_strict_stats =
            pqss_proxy_mesh::analyzePrimitiveMeshObj(
                spatial_group, root / "spatial_group_strict", group_strict);
        require(group_strict_stats.merged_spatial_primitive_groups == 0 &&
                    group_strict_stats.primitive_count > 6,
                "exact error must also preserve the surface-only candidate set");

        // Independently indexed solids may share an interface. That coincident
        // pair is internal to their union and must disappear without removing
        // either solid's exterior boundary.
        const auto touching_components = root / "touching_closed_components.obj";
        std::ostringstream touching_components_obj;
        std::size_t touching_vertex = 1;
        appendBoxObj(touching_components_obj, touching_vertex,
                     0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        appendBoxObj(touching_components_obj, touching_vertex,
                     1.0, 0.0, 0.0, 2.0, 1.0, 1.0);
        writeText(touching_components, touching_components_obj.str());
        const auto touching_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
            touching_components, root / "touching_closed_components", exact);
        require(touching_stats.containment_validation_passed &&
                    !hasAxisAlignedSurfaceAt(
                        root / "touching_closed_components" / "primitives.obj",
                        0, {1.0, 0.5, 0.5}),
                "a coincident interface between closed components must not survive");

        // A long, thin attachment may merge into the adjacent body when the
        // user-directed distance allows it. The result remains governed by
        // conservative coverage and the same global error limit.
        const auto nonlocal_attachment = root / "nonlocal_attachment.obj";
        std::ostringstream nonlocal_attachment_obj;
        std::size_t next_vertex = 1;
        appendBoxObj(nonlocal_attachment_obj, next_vertex,
                     0, 0, 0, 100, 100, 1);
        appendSubdividedBoxObj(nonlocal_attachment_obj, next_vertex,
                               35, 49.5, 1.1, 65, 50.5, 2.1, 2);
        writeText(nonlocal_attachment, nonlocal_attachment_obj.str());
        auto local_attachment_policy = exact;
        local_attachment_policy.maximum_open_error_distance = 100.0;
        local_attachment_policy.maximum_cavity_added_volume_ratio = 0.0;
        const auto local_attachment_stats =
            pqss_proxy_mesh::analyzePrimitiveMeshObj(
                nonlocal_attachment, root / "nonlocal_attachment_limited",
                local_attachment_policy);
        require(local_attachment_stats.containment_validation_passed &&
                    local_attachment_stats.open_max_distance <=
                        local_attachment_policy.maximum_open_error_distance +
                            1.0e-9,
                "an accepted attachment merge must remain conservative and error-bounded");

        // A typed multi-face attachment is one atomic enclosure candidate.
        // Pairwise convex fallback may still use it as an intrusion witness,
        // but must not consume one attachment at a time and grow a broad plane
        // across the empty support region between distinct attachments.
        const auto repeated_attachments = root / "repeated_attachments.obj";
        std::ostringstream repeated_attachment_obj;
        next_vertex = 1;
        appendBoxObj(repeated_attachment_obj, next_vertex,
                     0.0, 0.0, 0.0, 100.0, 10.0, 100.0);
        for (const auto& center : std::array<std::array<double, 2>, 4>{{
                 {{20.0, 20.0}}, {{80.0, 20.0}},
                 {{20.0, 80.0}}, {{80.0, 80.0}}}})
            appendNegativeYAttachmentObj(
                repeated_attachment_obj, next_vertex,
                center[0] - 5.0, center[1] - 5.0, -5.0,
                center[0] + 5.0, center[1] + 5.0, 3);
        repeated_attachment_obj
            << "v 45 30 45\nv 55 30 45\nv 50 30 55\n"
            << "f " << next_vertex << ' ' << next_vertex + 1 << ' '
            << next_vertex + 2 << '\n';
        next_vertex += 3;
        writeText(repeated_attachments, repeated_attachment_obj.str());
        auto repeated_attachment_policy = local_attachment_policy;
        repeated_attachment_policy.maximum_open_error_distance = 1.0;
        const auto repeated_attachment_stats =
            pqss_proxy_mesh::analyzePrimitiveMeshObj(
                repeated_attachments, root / "repeated_attachments",
                repeated_attachment_policy);
        const auto repeated_attachment_proxy =
            root / "repeated_attachments" / "primitives.obj";
        const std::string repeated_support_profile = readText(
            root / "repeated_attachments" /
                "support_protrusion_merge_profile.json");
        const std::string repeated_top_down_profile = readText(
            root / "repeated_attachments" /
                "top_down_surface_profile.json");
        require(repeated_attachment_stats.containment_validation_passed &&
                    repeated_support_profile.find(
                        "\"skipped\":true") != std::string::npos &&
                    repeated_support_profile.find(
                        "\"reason\":\"model_level_candidate_search\"") !=
                        std::string::npos &&
                    repeated_top_down_profile.find(
                        "\"complete\":true") != std::string::npos,
                "attachment recognition must not greedily consume responsibility before model-level selection");
        for (const auto& center : std::array<std::array<double, 2>, 4>{{
                 {{20.0, 20.0}}, {{80.0, 20.0}},
                 {{20.0, 80.0}}, {{80.0, 80.0}}}})
            require(hasAxisAlignedSurfaceAt(
                        repeated_attachment_proxy, 1,
                        {center[0], -5.0, center[1]}),
                    "each certified attachment must retain its outer surface");
        require(!hasAxisAlignedSurfaceAt(
                    repeated_attachment_proxy, 1, {50.0, -5.0, 50.0}),
                "atomic attachments must not grow a plane across exterior gaps");

        // Ordinary exterior gaps are not topological tunnels, regardless of
        // their width or the configured simplification distance.
        const auto narrow_gaps = root / "narrow_component_gaps.obj";
        std::ostringstream narrow_gap_obj;
        next_vertex = 1;
        appendBoxObj(narrow_gap_obj, next_vertex, 0.0, 0, 0, 1.0, 1, 1);
        appendBoxObj(narrow_gap_obj, next_vertex, 1.1, 0, 0, 2.1, 1, 1);
        appendBoxObj(narrow_gap_obj, next_vertex, 2.2, 0, 0, 3.2, 1, 1);
        writeText(narrow_gaps, narrow_gap_obj.str());
        auto gap_fill = exact;
        gap_fill.maximum_open_error_distance = 0.0;
        gap_fill.maximum_cavity_added_volume_ratio = 0.04;
        const auto gap_fill_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
            narrow_gaps, root / "narrow_component_gaps", gap_fill);
        require(gap_fill_stats.filled_intercomponent_gaps == 0 &&
                    !hasHorizontalSurfaceAt(
                        root / "narrow_component_gaps" /
                            "phase1_hole_filled.obj",
                        1.05, 0.5, 1.0),
                "an exterior component gap must remain connected to infinity");

        // Gap filling is not controlled by the simplification-error limit or
        // a local volume allowance.
        const auto large_gap = root / "large_component_gap.obj";
        std::ostringstream large_gap_obj;
        next_vertex = 1;
        appendBoxObj(large_gap_obj, next_vertex, 0, 0, 0, 1, 1, 1);
        appendBoxObj(large_gap_obj, next_vertex, 3, 0, 0, 4, 1, 1);
        writeText(large_gap, large_gap_obj.str());
        const auto large_gap_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
            large_gap, root / "large_component_gap", gap_fill);
        require(large_gap_stats.filled_intercomponent_gaps == 0 &&
                    !hasHorizontalSurfaceAt(
                        root / "large_component_gap" /
                            "phase1_hole_filled.obj",
                        2.0, 0.5, 1.0),
                "topology, not gap size, must reject an exterior gap");

        // Filling a recognized channel removes its two internal side walls;
        // their responsibility transfers to the retained fill surface.
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
        preserve_channel.maximum_open_error_distance = 1.0e-7;
        preserve_channel.maximum_cavity_added_volume_ratio = 0.08;
        const auto channel_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
            open_channel, root / "open_channel", preserve_channel);
        require(channel_stats.containment_validation_passed &&
                    channel_stats.filled_intercomponent_gaps == 0,
                "an open U-shaped channel must not be filled");

        // Two longitudinal solids and three transverse solids form two handles.
        // Filling is accepted only when the 3D cubical-complex certificate
        // reduces beta1 while preserving beta0 and beta2.
        const auto ladder_tunnels = root / "ladder_tunnels.obj";
        std::ostringstream ladder_obj;
        next_vertex = 1;
        appendBoxObj(ladder_obj, next_vertex, 0.0, 0.0, 0.0,
                     0.2, 3.0, 1.0);
        appendBoxObj(ladder_obj, next_vertex, 0.8, 0.0, 0.0,
                     1.0, 3.0, 1.0);
        appendBoxObj(ladder_obj, next_vertex, 0.0, 0.0, 0.0,
                     1.0, 0.2, 1.0);
        appendBoxObj(ladder_obj, next_vertex, 0.0, 1.4, 0.0,
                     1.0, 1.6, 1.0);
        appendBoxObj(ladder_obj, next_vertex, 0.0, 2.8, 0.0,
                     1.0, 3.0, 1.0);
        writeText(ladder_tunnels, ladder_obj.str());
        const auto ladder_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
            ladder_tunnels, root / "ladder_tunnels", gap_fill);
        const auto ladder_phase1 = root / "ladder_tunnels" /
            "phase1_hole_filled.obj";
        require(ladder_stats.containment_validation_passed &&
                    ladder_stats.filled_intercomponent_gaps == 2,
                "each independent ladder handle must receive one certified fill");
        require(hasAxisAlignedSurfaceAt(
                    ladder_phase1, 2, {0.5, 1.5, 0.0}) &&
                    hasAxisAlignedSurfaceAt(
                    ladder_phase1, 2, {0.5, 1.5, 1.0}),
                "filling a handle must retain the transverse solid's bottom and top");
        require(!hasAxisAlignedSurfaceAt(
                    ladder_phase1, 1, {0.5, 1.4, 0.5}) &&
                    !hasAxisAlignedSurfaceAt(
                    ladder_phase1, 1, {0.5, 1.6, 0.5}),
                "filled-volume coverage must remove both internal side walls");
        const std::string topology_profile = readText(
            root / "ladder_tunnels" / "intercomponent_gap_profile.json");
        require(topology_profile.find("\"initial_betti\":[1,2,0]") !=
                    std::string::npos &&
                topology_profile.find("\"final_betti\":[1,0,0]") !=
                    std::string::npos,
                "the fill profile must record the certified Betti-number change");

        std::cout << "primitive mesh staged-pipeline tests passed\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
