#include "pqss_proxy_mesh/primitive_mesh_analyzer.hpp"

#include <algorithm>
#include <cmath>
#include <array>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace
{

int check(const bool condition, const char* message)
{
    if (!condition) std::cerr << message << '\n';
    return condition ? 0 : 1;
}

void writeText(const std::filesystem::path& path, const std::string& text)
{
    std::ofstream stream(path);
    stream << text;
}

std::string readText(const std::filesystem::path& path)
{
    std::ifstream stream(path);
    return {std::istreambuf_iterator<char>(stream), std::istreambuf_iterator<char>()};
}

bool objIsClosedTriangleSoup(const std::filesystem::path& path)
{
    std::ifstream stream(path);
    std::vector<std::array<double, 3>> vertices;
    std::vector<std::array<std::size_t, 3>> faces;
    std::string line;
    while (std::getline(stream, line))
    {
        std::istringstream fields(line);
        std::string kind;
        fields >> kind;
        if (kind == "v")
        {
            std::array<double, 3> vertex{};
            fields >> vertex[0] >> vertex[1] >> vertex[2];
            vertices.push_back(vertex);
            continue;
        }
        if (kind != "f") continue;
        std::array<std::size_t, 3> face{};
        for (std::size_t& index : face)
        {
            std::string token;
            fields >> token;
            index = std::stoull(token.substr(0, token.find('/'))) - 1;
        }
        faces.push_back(face);
    }
    std::map<std::array<long long, 3>, std::size_t> welded_vertices;
    std::vector<std::size_t> welded(vertices.size());
    for (std::size_t index = 0; index < vertices.size(); ++index)
    {
        std::array<long long, 3> key{};
        for (int axis = 0; axis < 3; ++axis)
            key[axis] = std::llround(vertices[index][axis] * 1.0e9);
        const auto [iterator, inserted] =
            welded_vertices.emplace(key, welded_vertices.size());
        (void)inserted;
        welded[index] = iterator->second;
    }
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> edge_counts;
    for (const auto& face : faces)
    {
        for (int edge = 0; edge < 3; ++edge)
        {
            const auto first = welded[face[edge]];
            const auto second = welded[face[(edge + 1) % 3]];
            ++edge_counts[std::minmax(first, second)];
        }
    }
    return !edge_counts.empty() &&
        std::all_of(edge_counts.begin(), edge_counts.end(),
                    [](const auto& edge) { return edge.second == 2; });
}

std::string separatedBoxesObj(const double gap)
{
    std::ostringstream stream;
    constexpr std::array<std::array<int, 3>, 8> corners{{
        {{0, 0, 0}}, {{1, 0, 0}}, {{1, 1, 0}}, {{0, 1, 0}},
        {{0, 0, 1}}, {{1, 0, 1}}, {{1, 1, 1}}, {{0, 1, 1}},
    }};
    constexpr std::array<std::array<int, 3>, 12> faces{{
        {{0, 1, 2}}, {{0, 2, 3}}, {{4, 6, 5}}, {{4, 7, 6}},
        {{0, 4, 5}}, {{0, 5, 1}}, {{1, 5, 6}}, {{1, 6, 2}},
        {{2, 6, 7}}, {{2, 7, 3}}, {{3, 7, 4}}, {{3, 4, 0}},
    }};
    for (int box = 0; box < 4; ++box)
    {
        const double offset = box * (1.0 + gap);
        for (const auto& corner : corners)
            stream << "v " << offset + corner[0] << ' ' << corner[1] << ' '
                   << corner[2] << '\n';
    }
    for (int box = 0; box < 4; ++box)
    {
        const int offset = box * 8 + 1;
        for (const auto& face : faces)
            stream << "f " << offset + face[0] << ' ' << offset + face[1]
                   << ' ' << offset + face[2] << '\n';
    }
    return stream.str();
}

std::string intersectingBoxesObj()
{
    std::ostringstream stream;
    constexpr std::array<std::array<int, 3>, 12> faces{{
        {{0, 1, 2}}, {{0, 2, 3}}, {{4, 6, 5}}, {{4, 7, 6}},
        {{0, 4, 5}}, {{0, 5, 1}}, {{1, 5, 6}}, {{1, 6, 2}},
        {{2, 6, 7}}, {{2, 7, 3}}, {{3, 7, 4}}, {{3, 4, 0}},
    }};
    const auto appendBox = [&](const std::array<double, 3>& lower,
                               const std::array<double, 3>& upper,
                               const int vertex_offset)
    {
        constexpr std::array<std::array<int, 3>, 8> corners{{
            {{0, 0, 0}}, {{1, 0, 0}}, {{1, 1, 0}}, {{0, 1, 0}},
            {{0, 0, 1}}, {{1, 0, 1}}, {{1, 1, 1}}, {{0, 1, 1}},
        }};
        for (const auto& corner : corners)
            stream << "v " << (corner[0] ? upper[0] : lower[0]) << ' '
                   << (corner[1] ? upper[1] : lower[1]) << ' '
                   << (corner[2] ? upper[2] : lower[2]) << '\n';
        for (const auto& face : faces)
            stream << "f " << vertex_offset + face[0] + 1 << ' '
                   << vertex_offset + face[1] + 1 << ' '
                   << vertex_offset + face[2] + 1 << '\n';
    };
    appendBox({0.0, 0.0, 0.0}, {2.0, 1.0, 1.0}, 0);
    appendBox({1.0, 0.25, -0.5}, {3.0, 0.75, 1.5}, 8);
    return stream.str();
}

std::string partiallyTouchingBoxesObj()
{
    std::ostringstream stream;
    constexpr std::array<std::array<int, 3>, 12> faces{{
        {{0, 1, 2}}, {{0, 2, 3}}, {{4, 6, 5}}, {{4, 7, 6}},
        {{0, 4, 5}}, {{0, 5, 1}}, {{1, 5, 6}}, {{1, 6, 2}},
        {{2, 6, 7}}, {{2, 7, 3}}, {{3, 7, 4}}, {{3, 4, 0}},
    }};
    const auto appendBox = [&](const std::array<double, 3>& lower,
                               const std::array<double, 3>& upper,
                               const int vertex_offset)
    {
        constexpr std::array<std::array<int, 3>, 8> corners{{
            {{0, 0, 0}}, {{1, 0, 0}}, {{1, 1, 0}}, {{0, 1, 0}},
            {{0, 0, 1}}, {{1, 0, 1}}, {{1, 1, 1}}, {{0, 1, 1}},
        }};
        for (const auto& corner : corners)
            stream << "v " << (corner[0] ? upper[0] : lower[0]) << ' '
                   << (corner[1] ? upper[1] : lower[1]) << ' '
                   << (corner[2] ? upper[2] : lower[2]) << '\n';
        for (const auto& face : faces)
            stream << "f " << vertex_offset + face[0] + 1 << ' '
                   << vertex_offset + face[1] + 1 << ' '
                   << vertex_offset + face[2] + 1 << '\n';
    };
    appendBox({0.0, 0.0, 0.0}, {2.0, 2.0, 2.0}, 0);
    appendBox({2.0, 0.5, 0.5}, {3.0, 1.5, 1.5}, 8);
    return stream.str();
}

std::string repeatObjFaces(const std::string& obj, const std::size_t repeats)
{
    std::istringstream input(obj);
    std::ostringstream vertices;
    std::vector<std::string> faces;
    std::string line;
    while (std::getline(input, line))
    {
        if (line.starts_with("f ")) faces.push_back(line);
        else vertices << line << '\n';
    }
    for (std::size_t repeat = 0; repeat < repeats; ++repeat)
        for (const std::string& face : faces) vertices << face << '\n';
    return vertices.str();
}

std::string denselyNotchedPlanarObj(const std::size_t notch_count)
{
    std::vector<std::array<double, 2>> boundary;
    boundary.reserve(3 * notch_count + 4);
    boundary.push_back({0.0, 0.0});
    boundary.push_back({static_cast<double>(notch_count), 0.0});
    boundary.push_back({static_cast<double>(notch_count), 1.0});
    for (std::size_t reverse = notch_count; reverse-- > 0;)
    {
        const double x = static_cast<double>(reverse);
        boundary.push_back({x + 2.0 / 3.0, 1.0});
        boundary.push_back({x + 0.5, 0.995});
        boundary.push_back({x + 1.0 / 3.0, 1.0});
    }
    boundary.push_back({0.0, 1.0});

    std::ostringstream stream;
    for (const auto& point : boundary)
        stream << "v " << point[0] << ' ' << point[1] << " 0\n";
    for (std::size_t vertex = 1; vertex + 1 < boundary.size(); ++vertex)
        stream << "f 1 " << vertex + 1 << ' ' << vertex + 2 << '\n';
    return stream.str();
}

bool triangleMeshCoversPlanarPoint(const std::string& obj,
                                   const double x,
                                   const double y,
                                   const double z)
{
    struct Point { double x = 0.0; double y = 0.0; double z = 0.0; };
    std::vector<Point> vertices;
    std::istringstream stream(obj);
    std::string line;
    const auto side = [](const Point& first, const Point& second,
                         const double px, const double py)
    {
        return (second.x - first.x) * (py - first.y) -
               (second.y - first.y) * (px - first.x);
    };
    while (std::getline(stream, line))
    {
        std::istringstream fields(line);
        std::string kind;
        fields >> kind;
        if (kind == "v")
        {
            Point point;
            fields >> point.x >> point.y >> point.z;
            vertices.push_back(point);
            continue;
        }
        if (kind != "f") continue;
        std::array<std::size_t, 3> indices{};
        for (std::size_t& index : indices)
        {
            std::string token;
            fields >> token;
            index = static_cast<std::size_t>(std::stoul(token.substr(0, token.find('/'))) - 1);
        }
        const Point& first = vertices[indices[0]];
        const Point& second = vertices[indices[1]];
        const Point& third = vertices[indices[2]];
        if (std::abs(first.z - z) > 1.0e-9 ||
            std::abs(second.z - z) > 1.0e-9 ||
            std::abs(third.z - z) > 1.0e-9) continue;
        const double a = side(first, second, x, y);
        const double b = side(second, third, x, y);
        const double c = side(third, first, x, y);
        if ((a >= -1.0e-9 && b >= -1.0e-9 && c >= -1.0e-9) ||
            (a <= 1.0e-9 && b <= 1.0e-9 && c <= 1.0e-9)) return true;
    }
    return false;
}

bool triangleMeshContainsPoint(const std::string& obj,
                               const std::array<double, 3>& query)
{
    struct Point
    {
        double value[3]{};
    };
    const auto subtract = [](const Point& first, const Point& second)
    {
        return Point{{first.value[0] - second.value[0],
                      first.value[1] - second.value[1],
                      first.value[2] - second.value[2]}};
    };
    const auto dot = [](const Point& first, const Point& second)
    {
        return first.value[0] * second.value[0] +
               first.value[1] * second.value[1] +
               first.value[2] * second.value[2];
    };
    const auto cross = [](const Point& first, const Point& second)
    {
        return Point{{first.value[1] * second.value[2] -
                          first.value[2] * second.value[1],
                      first.value[2] * second.value[0] -
                          first.value[0] * second.value[2],
                      first.value[0] * second.value[1] -
                          first.value[1] * second.value[0]}};
    };
    std::vector<Point> vertices;
    std::vector<std::array<std::size_t, 3>> faces;
    std::istringstream stream(obj);
    std::string line;
    while (std::getline(stream, line))
    {
        std::istringstream fields(line);
        std::string kind;
        fields >> kind;
        if (kind == "v")
        {
            Point point;
            fields >> point.value[0] >> point.value[1] >> point.value[2];
            vertices.push_back(point);
        }
        else if (kind == "f")
        {
            std::array<std::size_t, 3> face{};
            for (auto& index : face)
            {
                std::string token;
                fields >> token;
                index = std::stoull(token.substr(0, token.find('/'))) - 1;
            }
            faces.push_back(face);
        }
    }
    const Point point{{query[0], query[1], query[2]}};
    for (const auto& face : faces)
    {
        const Point a = vertices[face[0]];
        const Point first = subtract(vertices[face[1]], a);
        const Point second = subtract(vertices[face[2]], a);
        const Point offset = subtract(point, a);
        const Point normal = cross(first, second);
        const double normal_squared = dot(normal, normal);
        if (normal_squared <= 1.0e-24 ||
            std::abs(dot(normal, offset)) > 1.0e-9 * std::sqrt(normal_squared))
            continue;
        const double first_squared = dot(first, first);
        const double second_squared = dot(second, second);
        const double mixed = dot(first, second);
        const double first_offset = dot(first, offset);
        const double second_offset = dot(second, offset);
        const double denominator = first_squared * second_squared - mixed * mixed;
        if (std::abs(denominator) <= 1.0e-24) continue;
        const double u = (second_squared * first_offset - mixed * second_offset) /
                         denominator;
        const double v = (first_squared * second_offset - mixed * first_offset) /
                         denominator;
        if (u >= -1.0e-9 && v >= -1.0e-9 && u + v <= 1.0 + 1.0e-9)
            return true;
    }
    return false;
}

} // namespace

int main()
{
    const double coarsest_threshold =
        pqss_proxy_mesh::analysisStrengthToAllowedExcessRatio(0.0);
    if (check(std::isinf(coarsest_threshold) && coarsest_threshold > 0.0,
              "strength zero must map to positive infinity without a special analysis path")) return 1;
    if (check(std::abs(pqss_proxy_mesh::analysisStrengthToAllowedExcessRatio(0.5) - 0.01) < 1.0e-12,
              "strength 0.5 must map to 0.01")) return 1;
    if (check(std::abs(pqss_proxy_mesh::analysisStrengthToAllowedExcessRatio(1.0) - 1.0e-5) < 1.0e-15,
              "strength one must map to 1e-5")) return 1;

    const auto root = std::filesystem::temp_directory_path() / "pqss_primitive_mesh_analyzer_test";
    std::filesystem::remove_all(root);
    std::filesystem::create_directories(root);

    const auto box = root / "2.obj";
    writeText(box,
        "v 0 0 0\nv 1 0 0\nv 1 1 0\nv 0 1 0\n"
        "v 0 0 1\nv 1 0 1\nv 1 1 1\nv 0 1 1\n"
        "f 1 2 3\nf 1 3 4\nf 5 7 6\nf 5 8 7\n"
        "f 1 5 6\nf 1 6 2\nf 2 6 7\nf 2 7 3\n"
        "f 3 7 8\nf 3 8 4\nf 4 8 5\nf 4 5 1\n");
    pqss_proxy_mesh::PrimitiveMeshAnalysisOptions coarse;
    coarse.analysis_strength = 0.0;
    coarse.allow_frustum = false;
    coarse.allow_polygon = true;
    const auto box_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(box, root / "box", coarse);
    if (check(box_stats.primitive_count == 6,
              "a box proxy must be exported as six rectangles")) return 1;
    if (check(box_stats.polygon_count == 6,
              "a box-shaped mesh must be emitted directly as six polygons")) return 1;
    if (check(box_stats.proxy_triangles == 12,
              "six box rectangles must triangulate to 12 faces")) return 1;

    // Surface simplification is a fixed-point process controlled only by the
    // directed proxy-to-source distance.  Three disconnected but coplanar
    // rectangles require two accepted merges; lowering the distance below the
    // half-gap must reject both without consulting triangle/BVH workload.
    const auto coplanar_gaps = root / "coplanar_gaps.obj";
    writeText(coplanar_gaps,
        "v 0 0 0\nv 1 0 0\nv 1 1 0\nv 0 1 0\n"
        "v 1.1 0 0\nv 2.1 0 0\nv 2.1 1 0\nv 1.1 1 0\n"
        "v 2.2 0 0\nv 3.2 0 0\nv 3.2 1 0\nv 2.2 1 0\n"
        "f 1 2 3\nf 1 3 4\nf 5 6 7\nf 5 7 8\n"
        "f 9 10 11\nf 9 11 12\n");
    auto staged_loose = coarse;
    staged_loose.maximum_open_error_distance = 0.051;
    const auto staged_loose_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        coplanar_gaps, root / "staged-coplanar-loose", staged_loose);
    if (check(staged_loose_stats.primitive_count == 1 &&
                  staged_loose_stats.proxy_triangles == 2 &&
                  staged_loose_stats.merged_local_planar_primitives == 2 &&
                  staged_loose_stats.open_max_distance <= 0.051 + 1.0e-9,
              "staged coplanar merging must iterate to a Hausdorff-compliant fixed point"))
        return 1;
    auto staged_strict = coarse;
    staged_strict.maximum_open_error_distance = 0.049;
    const auto staged_strict_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        coplanar_gaps, root / "staged-coplanar-strict", staged_strict);
    if (check(staged_strict_stats.primitive_count == 3 &&
                  staged_strict_stats.proxy_triangles == 6 &&
                  staged_strict_stats.merged_local_planar_primitives == 0,
              "a merge beyond the directed Hausdorff limit must remain split"))
        return 1;

    auto default_error_options = coarse;
    default_error_options.use_unified_candidate_optimizer = true;
    const auto default_error_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        box, root / "unified-default-error", default_error_options);
    if (check(std::abs(default_error_stats.maximum_open_error_distance_limit -
                       std::sqrt(3.0) * 0.08) < 1.0e-12,
              "the default maximum open error must be eight percent of the model diagonal"))
        return 1;

    // The unified optimizer accepts enclosing merges by one-sided
    // proxy-to-source distance, never by a watertight added-volume surrogate.
    pqss_proxy_mesh::PrimitiveMeshAnalysisOptions unified = coarse;
    unified.use_unified_candidate_optimizer = true;
    unified.use_staged_surface_pipeline = false;
    // Structural optimizer regressions below predate the user-facing distance
    // limit. Keep them focused on candidate construction; dedicated tests
    // below exercise strict maximum-error filtering.
    unified.maximum_open_error_distance = 1.0e6;
    const auto narrow_soup = root / "unified_narrow_soup.obj";
    writeText(narrow_soup, separatedBoxesObj(0.02));
    const auto narrow_soup_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        narrow_soup, root / "unified-narrow-soup", unified);
    if (check(narrow_soup_stats.proxy_triangles == 12,
              "the CPU watertight reference must close narrow model-scale cracks")) return 1;
    if (check(narrow_soup_stats.coverage_assigned_source_faces ==
                  narrow_soup_stats.source_triangles &&
              narrow_soup_stats.coverage_enclosure_source_faces ==
                  narrow_soup_stats.source_triangles &&
              narrow_soup_stats.coverage_unassigned_source_faces == 0 &&
              narrow_soup_stats.coverage_failed_source_faces == 0,
              "the final unified enclosure must certify every source triangle before export"))
        return 1;
    if (check(narrow_soup_stats.unified_selected_workload ==
                  static_cast<double>(narrow_soup_stats.proxy_triangles),
              "the unified optimizer must export the same post-union geometry whose workload it selected"))
        return 1;
    auto distance_compliant = unified;
    const auto distance_only_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        narrow_soup, root / "unified-distance-only-merge", distance_compliant);
    if (check(distance_only_stats.proxy_triangles == 12 &&
                  distance_only_stats.open_max_distance <=
                      distance_compliant.maximum_open_error_distance + 1.0e-8,
              "a distance-compliant enclosing merge must remain eligible")) return 1;
    auto open_gap_options = unified;
    const auto open_gap_soup = root / "unified_open_gap_soup.obj";
    writeText(open_gap_soup, separatedBoxesObj(0.4));
    const auto open_gap_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        open_gap_soup, root / "unified-open-gap-audit", open_gap_options);
    if (check(open_gap_stats.containment_validation_passed &&
                  open_gap_stats.open_max_distance > 0.0 &&
                  open_gap_stats.open_error_distance_sample_count > 0,
              "an exterior-accessible gap filled by a conservative proxy must produce a located open-region error"))
        return 1;
    if (check(std::filesystem::is_regular_file(
                  root / "unified-open-gap-audit" / "open_error.json") &&
                  readText(root / "unified-open-gap-audit" / "model.json")
                      .find("\"maximum_pair\"") != std::string::npos &&
                  readText(root / "unified-open-gap-audit" / "model.json")
                      .find("deterministic_area_surface_with_vertices_and_edges") !=
                      std::string::npos,
              "the final audit must persist its sampled maximum-error point pair and boundary-aware sampling method"))
        return 1;
    const auto wide_soup = root / "unified_wide_soup.obj";
    writeText(wide_soup, separatedBoxesObj(20.0));
    auto strict_wide_options = unified;
    strict_wide_options.maximum_open_error_distance = 0.1;
    const auto wide_soup_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        wide_soup, root / "unified-wide-soup", strict_wide_options);
    if (check(wide_soup_stats.proxy_triangles > 12,
              "a group envelope beyond the maximum one-sided distance must be rejected")) return 1;
    if (check(wide_soup_stats.open_max_distance <= 0.1 + 1.0e-8,
              "the final proxy must satisfy the requested maximum open-error distance")) return 1;
    if (check(wide_soup_stats.unified_adaptive_refinements > 0,
              "an over-limit enclosing candidate with finer children must be adaptively refined instead of discarded")) return 1;
    auto exact_wide_options = unified;
    exact_wide_options.maximum_open_error_distance = 0.0;
    const auto exact_wide_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        wide_soup, root / "unified-wide-soup-exact", exact_wide_options);
    const std::string exact_wide_metadata = readText(
        root / "unified-wide-soup-exact" / "model.json");
    if (check(exact_wide_stats.containment_validation_passed &&
                  exact_wide_stats.open_max_distance <= 1.0e-8 &&
                  exact_wide_metadata.find("inf") == std::string::npos &&
                  exact_wide_metadata.find("nan") == std::string::npos,
              "a zero maximum error must produce a finite exact conservative fallback")) return 1;

    // Intersecting conservative shells are converted to the boundary of their
    // solid union. Keeping the buried portions would make overlapping RSS
    // branches and duplicate triangle tests survive every near-contact query.
    const auto intersecting_boxes = root / "unified_intersecting_boxes.obj";
    writeText(intersecting_boxes, intersectingBoxesObj());
    auto intersecting_options = unified;
    intersecting_options.allow_frustum = false;
    intersecting_options.maximum_open_error_distance = 0.1;
    const auto intersecting_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        intersecting_boxes, root / "unified-intersecting-boxes",
        intersecting_options);
    const std::string intersecting_proxy = readText(
        root / "unified-intersecting-boxes" / "proxy.obj");
    if (check(intersecting_stats.unified_enclosure_extrusions == 2 &&
                  intersecting_stats.unified_occlusion_clipped_primitives > 0 &&
                  intersecting_stats.unified_occlusion_removed_area > 0.0,
              "duplicate extrusion axes must collapse before intersecting enclosures remove buried surface area"))
        return 1;
    if (check(!triangleMeshContainsPoint(
                  intersecting_proxy, {1.5, 0.25, 0.5}) &&
                  !triangleMeshContainsPoint(
                      intersecting_proxy, {1.5, 0.5, 1.0}),
              "the final proxy must not retain surfaces inside another enclosure"))
        return 1;
    if (check(triangleMeshContainsPoint(
                  intersecting_proxy, {1.5, 0.25, 1.25}),
              "union-boundary clipping must retain the exposed enclosure surface"))
        return 1;
    if (check(intersecting_stats.coverage_unassigned_source_faces == 0 &&
                  intersecting_stats.coverage_failed_source_faces == 0,
              "union-boundary clipping must preserve conservative source coverage"))
        return 1;

    const auto touching_boxes = root / "unified_partially_touching_boxes.obj";
    writeText(touching_boxes, partiallyTouchingBoxesObj());
    const auto touching_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        touching_boxes, root / "unified-partially-touching-boxes",
        intersecting_options);
    const std::string touching_proxy = readText(
        root / "unified-partially-touching-boxes" / "proxy.obj");
    if (check(touching_stats.unified_enclosure_extrusions == 2 &&
                  touching_stats.unified_occlusion_clipped_primitives > 0 &&
                  touching_stats.unified_occlusion_removed_area > 0.0 &&
                  touching_stats.removed_contained_primitives > 0 &&
                  !triangleMeshContainsPoint(touching_proxy, {2.0, 1.0, 1.0}) &&
                  triangleMeshContainsPoint(touching_proxy, {2.0, 0.25, 0.25}) &&
                  triangleMeshContainsPoint(touching_proxy, {3.0, 1.0, 1.0}),
              "pre-removal enclosure certificates must contract a shared interface after containment cleanup opens one output shell"))
        return 1;
    if (check(touching_stats.coverage_unassigned_source_faces == 0 &&
                  touching_stats.coverage_failed_source_faces == 0,
              "owned interface contraction must preserve conservative source coverage"))
        return 1;

    const auto edge_interface = root / "unified_edge_interface.obj";
    writeText(edge_interface,
        "v 0 0 0\nv 4 0 0\nv 4 4 0\nv 0 4 0\n"
        "v 0 0 4\nv 4 0 4\nv 4 4 4\nv 0 4 4\n"
        "f 1 2 3\nf 1 3 4\nf 5 7 6\nf 5 8 7\n"
        "f 1 5 6\nf 1 6 2\nf 2 6 7\nf 2 7 3\n"
        "f 3 7 8\nf 3 8 4\nf 4 8 5\nf 4 5 1\n"
        "v 4 0 0\nv 4 1 0\nv 4 1 4\nv 4 0 4\n"
        "v 5 0 0\nv 5 1 0\nv 5 1 4\nv 5 0 4\n"
        "f 9 10 11\nf 9 11 12\nf 13 15 14\nf 13 16 15\n"
        "f 9 13 14\nf 9 14 10\nf 10 14 15\nf 10 15 11\n"
        "f 11 15 16\nf 11 16 12\nf 12 16 13\nf 12 13 9\n");
    auto edge_interface_options = intersecting_options;
    edge_interface_options.maximum_open_error_distance = 0.1;
    const auto edge_interface_stats =
        pqss_proxy_mesh::analyzePrimitiveMeshObj(
            edge_interface, root / "unified-edge-interface",
            edge_interface_options);
    const std::string edge_interface_proxy = readText(
        root / "unified-edge-interface" / "proxy.obj");
    if (check(edge_interface_stats.unified_occlusion_clipped_primitives > 0 &&
                  !triangleMeshContainsPoint(
                      edge_interface_proxy, {4.0, 0.5, 2.0}) &&
                  triangleMeshContainsPoint(
                      edge_interface_proxy, {4.0, 2.0, 2.0}),
              "an adjacent closed volume must contract a covered edge strip without cutting a hole"))
        return 1;
    if (check(edge_interface_stats.coverage_unassigned_source_faces == 0 &&
                  edge_interface_stats.coverage_failed_source_faces == 0,
              "edge-interface contraction must preserve conservative source coverage"))
        return 1;

    const auto shallow_bumped_rectangle = root / "shallow_bumped_rectangle.obj";
    writeText(shallow_bumped_rectangle,
        "v 0 0 0\nv 100 0 0\nv 100 45 0\nv 100.05 50 0\n"
        "v 100 55 0\nv 100 100 0\nv 55 100 0\nv 50 100.05 0\n"
        "v 45 100 0\nv 0 100 0\nv 0 55 0\nv -0.05 50 0\n"
        "v 0 45 0\nv 50 50 0\n"
        "f 14 1 2\nf 14 2 3\nf 14 3 4\nf 14 4 5\n"
        "f 14 5 6\nf 14 6 7\nf 14 7 8\nf 14 8 9\n"
        "f 14 9 10\nf 14 10 11\nf 14 11 12\nf 14 12 13\n"
        "f 14 13 1\n");
    auto near_rectangle_options = unified;
    near_rectangle_options.tiny_planar_detail_area_ratio = 0.02;
    const auto shallow_bumped_stats =
        pqss_proxy_mesh::analyzePrimitiveMeshObj(
            shallow_bumped_rectangle, root / "shallow-bumped-rectangle",
            near_rectangle_options);
    const std::string shallow_bumped_proxy = readText(
        root / "shallow-bumped-rectangle" / "proxy.obj");
    if (check(shallow_bumped_stats.primitive_count == 1 &&
                  shallow_bumped_stats.proxy_triangles == 2 &&
                  triangleMeshCoversPlanarPoint(
                      shallow_bumped_proxy, 100.05, 50.0, 0.0) &&
                  triangleMeshCoversPlanarPoint(
                      shallow_bumped_proxy, -0.05, 50.0, 0.0),
              "a shallow near-rectangular outline must use its conservative enclosing rectangle"))
        return 1;

    const auto structural_hexagon = root / "unified_structural_hexagon.obj";
    writeText(structural_hexagon,
        "v 0 0 0\nv 2 0 0\nv 3 1 0\nv 2 2 0\nv 0 2 0\nv -1 1 0\n"
        "f 1 2 3\nf 1 3 4\nf 1 4 5\nf 1 5 6\n");
    const auto structural_hexagon_stats =
        pqss_proxy_mesh::analyzePrimitiveMeshObj(
            structural_hexagon, root / "unified-structural-hexagon",
            near_rectangle_options);
    if (check(structural_hexagon_stats.primitive_count == 1 &&
                  structural_hexagon_stats.proxy_triangles == 4,
              "a structurally non-rectangular silhouette must not pass shallow rectangle regularization"))
        return 1;

    const auto tessellated_medium_soup = root / "unified_tessellated_medium_soup.obj";
    writeText(tessellated_medium_soup,
              repeatObjFaces(separatedBoxesObj(0.5), 64));
    const auto tessellated_medium_stats =
        pqss_proxy_mesh::analyzePrimitiveMeshObj(
            tessellated_medium_soup,
            root / "unified-tessellated-medium-soup", unified);
    if (check(tessellated_medium_stats.proxy_triangles == 12,
              "source retessellation must not outweigh a feasible lower-work proxy")) return 1;

    const auto numerical_sliver = root / "numerical_sliver.obj";
    writeText(numerical_sliver,
        "v 0 0 0\nv 1 0 0\nv 0 1 0\n"
        "v 133.685 204.03 1321.68\n"
        "v 133.785 204.13 1321.68\n"
        "v 133.485 203.83 1321.68\n"
        "f 1 2 3\nf 4 5 6\n");
    const auto numerical_sliver_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        numerical_sliver, root / "numerical-sliver", coarse);
    if (check(numerical_sliver_stats.discarded_degenerate_triangles == 1 &&
              numerical_sliver_stats.proxy_triangles == 1,
              "scale-relative preprocessing must remove decimal-coordinate collinear slivers before primitive analysis")) return 1;

    const auto sparse_box_shell = root / "sparse_box_shell.obj";
    writeText(sparse_box_shell,
        // Two dominant supports plus a complex but narrow protruding strip. A
        // whole-box fit would create a mostly unsupported opposite side.
        "v 0 0 0\nv 0 10 0\nv 0 10 10\nv 0 0 10\n"
        "v 2 0 0\nv 2 10 0\nv 2 10 10\nv 2 0 10\n"
        "v 0 4 4\nv 4 4 4\nv 4 6 4\nv 0 6 4\n"
        "v 0 4 6\nv 4 4 6\nv 4 6 6\nv 0 6 6\n"
        "f 1 2 3\nf 1 3 4\nf 5 8 7\nf 5 7 6\n"
        "f 9 10 11\nf 9 11 12\nf 13 16 15\nf 13 15 14\n"
        "f 9 13 14\nf 9 14 10\nf 12 11 15\nf 12 15 16\n"
        "f 10 14 15\nf 10 15 11\n");
    pqss_proxy_mesh::PrimitiveMeshAnalysisOptions sparse_options = coarse;
    sparse_options.protrusion_max_area_excess_ratio = 0.03;
    const auto sparse_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        sparse_box_shell, root / "sparse-box-shell", sparse_options);
    if (check(sparse_stats.recognized_protrusion_box_shells == 0,
              "a protrusion box with a mostly unsupported output face must be rejected")) return 1;

    const auto disconnected_bounds = root / "disconnected_bounds.obj";
    writeText(disconnected_bounds,
        "v 0 0 0\nv 10 0 0\nv 10 10 0\nv 0 10 0\n"
        "v 0 0 1\nv 10 0 1\nv 10 10 1\nv 0 10 1\n"
        "v 4 -1 2\nv 6 -1 2\nv 6 0 2\nv 4 0 2\n"
        "v 4 -1 3\nv 6 -1 3\nv 6 0 3\nv 4 0 3\n"
        "f 1 2 3\nf 1 3 4\nf 5 7 6\nf 5 8 7\n"
        "f 1 5 6\nf 1 6 2\nf 2 6 7\nf 2 7 3\n"
        "f 3 7 8\nf 3 8 4\nf 4 8 5\nf 4 5 1\n"
        "f 9 10 11\nf 9 11 12\nf 13 15 14\nf 13 16 15\n"
        "f 9 13 14\nf 9 14 10\nf 10 14 15\nf 10 15 11\n"
        "f 11 15 16\nf 11 16 12\nf 12 16 13\nf 12 13 9\n");
    const auto disconnected_bounds_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        disconnected_bounds, root / "disconnected-bounds", coarse);
    const std::string disconnected_bounds_proxy =
        readText(root / "disconnected-bounds" / "proxy.obj");
    if (check(disconnected_bounds_stats.proxy_triangles > 0 &&
              triangleMeshCoversPlanarPoint(
                  disconnected_bounds_proxy, 5.0, 5.0, 0.0) &&
              !triangleMeshCoversPlanarPoint(
                  disconnected_bounds_proxy, 5.0, -0.5, 0.0),
              "a disconnected component must not stretch another component's local box envelope")) return 1;

    const auto close_component_group = root / "close_component_group.obj";
    writeText(close_component_group, separatedBoxesObj(0.02));
    const auto close_component_group_stats =
        pqss_proxy_mesh::analyzePrimitiveMeshObj(
            close_component_group, root / "close-component-group", coarse);
    if (check(close_component_group_stats.primitive_count == 6 &&
              close_component_group_stats.proxy_triangles == 12,
              "disconnected shells whose shared empty volume is globally in budget must merge into one envelope")) return 1;

    const auto distant_component_group = root / "distant_component_group.obj";
    writeText(distant_component_group, separatedBoxesObj(1.0));
    const auto distant_component_group_stats =
        pqss_proxy_mesh::analyzePrimitiveMeshObj(
            distant_component_group, root / "distant-component-group", coarse);
    if (check(distant_component_group_stats.recognized_closed_box_shells >= 4 &&
              distant_component_group_stats.proxy_triangles >= 48,
              "the same component hierarchy must reject a shared envelope when its empty volume is globally over budget")) return 1;

    const auto rectangle = root / "3.obj";
    writeText(rectangle, "v 0 0 0\nv 4 0 0\nv 4 2 0\nv 0 2 0\nf 1 2 3\nf 1 3 4\n");
    pqss_proxy_mesh::PrimitiveMeshAnalysisOptions planar;
    planar.analysis_strength = 1.0;
    planar.allow_frustum = false;
    planar.enable_volume_evaluated_envelope = false;
    const auto rectangle_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        rectangle, root / "rectangle", planar);
    if (check(rectangle_stats.polygon_count == 1, "a complete rectangle must remain one polygon")) return 1;
    if (check(rectangle_stats.proxy_triangles == 2, "a rectangle must triangulate to two faces")) return 1;

    const auto densely_notched = root / "densely_notched_planar.obj";
    writeText(densely_notched, denselyNotchedPlanarObj(256));
    const auto densely_notched_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        densely_notched, root / "densely-notched-planar", planar);
    const std::string densely_notched_proxy =
        readText(root / "densely-notched-planar" / "proxy.obj");
    if (check(densely_notched_stats.primitive_count == 1 &&
              triangleMeshCoversPlanarPoint(
                  densely_notched_proxy, 128.5, 1.0, 0.0),
              "dense shallow notches must simplify to one conservative planar primitive")) return 1;
    if (check(densely_notched_stats.analysis_seconds < 5.0,
              "dense planar boundary simplification must avoid all-edge quadratic validation per shortcut")) return 1;

    const auto ring = root / "ring.obj";
    writeText(ring,
        "v 0 0 0\nv 4 0 0\nv 4 4 0\nv 0 4 0\n"
        "v 1 1 0\nv 3 1 0\nv 3 3 0\nv 1 3 0\n"
        "f 1 2 6\nf 1 6 5\nf 2 3 7\nf 2 7 6\n"
        "f 3 4 8\nf 3 8 7\nf 4 1 5\nf 4 5 8\n");
    const auto ring_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        ring, root / "ring", planar);
    if (check(ring_stats.polygon_count == 1,
              "a rectangular patch with a hole must become one polygon")) return 1;
    if (check(ring_stats.filled_planar_holes == 1,
              "the inner boundary of a connected planar patch must be filled")) return 1;
    if (check(ring_stats.proxy_triangles == 2,
              "a filled rectangular patch must triangulate to two faces")) return 1;

    auto unified_zero_error_planar = planar;
    unified_zero_error_planar.use_unified_candidate_optimizer = true;
    unified_zero_error_planar.maximum_open_error_distance = 0.0;
    const auto unified_ring_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        ring, root / "unified-ring-zero-error", unified_zero_error_planar);
    if (check(unified_ring_stats.containment_validation_passed &&
                  !unified_ring_stats.unified_exact_fallback_selected &&
                  unified_ring_stats.proxy_triangles == 2 &&
                  unified_ring_stats.open_max_distance <= 1.0e-8,
              "a certified planar inner loop must be fillable without consuming the open-error budget")) return 1;

    const auto concave_l = root / "concave_l.obj";
    writeText(concave_l,
        "v 0 0 0\nv 4 0 0\nv 4 1 0\nv 1 1 0\nv 1 4 0\nv 0 4 0\n"
        "f 1 2 4\nf 2 3 4\nf 1 4 6\nf 4 5 6\n");
    const auto concave_l_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        concave_l, root / "unified-concave-l-zero-error",
        unified_zero_error_planar);
    const std::string concave_l_proxy = readText(
        root / "unified-concave-l-zero-error" / "proxy.obj");
    if (check(concave_l_stats.containment_validation_passed &&
                  concave_l_stats.proxy_triangles == 4 &&
                  !triangleMeshCoversPlanarPoint(
                      concave_l_proxy, 2.0, 2.0, 0.0),
              "an exterior concavity is not a planar hole and must remain chargeable")) return 1;

    const auto islands = root / "islands.obj";
    writeText(islands,
        "v 0 0 0\nv 1 0 0\nv 1 1 0\nv 0 1 0\n"
        "v 3 0 0\nv 4 0 0\nv 4 1 0\nv 3 1 0\n"
        "f 1 2 3\nf 1 3 4\nf 5 6 7\nf 5 7 8\n");
    const auto separate_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        islands, root / "islands-separate", planar);
    if (check(separate_stats.primitive_count == 2,
              "fine analysis must not merge disconnected coplanar objects")) return 1;

    const auto many_coplanar_islands = root / "many_coplanar_islands.obj";
    std::ostringstream many_coplanar_islands_obj;
    constexpr std::size_t island_count = 256;
    for (std::size_t island = 0; island < island_count; ++island)
    {
        const double x = static_cast<double>(island % 32) * 3.0;
        const double y = static_cast<double>(island / 32) * 3.0;
        many_coplanar_islands_obj
            << "v " << x << ' ' << y << " 0\n"
            << "v " << x + 1.0 << ' ' << y << " 0\n"
            << "v " << x << ' ' << y + 1.0 << " 0\n";
    }
    for (std::size_t island = 0; island < island_count; ++island)
    {
        const std::size_t first = island * 3 + 1;
        many_coplanar_islands_obj << "f " << first << ' ' << first + 1
                                  << ' ' << first + 2 << '\n';
    }
    writeText(many_coplanar_islands, many_coplanar_islands_obj.str());
    pqss_proxy_mesh::PrimitiveMeshAnalysisOptions disconnected_planar = planar;
    disconnected_planar.maximum_local_planar_fill_area_ratio = 0.0;
    const auto many_coplanar_islands_stats =
        pqss_proxy_mesh::analyzePrimitiveMeshObj(
            many_coplanar_islands, root / "many-coplanar-islands", disconnected_planar);
    if (check(many_coplanar_islands_stats.primitive_count == island_count &&
              many_coplanar_islands_stats.canonicalized_coplanar_groups == 0,
              "distant patches in one plane bucket must remain separate Boolean components")) return 1;

    pqss_proxy_mesh::PrimitiveMeshAnalysisOptions joined = planar;
    joined.analysis_strength = 0.0;
    joined.allow_frustum = true;
    joined.uniform_structure_policy = false;
    const auto joined_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        islands, root / "islands-joined", joined);
    if (check(joined_stats.frustum_count == 0,
              "disconnected planar objects must not be disguised as one circular primitive")) return 1;

    const auto globally_small_gap = root / "globally_small_gap.obj";
    writeText(globally_small_gap,
        "v 0 0 0\nv 1 0 0\nv 1 1 0\nv 0 1 0\n"
        "v 1.2 0 0\nv 2.2 0 0\nv 2.2 1 0\nv 1.2 1 0\n"
        "v 100 0 10\nv 110 0 10\nv 110 10 10\nv 100 10 10\n"
        "f 1 2 3\nf 1 3 4\nf 5 6 7\nf 5 7 8\nf 9 10 11\nf 9 11 12\n");
    pqss_proxy_mesh::PrimitiveMeshAnalysisOptions global_area = planar;
    global_area.maximum_local_planar_fill_area_ratio = 0.003;
    const auto global_area_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        globally_small_gap, root / "globally-small-gap", global_area);
    if (check(global_area_stats.merged_local_planar_primitives == 1,
              "planar merge error must be normalized by total model surface area")) return 1;
    if (check(global_area_stats.primitive_count == 2,
              "two small patches should merge while the unrelated model surface remains separate")) return 1;

    const auto different_depths = root / "different_depths.obj";
    writeText(different_depths,
        "v 0 0 0\nv 1 0 0\nv 1 1 0\nv 0 1 0\n"
        "v 1.2 0 1\nv 2.2 0 1\nv 2.2 1 1\nv 1.2 1 1\n"
        "f 1 2 3\nf 1 3 4\nf 5 6 7\nf 5 7 8\n");
    pqss_proxy_mesh::PrimitiveMeshAnalysisOptions unlimited_planar_merge = planar;
    unlimited_planar_merge.maximum_local_planar_fill_area_ratio =
        std::numeric_limits<double>::infinity();
    const auto different_depth_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        different_depths, root / "different-depths", unlimited_planar_merge);
    if (check(different_depth_stats.merged_local_planar_primitives == 0,
              "global error budget must never merge parallel patches at different depths")) return 1;

    const auto coplanar_overlap = root / "coplanar_overlap.obj";
    writeText(coplanar_overlap,
        "v 0 0 0\nv 4 0 0\nv 4 4 0\n"
        "v 1 1 0\nv 5 1 0\nv 5 5 0\nv 1 5 0\n"
        "f 1 2 3\nf 4 5 6\nf 4 6 7\n");
    pqss_proxy_mesh::PrimitiveMeshAnalysisOptions overlap_options = planar;
    overlap_options.maximum_local_planar_fill_area_ratio = 0.0;
    const auto overlap_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        coplanar_overlap, root / "coplanar-overlap", overlap_options);
    if (check(overlap_stats.canonicalized_coplanar_groups == 1,
              "geometrically coplanar overlapping regions must be unioned without shared OBJ edges")) return 1;
    if (check(overlap_stats.removed_coplanar_overlap_area > 4.49 &&
              overlap_stats.removed_coplanar_overlap_area < 4.51,
              "coplanar union must report the duplicated covered area")) return 1;

    const auto split_rectangle = root / "split_rectangle.obj";
    writeText(split_rectangle,
        "v 0 0 0\nv 2 0 0\nv 2 2 0\nv 0 2 0\n"
        "v 2 0 0\nv 4 0 0\nv 4 2 0\nv 2 2 0\n"
        "f 1 2 3\nf 1 3 4\nf 5 6 7\nf 5 7 8\n");
    const auto split_rectangle_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        split_rectangle, root / "split-rectangle", overlap_options);
    if (check(split_rectangle_stats.polygon_count == 1 &&
              split_rectangle_stats.primitive_count == 1 &&
              split_rectangle_stats.proxy_triangles == 2,
              "a rectangular planar union must be restored as one polygon primitive")) return 1;

    const auto parallel_overlap = root / "parallel_overlap.obj";
    writeText(parallel_overlap,
        "v 0 0 0\nv 4 0 0\nv 4 4 0\n"
        "v 1 1 1\nv 5 1 1\nv 5 5 1\nv 1 5 1\n"
        "f 1 2 3\nf 4 5 6\nf 4 6 7\n");
    const auto parallel_overlap_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        parallel_overlap, root / "parallel-overlap", overlap_options);
    if (check(parallel_overlap_stats.canonicalized_coplanar_groups == 0,
              "parallel regions at different depths must not enter one planar union")) return 1;

    const auto pentagon = root / "4.obj";
    writeText(pentagon,
        "v 0 0 0\nv 3 0 0\nv 3 2 0\nv 1 3 0\nv 0 2 0\n"
        "f 1 2 3\nf 1 3 4\nf 1 4 5\n");
    const auto pentagon_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        pentagon, root / "pentagon", planar);
    if (check(pentagon_stats.polygon_count == 1, "a planar pentagon must remain one polygon primitive")) return 1;
    if (check(pentagon_stats.proxy_triangles == 3, "pentagon proxy OBJ must have three faces")) return 1;
    const std::string pentagon_primitives = readText(root / "pentagon" / "primitives.obj");
    if (check(pentagon_primitives.find("g primitive_00000_polygon") != std::string::npos &&
              pentagon_primitives.find("f 1 2 3 4 5") != std::string::npos,
              "primitive analysis OBJ must preserve the pentagon as one five-vertex face")) return 1;

    const auto hexagon = root / "hexagon.obj";
    writeText(hexagon,
        "v 0 0 0\nv 2 0 0\nv 3 1 0\nv 2 2 0\nv 0 2 0\nv -1 1 0\n"
        "f 1 2 3\nf 1 3 4\nf 1 4 5\nf 1 5 6\n");
    const auto hexagon_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        hexagon, root / "hexagon", planar);
    if (check(hexagon_stats.polygon_count == 1 && hexagon_stats.proxy_triangles == 4,
              "one hexagonal primitive must triangulate to exactly four faces")) return 1;
    const std::string hexagon_metadata = readText(root / "hexagon" / "model.json");
    if (check(hexagon_metadata.find("\"type\":\"polygon\",\"vertex_count\":6,\"triangulated_face_count\":4") != std::string::npos,
              "metadata must expose the semantic hexagon and its four collision triangles")) return 1;

    const auto globally_small_hexagon_notch = root / "globally_small_hexagon_notch.obj";
    writeText(globally_small_hexagon_notch,
        "v 0 0 0\nv 4 0 0\nv 4 1 0\nv 3 2 0\nv 0 2 0\nv -1 1 0\n"
        "v 0 0 10\nv 100 0 10\nv 100 100 10\nv 0 100 10\n"
        "f 1 2 3\nf 1 3 4\nf 1 4 5\nf 1 5 6\n"
        "f 7 8 9\nf 7 9 10\n");
    const auto globally_small_hexagon_notch_stats =
        pqss_proxy_mesh::analyzePrimitiveMeshObj(
            globally_small_hexagon_notch, root / "globally-small-hexagon-notch", planar);
    const std::string globally_small_hexagon_notch_metadata =
        readText(root / "globally-small-hexagon-notch" / "model.json");
    if (check(globally_small_hexagon_notch_stats.primitive_count == 2 &&
              globally_small_hexagon_notch_metadata.find(
                  "\"type\":\"polygon\",\"vertex_count\":6,\"triangulated_face_count\":4") !=
                  std::string::npos,
              "global area normalization must not authorize replacing a structural hexagonal silhouette with its bounding rectangle")) return 1;

    const auto raised_parallel_patch = root / "raised_parallel_patch.obj";
    writeText(raised_parallel_patch,
        "v 0 0 0\nv 4 0 0\nv 4 4 0\nv 0 4 0\n"
        "v 1 1 0.1\nv 3 1 0.1\nv 3 3 0.1\nv 1 3 0.1\n"
        "v 0 0 10\nv 1 0 10\nv 0 1 10\n"
        "f 1 2 3\nf 1 3 4\n"
        "f 5 6 7\nf 5 7 8\n"
        "f 9 10 11\n");
    const auto raised_parallel_patch_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        raised_parallel_patch, root / "raised-parallel-patch", overlap_options);
    const std::string raised_parallel_patch_primitives =
        readText(root / "raised-parallel-patch" / "primitives.obj");
    if (check(raised_parallel_patch_stats.primitive_count == 3 &&
              raised_parallel_patch_primitives.find(" 0.10000000000000001") !=
                  std::string::npos,
              "projected coverage by a lower parallel plane must not delete a raised surface")) return 1;

    const auto outward_foot = root / "outward_foot.obj";
    writeText(outward_foot,
        "v 0 0 0\nv 10 0 0\nv 10 10 0\nv 0 10 0\n"
        "v 4 4 0\nv 6 4 0\nv 6 6 0\nv 4 6 0\n"
        "v 4 4 -1\nv 6 4 -1\nv 6 6 -1\nv 4 6 -1\n"
        "v 0 0 10\nv 1 0 10\nv 0 1 10\n"
        "f 1 2 3\nf 1 3 4\n"
        "f 5 9 10\nf 5 10 6\nf 6 10 11\nf 6 11 7\n"
        "f 7 11 12\nf 7 12 8\nf 8 12 9\nf 8 9 5\n"
        "f 9 12 11\nf 9 11 10\n"
        "f 13 14 15\n");
    const auto outward_foot_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        outward_foot, root / "outward-foot", overlap_options);
    const std::string outward_foot_source =
        readText(root / "outward-foot" / "source.obj");
    const std::string outward_foot_proxy =
        readText(root / "outward-foot" / "primitives.obj");
    if (check(outward_foot_stats.primitive_count >= 3 &&
              outward_foot_source.find(" -1\n") != std::string::npos &&
              outward_foot_proxy.find(" -1") != std::string::npos,
              "a shallow feature on the outward side of a support plane must not be removed as a blind cavity")) return 1;

    const auto shallow_box_with_foot = root / "shallow_box_with_foot.obj";
    const std::string shallow_box_vertices =
        "v 0 0 0\nv 10 0 0\nv 10 10 0\nv 0 10 0\n"
        "v 0 0 0.5\nv 10 0 0.5\nv 10 10 0.5\nv 0 10 0.5\n"
        "v 4 4 -0.2\nv 6 4 -0.2\nv 6 6 -0.2\nv 4 6 -0.2\n"
        "v 4 4 0\nv 6 4 0\nv 6 6 0\nv 4 6 0\n"
        "v 0 0 10\nv 1 0 10\nv 0 1 10\n";
    const std::string shallow_box_faces =
        "f 1 2 3\nf 1 3 4\nf 5 7 6\nf 5 8 7\n"
        "f 1 5 6\nf 1 6 2\nf 2 6 7\nf 2 7 3\n"
        "f 3 7 8\nf 3 8 4\nf 4 8 5\nf 4 5 1\n"
        "f 9 10 11\nf 9 11 12\nf 13 15 14\nf 13 16 15\n"
        "f 9 13 14\nf 9 14 10\nf 10 14 15\nf 10 15 11\n"
        "f 11 15 16\nf 11 16 12\nf 12 16 13\nf 12 13 9\n";
    writeText(shallow_box_with_foot,
        shallow_box_vertices + shallow_box_faces + shallow_box_faces +
        "f 17 18 19\n");
    const auto shallow_box_with_foot_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        shallow_box_with_foot, root / "shallow-box-with-foot", overlap_options);
    const std::string shallow_box_with_foot_proxy =
        readText(root / "shallow-box-with-foot" / "primitives.obj");
    if (check(shallow_box_with_foot_stats.recognized_closed_box_shells >= 1 &&
              shallow_box_with_foot_proxy.find("v 0 0 -0.199") !=
                  std::string::npos &&
              shallow_box_with_foot_proxy.find("v 10 10 0.5") != std::string::npos,
              "a base and shallow outward feet must merge into one complete conservative box envelope")) return 1;

    const auto shallow_hexagonal_shell = root / "shallow_hexagonal_shell.obj";
    writeText(shallow_hexagonal_shell,
        "v 0 0 0\nv 4 0 0\nv 4 0 1\nv 3 0 2\nv 0 0 2\nv -1 0 1\n"
        "v 0 0.01 0\nv 4 0.01 0\nv 4 0.01 1\nv 3 0.01 2\nv 0 0.01 2\nv -1 0.01 1\n"
        "v 0 -10 0\nv 1 -10 0\nv 0 -10 1\n"
        "f 1 2 3\nf 1 3 4\nf 1 4 5\nf 1 5 6\n"
        "f 7 9 8\nf 7 10 9\nf 7 11 10\nf 7 12 11\n"
        "f 1 7 8\nf 1 8 2\nf 2 8 9\nf 2 9 3\n"
        "f 3 9 10\nf 3 10 4\nf 4 10 11\nf 4 11 5\n"
        "f 5 11 12\nf 5 12 6\nf 6 12 7\nf 6 7 1\n"
        "f 13 14 15\n");
    const auto shallow_hexagonal_shell_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        shallow_hexagonal_shell, root / "shallow-hexagonal-shell", overlap_options);
    const std::string shallow_hexagonal_shell_metadata =
        readText(root / "shallow-hexagonal-shell" / "model.json");
    if (check(shallow_hexagonal_shell_stats.canonicalized_coplanar_groups >= 1 &&
              shallow_hexagonal_shell_metadata.find(
                  "\"type\":\"polygon\",\"vertex_count\":6,\"triangulated_face_count\":4") !=
                  std::string::npos,
              "a connected shallow parallel shell must coalesce onto one outward plane while preserving its six-vertex silhouette")) return 1;

    const auto short_bevel_shell = root / "short_bevel_shell.obj";
    writeText(short_bevel_shell,
        "v 0 0 0\nv 4 0 0\nv 4 0 1.95\nv 3.98 0 2.02\nv 2 0 4\nv 0 0 4\n"
        "v 0 0.01 0\nv 4 0.01 0\nv 4 0.01 1.95\nv 3.98 0.01 2.02\nv 2 0.01 4\nv 0 0.01 4\n"
        "f 1 2 3\nf 1 3 4\nf 1 4 5\nf 1 5 6\n"
        "f 7 9 8\nf 7 10 9\nf 7 11 10\nf 7 12 11\n"
        "f 1 7 8\nf 1 8 2\nf 2 8 9\nf 2 9 3\n"
        "f 3 9 10\nf 3 10 4\nf 4 10 11\nf 4 11 5\n"
        "f 5 11 12\nf 5 12 6\nf 6 12 7\nf 6 7 1\n");
    const auto short_bevel_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        short_bevel_shell, root / "short-bevel-shell", overlap_options);
    const std::string short_bevel_metadata =
        readText(root / "short-bevel-shell" / "model.json");
    if (check(short_bevel_stats.polygon_count == 7 &&
              short_bevel_metadata.find(
                  "\"type\":\"polygon\",\"vertex_count\":5,\"triangulated_face_count\":3") !=
                  std::string::npos,
              "a short bevel must simplify by extending its adjacent main edges, not by creating a box corner")) return 1;

    const auto terraced_recess = root / "terraced_recess.obj";
    writeText(terraced_recess,
        "v 0 0 0\nv 0 2 0\nv 0 2 4\nv 0 0 4\n"
        "v 0 4 0\nv 0 6 0\nv 0 6 4\nv 0 4 4\n"
        "v 0.5 0 0\nv 0.5 6 0\nv 0.5 6 4\nv 0.5 0 4\n"
        "f 1 2 3\nf 1 3 4\nf 5 6 7\nf 5 7 8\n"
        "f 9 10 11\nf 9 11 12\n"
        "f 2 5 8\nf 2 8 3\nf 1 9 12\nf 1 12 4\n"
        "f 6 10 11\nf 6 11 7\n");
    pqss_proxy_mesh::PrimitiveMeshAnalysisOptions terrace_options = overlap_options;
    terrace_options.maximum_cavity_added_volume_ratio = 1.0;
    terrace_options.shallow_parallel_merge_depth_relative = 0.01;
    const auto terrace_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        terraced_recess, root / "terraced-recess", terrace_options);
    if (check(terrace_stats.primitive_count < 6,
              "a connected terraced recess must be filled by globally evaluated layered sweep volume")) return 1;

    const auto synchronized_extrusion = root / "synchronized_extrusion.obj";
    writeText(synchronized_extrusion,
        // Identical conservative caps describe a chamfered extrusion.  The
        // source side responsibility is terraced and entirely inside that
        // extrusion; a parallel-shell merge may first broaden it to x=0.  The
        // semantic shell must then rebuild the side from the cap boundary,
        // trading one extra quad for removal of the empty triangular prism.
        "v 0.1 0 0\nv 0 0 1\nv 0 0 4\nv 4 0 4\nv 4 0 0\n"
        "v 0.1 10 0\nv 0 10 1\nv 0 10 4\nv 4 10 4\nv 4 10 0\n"
        "v 0.1 0 1\nv 0.1 10 1\n"
        "f 1 2 3\nf 1 3 4\nf 1 4 5\n"
        "f 6 8 7\nf 6 9 8\nf 6 10 9\n"
        "f 1 6 12\nf 1 12 11\n"
        "f 11 12 7\nf 11 7 2\n"
        "f 2 7 8\nf 2 8 3\n"
        "f 3 8 9\nf 3 9 4\n"
        "f 4 9 10\nf 4 10 5\n"
        "f 5 10 6\nf 5 6 1\n");
    pqss_proxy_mesh::PrimitiveMeshAnalysisOptions extrusion_options = overlap_options;
    extrusion_options.maximum_cavity_added_volume_ratio = 1.0;
    extrusion_options.shallow_parallel_merge_depth_relative = 0.02;
    const auto extrusion_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        synchronized_extrusion, root / "synchronized-extrusion", extrusion_options);
    const std::string extrusion_proxy =
        readText(root / "synchronized-extrusion" / "primitives.obj");
    if (check(extrusion_stats.proxy_triangles >= 14 &&
              extrusion_proxy.find("v 0.10000000000000001 0 0") != std::string::npos &&
              extrusion_proxy.find("v 0 0 1") != std::string::npos &&
              extrusion_proxy.find("v 0 10 1") != std::string::npos &&
              extrusion_proxy.find("v 0.10000000000000001 10 0") != std::string::npos &&
              extrusion_proxy.find("v 0 0 0\n") == std::string::npos,
              "synchronized opposing silhouettes must rebuild a tighter certified extruded side instead of retaining an empty triangular prism")) return 1;

    const auto chamfered_rectangle = root / "chamfered_rectangle.obj";
    writeText(chamfered_rectangle,
        "v 0 0 0\nv 4 0 0\nv 4 1.92 0\nv 3.92 2 0\nv 0 2 0\n"
        "f 1 2 3\nf 1 3 4\nf 1 4 5\n");
    pqss_proxy_mesh::PrimitiveMeshAnalysisOptions boundary_void = planar;
    boundary_void.analysis_strength = 0.5;
    const auto boundary_void_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        chamfered_rectangle, root / "chamfered-rectangle", boundary_void);
    if (check(boundary_void_stats.polygon_count == 1,
              "a small boundary chamfer must be closed before planar classification")) return 1;
    if (check(boundary_void_stats.filled_boundary_voids == 1,
              "the closed boundary chamfer must be reported")) return 1;
    const auto preserved_chamfer_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        chamfered_rectangle, root / "preserved-chamfer", planar);
    if (check(preserved_chamfer_stats.filled_boundary_voids == 1,
              "small structural boundary voids must also be closed at fine strength")) return 1;

    const auto chamfered_prism = root / "chamfered_prism.obj";
    writeText(chamfered_prism,
        "v 0 0 0\nv 4 0 0\nv 4 1.92 0\nv 3.92 2 0\nv 0 2 0\n"
        "v 0 0 1\nv 4 0 1\nv 4 1.92 1\nv 3.92 2 1\nv 0 2 1\n"
        "f 1 2 3\nf 1 3 4\nf 1 4 5\n"
        "f 6 8 7\nf 6 9 8\nf 6 10 9\n"
        "f 1 6 7\nf 1 7 2\nf 2 7 8\nf 2 8 3\n"
        "f 3 8 9\nf 3 9 4\nf 4 9 10\nf 4 10 5\n"
        "f 5 10 6\nf 5 6 1\n");
    const auto chamfered_prism_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        chamfered_prism, root / "chamfered-prism", planar);
    if (check(chamfered_prism_stats.filled_boundary_voids >= 2 &&
              chamfered_prism_stats.removed_sealed_void_wall_primitives >= 1 &&
              chamfered_prism_stats.primitive_count == 6 &&
              chamfered_prism_stats.proxy_triangles == 12,
              "filling opposing boundary voids must remove their internal bevel wall during structural cleanup")) return 1;

    auto unified_chamfer = planar;
    unified_chamfer.use_unified_candidate_optimizer = true;
    unified_chamfer.maximum_open_error_distance = 0.2;
    const auto unified_chamfer_stats =
        pqss_proxy_mesh::analyzePrimitiveMeshObj(
            chamfered_prism, root / "unified-chamfered-prism",
            unified_chamfer);
    const auto unified_chamfer_proxy_path =
        root / "unified-chamfered-prism" / "proxy.obj";
    const std::string unified_chamfer_proxy =
        readText(unified_chamfer_proxy_path);
    if (check(unified_chamfer_stats.coverage_unassigned_source_faces == 0 &&
                  unified_chamfer_stats.coverage_failed_source_faces == 0 &&
                  unified_chamfer_stats.proxy_triangles == 12 &&
                  objIsClosedTriangleSoup(unified_chamfer_proxy_path) &&
                  triangleMeshCoversPlanarPoint(
                      unified_chamfer_proxy, 4.0, 2.0, 0.0) &&
                  triangleMeshCoversPlanarPoint(
                      unified_chamfer_proxy, 4.0, 2.0, 1.0),
              "unified boundary-void cleanup must emit both caps before removing the bevel wall"))
        return 1;

    auto strict_unified_chamfer = unified_chamfer;
    strict_unified_chamfer.maximum_open_error_distance = 0.01;
    const auto strict_unified_chamfer_stats =
        pqss_proxy_mesh::analyzePrimitiveMeshObj(
            chamfered_prism, root / "strict-unified-chamfered-prism",
            strict_unified_chamfer);
    const auto strict_unified_chamfer_proxy_path =
        root / "strict-unified-chamfered-prism" / "proxy.obj";
    const std::string strict_unified_chamfer_proxy =
        readText(strict_unified_chamfer_proxy_path);
    if (check(strict_unified_chamfer_stats.coverage_unassigned_source_faces == 0 &&
                  strict_unified_chamfer_stats.coverage_failed_source_faces == 0 &&
                  objIsClosedTriangleSoup(strict_unified_chamfer_proxy_path) &&
                  !triangleMeshCoversPlanarPoint(
                      strict_unified_chamfer_proxy, 4.0, 2.0, 0.0) &&
                  !triangleMeshCoversPlanarPoint(
                      strict_unified_chamfer_proxy, 4.0, 2.0, 1.0),
              "a rejected cap must retain the bevel wall instead of deleting it independently"))
        return 1;

    const auto sealed_tube = root / "sealed_tube.obj";
    writeText(sealed_tube,
        "v 0 0 0\nv 4 0 0\nv 4 4 0\nv 0 4 0\n"
        "v 1 1 0\nv 3 1 0\nv 3 3 0\nv 1 3 0\n"
        "v 0 0 1\nv 4 0 1\nv 4 4 1\nv 0 4 1\n"
        "v 1 1 1\nv 3 1 1\nv 3 3 1\nv 1 3 1\n"
        "f 1 2 6\nf 1 6 5\nf 2 3 7\nf 2 7 6\nf 3 4 8\nf 3 8 7\nf 4 1 5\nf 4 5 8\n"
        "f 9 13 14\nf 9 14 10\nf 10 14 15\nf 10 15 11\nf 11 15 16\nf 11 16 12\nf 12 16 13\nf 12 13 9\n"
        "f 5 6 14\nf 5 14 13\nf 6 7 15\nf 6 15 14\nf 7 8 16\nf 7 16 15\nf 8 5 13\nf 8 13 16\n");
    const auto preserved_tube_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        sealed_tube, root / "preserved-tube", planar);
    if (check(preserved_tube_stats.removed_sealed_void_wall_primitives == 0,
              "a cavity whose individual swept volume exceeds the whole-model ratio must retain its walls")) return 1;
    auto fill_tube = planar;
    fill_tube.maximum_cavity_added_volume_ratio = 1.0;
    const auto sealed_tube_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        sealed_tube, root / "sealed-tube", fill_tube);
    if (check(sealed_tube_stats.removed_sealed_void_wall_primitives == 4,
              "four inner walls between filled hole caps must be removed")) return 1;
    if (check(sealed_tube_stats.excluded_sealed_void_wall_triangles == 8,
              "all source triangles of the sealed inner walls must be excluded")) return 1;
    const std::string sealed_tube_proxy = readText(root / "sealed-tube" / "proxy.obj");
    if (check(triangleMeshCoversPlanarPoint(sealed_tube_proxy, 2.0, 2.0, 0.0) &&
              triangleMeshCoversPlanarPoint(sealed_tube_proxy, 2.0, 2.0, 1.0),
              "an accepted through-hole fill must cover both sides of the local thin shell")) return 1;

    auto unified_fill_tube = fill_tube;
    unified_fill_tube.use_unified_candidate_optimizer = true;
    unified_fill_tube.maximum_open_error_distance = 1.0e6;
    const auto unified_sealed_tube_stats =
        pqss_proxy_mesh::analyzePrimitiveMeshObj(
            sealed_tube, root / "unified-sealed-tube", unified_fill_tube);
    if (check(unified_sealed_tube_stats.removed_sealed_void_wall_primitives == 4 &&
                  unified_sealed_tube_stats.excluded_sealed_void_wall_triangles == 8 &&
                  unified_sealed_tube_stats.primitive_count == 2 &&
                  unified_sealed_tube_stats.proxy_triangles == 4,
              "unified optimization must not refit the inner walls of an accepted hole fill as an enclosure"))
        return 1;
    if (check(unified_sealed_tube_stats.coverage_unassigned_source_faces == 0 &&
                  unified_sealed_tube_stats.coverage_failed_source_faces == 0,
              "excluding certified filled-hole walls must preserve unified conservative coverage"))
        return 1;

    // Two opposing, index-disconnected cavity walls exercise geometric void
    // grouping. OBJ triangle soups need not connect every boundary surface, so
    // preserving the opening cannot depend on wall-to-wall connectivity.
    const auto disconnected_tube = root / "disconnected_tube.obj";
    writeText(disconnected_tube,
        "v 0 0 0\nv 4 0 0\nv 4 4 0\nv 0 4 0\n"
        "v 1 1 0\nv 3 1 0\nv 3 3 0\nv 1 3 0\n"
        "v 0 0 1\nv 4 0 1\nv 4 4 1\nv 0 4 1\n"
        "v 1 1 1\nv 3 1 1\nv 3 3 1\nv 1 3 1\n"
        "f 1 2 6\nf 1 6 5\nf 2 3 7\nf 2 7 6\nf 3 4 8\nf 3 8 7\nf 4 1 5\nf 4 5 8\n"
        "f 9 13 14\nf 9 14 10\nf 10 14 15\nf 10 15 11\nf 11 15 16\nf 11 16 12\nf 12 16 13\nf 12 13 9\n"
        "f 5 6 14\nf 5 14 13\nf 7 8 16\nf 7 16 15\n");
    const auto disconnected_tube_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        disconnected_tube, root / "disconnected-tube", planar);
    const std::string disconnected_proxy =
        readText(root / "disconnected-tube" / "proxy.obj");
    if (check(disconnected_tube_stats.removed_sealed_void_wall_primitives == 0,
              "an over-budget cavity must retain independently indexed inner walls")) return 1;
    if (check(!triangleMeshCoversPlanarPoint(disconnected_proxy, 2.0, 2.0, 0.0) &&
              !triangleMeshCoversPlanarPoint(disconnected_proxy, 2.0, 2.0, 1.0),
              "geometrically grouped cavity walls must reopen both synthetic cap planes")) return 1;

    const auto rotated_disconnected_tube = root / "rotated_disconnected_tube.obj";
    writeText(rotated_disconnected_tube,
        "v 2 -0.8284271247461903 0\nv 4.82842712474619 2 0\n"
        "v 2 4.82842712474619 0\nv -0.8284271247461903 2 0\n"
        "v 2 0.5857864376269049 0\nv 3.414213562373095 2 0\n"
        "v 2 3.414213562373095 0\nv 0.5857864376269049 2 0\n"
        "v 2 -0.8284271247461903 1\nv 4.82842712474619 2 1\n"
        "v 2 4.82842712474619 1\nv -0.8284271247461903 2 1\n"
        "v 2 0.5857864376269049 1\nv 3.414213562373095 2 1\n"
        "v 2 3.414213562373095 1\nv 0.5857864376269049 2 1\n"
        "f 1 2 6\nf 1 6 5\nf 2 3 7\nf 2 7 6\nf 3 4 8\nf 3 8 7\nf 4 1 5\nf 4 5 8\n"
        "f 9 13 14\nf 9 14 10\nf 10 14 15\nf 10 15 11\nf 11 15 16\nf 11 16 12\nf 12 16 13\nf 12 13 9\n"
        "f 5 6 14\nf 5 14 13\nf 7 8 16\nf 7 16 15\n");
    const auto rotated_disconnected_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        rotated_disconnected_tube, root / "rotated-disconnected-tube", planar);
    const std::string rotated_disconnected_proxy =
        readText(root / "rotated-disconnected-tube" / "proxy.obj");
    if (check(rotated_disconnected_stats.primitive_count ==
                  disconnected_tube_stats.primitive_count &&
              rotated_disconnected_stats.proxy_triangles ==
                  disconnected_tube_stats.proxy_triangles &&
              !triangleMeshCoversPlanarPoint(
                  rotated_disconnected_proxy, 2.0, 2.0, 0.0) &&
              !triangleMeshCoversPlanarPoint(
                  rotated_disconnected_proxy, 2.0, 2.0, 1.0),
              "cavity grouping and planar reopening must be invariant under rigid rotation")) return 1;

    const auto circle = root / "5.obj";
    std::ofstream circle_stream(circle);
    circle_stream << "v 0 0 0\n";
    constexpr int segments = 12;
    for (int index = 0; index < segments; ++index)
    {
        const double angle = 2.0 * 3.14159265358979323846 * index / segments;
        circle_stream << "v " << std::cos(angle) << ' ' << std::sin(angle) << " 0\n";
    }
    for (int index = 0; index < segments; ++index)
        circle_stream << "f 1 " << index + 2 << ' ' << (index + 1) % segments + 2 << '\n';
    circle_stream.close();
    pqss_proxy_mesh::PrimitiveMeshAnalysisOptions circular = planar;
    circular.allow_frustum = true;
    circular.frustum_segments = 16;
    circular.maximum_added_volume_ratio = 0.1;
    const auto circle_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(circle, root / "circle", circular);
    if (check(circle_stats.disk_count == 1,
              "a circular planar patch must become one disk surface")) return 1;
    if (check(circle_stats.proxy_triangles == 16, "circle tessellation must honor frustum_segments")) return 1;

    const auto annulus = root / "annulus.obj";
    std::ofstream annulus_stream(annulus);
    for (int layer = 0; layer < 2; ++layer)
        for (const double radius : {1.0, 2.0})
            for (int index = 0; index < segments; ++index)
            {
                const double angle = 2.0 * 3.14159265358979323846 * index / segments;
                annulus_stream << "v " << radius * std::cos(angle) << ' '
                               << radius * std::sin(angle) << ' ' << -layer << '\n';
            }
    for (int layer = 0; layer < 2; ++layer)
    {
        const int offset = layer * 2 * segments;
        for (int index = 0; index < segments; ++index)
        {
            const int next = (index + 1) % segments;
            annulus_stream << "f " << offset + index + 1 << ' '
                           << offset + segments + index + 1 << ' '
                           << offset + segments + next + 1 << '\n';
            annulus_stream << "f " << offset + index + 1 << ' '
                           << offset + segments + next + 1 << ' '
                           << offset + next + 1 << '\n';
        }
    }
    annulus_stream.close();
    auto preserve_annulus = circular;
    preserve_annulus.maximum_cavity_added_volume_ratio = 0.0;
    const auto annulus_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        annulus, root / "annulus-preserved", preserve_annulus);
    if (check(annulus_stats.annulus_count == 2,
              "a circular inner boundary over budget must remain an annulus")) return 1;
    auto fill_annulus = circular;
    fill_annulus.maximum_cavity_added_volume_ratio = 0.4;
    const auto filled_annulus_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        annulus, root / "annulus-filled", fill_annulus);
    if (check(filled_annulus_stats.disk_count == 2 &&
              filled_annulus_stats.filled_planar_holes == 2,
              "the same circular holes must become disks when the global volume budget allows it")) return 1;

    const auto cylinder = root / "cylinder.obj";
    std::ofstream cylinder_stream(cylinder);
    for (int end = 0; end < 2; ++end)
        for (int index = 0; index < segments; ++index)
        {
            const double angle = 2.0 * 3.14159265358979323846 * index / segments;
            cylinder_stream << "v " << std::cos(angle) << ' ' << std::sin(angle)
                            << ' ' << 4.0 * end << '\n';
        }
    for (int index = 0; index < segments; ++index)
    {
        const int next = (index + 1) % segments;
        cylinder_stream << "f " << index + 1 << ' ' << next + 1 << ' '
                        << segments + next + 1 << '\n';
        cylinder_stream << "f " << index + 1 << ' ' << segments + next + 1 << ' '
                        << segments + index + 1 << '\n';
    }
    cylinder_stream.close();
    const auto cylinder_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        cylinder, root / "cylinder", circular);
    if (check(cylinder_stats.cylindrical_band_count == 1,
              "a certified cylindrical side must become one cylindrical band")) return 1;

    auto unified_lateral = circular;
    unified_lateral.use_unified_candidate_optimizer = true;
    unified_lateral.maximum_open_error_distance = 1.0e6;
    const auto unified_capless_cylinder_stats =
        pqss_proxy_mesh::analyzePrimitiveMeshObj(
            cylinder, root / "unified-capless-cylinder", unified_lateral);
    if (check(unified_capless_cylinder_stats.unified_lateral_sweep_certified == 1 &&
                  unified_capless_cylinder_stats.unified_final_analytic_choices > 0 &&
                  unified_capless_cylinder_stats.unified_selected_workload <
                      unified_capless_cylinder_stats.unified_best_analytic_workload &&
                  unified_capless_cylinder_stats.unified_selected_analytic_faces == 0,
              "a certified analytic candidate must compete by final workload instead of being forced into the output"))
        return 1;

    const auto unified_cylinder = root / "unified-capped-cylinder.obj";
    std::ofstream unified_cylinder_stream(unified_cylinder);
    for (int end = 0; end < 2; ++end)
        for (int index = 0; index < segments; ++index)
        {
            const double angle =
                2.0 * 3.14159265358979323846 * index / segments;
            unified_cylinder_stream << "v " << std::cos(angle) << ' '
                                    << std::sin(angle) << ' '
                                    << 6.0 * end << '\n';
        }
    unified_cylinder_stream << "v 0 0 0\nv 0 0 6\n";
    const int base_center = 2 * segments + 1;
    const int top_center = base_center + 1;
    for (int index = 0; index < segments; ++index)
    {
        const int next = (index + 1) % segments;
        unified_cylinder_stream << "f " << index + 1 << ' ' << next + 1
                                << ' ' << segments + next + 1 << '\n';
        unified_cylinder_stream << "f " << index + 1 << ' '
                                << segments + next + 1 << ' '
                                << segments + index + 1 << '\n';
        unified_cylinder_stream << "f " << base_center << ' ' << next + 1
                                << ' ' << index + 1 << '\n';
        unified_cylinder_stream << "f " << top_center << ' '
                                << segments + index + 1 << ' '
                                << segments + next + 1 << '\n';
    }
    unified_cylinder_stream.close();
    auto unified_circular = unified_lateral;
    const auto unified_cylinder_stats =
        pqss_proxy_mesh::analyzePrimitiveMeshObj(
            unified_cylinder, root / "unified-capped-cylinder",
            unified_circular);
    if (check(unified_cylinder_stats.unified_lateral_sweep_certified == 1 &&
                  unified_cylinder_stats.unified_final_analytic_choices > 0 &&
                  unified_cylinder_stats.unified_selected_workload <
                      unified_cylinder_stats.unified_best_analytic_workload,
              "a significant long analytic sweep must remain a candidate without overriding a cheaper feasible enclosure"))
        return 1;
    if (check(objIsClosedTriangleSoup(
                  root / "unified-capped-cylinder" / "proxy.obj"),
              "a triangulated analytic sweep must connect adjacent ring segments"))
        return 1;
    if (check(unified_cylinder_stats.coverage_assigned_source_faces ==
                  unified_cylinder_stats.source_triangles &&
              unified_cylinder_stats.coverage_enclosure_source_faces ==
                  unified_cylinder_stats.source_triangles &&
              unified_cylinder_stats.coverage_unassigned_source_faces == 0 &&
              unified_cylinder_stats.coverage_failed_source_faces == 0,
              "the final analytic enclosure must certify every source triangle before export"))
        return 1;

    const auto unified_short_cylinder = root / "unified-short-cylinder.obj";
    std::ofstream unified_short_cylinder_stream(unified_short_cylinder);
    for (int end = 0; end < 2; ++end)
        for (int index = 0; index < segments; ++index)
        {
            const double angle =
                2.0 * 3.14159265358979323846 * index / segments;
            unified_short_cylinder_stream << "v " << std::cos(angle) << ' '
                                          << std::sin(angle) << ' '
                                          << 0.25 * end << '\n';
        }
    unified_short_cylinder_stream << "v 0 0 0\nv 0 0 0.25\n";
    for (int index = 0; index < segments; ++index)
    {
        const int next = (index + 1) % segments;
        unified_short_cylinder_stream << "f " << index + 1 << ' '
            << next + 1 << ' ' << segments + next + 1 << '\n';
        unified_short_cylinder_stream << "f " << index + 1 << ' '
            << segments + next + 1 << ' ' << segments + index + 1 << '\n';
        unified_short_cylinder_stream << "f " << base_center << ' '
            << next + 1 << ' ' << index + 1 << '\n';
        unified_short_cylinder_stream << "f " << top_center << ' '
            << segments + index + 1 << ' ' << segments + next + 1 << '\n';
    }
    unified_short_cylinder_stream.close();
    const auto unified_short_cylinder_stats =
        pqss_proxy_mesh::analyzePrimitiveMeshObj(
            unified_short_cylinder, root / "unified-short-cylinder",
            unified_circular);
    if (check(unified_short_cylinder_stats.unified_lateral_sweep_certified == 1 &&
                  unified_short_cylinder_stats.unified_final_analytic_choices > 0 &&
                  unified_short_cylinder_stats.unified_selected_workload <= 12,
              "a certified model-significant disk must compete with a cheaper feasible enclosure"))
        return 1;

    const auto short_drum_with_attachments =
        root / "short-drum-with-attachments.obj";
    std::ofstream short_drum_with_attachments_stream(
        short_drum_with_attachments);
    short_drum_with_attachments_stream << readText(unified_short_cylinder);
    const auto append_box = [&](const double x0, const double y0,
                                const double z0, const double x1,
                                const double y1, const double z1,
                                const int first_vertex)
    {
        short_drum_with_attachments_stream
            << "v " << x0 << ' ' << y0 << ' ' << z0 << '\n'
            << "v " << x1 << ' ' << y0 << ' ' << z0 << '\n'
            << "v " << x1 << ' ' << y1 << ' ' << z0 << '\n'
            << "v " << x0 << ' ' << y1 << ' ' << z0 << '\n'
            << "v " << x0 << ' ' << y0 << ' ' << z1 << '\n'
            << "v " << x1 << ' ' << y0 << ' ' << z1 << '\n'
            << "v " << x1 << ' ' << y1 << ' ' << z1 << '\n'
            << "v " << x0 << ' ' << y1 << ' ' << z1 << '\n';
        const int a = first_vertex;
        short_drum_with_attachments_stream
            << "f " << a << ' ' << a + 1 << ' ' << a + 2 << '\n'
            << "f " << a << ' ' << a + 2 << ' ' << a + 3 << '\n'
            << "f " << a + 4 << ' ' << a + 6 << ' ' << a + 5 << '\n'
            << "f " << a + 4 << ' ' << a + 7 << ' ' << a + 6 << '\n'
            << "f " << a << ' ' << a + 4 << ' ' << a + 5 << '\n'
            << "f " << a << ' ' << a + 5 << ' ' << a + 1 << '\n'
            << "f " << a + 1 << ' ' << a + 5 << ' ' << a + 6 << '\n'
            << "f " << a + 1 << ' ' << a + 6 << ' ' << a + 2 << '\n'
            << "f " << a + 2 << ' ' << a + 6 << ' ' << a + 7 << '\n'
            << "f " << a + 2 << ' ' << a + 7 << ' ' << a + 3 << '\n'
            << "f " << a + 3 << ' ' << a + 7 << ' ' << a + 4 << '\n'
            << "f " << a + 3 << ' ' << a + 4 << ' ' << a << '\n';
    };
    append_box(1.4, -0.25, -0.125, 1.9, 0.25, 0.375, 27);
    append_box(-1.9, -0.25, -0.125, -1.4, 0.25, 0.375, 35);
    short_drum_with_attachments_stream.close();
    const auto short_drum_with_attachments_stats =
        pqss_proxy_mesh::analyzePrimitiveMeshObj(
            short_drum_with_attachments,
            root / "short-drum-with-attachments", unified_circular);
    if (check(
            short_drum_with_attachments_stats.unified_lateral_sweep_certified > 0 &&
                short_drum_with_attachments_stats
                        .unified_lateral_significant_candidates > 0 &&
                short_drum_with_attachments_stats
                        .unified_final_analytic_choices > 0,
            "a short analytic body with separate attachments must remain in the final model-level competition"))
        return 1;

    const auto dense_internal_cylinder = root / "dense-internal-cylinder.obj";
    std::ofstream dense_internal_stream(dense_internal_cylinder);
    dense_internal_stream << readText(unified_short_cylinder)
        << "v -5 -5 -5\nv 5 -5 -5\nv 5 5 -5\nv -5 5 -5\n"
        << "v -5 -5 5\nv 5 -5 5\nv 5 5 5\nv -5 5 5\n"
        << "f 27 28 29\nf 27 29 30\nf 31 33 32\nf 31 34 33\n"
        << "f 27 31 32\nf 27 32 28\nf 28 32 33\nf 28 33 29\n"
        << "f 29 33 34\nf 29 34 30\nf 30 34 31\nf 30 31 27\n";
    dense_internal_stream.close();
    const auto dense_internal_stats =
        pqss_proxy_mesh::analyzePrimitiveMeshObj(
            dense_internal_cylinder, root / "dense-internal-cylinder",
            unified_circular);
    if (check(dense_internal_stats.cylindrical_band_count == 0 &&
                  dense_internal_stats.proxy_triangles == 12,
              "a densely tessellated internal cylinder must not outweigh its enclosing shell"))
        return 1;
    if (check(dense_internal_stats.containment_validation_passed &&
                  dense_internal_stats.open_max_distance <= 1.0e-8,
              "filling geometry inside a closed source shell must not count as open-region error"))
        return 1;

    const auto cone = root / "cone.obj";
    std::ofstream cone_stream(cone);
    for (int end = 0; end < 2; ++end)
        for (int index = 0; index < segments; ++index)
        {
            const double radius = end == 0 ? 1.0 : 0.35;
            const double angle = 2.0 * 3.14159265358979323846 * index / segments;
            cone_stream << "v " << radius * std::cos(angle) << ' '
                        << radius * std::sin(angle) << ' ' << 4.0 * end << '\n';
        }
    for (int index = 0; index < segments; ++index)
    {
        const int next = (index + 1) % segments;
        cone_stream << "f " << index + 1 << ' ' << next + 1 << ' '
                    << segments + next + 1 << '\n';
        cone_stream << "f " << index + 1 << ' ' << segments + next + 1 << ' '
                    << segments + index + 1 << '\n';
    }
    cone_stream.close();
    const auto cone_stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(
        cone, root / "cone", circular);
    if (check(cone_stats.conical_band_count == 1,
              "a certified tapered side must become one conical band")) return 1;

    std::filesystem::remove_all(root);
    return 0;
}
