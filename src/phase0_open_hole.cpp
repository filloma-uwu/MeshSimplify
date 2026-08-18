#include "pqss_proxy_mesh/phase0_open_hole.hpp"

#include "pqss_proxy_mesh/topology_fill.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <numbers>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace pqss_proxy_mesh
{
namespace
{

struct Vec2
{
    double x = 0.0;
    double y = 0.0;
};

struct Vec3
{
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

Vec3 operator+(const Vec3 a, const Vec3 b)
{
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

Vec3 operator-(const Vec3 a, const Vec3 b)
{
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

Vec3 operator*(const Vec3 a, const double scale)
{
    return {a.x * scale, a.y * scale, a.z * scale};
}

double dot(const Vec3 a, const Vec3 b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vec3 cross(const Vec3 a, const Vec3 b)
{
    return {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x};
}

double norm(const Vec3 value)
{
    return std::sqrt(dot(value, value));
}

Vec3 normalized(const Vec3 value)
{
    const double length = norm(value);
    if (!(length > 0.0)) return {};
    return value * (1.0 / length);
}

Vec3 toVec(const Position3 point)
{
    return {point.x, point.y, point.z};
}

double orient2d(const Vec2 a, const Vec2 b, const Vec2 c)
{
    return (b.x - a.x) * (c.y - a.y) -
           (b.y - a.y) * (c.x - a.x);
}

bool pointInTriangle(
    const Vec2 point, const Vec2 a, const Vec2 b, const Vec2 c,
    const double orientation, const double tolerance)
{
    const double ab = orient2d(a, b, point) * orientation;
    const double bc = orient2d(b, c, point) * orientation;
    const double ca = orient2d(c, a, point) * orientation;
    return ab >= -tolerance && bc >= -tolerance && ca >= -tolerance;
}

int orientationSign(const Vec2 a, const Vec2 b, const Vec2 c, const double tolerance)
{
    const double value = orient2d(a, b, c);
    if (value > tolerance) return 1;
    if (value < -tolerance) return -1;
    return 0;
}

bool onSegment(const Vec2 a, const Vec2 b, const Vec2 point, const double tolerance)
{
    return point.x >= std::min(a.x, b.x) - tolerance &&
           point.x <= std::max(a.x, b.x) + tolerance &&
           point.y >= std::min(a.y, b.y) - tolerance &&
           point.y <= std::max(a.y, b.y) + tolerance &&
           std::abs(orient2d(a, b, point)) <= tolerance;
}

bool segmentsIntersect(
    const Vec2 a, const Vec2 b, const Vec2 c, const Vec2 d,
    const double tolerance)
{
    const int o1 = orientationSign(a, b, c, tolerance);
    const int o2 = orientationSign(a, b, d, tolerance);
    const int o3 = orientationSign(c, d, a, tolerance);
    const int o4 = orientationSign(c, d, b, tolerance);
    if (o1 != 0 && o2 != 0 && o3 != 0 && o4 != 0)
        return o1 != o2 && o3 != o4;
    return (o1 == 0 && onSegment(a, b, c, tolerance)) ||
           (o2 == 0 && onSegment(a, b, d, tolerance)) ||
           (o3 == 0 && onSegment(c, d, a, tolerance)) ||
           (o4 == 0 && onSegment(c, d, b, tolerance));
}

double signedArea(const std::vector<Vec2>& polygon)
{
    double twice_area = 0.0;
    for (std::size_t index = 0; index < polygon.size(); ++index)
    {
        const Vec2 a = polygon[index];
        const Vec2 b = polygon[(index + 1) % polygon.size()];
        twice_area += a.x * b.y - a.y * b.x;
    }
    return 0.5 * twice_area;
}

bool isSimplePolygon(const std::vector<Vec2>& polygon, const double tolerance)
{
    const std::size_t count = polygon.size();
    for (std::size_t first = 0; first < count; ++first)
    {
        const std::size_t first_next = (first + 1) % count;
        for (std::size_t second = first + 1; second < count; ++second)
        {
            const std::size_t second_next = (second + 1) % count;
            if (first == second || first_next == second ||
                second_next == first)
                continue;
            if (segmentsIntersect(
                    polygon[first], polygon[first_next],
                    polygon[second], polygon[second_next], tolerance))
                return false;
        }
    }
    return true;
}

std::vector<TriangleIndices> triangulateSimplePolygon(
    const std::vector<Vec2>& polygon,
    const std::vector<std::uint32_t>& vertex_ids,
    const double tolerance)
{
    if (polygon.size() != vertex_ids.size() || polygon.size() < 3) return {};
    const double area = signedArea(polygon);
    if (std::abs(area) <= tolerance) return {};
    const double winding = area > 0.0 ? 1.0 : -1.0;
    std::vector<std::size_t> remaining(polygon.size());
    for (std::size_t index = 0; index < remaining.size(); ++index)
        remaining[index] = index;
    std::vector<TriangleIndices> result;
    result.reserve(polygon.size() - 2);
    while (remaining.size() > 3)
    {
        bool clipped = false;
        for (std::size_t current = 0; current < remaining.size(); ++current)
        {
            const std::size_t previous = remaining[
                (current + remaining.size() - 1) % remaining.size()];
            const std::size_t center = remaining[current];
            const std::size_t next = remaining[(current + 1) % remaining.size()];
            if (orient2d(polygon[previous], polygon[center], polygon[next]) *
                    winding <= tolerance)
                continue;
            bool contains_vertex = false;
            for (const std::size_t candidate : remaining)
            {
                if (candidate == previous || candidate == center || candidate == next)
                    continue;
                if (pointInTriangle(
                        polygon[candidate], polygon[previous], polygon[center],
                        polygon[next], winding, tolerance))
                {
                    contains_vertex = true;
                    break;
                }
            }
            if (contains_vertex) continue;
            result.push_back({
                vertex_ids[previous], vertex_ids[center], vertex_ids[next]});
            remaining.erase(remaining.begin() + static_cast<std::ptrdiff_t>(current));
            clipped = true;
            break;
        }
        if (!clipped) return {};
    }
    result.push_back({
        vertex_ids[remaining[0]], vertex_ids[remaining[1]],
        vertex_ids[remaining[2]]});
    return result;
}

using PositionKey = std::array<double, 3>;
using EdgeKey = std::array<std::uint32_t, 2>;

PositionKey positionKey(const Position3 point)
{
    return {point.x, point.y, point.z};
}

struct BoundaryEdge
{
    std::uint32_t first = 0;
    std::uint32_t second = 0;
    std::uint32_t directed_first = 0;
    std::uint32_t directed_second = 0;
};

void writeCombinedObj(
    const std::filesystem::path& path,
    const Phase0OpenHoleResult& result)
{
    std::ofstream stream(path);
    if (!stream) throw std::runtime_error("failed to create phase-0 OBJ: " + path.string());
    stream << std::setprecision(17) << "o phase0_original_obj_open_hole_caps\n";
    for (const Position3 point : result.covered_mesh.vertices)
        stream << "v " << point.x << ' ' << point.y << ' ' << point.z << '\n';
    stream << "g phase0_source\n";
    for (std::size_t face = 0; face < result.stats.source_triangles; ++face)
    {
        const TriangleIndices triangle = result.covered_mesh.triangles[face];
        stream << "f " << triangle[0] + 1 << ' ' << triangle[1] + 1 << ' '
               << triangle[2] + 1 << '\n';
    }
    for (std::size_t cap = 0; cap < result.cap_records.size(); ++cap)
    {
        stream << "g phase0_cap_" << cap << '\n';
        const Phase0OpenHoleCap& record = result.cap_records[cap];
        for (std::size_t local = 0; local < record.cap_triangles; ++local)
        {
            const TriangleIndices triangle = result.covered_mesh.triangles[
                record.first_cap_triangle + local];
            stream << "f " << triangle[0] + 1 << ' ' << triangle[1] + 1 << ' '
                   << triangle[2] + 1 << '\n';
        }
    }
}

void writeCapsObj(
    const std::filesystem::path& path,
    const Phase0OpenHoleResult& result)
{
    std::ofstream stream(path);
    if (!stream) throw std::runtime_error("failed to create phase-0 cap OBJ: " + path.string());
    stream << std::setprecision(17) << "o phase0_open_hole_caps\n";
    for (const Position3 point : result.covered_mesh.vertices)
        stream << "v " << point.x << ' ' << point.y << ' ' << point.z << '\n';
    for (std::size_t cap = 0; cap < result.cap_records.size(); ++cap)
    {
        stream << "g phase0_cap_" << cap << '\n';
        const Phase0OpenHoleCap& record = result.cap_records[cap];
        for (std::size_t local = 0; local < record.cap_triangles; ++local)
        {
            const TriangleIndices triangle = result.covered_mesh.triangles[
                record.first_cap_triangle + local];
            stream << "f " << triangle[0] + 1 << ' ' << triangle[1] + 1 << ' '
                   << triangle[2] + 1 << '\n';
        }
    }
}

void writeStatsJson(
    const std::filesystem::path& path,
    const Phase0OpenHoleResult& result)
{
    std::ofstream stream(path);
    if (!stream) throw std::runtime_error("failed to create phase-0 statistics: " + path.string());
    const auto& stats = result.stats;
    stream << std::setprecision(17)
           << "{\n"
           << "  \"source_vertices\":" << stats.source_vertices << ",\n"
           << "  \"source_triangles\":" << stats.source_triangles << ",\n"
           << "  \"analysis_vertices\":" << stats.analysis_vertices << ",\n"
           << "  \"analysis_triangles\":" << stats.analysis_triangles << ",\n"
           << "  \"reality_boundary_edges\":" << stats.reality_boundary_edges << ",\n"
           << "  \"closed_boundary_loops\":" << stats.closed_boundary_loops << ",\n"
           << "  \"rejected_branched_boundaries\":" << stats.rejected_branched_boundaries << ",\n"
           << "  \"nonplanar_loops\":" << stats.nonplanar_loops << ",\n"
           << "  \"rejected_nonsimple_loops\":" << stats.rejected_nonsimple_loops << ",\n"
           << "  \"rejected_triangulation_failures\":" << stats.rejected_triangulation_failures << ",\n"
           << "  \"capped_loops\":" << stats.capped_loops << ",\n"
           << "  \"cap_triangles\":" << stats.cap_triangles << ",\n"
           << "  \"elapsed_seconds\":" << stats.elapsed_seconds << ",\n"
           << "  \"caps\":[\n";
    for (std::size_t index = 0; index < result.cap_records.size(); ++index)
    {
        const auto& cap = result.cap_records[index];
        stream << "    {\"id\":" << index
               << ",\"boundary_vertices\":" << cap.boundary_vertices
               << ",\"cap_triangles\":" << cap.cap_triangles
               << ",\"area\":" << cap.area
               << ",\"maximum_planarity_error\":"
               << cap.maximum_planarity_error << "}"
               << (index + 1 == result.cap_records.size() ? "\n" : ",\n");
    }
    stream << "  ]\n}\n";
}

} // namespace

Phase0OpenHoleResult coverPlanarOpenHoles(
    const MeshModel& source, const Phase0OpenHoleOptions& options)
{
    const auto started = std::chrono::steady_clock::now();
    if (!std::isfinite(options.projection_relative_tolerance) ||
        options.projection_relative_tolerance <= 0.0)
        throw std::invalid_argument("phase-0 projection tolerance must be positive and finite");

    AnalysisHalfedgeStats analysis_stats;
    const OrientedSurfaceMesh analysis = buildAnalysisHalfedgeMesh(source, &analysis_stats);
    Phase0OpenHoleResult result;
    result.covered_mesh = source;
    result.covered_mesh.name = source.name + "_phase0_open_holes_covered";
    result.caps.name = source.name + "_phase0_caps";
    result.stats.source_vertices = source.vertices.size();
    result.stats.source_triangles = source.triangles.size();
    result.stats.analysis_vertices = analysis.geometry.vertices.size();
    result.stats.analysis_triangles = analysis.geometry.triangles.size();

    Vec3 lower{
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity()};
    Vec3 upper{-lower.x, -lower.y, -lower.z};
    for (const Position3 point : source.vertices)
    {
        lower.x = std::min(lower.x, point.x);
        lower.y = std::min(lower.y, point.y);
        lower.z = std::min(lower.z, point.z);
        upper.x = std::max(upper.x, point.x);
        upper.y = std::max(upper.y, point.y);
        upper.z = std::max(upper.z, point.z);
    }
    const double scale = norm(upper - lower);
    const double planarity_diagnostic_tolerance = std::max(scale * 1.0e-7, 1.0e-12);

    std::map<PositionKey, std::uint32_t> node_by_position;
    std::vector<Position3> node_positions;
    std::vector<std::uint32_t> analysis_to_node(analysis.geometry.vertices.size());
    for (std::size_t vertex = 0; vertex < analysis.geometry.vertices.size(); ++vertex)
    {
        const Position3 point = analysis.geometry.vertices[vertex];
        const auto [iterator, inserted] = node_by_position.try_emplace(
            positionKey(point), static_cast<std::uint32_t>(node_positions.size()));
        if (inserted) node_positions.push_back(point);
        analysis_to_node[vertex] = iterator->second;
    }

    struct EdgeUse
    {
        std::uint32_t from = 0;
        std::uint32_t to = 0;
    };
    std::map<EdgeKey, std::vector<EdgeUse>> uses;
    for (const TriangleIndices triangle : analysis.geometry.triangles)
        for (std::size_t local = 0; local < 3; ++local)
        {
            const std::uint32_t from = analysis_to_node[triangle[local]];
            const std::uint32_t to = analysis_to_node[triangle[(local + 1) % 3]];
            if (from == to) continue;
            uses[{std::min(from, to), std::max(from, to)}].push_back({from, to});
        }

    std::vector<BoundaryEdge> boundary_edges;
    for (const auto& [edge, edge_uses] : uses)
    {
        if (edge_uses.size() != 1) continue;
        boundary_edges.push_back({
            edge[0], edge[1], edge_uses.front().from, edge_uses.front().to});
    }
    result.stats.reality_boundary_edges = boundary_edges.size();
    std::vector<std::vector<std::uint32_t>> incident(node_positions.size());
    for (std::uint32_t edge = 0; edge < boundary_edges.size(); ++edge)
    {
        incident[boundary_edges[edge].first].push_back(edge);
        incident[boundary_edges[edge].second].push_back(edge);
    }

    std::map<PositionKey, std::uint32_t> covered_vertex_by_position;
    for (std::uint32_t vertex = 0; vertex < result.covered_mesh.vertices.size(); ++vertex)
        covered_vertex_by_position.try_emplace(
            positionKey(result.covered_mesh.vertices[vertex]), vertex);

    std::vector<std::uint8_t> visited(boundary_edges.size(), 0);
    for (std::uint32_t seed = 0; seed < boundary_edges.size(); ++seed)
    {
        if (visited[seed]) continue;
        std::vector<std::uint32_t> component_edges{seed};
        visited[seed] = 1;
        std::set<std::uint32_t> component_nodes;
        for (std::size_t cursor = 0; cursor < component_edges.size(); ++cursor)
        {
            const BoundaryEdge edge = boundary_edges[component_edges[cursor]];
            for (const std::uint32_t node : {edge.first, edge.second})
            {
                component_nodes.insert(node);
                for (const std::uint32_t neighbor : incident[node])
                    if (!visited[neighbor])
                    {
                        visited[neighbor] = 1;
                        component_edges.push_back(neighbor);
                    }
            }
        }
        bool cycle = component_nodes.size() >= 3 &&
                     component_edges.size() == component_nodes.size();
        for (const std::uint32_t node : component_nodes)
            cycle = cycle && incident[node].size() == 2;
        if (!cycle)
        {
            ++result.stats.rejected_branched_boundaries;
            continue;
        }
        ++result.stats.closed_boundary_loops;

        std::vector<std::uint32_t> loop;
        loop.reserve(component_nodes.size());
        std::uint32_t previous_edge = std::numeric_limits<std::uint32_t>::max();
        std::uint32_t current_node = boundary_edges[component_edges.front()].directed_first;
        const std::uint32_t first_node = current_node;
        for (;;)
        {
            loop.push_back(current_node);
            const auto& adjacent = incident[current_node];
            const std::uint32_t next_edge = adjacent[0] == previous_edge
                ? adjacent[1] : adjacent[0];
            const BoundaryEdge edge = boundary_edges[next_edge];
            const std::uint32_t next_node = edge.first == current_node
                ? edge.second : edge.first;
            previous_edge = next_edge;
            current_node = next_node;
            if (current_node == first_node) break;
            if (loop.size() > component_nodes.size())
                throw std::runtime_error("phase-0 boundary cycle traversal overflow");
        }

        Vec3 center{};
        for (const std::uint32_t node : loop)
            center = center + toVec(node_positions[node]);
        center = center * (1.0 / static_cast<double>(loop.size()));
        Vec3 normal{};
        for (std::size_t index = 0; index < loop.size(); ++index)
        {
            const Vec3 a = toVec(node_positions[loop[index]]) - center;
            const Vec3 b = toVec(node_positions[loop[(index + 1) % loop.size()]]) - center;
            normal = normal + cross(a, b);
        }
        const double normal_length = norm(normal);
        if (!(normal_length > planarity_diagnostic_tolerance *
                                  planarity_diagnostic_tolerance))
        {
            ++result.stats.rejected_nonsimple_loops;
            continue;
        }
        normal = normal * (1.0 / normal_length);
        double maximum_planarity_error = 0.0;
        for (const std::uint32_t node : loop)
            maximum_planarity_error = std::max(
                maximum_planarity_error,
                std::abs(dot(toVec(node_positions[node]) - center, normal)));
        if (maximum_planarity_error > planarity_diagnostic_tolerance)
            ++result.stats.nonplanar_loops;

        Vec3 reference = std::abs(normal.x) <= std::abs(normal.y) &&
                                 std::abs(normal.x) <= std::abs(normal.z)
            ? Vec3{1.0, 0.0, 0.0}
            : std::abs(normal.y) <= std::abs(normal.z)
                ? Vec3{0.0, 1.0, 0.0}
                : Vec3{0.0, 0.0, 1.0};
        const Vec3 axis_u = normalized(cross(reference, normal));
        const Vec3 axis_v = cross(normal, axis_u);
        std::vector<Vec2> polygon;
        polygon.reserve(loop.size());
        double coordinate_scale = 0.0;
        for (const std::uint32_t node : loop)
        {
            const Vec3 relative = toVec(node_positions[node]) - center;
            polygon.push_back({dot(relative, axis_u), dot(relative, axis_v)});
            coordinate_scale = std::max(
                coordinate_scale,
                std::max(std::abs(polygon.back().x), std::abs(polygon.back().y)));
        }
        const double area_tolerance = std::max(
            coordinate_scale * coordinate_scale *
                options.projection_relative_tolerance,
            1.0e-18);
        const double polygon_area = std::abs(signedArea(polygon));
        if (!(polygon_area > area_tolerance) ||
            !isSimplePolygon(polygon, area_tolerance))
        {
            ++result.stats.rejected_nonsimple_loops;
            continue;
        }

        std::size_t direction_matches = 0;
        for (std::size_t index = 0; index < loop.size(); ++index)
        {
            const std::uint32_t from = loop[index];
            const std::uint32_t to = loop[(index + 1) % loop.size()];
            const EdgeKey key{std::min(from, to), std::max(from, to)};
            const EdgeUse use = uses.at(key).front();
            direction_matches += use.from == from && use.to == to;
        }
        // The cap boundary must oppose the existing surface boundary.
        if (direction_matches * 2 > loop.size())
        {
            std::reverse(loop.begin(), loop.end());
            std::reverse(polygon.begin(), polygon.end());
        }

        std::vector<std::uint32_t> cap_vertices;
        cap_vertices.reserve(loop.size());
        for (const std::uint32_t node : loop)
        {
            const Position3 point = node_positions[node];
            const auto [iterator, inserted] = covered_vertex_by_position.try_emplace(
                positionKey(point),
                static_cast<std::uint32_t>(result.covered_mesh.vertices.size()));
            if (inserted) result.covered_mesh.vertices.push_back(point);
            cap_vertices.push_back(iterator->second);
        }
        const std::vector<TriangleIndices> triangles = triangulateSimplePolygon(
            polygon, cap_vertices, area_tolerance);
        if (triangles.size() != loop.size() - 2)
        {
            ++result.stats.rejected_triangulation_failures;
            continue;
        }
        Phase0OpenHoleCap record;
        record.boundary_vertices = loop.size();
        record.first_cap_triangle = result.covered_mesh.triangles.size();
        record.cap_triangles = triangles.size();
        record.area = polygon_area;
        record.maximum_planarity_error = maximum_planarity_error;
        result.covered_mesh.triangles.insert(
            result.covered_mesh.triangles.end(), triangles.begin(), triangles.end());
        result.cap_records.push_back(record);
        ++result.stats.capped_loops;
        result.stats.cap_triangles += triangles.size();
    }

    result.caps.vertices = result.covered_mesh.vertices;
    for (const Phase0OpenHoleCap& record : result.cap_records)
        result.caps.triangles.insert(
            result.caps.triangles.end(),
            result.covered_mesh.triangles.begin() +
                static_cast<std::ptrdiff_t>(record.first_cap_triangle),
            result.covered_mesh.triangles.begin() +
                static_cast<std::ptrdiff_t>(record.first_cap_triangle + record.cap_triangles));
    result.stats.elapsed_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    return result;
}

Phase0OpenHoleStats generatePhase0OpenHoleModel(
    const std::filesystem::path& source_obj,
    const std::filesystem::path& output_directory,
    const Phase0OpenHoleOptions& options)
{
    const MeshModel source = readTriangleSoupObj(source_obj);
    const Phase0OpenHoleResult result = coverPlanarOpenHoles(source, options);
    std::filesystem::create_directories(output_directory);
    std::filesystem::copy_file(
        source_obj, output_directory / "source.obj",
        std::filesystem::copy_options::overwrite_existing);
    writeCombinedObj(output_directory / "phase0_open_holes.obj", result);
    writeCapsObj(output_directory / "phase0_caps.obj", result);
    writeStatsJson(output_directory / "phase0_stats.json", result);

    std::ofstream metadata(output_directory / "model.json");
    if (!metadata) throw std::runtime_error("failed to create phase-0 model metadata");
    metadata << std::setprecision(17)
             << "{\n"
             << "  \"stats\":{\"model\":\"" << source_obj.filename().string()
             << "\",\"source_triangles\":" << result.stats.source_triangles
             << ",\"primitive_count\":0,\"primitive_types\":{},"
             << "\"proxy_triangles\":" << result.stats.cap_triangles
             << ",\"timings_seconds\":{\"total\":"
             << result.stats.elapsed_seconds << "},"
             << "\"phase0_open_holes\":{\"reality_boundary_edges\":"
             << result.stats.reality_boundary_edges
             << ",\"closed_boundary_loops\":"
             << result.stats.closed_boundary_loops
             << ",\"capped_loops\":" << result.stats.capped_loops
             << ",\"cap_triangles\":" << result.stats.cap_triangles << "}},\n"
             << "  \"source\":\"source.obj\",\n"
             << "  \"phase0_caps\":\"phase0_caps.obj\",\n"
             << "  \"phase0_combined\":\"phase0_open_holes.obj\",\n"
             << "  \"phase0_stats\":\"phase0_stats.json\",\n"
             << "  \"proxy_components\":[],\n"
             << "  \"viewer_stages\":[\"source\",\"phase0\"]\n"
             << "}\n";
    return result.stats;
}

} // namespace pqss_proxy_mesh
