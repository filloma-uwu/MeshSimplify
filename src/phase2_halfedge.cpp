#include "pqss_proxy_mesh/phase2_halfedge.hpp"

#include "pqss_proxy_mesh/halfedge_validation.hpp"

#include <clipper2/clipper.h>

#include <array>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace pqss_proxy_mesh
{

WorkingHalfedgeValidationReport validateWorkingHalfedgeTopology(
    const OrientedSurfaceMesh& mesh, const double relative_tolerance)
{
    (void)relative_tolerance;
    if (mesh.geometry.vertices.empty() || mesh.geometry.triangles.empty())
        throw std::runtime_error("working halfedge mesh is empty");
    if (mesh.halfedges.size() != mesh.geometry.triangles.size() * 3 ||
        mesh.face_halfedges.size() != mesh.geometry.triangles.size() ||
        mesh.vertex_halfedges.size() != mesh.geometry.vertices.size())
        throw std::runtime_error("working halfedge arrays have inconsistent counts");

    WorkingHalfedgeValidationReport report;
    report.vertices = mesh.geometry.vertices.size();
    report.triangles = mesh.geometry.triangles.size();
    report.halfedges = mesh.halfedges.size();
    std::vector<std::uint8_t> visited(mesh.geometry.triangles.size(), 0);
    for (std::uint32_t face = 0; face < mesh.geometry.triangles.size(); ++face)
    {
        const std::uint32_t first = mesh.face_halfedges[face];
        if (first >= mesh.halfedges.size())
            throw std::runtime_error("working halfedge face has an invalid edge");
        std::uint32_t edge = first;
        for (std::uint32_t local = 0; local < 3; ++local)
        {
            const SurfaceHalfedge& halfedge = mesh.halfedges[edge];
            if (halfedge.face != face || halfedge.next >= mesh.halfedges.size() ||
                halfedge.origin != mesh.geometry.triangles[face][local])
                throw std::runtime_error("working halfedge face cycle is inconsistent");
            if (halfedge.opposite == invalid_surface_index)
                ++report.boundary_halfedges;
            else
            {
                if (halfedge.opposite >= mesh.halfedges.size() ||
                    mesh.halfedges[halfedge.opposite].opposite != edge)
                    throw std::runtime_error("working halfedge opposite is asymmetric");
                const std::uint32_t destination = mesh.halfedges[halfedge.next].origin;
                const SurfaceHalfedge& opposite = mesh.halfedges[halfedge.opposite];
                const std::uint32_t opposite_destination =
                    mesh.halfedges[opposite.next].origin;
                if (opposite.origin != destination ||
                    opposite_destination != halfedge.origin)
                    throw std::runtime_error(
                        "working halfedge opposite is not the reversed geometric edge");
                if (edge < halfedge.opposite) ++report.paired_edges;
            }
            edge = halfedge.next;
        }
        if (edge != first)
            throw std::runtime_error("working halfedge face is not triangular");
    }
    for (std::uint32_t seed = 0; seed < mesh.geometry.triangles.size(); ++seed)
    {
        if (visited[seed]) continue;
        ++report.face_components;
        std::vector<std::uint32_t> queue{seed};
        visited[seed] = 1;
        for (std::size_t cursor = 0; cursor < queue.size(); ++cursor)
        {
            std::uint32_t edge = mesh.face_halfedges[queue[cursor]];
            for (int local = 0; local < 3; ++local)
            {
                const std::uint32_t opposite = mesh.halfedges[edge].opposite;
                if (opposite < mesh.halfedges.size())
                {
                    const std::uint32_t neighbor = mesh.halfedges[opposite].face;
                    if (neighbor >= mesh.geometry.triangles.size())
                        throw std::runtime_error("working halfedge opposite has an invalid face");
                    if (!visited[neighbor])
                    {
                        visited[neighbor] = 1;
                        queue.push_back(neighbor);
                    }
                }
                edge = mesh.halfedges[edge].next;
            }
        }
    }
    return report;
}

namespace
{

struct Vec3
{
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

Vec3 toVec(const Position3 point) { return {point.x, point.y, point.z}; }
Vec3 operator-(const Vec3 first, const Vec3 second)
{ return {first.x - second.x, first.y - second.y, first.z - second.z}; }
Vec3 operator+(const Vec3 first, const Vec3 second)
{ return {first.x + second.x, first.y + second.y, first.z + second.z}; }
Vec3 operator*(const Vec3 point, const double scale)
{ return {point.x * scale, point.y * scale, point.z * scale}; }
double dot(const Vec3 first, const Vec3 second)
{ return first.x * second.x + first.y * second.y + first.z * second.z; }
Vec3 cross(const Vec3 first, const Vec3 second)
{
    return {first.y * second.z - first.z * second.y,
            first.z * second.x - first.x * second.z,
            first.x * second.y - first.y * second.x};
}
double norm(const Vec3 point) { return std::sqrt(dot(point, point)); }
Vec3 normalized(const Vec3 point)
{
    const double length = norm(point);
    return length > 0.0 ? point * (1.0 / length) : Vec3{};
}
Position3 toPosition(const Vec3 point) { return {point.x, point.y, point.z}; }

struct EdgeKey
{
    std::uint32_t first = 0;
    std::uint32_t second = 0;
    auto operator<=>(const EdgeKey&) const = default;
};

struct EdgeUse
{
    std::uint32_t face = 0;
    bool forward = false;
};

struct LineKey
{
    std::array<std::int64_t, 6> value{};
    auto operator<=>(const LineKey&) const = default;
};

LineKey lineKey(const Position3 first, const Position3 second, const double tolerance)
{
    Vec3 direction = toVec(second) - toVec(first);
    const double length = norm(direction);
    if (length > 0.0) direction = direction * (1.0 / length);
    if (direction.x < 0.0 ||
        (direction.x == 0.0 && direction.y < 0.0) ||
        (direction.x == 0.0 && direction.y == 0.0 && direction.z < 0.0))
        direction = direction * -1.0;
    const Vec3 moment = cross(direction, toVec(first));
    return {{{static_cast<std::int64_t>(std::llround(direction.x / 1.0e-10)),
              static_cast<std::int64_t>(std::llround(direction.y / 1.0e-10)),
              static_cast<std::int64_t>(std::llround(direction.z / 1.0e-10)),
              static_cast<std::int64_t>(std::llround(moment.x / tolerance)),
              static_cast<std::int64_t>(std::llround(moment.y / tolerance)),
              static_cast<std::int64_t>(std::llround(moment.z / tolerance))}}};
}

bool pointStrictlyInsideSegment(
    const Position3 point, const Position3 first, const Position3 second,
    const double tolerance)
{
    const Vec3 a = toVec(first);
    const Vec3 b = toVec(second);
    const Vec3 p = toVec(point);
    const Vec3 edge = b - a;
    const double length = norm(edge);
    if (length <= tolerance) return false;
    const double parameter = dot(p - a, edge) / (length * length);
    if (parameter <= tolerance / length || parameter >= 1.0 - tolerance / length)
        return false;
    return norm(cross(p - a, edge)) / length <= tolerance;
}

void splitUnmatchedEdgesAtCollinearVertices(MeshModel& mesh, const double tolerance)
{
    for (int pass = 0; pass < 8; ++pass)
    {
        std::map<EdgeKey, std::vector<std::uint32_t>> uses;
        for (std::uint32_t face = 0; face < mesh.triangles.size(); ++face)
        {
            const TriangleIndices triangle = mesh.triangles[face];
            for (int local = 0; local < 3; ++local)
            {
                const std::uint32_t first = triangle[local];
                const std::uint32_t second = triangle[(local + 1) % 3];
                uses[{std::min(first, second), std::max(first, second)}].push_back(face);
            }
        }
        std::map<LineKey, std::vector<std::uint32_t>> line_vertices;
        std::set<EdgeKey> unmatched;
        for (const auto& [edge, edge_uses] : uses)
        {
            if (edge_uses.size() == 2) continue;
            unmatched.insert(edge);
            const LineKey key = lineKey(
                mesh.vertices[edge.first], mesh.vertices[edge.second], tolerance);
            line_vertices[key].push_back(edge.first);
            line_vertices[key].push_back(edge.second);
        }
        for (auto& [key, vertices] : line_vertices)
        {
            (void)key;
            std::sort(vertices.begin(), vertices.end());
            vertices.erase(std::unique(vertices.begin(), vertices.end()), vertices.end());
        }

        bool changed = false;
        std::vector<TriangleIndices> rebuilt;
        rebuilt.reserve(mesh.triangles.size());
        for (const TriangleIndices triangle : mesh.triangles)
        {
            bool split = false;
            for (int local = 0; local < 3 && !split; ++local)
            {
                const std::uint32_t first = triangle[local];
                const std::uint32_t second = triangle[(local + 1) % 3];
                const std::uint32_t third = triangle[(local + 2) % 3];
                const EdgeKey edge{std::min(first, second), std::max(first, second)};
                if (!unmatched.contains(edge)) continue;
                const LineKey key = lineKey(
                    mesh.vertices[first], mesh.vertices[second], tolerance);
                const auto candidates = line_vertices.find(key);
                if (candidates == line_vertices.end()) continue;
                std::vector<std::pair<double, std::uint32_t>> cuts;
                const Vec3 a = toVec(mesh.vertices[first]);
                const Vec3 b = toVec(mesh.vertices[second]);
                const Vec3 segment = b - a;
                const double length2 = dot(segment, segment);
                for (const std::uint32_t vertex : candidates->second)
                {
                    if (vertex == first || vertex == second) continue;
                    if (!pointStrictlyInsideSegment(
                            mesh.vertices[vertex], mesh.vertices[first],
                            mesh.vertices[second], tolerance))
                        continue;
                    cuts.push_back({dot(toVec(mesh.vertices[vertex]) - a, segment) / length2,
                                    vertex});
                }
                if (cuts.empty()) continue;
                std::sort(cuts.begin(), cuts.end());
                std::uint32_t previous = first;
                for (const auto [parameter, vertex] : cuts)
                {
                    (void)parameter;
                    if (previous != vertex && vertex != third && previous != third)
                        rebuilt.push_back({previous, vertex, third});
                    previous = vertex;
                }
                if (previous != second && second != third && previous != third)
                    rebuilt.push_back({previous, second, third});
                split = true;
                changed = true;
            }
            if (!split) rebuilt.push_back(triangle);
        }
        mesh.triangles = std::move(rebuilt);
        if (!changed) return;
    }
}

void orientTriangleSoupConsistently(MeshModel& mesh)
{
    std::map<EdgeKey, std::vector<EdgeUse>> uses;
    for (std::uint32_t face = 0; face < mesh.triangles.size(); ++face)
    {
        const TriangleIndices triangle = mesh.triangles[face];
        for (int local = 0; local < 3; ++local)
        {
            const std::uint32_t first = triangle[local];
            const std::uint32_t second = triangle[(local + 1) % 3];
            const EdgeKey key{std::min(first, second), std::max(first, second)};
            uses[key].push_back({face, first == key.first});
        }
    }

    std::vector<std::vector<std::pair<std::uint32_t, bool>>> adjacency(
        mesh.triangles.size());
    for (const auto& [key, edge_uses] : uses)
    {
        (void)key;
        if (edge_uses.size() != 2) continue;
        const EdgeUse first = edge_uses[0];
        const EdgeUse second = edge_uses[1];
        const bool different_flip = first.forward == second.forward;
        adjacency[first.face].push_back({second.face, different_flip});
        adjacency[second.face].push_back({first.face, different_flip});
    }

    std::vector<std::uint8_t> visited(mesh.triangles.size(), 0);
    std::vector<std::uint8_t> flipped(mesh.triangles.size(), 0);
    for (std::uint32_t seed = 0; seed < mesh.triangles.size(); ++seed)
    {
        if (visited[seed]) continue;
        visited[seed] = 1;
        std::vector<std::uint32_t> queue{seed};
        for (std::size_t cursor = 0; cursor < queue.size(); ++cursor)
        {
            const std::uint32_t face = queue[cursor];
            for (const auto [neighbor, different_flip] : adjacency[face])
            {
                const bool face_flipped = flipped[face] != 0;
                const std::uint8_t neighbor_flip =
                    static_cast<std::uint8_t>(face_flipped != different_flip);
                if (!visited[neighbor])
                {
                    visited[neighbor] = 1;
                    flipped[neighbor] = neighbor_flip;
                    queue.push_back(neighbor);
                }
            }
        }
    }
    for (std::uint32_t face = 0; face < mesh.triangles.size(); ++face)
        if (flipped[face])
            std::swap(mesh.triangles[face][1], mesh.triangles[face][2]);
}

OrientedSurfaceMesh buildIndexedClosedHalfedge(MeshModel mesh)
{
    std::vector<std::uint8_t> used(mesh.vertices.size(), 0);
    for (const TriangleIndices face : mesh.triangles)
        for (const std::uint32_t vertex : face)
        {
            if (vertex >= mesh.vertices.size())
                throw std::runtime_error("phase-2 triangle has an invalid vertex");
            used[vertex] = 1;
        }
    std::vector<std::uint32_t> remap(mesh.vertices.size(), invalid_surface_index);
    MeshModel compact;
    compact.name = std::move(mesh.name);
    Vec3 lower = toVec(mesh.vertices.front());
    Vec3 upper = lower;
    for (const Position3 point : mesh.vertices)
    {
        const Vec3 value = toVec(point);
        lower.x = std::min(lower.x, value.x);
        lower.y = std::min(lower.y, value.y);
        lower.z = std::min(lower.z, value.z);
        upper.x = std::max(upper.x, value.x);
        upper.y = std::max(upper.y, value.y);
        upper.z = std::max(upper.z, value.z);
    }
    const double weld_tolerance = std::max(norm(upper - lower) * 1.0e-9, 1.0e-8);
    std::map<std::array<std::int64_t, 3>, std::uint32_t> welded_vertices;
    for (std::uint32_t vertex = 0; vertex < mesh.vertices.size(); ++vertex)
        if (used[vertex])
        {
            const Position3 point = mesh.vertices[vertex];
            const std::array<std::int64_t, 3> key{
                static_cast<std::int64_t>(std::llround(point.x / weld_tolerance)),
                static_cast<std::int64_t>(std::llround(point.y / weld_tolerance)),
                static_cast<std::int64_t>(std::llround(point.z / weld_tolerance))};
            const auto [iterator, inserted] = welded_vertices.try_emplace(
                key, static_cast<std::uint32_t>(compact.vertices.size()));
            if (inserted) compact.vertices.push_back(point);
            remap[vertex] = iterator->second;
        }
    compact.triangles.reserve(mesh.triangles.size());
    for (const TriangleIndices face : mesh.triangles)
    {
        const TriangleIndices triangle{remap[face[0]], remap[face[1]], remap[face[2]]};
        if (triangle[0] == triangle[1] || triangle[1] == triangle[2] ||
            triangle[2] == triangle[0])
            continue;
        compact.triangles.push_back(triangle);
    }
    if (compact.triangles.empty()) throw std::runtime_error("phase-2 mesh has no faces");
    splitUnmatchedEdgesAtCollinearVertices(
        compact, std::max(weld_tolerance * 1024.0, 1.0e-5));
    orientTriangleSoupConsistently(compact);

    OrientedSurfaceMesh result;
    result.geometry = std::move(compact);
    result.halfedges.resize(result.geometry.triangles.size() * 3);
    result.face_halfedges.resize(result.geometry.triangles.size());
    result.vertex_halfedges.assign(result.geometry.vertices.size(), invalid_surface_index);
    std::map<EdgeKey, std::vector<std::uint32_t>> uses;
    for (std::uint32_t face = 0; face < result.geometry.triangles.size(); ++face)
    {
        result.face_halfedges[face] = face * 3;
        const TriangleIndices triangle = result.geometry.triangles[face];
        if (triangle[0] == triangle[1] || triangle[1] == triangle[2] ||
            triangle[2] == triangle[0])
            throw std::runtime_error("phase-2 mesh has a degenerate index cycle");
        for (std::uint32_t local = 0; local < 3; ++local)
        {
            const std::uint32_t edge = face * 3 + local;
            const std::uint32_t origin = triangle[local];
            const std::uint32_t destination = triangle[(local + 1) % 3];
            result.halfedges[edge] = {origin, face, face * 3 + (local + 1) % 3,
                                      invalid_surface_index};
            if (result.vertex_halfedges[origin] == invalid_surface_index)
                result.vertex_halfedges[origin] = edge;
            uses[{std::min(origin, destination), std::max(origin, destination)}].push_back(edge);
        }
    }
    for (const auto& [key, edges] : uses)
    {
        (void)key;
        if (edges.size() != 2)
        {
            const Position3 first_point = result.geometry.vertices[key.first];
            const Position3 second_point = result.geometry.vertices[key.second];
            throw std::runtime_error(
                "phase-2 mesh is not edge-manifold and closed: edge=" +
                std::to_string(key.first) + "," + std::to_string(key.second) +
                " uses=" + std::to_string(edges.size()) +
                " first=(" + std::to_string(first_point.x) + "," +
                std::to_string(first_point.y) + "," + std::to_string(first_point.z) +
                ") second=(" + std::to_string(second_point.x) + "," +
                std::to_string(second_point.y) + "," + std::to_string(second_point.z) +
                ")");
        }
        const std::uint32_t first = edges[0];
        const std::uint32_t second = edges[1];
        const std::uint32_t first_end = result.halfedges[result.halfedges[first].next].origin;
        const std::uint32_t second_end = result.halfedges[result.halfedges[second].next].origin;
        if (result.halfedges[first].origin != second_end ||
            result.halfedges[second].origin != first_end)
            throw std::runtime_error("phase-2 mesh has an inconsistently oriented edge");
        result.halfedges[first].opposite = second;
        result.halfedges[second].opposite = first;
    }
    (void)validateClosedHalfedgeTopology(result);
    return result;
}

MeshModel canonicalizeOrientedCoplanarTriangleSoup(
    const MeshModel& source, const double relative_tolerance)
{
    Vec3 lower = toVec(source.vertices.front());
    Vec3 upper = lower;
    for (const Position3 point : source.vertices)
    {
        const Vec3 value = toVec(point);
        lower.x = std::min(lower.x, value.x);
        lower.y = std::min(lower.y, value.y);
        lower.z = std::min(lower.z, value.z);
        upper.x = std::max(upper.x, value.x);
        upper.y = std::max(upper.y, value.y);
        upper.z = std::max(upper.z, value.z);
    }
    const double scale = norm(upper - lower);
    const double distance_quantum = std::max(scale * relative_tolerance, 1.0e-12);
    constexpr double angular_quantum = 1.0e-10;
    constexpr int clipper_precision = 8;
    using PlaneKey = std::array<std::int64_t, 4>;
    struct PlaneGroup
    {
        Vec3 origin;
        Vec3 normal;
        Vec3 tangent;
        Vec3 bitangent;
        std::vector<std::uint32_t> faces;
    };
    std::map<PlaneKey, PlaneGroup> groups;

    for (std::uint32_t face_id = 0; face_id < source.triangles.size(); ++face_id)
    {
        const TriangleIndices face = source.triangles[face_id];
        const Vec3 a = toVec(source.vertices[face[0]]);
        const Vec3 b = toVec(source.vertices[face[1]]);
        const Vec3 c = toVec(source.vertices[face[2]]);
        Vec3 normal = cross(b - a, c - a);
        const double length = norm(normal);
        if (length <= std::max(scale * scale, 1.0) * 1.0e-24) continue;
        normal = normal * (1.0 / length);
        const double distance = dot(normal, a);
        PlaneKey key{
            static_cast<std::int64_t>(std::llround(normal.x / angular_quantum)),
            static_cast<std::int64_t>(std::llround(normal.y / angular_quantum)),
            static_cast<std::int64_t>(std::llround(normal.z / angular_quantum)),
            static_cast<std::int64_t>(std::llround(distance / distance_quantum))};
        auto [iterator, inserted] = groups.try_emplace(key);
        auto& group = iterator->second;
        if (inserted)
        {
            group.origin = a;
            group.normal = normal;
            const Vec3 helper = std::abs(normal.x) < 0.8 ? Vec3{1.0, 0.0, 0.0}
                                                        : Vec3{0.0, 1.0, 0.0};
            group.tangent = normalized(cross(helper, normal));
            group.bitangent = cross(normal, group.tangent);
        }
        const bool coplanar =
            std::abs(dot(group.normal, a - group.origin)) <= distance_quantum &&
            std::abs(dot(group.normal, b - group.origin)) <= distance_quantum &&
            std::abs(dot(group.normal, c - group.origin)) <= distance_quantum &&
            dot(group.normal, normal) >= 1.0 - 1.0e-12;
        if (coplanar) group.faces.push_back(face_id);
        else
        {
            key[3] ^= static_cast<std::int64_t>(face_id + 1) << 32;
            auto& unique = groups[key];
            unique.origin = a;
            unique.normal = normal;
            const Vec3 helper = std::abs(normal.x) < 0.8 ? Vec3{1.0, 0.0, 0.0}
                                                        : Vec3{0.0, 1.0, 0.0};
            unique.tangent = normalized(cross(helper, normal));
            unique.bitangent = cross(normal, unique.tangent);
            unique.faces.push_back(face_id);
        }
    }

    MeshModel result;
    result.name = source.name + "_oriented_coplanar_union";
    std::map<std::array<double, 3>, std::uint32_t> output_vertices;
    const auto appendVertex = [&](const Position3 point)
    {
        const std::array<double, 3> key{{point.x, point.y, point.z}};
        const auto [iterator, inserted] = output_vertices.try_emplace(
            key, static_cast<std::uint32_t>(result.vertices.size()));
        if (inserted) result.vertices.push_back(point);
        return iterator->second;
    };

    for (const auto& [key, group] : groups)
    {
        (void)key;
        Clipper2Lib::PathsD paths;
        for (const std::uint32_t face_id : group.faces)
        {
            Clipper2Lib::PathD path;
            for (const std::uint32_t vertex : source.triangles[face_id])
            {
                const Vec3 delta = toVec(source.vertices[vertex]) - group.origin;
                path.emplace_back(dot(delta, group.tangent), dot(delta, group.bitangent));
            }
            if (std::abs(Clipper2Lib::Area(path)) <= distance_quantum * distance_quantum)
                continue;
            if (Clipper2Lib::Area(path) < 0.0) std::reverse(path.begin(), path.end());
            paths.push_back(std::move(path));
        }
        Clipper2Lib::PathsD united;
        Clipper2Lib::ClipperD clipper(clipper_precision);
        clipper.PreserveCollinear(true);
        clipper.AddSubject(paths);
        if (!clipper.Execute(Clipper2Lib::ClipType::Union,
                             Clipper2Lib::FillRule::NonZero, united))
            united.clear();
        Clipper2Lib::PathsD triangles;
        if (united.empty() || Clipper2Lib::Triangulate(
                united, clipper_precision, triangles, false) !=
                Clipper2Lib::TriangulateResult::success)
        {
            for (const std::uint32_t face_id : group.faces)
            {
                TriangleIndices output_face{};
                for (int corner = 0; corner < 3; ++corner)
                    output_face[corner] =
                        appendVertex(source.vertices[source.triangles[face_id][corner]]);
                result.triangles.push_back(output_face);
            }
            continue;
        }
        for (const auto& triangle : triangles)
        {
            if (triangle.size() != 3) continue;
            std::array<std::uint32_t, 3> ids{};
            for (int corner = 0; corner < 3; ++corner)
            {
                const auto& point = triangle[corner];
                ids[corner] = appendVertex(toPosition(
                    group.origin + group.tangent * point.x + group.bitangent * point.y));
            }
            const Position3 a = result.vertices[ids[0]];
            const Position3 b = result.vertices[ids[1]];
            const Position3 c = result.vertices[ids[2]];
            if (dot(cross(toVec(b) - toVec(a), toVec(c) - toVec(a)), group.normal) < 0.0)
                std::swap(ids[1], ids[2]);
            result.triangles.push_back({ids[0], ids[1], ids[2]});
        }
    }
    if (result.triangles.empty())
        throw std::runtime_error("oriented coplanar union produced an empty mesh");
    return result;
}

MeshModel mergeCoplanarPatches(
    const OrientedSurfaceMesh& source, const double relative_tolerance,
    std::size_t* considered_faces)
{
    if (considered_faces) *considered_faces += source.geometry.triangles.size();
    MeshModel result = canonicalizeOrientedCoplanarTriangleSoup(
        source.geometry, relative_tolerance);
    result.name = "phase2_coplanar_union";
    orientTriangleSoupConsistently(result);
    return result;
}

struct OutwardSupportPlane
{
    Vec3 origin;
    Vec3 normal;
    double offset = 0.0;
    double area = 0.0;
    std::vector<std::uint32_t> faces;
};

struct OutwardEnvelopeResult
{
    MeshModel geometry;
    std::size_t moved_vertices = 0;
    std::size_t support_planes = 0;
    std::size_t scored_coplanar_triangles = 0;
    std::size_t selected_round = invalid_surface_index;
    std::vector<std::vector<Position3>> round_vertices;
    std::vector<std::size_t> round_cumulative_moved;
    std::vector<std::size_t> round_coplanar_scores;
};

double clamp01(const double value)
{
    return std::max(0.0, std::min(1.0, value));
}

Position3 closestPointOnTriangle(
    const Position3 point, const Position3 a, const Position3 b, const Position3 c)
{
    const Vec3 p = toVec(point);
    const Vec3 av = toVec(a);
    const Vec3 bv = toVec(b);
    const Vec3 cv = toVec(c);
    const Vec3 ab = bv - av;
    const Vec3 ac = cv - av;
    const Vec3 ap = p - av;
    const double d1 = dot(ab, ap);
    const double d2 = dot(ac, ap);
    if (d1 <= 0.0 && d2 <= 0.0) return a;
    const Vec3 bp = p - bv;
    const double d3 = dot(ab, bp);
    const double d4 = dot(ac, bp);
    if (d3 >= 0.0 && d4 <= d3) return b;
    const double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
    {
        const double v = d1 / (d1 - d3);
        return toPosition(av + ab * v);
    }
    const Vec3 cp = p - cv;
    const double d5 = dot(ab, cp);
    const double d6 = dot(ac, cp);
    if (d6 >= 0.0 && d5 <= d6) return c;
    const double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
    {
        const double w = d2 / (d2 - d6);
        return toPosition(av + ac * w);
    }
    const double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0)
    {
        const double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return toPosition(bv + (cv - bv) * w);
    }
    const double denominator = 1.0 / (va + vb + vc);
    const double v = vb * denominator;
    const double w = vc * denominator;
    return toPosition(av + ab * v + ac * w);
}

Vec3 faceNormal(const MeshModel& mesh, const std::uint32_t face_id)
{
    const TriangleIndices face = mesh.triangles[face_id];
    const Vec3 a = toVec(mesh.vertices[face[0]]);
    const Vec3 b = toVec(mesh.vertices[face[1]]);
    const Vec3 c = toVec(mesh.vertices[face[2]]);
    return normalized(cross(b - a, c - a));
}

struct LocalHalfspacePlane
{
    Vec3 normal;
    double offset = 0.0;
};

std::size_t estimateCoplanarMergeScore(
    const MeshModel& source, const double relative_tolerance)
{
    if (source.vertices.empty() || source.triangles.empty()) return 0;
    Vec3 lower = toVec(source.vertices.front());
    Vec3 upper = lower;
    for (const Position3 point : source.vertices)
    {
        const Vec3 value = toVec(point);
        lower.x = std::min(lower.x, value.x);
        lower.y = std::min(lower.y, value.y);
        lower.z = std::min(lower.z, value.z);
        upper.x = std::max(upper.x, value.x);
        upper.y = std::max(upper.y, value.y);
        upper.z = std::max(upper.z, value.z);
    }
    const double scale = std::max(norm(upper - lower), 1.0);
    const double distance_quantum = std::max(scale * relative_tolerance, 1.0e-12);
    constexpr double angular_quantum = 1.0e-10;
    using PlaneKey = std::array<std::int64_t, 4>;
    std::map<PlaneKey, std::size_t> groups;
    for (const TriangleIndices face : source.triangles)
    {
        const Vec3 a = toVec(source.vertices[face[0]]);
        const Vec3 b = toVec(source.vertices[face[1]]);
        const Vec3 c = toVec(source.vertices[face[2]]);
        Vec3 normal = cross(b - a, c - a);
        const double length = norm(normal);
        if (length <= scale * scale * 1.0e-24) continue;
        normal = normal * (1.0 / length);
        const double distance = dot(normal, a);
        PlaneKey key{
            static_cast<std::int64_t>(std::llround(normal.x / angular_quantum)),
            static_cast<std::int64_t>(std::llround(normal.y / angular_quantum)),
            static_cast<std::int64_t>(std::llround(normal.z / angular_quantum)),
            static_cast<std::int64_t>(std::llround(distance / distance_quantum))};
        ++groups[key];
    }
    std::size_t score = 0;
    for (const auto& [key, count] : groups)
    {
        (void)key;
        score += std::min<std::size_t>(count, 2);
    }
    return score;
}

OutwardEnvelopeResult makeConservativeOutwardEnvelope(
    const OrientedSurfaceMesh& source, const double maximum_distance,
    const double relative_tolerance)
{
    const MeshModel& mesh = source.geometry;
    if (mesh.vertices.empty() || mesh.triangles.empty()) return {mesh, 0, 0};

    Vec3 lower = toVec(mesh.vertices.front());
    Vec3 upper = lower;
    for (const Position3 point : mesh.vertices)
    {
        const Vec3 value = toVec(point);
        lower.x = std::min(lower.x, value.x);
        lower.y = std::min(lower.y, value.y);
        lower.z = std::min(lower.z, value.z);
        upper.x = std::max(upper.x, value.x);
        upper.y = std::max(upper.y, value.y);
        upper.z = std::max(upper.z, value.z);
    }
    const double diagonal = std::max(norm(upper - lower), 1.0);
    std::vector<std::vector<std::uint32_t>> vertex_faces(mesh.vertices.size());
    for (std::uint32_t face_id = 0; face_id < mesh.triangles.size(); ++face_id)
        for (const std::uint32_t vertex : mesh.triangles[face_id])
            if (vertex < vertex_faces.size()) vertex_faces[vertex].push_back(face_id);

    std::vector<std::vector<std::uint32_t>> vertex_neighbors(mesh.vertices.size());
    std::vector<std::vector<std::uint32_t>> face_neighbors(mesh.triangles.size());
    if (source.face_halfedges.size() == mesh.triangles.size())
    {
        for (std::uint32_t face = 0; face < mesh.triangles.size(); ++face)
        {
            std::uint32_t edge = source.face_halfedges[face];
            for (int local = 0; local < 3 && edge < source.halfedges.size(); ++local)
            {
                const std::uint32_t opposite = source.halfedges[edge].opposite;
                if (opposite < source.halfedges.size())
                {
                    const std::uint32_t neighbor = source.halfedges[opposite].face;
                    if (neighbor < mesh.triangles.size())
                        face_neighbors[face].push_back(neighbor);
                }
                const std::uint32_t origin = source.halfedges[edge].origin;
                const std::uint32_t destination =
                    source.halfedges[source.halfedges[edge].next].origin;
                if (origin < vertex_neighbors.size() && destination < vertex_neighbors.size())
                {
                    vertex_neighbors[origin].push_back(destination);
                    vertex_neighbors[destination].push_back(origin);
                }
                edge = source.halfedges[edge].next;
            }
        }
    }
    for (auto& neighbors : vertex_neighbors)
    {
        std::sort(neighbors.begin(), neighbors.end());
        neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
    }

    std::vector<Vec3> original_normals(mesh.triangles.size());
    for (std::uint32_t face = 0; face < mesh.triangles.size(); ++face)
        original_normals[face] = faceNormal(mesh, face);

    MeshModel result = mesh;
    result.name = "phase2_conservative_outward_envelope";
    MeshModel best = result;
    std::size_t best_score = std::numeric_limits<std::size_t>::max();
    std::size_t best_moved = 0;
    std::size_t selected_round = invalid_surface_index;
    std::vector<std::vector<Position3>> round_vertices;
    std::vector<std::size_t> round_cumulative_moved;
    std::vector<std::size_t> round_coplanar_scores;
    const std::uint32_t neighborhood_depth = 4;
    const std::uint32_t maximum_rounds = 32;
    const double concavity_tolerance = std::max(diagonal * 1.0e-8, 1.0e-5);
    std::size_t total_moved = 0;
    std::vector<std::uint32_t> visit_stamp(mesh.triangles.size(), 0);
    std::vector<std::uint32_t> vertex_stamp(mesh.vertices.size(), 0);
    std::uint32_t stamp = 1;

    const auto buildFacePlanes = [&](const MeshModel& current)
    {
        std::vector<LocalHalfspacePlane> planes(current.triangles.size());
        for (std::uint32_t face = 0; face < current.triangles.size(); ++face)
        {
            const TriangleIndices triangle = current.triangles[face];
            const Vec3 normal = faceNormal(current, face);
            planes[face] = {normal, dot(normal, toVec(current.vertices[triangle[0]]))};
        }
        return planes;
    };

    const auto gatherNeighborhoodFaces = [&](
        const std::uint32_t vertex, const std::uint32_t depth_limit)
    {
        std::vector<std::uint32_t> faces;
        std::vector<std::uint32_t> depths;
        if (vertex >= vertex_faces.size()) return faces;
        for (const std::uint32_t face : vertex_faces[vertex])
            if (visit_stamp[face] != stamp)
            {
                visit_stamp[face] = stamp;
                faces.push_back(face);
                depths.push_back(0);
            }
        for (std::size_t cursor = 0; cursor < faces.size(); ++cursor)
        {
            if (depths[cursor] >= depth_limit) continue;
            for (const std::uint32_t neighbor : face_neighbors[faces[cursor]])
                if (visit_stamp[neighbor] != stamp)
                {
                    visit_stamp[neighbor] = stamp;
                    faces.push_back(neighbor);
                    depths.push_back(depths[cursor] + 1);
                }
        }
        ++stamp;
        if (stamp == 0)
        {
            std::fill(visit_stamp.begin(), visit_stamp.end(), 0);
            stamp = 1;
        }
        return faces;
    };

    const auto concaveVertices = [&](
        const MeshModel& current, const std::vector<LocalHalfspacePlane>& planes)
    {
        std::vector<std::uint8_t> concave(current.vertices.size(), 0);
        std::size_t concave_vertices = 0;
        for (std::uint32_t vertex = 0; vertex < current.vertices.size(); ++vertex)
        {
            const auto faces = gatherNeighborhoodFaces(vertex, neighborhood_depth);
            double worst_violation = 0.0;
            const Vec3 point = toVec(current.vertices[vertex]);
            for (const std::uint32_t face : faces)
            {
                const double signed_distance =
                    dot(planes[face].normal, point) - planes[face].offset;
                worst_violation = std::min(worst_violation, signed_distance);
            }
            if (worst_violation < -concavity_tolerance)
            {
                concave[vertex] = 1;
                ++concave_vertices;
            }
        }
        return std::pair<std::vector<std::uint8_t>, std::size_t>{
            std::move(concave), concave_vertices};
    };

    const auto acceptsVertexMove = [&](const MeshModel& current,
                                       const std::uint32_t vertex,
                                       const Position3 candidate)
    {
        const Vec3 original_delta = toVec(candidate) - toVec(mesh.vertices[vertex]);
        if (norm(original_delta) > maximum_distance * 0.999) return false;
        for (const std::uint32_t face : vertex_faces[vertex])
        {
            const TriangleIndices triangle = current.triangles[face];
            std::array<Position3, 3> points{
                current.vertices[triangle[0]],
                current.vertices[triangle[1]],
                current.vertices[triangle[2]]};
            for (int corner = 0; corner < 3; ++corner)
                if (triangle[corner] == vertex) points[corner] = candidate;
            const Vec3 raw = cross(
                toVec(points[1]) - toVec(points[0]),
                toVec(points[2]) - toVec(points[0]));
            if (norm(raw) <= diagonal * diagonal * 1.0e-18) return false;
            if (dot(raw, original_normals[face]) <= 0.0) return false;
        }
        return true;
    };

    const auto acceptsComponentMove = [&](
        const MeshModel& current,
        const std::vector<std::pair<std::uint32_t, Position3>>& moves)
    {
        if (moves.empty()) return false;
        std::map<std::uint32_t, Position3> replacements;
        for (const auto& [vertex, point] : moves)
        {
            const Vec3 original_delta = toVec(point) - toVec(mesh.vertices[vertex]);
            if (norm(original_delta) > maximum_distance * 0.999) return false;
            replacements[vertex] = point;
        }
        std::vector<std::uint32_t> affected_faces;
        for (const auto& [vertex, point] : moves)
        {
            (void)point;
            for (const std::uint32_t face : vertex_faces[vertex])
                if (visit_stamp[face] != stamp)
                {
                    visit_stamp[face] = stamp;
                    affected_faces.push_back(face);
                }
        }
        ++stamp;
        if (stamp == 0)
        {
            std::fill(visit_stamp.begin(), visit_stamp.end(), 0);
            stamp = 1;
        }
        for (const std::uint32_t face : affected_faces)
        {
            const TriangleIndices triangle = current.triangles[face];
            std::array<Position3, 3> points{
                current.vertices[triangle[0]],
                current.vertices[triangle[1]],
                current.vertices[triangle[2]]};
            for (int corner = 0; corner < 3; ++corner)
            {
                const auto replacement = replacements.find(triangle[corner]);
                if (replacement != replacements.end())
                    points[corner] = replacement->second;
            }
            const Vec3 raw = cross(
                toVec(points[1]) - toVec(points[0]),
                toVec(points[2]) - toVec(points[0]));
            if (norm(raw) <= diagonal * diagonal * 1.0e-18) return false;
            if (dot(raw, original_normals[face]) <= 0.0) return false;
        }
        return true;
    };

    for (std::uint32_t round = 0; round < maximum_rounds; ++round)
    {
        const std::vector<LocalHalfspacePlane> current_planes = buildFacePlanes(result);
        const auto [is_concave, concave_vertex_count] =
            concaveVertices(result, current_planes);
        if (concave_vertex_count == 0) break;
        std::vector<std::uint8_t> visited_vertex(mesh.vertices.size(), 0);
        std::size_t round_moved = 0;
        for (std::uint32_t seed = 0; seed < mesh.vertices.size(); ++seed)
        {
            if (!is_concave[seed] || visited_vertex[seed]) continue;
            std::vector<std::uint32_t> component{seed};
            visited_vertex[seed] = 1;
            for (std::size_t cursor = 0; cursor < component.size(); ++cursor)
            {
                for (const std::uint32_t neighbor : vertex_neighbors[component[cursor]])
                    if (neighbor < is_concave.size() && is_concave[neighbor] &&
                        !visited_vertex[neighbor])
                    {
                        visited_vertex[neighbor] = 1;
                        component.push_back(neighbor);
                    }
            }
            if (component.size() < 2) continue;

            std::vector<std::uint32_t> face_queue;
            std::vector<std::uint32_t> depth;
            for (const std::uint32_t vertex : component)
                for (const std::uint32_t face : vertex_faces[vertex])
                    if (visit_stamp[face] != stamp)
                    {
                        visit_stamp[face] = stamp;
                        face_queue.push_back(face);
                        depth.push_back(0);
                    }
            for (std::size_t cursor = 0; cursor < face_queue.size(); ++cursor)
            {
                if (depth[cursor] >= neighborhood_depth) continue;
                for (const std::uint32_t neighbor : face_neighbors[face_queue[cursor]])
                    if (visit_stamp[neighbor] != stamp)
                    {
                        visit_stamp[neighbor] = stamp;
                        face_queue.push_back(neighbor);
                        depth.push_back(depth[cursor] + 1);
                    }
            }
            ++stamp;
            if (stamp == 0)
            {
                std::fill(visit_stamp.begin(), visit_stamp.end(), 0);
                stamp = 1;
            }

            std::vector<std::pair<std::uint32_t, Position3>> moves;
            for (const std::uint32_t vertex : component)
            {
                Vec3 candidate = toVec(result.vertices[vertex]);
                const Vec3 original = toVec(mesh.vertices[vertex]);
                for (int projection = 0; projection < 16; ++projection)
                {
                    double worst_signed_distance = 0.0;
                    std::uint32_t worst_face = invalid_surface_index;
                    for (const std::uint32_t face : face_queue)
                    {
                        const double signed_distance =
                            dot(current_planes[face].normal, candidate) -
                            current_planes[face].offset;
                        if (signed_distance < worst_signed_distance)
                        {
                            worst_signed_distance = signed_distance;
                            worst_face = face;
                        }
                    }
                    if (worst_face == invalid_surface_index ||
                        worst_signed_distance >= -concavity_tolerance)
                        break;
                    candidate = candidate +
                        current_planes[worst_face].normal * (-worst_signed_distance);
                    const Vec3 from_original = candidate - original;
                    const double moved_from_original = norm(from_original);
                    if (moved_from_original > maximum_distance * 0.999)
                    {
                        candidate = original +
                            from_original * ((maximum_distance * 0.999) /
                                             moved_from_original);
                        break;
                    }
                }
                if (norm(candidate - toVec(result.vertices[vertex])) <=
                    std::max(diagonal * 1.0e-10, 1.0e-6))
                    continue;
                const Position3 candidate_position = toPosition(candidate);
                moves.push_back({vertex, candidate_position});
            }
            if (!acceptsComponentMove(result, moves))
            {
                std::vector<std::pair<std::uint32_t, Position3>> scaled_moves;
                bool accepted_scaled = false;
                for (const double scale_factor : {0.5, 0.25, 0.125, 0.0625})
                {
                    scaled_moves.clear();
                    scaled_moves.reserve(moves.size());
                    for (const auto& [vertex, point] : moves)
                    {
                        const Vec3 current = toVec(result.vertices[vertex]);
                        const Vec3 target = toVec(point);
                        scaled_moves.push_back({
                            vertex, toPosition(current + (target - current) * scale_factor)});
                    }
                    if (acceptsComponentMove(result, scaled_moves))
                    {
                        moves = scaled_moves;
                        accepted_scaled = true;
                        break;
                    }
                }
                if (!accepted_scaled) continue;
            }
            for (const auto& [vertex, point] : moves)
            {
                result.vertices[vertex] = point;
                ++round_moved;
            }
        }
        if (round_moved == 0) break;
        total_moved += round_moved;
        const std::size_t score = estimateCoplanarMergeScore(result, relative_tolerance);
        std::cerr << "monitor: stage=phase2_outward_round round=" << round
                  << " concave_vertices=" << concave_vertex_count
                  << " moved_vertices=" << round_moved
                  << " scored_coplanar_triangles=" << score << '\n';
        const std::size_t round_index = round_vertices.size();
        round_vertices.push_back(result.vertices);
        round_cumulative_moved.push_back(total_moved);
        round_coplanar_scores.push_back(score);
        if (score < best_score)
        {
            best_score = score;
            best = result;
            best_moved = total_moved;
            selected_round = round_index;
        }
    }
    if (best_score == std::numeric_limits<std::size_t>::max())
        return {mesh, 0, 0, 0};
    OutwardEnvelopeResult output{std::move(best), best_moved, 0, best_score};
    output.selected_round = selected_round;
    output.round_vertices = std::move(round_vertices);
    output.round_cumulative_moved = std::move(round_cumulative_moved);
    output.round_coplanar_scores = std::move(round_coplanar_scores);
    return output;
}

struct ViewerLinks
{
    std::string source = "phase1.obj";
    std::string phase0;
};

std::optional<std::string> jsonStringValue(const std::string& text, const std::string& key)
{
    const std::string needle = "\"" + key + "\"";
    std::size_t position = text.find(needle);
    if (position == std::string::npos) return std::nullopt;
    position = text.find(':', position + needle.size());
    if (position == std::string::npos) return std::nullopt;
    position = text.find('"', position + 1);
    if (position == std::string::npos) return std::nullopt;
    std::string value;
    for (++position; position < text.size(); ++position)
    {
        const char character = text[position];
        if (character == '\\' && position + 1 < text.size())
        {
            value.push_back(text[++position]);
            continue;
        }
        if (character == '"') return value;
        value.push_back(character);
    }
    return std::nullopt;
}

ViewerLinks readPhase1ViewerLinks(const std::filesystem::path& phase1_halfedge)
{
    ViewerLinks links;
    std::ifstream stream(phase1_halfedge.parent_path() / "model.json");
    if (!stream) return links;
    std::ostringstream buffer;
    buffer << stream.rdbuf();
    const std::string text = buffer.str();
    if (auto source = jsonStringValue(text, "source")) links.source = *source;
    if (auto phase0 = jsonStringValue(text, "phase0_combined")) links.phase0 = *phase0;
    else if (auto caps = jsonStringValue(text, "phase0_caps")) links.phase0 = *caps;
    return links;
}

void writeObj(const std::filesystem::path& path, const MeshModel& mesh)
{
    std::ofstream stream(path);
    if (!stream) throw std::runtime_error("failed to create phase-2 OBJ: " + path.string());
    stream << std::setprecision(17) << "o " << mesh.name << '\n';
    for (const Position3 point : mesh.vertices)
        stream << "v " << point.x << ' ' << point.y << ' ' << point.z << '\n';
    for (const TriangleIndices face : mesh.triangles)
        stream << "f " << face[0] + 1 << ' ' << face[1] + 1 << ' ' << face[2] + 1 << '\n';
}

void writeErrorVisualization(
    const std::filesystem::path& path, const DirectedHausdorffCertificate& certificate)
{
    std::ofstream stream(path);
    if (!stream) throw std::runtime_error("failed to create phase-2 error visualization");
    stream << std::setprecision(17)
           << "{\n  \"maximum_pair\":{\"distance\":" << certificate.lower_bound
           << ",\"proxy\":[" << certificate.maximum_proxy_point[0] << ','
           << certificate.maximum_proxy_point[1] << ','
           << certificate.maximum_proxy_point[2] << "],\"source\":["
           << certificate.nearest_reference_point[0] << ','
           << certificate.nearest_reference_point[1] << ','
           << certificate.nearest_reference_point[2] << "]}\n}\n";
}

void writeMetadata(
    const std::filesystem::path& path, const Phase2HalfedgeStats& stats,
    const Phase2HalfedgeOptions& options, const ViewerLinks& links)
{
    std::ofstream stream(path);
    if (!stream) throw std::runtime_error("failed to create phase-2 metadata");
    stream << std::setprecision(17)
           << "{\n  \"stats\":{\"source_triangles\":" << stats.input_triangles
           << ",\"proxy_triangles\":" << stats.output_triangles
           << ",\"primitive_count\":0,\"primitive_types\":{}"
           << ",\"timings_seconds\":{\"total\":" << stats.elapsed_seconds << '}'
           << ",\"simplification_error\":{\"direction\":\"proxy_to_phase1\""
           << ",\"method\":\"triangle_bvh_1_lipschitz_adaptive_upper_bound\""
           << ",\"maximum_is_certified_upper_bound\":true"
           << ",\"maximum_distance_limit\":" << options.maximum_directed_hausdorff
           << ",\"maximum_distance\":" << stats.directed_hausdorff.upper_bound
           << ",\"observed_lower_bound\":" << stats.directed_hausdorff.lower_bound
           << ",\"reference_queries\":" << stats.directed_hausdorff.reference_queries
           << ",\"subdivision_nodes\":" << stats.directed_hausdorff.subdivision_nodes
           << "},\"phase2\":{\"selected_candidate\":\""
           << stats.selected_candidate << "\",\"outward_vertices_moved\":"
           << stats.outward_vertices_moved << ",\"outward_support_planes\":"
           << stats.outward_support_planes << "}},\n"
           << "  \"source\":\"phase1.obj\",\n"
           << "  \"comparison_source_label\":\"阶段 1\",\n"
           << "  \"phase1_halfedge\":\"phase1_halfedge.bin\",\n"
           << "  \"phase2_halfedge\":\"phase2_halfedge.bin\",\n"
           << "  \"phase4_triangulated\":\"phase2.obj\",\n"
           << "  \"proxy\":\"phase2.obj\",\n"
           << "  \"open_error_visualization\":\"phase2_error.json\",\n"
           << "  \"proxy_components\":[],\n"
           << "  \"viewer_stages\":[\"source\",\"phase1\",\"phase2\",\"error\",\"split\"]\n}\n";
}

void writeViewerMetadata(
    const std::filesystem::path& path, const Phase2HalfedgeStats& stats,
    const Phase2HalfedgeOptions& options, const ViewerLinks& links,
    const bool has_phase2_halfedge, const OutwardEnvelopeResult* rounds)
{
    std::ofstream stream(path);
    if (!stream) throw std::runtime_error("failed to create phase-2 metadata");
    stream << std::setprecision(17)
           << "{\n  \"stats\":{\"source_triangles\":" << stats.input_triangles
           << ",\"proxy_triangles\":" << stats.output_triangles
           << ",\"primitive_count\":0,\"primitive_types\":{}"
           << ",\"timings_seconds\":{\"total\":" << stats.elapsed_seconds << '}'
           << ",\"simplification_error\":{\"direction\":\"proxy_to_phase1\""
           << ",\"method\":\"triangle_bvh_1_lipschitz_adaptive_upper_bound\""
           << ",\"maximum_is_certified_upper_bound\":true"
           << ",\"maximum_distance_limit\":" << options.maximum_directed_hausdorff
           << ",\"maximum_distance\":" << stats.directed_hausdorff.upper_bound
           << ",\"observed_lower_bound\":" << stats.directed_hausdorff.lower_bound
           << ",\"reference_queries\":" << stats.directed_hausdorff.reference_queries
           << ",\"subdivision_nodes\":" << stats.directed_hausdorff.subdivision_nodes
           << "},\"phase2\":{\"selected_candidate\":\""
           << stats.selected_candidate << "\",\"outward_vertices_moved\":"
           << stats.outward_vertices_moved << ",\"outward_support_planes\":"
           << stats.outward_support_planes;
    if (rounds && !rounds->round_vertices.empty())
    {
        stream << ",\"selected_round\":" << rounds->selected_round
               << ",\"round_count\":" << rounds->round_vertices.size()
               << ",\"rounds\":[";
        for (std::size_t index = 0; index < rounds->round_vertices.size(); ++index)
        {
            if (index) stream << ',';
            stream << "{\"round\":" << index
                   << ",\"halfedge\":\"phase2_round_" << index << "_halfedge.bin\""
                   << ",\"cumulative_moved_vertices\":"
                   << rounds->round_cumulative_moved[index]
                   << ",\"scored_coplanar_triangles\":"
                   << rounds->round_coplanar_scores[index]
                   << ",\"selected\":"
                   << (index == rounds->selected_round ? "true" : "false") << '}';
        }
        stream << ']';
    }
    stream << "}},\n"
           << "  \"source\":\"" << links.source << "\",\n"
           << "  \"comparison_source_label\":\"original\",\n";
    if (!links.phase0.empty())
        stream << "  \"phase0_combined\":\"" << links.phase0 << "\",\n";
    stream << "  \"phase1_halfedge\":\"phase1_halfedge.bin\",\n";
    if (has_phase2_halfedge)
        stream << "  \"phase2_halfedge\":\"phase2_halfedge.bin\",\n";
    if (rounds && !rounds->round_vertices.empty())
    {
        stream << "  \"phase2_round_halfedges\":[";
        for (std::size_t index = 0; index < rounds->round_vertices.size(); ++index)
        {
            if (index) stream << ',';
            stream << "\"phase2_round_" << index << "_halfedge.bin\"";
        }
        stream << "],\n";
    }
    stream << "  \"phase4_triangulated\":\"phase2.obj\",\n"
           << "  \"proxy\":\"phase2.obj\",\n"
           << "  \"open_error_visualization\":\"phase2_error.json\",\n"
           << "  \"proxy_components\":[],\n"
           << "  \"viewer_stages\":[\"source\"";
    if (!links.phase0.empty()) stream << ",\"phase0\"";
    stream << ",\"phase1\",\"phase2\",\"error\",\"split\"]\n}\n";
}

} // namespace

Phase2HalfedgeStats generatePhase2Halfedge(
    const std::filesystem::path& phase1_halfedge,
    const std::filesystem::path& output_directory,
    const Phase2HalfedgeOptions& options)
{
    const auto started = std::chrono::steady_clock::now();
    if (!std::isfinite(options.maximum_directed_hausdorff) ||
        options.maximum_directed_hausdorff < 0.0)
        throw std::invalid_argument("maximum directed Hausdorff distance must be finite and non-negative");
    const OrientedSurfaceMesh phase1 = readAnalysisHalfedgeMesh(phase1_halfedge);
    const ViewerLinks viewer_links = readPhase1ViewerLinks(phase1_halfedge);
    (void)validatePhase1HalfedgeTopology(phase1);
    const MeshModel& reference = phase1.geometry;

    Phase2HalfedgeStats stats;
    stats.input_vertices = reference.vertices.size();
    stats.input_triangles = reference.triangles.size();
    stats.output_vertices = stats.input_vertices;
    stats.output_triangles = stats.input_triangles;
    stats.selected_candidate = "phase1_exact_fallback";
    stats.directed_hausdorff.passed = true;
    stats.source_coverage_certified = true;
    MeshModel selected = reference;
    selected.name = "phase2_exact_fallback";
    OrientedSurfaceMesh selected_halfedge = phase1;
    OutwardEnvelopeResult selected_outward_rounds;
    bool has_selected_outward_rounds = false;

    const auto considerOutwardEnvelope = [&](OutwardEnvelopeResult candidate,
                                             const std::string& name)
    {
        if (candidate.geometry.triangles.empty() || candidate.moved_vertices == 0)
            return;
        OrientedSurfaceMesh candidate_halfedge = phase1;
        candidate_halfedge.geometry = candidate.geometry;
        candidate_halfedge.geometry.name = name;
        (void)validateWorkingHalfedgeTopology(candidate_halfedge);
        const DirectedHausdorffCertificate certificate = certifyDirectedHausdorff(
            candidate.geometry, reference, options.maximum_directed_hausdorff,
            options.maximum_certificate_depth, options.numerical_tolerance);
        if (!certificate.passed || certificate.upper_bound >
            options.maximum_directed_hausdorff + options.numerical_tolerance)
            return;
        selected = std::move(candidate.geometry);
        selected.name = name;
        stats.selected_candidate = name;
        stats.directed_hausdorff = certificate;
        stats.output_vertices = selected.vertices.size();
        stats.output_triangles = selected.triangles.size();
        stats.outward_vertices_moved = candidate.moved_vertices;
        stats.outward_support_planes = candidate.support_planes;
        ++stats.candidate_faces_accepted;
        selected_halfedge = std::move(candidate_halfedge);
        selected_outward_rounds = std::move(candidate);
        has_selected_outward_rounds = true;
        std::cerr << std::setprecision(17)
                  << "monitor: stage=phase2_candidate_accepted candidate=" << name
                  << " triangles=" << selected.triangles.size()
                  << " moved_vertices=" << stats.outward_vertices_moved
                  << " support_planes=" << stats.outward_support_planes
                  << " lower_bound=" << certificate.lower_bound
                  << " upper_bound=" << certificate.upper_bound << '\n';
    };

    OutwardEnvelopeResult outward = makeConservativeOutwardEnvelope(
        phase1, options.maximum_directed_hausdorff,
        options.coplanar_relative_tolerance);
    stats.candidate_faces_considered = reference.triangles.size();
    stats.outward_support_planes = outward.support_planes;
    considerOutwardEnvelope(
        std::move(outward), "conservative_outward_envelope");

    stats.output_halfedges = selected_halfedge.halfedges.size();
    const WorkingHalfedgeValidationReport phase2_validation =
        validateWorkingHalfedgeTopology(selected_halfedge);
    stats.paired_edges = phase2_validation.paired_edges;
    stats.boundary_halfedges = phase2_validation.boundary_halfedges;
    stats.face_components = phase2_validation.face_components;
    std::filesystem::create_directories(output_directory);
    writeAnalysisHalfedgeMesh(output_directory / "phase2_halfedge.bin", selected_halfedge);
    if (has_selected_outward_rounds)
    {
        for (std::size_t index = 0; index < selected_outward_rounds.round_vertices.size(); ++index)
        {
            OrientedSurfaceMesh round_halfedge = phase1;
            round_halfedge.geometry.name =
                "phase2_conservative_outward_envelope_round_" + std::to_string(index);
            round_halfedge.geometry.vertices = selected_outward_rounds.round_vertices[index];
            (void)validateWorkingHalfedgeTopology(round_halfedge);
            writeAnalysisHalfedgeMesh(
                output_directory /
                    ("phase2_round_" + std::to_string(index) + "_halfedge.bin"),
                round_halfedge);
        }
    }
    writeAnalysisHalfedgeMesh(output_directory / "phase1_halfedge.bin", phase1);
    writeObj(output_directory / "phase1.obj", reference);
    writeObj(output_directory / "phase2.obj", selected);
    writeErrorVisualization(output_directory / "phase2_error.json", stats.directed_hausdorff);
    stats.elapsed_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    writeViewerMetadata(
        output_directory / "model.json", stats, options, viewer_links,
        true, has_selected_outward_rounds ? &selected_outward_rounds : nullptr);
    std::ofstream manifest(output_directory / "viewer_manifest.json");
    if (!manifest) throw std::runtime_error("failed to create phase-2 viewer manifest");
    manifest << "{\n  \"algorithm\":\"CppPhase2ConservativeOutwardEnvelopeV1\",\n"
             << "  \"complete\":true,\n  \"model_count\":1,\n  \"models\":[{\"id\":\""
             << (options.model_id.empty() ? "0" : options.model_id)
             << "\",\"metadata\":\"model.json\"}],\n"
             << "  \"options\":{\"maximum_directed_hausdorff\":"
             << std::setprecision(17) << options.maximum_directed_hausdorff
             << ",\"direction\":\"proxy_to_phase1\"}\n}\n";
    return stats;
}

} // namespace pqss_proxy_mesh
