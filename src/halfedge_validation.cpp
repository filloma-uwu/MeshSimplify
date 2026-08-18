#include "pqss_proxy_mesh/halfedge_validation.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>

namespace pqss_proxy_mesh
{
namespace
{

using EdgeKey = std::array<std::uint32_t, 2>;

EdgeKey edgeKey(const std::uint32_t first, const std::uint32_t second)
{
    return {std::min(first, second), std::max(first, second)};
}

double signedVolume(const MeshModel& geometry)
{
    double volume = 0.0;
    for (const TriangleIndices triangle : geometry.triangles)
    {
        const Position3 a = geometry.vertices[triangle[0]];
        const Position3 b = geometry.vertices[triangle[1]];
        const Position3 c = geometry.vertices[triangle[2]];
        volume += a.x * (b.y * c.z - b.z * c.y) +
                  a.y * (b.z * c.x - b.x * c.z) +
                  a.z * (b.x * c.y - b.y * c.x);
    }
    return volume / 6.0;
}

} // namespace

HalfedgeTopologyValidationReport validateClosedHalfedgeTopology(
    const OrientedSurfaceMesh& mesh)
{
    const std::size_t vertex_count = mesh.geometry.vertices.size();
    const std::size_t face_count = mesh.geometry.triangles.size();
    const std::size_t halfedge_count = mesh.halfedges.size();
    if (vertex_count == 0 || face_count == 0)
        throw std::runtime_error("phase-1 halfedge mesh is empty");
    if (halfedge_count != face_count * 3 ||
        mesh.face_halfedges.size() != face_count ||
        mesh.vertex_halfedges.size() != vertex_count)
        throw std::runtime_error("phase-1 halfedge arrays have inconsistent counts");

    std::map<EdgeKey, std::vector<std::uint32_t>> edge_uses;
    std::set<TriangleIndices> unique_faces;
    std::vector<std::uint32_t> vertex_valence(vertex_count, 0);
    for (std::uint32_t face = 0; face < face_count; ++face)
    {
        const TriangleIndices triangle = mesh.geometry.triangles[face];
        for (const std::uint32_t vertex : triangle)
            if (vertex >= vertex_count)
                throw std::runtime_error("phase-1 triangle has an invalid vertex");
        if (triangle[0] == triangle[1] || triangle[1] == triangle[2] ||
            triangle[2] == triangle[0])
            throw std::runtime_error("phase-1 mesh has a degenerate triangle index cycle");
        for (const std::uint32_t vertex : triangle)
        {
            const Position3 point = mesh.geometry.vertices[vertex];
            if (!std::isfinite(point.x) || !std::isfinite(point.y) ||
                !std::isfinite(point.z))
                throw std::runtime_error("phase-1 mesh has a non-finite vertex");
        }
        TriangleIndices signature = triangle;
        std::sort(signature.begin(), signature.end());
        if (!unique_faces.insert(signature).second)
            throw std::runtime_error("phase-1 mesh has a duplicate triangle");

        const std::uint32_t first = mesh.face_halfedges[face];
        if (first >= halfedge_count)
            throw std::runtime_error("phase-1 face has an invalid halfedge");
        std::uint32_t edge = first;
        for (std::uint32_t local = 0; local < 3; ++local)
        {
            if (edge >= halfedge_count)
                throw std::runtime_error("phase-1 face cycle leaves the halfedge array");
            const SurfaceHalfedge& halfedge = mesh.halfedges[edge];
            if (halfedge.face != face || halfedge.origin != triangle[local] ||
                halfedge.next >= halfedge_count)
                throw std::runtime_error("phase-1 face cycle is inconsistent");
            ++vertex_valence[halfedge.origin];
            edge_uses[edgeKey(halfedge.origin, triangle[(local + 1) % 3])]
                .push_back(edge);
            edge = halfedge.next;
        }
        if (edge != first)
            throw std::runtime_error("phase-1 face is not a three-halfedge cycle");
    }

    for (const auto& [key, uses] : edge_uses)
    {
        (void)key;
        if (uses.size() != 2)
            throw std::runtime_error("phase-1 mesh is not edge-manifold and closed");
        const std::uint32_t first = uses[0];
        const std::uint32_t second = uses[1];
        if (mesh.halfedges[first].opposite != second ||
            mesh.halfedges[second].opposite != first)
            throw std::runtime_error("phase-1 opposite relation is missing or asymmetric");
        const std::uint32_t first_end =
            mesh.halfedges[mesh.halfedges[first].next].origin;
        const std::uint32_t second_end =
            mesh.halfedges[mesh.halfedges[second].next].origin;
        if (mesh.halfedges[first].origin != second_end ||
            mesh.halfedges[second].origin != first_end)
            throw std::runtime_error("phase-1 adjacent faces have inconsistent orientation");
    }

    for (std::uint32_t vertex = 0; vertex < vertex_count; ++vertex)
    {
        const std::uint32_t start = mesh.vertex_halfedges[vertex];
        if (start >= halfedge_count || mesh.halfedges[start].origin != vertex)
            throw std::runtime_error("phase-1 vertex has no valid outgoing halfedge");
        std::uint32_t current = start;
        std::uint32_t fan_size = 0;
        do
        {
            if (current >= halfedge_count ||
                mesh.halfedges[current].origin != vertex ||
                ++fan_size > vertex_valence[vertex])
                throw std::runtime_error("phase-1 vertex fan is invalid");
            const std::uint32_t next = mesh.halfedges[current].next;
            const std::uint32_t previous = mesh.halfedges[next].next;
            current = mesh.halfedges[previous].opposite;
        }
        while (current != start);
        if (fan_size != vertex_valence[vertex])
            throw std::runtime_error("phase-1 vertex has multiple disconnected fans");
    }

    std::vector<std::uint8_t> visited(face_count, 0);
    std::size_t component_count = 0;
    for (std::uint32_t seed = 0; seed < face_count; ++seed)
    {
        if (visited[seed]) continue;
        ++component_count;
        std::vector<std::uint32_t> queue{seed};
        visited[seed] = 1;
        for (std::size_t cursor = 0; cursor < queue.size(); ++cursor)
        {
            std::uint32_t edge = mesh.face_halfedges[queue[cursor]];
            for (int local = 0; local < 3; ++local)
            {
                const std::uint32_t opposite = mesh.halfedges[edge].opposite;
                if (opposite >= halfedge_count)
                    throw std::runtime_error("phase-1 mesh has a boundary halfedge");
                const std::uint32_t neighbor = mesh.halfedges[opposite].face;
                if (neighbor >= face_count)
                    throw std::runtime_error("phase-1 opposite has an invalid face");
                if (!visited[neighbor])
                {
                    visited[neighbor] = 1;
                    queue.push_back(neighbor);
                }
                edge = mesh.halfedges[edge].next;
            }
        }
    }
    const int euler = static_cast<int>(vertex_count) -
        static_cast<int>(edge_uses.size()) + static_cast<int>(face_count);
    if (euler != static_cast<int>(2 * component_count))
        throw std::runtime_error("halfedge components are not closed genus-zero surfaces");
    const double volume = signedVolume(mesh.geometry);
    if (!std::isfinite(volume))
        throw std::runtime_error("phase-1 halfedge mesh has non-finite oriented volume");

    HalfedgeTopologyValidationReport report;
    report.vertices = vertex_count;
    report.triangles = face_count;
    report.halfedges = halfedge_count;
    report.boundary_halfedges = 0;
    report.face_components = component_count;
    report.euler_characteristic = euler;
    report.signed_volume = volume;
    return report;
}

HalfedgeTopologyValidationReport validatePhase1HalfedgeTopology(
    const OrientedSurfaceMesh& mesh)
{
    HalfedgeTopologyValidationReport report = validateClosedHalfedgeTopology(mesh);
    if (report.face_components != 1)
        throw std::runtime_error("phase-1 halfedge mesh has multiple face components");
    return report;
}

} // namespace pqss_proxy_mesh
