#include "pqss_proxy_mesh/hausdorff_simplifier.hpp"
#include "pqss_proxy_mesh/halfedge_validation.hpp"

#include "QuickHull.hpp"
#include "clipper2/clipper.h"
#include "clipper2/clipper.triangulation.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <unordered_map>

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#include <psapi.h>
#pragma comment(lib, "psapi.lib")
#endif

namespace pqss_proxy_mesh
{
namespace
{

double workingSetMiB()
{
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS counters{};
    counters.cb = sizeof(counters);
    if (GetProcessMemoryInfo(GetCurrentProcess(), &counters, sizeof(counters)))
        return static_cast<double>(counters.WorkingSetSize) / (1024.0 * 1024.0);
#endif
    return -1.0;
}

void logMonitor(const char* stage, const std::string& candidate = {})
{
    std::cerr << "monitor: stage=" << stage;
    if (!candidate.empty()) std::cerr << " candidate=" << candidate;
    const double working_set = workingSetMiB();
    if (working_set >= 0.0)
        std::cerr << " working_set_mib=" << working_set;
    std::cerr << '\n';
}

struct Vec3
{
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

Vec3 toVec(const Position3 p) { return {p.x, p.y, p.z}; }
Position3 toPosition(const Vec3 p) { return {p.x, p.y, p.z}; }
Vec3 operator+(const Vec3 a, const Vec3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
Vec3 operator-(const Vec3 a, const Vec3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
Vec3 operator*(const Vec3 a, const double s) { return {a.x * s, a.y * s, a.z * s}; }
double dot(const Vec3 a, const Vec3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
Vec3 cross(const Vec3 a, const Vec3 b)
{
    return {a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x};
}
double normSquared(const Vec3 a) { return dot(a, a); }
double norm(const Vec3 a) { return std::sqrt(normSquared(a)); }

struct Bounds
{
    Vec3 lower{std::numeric_limits<double>::infinity(),
               std::numeric_limits<double>::infinity(),
               std::numeric_limits<double>::infinity()};
    Vec3 upper{-std::numeric_limits<double>::infinity(),
               -std::numeric_limits<double>::infinity(),
               -std::numeric_limits<double>::infinity()};
};

void include(Bounds& bounds, const Vec3 p)
{
    bounds.lower.x = std::min(bounds.lower.x, p.x);
    bounds.lower.y = std::min(bounds.lower.y, p.y);
    bounds.lower.z = std::min(bounds.lower.z, p.z);
    bounds.upper.x = std::max(bounds.upper.x, p.x);
    bounds.upper.y = std::max(bounds.upper.y, p.y);
    bounds.upper.z = std::max(bounds.upper.z, p.z);
}

double boundsDistanceSquared(const Bounds& bounds, const Vec3 p)
{
    const auto axis = [](const double value, const double lower, const double upper)
    {
        if (value < lower) return lower - value;
        if (value > upper) return value - upper;
        return 0.0;
    };
    const double dx = axis(p.x, bounds.lower.x, bounds.upper.x);
    const double dy = axis(p.y, bounds.lower.y, bounds.upper.y);
    const double dz = axis(p.z, bounds.lower.z, bounds.upper.z);
    return dx*dx + dy*dy + dz*dz;
}

Vec3 closestPointOnTriangle(const Vec3 p, const Vec3 a, const Vec3 b, const Vec3 c)
{
    const Vec3 ab = b - a;
    const Vec3 ac = c - a;
    const Vec3 ap = p - a;
    const double d1 = dot(ab, ap);
    const double d2 = dot(ac, ap);
    if (d1 <= 0.0 && d2 <= 0.0) return a;

    const Vec3 bp = p - b;
    const double d3 = dot(ab, bp);
    const double d4 = dot(ac, bp);
    if (d3 >= 0.0 && d4 <= d3) return b;

    const double vc = d1*d4 - d3*d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
        return a + ab * (d1 / (d1 - d3));

    const Vec3 cp = p - c;
    const double d5 = dot(ab, cp);
    const double d6 = dot(ac, cp);
    if (d6 >= 0.0 && d5 <= d6) return c;

    const double vb = d5*d2 - d1*d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
        return a + ac * (d2 / (d2 - d6));

    const double va = d3*d6 - d5*d4;
    if (va <= 0.0 && d4 - d3 >= 0.0 && d5 - d6 >= 0.0)
        return b + (c - b) * ((d4 - d3) / ((d4 - d3) + (d5 - d6)));

    const double denominator = 1.0 / (va + vb + vc);
    return a + ab * (vb * denominator) + ac * (vc * denominator);
}

class TriangleBvh
{
public:
    explicit TriangleBvh(const MeshModel& mesh) : mesh_(mesh)
    {
        if (mesh.triangles.empty()) throw std::invalid_argument("reference mesh is empty");
        ids_.resize(mesh.triangles.size());
        std::iota(ids_.begin(), ids_.end(), std::uint32_t{0});
        // Reserve the exact node count produced by the median split. A loose
        // per-triangle bound needlessly consumes hundreds of MiB on large
        // phase-1 meshes.
        nodes_.reserve(requiredNodeCount(mesh.triangles.size()));
        build(0, static_cast<std::uint32_t>(ids_.size()));
    }

    std::pair<double, Vec3> nearest(const Vec3 point) const
    {
        struct Visit { double lower; std::uint32_t node; };
        struct Greater { bool operator()(const Visit& a, const Visit& b) const { return a.lower > b.lower; } };
        std::priority_queue<Visit, std::vector<Visit>, Greater> queue;
        queue.push({boundsDistanceSquared(nodes_[0].bounds, point), 0});
        double best = std::numeric_limits<double>::infinity();
        Vec3 nearest{};
        while (!queue.empty())
        {
            const Visit visit = queue.top();
            queue.pop();
            if (visit.lower >= best) continue;
            const Node& node = nodes_[visit.node];
            if (node.count != 0)
            {
                for (std::uint32_t offset = 0; offset < node.count; ++offset)
                {
                    const auto face = mesh_.triangles[ids_[node.begin + offset]];
                    const Vec3 candidate = closestPointOnTriangle(
                        point, toVec(mesh_.vertices[face[0]]),
                        toVec(mesh_.vertices[face[1]]), toVec(mesh_.vertices[face[2]]));
                    const double distance = normSquared(point - candidate);
                    if (distance < best) { best = distance; nearest = candidate; }
                }
            }
            else
            {
                for (const std::uint32_t child : {node.left, node.right})
                {
                    const double lower = boundsDistanceSquared(nodes_[child].bounds, point);
                    if (lower < best) queue.push({lower, child});
                }
            }
        }
        return {std::sqrt(best), nearest};
    }

private:
    static std::size_t requiredNodeCount(const std::size_t triangle_count)
    {
        if (triangle_count <= 8) return 1;
        const std::size_t left = triangle_count / 2;
        return 1 + requiredNodeCount(left) +
               requiredNodeCount(triangle_count - left);
    }

    struct Node
    {
        Bounds bounds;
        std::uint32_t begin = 0;
        std::uint32_t count = 0;
        std::uint32_t left = 0;
        std::uint32_t right = 0;
    };

    std::uint32_t build(const std::uint32_t begin, const std::uint32_t end)
    {
        const std::uint32_t id = static_cast<std::uint32_t>(nodes_.size());
        nodes_.push_back({});
        Bounds bounds;
        Bounds centers;
        for (std::uint32_t i = begin; i < end; ++i)
        {
            const auto face = mesh_.triangles[ids_[i]];
            Vec3 center{};
            for (const auto vertex : face)
            {
                const Vec3 p = toVec(mesh_.vertices[vertex]);
                include(bounds, p);
                center = center + p * (1.0 / 3.0);
            }
            include(centers, center);
        }
        nodes_[id].bounds = bounds;
        const std::uint32_t count = end - begin;
        if (count <= 8)
        {
            nodes_[id].begin = begin;
            nodes_[id].count = count;
            return id;
        }
        const Vec3 extent = centers.upper - centers.lower;
        const int axis = extent.y > extent.x ? (extent.z > extent.y ? 2 : 1)
                                             : (extent.z > extent.x ? 2 : 0);
        const auto coordinate = [axis](const Vec3 p) { return axis == 0 ? p.x : axis == 1 ? p.y : p.z; };
        const std::uint32_t middle = begin + count / 2;
        std::nth_element(ids_.begin() + begin, ids_.begin() + middle, ids_.begin() + end,
            [&](const std::uint32_t first, const std::uint32_t second)
            {
                const auto centroid = [&](const std::uint32_t face_id)
                {
                    const auto f = mesh_.triangles[face_id];
                    return (toVec(mesh_.vertices[f[0]]) + toVec(mesh_.vertices[f[1]]) +
                            toVec(mesh_.vertices[f[2]])) * (1.0 / 3.0);
                };
                return coordinate(centroid(first)) < coordinate(centroid(second));
            });
        const std::uint32_t left = build(begin, middle);
        const std::uint32_t right = build(middle, end);
        nodes_[id].left = left;
        nodes_[id].right = right;
        return id;
    }

    const MeshModel& mesh_;
    std::vector<std::uint32_t> ids_;
    std::vector<Node> nodes_;
};

DirectedHausdorffCertificate certifyWithBvh(
    const MeshModel& proxy, const TriangleBvh& reference_bvh,
    const double maximum_distance, const std::uint32_t maximum_depth,
    const double numerical_tolerance, const bool tight_triangle_radius = true,
    const bool collect_global_violation = false)
{
    DirectedHausdorffCertificate result;
    result.passed = true;
    struct Node { Vec3 a,b,c; std::uint32_t depth; };
    std::vector<Node> stack;
    const double tolerance = numerical_tolerance * std::max(1.0, maximum_distance);
    for (std::size_t proxy_face = 0; proxy_face < proxy.triangles.size(); ++proxy_face)
    {
        const TriangleIndices face = proxy.triangles[proxy_face];
        stack.push_back({toVec(proxy.vertices[face[0]]), toVec(proxy.vertices[face[1]]),
                         toVec(proxy.vertices[face[2]]), 0});
        while (!stack.empty())
        {
            const Node node = stack.back();
            stack.pop_back();
            ++result.subdivision_nodes;
            bool node_exceeds_limit = false;
            std::array<double, 4> distances{};
            std::array<Vec3, 4> nearest{};
            const std::array<Vec3, 4> samples{node.a, node.b, node.c,
                                              (node.a + node.b + node.c) * (1.0 / 3.0)};
            for (int i = 0; i < 4; ++i)
            {
                std::tie(distances[i], nearest[i]) = reference_bvh.nearest(samples[i]);
                ++result.reference_queries;
                if (distances[i] > result.lower_bound)
                {
                    result.lower_bound = distances[i];
                    result.maximum_proxy_point = {samples[i].x, samples[i].y, samples[i].z};
                    result.nearest_reference_point = {nearest[i].x, nearest[i].y, nearest[i].z};
                    result.maximum_proxy_face = proxy_face;
                }
                if (distances[i] > maximum_distance + tolerance)
                {
                    result.passed = false;
                    result.failure_reason = "witness point exceeds the directed Hausdorff limit";
                    result.upper_bound = std::numeric_limits<double>::infinity();
                    node_exceeds_limit = true;
                    if (!collect_global_violation) return result;
                }
            }
            if (node_exceeds_limit) continue;
            const Vec3 centroid = (node.a + node.b + node.c) * (1.0 / 3.0);
            const double radius = tight_triangle_radius
                ? std::max({norm(node.a - centroid), norm(node.b - centroid),
                            norm(node.c - centroid)})
                : std::max({norm(node.a-node.b), norm(node.b-node.c),
                            norm(node.c-node.a)}) * (2.0 / 3.0);
            const double local_upper = distances[3] + radius;
            if (local_upper <= maximum_distance + tolerance)
            {
                result.upper_bound = std::max(result.upper_bound, local_upper);
                continue;
            }
            if (node.depth >= maximum_depth)
            {
                result.passed = false;
                result.failure_reason = "certificate depth exhausted before proving the upper bound";
                return result;
            }
            const Vec3 ab = (node.a + node.b) * 0.5;
            const Vec3 bc = (node.b + node.c) * 0.5;
            const Vec3 ca = (node.c + node.a) * 0.5;
            const auto depth = node.depth + 1;
            stack.push_back({node.a,ab,ca,depth});
            stack.push_back({ab,node.b,bc,depth});
            stack.push_back({ca,bc,node.c,depth});
            stack.push_back({ab,bc,ca,depth});
        }
    }
    return result;
}

MeshModel makeBox(const MeshModel& reference)
{
    Bounds bounds;
    for (const Position3 p : reference.vertices) include(bounds, toVec(p));
    MeshModel box;
    box.name = "axis_aligned_enclosing_box";
    for (int x = 0; x < 2; ++x)
        for (int y = 0; y < 2; ++y)
            for (int z = 0; z < 2; ++z)
                box.vertices.push_back({x ? bounds.upper.x : bounds.lower.x,
                                        y ? bounds.upper.y : bounds.lower.y,
                                        z ? bounds.upper.z : bounds.lower.z});
    box.triangles = {
        {0,2,3},{0,3,1}, {4,5,7},{4,7,6},
        {0,1,5},{0,5,4}, {2,6,7},{2,7,3},
        {0,4,6},{0,6,2}, {1,3,7},{1,7,5}};
    return box;
}

void orientPositiveVolume(MeshModel& mesh)
{
    double six_volume = 0.0;
    for (const TriangleIndices face : mesh.triangles)
    {
        const Vec3 a = toVec(mesh.vertices[face[0]]);
        const Vec3 b = toVec(mesh.vertices[face[1]]);
        const Vec3 c = toVec(mesh.vertices[face[2]]);
        six_volume += dot(a, cross(b, c));
    }
    if (six_volume < 0.0)
        for (TriangleIndices& face : mesh.triangles) std::swap(face[1], face[2]);
}

void orientConvexFacesOutward(MeshModel& mesh)
{
    Vec3 center{};
    for (const Position3 point : mesh.vertices) center = center + toVec(point);
    center = center * (1.0 / mesh.vertices.size());
    for (TriangleIndices& face : mesh.triangles)
    {
        const Vec3 a = toVec(mesh.vertices[face[0]]);
        const Vec3 normal = cross(toVec(mesh.vertices[face[1]]) - a,
                                  toVec(mesh.vertices[face[2]]) - a);
        if (dot(normal, center - a) > 0.0) std::swap(face[1], face[2]);
    }
}

struct HalfspacePolytope
{
    std::vector<Vec3> normals;
    std::vector<double> offsets;
};

HalfspacePolytope extractHalfspaces(const MeshModel& mesh)
{
    HalfspacePolytope result;
    Vec3 center{};
    for (const Position3 point : mesh.vertices) center = center + toVec(point);
    center = center * (1.0 / mesh.vertices.size());
    for (const TriangleIndices face : mesh.triangles)
    {
        const Vec3 a = toVec(mesh.vertices[face[0]]);
        Vec3 normal = cross(toVec(mesh.vertices[face[1]]) - a,
                            toVec(mesh.vertices[face[2]]) - a);
        const double length = norm(normal);
        if (!(length > 0.0)) continue;
        normal = normal * (1.0 / length);
        if (dot(normal, center - a) > 0.0) normal = normal * -1.0;
        const double offset = dot(normal, a);
        bool duplicate = false;
        for (std::size_t index = 0; index < result.normals.size(); ++index)
            duplicate |= dot(normal, result.normals[index]) > 1.0 - 1.0e-10 &&
                         std::abs(offset - result.offsets[index]) <= 1.0e-8;
        if (!duplicate)
        {
            result.normals.push_back(normal);
            result.offsets.push_back(offset);
        }
    }
    return result;
}

MeshModel makeConvexHull(const MeshModel& reference);

MeshModel buildHalfspacePolytope(const HalfspacePolytope& planes)
{
    MeshModel result;
    result.name = "closed_halfspace_polytope";
    if (planes.normals.size() < 4) throw std::runtime_error("too few halfspaces");
    double scale = 1.0;
    for (const double offset : planes.offsets) scale = std::max(scale, std::abs(offset));
    const double tolerance = std::max(1.0e-9, scale * 1.0e-12);
    for (std::size_t first = 0; first < planes.normals.size(); ++first)
        for (std::size_t second = first + 1; second < planes.normals.size(); ++second)
            for (std::size_t third = second + 1; third < planes.normals.size(); ++third)
            {
                const Vec3 cross_second_third = cross(
                    planes.normals[second], planes.normals[third]);
                const double determinant = dot(planes.normals[first], cross_second_third);
                if (std::abs(determinant) <= 1.0e-12) continue;
                const Vec3 point = (cross_second_third * planes.offsets[first] +
                    cross(planes.normals[third], planes.normals[first]) * planes.offsets[second] +
                    cross(planes.normals[first], planes.normals[second]) * planes.offsets[third]) *
                    (1.0 / determinant);
                bool inside = true;
                for (std::size_t plane = 0; plane < planes.normals.size(); ++plane)
                    inside &= dot(planes.normals[plane], point) <=
                              planes.offsets[plane] + tolerance;
                if (!inside) continue;
                if (std::any_of(result.vertices.begin(), result.vertices.end(),
                    [&](const Position3 existing)
                    { return norm(toVec(existing) - point) <= tolerance; }))
                    continue;
                result.vertices.push_back(toPosition(point));
            }
    if (result.vertices.size() < 4)
        throw std::runtime_error("halfspace intersection is lower dimensional");
    return makeConvexHull(result);
}

MeshModel clipClosedConvexMesh(const MeshModel& input, const Vec3 normal,
                               const double offset, const double tolerance)
{
    MeshModel result;
    result.name = "support_clipped_closed_convex_mesh";
    using VertexKey = std::array<std::int64_t, 3>;
    std::map<VertexKey, std::uint32_t> vertex_map;
    const double quantum = std::max(tolerance, 1.0e-10);
    const auto appendVertex = [&](const Vec3 point)
    {
        const VertexKey key{
            static_cast<std::int64_t>(std::llround(point.x / quantum)),
            static_cast<std::int64_t>(std::llround(point.y / quantum)),
            static_cast<std::int64_t>(std::llround(point.z / quantum))};
        const auto [iterator, inserted] = vertex_map.try_emplace(
            key, static_cast<std::uint32_t>(result.vertices.size()));
        if (inserted) result.vertices.push_back(toPosition(point));
        return iterator->second;
    };
    std::vector<Vec3> cap_points;
    for (const TriangleIndices face : input.triangles)
    {
        std::vector<Vec3> polygon{
            toVec(input.vertices[face[0]]), toVec(input.vertices[face[1]]),
            toVec(input.vertices[face[2]])};
        std::vector<Vec3> clipped;
        for (std::size_t index = 0; index < polygon.size(); ++index)
        {
            const Vec3 first = polygon[index];
            const Vec3 second = polygon[(index + 1) % polygon.size()];
            const double first_distance = dot(normal, first) - offset;
            const double second_distance = dot(normal, second) - offset;
            const bool first_inside = first_distance <= tolerance;
            const bool second_inside = second_distance <= tolerance;
            if (first_inside) clipped.push_back(first);
            if (first_inside == second_inside) continue;
            const Vec3 intersection = first + (second - first) *
                (first_distance / (first_distance - second_distance));
            clipped.push_back(intersection);
            cap_points.push_back(intersection);
        }
        if (clipped.size() < 3) continue;
        const std::uint32_t first = appendVertex(clipped[0]);
        for (std::size_t index = 1; index + 1 < clipped.size(); ++index)
            result.triangles.push_back(
                {first, appendVertex(clipped[index]), appendVertex(clipped[index + 1])});
    }
    std::vector<std::uint32_t> cap_vertices;
    for (const Vec3 point : cap_points)
    {
        const std::uint32_t vertex = appendVertex(point);
        if (std::find(cap_vertices.begin(), cap_vertices.end(), vertex) == cap_vertices.end())
            cap_vertices.push_back(vertex);
    }
    if (cap_vertices.size() >= 3)
    {
        Vec3 center{};
        for (const std::uint32_t vertex : cap_vertices)
            center = center + toVec(result.vertices[vertex]);
        center = center * (1.0 / cap_vertices.size());
        const Vec3 helper = std::abs(normal.x) < 0.8 ? Vec3{1,0,0} : Vec3{0,1,0};
        const Vec3 raw_tangent = cross(helper, normal);
        const Vec3 tangent = raw_tangent * (1.0 / norm(raw_tangent));
        const Vec3 bitangent = cross(normal, tangent);
        std::sort(cap_vertices.begin(), cap_vertices.end(), [&](const std::uint32_t first,
                                                                 const std::uint32_t second)
        {
            const Vec3 a = toVec(result.vertices[first]) - center;
            const Vec3 b = toVec(result.vertices[second]) - center;
            return std::atan2(dot(a, bitangent), dot(a, tangent)) <
                   std::atan2(dot(b, bitangent), dot(b, tangent));
        });
        for (std::size_t index = 1; index + 1 < cap_vertices.size(); ++index)
            result.triangles.push_back(
                {cap_vertices[0], cap_vertices[index], cap_vertices[index + 1]});
    }
    orientPositiveVolume(result);
    return result;
}

MeshModel makeConvexHull(const MeshModel& reference)
{
    std::vector<quickhull::Vector3<double>> points;
    points.reserve(reference.vertices.size());
    for (const Position3 p : reference.vertices) points.emplace_back(p.x, p.y, p.z);
    std::sort(points.begin(),points.end(),[](const auto& first,const auto& second)
    { return std::tie(first.x,first.y,first.z)<std::tie(second.x,second.y,second.z); });
    points.erase(std::unique(points.begin(),points.end(),[](const auto& first,const auto& second)
    { return first.x==second.x && first.y==second.y && first.z==second.z; }),points.end());
    if (points.size()<4) throw std::runtime_error("convex hull input is lower dimensional");
    quickhull::QuickHull<double> builder;
    const auto hull = builder.getConvexHull(points, true, false, 0.0);
    if (builder.getDiagnostics().m_failedHorizonEdges != 0)
        throw std::runtime_error("QuickHull reported an unresolved horizon edge");
    MeshModel result;
    result.name = "convex_hull_of_all_phase1_vertices";
    for (const auto& p : hull.getVertexBuffer()) result.vertices.push_back({p.x, p.y, p.z});
    const auto& indices = hull.getIndexBuffer();
    for (std::size_t i = 0; i + 2 < indices.size(); i += 3)
        result.triangles.push_back({static_cast<std::uint32_t>(indices[i]),
                                    static_cast<std::uint32_t>(indices[i + 1]),
                                    static_cast<std::uint32_t>(indices[i + 2])});
    if (result.triangles.empty()) throw std::runtime_error("convex hull construction failed");
    Bounds hull_bounds;
    for (const Position3 point : result.vertices) include(hull_bounds,toVec(point));
    const double scale=norm(hull_bounds.upper-hull_bounds.lower);
    const double area_tolerance = std::max(1.0e-20,scale*scale*1.0e-14);
    Vec3 center{};
    for (const Position3 point : result.vertices) center = center + toVec(point);
    center = center * (1.0 / result.vertices.size());
    std::vector<TriangleIndices> valid_faces;
    std::set<TriangleIndices> signatures;
    for (TriangleIndices face : result.triangles)
    {
        const Vec3 a = toVec(result.vertices[face[0]]);
        Vec3 normal = cross(toVec(result.vertices[face[1]]) - a,
                            toVec(result.vertices[face[2]]) - a);
        if (norm(normal) <= area_tolerance) continue;
        if (dot(normal, center - a) > 0.0) std::swap(face[1], face[2]);
        TriangleIndices signature = face;
        std::sort(signature.begin(), signature.end());
        if (signatures.insert(signature).second) valid_faces.push_back(face);
    }
    result.triangles = std::move(valid_faces);
    if (result.triangles.empty())
        throw std::runtime_error("convex hull contains no non-degenerate faces");
    return result;
}

Vec3 normalized(Vec3 value);

MeshModel makeConvexHull(const MeshModel& reference,
                         const std::vector<std::uint32_t>& faces)
{
    MeshModel subset;
    subset.name = "convex_face_group";
    subset.vertices.reserve(faces.size() * 3);
    for (const std::uint32_t face_id : faces)
        for (const std::uint32_t vertex : reference.triangles[face_id])
            subset.vertices.push_back(reference.vertices[vertex]);
    return makeConvexHull(subset);
}

struct ConvexPiece
{
    MeshModel mesh;
    Bounds bounds;
    std::vector<Vec3> normals;
    std::vector<double> offsets;
    std::vector<std::uint32_t> responsibility;
};

ConvexPiece makeConvexPiece(MeshModel mesh, std::vector<std::uint32_t> responsibility)
{
    ConvexPiece piece;
    piece.mesh = std::move(mesh);
    piece.responsibility = std::move(responsibility);
    Vec3 center{};
    for (const Position3 point : piece.mesh.vertices)
    {
        include(piece.bounds, toVec(point));
        center = center + toVec(point);
    }
    center = center * (1.0 / piece.mesh.vertices.size());
    piece.normals.reserve(piece.mesh.triangles.size());
    piece.offsets.reserve(piece.mesh.triangles.size());
    for (TriangleIndices& face : piece.mesh.triangles)
    {
        const Vec3 a = toVec(piece.mesh.vertices[face[0]]);
        const Vec3 b = toVec(piece.mesh.vertices[face[1]]);
        const Vec3 c = toVec(piece.mesh.vertices[face[2]]);
        Vec3 normal = normalized(cross(b-a, c-a));
        if (dot(normal, center-a) > 0.0)
        {
            std::swap(face[1], face[2]);
            normal = normal * -1.0;
        }
        piece.normals.push_back(normal);
        piece.offsets.push_back(dot(normal, a));
    }
    return piece;
}

bool boundsOverlap(const Bounds& first, const Bounds& second, const double tolerance)
{
    return first.upper.x >= second.lower.x-tolerance &&
           second.upper.x >= first.lower.x-tolerance &&
           first.upper.y >= second.lower.y-tolerance &&
           second.upper.y >= first.lower.y-tolerance &&
           first.upper.z >= second.lower.z-tolerance &&
           second.upper.z >= first.lower.z-tolerance;
}

std::vector<Vec3> clipPolygonInsideConvexPiece(
    std::vector<Vec3> polygon, const ConvexPiece& piece, const double tolerance)
{
    for (std::size_t plane = 0; plane < piece.normals.size() && !polygon.empty(); ++plane)
    {
        std::vector<Vec3> clipped;
        clipped.reserve(polygon.size()+2);
        for (std::size_t index = 0; index < polygon.size(); ++index)
        {
            const Vec3 first = polygon[index];
            const Vec3 second = polygon[(index+1)%polygon.size()];
            const double first_distance = dot(piece.normals[plane], first)-piece.offsets[plane];
            const double second_distance = dot(piece.normals[plane], second)-piece.offsets[plane];
            const bool first_inside = first_distance <= tolerance;
            const bool second_inside = second_distance <= tolerance;
            if (first_inside) clipped.push_back(first);
            if (first_inside == second_inside) continue;
            const double denominator = first_distance-second_distance;
            if (std::abs(denominator) <= 1.0e-30) continue;
            clipped.push_back(first+(second-first)*
                std::clamp(first_distance/denominator, 0.0, 1.0));
        }
        polygon = std::move(clipped);
    }
    return polygon;
}

MeshModel convexUnionOuterSurface(
    const std::vector<ConvexPiece>& pieces, const MeshModel& exact_fallback,
    const double scale)
{
    constexpr int precision = 8;
    const double tolerance = std::max(1.0e-10, scale*1.0e-10);
    MeshModel result;
    result.name = "adaptive_certified_convex_cover_exposed_boundary";
    const auto appendTriangle = [&](const Vec3 a, const Vec3 b, const Vec3 c)
    {
        if (normSquared(cross(b-a,c-a)) <= tolerance*tolerance) return;
        const std::uint32_t base = static_cast<std::uint32_t>(result.vertices.size());
        result.vertices.push_back(toPosition(a));
        result.vertices.push_back(toPosition(b));
        result.vertices.push_back(toPosition(c));
        result.triangles.push_back({base,base+1,base+2});
    };
    for (std::size_t piece_id = 0; piece_id < pieces.size(); ++piece_id)
    {
        const ConvexPiece& piece = pieces[piece_id];
        for (std::size_t face_id = 0; face_id < piece.mesh.triangles.size(); ++face_id)
        {
            const TriangleIndices face = piece.mesh.triangles[face_id];
            const Vec3 origin = toVec(piece.mesh.vertices[face[0]]);
            const Vec3 normal = piece.normals[face_id];
            const Vec3 basis_u = std::abs(normal.x) < 0.8
                ? normalized(cross(normal,{1,0,0}))
                : normalized(cross(normal,{0,1,0}));
            const Vec3 basis_v = cross(normal,basis_u);
            const auto project = [&](const Vec3 point)
            {
                const Vec3 relative = point-origin;
                return Clipper2Lib::PointD(dot(relative,basis_u),dot(relative,basis_v));
            };
            Clipper2Lib::PathD original;
            for (const std::uint32_t vertex : face)
                original.push_back(project(toVec(piece.mesh.vertices[vertex])));
            if (Clipper2Lib::Area(original)<0.0) std::reverse(original.begin(),original.end());
            Clipper2Lib::PathsD remaining{original};
            for (std::size_t other_id = 0; other_id < pieces.size() && !remaining.empty(); ++other_id)
            {
                if (piece_id == other_id || !boundsOverlap(piece.bounds,pieces[other_id].bounds,tolerance))
                    continue;
                const ConvexPiece& other = pieces[other_id];
                bool same_facing_coplanar = false;
                for (std::size_t other_face = 0; other_face < other.normals.size(); ++other_face)
                    if (dot(normal,other.normals[other_face]) > 1.0-1.0e-10 &&
                        std::abs(dot(normal,toVec(other.mesh.vertices[
                            other.mesh.triangles[other_face][0]]))-piece.offsets[face_id]) <= tolerance*8.0)
                    {
                        same_facing_coplanar = true;
                        break;
                    }
                if (same_facing_coplanar && piece_id < other_id) continue;
                std::vector<Vec3> intersection = clipPolygonInsideConvexPiece(
                    {toVec(piece.mesh.vertices[face[0]]),
                     toVec(piece.mesh.vertices[face[1]]),
                     toVec(piece.mesh.vertices[face[2]])}, other, tolerance*8.0);
                if (intersection.size()<3) continue;
                Clipper2Lib::PathD clip;
                for (const Vec3 point : intersection) clip.push_back(project(point));
                if (std::abs(Clipper2Lib::Area(clip))<=tolerance*tolerance) continue;
                if (Clipper2Lib::Area(clip)<0.0) std::reverse(clip.begin(),clip.end());
                remaining = Clipper2Lib::Difference(
                    remaining,Clipper2Lib::PathsD{clip},Clipper2Lib::FillRule::NonZero,precision);
            }
            Clipper2Lib::PathsD triangles;
            if (!remaining.empty() && Clipper2Lib::Triangulate(
                    remaining,precision,triangles,false)==Clipper2Lib::TriangulateResult::success)
                for (const auto& triangle : triangles)
                    if (triangle.size()==3)
                        appendTriangle(origin+basis_u*triangle[0].x+basis_v*triangle[0].y,
                                       origin+basis_u*triangle[1].x+basis_v*triangle[1].y,
                                       origin+basis_u*triangle[2].x+basis_v*triangle[2].y);
        }
    }
    // Exact fallback triangles are zero-error surface responsibility. Remove
    // only the portions certified inside an accepted closed convex piece.
    for (const TriangleIndices face : exact_fallback.triangles)
    {
        const Vec3 origin=toVec(exact_fallback.vertices[face[0]]);
        Vec3 normal=cross(toVec(exact_fallback.vertices[face[1]])-origin,
                          toVec(exact_fallback.vertices[face[2]])-origin);
        if (norm(normal)<=tolerance*tolerance) continue;
        normal=normalized(normal);
        const Vec3 basis_u=std::abs(normal.x)<0.8
            ? normalized(cross(normal,{1,0,0}))
            : normalized(cross(normal,{0,1,0}));
        const Vec3 basis_v=cross(normal,basis_u);
        const auto project=[&](const Vec3 point)
        {
            const Vec3 relative=point-origin;
            return Clipper2Lib::PointD(dot(relative,basis_u),dot(relative,basis_v));
        };
        Clipper2Lib::PathD original;
        for (const std::uint32_t vertex:face)
            original.push_back(project(toVec(exact_fallback.vertices[vertex])));
        if (Clipper2Lib::Area(original)<0.0) std::reverse(original.begin(),original.end());
        Clipper2Lib::PathsD remaining{original};
        Bounds face_bounds;
        for (const std::uint32_t vertex:face) include(face_bounds,toVec(exact_fallback.vertices[vertex]));
        for (const ConvexPiece& piece:pieces)
        {
            if (!boundsOverlap(face_bounds,piece.bounds,tolerance)||remaining.empty()) continue;
            std::vector<Vec3> intersection=clipPolygonInsideConvexPiece(
                {toVec(exact_fallback.vertices[face[0]]),
                 toVec(exact_fallback.vertices[face[1]]),
                 toVec(exact_fallback.vertices[face[2]])},piece,tolerance*8.0);
            if (intersection.size()<3) continue;
            Clipper2Lib::PathD clip;
            for (const Vec3 point:intersection) clip.push_back(project(point));
            if (std::abs(Clipper2Lib::Area(clip))<=tolerance*tolerance) continue;
            if (Clipper2Lib::Area(clip)<0.0) std::reverse(clip.begin(),clip.end());
            remaining=Clipper2Lib::Difference(
                remaining,Clipper2Lib::PathsD{clip},Clipper2Lib::FillRule::NonZero,precision);
        }
        Clipper2Lib::PathsD triangles;
        if (!remaining.empty()&&Clipper2Lib::Triangulate(
                remaining,precision,triangles,false)==Clipper2Lib::TriangulateResult::success)
            for (const auto& triangle:triangles)
                if (triangle.size()==3)
                    appendTriangle(origin+basis_u*triangle[0].x+basis_v*triangle[0].y,
                                   origin+basis_u*triangle[1].x+basis_v*triangle[1].y,
                                   origin+basis_u*triangle[2].x+basis_v*triangle[2].y);
    }
    return canonicalizeCoplanarTriangleSoup(result);
}

MeshModel makeOrientedBoxFromPoints(const std::vector<Vec3>& points)
{
    Vec3 center{};
    if (points.empty()) throw std::runtime_error("oriented box has no responsibility");
    for (const Vec3 point : points) center = center + point;
    center = center * (1.0 / points.size());
    double covariance[3][3]{};
    for (const Vec3 point : points)
    {
        const Vec3 delta = point - center;
        const double value[3]{delta.x, delta.y, delta.z};
        for (int row = 0; row < 3; ++row)
            for (int column = 0; column < 3; ++column)
                covariance[row][column] += value[row] * value[column];
    }
    const auto multiply = [&](const Vec3 value)
    {
        return Vec3{
            covariance[0][0] * value.x + covariance[0][1] * value.y + covariance[0][2] * value.z,
            covariance[1][0] * value.x + covariance[1][1] * value.y + covariance[1][2] * value.z,
            covariance[2][0] * value.x + covariance[2][1] * value.y + covariance[2][2] * value.z};
    };
    const auto powerAxis = [&](Vec3 axis, const Vec3 excluded)
    {
        for (int iteration = 0; iteration < 32; ++iteration)
        {
            Vec3 next = multiply(axis);
            next = next - excluded * dot(next, excluded);
            if (!(norm(next) > 1.0e-30)) break;
            axis = normalized(next);
        }
        return normalized(axis - excluded * dot(axis, excluded));
    };
    int first_dimension = 0;
    if (covariance[1][1] > covariance[first_dimension][first_dimension]) first_dimension = 1;
    if (covariance[2][2] > covariance[first_dimension][first_dimension]) first_dimension = 2;
    Vec3 first = first_dimension == 0 ? Vec3{1,0,0}
               : first_dimension == 1 ? Vec3{0,1,0} : Vec3{0,0,1};
    first = powerAxis(first, {});
    Vec3 second_seed = std::abs(first.x) < 0.8 ? Vec3{1,0,0} : Vec3{0,1,0};
    Vec3 second = powerAxis(second_seed, first);
    Vec3 third = normalized(cross(first, second));
    second = cross(third, first);
    const Vec3 axes[3]{first, second, third};
    double lower[3]{std::numeric_limits<double>::infinity(),
                    std::numeric_limits<double>::infinity(),
                    std::numeric_limits<double>::infinity()};
    double upper[3]{-std::numeric_limits<double>::infinity(),
                    -std::numeric_limits<double>::infinity(),
                    -std::numeric_limits<double>::infinity()};
    for (const Vec3 point : points)
        for (int axis = 0; axis < 3; ++axis)
        {
            const double coordinate = dot(axes[axis], point);
            lower[axis] = std::min(lower[axis], coordinate);
            upper[axis] = std::max(upper[axis], coordinate);
        }
    MeshModel box;
    box.name = "oriented_responsibility_box";
    for (int x = 0; x < 2; ++x)
        for (int y = 0; y < 2; ++y)
            for (int z = 0; z < 2; ++z)
            {
                const Vec3 point = axes[0] * (x ? upper[0] : lower[0]) +
                                   axes[1] * (y ? upper[1] : lower[1]) +
                                   axes[2] * (z ? upper[2] : lower[2]);
                box.vertices.push_back(toPosition(point));
            }
    box.triangles = {
        {0,2,3},{0,3,1}, {4,5,7},{4,7,6},
        {0,1,5},{0,5,4}, {2,6,7},{2,7,3},
        {0,4,6},{0,6,2}, {1,3,7},{1,7,5}};
    return box;
}

MeshModel makeOrientedBox(const MeshModel& reference,
                          const std::vector<std::uint32_t>& faces)
{
    std::vector<Vec3> points;
    points.reserve(faces.size() * 3);
    for (const std::uint32_t face_id : faces)
        for (const std::uint32_t vertex : reference.triangles[face_id])
            points.push_back(toVec(reference.vertices[vertex]));
    return makeOrientedBoxFromPoints(points);
}

MeshModel makeOrientedTriangularPrism(const MeshModel& reference,
                                      const std::vector<std::uint32_t>& faces,
                                      const bool swap_cross_axes)
{
    const MeshModel box = makeOrientedBox(reference, faces);
    const Vec3 corner = toVec(box.vertices[0]);
    std::array<Vec3, 3> axes{
        normalized(toVec(box.vertices[4]) - corner),
        normalized(toVec(box.vertices[2]) - corner),
        normalized(toVec(box.vertices[1]) - corner)};
    std::array<double, 3> lower{};
    std::array<double, 3> upper{};
    for (int axis = 0; axis < 3; ++axis)
    {
        lower[axis] = std::numeric_limits<double>::infinity();
        upper[axis] = -std::numeric_limits<double>::infinity();
        for (const Position3 point : box.vertices)
        {
            const double coordinate = dot(axes[axis], toVec(point));
            lower[axis] = std::min(lower[axis], coordinate);
            upper[axis] = std::max(upper[axis], coordinate);
        }
    }
    int extrusion = 0;
    if (upper[1] - lower[1] > upper[extrusion] - lower[extrusion]) extrusion = 1;
    if (upper[2] - lower[2] > upper[extrusion] - lower[extrusion]) extrusion = 2;
    std::array<int, 2> cross_axes{};
    int cross_count = 0;
    for (int axis = 0; axis < 3; ++axis)
        if (axis != extrusion) cross_axes[cross_count++] = axis;
    if (swap_cross_axes) std::swap(cross_axes[0], cross_axes[1]);
    const int u = cross_axes[0];
    const int v = cross_axes[1];
    const double width = upper[u] - lower[u];
    const double height = upper[v] - lower[v];
    const std::array<std::array<double, 2>, 3> triangle{{
        {lower[u] - 0.5 * width, lower[v]},
        {upper[u] + 0.5 * width, lower[v]},
        {0.5 * (lower[u] + upper[u]), lower[v] + 2.0 * height}}};
    MeshModel prism;
    prism.name = "oriented_responsibility_triangular_prism";
    for (int side = 0; side < 2; ++side)
        for (const auto coordinate : triangle)
            prism.vertices.push_back(toPosition(
                axes[u] * coordinate[0] + axes[v] * coordinate[1] +
                axes[extrusion] * (side ? upper[extrusion] : lower[extrusion])));
    prism.triangles = {
        {0,2,1}, {3,4,5},
        {0,1,4},{0,4,3}, {1,2,5},{1,5,4}, {2,0,3},{2,3,5}};
    return prism;
}

MeshModel makeOrientedTriangularPrismFromPoints(
    const std::vector<Vec3>& points, const bool swap_cross_axes)
{
    const MeshModel box = makeOrientedBoxFromPoints(points);
    MeshModel box_reference;
    box_reference.vertices = box.vertices;
    box_reference.triangles = box.triangles;
    std::vector<std::uint32_t> faces(box.triangles.size());
    std::iota(faces.begin(), faces.end(), std::uint32_t{0});
    return makeOrientedTriangularPrism(box_reference, faces, swap_cross_axes);
}

MeshModel makeOrientedTetrahedron(const MeshModel& reference,
                                  const std::vector<std::uint32_t>& faces)
{
    const MeshModel box = makeOrientedBox(reference, faces);
    const Vec3 corner = toVec(box.vertices[0]);
    const std::array<Vec3, 3> axes{
        normalized(toVec(box.vertices[4]) - corner),
        normalized(toVec(box.vertices[2]) - corner),
        normalized(toVec(box.vertices[1]) - corner)};
    Vec3 center{};
    for (const Position3 point : box.vertices) center = center + toVec(point);
    center = center * (1.0 / box.vertices.size());
    const std::array<double, 3> half_extent{
        0.5 * norm(toVec(box.vertices[4]) - corner),
        0.5 * norm(toVec(box.vertices[2]) - corner),
        0.5 * norm(toVec(box.vertices[1]) - corner)};
    MeshModel tetrahedron;
    tetrahedron.name = "oriented_responsibility_tetrahedron";
    for (const std::array<double, 3> signs :
         {std::array<double,3>{1,1,1}, {1,-1,-1}, {-1,1,-1}, {-1,-1,1}})
        tetrahedron.vertices.push_back(toPosition(
            center + axes[0] * (3.0 * signs[0] * half_extent[0]) +
                     axes[1] * (3.0 * signs[1] * half_extent[1]) +
                     axes[2] * (3.0 * signs[2] * half_extent[2])));
    tetrahedron.triangles = {{0,2,1},{0,1,3},{0,3,2},{1,2,3}};
    return tetrahedron;
}

MeshModel boundedBoxCover(const MeshModel& reference, const std::size_t target_pieces,
                          const double split_fraction, const std::uint32_t clip_iterations = 0,
                          const double maximum_distance = 0.0,
                          const std::uint32_t maximum_depth = 18,
                          const double numerical_tolerance = 1.0e-9,
                          const bool triangular_prisms = false,
                          const bool swap_prism_axes = false,
                          const double retreat_factor = 1.75,
                          const std::uint32_t prism_group_mask = 0,
                          const bool tetrahedra = false,
                          const double retreat_margin = 0.0,
                          const std::uint32_t final_clip_depth = 0,
                          const std::uint32_t final_clip_iterations = 0,
                          const std::uint32_t witness_split_budget = 0)
{
    struct Group
    {
        std::vector<std::uint32_t> faces;
        std::vector<Vec3> responsibility_points;
    };
    Group root;
    root.faces.resize(reference.triangles.size());
    std::iota(root.faces.begin(), root.faces.end(), std::uint32_t{0});
    std::vector<Group> groups;
    groups.push_back(std::move(root));
    const auto responsibilityPoints = [&](const Group& group)
    {
        if (!group.responsibility_points.empty()) return group.responsibility_points;
        std::vector<Vec3> points;
        points.reserve(group.faces.size() * 3);
        for (const std::uint32_t face_id : group.faces)
            for (const std::uint32_t vertex : reference.triangles[face_id])
                points.push_back(toVec(reference.vertices[vertex]));
        return points;
    };
    while (groups.size() < target_pieces)
    {
        std::size_t selected = groups.size();
        double selected_extent = -1.0;
        Vec3 selected_axis{};
        bool witness_split = false;
        double witness_split_coordinate = 0.0;
        if (target_pieces == 9 && groups.size() == 8 && maximum_distance > 0.0)
        {
            MeshModel eight_boxes;
            std::vector<MeshModel> group_boxes;
            group_boxes.reserve(groups.size());
            for (const Group& group : groups)
            {
                group_boxes.push_back(makeOrientedBox(reference, group.faces));
                const MeshModel& box = group_boxes.back();
                const std::uint32_t base = static_cast<std::uint32_t>(eight_boxes.vertices.size());
                eight_boxes.vertices.insert(
                    eight_boxes.vertices.end(), box.vertices.begin(), box.vertices.end());
                for (const TriangleIndices face : box.triangles)
                    eight_boxes.triangles.push_back(
                        {base + face[0], base + face[1], base + face[2]});
            }
            const TriangleBvh reference_bvh(reference);
            const DirectedHausdorffCertificate certificate = certifyWithBvh(
                eight_boxes, reference_bvh, maximum_distance, maximum_depth,
                numerical_tolerance, false);
            const Vec3 witness{certificate.maximum_proxy_point[0],
                               certificate.maximum_proxy_point[1],
                               certificate.maximum_proxy_point[2]};
            const Vec3 nearest_reference{certificate.nearest_reference_point[0],
                                         certificate.nearest_reference_point[1],
                                         certificate.nearest_reference_point[2]};
            double nearest = std::numeric_limits<double>::infinity();
            for (std::size_t index = 0; index < group_boxes.size(); ++index)
            {
                const TriangleBvh box_bvh(group_boxes[index]);
                const double distance = box_bvh.nearest(witness).first;
                if (distance < nearest)
                {
                    nearest = distance;
                    selected = index;
                }
            }
            const Vec3 separation = witness - nearest_reference;
            if (selected != groups.size() && norm(separation) > 1.0e-12)
            {
                selected_axis = normalized(separation);
                witness_split_coordinate = dot(selected_axis, witness);
                witness_split = true;
            }
        }
        for (std::size_t index = 0; index < groups.size(); ++index)
        {
            if (selected != groups.size() && index != selected) continue;
            if (groups[index].faces.size() < 2) continue;
            if (witness_split)
            {
                selected_extent = 0.0;
                break;
            }
            Bounds centers;
            for (const std::uint32_t face_id : groups[index].faces)
            {
                const auto face = reference.triangles[face_id];
                const Vec3 center = (toVec(reference.vertices[face[0]]) +
                                     toVec(reference.vertices[face[1]]) +
                                     toVec(reference.vertices[face[2]])) * (1.0 / 3.0);
                include(centers, center);
            }
            const Vec3 extent = centers.upper - centers.lower;
            const int axis = extent.y > extent.x ? (extent.z > extent.y ? 2 : 1)
                                                 : (extent.z > extent.x ? 2 : 0);
            const double value = axis == 0 ? extent.x : axis == 1 ? extent.y : extent.z;
            if (value > selected_extent)
            {
                selected = index;
                selected_extent = value;
                selected_axis = axis == 0 ? Vec3{1,0,0}
                              : axis == 1 ? Vec3{0,1,0} : Vec3{0,0,1};
            }
        }
        if (selected == groups.size()) break;
        Group group = std::move(groups[selected]);
        if (witness_split)
        {
            Group first, second;
            const auto appendClipped = [&](std::vector<Vec3> polygon, const bool lower)
            {
                std::vector<Vec3> clipped;
                for (std::size_t index = 0; index < polygon.size(); ++index)
                {
                    const Vec3 a = polygon[index];
                    const Vec3 b = polygon[(index + 1) % polygon.size()];
                    const double da = dot(selected_axis, a) - witness_split_coordinate;
                    const double db = dot(selected_axis, b) - witness_split_coordinate;
                    const bool a_inside = lower ? da <= 0.0 : da >= 0.0;
                    const bool b_inside = lower ? db <= 0.0 : db >= 0.0;
                    if (a_inside) clipped.push_back(a);
                    if (a_inside != b_inside)
                        clipped.push_back(a + (b - a) * (da / (da - db)));
                }
                auto& target = lower ? first.responsibility_points
                                     : second.responsibility_points;
                target.insert(target.end(), clipped.begin(), clipped.end());
            };
            for (const std::uint32_t face_id : group.faces)
            {
                const TriangleIndices face = reference.triangles[face_id];
                const std::vector<Vec3> triangle{
                    toVec(reference.vertices[face[0]]),
                    toVec(reference.vertices[face[1]]),
                    toVec(reference.vertices[face[2]])};
                appendClipped(triangle, true);
                appendClipped(triangle, false);
            }
            if (!first.responsibility_points.empty() &&
                !second.responsibility_points.empty())
            {
                groups[selected] = std::move(first);
                groups.push_back(std::move(second));
                continue;
            }
        }
        const auto coordinate = [&](const std::uint32_t face_id)
        {
            const auto face = reference.triangles[face_id];
            if (witness_split)
                return std::max({dot(selected_axis, toVec(reference.vertices[face[0]])),
                                 dot(selected_axis, toVec(reference.vertices[face[1]])),
                                 dot(selected_axis, toVec(reference.vertices[face[2]]))});
            const Vec3 center = (toVec(reference.vertices[face[0]]) +
                                 toVec(reference.vertices[face[1]]) +
                                 toVec(reference.vertices[face[2]])) * (1.0 / 3.0);
            return dot(selected_axis, center);
        };
        std::size_t middle = 0;
        if (witness_split)
        {
            const auto split = std::stable_partition(
                group.faces.begin(), group.faces.end(), [&](const std::uint32_t face_id)
                { return coordinate(face_id) < witness_split_coordinate; });
            middle = static_cast<std::size_t>(split - group.faces.begin());
        }
        if (middle == 0 || middle == group.faces.size())
        {
            middle = std::clamp<std::size_t>(
                static_cast<std::size_t>(group.faces.size() * split_fraction),
                1, group.faces.size() - 1);
            std::nth_element(group.faces.begin(), group.faces.begin() + middle,
                             group.faces.end(), [&](const std::uint32_t first,
                                                    const std::uint32_t second)
                             { return coordinate(first) < coordinate(second); });
        }
        Group first, second;
        first.faces.assign(group.faces.begin(), group.faces.begin() + middle);
        second.faces.assign(group.faces.begin() + middle, group.faces.end());
        groups[selected] = std::move(first);
        groups.push_back(std::move(second));
    }

    std::vector<MeshModel> boxes;
    std::vector<HalfspacePolytope> halfspaces;
    struct InternalPlane { Vec3 normal{}; double offset = 0.0; bool enabled = false; };
    std::vector<InternalPlane> internal_planes;
    boxes.reserve(groups.size());
    halfspaces.reserve(groups.size());
    internal_planes.reserve(groups.size());
    for (std::size_t group_id = 0; group_id < groups.size(); ++group_id)
    {
        const std::vector<Vec3> points = responsibilityPoints(groups[group_id]);
        boxes.push_back(tetrahedra
            ? makeOrientedTetrahedron(reference, groups[group_id].faces)
            : triangular_prisms ||
                        (prism_group_mask & (std::uint32_t{1} << group_id)) != 0
            ? makeOrientedTriangularPrismFromPoints(points, swap_prism_axes)
            : makeOrientedBoxFromPoints(points));
        orientPositiveVolume(boxes.back());
        halfspaces.push_back(extractHalfspaces(boxes.back()));
        internal_planes.push_back({});
    }
    const auto assemble = [&]()
    {
        MeshModel result;
        result.name = "bounded_partition_oriented_box_surfaces";
        for (std::size_t box_id = 0; box_id < boxes.size(); ++box_id)
        {
            const MeshModel& box = boxes[box_id];
            const std::uint32_t base = static_cast<std::uint32_t>(result.vertices.size());
            result.vertices.insert(result.vertices.end(), box.vertices.begin(), box.vertices.end());
            for (const TriangleIndices face : box.triangles)
            {
                if (internal_planes[box_id].enabled)
                {
                    bool internal = true;
                    for (const std::uint32_t vertex : face)
                        internal &= std::abs(dot(internal_planes[box_id].normal,
                            toVec(box.vertices[vertex])) -
                            internal_planes[box_id].offset) <= 1.0e-7;
                    if (internal) continue;
                }
                result.triangles.push_back({base + face[0], base + face[1], base + face[2]});
            }
        }
        return result;
    };
    std::uint32_t witness_splits = 0;
    std::vector<std::size_t> strict_plane_starts(boxes.size(),
        std::numeric_limits<std::size_t>::max());
    for (std::uint32_t clip_iteration = 0;
         clip_iteration < clip_iterations && !boxes.empty(); ++clip_iteration)
    {
        const TriangleBvh reference_bvh(reference);
        const std::uint32_t certificate_depth =
            final_clip_depth != 0 &&
                    clip_iteration + final_clip_iterations >= clip_iterations
                ? final_clip_depth : maximum_depth;
        const DirectedHausdorffCertificate raw_certificate = certifyWithBvh(
            assemble(), reference_bvh, maximum_distance, certificate_depth,
            numerical_tolerance, true, certificate_depth == 0);
        if (raw_certificate.passed) break;
        const Vec3 witness{raw_certificate.maximum_proxy_point[0],
                           raw_certificate.maximum_proxy_point[1],
                           raw_certificate.maximum_proxy_point[2]};
        std::size_t worst_box = boxes.size();
        std::size_t first_face = 0;
        for (std::size_t box_id = 0; box_id < boxes.size(); ++box_id)
        {
            const std::size_t end_face = first_face + boxes[box_id].triangles.size();
            if (raw_certificate.maximum_proxy_face < end_face)
            {
                worst_box = box_id;
                break;
            }
            first_face = end_face;
        }
        if (worst_box == boxes.size())
            throw std::runtime_error("Hausdorff witness has no proxy component owner");
        MeshModel& box = boxes[worst_box];
        const Vec3 nearest_reference{raw_certificate.nearest_reference_point[0],
                                     raw_certificate.nearest_reference_point[1],
                                     raw_certificate.nearest_reference_point[2]};
        Vec3 center{};
        for (const Position3 point : box.vertices) center = center + toVec(point);
        center = center * (1.0 / box.vertices.size());
        const Vec3 direction = normalized(witness - nearest_reference);
        double support = -std::numeric_limits<double>::infinity();
        for (const Vec3 point : responsibilityPoints(groups[worst_box]))
            support = std::max(support, dot(direction, point));
        const double witness_offset = dot(direction, witness);
        const double effective_retreat_factor =
            certificate_depth == final_clip_depth && final_clip_depth != 0
                ? retreat_factor * 1.7 : retreat_factor;
        const double effective_retreat_margin =
            certificate_depth == final_clip_depth && final_clip_depth != 0
                ? retreat_margin : 0.0;
        const double retreat = std::max(0.0, raw_certificate.lower_bound - maximum_distance) *
                               effective_retreat_factor + effective_retreat_margin;
        const double clip_offset = std::max(support, witness_offset - retreat);
        if (certificate_depth == final_clip_depth && final_clip_depth != 0)
            std::cerr << std::setprecision(17)
                      << "monitor: stage=deep_witness iteration=" << clip_iteration
                      << " box=" << worst_box
                      << " triangles=" << assemble().triangles.size()
                      << " lower_bound=" << raw_certificate.lower_bound
                      << " witness_offset=" << witness_offset
                      << " responsibility_support=" << support
                      << " clip_offset=" << clip_offset
                      << " retreat=" << retreat << '\n';
        if (clip_offset >= witness_offset &&
            certificate_depth == final_clip_depth &&
            witness_splits < witness_split_budget &&
            groups[worst_box].faces.size() >= 2)
        {
            const Vec3 split_axis = direction;
            const double split_coordinate = witness_offset;
            const auto faceCoordinate = [&](const std::uint32_t face_id)
            {
                const TriangleIndices face = reference.triangles[face_id];
                return dot(split_axis,
                    (toVec(reference.vertices[face[0]]) +
                     toVec(reference.vertices[face[1]]) +
                     toVec(reference.vertices[face[2]])) * (1.0 / 3.0));
            };
            auto& responsibility = groups[worst_box].faces;
            auto split = std::stable_partition(
                responsibility.begin(), responsibility.end(),
                [&](const std::uint32_t face_id)
                { return faceCoordinate(face_id) < split_coordinate; });
            std::size_t middle = static_cast<std::size_t>(split - responsibility.begin());
            if (middle == 0 || middle == responsibility.size())
            {
                middle = responsibility.size() / 2;
                std::nth_element(responsibility.begin(), responsibility.begin() + middle,
                                 responsibility.end(), [&](const std::uint32_t first,
                                                            const std::uint32_t second)
                                 { return faceCoordinate(first) < faceCoordinate(second); });
            }
            Group second_group;
            second_group.faces.assign(responsibility.begin() + middle, responsibility.end());
            responsibility.erase(responsibility.begin() + middle, responsibility.end());
            groups.push_back(std::move(second_group));
            double first_support = -std::numeric_limits<double>::infinity();
            double second_lower = std::numeric_limits<double>::infinity();
            for (const Vec3 point : responsibilityPoints(groups[worst_box]))
                first_support = std::max(first_support, dot(split_axis, point));
            for (const Vec3 point : responsibilityPoints(groups.back()))
                second_lower = std::min(second_lower, dot(split_axis, point));
            const double split_overlap = std::max(1.0e-8,
                maximum_distance * 1.0e-6);
            HalfspacePolytope second_halfspaces = halfspaces[worst_box];
            halfspaces[worst_box].normals.push_back(split_axis);
            halfspaces[worst_box].offsets.push_back(first_support + split_overlap);
            second_halfspaces.normals.push_back(split_axis * -1.0);
            second_halfspaces.offsets.push_back(-second_lower + split_overlap);
            boxes[worst_box] = buildHalfspacePolytope(halfspaces[worst_box]);
            boxes.push_back(buildHalfspacePolytope(second_halfspaces));
            halfspaces.push_back(std::move(second_halfspaces));
            internal_planes[worst_box] = {};
            internal_planes.push_back({});
            ++witness_splits;
            continue;
        }
        if (clip_offset < witness_offset)
        {
            if (certificate_depth == final_clip_depth && final_clip_depth != 0 &&
                strict_plane_starts[worst_box] == std::numeric_limits<std::size_t>::max())
                strict_plane_starts[worst_box] = halfspaces[worst_box].normals.size();
            HalfspacePolytope trial_halfspaces = halfspaces[worst_box];
            trial_halfspaces.normals.push_back(direction);
            trial_halfspaces.offsets.push_back(clip_offset);
            MeshModel trial_box = buildHalfspacePolytope(trial_halfspaces);
            MeshModel previous_box = box;
            box = std::move(trial_box);
            const MeshModel trial = assemble();
            bool accept = trial.triangles.size() <= 100;
            if (certificate_depth == final_clip_depth && final_clip_depth != 0)
            {
                HalfspacePolytope tightened = halfspaces[worst_box];
                const double tightening =
                    std::max(0.0, raw_certificate.lower_bound - maximum_distance) * 4.0 +
                    effective_retreat_margin;
                for (std::size_t plane = strict_plane_starts[worst_box];
                     plane < tightened.normals.size(); ++plane)
                {
                    double plane_support = -std::numeric_limits<double>::infinity();
                    for (const Vec3 point : responsibilityPoints(groups[worst_box]))
                        plane_support = std::max(plane_support,
                            dot(tightened.normals[plane], point));
                    tightened.offsets[plane] = std::max(
                        plane_support, tightened.offsets[plane] - tightening);
                }
                MeshModel tightened_box = buildHalfspacePolytope(tightened);
                box = tightened_box;
                const MeshModel tightened_mesh = assemble();
                if (tightened_mesh.triangles.size() <= 100)
                {
                    const DirectedHausdorffCertificate tightened_certificate = certifyWithBvh(
                        tightened_mesh, reference_bvh, maximum_distance,
                        certificate_depth, numerical_tolerance);
                    if (tightened_certificate.passed)
                    {
                        trial_halfspaces = std::move(tightened);
                        box = std::move(tightened_box);
                        accept = true;
                    }
                }
            }
            if (!accept && certificate_depth == final_clip_depth && final_clip_depth != 0)
            {
                DirectedHausdorffCertificate best_certificate;
                best_certificate.lower_bound = std::numeric_limits<double>::infinity();
                std::optional<HalfspacePolytope> best_halfspaces;
                MeshModel best_box;
                for (std::size_t remove = strict_plane_starts[worst_box];
                     remove < halfspaces[worst_box].normals.size(); ++remove)
                {
                    HalfspacePolytope replacement = halfspaces[worst_box];
                    replacement.normals.erase(replacement.normals.begin() + remove);
                    replacement.offsets.erase(replacement.offsets.begin() + remove);
                    replacement.normals.push_back(direction);
                    replacement.offsets.push_back(clip_offset);
                    MeshModel replacement_box = buildHalfspacePolytope(replacement);
                    box = replacement_box;
                    const MeshModel replacement_mesh = assemble();
                    if (replacement_mesh.triangles.size() > 100) continue;
                    const DirectedHausdorffCertificate replacement_certificate = certifyWithBvh(
                        replacement_mesh, reference_bvh, maximum_distance,
                        certificate_depth, numerical_tolerance);
                    if (replacement_certificate.passed ||
                        replacement_certificate.lower_bound < best_certificate.lower_bound)
                    {
                        best_certificate = replacement_certificate;
                        best_halfspaces = std::move(replacement);
                        best_box = std::move(replacement_box);
                        if (replacement_certificate.passed) break;
                    }
                }
                if (best_halfspaces)
                {
                    trial_halfspaces = std::move(*best_halfspaces);
                    box = std::move(best_box);
                    accept = true;
                }
            }
            if (accept)
                halfspaces[worst_box] = std::move(trial_halfspaces);
            else
                box = std::move(previous_box);
        }
    }

    MeshModel result = assemble();
    if (witness_splits != 0)
        result = canonicalizeCoplanarTriangleSoup(result);
    return result;
}

MeshModel adaptiveConvexCover(
    const MeshModel& reference, const double maximum_distance,
    const std::uint32_t maximum_depth, const double numerical_tolerance)
{
    const TriangleBvh reference_bvh(reference);
    std::vector<ConvexPiece> accepted_pieces;
    MeshModel exact_fallback;
    exact_fallback.name = "adaptive_convex_exact_fallback";
    struct Group { std::vector<std::uint32_t> faces; std::uint32_t depth = 0; };
    Group root;
    root.faces.resize(reference.triangles.size());
    std::iota(root.faces.begin(), root.faces.end(), std::uint32_t{0});
    std::vector<Group> stack;
    stack.push_back(std::move(root));
    const auto appendFallback = [&](const MeshModel& patch)
    {
        const std::uint32_t base = static_cast<std::uint32_t>(exact_fallback.vertices.size());
        exact_fallback.vertices.insert(exact_fallback.vertices.end(),patch.vertices.begin(),patch.vertices.end());
        for (const TriangleIndices face : patch.triangles)
            exact_fallback.triangles.push_back({base+face[0],base+face[1],base+face[2]});
    };
    while (!stack.empty())
    {
        Group group = std::move(stack.back());
        stack.pop_back();
        MeshModel hull;
        bool hull_built = false;
        if (group.faces.size() >= 4)
        {
            try
            {
                hull = makeConvexHull(reference, group.faces);
                hull_built = !hull.triangles.empty() && hull.triangles.size() < group.faces.size();
            }
            catch (const std::exception&) {}
        }
        if (hull_built)
        {
            const auto certificate = certifyWithBvh(
                hull, reference_bvh, maximum_distance, maximum_depth, numerical_tolerance);
            if (certificate.passed)
            {
                ConvexPiece piece=makeConvexPiece(std::move(hull),group.faces);
                Bounds piece_bounds;
                for (const Position3 point : piece.mesh.vertices) include(piece_bounds,toVec(point));
                const double piece_scale=norm(piece_bounds.upper-piece_bounds.lower);
                const double coverage_tolerance=std::max(1.0e-9,piece_scale*1.0e-9);
                bool covers=true;
                for (const std::uint32_t face_id : group.faces)
                    for (const std::uint32_t vertex : reference.triangles[face_id])
                        for (std::size_t plane=0;plane<piece.normals.size();++plane)
                            if (dot(piece.normals[plane],toVec(reference.vertices[vertex])) >
                                piece.offsets[plane]+coverage_tolerance)
                            {
                                covers=false;
                                break;
                            }
                if (covers)
                {
                    accepted_pieces.push_back(std::move(piece));
                    continue;
                }
            }
        }
        if (group.faces.size() <= 8 || group.depth >= 32)
        {
            MeshModel fallback;
            fallback.name = "exact_face_group_fallback";
            std::unordered_map<std::uint32_t,std::uint32_t> vertex_map;
            fallback.triangles.reserve(group.faces.size());
            for (const std::uint32_t face_id : group.faces)
            {
                TriangleIndices face{};
                for (int local = 0; local < 3; ++local)
                {
                    const std::uint32_t source_vertex = reference.triangles[face_id][local];
                    const auto [iterator, inserted] = vertex_map.try_emplace(
                        source_vertex, static_cast<std::uint32_t>(fallback.vertices.size()));
                    if (inserted) fallback.vertices.push_back(reference.vertices[source_vertex]);
                    face[local] = iterator->second;
                }
                fallback.triangles.push_back(face);
            }
            appendFallback(fallback);
            continue;
        }
        Bounds centers;
        std::vector<std::pair<double,std::uint32_t>> ordered;
        ordered.reserve(group.faces.size());
        for (const std::uint32_t face_id : group.faces)
        {
            const auto face = reference.triangles[face_id];
            const Vec3 center = (toVec(reference.vertices[face[0]]) +
                                 toVec(reference.vertices[face[1]]) +
                                 toVec(reference.vertices[face[2]])) * (1.0/3.0);
            include(centers, center);
        }
        const Vec3 extent = centers.upper - centers.lower;
        const int axis = extent.y > extent.x ? (extent.z > extent.y ? 2 : 1)
                                             : (extent.z > extent.x ? 2 : 0);
        const auto coordinate = [axis](const Vec3 point)
        { return axis == 0 ? point.x : axis == 1 ? point.y : point.z; };
        for (const std::uint32_t face_id : group.faces)
        {
            const auto face = reference.triangles[face_id];
            const Vec3 center = (toVec(reference.vertices[face[0]]) +
                                 toVec(reference.vertices[face[1]]) +
                                 toVec(reference.vertices[face[2]])) * (1.0/3.0);
            ordered.emplace_back(coordinate(center), face_id);
        }
        const std::size_t middle = ordered.size()/2;
        std::nth_element(ordered.begin(), ordered.begin()+middle, ordered.end());
        Group first, second;
        first.depth = second.depth = group.depth + 1;
        first.faces.reserve(middle);
        second.faces.reserve(ordered.size()-middle);
        for (std::size_t i = 0; i < ordered.size(); ++i)
            (i < middle ? first.faces : second.faces).push_back(ordered[i].second);
        stack.push_back(std::move(second));
        stack.push_back(std::move(first));
    }
    Bounds reference_bounds;
    for (const Position3 point : reference.vertices) include(reference_bounds,toVec(point));
    const double scale = norm(reference_bounds.upper-reference_bounds.lower);
    // Every responsibility triangle is either certified inside one accepted
    // convex piece above or copied exactly. Removing faces inside another
    // accepted convex solid does not change that union coverage.
    return convexUnionOuterSurface(accepted_pieces,exact_fallback,scale);
}

Vec3 normalized(const Vec3 value)
{
    const double length = norm(value);
    if (!(length > 0.0)) throw std::invalid_argument("zero polytope direction");
    return value * (1.0 / length);
}

std::vector<Vec3> kdopDirections(const int family)
{
    std::vector<Vec3> result{{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
    const auto addSigned = [&](const Vec3 base)
    {
        for (int sx : {-1, 1}) for (int sy : {-1, 1}) for (int sz : {-1, 1})
        {
            Vec3 direction{base.x*sx, base.y*sy, base.z*sz};
            if (normSquared(direction) == 0.0) continue;
            direction = normalized(direction);
            if (std::none_of(result.begin(), result.end(), [&](const Vec3 existing)
                { return normSquared(existing - direction) < 1.0e-20; }))
                result.push_back(direction);
        }
    };
    if (family >= 18)
    {
        addSigned({1,1,0});
        addSigned({1,0,1});
        addSigned({0,1,1});
    }
    if (family == 14 || family >= 26) addSigned({1,1,1});
    return result;
}

MeshModel makeSupportPolytope(const MeshModel& reference, const int family)
{
    const std::vector<Vec3> directions = kdopDirections(family);
    std::vector<double> supports(directions.size(), -std::numeric_limits<double>::infinity());
    Bounds bounds;
    for (const Position3 point : reference.vertices)
    {
        const Vec3 p = toVec(point);
        include(bounds, p);
        for (std::size_t i = 0; i < directions.size(); ++i)
            supports[i] = std::max(supports[i], dot(directions[i], p));
    }
    const double scale = norm(bounds.upper - bounds.lower);
    const double tolerance = std::max(1.0e-10, scale * 1.0e-10);
    MeshModel result;
    result.name = std::to_string(family) + "_dop_support_polytope";
    for (std::size_t i = 0; i < directions.size(); ++i)
        for (std::size_t j = i + 1; j < directions.size(); ++j)
            for (std::size_t k = j + 1; k < directions.size(); ++k)
            {
                const Vec3 jk = cross(directions[j], directions[k]);
                const double determinant = dot(directions[i], jk);
                if (std::abs(determinant) <= 1.0e-12) continue;
                const Vec3 point = (jk * supports[i] +
                    cross(directions[k], directions[i]) * supports[j] +
                    cross(directions[i], directions[j]) * supports[k]) * (1.0 / determinant);
                bool inside = true;
                for (std::size_t plane = 0; plane < directions.size(); ++plane)
                    inside &= dot(directions[plane], point) <= supports[plane] + tolerance;
                if (!inside) continue;
                if (std::none_of(result.vertices.begin(), result.vertices.end(), [&](const Position3 existing)
                    { return normSquared(toVec(existing) - point) <= tolerance*tolerance; }))
                    result.vertices.push_back(toPosition(point));
            }
    for (std::size_t plane = 0; plane < directions.size(); ++plane)
    {
        std::vector<std::uint32_t> face_vertices;
        for (std::uint32_t vertex = 0; vertex < result.vertices.size(); ++vertex)
            if (std::abs(dot(directions[plane], toVec(result.vertices[vertex])) - supports[plane]) <=
                tolerance * 4.0)
                face_vertices.push_back(vertex);
        if (face_vertices.size() < 3) continue;
        Vec3 basis_u = std::abs(directions[plane].x) < 0.8
            ? normalized(cross(directions[plane], {1,0,0}))
            : normalized(cross(directions[plane], {0,1,0}));
        const Vec3 basis_v = cross(directions[plane], basis_u);
        Vec3 center{};
        for (const auto vertex : face_vertices) center = center + toVec(result.vertices[vertex]);
        center = center * (1.0 / face_vertices.size());
        std::sort(face_vertices.begin(), face_vertices.end(), [&](const std::uint32_t first,
                                                                  const std::uint32_t second)
        {
            const Vec3 a = toVec(result.vertices[first]) - center;
            const Vec3 b = toVec(result.vertices[second]) - center;
            return std::atan2(dot(a,basis_v), dot(a,basis_u)) <
                   std::atan2(dot(b,basis_v), dot(b,basis_u));
        });
        for (std::size_t i = 1; i + 1 < face_vertices.size(); ++i)
            result.triangles.push_back({face_vertices[0], face_vertices[i], face_vertices[i+1]});
    }
    if (result.triangles.empty()) throw std::runtime_error("support polytope construction failed");
    return result;
}

struct EdgeKey
{
    std::uint32_t first = 0;
    std::uint32_t second = 0;
    auto operator<=>(const EdgeKey&) const = default;
};

struct EdgeHash
{
    std::size_t operator()(const EdgeKey key) const
    {
        return (static_cast<std::size_t>(key.first) << 32) ^ key.second;
    }
};

MeshModel convexifyCoplanarComponents(
    const MeshModel& input, const MeshModel& reference,
    const double maximum_distance, const std::uint32_t maximum_depth,
    const double numerical_tolerance)
{
    if (input.triangles.empty()) return input;
    Bounds bounds;
    for (const Position3 point : input.vertices) include(bounds, toVec(point));
    const double scale = norm(bounds.upper - bounds.lower);
    const double planar_tolerance = std::max(1.0e-10, scale * 1.0e-11);
    std::vector<std::uint32_t> parent(input.triangles.size());
    std::iota(parent.begin(), parent.end(), std::uint32_t{0});
    const auto root = [&](std::uint32_t value)
    {
        while (parent[value] != value)
        {
            parent[value] = parent[parent[value]];
            value = parent[value];
        }
        return value;
    };
    const auto unite = [&](std::uint32_t first, std::uint32_t second)
    {
        first = root(first);
        second = root(second);
        if (first != second) parent[second] = first;
    };
    struct Plane { Vec3 normal; double offset = 0.0; bool valid = false; };
    std::vector<Plane> planes(input.triangles.size());
    for (std::size_t face_id = 0; face_id < input.triangles.size(); ++face_id)
    {
        const auto face = input.triangles[face_id];
        const Vec3 a = toVec(input.vertices[face[0]]);
        const Vec3 normal_value = cross(toVec(input.vertices[face[1]]) - a,
                                        toVec(input.vertices[face[2]]) - a);
        const double length = norm(normal_value);
        if (length > planar_tolerance * planar_tolerance)
        {
            planes[face_id] = {normal_value * (1.0 / length), 0.0, true};
            planes[face_id].offset = dot(planes[face_id].normal, a);
        }
    }
    logMonitor("planar_convexification_planes_built");
    std::unordered_map<EdgeKey, std::uint32_t, EdgeHash> edges;
    for (std::uint32_t face_id = 0; face_id < input.triangles.size(); ++face_id)
    {
        const auto face = input.triangles[face_id];
        for (int local = 0; local < 3; ++local)
        {
            const EdgeKey edge{std::min(face[local], face[(local+1)%3]),
                               std::max(face[local], face[(local+1)%3])};
            const auto [iterator, inserted] = edges.emplace(edge, face_id);
            if (inserted) continue;
            const std::uint32_t other = iterator->second;
            if (!planes[face_id].valid || !planes[other].valid) continue;
            const double alignment = dot(planes[face_id].normal, planes[other].normal);
            if (alignment < 1.0 - 1.0e-12) continue;
            const Position3 point = input.vertices[input.triangles[other][0]];
            if (std::abs(dot(planes[face_id].normal, toVec(point)) - planes[face_id].offset) >
                planar_tolerance) continue;
            unite(face_id, other);
        }
    }
    logMonitor("planar_convexification_edges_built");
    std::unordered_map<std::uint32_t, std::vector<std::uint32_t>> groups;
    for (std::uint32_t face = 0; face < input.triangles.size(); ++face)
        groups[root(face)].push_back(face);
    logMonitor("planar_convexification_groups_built");

    std::vector<std::uint8_t> replaced(input.triangles.size(), 0);
    const TriangleBvh reference_bvh(reference);
    logMonitor("planar_convexification_reference_bvh_built");
    MeshModel result;
    result.name = "certified_planar_component_convexification";
    result.vertices = input.vertices;
    for (const auto& [representative, faces] : groups)
    {
        if (faces.size() <= 1 || !planes[representative].valid) continue;
        const Vec3 normal = planes[representative].normal;
        const Vec3 basis_u = std::abs(normal.x) < 0.8
            ? normalized(cross(normal, {1,0,0}))
            : normalized(cross(normal, {0,1,0}));
        const Vec3 basis_v = cross(normal, basis_u);
        struct Point2 { double x,y; std::uint32_t vertex; };
        std::vector<Point2> points;
        for (const std::uint32_t face_id : faces)
            for (const std::uint32_t vertex : input.triangles[face_id])
                points.push_back({dot(toVec(input.vertices[vertex]), basis_u),
                                  dot(toVec(input.vertices[vertex]), basis_v), vertex});
        std::sort(points.begin(), points.end(), [](const Point2& a, const Point2& b)
        { return std::tie(a.x,a.y,a.vertex) < std::tie(b.x,b.y,b.vertex); });
        points.erase(std::unique(points.begin(), points.end(), [&](const Point2& a, const Point2& b)
        { return std::abs(a.x-b.x) <= planar_tolerance && std::abs(a.y-b.y) <= planar_tolerance; }),
                     points.end());
        if (points.size() < 3) continue;
        const auto turn = [](const Point2& a, const Point2& b, const Point2& c)
        { return (b.x-a.x)*(c.y-a.y) - (b.y-a.y)*(c.x-a.x); };
        std::vector<Point2> hull(points.size()*2);
        std::size_t count = 0;
        for (const Point2 point : points)
        {
            while (count >= 2 && turn(hull[count-2], hull[count-1], point) <= 0.0) --count;
            hull[count++] = point;
        }
        const std::size_t lower_count = count;
        for (auto iterator = points.rbegin() + 1; iterator != points.rend(); ++iterator)
        {
            while (count > lower_count && turn(hull[count-2], hull[count-1], *iterator) <= 0.0) --count;
            hull[count++] = *iterator;
        }
        if (count > 1) --count;
        hull.resize(count);
        if (hull.size() < 3 || hull.size() - 2 >= faces.size()) continue;
        MeshModel patch;
        patch.name = "planar_component_candidate";
        for (const Point2 point : hull)
        {
            const Vec3 projected = basis_u * point.x + basis_v * point.y +
                normal * planes[representative].offset;
            patch.vertices.push_back(toPosition(projected));
        }
        for (std::size_t index = 1; index + 1 < patch.vertices.size(); ++index)
            patch.triangles.push_back({0, static_cast<std::uint32_t>(index),
                                      static_cast<std::uint32_t>(index+1)});
        const DirectedHausdorffCertificate certificate = certifyWithBvh(
            patch, reference_bvh, maximum_distance, maximum_depth, numerical_tolerance);
        if (!certificate.passed) continue;
        const std::uint32_t base = static_cast<std::uint32_t>(result.vertices.size());
        result.vertices.insert(result.vertices.end(), patch.vertices.begin(), patch.vertices.end());
        for (const TriangleIndices face : patch.triangles)
            result.triangles.push_back({base+face[0], base+face[1], base+face[2]});
        for (const std::uint32_t face : faces) replaced[face] = 1;
    }
    for (std::size_t face = 0; face < input.triangles.size(); ++face)
        if (!replaced[face]) result.triangles.push_back(input.triangles[face]);
    logMonitor("planar_convexification_patches_built");
    logMonitor("planar_convexification_final_union_started");
    MeshModel canonical = canonicalizeCoplanarTriangleSoup(result);
    logMonitor("planar_convexification_final_union_finished");
    return canonical;
}

void writeObj(const std::filesystem::path& path, const MeshModel& mesh)
{
    std::ofstream stream(path);
    if (!stream) throw std::runtime_error("failed to create OBJ: " + path.string());
    stream << std::setprecision(17) << "o " << mesh.name << '\n';
    for (const Position3 p : mesh.vertices) stream << "v " << p.x << ' ' << p.y << ' ' << p.z << '\n';
    for (const TriangleIndices f : mesh.triangles)
        stream << "f " << f[0] + 1 << ' ' << f[1] + 1 << ' ' << f[2] + 1 << '\n';
}

template<class T>
void readExact(std::ifstream& stream, T* data, const std::size_t count, const char* label)
{
    stream.read(reinterpret_cast<char*>(data), static_cast<std::streamsize>(sizeof(T) * count));
    if (!stream) throw std::runtime_error(std::string("truncated PQSSHED1 ") + label);
}

} // namespace

struct DirectedHausdorffCertifier::Implementation
{
    explicit Implementation(const MeshModel& reference) : bvh(reference) {}
    TriangleBvh bvh;
};

DirectedHausdorffCertifier::DirectedHausdorffCertifier(const MeshModel& reference)
    : implementation_(std::make_unique<Implementation>(reference))
{
}

DirectedHausdorffCertifier::~DirectedHausdorffCertifier() = default;
DirectedHausdorffCertifier::DirectedHausdorffCertifier(
    DirectedHausdorffCertifier&&) noexcept = default;
DirectedHausdorffCertifier& DirectedHausdorffCertifier::operator=(
    DirectedHausdorffCertifier&&) noexcept = default;

DirectedHausdorffCertificate DirectedHausdorffCertifier::certify(
    const MeshModel& proxy, const double maximum_distance,
    const std::uint32_t maximum_depth, const double numerical_tolerance) const
{
    if (!std::isfinite(maximum_distance) || maximum_distance < 0.0)
        throw std::invalid_argument(
            "maximum directed Hausdorff distance must be finite and non-negative");
    return certifyWithBvh(proxy, implementation_->bvh, maximum_distance,
                          maximum_depth, numerical_tolerance);
}

OrientedSurfaceMesh readAnalysisHalfedgeMesh(const std::filesystem::path& path)
{
    std::ifstream stream(path, std::ios::binary);
    if (!stream) throw std::runtime_error("failed to open halfedge mesh: " + path.string());
    std::array<char, 8> magic{};
    readExact(stream, magic.data(), magic.size(), "magic");
    if (std::memcmp(magic.data(), "PQSSHED1", 8) != 0)
        throw std::runtime_error("invalid PQSSHED1 magic");
    std::array<std::uint32_t, 5> counts{};
    readExact(stream, counts.data(), counts.size(), "counts");
    OrientedSurfaceMesh mesh;
    mesh.geometry.name = path.stem().string();
    mesh.geometry.vertices.resize(counts[0]);
    mesh.geometry.triangles.resize(counts[1]);
    mesh.halfedges.resize(counts[2]);
    mesh.vertex_halfedges.resize(counts[3]);
    mesh.face_halfedges.resize(counts[4]);
    for (Position3& p : mesh.geometry.vertices)
    {
        std::array<double, 3> value{};
        readExact(stream, value.data(), value.size(), "vertices");
        p = {value[0], value[1], value[2]};
    }
    for (auto& face : mesh.geometry.triangles) readExact(stream, face.data(), face.size(), "faces");
    for (auto& edge : mesh.halfedges)
    {
        std::array<std::uint32_t, 4> value{};
        readExact(stream, value.data(), value.size(), "halfedges");
        edge = {value[0], value[1], value[2], value[3]};
    }
    readExact(stream, mesh.vertex_halfedges.data(), mesh.vertex_halfedges.size(), "vertex halfedges");
    readExact(stream, mesh.face_halfedges.data(), mesh.face_halfedges.size(), "face halfedges");
    if (stream.peek() != std::ifstream::traits_type::eof())
        throw std::runtime_error("PQSSHED1 contains trailing bytes");
    return mesh;
}

DirectedHausdorffCertificate certifyDirectedHausdorff(
    const MeshModel& proxy, const MeshModel& reference, const double maximum_distance,
    const std::uint32_t maximum_depth, const double numerical_tolerance)
{
    if (!std::isfinite(maximum_distance) || maximum_distance < 0.0)
        throw std::invalid_argument("maximum directed Hausdorff distance must be finite and non-negative");
    DirectedHausdorffCertificate result;
    result.passed = true;
    const auto samePoint = [&](const Position3 first, const Position3 second)
    {
        return std::abs(first.x - second.x) <= numerical_tolerance &&
               std::abs(first.y - second.y) <= numerical_tolerance &&
               std::abs(first.z - second.z) <= numerical_tolerance;
    };
    const bool identical = proxy.vertices.size() == reference.vertices.size() &&
        proxy.triangles == reference.triangles &&
        std::equal(proxy.vertices.begin(), proxy.vertices.end(), reference.vertices.begin(),
            [&](const Position3 first, const Position3 second)
            { return samePoint(first, second); });
    if (identical)
    {
        result.upper_bound = 0.0;
        return result;
    }
    DirectedHausdorffCertifier certifier(reference);
    return certifier.certify(proxy, maximum_distance, maximum_depth,
                             numerical_tolerance);
}

MeshModel convexifyCertifiedPlanarComponents(
    const MeshModel& proxy, const MeshModel& reference,
    const double maximum_distance, const std::uint32_t maximum_depth,
    const double numerical_tolerance)
{
    return convexifyCoplanarComponents(
        canonicalizeCoplanarTriangleSoup(proxy), reference,
        maximum_distance, maximum_depth, numerical_tolerance);
}

HausdorffSimplificationStats simplifyPhase1Halfedge(
    const std::filesystem::path& phase1_halfedge, const std::filesystem::path& source_obj,
    const std::filesystem::path& output_directory, const HausdorffSimplificationOptions& options)
{
    const auto started = std::chrono::steady_clock::now();
    if (!std::isfinite(options.maximum_directed_hausdorff) ||
        options.maximum_directed_hausdorff < 0.0)
        throw std::invalid_argument("maximum directed Hausdorff distance must be finite and non-negative");
    std::filesystem::create_directories(output_directory);
    const OrientedSurfaceMesh phase1 = readAnalysisHalfedgeMesh(phase1_halfedge);
    logMonitor("phase1_loaded");
    const MeshModel& reference = phase1.geometry;
    HausdorffSimplificationStats stats;
    stats.phase1_vertices = reference.vertices.size();
    stats.phase1_triangles = reference.triangles.size();
    stats.maximum_directed_hausdorff = options.maximum_directed_hausdorff;

    MeshModel selected = reference;
    selected.name = "exact_phase1_fallback";
    stats.selected_candidate = selected.name;
    HausdorffCandidateStats fallback;
    fallback.name = selected.name;
    fallback.triangles = selected.triangles.size();
    fallback.conservative_coverage = true;
    fallback.hausdorff.passed = true;
    fallback.hausdorff.upper_bound = 0.0;
    fallback.selected = true;
    stats.candidates.push_back(fallback);

    const auto recordGenerationFailure = [&](const std::string& name,
                                             const std::string& reason)
    {
        HausdorffCandidateStats record;
        record.name = name;
        record.hausdorff.failure_reason = "candidate generation failed: " + reason;
        stats.candidates.push_back(std::move(record));
        std::cerr << "monitor: stage=candidate_generation_failed candidate="
                  << name << " reason=" << reason << '\n';
    };

    const auto consider = [&](MeshModel candidate, const std::string& name,
                              const bool conservative_coverage,
                              const bool exact_surface_equivalence = false)
    {
        logMonitor("candidate_evaluation_started", name);
        HausdorffCandidateStats record;
        record.name = name;
        record.triangles = candidate.triangles.size();
        record.conservative_coverage = conservative_coverage;
        if (conservative_coverage && candidate.triangles.size() < selected.triangles.size() &&
            exact_surface_equivalence)
        {
            record.hausdorff.passed = true;
            record.hausdorff.upper_bound = 0.0;
        }
        else if (conservative_coverage && candidate.triangles.size() < selected.triangles.size())
        {
            logMonitor("candidate_certificate_started", name);
            record.hausdorff = certifyDirectedHausdorff(
                candidate, reference, options.maximum_directed_hausdorff,
                options.maximum_certificate_depth, options.numerical_tolerance);
            std::cerr << "monitor: stage=candidate_certificate_result candidate=" << name
                      << " passed=" << (record.hausdorff.passed ? "true" : "false")
                      << " lower_bound=" << record.hausdorff.lower_bound
                      << " upper_bound=" << record.hausdorff.upper_bound
                      << " reason=" << record.hausdorff.failure_reason << '\n';
            logMonitor("candidate_certificate_finished", name);
        }
        else
        {
            record.hausdorff.passed = false;
            record.hausdorff.failure_reason = conservative_coverage
                ? "candidate does not reduce final triangle count"
                : "candidate lacks a conservative coverage certificate";
        }
        if (record.hausdorff.passed)
        {
            for (auto& existing : stats.candidates) existing.selected = false;
            record.selected = true;
            selected = std::move(candidate);
            selected.name = name;
            stats.selected_candidate = name;
        }
        stats.candidates.push_back(std::move(record));
        logMonitor("candidate_evaluation_finished", name);
    };

    if (options.enable_exact_coplanar_union)
    {
        try
        {
            logMonitor("candidate_generation_started", "exact_coplanar_union");
            MeshModel coplanar = canonicalizeCoplanarTriangleSoup(reference);
            logMonitor("candidate_generation_finished", "exact_coplanar_union");
            if (options.enable_planar_component_convexification)
            {
                try
                {
                    logMonitor("candidate_generation_started", "planar_component_convexification");
                    MeshModel convexified = convexifyCoplanarComponents(
                        coplanar, reference, options.maximum_directed_hausdorff,
                        options.maximum_certificate_depth, options.numerical_tolerance);
                    logMonitor("candidate_generation_finished", "planar_component_convexification");
                    consider(std::move(convexified), "planar_component_convexification", true);
                }
                catch (const std::exception& error)
                {
                    recordGenerationFailure("planar_component_convexification", error.what());
                }
            }
            consider(std::move(coplanar), "exact_coplanar_union", true, true);
        }
        catch (const std::exception& error)
        {
            recordGenerationFailure("exact_coplanar_union", error.what());
            recordGenerationFailure("planar_component_convexification",
                "exact coplanar input unavailable");
        }
    }
    if (options.enable_convex_hull)
    {
        const std::string name = "convex_hull_of_all_phase1_vertices";
        try
        {
            logMonitor("candidate_generation_started", name);
            MeshModel candidate = makeConvexHull(reference);
            logMonitor("candidate_generation_finished", name);
            consider(std::move(candidate), name, true);
        }
        catch (const std::exception& error)
        {
            recordGenerationFailure(name, error.what());
        }
    }
    if (options.enable_adaptive_convex_cover)
    {
        const std::string name = "adaptive_certified_convex_cover";
        try
        {
            logMonitor("candidate_generation_started", name);
            MeshModel candidate = adaptiveConvexCover(
                reference, options.maximum_directed_hausdorff,
                options.maximum_certificate_depth, options.numerical_tolerance);
            logMonitor("candidate_generation_finished", name);
            consider(std::move(candidate), name, true);
        }
        catch (const std::exception& error)
        {
            recordGenerationFailure(name, error.what());
        }
    }
    if (options.enable_bounded_box_cover)
    {
        for (const int retreat_ten_thousandths : {10000})
        {
            const std::string name = "closed_4_obb_q50_clip11_final_strict4_margin5e5_cert28";
            try
            {
                logMonitor("candidate_generation_started", name);
                MeshModel candidate = boundedBoxCover(
                    reference, 4, 0.5, 15,
                    options.maximum_directed_hausdorff, 0,
                    0.0, false, false,
                    retreat_ten_thousandths / 10000.0, 0, false,
                    options.maximum_directed_hausdorff * 5.0e-5, 28, 4);
                if (candidate.triangles.size() > 100)
                    throw std::runtime_error("candidate exceeds 100 triangles");
                AnalysisHalfedgeStats topology;
                const OrientedSurfaceMesh candidate_halfedge =
                    buildAnalysisHalfedgeMesh(candidate, &topology);
                if (topology.boundary_halfedges != 0 ||
                    topology.nonmanifold_edge_groups != 0 ||
                    topology.inconsistent_orientation_edges != 0)
                    throw std::runtime_error(
                        "candidate halfedge invalid: boundary=" +
                        std::to_string(topology.boundary_halfedges) +
                        " nonmanifold=" + std::to_string(topology.nonmanifold_edge_groups) +
                        " orientation=" +
                        std::to_string(topology.inconsistent_orientation_edges) +
                        " components=" + std::to_string(topology.face_components));
                (void)validateClosedHalfedgeTopology(candidate_halfedge);
                logMonitor("candidate_generation_finished", name);
                logMonitor("candidate_evaluation_started", name);
                HausdorffCandidateStats record;
                record.name = name;
                record.triangles = candidate.triangles.size();
                record.conservative_coverage = true;
                logMonitor("candidate_certificate_started", name);
                record.hausdorff = certifyDirectedHausdorff(
                    candidate, reference, options.maximum_directed_hausdorff,
                    28, options.numerical_tolerance);
                std::cerr << std::setprecision(17)
                          << "monitor: stage=candidate_certificate_result candidate=" << name
                          << " passed=" << (record.hausdorff.passed ? "true" : "false")
                          << " lower_bound=" << record.hausdorff.lower_bound
                          << " upper_bound=" << record.hausdorff.upper_bound
                          << " reason=" << record.hausdorff.failure_reason << '\n';
                logMonitor("candidate_certificate_finished", name);
                if (record.hausdorff.passed &&
                    record.hausdorff.lower_bound <= options.maximum_directed_hausdorff &&
                    record.hausdorff.upper_bound <= options.maximum_directed_hausdorff &&
                    candidate.triangles.size() < selected.triangles.size())
                {
                    for (auto& existing : stats.candidates) existing.selected = false;
                    record.selected = true;
                    selected = std::move(candidate);
                    selected.name = name;
                    stats.selected_candidate = name;
                }
                stats.candidates.push_back(std::move(record));
                logMonitor("candidate_evaluation_finished", name);
            }
            catch (const std::exception& error)
            {
                recordGenerationFailure(name, error.what());
            }
        }
    }
    if (options.enable_discrete_orientation_polytopes)
        for (const int family : {14, 18, 26})
        {
            const std::string name =
                std::to_string(family) + "_dop_support_polytope";
            try
            {
                logMonitor("candidate_generation_started", name);
                MeshModel candidate = makeSupportPolytope(reference, family);
                logMonitor("candidate_generation_finished", name);
                consider(std::move(candidate), name, true);
            }
            catch (const std::exception& error)
            {
                recordGenerationFailure(name, error.what());
            }
        }
    if (options.enable_axis_aligned_box)
    {
        const std::string name = "axis_aligned_enclosing_box";
        try
        {
            logMonitor("candidate_generation_started", name);
            MeshModel candidate = makeBox(reference);
            logMonitor("candidate_generation_finished", name);
            consider(std::move(candidate), name, true);
        }
        catch (const std::exception& error)
        {
            recordGenerationFailure(name, error.what());
        }
    }

    stats.final_vertices = selected.vertices.size();
    stats.final_triangles = selected.triangles.size();
    stats.elapsed_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    const bool has_phase2 = stats.selected_candidate != "exact_phase1_fallback";
    if (has_phase2)
    {
        AnalysisHalfedgeStats selected_topology;
        const OrientedSurfaceMesh selected_halfedge =
            buildAnalysisHalfedgeMesh(selected, &selected_topology);
        (void)validateClosedHalfedgeTopology(selected_halfedge);
        writeAnalysisHalfedgeMesh(output_directory / "phase2_halfedge.bin", selected_halfedge);
    }
    writeObj(output_directory / "proxy.obj", selected);
    std::filesystem::copy_file(phase1_halfedge, output_directory / "phase1_halfedge.bin",
                               std::filesystem::copy_options::overwrite_existing);
    std::filesystem::copy_file(source_obj, output_directory / "source.obj",
                               std::filesystem::copy_options::overwrite_existing);

    std::ofstream model(output_directory / "model.json");
    const auto jsonNumber = [](const double value)
    {
        if (!std::isfinite(value)) return std::string("null");
        std::ostringstream stream;
        stream << std::setprecision(17) << value;
        return stream.str();
    };
    model << std::setprecision(17)
          << "{\n  \"stats\":{\"source_triangles\":" << stats.phase1_triangles
          << ",\"proxy_triangles\":" << stats.final_triangles
          << ",\"primitive_count\":1,\"primitive_types\":{\"polygon\":1}"
          << ",\"timings_seconds\":{\"total\":" << stats.elapsed_seconds << "}"
          << ",\"simplification_error\":{\"direction\":\"proxy_to_phase1\""
          << ",\"method\":\"triangle_bvh_1_lipschitz_adaptive_upper_bound\""
          << ",\"maximum_is_certified_upper_bound\":true"
          << ",\"maximum_distance_limit\":" << options.maximum_directed_hausdorff;
    const auto selected_record = std::find_if(stats.candidates.begin(), stats.candidates.end(),
        [](const auto& candidate) { return candidate.selected; });
    model << ",\"maximum_distance\":" << selected_record->hausdorff.upper_bound
          << ",\"observed_lower_bound\":" << selected_record->hausdorff.lower_bound
          << ",\"reference_queries\":" << selected_record->hausdorff.reference_queries
          << ",\"subdivision_nodes\":" << selected_record->hausdorff.subdivision_nodes << '}'
          << ",\"optimization\":{\"objective\":\"minimum_final_obj_triangle_count\""
          << ",\"status\":\"best-known feasible over generated certified candidates\""
          << ",\"selected_candidate\":\"" << stats.selected_candidate << "\""
          << ",\"candidates\":[";
    for (std::size_t i = 0; i < stats.candidates.size(); ++i)
    {
        if (i) model << ',';
        const auto& c = stats.candidates[i];
        model << "{\"name\":\"" << c.name << "\",\"triangles\":" << c.triangles
              << ",\"coverage_certified\":" << (c.conservative_coverage ? "true" : "false")
              << ",\"hausdorff_certified\":" << (c.hausdorff.passed ? "true" : "false")
              << ",\"upper_bound\":" << jsonNumber(c.hausdorff.upper_bound)
              << ",\"selected\":" << (c.selected ? "true" : "false") << '}';
    }
    model << "]}},\n  \"source\":\"source.obj\",\n  \"phase1_halfedge\":\"phase1_halfedge.bin\","
          << (has_phase2 ? "\n  \"phase2_halfedge\":\"phase2_halfedge.bin\"," : "")
          << "\n  \"phase4_triangulated\":\"proxy.obj\",\n  \"proxy\":\"proxy.obj\","
          << "\n  \"proxy_components\":[{\"id\":0,\"type\":\"polygon\",\"triangulated_face_count\":"
          << stats.final_triangles << "}],\n  \"viewer_stages\":[\"source\",\"phase1\",\"phase4\",\"split\"]\n}\n";
    std::ofstream manifest(output_directory / "viewer_manifest.json");
    const std::string model_id = options.model_id.empty()
        ? source_obj.stem().string() : options.model_id;
    manifest << "{\"algorithm\":\"CertifiedDirectedHausdorffSimplifierV1\""
             << ",\"complete\":" << (has_phase2 ? "true" : "false")
             << ",\"model_count\":1,\"models\":[{\"id\":\""
             << model_id << "\",\"metadata\":\"model.json\"}]"
             << ",\"options\":{\"maximum_directed_hausdorff\":"
             << std::setprecision(17) << options.maximum_directed_hausdorff
             << ",\"direction\":\"proxy_to_phase1\"}}\n";
    return stats;
}

} // namespace pqss_proxy_mesh
