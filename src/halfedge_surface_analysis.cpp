#include "pqss_proxy_mesh/halfedge_surface_analysis.hpp"
#include "pqss_proxy_mesh/hausdorff_simplifier.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <deque>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

namespace pqss_proxy_mesh
{
namespace
{

struct Vec3
{
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct Vec2
{
    double x = 0.0;
    double y = 0.0;
    std::uint32_t source_vertex = invalid_surface_index;
};

Vec3 toVec(const Position3 p) { return {p.x, p.y, p.z}; }
Position3 toPosition(const Vec3 p) { return {p.x, p.y, p.z}; }
Vec3 operator+(const Vec3 a, const Vec3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
Vec3 operator-(const Vec3 a, const Vec3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
Vec3 operator*(const Vec3 a, const double s) { return {a.x * s, a.y * s, a.z * s}; }
double dot(const Vec3 a, const Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
Vec3 cross(const Vec3 a, const Vec3 b)
{
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x};
}
double norm(const Vec3 a) { return std::sqrt(dot(a, a)); }
Vec3 normalized(const Vec3 a)
{
    const double length = norm(a);
    if (!(length > 0.0)) throw std::invalid_argument("zero-length direction");
    return a * (1.0 / length);
}
Vec2 operator-(const Vec2 a, const Vec2 b) { return {a.x - b.x, a.y - b.y, 0}; }
double cross2(const Vec2 a, const Vec2 b) { return a.x * b.y - a.y * b.x; }

Vec3 faceRawNormal(const OrientedSurfaceMesh& mesh, const std::uint32_t face)
{
    const auto triangle = mesh.geometry.triangles[face];
    const Vec3 a = toVec(mesh.geometry.vertices[triangle[0]]);
    const Vec3 b = toVec(mesh.geometry.vertices[triangle[1]]);
    const Vec3 c = toVec(mesh.geometry.vertices[triangle[2]]);
    const Vec3 raw = cross(b - a, c - a);
    if (!(norm(raw) > 0.0))
        throw std::runtime_error("halfedge mesh has a degenerate face");
    return raw;
}

struct PlaneFrame
{
    Vec3 normal;
    Vec3 tangent;
    Vec3 bitangent;
    double offset = 0.0;
};

PlaneFrame makePlaneFrame(const Vec3 normal, const double offset)
{
    const Vec3 helper = std::abs(normal.x) < 0.8 ? Vec3{1.0, 0.0, 0.0}
                                                  : Vec3{0.0, 1.0, 0.0};
    const Vec3 tangent = normalized(cross(helper, normal));
    return {normal, tangent, cross(normal, tangent), offset};
}

Vec2 project(const PlaneFrame& frame, const Vec3 point, const std::uint32_t source_vertex)
{
    return {dot(point, frame.tangent), dot(point, frame.bitangent), source_vertex};
}

Vec3 lift(const PlaneFrame& frame, const Vec2 point)
{
    return frame.normal * frame.offset + frame.tangent * point.x +
           frame.bitangent * point.y;
}

std::vector<Vec2> convexHull(std::vector<Vec2> points, const double tolerance)
{
    std::sort(points.begin(), points.end(), [](const Vec2 a, const Vec2 b)
    {
        if (a.x != b.x) return a.x < b.x;
        if (a.y != b.y) return a.y < b.y;
        return a.source_vertex < b.source_vertex;
    });
    points.erase(std::unique(points.begin(), points.end(), [&](const Vec2 a, const Vec2 b)
    {
        return std::abs(a.x - b.x) <= tolerance &&
               std::abs(a.y - b.y) <= tolerance;
    }), points.end());
    if (points.size() < 3) return {};

    const auto chain = [&](auto begin, auto end)
    {
        std::vector<Vec2> result;
        for (auto iterator = begin; iterator != end; ++iterator)
        {
            const Vec2 point = *iterator;
            while (result.size() >= 2 &&
                   cross2(result.back() - result[result.size() - 2],
                          point - result.back()) <= tolerance)
                result.pop_back();
            result.push_back(point);
        }
        return result;
    };
    std::vector<Vec2> hull = chain(points.begin(), points.end());
    std::vector<Vec2> upper = chain(points.rbegin(), points.rend());
    hull.pop_back();
    upper.pop_back();
    hull.insert(hull.end(), upper.begin(), upper.end());
    return hull.size() >= 3 ? hull : std::vector<Vec2>{};
}

template <typename T>
std::vector<T> mergedUnique(const std::vector<T>& first, const std::vector<T>& second)
{
    std::vector<T> result;
    result.reserve(first.size() + second.size());
    std::set_union(first.begin(), first.end(), second.begin(), second.end(),
                   std::back_inserter(result));
    return result;
}

struct MergedRegionData
{
    double certified_upper_bound = 0.0;
    std::vector<std::uint32_t> faces;
    std::vector<std::uint32_t> responsibility_vertices;
    std::vector<std::uint32_t> neighbors;
    std::vector<Position3> polygon_vertices;
};

// Leaves have no per-region containers. Inactive merged payloads are released
// immediately after their replacement parent has taken responsibility.
struct Region
{
    std::unique_ptr<MergedRegionData> merged;
    std::uint32_t replacement = invalid_surface_index;
    std::uint32_t triangle_count = 1;
};

static_assert(sizeof(Region) <= 16, "leaf region layout unexpectedly grew");

struct MergeCandidate
{
    bool orientation_consistent = false;
    std::vector<std::uint32_t> faces;
    std::vector<std::uint32_t> responsibility_vertices;
    MeshModel polygon;
};

MergeCandidate makeMergeCandidate(const OrientedSurfaceMesh& mesh,
                                   const std::uint32_t first_id, const Region& first,
                                   const std::uint32_t second_id, const Region& second,
                                   const double tolerance)
{
    MergeCandidate candidate;
    const auto appendResponsibility = [&](const std::uint32_t id, const Region& region)
    {
        if (region.merged)
        {
            candidate.faces.insert(candidate.faces.end(), region.merged->faces.begin(),
                                   region.merged->faces.end());
            candidate.responsibility_vertices.insert(
                candidate.responsibility_vertices.end(),
                region.merged->responsibility_vertices.begin(),
                region.merged->responsibility_vertices.end());
            return;
        }
        candidate.faces.push_back(id);
        const auto face = mesh.geometry.triangles[id];
        candidate.responsibility_vertices.insert(
            candidate.responsibility_vertices.end(), face.begin(), face.end());
    };
    appendResponsibility(first_id, first);
    appendResponsibility(second_id, second);
    std::sort(candidate.faces.begin(), candidate.faces.end());
    candidate.faces.erase(std::unique(candidate.faces.begin(), candidate.faces.end()),
                          candidate.faces.end());
    std::sort(candidate.responsibility_vertices.begin(),
              candidate.responsibility_vertices.end());
    candidate.responsibility_vertices.erase(
        std::unique(candidate.responsibility_vertices.begin(),
                    candidate.responsibility_vertices.end()),
        candidate.responsibility_vertices.end());

    Vec3 weighted_normal{};
    for (const std::uint32_t face : candidate.faces)
        weighted_normal = weighted_normal + faceRawNormal(mesh, face);
    if (!(norm(weighted_normal) > tolerance)) return candidate;
    const Vec3 normal = normalized(weighted_normal);
    for (const std::uint32_t face : candidate.faces)
    {
        const Vec3 raw_normal = faceRawNormal(mesh, face);
        if (dot(normal, raw_normal) <= tolerance * norm(raw_normal)) return candidate;
    }
    candidate.orientation_consistent = true;

    double support = -std::numeric_limits<double>::infinity();
    for (const std::uint32_t vertex : candidate.responsibility_vertices)
        support = std::max(support, dot(normal, toVec(mesh.geometry.vertices[vertex])));
    const PlaneFrame frame = makePlaneFrame(normal, support);
    std::vector<Vec2> projected;
    projected.reserve(candidate.responsibility_vertices.size());
    for (const std::uint32_t vertex : candidate.responsibility_vertices)
        projected.push_back(project(frame, toVec(mesh.geometry.vertices[vertex]), vertex));
    const std::vector<Vec2> hull = convexHull(std::move(projected), tolerance);
    if (hull.size() < 3)
    {
        candidate.orientation_consistent = false;
        return candidate;
    }

    candidate.polygon.name = "support_plane_region";
    candidate.polygon.vertices.reserve(hull.size());
    for (const Vec2 point : hull)
        candidate.polygon.vertices.push_back(toPosition(lift(frame, point)));
    for (std::uint32_t index = 1; index + 1 < hull.size(); ++index)
        candidate.polygon.triangles.push_back({0, index, index + 1});
    return candidate;
}

void appendRegion(const OrientedSurfaceMesh& mesh, const std::uint32_t region_id,
                   const Region& region, MeshModel& output,
                   std::map<std::array<double, 3>, std::uint32_t>& vertex_map)
{
    const auto appendTriangle = [&](const std::array<Position3, 3>& points)
    {
        std::array<std::uint32_t, 3> mapped{};
        for (int corner = 0; corner < 3; ++corner)
        {
            const Position3 point = points[corner];
            const std::array<double, 3> key{point.x, point.y, point.z};
            const auto [iterator, inserted] = vertex_map.try_emplace(
                key, static_cast<std::uint32_t>(output.vertices.size()));
            if (inserted) output.vertices.push_back(point);
            mapped[corner] = iterator->second;
        }
        output.triangles.push_back(mapped);
    };
    if (!region.merged)
    {
        const auto face = mesh.geometry.triangles[region_id];
        appendTriangle({mesh.geometry.vertices[face[0]], mesh.geometry.vertices[face[1]],
                        mesh.geometry.vertices[face[2]]});
        return;
    }
    for (std::uint32_t index = 1;
         index + 1 < region.merged->polygon_vertices.size(); ++index)
        appendTriangle({region.merged->polygon_vertices[0],
                        region.merged->polygon_vertices[index],
                        region.merged->polygon_vertices[index + 1]});
}

void writeObj(const std::filesystem::path& path, const MeshModel& mesh)
{
    std::ofstream stream(path);
    if (!stream) throw std::runtime_error("failed to create OBJ: " + path.string());
    stream << std::setprecision(17) << "o " << mesh.name << '\n';
    for (const Position3 point : mesh.vertices)
        stream << "v " << point.x << ' ' << point.y << ' ' << point.z << '\n';
    for (const auto face : mesh.triangles)
        stream << "f " << face[0] + 1 << ' ' << face[1] + 1 << ' ' << face[2] + 1 << '\n';
}

} // namespace

HalfedgeSurfaceAnalysisStats analyzeHalfedgeSurface(
    const OrientedSurfaceMesh& mesh, const std::filesystem::path& output_directory,
    const HalfedgeSurfaceAnalysisOptions& options)
{
    const auto started = std::chrono::steady_clock::now();
    auto last_progress = started;
    if (!std::isfinite(options.maximum_directed_hausdorff) ||
        options.maximum_directed_hausdorff < 0.0)
        throw std::invalid_argument("maximum directed Hausdorff distance must be non-negative");
    if (mesh.halfedges.size() != mesh.geometry.triangles.size() * 3 ||
        mesh.face_halfedges.size() != mesh.geometry.triangles.size())
        throw std::invalid_argument("halfedge mesh has inconsistent array counts");
    std::filesystem::create_directories(output_directory);

    HalfedgeSurfaceAnalysisStats stats;
    stats.phase1_vertices = mesh.geometry.vertices.size();
    stats.phase1_triangles = mesh.geometry.triangles.size();
    stats.maximum_directed_hausdorff = options.maximum_directed_hausdorff;
    const DirectedHausdorffCertifier certifier(mesh.geometry);
    std::cerr << "[stage] phase1 BVH built; initializing triangle regions\n";

    using Pair = std::pair<std::uint32_t, std::uint32_t>;
    const std::size_t face_count = mesh.geometry.triangles.size();
    if (face_count > std::numeric_limits<std::uint32_t>::max() / 2)
        throw std::invalid_argument("too many faces for compact region identifiers");
    std::vector<Region> regions;
    regions.reserve(face_count * 2);
    for (std::uint32_t face = 0; face < face_count; ++face)
    {
        Region leaf;
        leaf.replacement = face;
        regions.push_back(std::move(leaf));
    }

    // Initial leaf adjacency is CSR: one offset per face and one uint32 per
    // directed neighboring face, with no allocator state per leaf.
    std::vector<std::uint32_t> adjacency_offsets(face_count + 1, 0);
    std::deque<Pair> queue;
    const auto enqueue = [&](const std::uint32_t first, const std::uint32_t second)
    {
        if (first == second) return;
        queue.emplace_back(std::min(first, second), std::max(first, second));
    };
    for (std::size_t edge = 0; edge < mesh.halfedges.size(); ++edge)
    {
        const std::uint32_t opposite = mesh.halfedges[edge].opposite;
        if (opposite == invalid_surface_index || opposite <= edge) continue;
        if (opposite >= mesh.halfedges.size())
            throw std::invalid_argument("halfedge mesh has an invalid opposite index");
        const std::uint32_t first = mesh.halfedges[edge].face;
        const std::uint32_t second = mesh.halfedges[opposite].face;
        if (first >= face_count || second >= face_count)
            throw std::invalid_argument("halfedge mesh has an invalid face index");
        if (first == second) continue;
        ++adjacency_offsets[first + 1];
        ++adjacency_offsets[second + 1];
        enqueue(first, second);
    }
    for (std::size_t face = 1; face < adjacency_offsets.size(); ++face)
    {
        const std::uint64_t total = static_cast<std::uint64_t>(adjacency_offsets[face]) +
                                    adjacency_offsets[face - 1];
        if (total > std::numeric_limits<std::uint32_t>::max())
            throw std::invalid_argument("halfedge adjacency exceeds compact CSR capacity");
        adjacency_offsets[face] = static_cast<std::uint32_t>(total);
    }
    std::vector<std::uint32_t> adjacency(adjacency_offsets.back());
    std::vector<std::uint32_t> adjacency_cursor = adjacency_offsets;
    for (std::size_t edge = 0; edge < mesh.halfedges.size(); ++edge)
    {
        const std::uint32_t opposite = mesh.halfedges[edge].opposite;
        if (opposite == invalid_surface_index || opposite <= edge) continue;
        const std::uint32_t first = mesh.halfedges[edge].face;
        const std::uint32_t second = mesh.halfedges[opposite].face;
        if (first == second) continue;
        adjacency[adjacency_cursor[first]++] = second;
        adjacency[adjacency_cursor[second]++] = first;
    }
    adjacency_cursor.clear();
    adjacency_cursor.shrink_to_fit();

    const auto findActive = [&](std::uint32_t id)
    {
        std::uint32_t root = id;
        while (regions[root].replacement != root)
            root = regions[root].replacement;
        while (regions[id].replacement != id)
        {
            const std::uint32_t next = regions[id].replacement;
            regions[id].replacement = root;
            id = next;
        }
        return root;
    };
    const auto materializeNeighbors = [&](const std::uint32_t id)
    {
        std::vector<std::uint32_t> result;
        const Region& region = regions[id];
        if (region.merged)
        {
            result.reserve(region.merged->neighbors.size());
            for (const std::uint32_t neighbor : region.merged->neighbors)
            {
                const std::uint32_t active = findActive(neighbor);
                if (active != id) result.push_back(active);
            }
        }
        else
        {
            result.reserve(adjacency_offsets[id + 1] - adjacency_offsets[id]);
            for (std::uint32_t index = adjacency_offsets[id];
                 index < adjacency_offsets[id + 1]; ++index)
            {
                const std::uint32_t active = findActive(adjacency[index]);
                if (active != id) result.push_back(active);
            }
        }
        std::sort(result.begin(), result.end());
        result.erase(std::unique(result.begin(), result.end()), result.end());
        return result;
    };

    std::size_t maximum_queue_size = queue.size();
    std::cerr << "[memory] leaf_region_bytes=" << sizeof(Region)
              << " regions=" << regions.size()
              << " region_capacity=" << regions.capacity()
              << " csr_offsets=" << adjacency_offsets.size()
              << " csr_neighbors=" << adjacency.size()
              << " queue=" << queue.size() << '\n';

    std::size_t attempts = 0;
    while (!queue.empty())
    {
        const Pair pair = queue.front();
        queue.pop_front();
        Region& first = regions[pair.first];
        Region& second = regions[pair.second];
        // A replacement parent explicitly queues every unique current neighbor,
        // so stale child pairs are discarded rather than resolved and retried.
        if (first.replacement != pair.first || second.replacement != pair.second) continue;
        ++attempts;
        const auto attempt_time = std::chrono::steady_clock::now();
        if (attempt_time - last_progress >= std::chrono::seconds(10))
        {
            std::cerr << "[progress] attempts=" << attempts
                      << " accepted=" << stats.accepted_region_merges
                      << " active_regions=" << face_count - stats.accepted_region_merges
                      << " region_capacity=" << regions.capacity()
                      << " pending=" << queue.size()
                      << " max_pending=" << maximum_queue_size << '\n';
            last_progress = attempt_time;
        }

        const std::vector<std::uint32_t> first_neighbors = materializeNeighbors(pair.first);
        const std::vector<std::uint32_t> second_neighbors = materializeNeighbors(pair.second);
        std::vector<std::uint32_t> new_neighbors = mergedUnique(
            first_neighbors, second_neighbors);
        new_neighbors.erase(std::remove(new_neighbors.begin(), new_neighbors.end(), pair.first),
                            new_neighbors.end());
        new_neighbors.erase(std::remove(new_neighbors.begin(), new_neighbors.end(), pair.second),
                            new_neighbors.end());
        MergeCandidate candidate = makeMergeCandidate(
            mesh, pair.first, first, pair.second, second,
            options.numerical_tolerance);
        if (!candidate.orientation_consistent)
        {
            ++stats.rejected_orientation_merges;
            continue;
        }
        const std::size_t child_triangles =
            static_cast<std::size_t>(first.triangle_count) + second.triangle_count;
        const std::size_t candidate_triangles = candidate.polygon.triangles.size();
        if (candidate_triangles > child_triangles ||
            (candidate_triangles == child_triangles && new_neighbors.empty()))
        {
            ++stats.rejected_triangle_count_merges;
            continue;
        }

        const DirectedHausdorffCertificate certificate = certifier.certify(
            candidate.polygon, options.maximum_directed_hausdorff,
            options.maximum_certificate_depth, options.numerical_tolerance);
        if (!certificate.passed)
        {
            // The failed parent is discarded; its two certified children remain
            // active and therefore retain all source responsibility.
            ++stats.rejected_hausdorff_merges;
            continue;
        }

        const std::uint32_t parent_id = static_cast<std::uint32_t>(regions.size());
        Region parent;
        parent.replacement = parent_id;
        parent.triangle_count = static_cast<std::uint32_t>(candidate_triangles);
        parent.merged = std::make_unique<MergedRegionData>();
        parent.merged->certified_upper_bound = certificate.upper_bound;
        parent.merged->faces = std::move(candidate.faces);
        parent.merged->responsibility_vertices =
            std::move(candidate.responsibility_vertices);
        parent.merged->neighbors = std::move(new_neighbors);
        parent.merged->polygon_vertices = std::move(candidate.polygon.vertices);
        regions.push_back(std::move(parent));

        for (const std::uint32_t neighbor : regions[parent_id].merged->neighbors)
            enqueue(parent_id, neighbor);
        regions[pair.first].replacement = parent_id;
        regions[pair.second].replacement = parent_id;
        regions[pair.first].merged.reset();
        regions[pair.second].merged.reset();
        ++stats.accepted_region_merges;
        maximum_queue_size = std::max(maximum_queue_size, queue.size());
    }
    std::cerr << "[stage] region merge fixed point reached after " << attempts
              << " candidate attempts; regions=" << regions.size()
              << " active=" << face_count - stats.accepted_region_merges
              << " region_capacity=" << regions.capacity()
              << " max_queue=" << maximum_queue_size << "\n";

    MeshModel proxy;
    proxy.name = "phase2_halfedge_support_plane_proxy";
    std::map<std::array<double, 3>, std::uint32_t> vertex_map;
    for (std::uint32_t region_id = 0; region_id < regions.size(); ++region_id)
    {
        const Region& region = regions[region_id];
        if (region.replacement != region_id) continue;
        ++stats.planar_regions;
        if (!region.merged)
            ++stats.fallback_triangles;
        else
        {
            ++stats.emitted_planar_polygons;
            stats.certified_directed_hausdorff_upper_bound = std::max(
                stats.certified_directed_hausdorff_upper_bound,
                region.merged->certified_upper_bound);
        }
        appendRegion(mesh, region_id, region, proxy, vertex_map);
    }

    stats.final_vertices = proxy.vertices.size();
    stats.final_triangles = proxy.triangles.size();
    stats.coverage_certified = true;
    const DirectedHausdorffCertificate final_certificate = certifier.certify(
        proxy, options.maximum_directed_hausdorff,
        options.maximum_certificate_depth, options.numerical_tolerance);
    if (!final_certificate.passed)
        throw std::runtime_error(
            "assembled proxy failed its global directed Hausdorff certificate: " +
            final_certificate.failure_reason);
    stats.global_hausdorff_certified = true;
    stats.certified_directed_hausdorff_upper_bound = final_certificate.upper_bound;
    stats.global_hausdorff_reference_queries = final_certificate.reference_queries;
    stats.elapsed_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
    writeObj(output_directory / "proxy.obj", proxy);

    std::ofstream model(output_directory / "model.json");
    model << std::setprecision(17)
          << "{\"stats\":{\"phase1_triangles\":" << stats.phase1_triangles
          << ",\"proxy_triangles\":" << stats.final_triangles
          << ",\"regions\":" << stats.planar_regions
          << ",\"support_plane_polygons\":" << stats.emitted_planar_polygons
          << ",\"fallback_triangles\":" << stats.fallback_triangles
          << ",\"accepted_region_merges\":" << stats.accepted_region_merges
          << ",\"rejected_orientation_merges\":" << stats.rejected_orientation_merges
          << ",\"rejected_triangle_count_merges\":" << stats.rejected_triangle_count_merges
           << ",\"rejected_hausdorff_merges\":" << stats.rejected_hausdorff_merges
           << ",\"coverage_certified\":" << (stats.coverage_certified ? "true" : "false")
           << ",\"global_hausdorff_certified\":"
           << (stats.global_hausdorff_certified ? "true" : "false")
           << ",\"maximum_directed_hausdorff\":" << stats.maximum_directed_hausdorff
           << ",\"certified_directed_hausdorff_upper_bound\":"
           << stats.certified_directed_hausdorff_upper_bound
           << ",\"global_hausdorff_reference_queries\":"
           << stats.global_hausdorff_reference_queries
           << ",\"elapsed_seconds\":" << stats.elapsed_seconds << "}}\n";
    return stats;
}

} // namespace pqss_proxy_mesh
