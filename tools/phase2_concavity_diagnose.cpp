#include "pqss_proxy_mesh/hausdorff_simplifier.hpp"
#include "pqss_proxy_mesh/topology_fill.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

struct Vec3
{
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

Vec3 toVec(const pqss_proxy_mesh::Position3 p) { return {p.x, p.y, p.z}; }
pqss_proxy_mesh::Position3 toPosition(const Vec3 p) { return {p.x, p.y, p.z}; }
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
    return length > 0.0 ? a * (1.0 / length) : Vec3{};
}

struct FacePlane
{
    Vec3 normal;
    double offset = 0.0;
};

struct VertexConcavityExplanation
{
    std::uint32_t vertex = 0;
    std::uint32_t worst_face = 0;
    double worst_signed_distance = 0.0;
    double incident_min_dot = 1.0;
    std::uint32_t gathered_faces = 0;
};

Vec3 faceNormal(const pqss_proxy_mesh::MeshModel& mesh, const std::uint32_t face_id)
{
    const auto face = mesh.triangles[face_id];
    const Vec3 a = toVec(mesh.vertices[face[0]]);
    const Vec3 b = toVec(mesh.vertices[face[1]]);
    const Vec3 c = toVec(mesh.vertices[face[2]]);
    return normalized(cross(b - a, c - a));
}

void explainFaceDihedral(
    const pqss_proxy_mesh::OrientedSurfaceMesh& surface,
    const std::uint32_t face_id)
{
    const auto& mesh = surface.geometry;
    if (face_id >= mesh.triangles.size())
        throw std::out_of_range("explain face ID exceeds face count");
    std::vector<FacePlane> planes(mesh.triangles.size());
    for (std::uint32_t face = 0; face < mesh.triangles.size(); ++face)
    {
        const Vec3 normal = faceNormal(mesh, face);
        planes[face] = {normal, dot(normal, toVec(mesh.vertices[mesh.triangles[face][0]]))};
    }
    std::cout << std::setprecision(17) << "explain_face=" << face_id << '\n';
    std::uint32_t edge = surface.face_halfedges[face_id];
    for (int local = 0; local < 3; ++local)
    {
        const auto& halfedge = surface.halfedges[edge];
        const std::uint32_t next = halfedge.next;
        const std::uint32_t first = halfedge.origin;
        const std::uint32_t second = surface.halfedges[next].origin;
        const std::uint32_t opposite = halfedge.opposite;
        std::cout << "  halfedge=" << edge << " vertices=" << first << "->"
                  << second << " opposite=" << opposite;
        if (opposite < surface.halfedges.size())
        {
            const std::uint32_t opposite_face = surface.halfedges[opposite].face;
            const Vec3 direction = normalized(toVec(mesh.vertices[second]) -
                                              toVec(mesh.vertices[first]));
            const double signed_sine = dot(
                direction,
                cross(planes[face_id].normal, planes[opposite_face].normal));
            const double normal_dot =
                dot(planes[face_id].normal, planes[opposite_face].normal);
            std::cout << " opposite_face=" << opposite_face
                      << " signed_sine=" << signed_sine
                      << " normal_dot=" << normal_dot;
        }
        std::cout << '\n';
        edge = next;
    }
}

void writeObj(const std::filesystem::path& path, const pqss_proxy_mesh::MeshModel& mesh)
{
    std::ofstream stream(path);
    if (!stream) throw std::runtime_error("failed to write OBJ: " + path.string());
    stream << std::setprecision(17) << "o " << mesh.name << '\n';
    for (const auto point : mesh.vertices)
        stream << "v " << point.x << ' ' << point.y << ' ' << point.z << '\n';
    for (const auto face : mesh.triangles)
        stream << "f " << face[0] + 1 << ' ' << face[1] + 1 << ' ' << face[2] + 1 << '\n';
}

struct Diagnosis
{
    std::vector<std::uint8_t> concave_vertices;
    std::vector<std::pair<std::uint32_t, std::uint32_t>> concave_edges;
    std::vector<std::vector<std::uint32_t>> components;
    std::vector<std::vector<std::uint32_t>> component_faces;
    std::size_t highlighted_faces = 0;
    double seconds = 0.0;
};

Diagnosis diagnoseConcavity(
    const pqss_proxy_mesh::OrientedSurfaceMesh& surface,
    const std::uint32_t neighborhood_depth)
{
    const auto started = std::chrono::steady_clock::now();
    const auto& mesh = surface.geometry;
    if (mesh.vertices.empty() || mesh.triangles.empty())
        throw std::runtime_error("empty halfedge mesh");

    Vec3 lower = toVec(mesh.vertices.front());
    Vec3 upper = lower;
    for (const auto point : mesh.vertices)
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
    const double tolerance = std::max(diagonal * 1.0e-8, 1.0e-5);

    std::vector<std::vector<std::uint32_t>> vertex_faces(mesh.vertices.size());
    for (std::uint32_t face = 0; face < mesh.triangles.size(); ++face)
        for (const std::uint32_t vertex : mesh.triangles[face])
            vertex_faces[vertex].push_back(face);

    std::vector<std::vector<std::uint32_t>> vertex_neighbors(mesh.vertices.size());
    for (std::uint32_t face = 0; face < mesh.triangles.size(); ++face)
    {
        std::uint32_t edge = surface.face_halfedges[face];
        for (int local = 0; local < 3; ++local)
        {
            const auto& halfedge = surface.halfedges[edge];
            const std::uint32_t next = halfedge.next;
            const std::uint32_t destination = surface.halfedges[next].origin;
            vertex_neighbors[halfedge.origin].push_back(destination);
            vertex_neighbors[destination].push_back(halfedge.origin);
            edge = next;
        }
    }
    for (auto& neighbors : vertex_neighbors)
    {
        std::sort(neighbors.begin(), neighbors.end());
        neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
    }

    std::vector<FacePlane> planes(mesh.triangles.size());
    for (std::uint32_t face = 0; face < mesh.triangles.size(); ++face)
    {
        const Vec3 normal = faceNormal(mesh, face);
        planes[face] = {normal, dot(normal, toVec(mesh.vertices[mesh.triangles[face][0]]))};
    }

    std::vector<std::uint8_t> concave(mesh.vertices.size(), 0);
    std::vector<std::pair<std::uint32_t, std::uint32_t>> concave_edges;
    std::map<std::pair<std::uint32_t, std::uint32_t>,
             std::array<std::uint32_t, 2>>
        edge_faces;
    const double angle_tolerance = 1.0e-6;
    for (std::uint32_t edge = 0; edge < surface.halfedges.size(); ++edge)
    {
        const auto& halfedge = surface.halfedges[edge];
        const std::uint32_t opposite = halfedge.opposite;
        if (opposite >= surface.halfedges.size() || edge > opposite) continue;
        const auto& opposite_halfedge = surface.halfedges[opposite];
        const std::uint32_t face = halfedge.face;
        const std::uint32_t opposite_face = opposite_halfedge.face;
        if (face >= mesh.triangles.size() || opposite_face >= mesh.triangles.size())
            continue;
        const std::uint32_t first = halfedge.origin;
        const std::uint32_t second = surface.halfedges[halfedge.next].origin;
        if (first >= mesh.vertices.size() || second >= mesh.vertices.size() ||
            first == second)
            continue;
        const Vec3 direction = normalized(toVec(mesh.vertices[second]) -
                                          toVec(mesh.vertices[first]));
        const double signed_sine = dot(
            direction, cross(planes[face].normal, planes[opposite_face].normal));
        if (signed_sine < -angle_tolerance)
        {
            const auto key = std::minmax(first, second);
            concave[first] = 1;
            concave[second] = 1;
            concave_edges.push_back(key);
            edge_faces[key] = {face, opposite_face};
        }
    }
    std::sort(concave_edges.begin(), concave_edges.end());
    concave_edges.erase(std::unique(concave_edges.begin(), concave_edges.end()),
                        concave_edges.end());

    std::vector<std::uint8_t> visited(mesh.vertices.size(), 0);
    Diagnosis result;
    result.concave_vertices = std::move(concave);
    result.concave_edges = std::move(concave_edges);
    vertex_neighbors.assign(mesh.vertices.size(), {});
    for (const auto [first, second] : result.concave_edges)
    {
        vertex_neighbors[first].push_back(second);
        vertex_neighbors[second].push_back(first);
    }
    for (std::uint32_t seed = 0; seed < mesh.vertices.size(); ++seed)
    {
        if (!result.concave_vertices[seed] || visited[seed]) continue;
        std::vector<std::uint32_t> component{seed};
        visited[seed] = 1;
        for (std::size_t cursor = 0; cursor < component.size(); ++cursor)
            for (const std::uint32_t neighbor : vertex_neighbors[component[cursor]])
                if (result.concave_vertices[neighbor] && !visited[neighbor])
                {
                    visited[neighbor] = 1;
                    component.push_back(neighbor);
                }
        result.components.push_back(std::move(component));
    }

    result.component_faces.resize(result.components.size());
    for (std::size_t id = 0; id < result.components.size(); ++id)
    {
        const std::set<std::uint32_t> component_vertices(
            result.components[id].begin(), result.components[id].end());
        std::set<std::uint32_t> faces;
        for (const auto [first, second] : result.concave_edges)
            if (component_vertices.contains(first) && component_vertices.contains(second))
            {
                const auto iterator = edge_faces.find({first, second});
                if (iterator != edge_faces.end())
                {
                    faces.insert(iterator->second[0]);
                    faces.insert(iterator->second[1]);
                }
            }
        result.component_faces[id].assign(faces.begin(), faces.end());
        result.highlighted_faces += result.component_faces[id].size();
    }
    result.seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    return result;
}

VertexConcavityExplanation explainVertexConcavity(
    const pqss_proxy_mesh::OrientedSurfaceMesh& surface,
    const std::uint32_t vertex,
    const std::uint32_t neighborhood_depth)
{
    const auto& mesh = surface.geometry;
    if (vertex >= mesh.vertices.size())
        throw std::out_of_range("explain vertex ID exceeds vertex count");
    std::vector<std::vector<std::uint32_t>> vertex_faces(mesh.vertices.size());
    for (std::uint32_t face = 0; face < mesh.triangles.size(); ++face)
        for (const std::uint32_t corner : mesh.triangles[face])
            vertex_faces[corner].push_back(face);

    std::vector<std::vector<std::uint32_t>> face_neighbors(mesh.triangles.size());
    for (std::uint32_t face = 0; face < mesh.triangles.size(); ++face)
    {
        std::uint32_t edge = surface.face_halfedges[face];
        for (int local = 0; local < 3; ++local)
        {
            const auto& halfedge = surface.halfedges[edge];
            const std::uint32_t opposite = halfedge.opposite;
            if (opposite < surface.halfedges.size())
            {
                const std::uint32_t neighbor = surface.halfedges[opposite].face;
                if (neighbor < mesh.triangles.size())
                    face_neighbors[face].push_back(neighbor);
            }
            edge = halfedge.next;
        }
    }

    std::vector<FacePlane> planes(mesh.triangles.size());
    for (std::uint32_t face = 0; face < mesh.triangles.size(); ++face)
    {
        const Vec3 normal = faceNormal(mesh, face);
        planes[face] = {normal, dot(normal, toVec(mesh.vertices[mesh.triangles[face][0]]))};
    }

    std::vector<std::uint32_t> faces;
    std::vector<std::uint32_t> depths;
    std::vector<std::uint8_t> visited(mesh.triangles.size(), 0);
    for (const std::uint32_t face : vertex_faces[vertex])
        if (!visited[face])
        {
            visited[face] = 1;
            faces.push_back(face);
            depths.push_back(0);
        }
    for (std::size_t cursor = 0; cursor < faces.size(); ++cursor)
    {
        if (depths[cursor] >= neighborhood_depth) continue;
        for (const std::uint32_t neighbor : face_neighbors[faces[cursor]])
            if (!visited[neighbor])
            {
                visited[neighbor] = 1;
                faces.push_back(neighbor);
                depths.push_back(depths[cursor] + 1);
            }
    }

    const Vec3 point = toVec(mesh.vertices[vertex]);
    VertexConcavityExplanation result;
    result.vertex = vertex;
    result.gathered_faces = static_cast<std::uint32_t>(faces.size());
    result.worst_signed_distance = std::numeric_limits<double>::infinity();
    for (const std::uint32_t face : faces)
    {
        const double signed_distance =
            dot(planes[face].normal, point) - planes[face].offset;
        if (signed_distance < result.worst_signed_distance)
        {
            result.worst_signed_distance = signed_distance;
            result.worst_face = face;
        }
    }
    for (const std::uint32_t incident : vertex_faces[vertex])
        result.incident_min_dot = std::min(
            result.incident_min_dot,
            dot(planes[incident].normal, planes[result.worst_face].normal));
    return result;
}

void writeConcaveOverlayObj(
    const std::filesystem::path& path,
    const pqss_proxy_mesh::OrientedSurfaceMesh& surface,
    const Diagnosis& diagnosis)
{
    const auto& mesh = surface.geometry;
    std::ofstream stream(path);
    if (!stream) throw std::runtime_error("failed to write concavity overlay OBJ");
    stream << std::setprecision(17) << "o phase2_concave_components\n";
    std::uint32_t next_vertex = 1;
    for (std::size_t component = 0; component < diagnosis.component_faces.size(); ++component)
    {
        if (diagnosis.component_faces[component].empty()) continue;
        stream << "g region_" << component << "_component\n";
        for (const std::uint32_t face_id : diagnosis.component_faces[component])
        {
            const auto triangle = mesh.triangles[face_id];
            const Vec3 normal = faceNormal(mesh, face_id);
            for (const std::uint32_t vertex : triangle)
            {
                const Vec3 shifted = toVec(mesh.vertices[vertex]) + normal * 0.35;
                stream << "v " << shifted.x << ' ' << shifted.y << ' ' << shifted.z << '\n';
            }
            stream << "f " << next_vertex << ' ' << next_vertex + 1 << ' '
                   << next_vertex + 2 << '\n';
            next_vertex += 3;
        }
    }
}

void writeConcavePolylineObj(
    const std::filesystem::path& path,
    const pqss_proxy_mesh::OrientedSurfaceMesh& surface,
    const Diagnosis& diagnosis)
{
    const auto& mesh = surface.geometry;
    std::ofstream stream(path);
    if (!stream) throw std::runtime_error("failed to write concavity polyline OBJ");
    stream << std::setprecision(17) << "o phase2_concave_point_components\n";

    std::uint32_t next_vertex = 1;
    for (std::size_t component = 0; component < diagnosis.components.size(); ++component)
    {
        if (diagnosis.components[component].empty()) continue;
        stream << "g region_" << component << "_component\n";
        std::map<std::uint32_t, std::uint32_t> output_indices;
        for (const std::uint32_t vertex : diagnosis.components[component])
        {
            output_indices[vertex] = next_vertex++;
            const auto point = mesh.vertices[vertex];
            stream << "v " << point.x << ' ' << point.y << ' ' << point.z << '\n';
        }
        for (const std::uint32_t vertex : diagnosis.components[component])
            stream << "p " << output_indices[vertex] << '\n';
        std::set<std::pair<std::uint32_t, std::uint32_t>> lines;
        for (const auto [first, second] : diagnosis.concave_edges)
            if (output_indices.contains(first) && output_indices.contains(second))
                lines.insert({first, second});
        for (const auto [first, second] : lines)
            stream << "l " << output_indices[first] << ' ' << output_indices[second] << '\n';
    }
}

void writeModelJson(
    const std::filesystem::path& path, const Diagnosis& diagnosis,
    const std::uint32_t neighborhood_depth, const std::size_t source_triangles)
{
    std::ofstream stream(path);
    if (!stream) throw std::runtime_error("failed to write model.json");
    std::size_t concave_count = 0;
    for (const auto flag : diagnosis.concave_vertices) if (flag) ++concave_count;
    stream << std::setprecision(17)
           << "{\n  \"stats\":{\"source_triangles\":" << source_triangles
           << ",\"proxy_triangles\":" << diagnosis.highlighted_faces
           << ",\"primitive_count\":" << diagnosis.components.size()
           << ",\"primitive_types\":{\"component\":" << diagnosis.components.size()
           << "},\"timings_seconds\":{\"total\":" << diagnosis.seconds
           << "},\"phase2_concavity\":{\"neighborhood_depth\":"
           << neighborhood_depth << ",\"concave_vertices\":" << concave_count
           << ",\"concave_edges\":" << diagnosis.concave_edges.size()
           << ",\"components\":" << diagnosis.components.size()
           << ",\"highlighted_faces\":" << diagnosis.highlighted_faces << "}},\n"
           << "  \"source\":\"phase1.obj\",\n"
           << "  \"comparison_source_label\":\"阶段 1\",\n"
           << "  \"phase1_halfedge\":\"phase1_halfedge.bin\",\n"
           << "  \"phase3_simplified_surfaces\":\"phase2_concave_points.obj\",\n"
           << "  \"proxy\":\"phase2_concave_points.obj\",\n"
           << "  \"proxy_components\":[],\n"
           << "  \"viewer_stages\":[\"source\",\"phase1\",\"phase3\",\"split\"]\n}\n";
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        std::filesystem::path phase1_path;
        std::filesystem::path output_root;
        std::string model_id = "0";
        std::uint32_t neighborhood_depth = 4;
        std::vector<std::uint32_t> explain_vertices;
        std::vector<std::uint32_t> explain_faces;
        for (int index = 1; index < argc; ++index)
        {
            const std::string argument = argv[index];
            const auto value = [&]() -> std::string
            {
                if (++index >= argc) throw std::invalid_argument("missing value for " + argument);
                return argv[index];
            };
            if (argument == "--phase1-halfedge") phase1_path = value();
            else if (argument == "--output-root") output_root = value();
            else if (argument == "--model-id") model_id = value();
            else if (argument == "--neighborhood-depth")
                neighborhood_depth = static_cast<std::uint32_t>(std::stoul(value()));
            else if (argument == "--explain-vertex")
                explain_vertices.push_back(static_cast<std::uint32_t>(std::stoul(value())));
            else if (argument == "--explain-face")
                explain_faces.push_back(static_cast<std::uint32_t>(std::stoul(value())));
            else throw std::invalid_argument("unknown argument: " + argument);
        }
        if (phase1_path.empty() || output_root.empty())
            throw std::invalid_argument(
                "usage: pqss-phase2-concavity-diagnose --phase1-halfedge FILE "
                "--output-root DIR --model-id ID [--neighborhood-depth N]");
        const auto surface = pqss_proxy_mesh::readAnalysisHalfedgeMesh(phase1_path);
        for (const std::uint32_t face : explain_faces)
            explainFaceDihedral(surface, face);
        for (const std::uint32_t vertex : explain_vertices)
        {
            const auto explanation =
                explainVertexConcavity(surface, vertex, neighborhood_depth);
            std::cout << std::setprecision(17)
                      << "explain_vertex=" << explanation.vertex << '\n'
                      << "  gathered_faces=" << explanation.gathered_faces << '\n'
                      << "  worst_face=" << explanation.worst_face << '\n'
                      << "  worst_signed_distance="
                      << explanation.worst_signed_distance << '\n'
                      << "  incident_min_normal_dot_to_worst="
                      << explanation.incident_min_dot << '\n';
        }
        const auto diagnosis = diagnoseConcavity(surface, neighborhood_depth);
        const auto output_dir = output_root / "models" / model_id;
        std::filesystem::create_directories(output_dir);
        pqss_proxy_mesh::writeAnalysisHalfedgeMesh(
            output_dir / "phase1_halfedge.bin", surface);
        writeObj(output_dir / "phase1.obj", surface.geometry);
        writeConcaveOverlayObj(
            output_dir / "phase2_concave_components.obj", surface, diagnosis);
        writeConcavePolylineObj(
            output_dir / "phase2_concave_points.obj", surface, diagnosis);
        writeModelJson(
            output_dir / "model.json", diagnosis, neighborhood_depth,
            surface.geometry.triangles.size());
        std::ofstream manifest(output_root / "viewer_manifest.json");
        if (!manifest) throw std::runtime_error("failed to write viewer manifest");
        manifest << "{\n  \"algorithm\":\"CppPhase2ConcavityDiagnosisV1\",\n"
                 << "  \"complete\":true,\n  \"model_count\":1,\n"
                 << "  \"models\":[{\"id\":" << model_id
                 << ",\"metadata\":\"models/" << model_id << "/model.json\"}],\n"
                 << "  \"options\":{\"neighborhood_depth\":" << neighborhood_depth
                 << "}\n}\n";
        std::size_t concave_count = 0;
        for (const auto flag : diagnosis.concave_vertices) if (flag) ++concave_count;
        std::cout << "concave_vertices=" << concave_count << '\n'
                  << "components=" << diagnosis.components.size() << '\n'
                  << "highlighted_faces=" << diagnosis.highlighted_faces << '\n'
                  << "elapsed_seconds=" << diagnosis.seconds << '\n';
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
