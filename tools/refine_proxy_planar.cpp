#include "pqss_proxy_mesh/hausdorff_simplifier.hpp"
#include "pqss_proxy_mesh/topology_fill.hpp"

#include <exception>
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <map>
#include <vector>

namespace
{
using Edge = std::pair<std::uint32_t, std::uint32_t>;

pqss_proxy_mesh::MeshModel capPlanarConvexBoundaryLoops(
    pqss_proxy_mesh::MeshModel mesh, std::size_t& capped_loops)
{
    std::map<Edge, std::size_t> edge_counts;
    std::map<Edge, Edge> directed_edges;
    for (const auto face : mesh.triangles)
        for (int local = 0; local < 3; ++local)
        {
            const std::uint32_t a = face[local];
            const std::uint32_t b = face[(local + 1) % 3];
            const Edge edge{std::min(a, b), std::max(a, b)};
            ++edge_counts[edge];
            directed_edges.try_emplace(edge, Edge{a, b});
        }
    std::map<std::uint32_t, std::vector<std::uint32_t>> adjacency;
    for (const auto& [edge, count] : edge_counts)
        if (count == 1)
        {
            adjacency[edge.first].push_back(edge.second);
            adjacency[edge.second].push_back(edge.first);
        }
    std::map<Edge, bool> visited;
    for (const auto& [start, neighbors] : adjacency)
    {
        for (const std::uint32_t initial_next : neighbors)
        {
            const Edge initial{std::min(start, initial_next), std::max(start, initial_next)};
            if (visited[initial]) continue;
            std::vector<std::uint32_t> loop{start};
            std::uint32_t previous = start;
            std::uint32_t current = initial_next;
            while (current != start && loop.size() <= adjacency.size())
            {
                loop.push_back(current);
                visited[{std::min(previous, current), std::max(previous, current)}] = true;
                const auto iterator = adjacency.find(current);
                if (iterator == adjacency.end() || iterator->second.size() != 2) break;
                const std::uint32_t next = iterator->second[0] == previous
                    ? iterator->second[1] : iterator->second[0];
                previous = current;
                current = next;
            }
            if (current != start || loop.size() < 3) continue;
            visited[{std::min(previous, start), std::max(previous, start)}] = true;
            std::size_t same_direction = 0;
            for (std::size_t index = 0; index < loop.size(); ++index)
            {
                const std::uint32_t a = loop[index];
                const std::uint32_t b = loop[(index + 1) % loop.size()];
                const Edge edge{std::min(a, b), std::max(a, b)};
                const auto directed = directed_edges.find(edge);
                if (directed != directed_edges.end() && directed->second == Edge{a, b})
                    ++same_direction;
            }
            if (same_direction * 2 > loop.size()) std::reverse(loop.begin(), loop.end());
            std::array<double, 3> normal{};
            for (std::size_t index = 0; index < loop.size(); ++index)
            {
                const auto a = mesh.vertices[loop[index]];
                const auto b = mesh.vertices[loop[(index + 1) % loop.size()]];
                normal[0] += (a.y - b.y) * (a.z + b.z);
                normal[1] += (a.z - b.z) * (a.x + b.x);
                normal[2] += (a.x - b.x) * (a.y + b.y);
            }
            const double length = std::sqrt(
                normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
            if (!(length > 1.0e-12)) continue;
            for (double& value : normal) value /= length;
            const auto origin = mesh.vertices[loop[0]];
            bool planar = true;
            bool convex = true;
            double turn_sign = 0.0;
            for (std::size_t index = 0; index < loop.size(); ++index)
            {
                const auto point = mesh.vertices[loop[index]];
                const double plane_distance = normal[0] * (point.x - origin.x) +
                    normal[1] * (point.y - origin.y) + normal[2] * (point.z - origin.z);
                planar &= std::abs(plane_distance) <= 1.0e-7;
                const auto a = mesh.vertices[loop[index]];
                const auto b = mesh.vertices[loop[(index + 1) % loop.size()]];
                const auto c = mesh.vertices[loop[(index + 2) % loop.size()]];
                const double ab[3]{b.x-a.x,b.y-a.y,b.z-a.z};
                const double bc[3]{c.x-b.x,c.y-b.y,c.z-b.z};
                const double turn = normal[0] * (ab[1]*bc[2]-ab[2]*bc[1]) +
                    normal[1] * (ab[2]*bc[0]-ab[0]*bc[2]) +
                    normal[2] * (ab[0]*bc[1]-ab[1]*bc[0]);
                if (std::abs(turn) > 1.0e-12)
                {
                    if (turn_sign == 0.0) turn_sign = turn;
                    else convex &= turn * turn_sign > 0.0;
                }
            }
            if (!planar || !convex) continue;
            for (std::uint32_t index = 1; index + 1 < loop.size(); ++index)
                mesh.triangles.push_back({loop[0], loop[index], loop[index + 1]});
            ++capped_loops;
        }
    }
    return mesh;
}
}

int main(int argc, char** argv)
{
    try
    {
        std::filesystem::path proxy_path;
        std::filesystem::path phase1_path;
        std::filesystem::path output_path;
        double limit = -1.0;
        for (int index = 1; index < argc; ++index)
        {
            const std::string argument = argv[index];
            const auto value = [&]() -> std::string
            {
                if (++index >= argc) throw std::invalid_argument("missing value for " + argument);
                return argv[index];
            };
            if (argument == "--proxy") proxy_path = value();
            else if (argument == "--phase1-halfedge") phase1_path = value();
            else if (argument == "--output") output_path = value();
            else if (argument == "--maximum-distance") limit = std::stod(value());
            else throw std::invalid_argument("unknown argument: " + argument);
        }
        if (proxy_path.empty() || phase1_path.empty() || output_path.empty() || limit < 0.0)
            throw std::invalid_argument(
                "usage: pqss-refine-proxy-planar --proxy FILE --phase1-halfedge FILE "
                "--output FILE --maximum-distance VALUE");
        const pqss_proxy_mesh::MeshModel proxy =
            pqss_proxy_mesh::readTriangleSoupObj(proxy_path);
        const pqss_proxy_mesh::OrientedSurfaceMesh phase1 =
            pqss_proxy_mesh::readAnalysisHalfedgeMesh(phase1_path);
        pqss_proxy_mesh::MeshModel refined =
            pqss_proxy_mesh::convexifyCertifiedPlanarComponents(
                proxy, phase1.geometry, limit);
        std::size_t capped_loops = 0;
        refined = capPlanarConvexBoundaryLoops(std::move(refined), capped_loops);
        const auto certificate = pqss_proxy_mesh::certifyDirectedHausdorff(
            refined, phase1.geometry, limit);
        if (!certificate.passed)
            throw std::runtime_error("refined proxy failed global Hausdorff certification");
        std::ofstream output(output_path);
        if (!output) throw std::runtime_error("failed to create output OBJ");
        output << std::setprecision(17) << "o certified_planar_refined_proxy\n";
        for (const auto point : refined.vertices)
            output << "v " << point.x << ' ' << point.y << ' ' << point.z << '\n';
        for (const auto face : refined.triangles)
            output << "f " << face[0] + 1 << ' ' << face[1] + 1 << ' ' << face[2] + 1 << '\n';
        std::cout << std::setprecision(17)
                  << "input_triangles=" << proxy.triangles.size() << '\n'
                  << "output_triangles=" << refined.triangles.size() << '\n'
                  << "capped_planar_convex_loops=" << capped_loops << '\n'
                  << "certified_upper_bound=" << certificate.upper_bound << '\n';
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
