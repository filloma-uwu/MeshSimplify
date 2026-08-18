#include "pqss_proxy_mesh/hausdorff_simplifier.hpp"
#include "pqss_proxy_mesh/topology_fill.hpp"

#include <exception>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv)
{
    try
    {
        std::filesystem::path proxy_path;
        std::filesystem::path phase1_path;
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
            else if (argument == "--maximum-distance") limit = std::stod(value());
            else throw std::invalid_argument("unknown argument: " + argument);
        }
        if (proxy_path.empty() || phase1_path.empty() || limit < 0.0)
            throw std::invalid_argument(
                "usage: pqss-certify-directed-hausdorff --proxy FILE "
                "--phase1-halfedge FILE --maximum-distance VALUE");
        const pqss_proxy_mesh::MeshModel proxy =
            pqss_proxy_mesh::readTriangleSoupObj(proxy_path);
        pqss_proxy_mesh::AnalysisHalfedgeStats topology;
        (void)pqss_proxy_mesh::buildAnalysisHalfedgeMesh(proxy, &topology);
        const pqss_proxy_mesh::OrientedSurfaceMesh phase1 =
            pqss_proxy_mesh::readAnalysisHalfedgeMesh(phase1_path);
        const auto certificate = pqss_proxy_mesh::certifyDirectedHausdorff(
            proxy, phase1.geometry, limit);
        std::cout << std::setprecision(17)
                  << "{\"passed\":" << (certificate.passed ? "true" : "false")
                  << ",\"proxy_triangles\":" << proxy.triangles.size()
                  << ",\"boundary_halfedges\":" << topology.boundary_halfedges
                  << ",\"nonmanifold_edge_groups\":" << topology.nonmanifold_edge_groups
                  << ",\"inconsistent_orientation_edges\":"
                  << topology.inconsistent_orientation_edges
                  << ",\"face_components\":" << topology.face_components
                  << ",\"lower_bound\":" << certificate.lower_bound
                  << ",\"upper_bound\":" << certificate.upper_bound
                  << ",\"reason\":\"" << certificate.failure_reason << "\"}\n";
        return certificate.passed ? 0 : 1;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 2;
    }
}
