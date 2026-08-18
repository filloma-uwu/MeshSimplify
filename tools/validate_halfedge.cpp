#include "pqss_proxy_mesh/halfedge_validation.hpp"
#include "pqss_proxy_mesh/hausdorff_simplifier.hpp"

#include <exception>
#include <iomanip>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv)
{
    try
    {
        if (argc != 2)
            throw std::invalid_argument(
                "usage: pqss-validate-halfedge phase1_halfedge.bin");
        const auto mesh = pqss_proxy_mesh::readAnalysisHalfedgeMesh(argv[1]);
        const auto report =
            pqss_proxy_mesh::validatePhase1HalfedgeTopology(mesh);
        std::cout << std::setprecision(17)
                  << "{\"valid\":true"
                  << ",\"validation\":\"combinatorial_halfedge_topology\""
                  << ",\"vertices\":" << report.vertices
                  << ",\"triangles\":" << report.triangles
                  << ",\"halfedges\":" << report.halfedges
                  << ",\"boundary_halfedges\":" << report.boundary_halfedges
                  << ",\"face_components\":" << report.face_components
                  << ",\"euler_characteristic\":"
                  << report.euler_characteristic
                  << ",\"signed_volume\":" << report.signed_volume
                  << ",\"geometric_intersections_tested\":false}\n";
        return 0;
    }
    catch (const std::exception& exception)
    {
        std::cerr << exception.what() << '\n';
        return 1;
    }
}
