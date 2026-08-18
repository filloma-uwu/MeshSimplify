#include "pqss_proxy_mesh/topology_fill.hpp"

#include <exception>
#include <filesystem>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv)
{
    try
    {
        if (argc != 3)
            throw std::invalid_argument("usage: pqss-obj-to-halfedge input.obj output.bin");
        const pqss_proxy_mesh::MeshModel mesh =
            pqss_proxy_mesh::readTriangleSoupObj(argv[1]);
        pqss_proxy_mesh::AnalysisHalfedgeStats stats;
        const pqss_proxy_mesh::OrientedSurfaceMesh halfedge =
            pqss_proxy_mesh::buildAnalysisHalfedgeMesh(mesh, &stats);
        pqss_proxy_mesh::writeAnalysisHalfedgeMesh(argv[2], halfedge);
        std::cout << "vertices=" << halfedge.geometry.vertices.size() << '\n'
                  << "faces=" << halfedge.geometry.triangles.size() << '\n'
                  << "halfedges=" << halfedge.halfedges.size() << '\n'
                  << "boundary_halfedges=" << stats.boundary_halfedges << '\n'
                  << "nonmanifold_edge_groups=" << stats.nonmanifold_edge_groups << '\n'
                  << "inconsistent_orientation_edges="
                  << stats.inconsistent_orientation_edges << '\n'
                  << "face_components=" << stats.face_components << '\n';
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
