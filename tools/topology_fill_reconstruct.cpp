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
        std::filesystem::path source_path;
        std::filesystem::path voxel_path;
        std::filesystem::path output_path;
        for (int index = 1; index < argc; ++index)
        {
            const std::string argument = argv[index];
            const auto value = [&]() -> std::filesystem::path
            {
                if (++index >= argc)
                    throw std::invalid_argument("missing value for " + argument);
                return argv[index];
            };
            if (argument == "--source") source_path = value();
            else if (argument == "--vox") voxel_path = value();
            else if (argument == "--output") output_path = value();
            else throw std::invalid_argument("unknown argument: " + argument);
        }
        if (source_path.empty() || voxel_path.empty() || output_path.empty())
            throw std::invalid_argument(
                "usage: pqss-topology-reconstruct --source source.obj "
                "--vox phase1.vox --output boundary.obj");

        auto occupancy = pqss_proxy_mesh::readVoxelGrid(voxel_path);
        const auto source = pqss_proxy_mesh::readTriangleSoupObj(source_path);
        const auto solid = pqss_proxy_mesh::buildPhase1Solid(
            std::move(occupancy), source);
        const auto source_occupancy = pqss_proxy_mesh::voxelizeTriangleSoupOnGrid(
            source, solid.occupancy);
        const auto projected = pqss_proxy_mesh::projectPhase1BoundaryToSource(
            solid.boundary, solid.occupancy, source_occupancy, source);
        pqss_proxy_mesh::writePhase1BoundaryObj(output_path, projected);
        std::cout << std::setprecision(17)
                  << "{\"mesh_vertices\":"
                  << projected.geometry.vertices.size()
                  << ",\"mesh_triangles\":"
                  << projected.geometry.triangles.size()
                  << ",\"mesh_euler_characteristic\":"
                  << projected.euler_characteristic
                  << ",\"mesh_signed_volume\":"
                  << projected.signed_volume << "}\n";
        return 0;
    }
    catch (const std::exception& exception)
    {
        std::cerr << exception.what() << '\n';
        return 1;
    }
}
