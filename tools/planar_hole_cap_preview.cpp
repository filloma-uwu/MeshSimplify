#include "pqss_proxy_mesh/topology_fill.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv)
{
    try
    {
        if (argc != 3)
            throw std::invalid_argument(
                "usage: pqss-planar-hole-cap-preview input.obj output.obj");
        const auto source = pqss_proxy_mesh::readTriangleSoupObj(argv[1]);
        const auto occupancy = pqss_proxy_mesh::voxelizeTriangleSoup(
            source, 12'000'000, 4);
        const auto filled = pqss_proxy_mesh::enclosingTopologyFill(
            occupancy, 32);
        const auto preview = pqss_proxy_mesh::buildPlanarHoleCapPreview(
            filled, source);

        const std::filesystem::path output_path = argv[2];
        std::ofstream output(output_path);
        if (!output)
            throw std::runtime_error(
                "failed to create preview OBJ: " + output_path.string());
        output << std::setprecision(17) << "o planar_hole_cap_preview\n";
        for (const auto point : preview.vertices)
            output << "v " << point.x << ' ' << point.y << ' ' << point.z << '\n';
        for (const auto face : preview.triangles)
            output << "f " << face[0] + 1 << ' ' << face[1] + 1 << ' '
                   << face[2] + 1 << '\n';
        std::cout << "source_triangles=" << source.triangles.size() << '\n'
                  << "preview_triangles=" << preview.triangles.size() << '\n';
        return 0;
    }
    catch (const std::exception& exception)
    {
        std::cerr << exception.what() << '\n';
        return 1;
    }
}
