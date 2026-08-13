#include "pqss_proxy_mesh/topology_fill.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv)
{
    try
    {
        if (argc != 4)
            throw std::invalid_argument(
                "usage: pqss-source-dominated-preview source.obj phase1.vox output.obj");
        const auto source = pqss_proxy_mesh::readTriangleSoupObj(argv[1]);
        const auto occupancy = pqss_proxy_mesh::readVoxelGrid(argv[2]);
        const auto result = pqss_proxy_mesh::buildSourceDominatedPhase1Preview(
            occupancy, source);
        std::ofstream output(argv[3]);
        if (!output) throw std::runtime_error("failed to create output OBJ");
        output << std::setprecision(17) << "o source_dominated_phase1_preview\n";
        for (const auto point : result.vertices)
            output << "v " << point.x << ' ' << point.y << ' ' << point.z << '\n';
        for (const auto face : result.triangles)
            output << "f " << face[0] + 1 << ' ' << face[1] + 1 << ' '
                   << face[2] + 1 << '\n';
        std::cout << "vertices=" << result.vertices.size()
                  << " triangles=" << result.triangles.size() << '\n';
        return 0;
    }
    catch (const std::exception& exception)
    {
        std::cerr << exception.what() << '\n';
        return 1;
    }
}
